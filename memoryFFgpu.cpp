// memory model
// field-free case

#define EIGEN_USE_BLAS 
#define EIGEN_USE_LAPACKE

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include <cxxopts.hpp>
#include "cnpy.h"
#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>

// full k-step propagator with time step h,
// using a strictly diagonal and real hamiltonian (hence VectorXd)
// this is a "full TDCI density matrix" propagator
Eigen::MatrixXcd prop(int k, double h, const Eigen::VectorXd& hamiltonian)
{
  using namespace std::complex_literals;
  std::complex<double> coeff = -1.0i * static_cast<double>(k) * h;
  Eigen::VectorXcd factor = coeff*hamiltonian.array();
  Eigen::VectorXcd diagvals = factor.array().exp();
  return diagvals.asDiagonal();
}

// reduced k-step propagator with time step h,
// using a strictly diagonal and real hamiltonian (hence VectorXd)
// by "reduced" what is meant is that only goodCols are retained
// however this is still a "full TDCI density matrix" propagator
Eigen::MatrixXcd redprop(int k, double h, const Eigen::VectorXd& hamiltonian, const std::vector<int>& redCols)
{
  Eigen::MatrixXcd thisprop = prop(k, h, hamiltonian);
  Eigen::MatrixXcd hamkron = Eigen::kroneckerProduct(thisprop, thisprop.conjugate());
  Eigen::VectorXcd redpropdiag(redCols.size());
  for (size_t j=0; j<redCols.size(); ++j) 
    redpropdiag(j) = hamkron(redCols[j], redCols[j]);
  return redpropdiag.asDiagonal();
}

// Function to compute the Moore–Penrose pseudoinverse
Eigen::MatrixXcd pseudoInverse(const Eigen::MatrixXcd& input, double tol) 
{
    int m = input.rows();
    int n = input.cols();
    int lda = m;
    int min_mn = std::min(m, n);

    // 1. Initialize cuSOLVER
    cusolverDnHandle_t cusolverH = nullptr;
    cusolverDnCreate(&cusolverH);

    // 2. Allocate Managed Memory (accessible by both CPU and GPU)
    // cuDoubleComplex is CUDA's version of std::complex<double>
    cuDoubleComplex *d_A, *d_U, *d_VT;
    double *d_S, *d_rwork;
    int *devInfo;

    cudaMallocManaged(&d_A, m * n * sizeof(cuDoubleComplex));
    cudaMallocManaged(&d_U, m * min_mn * sizeof(cuDoubleComplex));   // Thin U
    cudaMallocManaged(&d_VT, min_mn * n * sizeof(cuDoubleComplex));  // Thin V^T
    cudaMallocManaged(&d_S, min_mn * sizeof(double));                // Singular values
    cudaMallocManaged(&devInfo, sizeof(int));
    cudaMallocManaged(&d_rwork, 5 * min_mn * sizeof(double));        // Real workspace for ZGESVD

    // Copy input data directly into the managed buffer (zero translation needed)
    memcpy(d_A, input.data(), m * n * sizeof(cuDoubleComplex));

    // 3. Query cuSOLVER for required workspace size
    int lwork = 0;
    // 'S' requests the "Thin" SVD, which is vastly faster for tall/skinny matrices
    cusolverDnZgesvd_bufferSize(cusolverH, m, n, &lwork);
  
    cuDoubleComplex *d_work;
    cudaMallocManaged(&d_work, lwork * sizeof(cuDoubleComplex));

    // 4. Fire off the SVD on the GPU!
    cusolverDnZgesvd(
        cusolverH, 'S', 'S', 
        m, n, 
        d_A, lda, 
        d_S, 
        d_U, lda, 
        d_VT, min_mn, 
        d_work, lwork, 
        d_rwork, devInfo
    );
    
    // Wait for the GPU to finish before the CPU reads the results
    cudaDeviceSynchronize();

    std::cout << "Done with SVD on GPU!\n";

    // 5. Wrap the managed pointers in Eigen Maps for easy math
    // reinterpret_cast is perfectly safe here because std::complex<double> and cuDoubleComplex have identical memory layouts
    Eigen::Map<Eigen::MatrixXcd> VT(reinterpret_cast<std::complex<double>*>(d_VT), min_mn, n);

    // 6. Scale U by the inverse singular values
    Eigen::Map<Eigen::MatrixXcd> U(reinterpret_cast<std::complex<double>*>(d_U), m, min_mn);
    Eigen::Map<Eigen::VectorXd> S(d_S, min_mn);
    
    Eigen::VectorXd S_inv(min_mn);
    for (int i = 0; i < min_mn; ++i) {
        S_inv(i) = (S(i) > tol) ? 1.0 / S(i) : 0.0;
    }

    // Scale U directly in Unified Memory (This is an O(N^2) memory operation, so the CPU handles it in <1 second)
    U = U * S_inv.asDiagonal();

    // 7. Reconstruct Pseudoinverse on GPU using cuBLAS
    // We compute P^H = U_scaled * d_VT
    cuDoubleComplex *d_Ph;
    cudaMallocManaged(&d_Ph, m * n * sizeof(cuDoubleComplex));

    cublasHandle_t cublasH;
    cublasCreate(&cublasH);

    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta = {0.0, 0.0};

    // cublasZgemm performs C = alpha * A * B + beta * C
    // A = d_U (m x min_mn), B = d_VT (min_mn x n), C = d_Ph (m x n)
    cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_N, 
                m, n, min_mn, 
                &alpha, 
                d_U, m, 
                d_VT, min_mn, 
                &beta, 
                d_Ph, m);

    cudaDeviceSynchronize();

    std::cout << "Done with cublasZgemm!\n";

    // Map the GPU result back to Eigen and take the adjoint to get P
    Eigen::Map<Eigen::MatrixXcd> Ph(reinterpret_cast<std::complex<double>*>(d_Ph), m, n);
    Eigen::MatrixXcd Pinv = Ph.adjoint();

    // 8. Clean up VRAM
    cudaFree(d_A); cudaFree(d_U); cudaFree(d_VT); cudaFree(d_S); 
    cudaFree(d_work); cudaFree(d_rwork); cudaFree(devInfo); cudaFree(d_Ph);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);
    
    // NRVO kicks in, no copy is made on return!
    return Pinv;
}

// 1RDM propagator with memory length n (time steps) and time step h
Eigen::MatrixXcd qprop(int n, double h, const Eigen::VectorXd& hamiltonian, const std::vector<int>& redCols, const Eigen::MatrixXd& BmatT, double tol=1e-6)
{
  int drc2 = BmatT.rows();

  Eigen::MatrixXcd bigmat((n+1)*drc2, redCols.size());
  for (int j=0; j<=n; ++j)
    bigmat.block(j*drc2, 0, drc2, redCols.size()) = BmatT * redprop(-j, h, hamiltonian, redCols);
    
  Eigen::MatrixXcd bigmatpinv = pseudoInverse(bigmat, tol);
  std::cout << "Computed pseudoInverse!\n";
  Eigen::MatrixXcd thisredprop = redprop(1, h, hamiltonian, redCols);
  std::cout << "Computed thisredprop!\n";
  return BmatT * thisredprop * bigmatpinv;
}

int main(int argc, char** argv)
{
  // omp_set_max_active_levels(1); 
  int num_threads = omp_get_max_threads();
  std::cout << "num_threads = " << num_threads << "\n";
  // omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  
  cxxopts::Options options("memoryFF", "Field-free memory model for 1RDM propagation");

  options.add_options()
  ("h,help", "Print usage")
  ("dt", "Time step size", cxxopts::value<double>())
  ("delay", "Integer value of delay", cxxopts::value<int>())
  ("infile", "Input file name", cxxopts::value<std::string>())
  ("savemae", "Path at which to save MAE to disk", cxxopts::value<std::string>())
  ("verbose", "Print lots of information to screen", cxxopts::value<bool>()->default_value("false"))
  ("saveqprop", "Save qprop matrix to disk", cxxopts::value<bool>()->default_value("false"))
  ("tol", "SVD tolerance", cxxopts::value<double>()->default_value("1e-6"));
  
  auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0;
  }

  double dt, tol;
  int delay;
  std::string infile, outpath;
  bool verbose, savemae, saveqprop;

  if (result.count("dt")==0)
  {
    std::cout << "Must specify dt!\n";
    return 1;
  }
  else
    dt = result["dt"].as<double>();
  if (result.count("delay")==0)
  {
    std::cout << "Must specify delay!\n";
    return 1;
  }
  else 
    delay = result["delay"].as<int>();
  if (result.count("infile")==0)
  {
    std::cout << "Must specify input file!\n";
    return 1;
  }
  else
    infile = result["infile"].as<std::string>();
  if (result.count("tol")==0)
    tol = 1e-6;
  else
    tol = result["tol"].as<double>();
  if (result.count("verbose")==0)
    verbose = false;
  else
    verbose = true;    
  if (result.count("savemae")==0)
  {
    savemae = false;
  }
  else
  {
    savemae = true;
    outpath = result["savemae"].as<std::string>();
  }
  if (result.count("saveqprop")==0)
    saveqprop = false;
  else
    saveqprop = true;

  // Load the entire .npz file into a map-like structure
  cnpy::npz_t my_npz = cnpy::npz_load(infile);

  // Load ham
  cnpy::NpyArray arr = my_npz["ham"];
  double* hamdata = arr.data<double>();
  size_t length = arr.shape[0];
  Eigen::Map<Eigen::VectorXd> ham(hamdata, length);
  int drcCI = (int) length;
  int drcCI2 = drcCI*drcCI;
  if (verbose)
    std::cout << "drcCI = " << drcCI << "\n";

  // Load Bten  
  arr = my_npz["Bten"];
  double* Btendata = arr.data<double>();
  length = arr.shape[0]*arr.shape[1]*arr.shape[2]*arr.shape[3];
  int drc = (int) arr.shape[2];
  int drc2 = drc*drc;
  if (((verbose)))
    std::cout << "drc = " << drc << "\n";
  Eigen::Map<Eigen::VectorXd> Bten(Btendata, length);
  Eigen::Map<const Eigen::MatrixXd> BmatT(Bten.data(), drc2, drcCI2);

  Eigen::VectorXi diagcols(drcCI);
  for (int k=0; k<drcCI; ++k) 
    diagcols(k) = k * drcCI + k; 

  // reduce the number of columns of BmatT
  Eigen::VectorXd colnorms(drcCI2);
  std::vector<int> goodCols;
  goodCols.reserve(drcCI2);
  int diagcheck = 0;
  Eigen::VectorXd bvec = Eigen::VectorXd::Zero(drc2);
  for (int j=0; j<drcCI2; ++j)
  {
    // skip column whose index corresponds to diagonal entry of full density
    if (j==diagcols(diagcheck))
    {
      bvec += BmatT.col(j)/drcCI;
      diagcheck++;
      continue;
    }
    // retain columns only if their norm exceeds a threshold
    colnorms(j) = BmatT.col(j).norm();
    if (colnorms[j] > 1e-14)
      goodCols.push_back(j);
  }
  if (verbose)
    std::cout << "number of good cols = " << goodCols.size() << "\n";

  // create new matrix with selected columns
  Eigen::MatrixXd BmatTgood(drc2, goodCols.size());
  for (size_t j=0; j<goodCols.size(); ++j)
    BmatTgood.col(j) = BmatT.col(goodCols[j]);

  if (verbose)
  {
    std::cout << "BmatTgood # of rows = " << BmatTgood.rows() << "\n";
    std::cout << "BmatTgood # of cols = " << BmatTgood.cols() << "\n";
  }
  Eigen::MatrixXcd thisqprop = qprop(delay, dt, ham, goodCols, BmatTgood, tol);
  if (verbose)
  {
    std::cout << "thisqprop # of rows = " << thisqprop.rows() << "\n";
    std::cout << "thisqprop # of cols = " << thisqprop.cols() << "\n";
  }
  // save qprop to disk
  if (saveqprop)
  {
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"); 
    std::string name = "thisqpropR.csv";
    std::ofstream fileR(name.c_str());
    fileR << thisqprop.real().matrix().format(CSVFormat);
    name = "thisqpropI.csv";
    std::ofstream fileI(name.c_str());
    fileI << thisqprop.imag().matrix().format(CSVFormat);
  }

  // initialize coeff matrix and set initial condition
  double T = 200.0;
  int nsteps = static_cast<int>(std::ceil(T/dt));
  if (verbose)
    std::cout << "About to propagate full TDCI coefficients for " << nsteps << " steps\n";

  Eigen::MatrixXcd coeffs(drcCI, nsteps+1);
  for (int j=0; j<drcCI; ++j)
    coeffs(j, 0) = 1.0/std::sqrt(drcCI);

  // propagate TDCI coefficients
  Eigen::MatrixXcd fullprop = prop(1, dt, ham);
  for (int k=0; k<nsteps; ++k)
  {
    coeffs.col(k+1) = fullprop * coeffs.col(k);
  }

  if (verbose)
    std::cout << "Norm of solution at final time = " << coeffs.col(nsteps).norm() << "\n";

  // compute ground truth 1RDMs
  Eigen::MatrixXcd true1rdms(drc2, nsteps+1);
  Eigen::MatrixXcd op;
  for (int k=0; k<=nsteps; ++k)
  {
    // outer product
    op = coeffs.col(k) * coeffs.col(k).adjoint();
    // reduction
    true1rdms.col(k) = BmatT * op.transpose().reshaped();
  }

  // set up predicted 1RDMs
  Eigen::MatrixXcd pred1rdms(drc2, nsteps+1);
  if (verbose)
    std::cout << "About to propagate 1RDMs for " << nsteps << " steps\n";
  for (int k=0; k<=delay; ++k)
    pred1rdms.col(k) = true1rdms.col(k);
  
  for (int k=delay; k<nsteps; ++k)
  {
    Eigen::MatrixXcd temp = pred1rdms.block(0, k-delay, drc2, delay+1).rowwise().reverse();
    Eigen::MatrixXcd temp2 = (temp.colwise() - bvec).reshaped();
    pred1rdms.col(k+1) = bvec + thisqprop * temp2;
  }
  
  double mae = (true1rdms - pred1rdms).array().abs().mean();
  if (verbose)
    std::cout << "Mean Absolute Error: " << mae << "\n";

  std::ostringstream oss;
  oss << std::setprecision(17) << mae;
  std::string maestr = oss.str();
  if (savemae)
  {
    std::filesystem::path p(infile);
    std::string stem = p.stem().string();
    std::string filename = stem + "_" + std::to_string(dt) + "_" + std::to_string(delay) + ".txt";
    std::filesystem::path dir(outpath);
    std::filesystem::path outfile = dir / filename;
    std::ofstream out(outfile);
    out << maestr << "\n";
  }
  else
    std::cout << maestr << "\n";

  return 0;
}


