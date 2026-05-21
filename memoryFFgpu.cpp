// memory model
// field-free case

#define EIGEN_USE_MKL

#include <omp.h>
#include <mkl.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cuComplex.h>

#include <Eigen/Core>
#include <Eigen/Dense>
// #include <unsupported/Eigen/KroneckerProduct>

#include <cxxopts.hpp>
#include "cnpy.h"
#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <utility>

using namespace std::complex_literals;

__global__ void scale_U_kernel(cuDoubleComplex* U, const double* S, int m, int min_mn, double tol) 
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = m * min_mn;
    
    if (tid < total_elements) 
    {
        int c = tid / m; // Column index
        double s_val = S[c];
        
        // Use the tolerance to find the inverse
        double sinv_val = (s_val > tol) ? 1.0 / s_val : 0.0;
        
        // Scale the complex value
        cuDoubleComplex u_val = U[tid];
        u_val.x *= sinv_val;
        u_val.y *= sinv_val;
        U[tid] = u_val;
    }
}

__global__ void add_bvec_kernel(cuDoubleComplex* d_target_col, const cuDoubleComplex* d_bvec, int drc2) 
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < drc2) 
    {
        cuDoubleComplex target_val = d_target_col[tid];
        cuDoubleComplex b_val = d_bvec[tid];
        
        // Add bvec to the matrix multiplication result
        target_val.x += b_val.x;
        target_val.y += b_val.y;
        
        d_target_col[tid] = target_val;
    }
}

__global__ void prep_temp2_kernel(
    const cuDoubleComplex* __restrict__ d_pred1rdms,
    const cuDoubleComplex* __restrict__ d_bvec,
    cuDoubleComplex* __restrict__ d_temp2,
    int k, int delay, int drc2)
{
    // Total number of elements we need to process
    int total_elements = drc2 * (delay + 1);
    
    // Get the global thread ID
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid < total_elements)
    {
        // Because Eigen is Column-Major, we can extract the local 2D coordinates
        int c = tid / drc2; // Integer division gives the "column" chunk
        int r = tid % drc2; // Modulo gives the row index
        
        // The reverse() logic: map to the source column in pred1rdms
        int source_col = k - c;
        
        // Calculate flat indices
        int source_idx = source_col * drc2 + r; // pred1rdms has a leading dimension of drc2
        
        // Fetch values
        cuDoubleComplex val = d_pred1rdms[source_idx];
        cuDoubleComplex bval = d_bvec[r];
        
        // cuDoubleComplex doesn't support the raw '-' operator in device code, so we use cuCsub
        d_temp2[tid] = cuCsub(val, bval);
    }
}

// full k-step propagator with time step h,
// using a strictly diagonal and real hamiltonian (hence VectorXd)
// this is a "full TDCI density matrix" propagator
Eigen::MatrixXcd prop(int k, double h, const Eigen::VectorXd& hamiltonian)
{
  std::complex<double> coeff = -1.0i * static_cast<double>(k) * h;
  Eigen::VectorXcd factor = coeff*hamiltonian.array();
  Eigen::VectorXcd diagvals = factor.array().exp();
  return diagvals.asDiagonal();
}

// reduced k-step propagator with time step h,
// using a strictly diagonal and real hamiltonian (hence VectorXd)
// by "reduced" what is meant is that only goodCols are retained
// however this is still a "full TDCI density matrix" propagator
Eigen::VectorXcd redprop(int k, double h, const Eigen::VectorXd& hamiltonian, const std::vector<int>& redCols)
{
  std::complex<double> coeff = -1.0i * static_cast<double>(k) * h;
  Eigen::VectorXcd factor = coeff*hamiltonian.array();
  Eigen::VectorXcd diagvals = factor.array().exp();
  Eigen::MatrixXcd outer_product = diagvals * diagvals.adjoint();
  Eigen::VectorXcd kronprop = outer_product.reshaped<Eigen::RowMajor>();
  return kronprop(redCols);
}

// Function to compute the bare metal ingredients to compute the Moore–Penrose pseudoinverse
std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> pseudoInverse(const Eigen::MatrixXcd& input, double tol) 
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
    size_t m_sz = m;
    size_t n_sz = n;
    size_t min_mn_sz = min_mn;

    cuDoubleComplex *d_A, *d_U, *d_VT;
    double *d_S, *d_rwork;
    int *devInfo;

    cudaMallocManaged(&d_A, m_sz * n_sz * sizeof(cuDoubleComplex));
    cudaMallocManaged(&d_U, m_sz * min_mn_sz * sizeof(cuDoubleComplex));   
    cudaMallocManaged(&d_VT, min_mn_sz * n_sz * sizeof(cuDoubleComplex));  
    cudaMallocManaged(&d_S, min_mn_sz * sizeof(double));                
    cudaMallocManaged(&devInfo, sizeof(int));
    cudaMallocManaged(&d_rwork, 5 * min_mn * sizeof(double));        // Real workspace for ZGESVD
    
    cudaMemcpy(d_A, input.data(), m_sz * n_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

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

    // Safely scale the columns of U in-place using the custom GPU kernel we discussed earlier
    int threads = 256;
    int blocks = (m_sz * min_mn_sz + threads - 1) / threads;
    scale_U_kernel<<<blocks, threads>>>(d_U, d_S, m, min_mn, tol);
    cudaDeviceSynchronize();
    std::cout << "Columns of U have been scaled in-place on GPU!\n";

    Eigen::MatrixXcd U_ret(m, min_mn);
    Eigen::MatrixXcd VT_ret(min_mn, n);

    // Explicitly copy the data from the GPU directly into the Eigen buffers
    cudaMemcpy(U_ret.data(), d_U, m * min_mn * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(VT_ret.data(), d_VT, min_mn * n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    // Clean up VRAM safely
    cudaFree(d_A); cudaFree(d_U); cudaFree(d_VT); cudaFree(d_S); 
    cudaFree(d_work); cudaFree(d_rwork); cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    
    return {U_ret, VT_ret};
}

// 1RDM propagator with memory length n (time steps) and time step h
Eigen::MatrixXcd qprop(int n, double h, const Eigen::VectorXd& hamiltonian, const std::vector<int>& redCols, const Eigen::MatrixXd& Bmat, double tol=1e-6)
{
  int drc2 = Bmat.rows();
  int rCs = redCols.size();
  Eigen::MatrixXcd bigmat((n+1)*drc2, rCs);
  Eigen::VectorXcd base_prop_vec = redprop(-1, h, hamiltonian, redCols);
  Eigen::VectorXcd current_prop_vec = Eigen::VectorXcd::Ones(rCs);
  for (int j=0; j<=n; ++j)
  {
    bigmat.block(j * drc2, 0, drc2, rCs).noalias() = Bmat * current_prop_vec.asDiagonal();
    current_prop_vec.array() *= base_prop_vec.array();
  }
  auto [U_scaled, VT] = pseudoInverse(bigmat, tol);
  std::cout << "Computed pseudoinverse!\n";
  Eigen::MatrixXcd thisredprop = redprop(1, h, hamiltonian, redCols).asDiagonal();
  Eigen::MatrixXcd thisqprop = (((Bmat * thisredprop) * VT.adjoint()) * U_scaled.adjoint());
  std::cout << "Computed thisqprop!\n";
  return thisqprop;
}

int main(int argc, char** argv)
{
  // omp_set_max_active_levels(1); 
  int num_threads = 56; // omp_get_max_threads();
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
  Eigen::Map<const Eigen::MatrixXd> Bmat(Bten.data(), drc2, drcCI2);

  Eigen::VectorXi diagcols(drcCI);
  for (int k=0; k<drcCI; ++k) 
    diagcols(k) = k * drcCI + k; 

  // reduce the number of columns of Bmat
  Eigen::VectorXd colnorms(drcCI2);
  std::vector<int> goodCols;
  goodCols.reserve(drcCI2);
  int diagcheck = 0;
  Eigen::VectorXcd bvec = Eigen::VectorXd::Zero(drc2);
  for (int j=0; j<drcCI2; ++j)
  {
    // skip column whose index corresponds to diagonal entry of full density
    if (j==diagcols(diagcheck))
    {
      bvec += Bmat.col(j)/drcCI;
      diagcheck++;
      continue;
    }
    // retain columns only if their norm exceeds a threshold
    colnorms(j) = Bmat.col(j).norm();
    if (colnorms[j] > 1e-14)
      goodCols.push_back(j);
  }
  if (verbose)
    std::cout << "number of good cols = " << goodCols.size() << "\n";

  // create new matrix with selected columns
  Eigen::MatrixXd BmatR = Bmat(Eigen::placeholders::all, goodCols);
  if (verbose)
  {
    std::cout << "BmatR # of rows = " << BmatR.rows() << "\n";
    std::cout << "BmatR # of cols = " << BmatR.cols() << "\n";
  }
  Eigen::MatrixXcd thisqprop = qprop(delay, dt, ham, goodCols, BmatR, tol);
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
  #pragma omp parallel for 
  for (int k=0; k<=nsteps; ++k)
  {
    // outer product
    Eigen::MatrixXcd op = coeffs.col(k) * coeffs.col(k).adjoint();
    // reduction
    true1rdms.col(k) = Bmat * op.transpose().reshaped();
  }

  // set up predicted 1RDMs
  Eigen::MatrixXcd pred1rdms(drc2, nsteps+1);
  if (verbose)
    std::cout << "About to propagate 1RDMs for " << nsteps << " steps\n";
  
  for (int k=0; k<=delay; ++k)
    pred1rdms.col(k) = true1rdms.col(k);

  // (Assuming d_pred1rdms, d_bvec, d_temp2, d_thisqprop are already allocated and populated on the GPU)
  int temp2_size = drc2 * (delay + 1);
  int threadsPerBlock = 256;
  int blocksPerGrid = (temp2_size + threadsPerBlock - 1) / threadsPerBlock;

  cuDoubleComplex *d_thisqprop, *d_pred1rdms, *d_bvec, *d_temp2;
  
  size_t bytes = thisqprop.size() * sizeof(cuDoubleComplex);
  cudaMallocManaged(&d_thisqprop, bytes);
  cudaMemcpy(d_thisqprop, thisqprop.data(), bytes, cudaMemcpyHostToDevice);

  bytes = pred1rdms.size() * sizeof(cuDoubleComplex);
  cudaMallocManaged(&d_pred1rdms, bytes);
  cudaMemcpy(d_pred1rdms, pred1rdms.data(), bytes, cudaMemcpyHostToDevice);

  bytes = bvec.size() * sizeof(cuDoubleComplex);
  cudaMallocManaged(&d_bvec, bytes);
  cudaMemcpy(d_bvec, bvec.data(), bytes, cudaMemcpyHostToDevice);

  bytes = temp2_size * sizeof(cuDoubleComplex);
  cudaMallocManaged(&d_temp2, bytes);
  
  cublasHandle_t cublasH;
  cublasCreate(&cublasH);
  cuDoubleComplex alpha = {1.0, 0.0};
  cuDoubleComplex beta_zero = {0.0, 0.0}; // Beta is 0.0 so cuBLAS completely OVERWRITES the target

  int bvec_blocks = (drc2 + threadsPerBlock - 1) / threadsPerBlock;

  for (int k = delay; k < nsteps; ++k)
  {
    prep_temp2_kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_pred1rdms, d_bvec, d_temp2, k, delay, drc2
    );

    // Check for synchronous launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      printf("CUDA Error in prep_temp2_kernel: %s\n", cudaGetErrorString(err));
      exit(-1);
    }

    cuDoubleComplex* d_target_col = d_pred1rdms + ((k + 1) * drc2);

    // 2. cuBLAS computes: target = thisqprop * temp2 (overwriting whatever was there)
    cublasZgemv(cublasH, CUBLAS_OP_N,
                drc2, temp2_size,
                &alpha, 
                d_thisqprop, drc2, 
                d_temp2, 1,        
                &beta_zero,         // <--- STRICT ZERO
                d_target_col, 1);  

    // 3. Add bvec explicitly via kernel (target = target + bvec)
    add_bvec_kernel<<<bvec_blocks, threadsPerBlock>>>(d_target_col, d_bvec, drc2);
  }
  cudaDeviceSynchronize();

  cudaMemcpy(
    pred1rdms.data(), 
    d_pred1rdms, 
    pred1rdms.size() * sizeof(cuDoubleComplex), 
    cudaMemcpyDeviceToHost
  );
  
  /*
  for (int k=delay; k<nsteps; ++k)
  {
    Eigen::MatrixXcd temp = pred1rdms.block(0, k-delay, drc2, delay+1).rowwise().reverse();
    Eigen::MatrixXcd temp2 = (temp.colwise() - bvec).reshaped();
    pred1rdms.col(k+1) = bvec + thisqprop * temp2;
  }
  */
  
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


