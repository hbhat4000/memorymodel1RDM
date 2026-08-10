// memory model
// field-free case

#define EIGEN_USE_MKL_ALL

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
#include <queue>
#include <set>
#include <sstream>



#define CHECK_CUDA(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

using namespace std::complex_literals;

__global__ void scale_U_kernel(cuDoubleComplex* U, const double* S, size_t m, size_t min_mn, double tol) 
{
    // Cast blockIdx.x to size_t BEFORE multiplication to prevent 32-bit wrap-around
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t total_elements = m * min_mn;
    
    if (tid < total_elements) 
    {
        size_t c = tid / m; // Column index
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

        /*
        // Assuming max_allocated_cols is the true number of columns allocated for d_pred1rdms
        int max_allocated_cols = 20001; 
        int max_valid_idx = max_allocated_cols * drc2;
        
        if (source_idx < 0 || source_idx >= max_valid_idx) {
            printf("CRASH PREVENTED: k=%d, c=%d, source_idx=%d is out of bounds!\n", k, c, source_idx);
            return; // Exit thread early to prevent the actual segfault
        }
        */
        
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

Eigen::MatrixXcd compute_qprop_gpu(const Eigen::MatrixXcd& input, const Eigen::MatrixXcd& M1, double tol) 
{
    int m = input.rows();
    int n = input.cols();
    int lda = m;
    int min_mn = std::min(m, n);

    // 1. Initialize Handles
    cusolverDnHandle_t cusolverH = nullptr;
    cusolverDnCreate(&cusolverH);
    cusolverDnParams_t params;
    cusolverDnCreateParams(&params);
    
    cublasHandle_t cublasH = nullptr;
    cublasCreate(&cublasH);

    // 2. Allocate SVD Device Memory
    size_t m_sz = m;
    size_t n_sz = n;
    size_t min_mn_sz = min_mn;

    cuDoubleComplex *d_A, *d_U, *d_VT;
    double *d_S;
    int *devInfo;

    CHECK_CUDA(cudaMalloc(&d_A, m_sz * n_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_U, m_sz * min_mn_sz * sizeof(cuDoubleComplex)));   
    CHECK_CUDA(cudaMalloc(&d_VT, min_mn_sz * n_sz * sizeof(cuDoubleComplex)));  
    CHECK_CUDA(cudaMalloc(&d_S, min_mn_sz * sizeof(double)));                
    CHECK_CUDA(cudaMalloc(&devInfo, sizeof(int)));

    CHECK_CUDA(cudaMemcpy(d_A, input.data(), m_sz * n_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // 3. 64-bit SVD Workspace
    size_t workspaceInBytesOnDevice = 0;
    size_t workspaceInBytesOnHost = 0;

    cusolverDnXgesvd_bufferSize(
        cusolverH, params, 'S', 'S', 
        static_cast<int64_t>(m), static_cast<int64_t>(n), 
        CUDA_C_64F, d_A, static_cast<int64_t>(lda), 
        CUDA_R_64F, d_S, 
        CUDA_C_64F, d_U, static_cast<int64_t>(lda), 
        CUDA_C_64F, d_VT, static_cast<int64_t>(min_mn), 
        CUDA_C_64F, 
        &workspaceInBytesOnDevice, 
        &workspaceInBytesOnHost
    );

    void *d_work = nullptr;
    if (workspaceInBytesOnDevice > 0) CHECK_CUDA(cudaMalloc(&d_work, workspaceInBytesOnDevice));
    void *h_work = nullptr;
    if (workspaceInBytesOnHost > 0) h_work = malloc(workspaceInBytesOnHost);

    // 4. Run SVD
    cusolverDnXgesvd(
        cusolverH, params, 'S', 'S', 
        static_cast<int64_t>(m), static_cast<int64_t>(n), 
        CUDA_C_64F, d_A, static_cast<int64_t>(lda), 
        CUDA_R_64F, d_S, 
        CUDA_C_64F, d_U, static_cast<int64_t>(lda), 
        CUDA_C_64F, d_VT, static_cast<int64_t>(min_mn), 
        CUDA_C_64F, 
        d_work, workspaceInBytesOnDevice, h_work, workspaceInBytesOnHost, devInfo
    );
    CHECK_CUDA(cudaDeviceSynchronize());
    std::cout << "Done with SVD on GPU!\n";

    // 5. Scale U
    size_t threads = 256;
    size_t blocks = (m_sz * min_mn_sz + threads - 1) / threads;
    
    // Pass m_sz and min_mn_sz directly, as they are already size_t
    scale_U_kernel<<<blocks, threads>>>(d_U, d_S, m_sz, min_mn_sz, tol);
    
    CHECK_CUDA(cudaDeviceSynchronize());
    std::cout << "Columns of U scaled on GPU!\n";
    
    // =========================================================================
    // NEW MULTIPLICATION CHAIN (Avoids CPU download)
    // =========================================================================
    
    int m1_rows = M1.rows();
    int m1_cols = M1.cols();
    
    cuDoubleComplex *d_M1, *d_M2, *d_thisqprop;
    CHECK_CUDA(cudaMalloc(&d_M1, m1_rows * m1_cols * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_M2, m1_rows * min_mn_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_thisqprop, m1_rows * m_sz * sizeof(cuDoubleComplex)));

    // Upload tiny M1 to GPU
    CHECK_CUDA(cudaMemcpy(d_M1, M1.data(), m1_rows * m1_cols * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    
    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta = {0.0, 0.0};
    
    // Step A: M2 = M1 * d_VT^H 
    // M1 is (100 x 6530), d_VT^H is (6530 x 6530)
    cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_C, 
                m1_rows, min_mn_sz, m1_cols, 
                &alpha, 
                d_M1, m1_rows, 
                d_VT, min_mn_sz, 
                &beta, 
                d_M2, m1_rows);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step B: thisqprop = M2 * d_U^H
    // M2 is (100 x 6530), d_U^H is (6530 x 300100)
    cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_C, 
                m1_rows, m_sz, min_mn_sz, 
                &alpha, 
                d_M2, m1_rows, 
                d_U, m_sz, 
                &beta, 
                d_thisqprop, m1_rows);
    CHECK_CUDA(cudaDeviceSynchronize());
    std::cout << "Computed thisqprop on GPU!\n";
    
    // Download ONLY the final 480 MB result
    Eigen::MatrixXcd thisqprop_ret(m1_rows, m_sz);
    CHECK_CUDA(cudaMemcpy(thisqprop_ret.data(), d_thisqprop, m1_rows * m_sz * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    
    // 6. Clean up everything
    cudaFree(d_A); cudaFree(d_U); cudaFree(d_VT); cudaFree(d_S); 
    cudaFree(d_M1); cudaFree(d_M2); cudaFree(d_thisqprop);
    if (d_work) cudaFree(d_work); 
    cudaFree(devInfo);
    if (h_work) free(h_work);
    cusolverDnDestroyParams(params);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);
    
    return thisqprop_ret;
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
  
  // Compute M1 on the CPU (it's extremely small and fast)
  Eigen::MatrixXcd thisredprop = redprop(1, h, hamiltonian, redCols).asDiagonal();
  Eigen::MatrixXcd M1 = Bmat * thisredprop;
  
  // Hand bigmat and M1 to the GPU. The GPU will do the SVD, scaling, and all GEMMs.
  Eigen::MatrixXcd thisqprop = compute_qprop_gpu(bigmat, M1, tol);
  
  return thisqprop;
}

std::set<int> parseInts(const std::string& s)
{
  std::set<int> result;
  std::stringstream ss(s);
  std::string token;

  while (std::getline(ss, token, ','))
    result.insert(std::stoi(token));

  return result;
}

int main(int argc, char** argv)
{
  omp_set_max_active_levels(1); 
  int num_threads = 72;
  std::cout << "num_threads = " << num_threads << "\n";
  omp_set_num_threads(num_threads);
  Eigen::setNbThreads(num_threads);
  
  cxxopts::Options options("memoryFF", "Field-free memory model for 1RDM propagation");

  options.add_options()
  ("h,help", "Print usage")
  ("states", "String of states to put in initial superposition", cxxopts::value<std::string>())
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
  std::string infile, outpath, states;
  bool verbose, savemae, saveqprop;

  if (result.count("states")==0)
  {
    std::cout << "Must specify states!\n";
    return 1;
  }
  else
    states = result["states"].as<std::string>();
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

  // parse string of states into a set of ints
  std::set<int> statesset = parseInts(states);

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
  std::vector<int> goodCols, spingoodCols;
  goodCols.reserve(drcCI2);
  int diagcheck = 0;
  Eigen::VectorXcd bvec = Eigen::VectorXd::Zero(drc2);
  for (int j=0; j<drcCI2; ++j)
  {
    // skip column whose index corresponds to diagonal entry of full density
    if (j==diagcols(diagcheck))
    {
      diagcheck++;
      continue;
    }
    // retain columns only if their norm exceeds a threshold
    // and if the column corresponds to a state in our initial superposition
    colnorms(j) = Bmat.col(j).norm();
    if (colnorms[j] > 1e-10)
      spingoodCols.push_back(j);
    
    int jr = j / drcCI;
    int jc = j % drcCI;
    bool membership = (statesset.contains(jr) || statesset.contains(jc));
    if ((colnorms[j] > 1e-10) && membership)
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

  // added on Tues August 4th --> more physical initial condition
  // let us first make a binary matrix of size drcCI x drcCI
  // this matrix corresponds to the good columns in spingoodCols!
  Eigen::MatrixXi binmat = Eigen::MatrixXi::Zero(drcCI, drcCI);
  for (int ij=0; ij<spingoodCols.size(); ++ij)
  {
    int trow = spingoodCols[ij] / drcCI;
    int tcol = spingoodCols[ij] % drcCI;
    binmat(trow, tcol) = 1;
  }

  std::vector<bool> visited(drcCI, false);
  
  // This will hold the new order of our rows/columns
  std::vector<int> new_order; 
  
  // Optional: Keep track of the distinct families for analysis
  std::vector<std::vector<int>> spin_families; 

  // ---------------------------------------------------------
  // 2. Graph Traversal (Breadth-First Search)
  // ---------------------------------------------------------
  for (int i = 0; i < drcCI; ++i) {
      if (!visited[i]) {
          std::queue<int> q;
          std::vector<int> current_family;
          
          q.push(i);
          visited[i] = true;

          while (!q.empty()) {
              int curr = q.front();
              q.pop();
              
              // Add this state to our current block
              current_family.push_back(curr);
              new_order.push_back(curr); 

              // Look for all states connected to 'curr'
              for (int j = 0; j < drcCI; ++j) {
                  if (binmat(curr, j) != 0 && !visited[j]) {
                      visited[j] = true;
                      q.push(j);
                  }
              }
          }
          spin_families.push_back(current_family);
      }
  }

  // ---------------------------------------------------------
  // 3. Reorder the Matrix into Block-Diagonal Form
  // ---------------------------------------------------------
  Eigen::MatrixXi block_mat(drcCI, drcCI);
  for (int i = 0; i < drcCI; ++i) {
      for (int j = 0; j < drcCI; ++j) {
          // Pull from the original matrix using our grouped indices
          block_mat(i, j) = binmat(new_order[i], new_order[j]);
      }
  }

  // ---------------------------------------------------------
  // 4. Print the Results
  // ---------------------------------------------------------
  std::cout << "Original Scrambled Matrix:\n" << binmat << "\n\n";
  std::cout << "Reordered Block-Diagonal Matrix:\n" << block_mat << "\n\n";
  std::cout << "Identified Spin Families:\n";
  for (size_t f = 0; f < spin_families.size(); ++f) {
      std::cout << "Family " << f << " contains original states: ";
      // std::cout << "Size of family: " << spin_families[f].size() << "\n";
      for (int state : spin_families[f]) {
          std::cout << state << " ";
      }
      std::cout << "\n";
  }
  
  Eigen::MatrixXcd coeffs = Eigen::MatrixXcd::Zero(drcCI, nsteps+1);
  int smax = statesset.size();
  for (int x : statesset)
  {
    int goodstate = x; // spin_families[0][x];
    std::cout << "Superposition contains state " << goodstate << "\n";
    coeffs(goodstate, 0) = 1.0/std::sqrt(smax);
    bvec += Bmat.col(diagcols(goodstate))/smax;
  }
  
  // propagate TDCI coefficients
  Eigen::MatrixXcd fullprop = prop(1, dt, ham);
  for (int k=0; k<nsteps; ++k)
  {
    coeffs.col(k+1) = fullprop * coeffs.col(k);
  }
  
  if (verbose)
    std::cout << "Norm of solution at final time = " << coeffs.col(nsteps).norm() << "\n";

  /* UNCOMMENT THIS BLOCK IF YOU WANT TO SAVE THE TDCI TRAJECTORY TO DISK
  if (verbose)
  {
    std::filesystem::path p(infile);
    std::string stem = p.stem().string();
    std::string filename = stem + "_" + std::to_string(dt) + "_" + std::to_string(delay) + "_coeffs.txt";
    std::filesystem::path dir(outpath);
    std::filesystem::path outfile = dir / filename;
    std::ofstream out(outfile);
    for (int k=0; k<=nsteps; ++k)
    {
      for (int l=0; l<drcCI; ++l)
        out << coeffs(l, k).real() << "+" << coeffs(l, k).imag() << "j" << (l<(drcCI-1) ? "," : "");
      out << "\n";
    }
  }
  */

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
  CHECK_CUDA(cudaMalloc(&d_thisqprop, bytes));
  CHECK_CUDA(cudaMemcpy(d_thisqprop, thisqprop.data(), bytes, cudaMemcpyHostToDevice));

  bytes = pred1rdms.size() * sizeof(cuDoubleComplex);
  CHECK_CUDA(cudaMalloc(&d_pred1rdms, bytes));
  CHECK_CUDA(cudaMemcpy(d_pred1rdms, pred1rdms.data(), bytes, cudaMemcpyHostToDevice));

  bytes = bvec.size() * sizeof(cuDoubleComplex);
  CHECK_CUDA(cudaMalloc(&d_bvec, bytes));
  CHECK_CUDA(cudaMemcpy(d_bvec, bvec.data(), bytes, cudaMemcpyHostToDevice));

  bytes = temp2_size * sizeof(cuDoubleComplex);
  CHECK_CUDA(cudaMalloc(&d_temp2, bytes));
  
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
  CHECK_CUDA(cudaDeviceSynchronize());

  CHECK_CUDA(cudaMemcpy(
    pred1rdms.data(), 
    d_pred1rdms, 
    pred1rdms.size() * sizeof(cuDoubleComplex), 
    cudaMemcpyDeviceToHost
  ));
  
  /*
  for (int k=delay; k<nsteps; ++k)
  {
    Eigen::MatrixXcd temp = pred1rdms.block(0, k-delay, drc2, delay+1).rowwise().reverse();
    Eigen::MatrixXcd temp2 = (temp.colwise() - bvec).reshaped();
    pred1rdms.col(k+1) = bvec + thisqprop * temp2;
  }
  */

  double mae = (true1rdms(Eigen::placeholders::all,Eigen::seq(delay+1,Eigen::placeholders::last)) - pred1rdms(Eigen::placeholders::all,Eigen::seq(delay+1,Eigen::placeholders::last))).array().abs().mean();
  if (verbose)
    std::cout << "Mean Absolute Error: " << mae << "\n";

  std::ostringstream oss;
  oss << std::setprecision(17) << mae;
  std::string maestr = oss.str();
  if (savemae)
  {
    std::filesystem::path p(infile);
    std::string stem = p.stem().string();
    std::string filename = stem + "_" + std::to_string(dt) + "_" + std::to_string(delay) + "_" + states + ".txt";
    std::filesystem::path dir(outpath);
    std::filesystem::path outfile = dir / filename;
    std::ofstream out(outfile);
    out << maestr << "\n";
  }
  else
    std::cout << maestr << "\n";

  return 0;
}


