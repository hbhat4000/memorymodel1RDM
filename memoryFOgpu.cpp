// memory model
// field-on case
// redo the code to cache products
// also, no more incremental SVD

// #define EIGEN_USE_BLAS 
// #define EIGEN_USE_LAPACKE

#define EIGEN_USE_MKL_ALL

#include <mkl.h>
#include <omp.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cuComplex.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cxxopts.hpp>
#include "cnpy.h"

#include <utility>
#include <filesystem>
#include <string>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <chrono>

#define CHECK_CUDA(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

using namespace std::complex_literals;
using MatrixXcdRowMajor = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

__global__ void extract_R_kernel(cuDoubleComplex* d_R, const cuDoubleComplex* d_A, int m, int n) 
{
    int r = blockIdx.y * blockDim.y + threadIdx.y;
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (r < n && c < n) {
        if (r <= c) { // Upper triangular part
            d_R[c * n + r] = d_A[c * m + r];
        } else {      // Zero out the lower part
            d_R[c * n + r] = {0.0, 0.0};
        }
    }
}

__global__ void extract_reverse_qhist_kernel(
    cuDoubleComplex* d_qhist, 
    const cuDoubleComplex* d_pred1rdms,
    int J, int delay, int drc2, int lda)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = drc2 * (delay + 1);
    
    if (tid < total_elements) {
        int r = tid % drc2;             // Row index
        int col_offset = tid / drc2;    // Column offset (0 to delay)
        
        // Reverse the columns: offset 0 gets column J, offset delay gets J-delay
        int source_col = J - col_offset;
        
        d_qhist[tid] = d_pred1rdms[source_col * lda + r];
    }
}

// Kernel to handle the Eigen RowMajor -> ColMajor block assignment
// Equivalent to: bigmat.block(0, 0, drc2, N2) = BmatR;
__global__ void copy_bmatr_to_bigmat_kernel(
    cuDoubleComplex* d_bigmat,
    const cuDoubleComplex* d_BmatR,
    int delay, int drc2, int N2)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_elements = drc2 * N2;
    
    if (idx < total_elements) {
        int r = idx / N2; // row in BmatR
        int c = idx % N2; // col in BmatR
        
        int bigmat_lda = (delay + 1) * drc2;
        d_bigmat[c * bigmat_lda + r] = d_BmatR[idx];
    }
}

// Kernel to handle BmatMultNP and the strided assignment into bigmat
__global__ void bmat_mult_np_kernel(
    cuDoubleComplex* d_bigmat, 
    const cuDoubleComplex* d_BmatR, 
    const cuDoubleComplex* d_pcc,
    int J, int j, int delay, int N, int drc2, int N2)
{
    // Each block handles one 'k' (one row of BmatR)
    int k = blockIdx.x; 
    if (k >= drc2) return;

    // Thread handles specific output elements in the N x N grid
    int tid = threadIdx.x;
    
    // BmatR is RowMajor (drc2 x N2)
    const cuDoubleComplex* Bk_ptr = d_BmatR + k * N2; 
    
    // pcc is uploaded as a flattened 1D array of ColMajor matrices
    const cuDoubleComplex* pcc_J_ptr = d_pcc + J * (delay * N * N);
    int F_lda = delay * N;
    int F_row_offset = (j - 1) * N;

    for (int i = tid; i < N2; i += blockDim.x) {
        int u = i / N; // row index of R_k
        int v = i % N; // col index of R_k
        
        cuDoubleComplex sum = {0.0, 0.0};
        
        for (int x = 0; x < N; ++x) {
            // F(u, x) -> col x, row u + F_row_offset
            cuDoubleComplex F_ux = pcc_J_ptr[x * F_lda + (u + F_row_offset)];
            cuDoubleComplex L_ux = {F_ux.x, -F_ux.y}; // Conjugate
            
            for (int y = 0; y < N; ++y) {
                // Bk is RowMajor
                cuDoubleComplex B_xy = Bk_ptr[x * N + y];
                
                // F^T(y, v) = F(v, y) -> col y, row v + F_row_offset
                cuDoubleComplex F_vy = pcc_J_ptr[y * F_lda + (v + F_row_offset)];
                
                // Complex multiplication: B_xy * F_vy
                double r_temp = B_xy.x * F_vy.x - B_xy.y * F_vy.y;
                double i_temp = B_xy.x * F_vy.y + B_xy.y * F_vy.x;
                
                // Multiply by L_ux and accumulate
                sum.x += L_ux.x * r_temp - L_ux.y * i_temp;
                sum.y += L_ux.x * i_temp + L_ux.y * r_temp;
            }
        }
        
        // Write out to d_bigmat (ColMajor)
        int global_row = j * drc2 + k;
        int global_col = i;
        int bigmat_lda = (delay + 1) * drc2;
        
        d_bigmat[global_col * bigmat_lda + global_row] = sum;
    }
}

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

__global__ void scale_temp_kernel(cuDoubleComplex* temp, const double* S, size_t min_mn, double tol) 
{
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < min_mn) 
    {
        double s_val = S[tid];
        // Hard thresholding pseudoinverse
        double sinv_val = (s_val > tol) ? 1.0 / s_val : 0.0;
        
        // Scale the vector element in-place
        cuDoubleComplex t_val = temp[tid];
        t_val.x *= sinv_val;
        t_val.y *= sinv_val;
        temp[tid] = t_val;
    }
}

void bigmatFromCache_gpu(
    cuDoubleComplex* d_bigmat,
    const cuDoubleComplex* d_BmatR,
    const cuDoubleComplex* d_pcc,
    int J, int delay, int drc2, int N, int N2)
{
    // 1. Copy BmatR into the top block
    int threads_copy = 256;
    int blocks_copy = (drc2 * N2 + threads_copy - 1) / threads_copy;
    copy_bmatr_to_bigmat_kernel<<<blocks_copy, threads_copy>>>(d_bigmat, d_BmatR, delay, drc2, N2);
    
    // 2. Compute and assign the BmatMultNP blocks
    int threads_mult = 256; // 256 threads per row k
    int blocks_mult = drc2; // One block per row k
    
    for (int j = 1; j <= delay; ++j) {
        bmat_mult_np_kernel<<<blocks_mult, threads_mult>>>(
            d_bigmat, d_BmatR, d_pcc, J, j, delay, N, drc2, N2
        );
    }
    CHECK_CUDA(cudaDeviceSynchronize());
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

Eigen::VectorXcd compute_qprop_step_gpu(
    const Eigen::MatrixXcd& bigmat, 
    const Eigen::VectorXcd& qhist_vec,
    const Eigen::MatrixXcd& F,
    const MatrixXcdRowMajor& BmatR,
    int N,
    int drc2,
    double tol) 
{
    int m = bigmat.rows();
    int n = bigmat.cols(); // Note: n == N * N (or N2)
    int lda = m;
    int min_mn = std::min(m, n);

    // 1. Initialize Handles
    cusolverDnHandle_t cusolverH = nullptr;
    cusolverDnCreate(&cusolverH);
    cusolverDnParams_t params;
    cusolverDnCreateParams(&params);
    
    cublasHandle_t cublasH = nullptr;
    cublasCreate(&cublasH);

    // 2. Allocate Device Memory
    size_t m_sz = m;
    size_t n_sz = n;
    size_t min_mn_sz = min_mn;

    cuDoubleComplex *d_A, *d_U, *d_VT;
    double *d_S;
    int *devInfo;
    
    cuDoubleComplex *d_qhist, *d_temp, *d_PreconVec;
    cuDoubleComplex *d_F, *d_TempMat, *d_InnerMat;
    cuDoubleComplex *d_BmatR, *d_Result;

    CHECK_CUDA(cudaMalloc(&d_A, m_sz * n_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_U, m_sz * min_mn_sz * sizeof(cuDoubleComplex)));   
    CHECK_CUDA(cudaMalloc(&d_VT, min_mn_sz * n_sz * sizeof(cuDoubleComplex)));  
    CHECK_CUDA(cudaMalloc(&d_S, min_mn_sz * sizeof(double)));                
    CHECK_CUDA(cudaMalloc(&devInfo, sizeof(int)));

    // Allocations for the multiplication chain
    CHECK_CUDA(cudaMalloc(&d_qhist, m_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_temp, min_mn_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_PreconVec, n_sz * sizeof(cuDoubleComplex)));
    
    CHECK_CUDA(cudaMalloc(&d_F, n_sz * sizeof(cuDoubleComplex))); // n_sz == N * N
    CHECK_CUDA(cudaMalloc(&d_TempMat, n_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_InnerMat, n_sz * sizeof(cuDoubleComplex)));
    
    CHECK_CUDA(cudaMalloc(&d_BmatR, drc2 * n_sz * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_Result, drc2 * sizeof(cuDoubleComplex)));

    // Upload Data
    CHECK_CUDA(cudaMemcpy(d_A, bigmat.data(), m_sz * n_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_qhist, qhist_vec.data(), m_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_F, F.data(), n_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_BmatR, BmatR.data(), drc2 * n_sz * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // 3. Run SVD Workspace & Execution (Matching your original working code)
    size_t workspaceInBytesOnDevice = 0;
    size_t workspaceInBytesOnHost = 0;
    cusolverDnXgesvd_bufferSize(
        cusolverH, params, 'S', 'S', 
        static_cast<int64_t>(m), static_cast<int64_t>(n), 
        CUDA_C_64F, d_A, static_cast<int64_t>(lda), 
        CUDA_R_64F, d_S, 
        CUDA_C_64F, d_U, static_cast<int64_t>(lda), 
        CUDA_C_64F, d_VT, static_cast<int64_t>(min_mn), 
        CUDA_C_64F, &workspaceInBytesOnDevice, &workspaceInBytesOnHost
    );

    void *d_work = nullptr;
    if (workspaceInBytesOnDevice > 0) CHECK_CUDA(cudaMalloc(&d_work, workspaceInBytesOnDevice));
    void *h_work = nullptr;
    if (workspaceInBytesOnHost > 0) h_work = malloc(workspaceInBytesOnHost);

    cusolverDnXgesvd(
        cusolverH, params, 'S', 'S', 
        static_cast<int64_t>(m), static_cast<int64_t>(n), 
        CUDA_C_64F, d_A, static_cast<int64_t>(lda), 
        CUDA_R_64F, d_S, 
        CUDA_C_64F, d_U, static_cast<int64_t>(lda), 
        CUDA_C_64F, d_VT, static_cast<int64_t>(min_mn), 
        CUDA_C_64F, d_work, workspaceInBytesOnDevice, h_work, workspaceInBytesOnHost, devInfo
    );
    CHECK_CUDA(cudaDeviceSynchronize());

    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta = {0.0, 0.0};

    // =========================================================================
    // MULTIPLICATION CHAIN 
    // =========================================================================

    // Step A: temp = U^H * qhist
    // d_U is (m x min_mn), CUBLAS_OP_C transposes it to (min_mn x m)
    cublasZgemv(cublasH, CUBLAS_OP_C, 
                m, min_mn, 
                &alpha, d_U, m, 
                d_qhist, 1, 
                &beta, d_temp, 1);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step B: Scale temp by inverse singular values
    size_t threads = 256;
    size_t blocks = (min_mn_sz + threads - 1) / threads;
    scale_temp_kernel<<<blocks, threads>>>(d_temp, d_S, min_mn_sz, tol);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step C: PreconVec = VT^H * temp
    // d_VT is (min_mn x n), CUBLAS_OP_C transposes it to (n x min_mn)
    cublasZgemv(cublasH, CUBLAS_OP_C, 
                min_mn, n, 
                &alpha, d_VT, min_mn, 
                d_temp, 1, 
                &beta, d_PreconVec, 1);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step D: TempMat = Precon * F^T
    // d_PreconVec acts as an N x N column-major matrix
    cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_T, 
                N, N, N, 
                &alpha, d_PreconVec, N, 
                d_F, N, 
                &beta, d_TempMat, N);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step E: InnerMat = F^H * TempMat
    cublasZgemm(cublasH, CUBLAS_OP_C, CUBLAS_OP_N, 
                N, N, N, 
                &alpha, d_F, N, 
                d_TempMat, N, 
                &beta, d_InnerMat, N);
    CHECK_CUDA(cudaDeviceSynchronize());

    // Step F: Result = BmatR * InnerMat_vec
    // BmatR is stored RowMajor in host memory. When copied to device, 
    // it perfectly maps to an (n x drc2) ColumnMajor matrix.
    // Transposing it with CUBLAS_OP_T gives us the exact (drc2 x n) multiplication we need.
    cublasZgemv(cublasH, CUBLAS_OP_T, 
                n, drc2, 
                &alpha, d_BmatR, n, 
                d_InnerMat, 1, 
                &beta, d_Result, 1);
    CHECK_CUDA(cudaDeviceSynchronize());

    // 4. Download and Clean Up
    Eigen::VectorXcd result_vec(drc2);
    CHECK_CUDA(cudaMemcpy(result_vec.data(), d_Result, drc2 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    cudaFree(d_A); cudaFree(d_U); cudaFree(d_VT); cudaFree(d_S); cudaFree(devInfo);
    cudaFree(d_qhist); cudaFree(d_temp); cudaFree(d_PreconVec);
    cudaFree(d_F); cudaFree(d_TempMat); cudaFree(d_InnerMat);
    cudaFree(d_BmatR); cudaFree(d_Result);
    if (d_work) cudaFree(d_work); 
    if (h_work) free(h_work);
    cusolverDnDestroyParams(params);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);
    
    return result_vec;
}

class memoryModel
{
  const double h, T, tol, freq, amp;
  const int ncyc, nsteps, delay, numthreads;
  const std::string inpath, outpath;
  const bool verbose;
  int drc, drc2, drcCI, drcCI2;
  int N, N2;                               // N = number of good states
  int offstep;                             // time step at which the field is switched off
  bool havecoeffs = false;                 // have we computed time-dependent coefficients yet?
  bool have1rdms = false;                  // have we computed 1-RDMS yet?
  bool filtered = false;                   // have we filtered out bad indices yet?
  bool builtpcc = false;                   // have we built the propagator chain cache yet?
  Eigen::VectorXd H0;                      // core Hamiltonian
  Eigen::MatrixXd dimatz;                  // dipole moment matrix in the z-direction
  Eigen::MatrixXd Bmat;                    // B tensor, in matrix form, reshaped to size (drc2 x drcCI2)
  MatrixXcdRowMajor BmatR;                 // reduced B tensor, in matrix form, reshaped to size (drc2 x N2)
  Eigen::MatrixXcd coeffs;                 // time-dependent coefficients
  std::vector<Eigen::MatrixXcd> props;     // all time-dependent propagators "exp(-i H_n dt)"
  std::vector<Eigen::MatrixXcd> fprops;    // all time-dependent propagators "exp(-i H_n dt)"
  Eigen::MatrixXcd true1rdms;              // ground truth 1-RDMs
  Eigen::MatrixXcd pred1rdms;              // predicted 1-RDMs (at each of the delay values in delayrange)s
  std::unordered_set<int> goodStates;      // non-trivial indices of coefficient vector "coeffs"
  std::vector<int> goodCols;               // columns in "drcCI**2" space to retain
  Eigen::MatrixXcd kronprop;               // one-step field-free Kronecker product propagator
  // propagator chain cache (pcc)
  // the index for the std::vector is J, a discrete time step
  // each complex matrix in the cache is of size (ell*N) x N
  std::vector<Eigen::MatrixXcd> pcc;
 
  public:
    //constructor
    memoryModel(double dt, double T, double freq, double amp, int ncyc, int delay, double svdtol, int numthreads, std::string infile, std::string outpath, bool verbose);

    int getdrc(void) { return drc; }
    int getdrc2(void) { return drc2; }
    int getdrcCI(void) { return drcCI; }
    int getN(void) { return N; }
    int getnsteps(void) { return nsteps; }
    Eigen::MatrixXcd getTrue1rdms(void) { return true1rdms; }
    Eigen::MatrixXcd getPred1rdms(void) { return pred1rdms; }

    // solve TDSE forward in time starting from a particular initial condition,
    // saving propagators as we go
    int tdseProp(const Eigen::VectorXcd& ic);

    // compute and store all exact 1RDMS
    int exact1RDMS(void);

    // figure out which indices (among the drcCI**2 linear indices) we should retain/delete
    int filterIndices(void);

    MatrixXcdRowMajor BmatMult(const Eigen::MatrixXcd& leftmat, const Eigen::MatrixXcd& rightmat);
    MatrixXcdRowMajor BmatMultNP(const Eigen::MatrixXcd& leftmat, const Eigen::MatrixXcd& rightmat);

    // build propagator chain cache
    int buildPCC(void);

    // form bigmat at particular point in time for particular delay value
    int bigmatFromCache(const int J, Eigen::MatrixXcd& bigmat);
    int bigmatBuildLocal(const int J, Eigen::MatrixXcd& bigmat);

    // print bigmat to screen
    int bigmatPrint(Eigen::MatrixXcd bigmat);

    // propagate 1RDM with memory model for all delay values at once
    int qpropALLV2(void);
    int saveResults(void);
};

memoryModel::memoryModel(double dt, double T, double freq, double amp, int ncyc, int delay, double svdtol, int numthreads, std::string infile, std::string outpath, bool verbose)
   : h(dt), T(T), freq(freq), amp(amp), ncyc(ncyc), 
     delay(delay),
     tol(svdtol), numthreads(numthreads), 
     inpath(std::move(infile)), outpath(std::move(outpath)), 
     nsteps(static_cast<int>(std::ceil(T/h))), verbose(verbose)
{
  offstep = static_cast<int>(std::ceil(ncyc / (dt * freq)));
  if (verbose)
    std::cout << "Field will be on for " << offstep << " time steps\n";
  // Load the entire .npz file into a map-like structure
  cnpy::npz_t my_npz = cnpy::npz_load(inpath);

  // Load ham
  cnpy::NpyArray arr = my_npz["ham"];
  double* hamdata = arr.data<double>();
  size_t length = arr.shape[0];
  Eigen::Map<Eigen::VectorXd> H0_map(hamdata, length);
  H0 = H0_map;
  drcCI = (int) length;
  drcCI2 = drcCI*drcCI;
  if (verbose)
    std::cout << "drcCI = " << drcCI << "\n";

  // Load dipole moment matrix (in z direction)
  arr = my_npz["CIdimatz"];
  double* cidimatz = arr.data<double>();
  Eigen::Map<Eigen::VectorXd> CIdimatz(cidimatz, drcCI2);
  Eigen::Map<const Eigen::MatrixXd> dimatz_map(CIdimatz.data(), drcCI, drcCI);
  dimatz = dimatz_map;
    
  // Load Bten  
  arr = my_npz["Bten"];
  double* Btendata = arr.data<double>();
  length = arr.shape[0]*arr.shape[1]*arr.shape[2]*arr.shape[3];
  drc = (int) arr.shape[2];
  drc2 = drc*drc;
  if (verbose)
    std::cout << "drc = " << drc << "\n";
  Eigen::Map<Eigen::VectorXd> Bten(Btendata, length);
  Eigen::Map<const Eigen::MatrixXd> Bmat_map(Bten.data(), drc2, drcCI2);
  Bmat = Bmat_map;

  // initialize various objects
  coeffs.setZero(drcCI, nsteps+1);
  true1rdms.setZero(drc2, nsteps+1);
  pred1rdms.setZero(drc2, nsteps+1);
}

int memoryModel::tdseProp(const Eigen::VectorXcd& ic)
{
  double field = 0.0;
  Eigen::MatrixXcd H;
  Eigen::MatrixXcd prop;
  Eigen::VectorXcd D, Dexp;
  std::complex<double> scalarfac = -1.0i * h;

  // copy initial condition into coeffs
  for (int j=0; j<drcCI; ++j)
    coeffs(j, 0) = ic(j);
    
  // initialize propagators to be the identity
  props.resize(nsteps, Eigen::MatrixXcd::Identity(drcCI, drcCI));

  for (int k=0; k<nsteps; ++k)
  {
    if (k < offstep)
    {
      H = H0.asDiagonal();
      field = amp * std::sin(2 * EIGEN_PI * freq * k * h);
      H += field * dimatz;
      if (! H.isApprox(H.adjoint(), 1e-12)) std::cout << "H is not Hermitian at step " << k << "\n";
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);
      D = solver.eigenvalues();
      D = scalarfac * D;
      Dexp = D.array().exp();
      prop = solver.eigenvectors() * Dexp.asDiagonal() * solver.eigenvectors().adjoint();
      props[k] = prop;
    }
    else
    {
      D = scalarfac * H0;
      Dexp = D.array().exp();
      props[k] = Dexp.asDiagonal();
    }
    coeffs.col(k+1) = props[k] * coeffs.col(k);
  }
  havecoeffs = true;
  return 0;
}

int memoryModel::exact1RDMS(void)
{
  if (!havecoeffs)
  {
    std::cout << "TDCI coefficients have not been computed yet!\n";
    return 1;
  }
  // compute ground truth 1RDMs
  #pragma omp parallel for
  for (int k=0; k<=nsteps; ++k)
  {
    // outer product
    Eigen::MatrixXcd op = coeffs.col(k) * coeffs.col(k).adjoint();
    // reduction
    true1rdms.col(k).noalias() = Bmat * op.transpose().reshaped();
  }
  have1rdms = true;
  return 0;
}

int memoryModel::filterIndices(void)
{
  if (!havecoeffs)
  {
    std::cout << "TDCI coefficients have not been computed yet!\n";
    return 1;
  }
  const double thresh = 1e-10;
  double rownorm;
  for (int j=0; j<drcCI; ++j)
  {
    rownorm = coeffs.row(j).norm();
    if (rownorm >= thresh) goodStates.insert(j);
  }

  // now figure out which entries in drcCI**2 space we should retain
  int cnt = 0;
  for (int row=0; row<drcCI; ++row)
  {
    for (int col=0; col<drcCI; ++col)
    {
      if (goodStates.find(row) != goodStates.end())
        if (goodStates.find(col) != goodStates.end())
          goodCols.push_back(cnt);
      cnt++;
    }    
  }

  N = goodStates.size();
  N2 = N*N;
  if (verbose)
    std::cout << "Retaining " << N2 << " or " << goodCols.size() << " entries\n";

  // reduce the B matrix
  BmatR = Bmat(Eigen::placeholders::all, goodCols);

  // filter all the propagators
  fprops.resize(nsteps, Eigen::MatrixXcd::Identity(N, N));
  std::vector<int> goodStatesVec(goodStates.begin(), goodStates.end());
  std::sort(goodStatesVec.begin(), goodStatesVec.end());

  #pragma omp parallel for
  for (int k=0; k<nsteps; ++k)
    fprops[k].noalias() = props[k](goodStatesVec, goodStatesVec);

  // we need this field-free Kronecker propagator below
  // we compute it here because we have goodStatesVec handy
  std::complex<double> coeff = -1.0i * h;
  Eigen::VectorXcd factor = coeff*H0.array();
  Eigen::VectorXcd diagvals = factor.array().exp();
  diagvals = diagvals(goodStatesVec);
  Eigen::MatrixXcd outer_product = diagvals * diagvals.adjoint();
  Eigen::VectorXcd kron_diag = outer_product.reshaped<Eigen::RowMajor>();
  kronprop = kron_diag.asDiagonal();
  
  filtered = true;
  return 0;
}

MatrixXcdRowMajor memoryModel::BmatMultNP(const Eigen::MatrixXcd& leftmat, const Eigen::MatrixXcd& rightmat)
{
  MatrixXcdRowMajor Result(drc2, N2);
  for (int k=0; k<drc2; ++k)
  {
    Eigen::Map<MatrixXcdRowMajor> BmatRk(BmatR.row(k).data(), N, N);
    Eigen::Map<MatrixXcdRowMajor> Rk(Result.row(k).data(), N, N);
    Rk.noalias() = (leftmat * BmatRk) * rightmat;
  }
  return Result;
}

MatrixXcdRowMajor memoryModel::BmatMult(const Eigen::MatrixXcd& leftmat, const Eigen::MatrixXcd& rightmat)
{
  MatrixXcdRowMajor Result(drc2, N2);
  #pragma omp parallel for schedule(dynamic)
  for (int k=0; k<drc2; ++k)
  {
    // Get the ID of the thread currently running this iteration
    Eigen::Map<MatrixXcdRowMajor> BmatRk(BmatR.row(k).data(), N, N);
    Eigen::Map<MatrixXcdRowMajor> Rk(Result.row(k).data(), N, N);
    Rk.noalias() = (leftmat * BmatRk) * rightmat;
  }
  return Result;
}

int memoryModel::buildPCC(void)
{
  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return 1;
  }
  // let the major axis be nsteps to avoid index madness later
  // essentially, we want to be able to refer to this things using the absolute time step index J
  // pcc[J] should be of size (ell*N) x N
  int upper = delay + offstep;
  int pccsize = (nsteps < upper) ? nsteps : upper;
  pcc.resize(pccsize + 1);
  for (int J=delay; J<=(delay+offstep); ++J)
  {
    if (J==nsteps) break;
    pcc[J].setZero(delay*N, N);
    pcc[J].block(0, 0, N, N) = fprops[J-1];
    if (J==delay) // construct a propagator chain that goes back delay steps
    {
      for (int j=2; j<=delay; ++j)
        pcc[J].block((j-1)*N, 0, N, N) = pcc[J].block((j-2)*N, 0, N, N) * fprops[J-j];
    }
    else 
    {
      #pragma omp parallel for schedule(dynamic)
      for (int j=2; j<=delay; ++j)
      {
        pcc[J].block((j-1)*N, 0, N, N).noalias() = fprops[J-1] * pcc[J-1].block((j-2)*N, 0, N, N);
      }
    }
  }
  builtpcc = true;
  if (verbose)
    std::cout << "Built propagator chain cache!\n";
  return 0;
}

int memoryModel::bigmatFromCache(const int J, Eigen::MatrixXcd& bigmat)
{
  Eigen::MatrixXcd thisblock;
  if (!builtpcc)
  {
    std::cout << "Propagator chain cache must be computed first!\n";
    return 1;
  }
  bigmat.block(0, 0, drc2, N2) = BmatR;
  for (int j=1; j<=delay; ++j)
  {
    thisblock = pcc[J].block((j-1)*N,0,N,N);
    bigmat.block(j*drc2, 0, drc2, N2).noalias() = BmatMultNP(thisblock.conjugate(), thisblock.transpose());
  }
  return 0;
}

// form bigmat at particular point in time for particular delay value
int memoryModel::bigmatBuildLocal(const int J, Eigen::MatrixXcd& bigmat)
{
  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return 1;
  }
  bigmat.block(0, 0, drc2, N2) = BmatR;
  // this will hold our propagator chain
  Eigen::MatrixXcd Amat(N, N);
  Amat.setIdentity( N, N );
  for (int j=1; j<=delay; ++j)
  {
    Amat = Amat * fprops[J - j];
    bigmat.block(j*drc2, 0, drc2, N2) = BmatMultNP(Amat.conjugate(), Amat.transpose());
  }
  return 0;
}

int memoryModel::bigmatPrint(Eigen::MatrixXcd bigmat)
{
  for (int j=0; j<bigmat.rows(); ++j)
  {
    for (int k=0; k<N2; ++k)
    {
      std::cout << std::setprecision(15) << bigmat(j,k).real() << "+" << bigmat(j,k).imag() << "j";
      if (k<(N2-1)) std::cout << ",";
    }
    std::cout << "\n";
  }
  return 0;
}

int memoryModel::qpropALLV2(void)
{
    if (!builtpcc) {
        std::cout << "Propagator cache must be computed first!\n";
        return 1;
    }

    // Matrix dimensions
    int m = (delay + 1) * drc2;
    int n = N2;

    // -------------------------------------------------------------------------
    // 1. INITIALIZE HANDLES & PRE-COMPUTE CPU MATRICES
    // -------------------------------------------------------------------------
    cusolverDnHandle_t cusolverH = nullptr;
    cusolverDnCreate(&cusolverH);
    cublasHandle_t cublasH = nullptr;
    cublasCreate(&cublasH);

    Eigen::MatrixXcd M1 = BmatR * kronprop;
    int m1_rows = M1.rows();
    int m1_cols = M1.cols(); // Note: m1_cols is exactly equal to n (N2)

    // Flatten std::vectors for contiguous GPU upload
    std::vector<cuDoubleComplex> pcc_flat(pcc.size() * delay * N * N);
    for (size_t J = 0; J < pcc.size(); ++J) {
        if (pcc[J].size() > 0)
            memcpy(&pcc_flat[J * delay * N * N], pcc[J].data(), delay * N * N * sizeof(cuDoubleComplex));
    }
    std::vector<cuDoubleComplex> fprops_flat(nsteps * N * N);
    for (int J = 0; J < nsteps; ++J) {
        memcpy(&fprops_flat[J * N * N], fprops[J].data(), N * N * sizeof(cuDoubleComplex));
    }

    // Copy in ground truth 1RDMs for the initial delay period
    for (int k = 0; k <= delay; ++k)
        pred1rdms.col(k) = true1rdms.col(k);

    // -------------------------------------------------------------------------
    // 2. ALLOCATE ALL DEVICE MEMORY (ONCE)
    // -------------------------------------------------------------------------
    cuDoubleComplex *d_pred1rdms, *d_BmatR, *d_pcc, *d_fprops;
    CHECK_CUDA(cudaMalloc(&d_pred1rdms, drc2 * (nsteps + 1) * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_BmatR, drc2 * N2 * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_pcc, pcc_flat.size() * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_fprops, fprops_flat.size() * sizeof(cuDoubleComplex)));

    // TSQR Matrices (Notice everything except bigmat is now tiny n x n)
    cuDoubleComplex *d_bigmat, *d_R, *d_UR, *d_VT, *d_tau;
    double *d_S, *d_rwork;
    int *devInfo;
    CHECK_CUDA(cudaMalloc(&d_bigmat, m * n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_R, n * n * sizeof(cuDoubleComplex)));        
    CHECK_CUDA(cudaMalloc(&d_UR, n * n * sizeof(cuDoubleComplex)));       
    CHECK_CUDA(cudaMalloc(&d_VT, n * n * sizeof(cuDoubleComplex)));  // V^H
    CHECK_CUDA(cudaMalloc(&d_tau, n * sizeof(cuDoubleComplex)));          
    CHECK_CUDA(cudaMalloc(&d_S, n * sizeof(double)));                     
    CHECK_CUDA(cudaMalloc(&d_rwork, 5 * n * sizeof(double)));             
    CHECK_CUDA(cudaMalloc(&devInfo, sizeof(int)));

    // Memory for Field-On Steps
    cuDoubleComplex *d_qhist, *d_temp, *d_temp2, *d_PreconVec, *d_TempMat, *d_InnerMat;
    CHECK_CUDA(cudaMalloc(&d_qhist, m * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_temp, n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_temp2, n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_PreconVec, n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_TempMat, n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_InnerMat, n * sizeof(cuDoubleComplex)));

    // Memory for Field-Off Steps
    cuDoubleComplex *d_M1, *d_M2, *d_M3, *d_rdmprop;
    CHECK_CUDA(cudaMalloc(&d_M1, m1_rows * n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_M2, m1_rows * n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_M3, m1_rows * n * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_rdmprop, m1_rows * m * sizeof(cuDoubleComplex)));

    // -------------------------------------------------------------------------
    // 3. UPLOAD DATA & SVD WORKSPACE SIZING
    // -------------------------------------------------------------------------
    CHECK_CUDA(cudaMemcpy(d_pred1rdms, pred1rdms.data(), drc2 * (nsteps + 1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_BmatR, BmatR.data(), drc2 * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_pcc, pcc_flat.data(), pcc_flat.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_fprops, fprops_flat.data(), fprops_flat.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_M1, M1.data(), m1_rows * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // Calculate maximum workspace required across all TSQR operations
    int lwork_geqrf = 0, lwork_ungqr = 0, lwork_svd = 0;
    cusolverDnZgeqrf_bufferSize(cusolverH, m, n, d_bigmat, m, &lwork_geqrf);
    cusolverDnZungqr_bufferSize(cusolverH, m, n, n, d_bigmat, m, d_tau, &lwork_ungqr);
    cusolverDnZgesvd_bufferSize(cusolverH, n, n, &lwork_svd);
    
    int lwork = std::max({lwork_geqrf, lwork_ungqr, lwork_svd});
    cuDoubleComplex *d_work = nullptr;
    CHECK_CUDA(cudaMalloc(&d_work, lwork * sizeof(cuDoubleComplex)));

    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta = {0.0, 0.0};

    // -------------------------------------------------------------------------
    // 4. THE TIME-STEPPING LOOP (DEVICE RESIDENT, TSQR)
    // -------------------------------------------------------------------------
    for (int J = delay; J < nsteps; ++J)
    {
        if (verbose) std::cout << "Propagating 1RDM at time step " << J << "\n";

        if ((J - delay) >= offstep) 
        {
            if ((J - delay) == offstep) 
            {
                // Field just turned off: Compute static rdmprop ONCE
                bigmatFromCache_gpu(d_bigmat, d_BmatR, d_pcc, J, delay, drc2, N, N2);
                
                // 1. QR Decomposition
                cusolverDnZgeqrf(cusolverH, m, n, d_bigmat, m, d_tau, d_work, lwork, devInfo);
                
                // 2. Extract R
                dim3 threads(16, 16);
                dim3 blocks((n + 15) / 16, (n + 15) / 16);
                extract_R_kernel<<<blocks, threads>>>(d_R, d_bigmat, m, n);
                
                // 3. Form Q explicitly in d_bigmat
                cusolverDnZungqr(cusolverH, m, n, n, d_bigmat, m, d_tau, d_work, lwork, devInfo);
                
                // 4. Tiny SVD on R
                cusolverDnZgesvd(cusolverH, 'S', 'S', n, n, d_R, n, d_S, d_UR, n, d_VT, n, d_work, lwork, d_rwork, devInfo);
                
                // Scale U_R
                size_t threads_scale = 256;
                size_t blocks_scale = (n * n + threads_scale - 1) / threads_scale;
                scale_U_kernel<<<blocks_scale, threads_scale>>>(d_UR, d_S, n, n, tol);
                
                // M2 = M1 * V (d_VT is V^H)
                cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_C, 
                            m1_rows, n, n, &alpha, 
                            d_M1, m1_rows, d_VT, n, &beta, d_M2, m1_rows);

                // M3 = M2 * U_R^H
                cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_C, 
                            m1_rows, n, n, &alpha, 
                            d_M2, m1_rows, d_UR, n, &beta, d_M3, m1_rows);
                
                // rdmprop = M3 * Q^H (d_bigmat now contains Q)
                cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_C, 
                            m1_rows, m, n, &alpha, 
                            d_M3, m1_rows, d_bigmat, m, &beta, d_rdmprop, m1_rows);
            }
            
            // Extract reversed history
            int threads_q = 256;
            int blocks_q = (m + threads_q - 1) / threads_q;
            extract_reverse_qhist_kernel<<<blocks_q, threads_q>>>(d_qhist, d_pred1rdms, J, delay, drc2, drc2);

            // pred1rdms.col(J+1) = rdmprop * qhist
            cublasZgemv(cublasH, CUBLAS_OP_N, drc2, m, 
                        &alpha, d_rdmprop, drc2, d_qhist, 1, 
                        &beta, d_pred1rdms + (J + 1) * drc2, 1);
        }
        else 
        {
            // Field is ON: Dynamic Memory Propagation
            bigmatFromCache_gpu(d_bigmat, d_BmatR, d_pcc, J, delay, drc2, N, N2);
            
            // 1. QR Decomposition
            cusolverDnZgeqrf(cusolverH, m, n, d_bigmat, m, d_tau, d_work, lwork, devInfo);
            
            // 2. Extract R
            dim3 threads_ext(16, 16);
            dim3 blocks_ext((n + 15) / 16, (n + 15) / 16);
            extract_R_kernel<<<blocks_ext, threads_ext>>>(d_R, d_bigmat, m, n);
            
            // 3. Form Q explicitly in d_bigmat
            cusolverDnZungqr(cusolverH, m, n, n, d_bigmat, m, d_tau, d_work, lwork, devInfo);
            
            // 4. Tiny SVD on R
            cusolverDnZgesvd(cusolverH, 'S', 'S', n, n, d_R, n, d_S, d_UR, n, d_VT, n, d_work, lwork, d_rwork, devInfo);
            
            // Extract history
            int threads_q = 256;
            int blocks_q = (m + threads_q - 1) / threads_q;
            extract_reverse_qhist_kernel<<<blocks_q, threads_q>>>(d_qhist, d_pred1rdms, J, delay, drc2, drc2);
            
            // b = Q^H * qhist (d_bigmat is Q)
            cublasZgemv(cublasH, CUBLAS_OP_C, m, n, 
                        &alpha, d_bigmat, m, d_qhist, 1, 
                        &beta, d_temp, 1);

            // temp2 = U_R^H * b
            cublasZgemv(cublasH, CUBLAS_OP_C, n, n, 
                        &alpha, d_UR, n, d_temp, 1, 
                        &beta, d_temp2, 1);
            
            // Scale temp2
            int threads_scale = 256;
            int blocks_scale = (n + threads_scale - 1) / threads_scale;
            scale_temp_kernel<<<blocks_scale, threads_scale>>>(d_temp2, d_S, n, tol);
            
            // PreconVec = V * temp2 (d_VT is V^H, OP_C makes it V)
            cublasZgemv(cublasH, CUBLAS_OP_C, n, n, 
                        &alpha, d_VT, n, d_temp2, 1, 
                        &beta, d_PreconVec, 1);
            
            // Pointer to fprops[J] on device
            cuDoubleComplex* d_F_J = d_fprops + J * N * N;
            
            // TempMat = Precon * F^T
            cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_T, 
                        N, N, N, &alpha, d_PreconVec, N, d_F_J, N, &beta, d_TempMat, N);
            
            // InnerMat = F^H * TempMat
            cublasZgemm(cublasH, CUBLAS_OP_C, CUBLAS_OP_N, 
                        N, N, N, &alpha, d_F_J, N, d_TempMat, N, &beta, d_InnerMat, N);
            
            // Result = BmatR * InnerMat_vec
            cublasZgemv(cublasH, CUBLAS_OP_T, n, drc2, 
                        &alpha, d_BmatR, n, d_InnerMat, 1, 
                        &beta, d_pred1rdms + (J + 1) * drc2, 1);
        }
    }

    // -------------------------------------------------------------------------
    // 5. DOWNLOAD RESULTS & CLEANUP
    // -------------------------------------------------------------------------
    CHECK_CUDA(cudaDeviceSynchronize());
    CHECK_CUDA(cudaMemcpy(pred1rdms.data(), d_pred1rdms, drc2 * (nsteps + 1) * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    cudaFree(d_pred1rdms); cudaFree(d_BmatR); cudaFree(d_pcc); cudaFree(d_fprops);
    cudaFree(d_bigmat); cudaFree(d_R); cudaFree(d_UR); cudaFree(d_VT); cudaFree(d_tau);
    cudaFree(d_S); cudaFree(d_rwork); cudaFree(devInfo);
    cudaFree(d_qhist); cudaFree(d_temp); cudaFree(d_temp2); cudaFree(d_PreconVec); 
    cudaFree(d_TempMat); cudaFree(d_InnerMat);
    cudaFree(d_M1); cudaFree(d_M2); cudaFree(d_M3); cudaFree(d_rdmprop); 
    if (d_work) cudaFree(d_work);
    
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);

    return 0;
}
  
int memoryModel::saveResults(void)
{
  // save to outfile
  std::filesystem::path p(inpath);
  std::string stem = p.stem().string();
  std::string filename = stem + "_" + std::to_string(h) + ".txt";
  // outpath comes from command-line argument
  std::filesystem::path dir(outpath);
  // this way we don't have to worry about trailing slashes
  // outfile is defined right here for the first time
  std::filesystem::path outfile = dir / filename;
  // out is the actual thing we write to (below)
  std::ofstream out(outfile);
  // save parameters
  out << "h = " << h << ", T = " << T << ", nsteps = " << nsteps << "\n";
  out << "freq = " << freq << ", amp = " << amp << ", ncyc = " << ncyc << "\n";
  Eigen::MatrixXcd residual = getTrue1rdms() - getPred1rdms();
  double maeTruth = (residual(Eigen::placeholders::all,Eigen::seq(delay+1,Eigen::placeholders::last))).array().abs().mean();
  std::ostringstream oss;
  oss << std::setprecision(17) << maeTruth;
  out << delay << ", " << oss.str() << "\n";
  return 0;
}

// 1. For REQUIRED parameters
template <typename T>
void getRequired(const cxxopts::ParseResult& res, const std::string& key, T& var, const std::string& errMsg)
{
    if (res.count(key) == 0)
    {
        std::cerr << errMsg << "\n";
        exit(1);
    }
    var = res[key].as<T>();
}

// 2. For OPTIONAL parameters with a default
template <typename T>
void getOptional(const cxxopts::ParseResult& res, const std::string& key, T& var, const T& defaultValue)
{
    if (res.count(key) == 0)
        var = defaultValue;
    else
        var = res[key].as<T>();
}

int main(int argc, char** argv)
{
  omp_set_max_active_levels(1); 
  int num_threads = omp_get_max_threads();
  std::cout << "num_threads = " << num_threads << "\n";
  omp_set_num_threads(num_threads);
  
  cxxopts::Options options("memoryFO", "Field-on memory model for 1RDM propagation");
  options.add_options()
  ("h,help", "Print usage")
  ("time", "Time-stepping parameters dt,T", cxxopts::value<std::vector<double>>())
  ("field", "Freq,amp,ncyc of applied electric field in z direction", cxxopts::value<std::vector<double>>())
  ("delay", "Integer value of delay/memory in time steps", cxxopts::value<int>())
  ("infile", "Input file path", cxxopts::value<std::string>())
  ("outpath", "Output file path (for MAEs)", cxxopts::value<std::string>())
  ("tol", "SVD tolerance", cxxopts::value<double>()->default_value("1e-6"))
  ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
  
  auto result = options.parse(argc, argv);
  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0;
  }
  double tol;
  std::vector<double> timeparams, fieldparams;
  int delay;
  std::string inpath, outpath;
  bool verbose;
  
  // Required parameters
  getRequired(result, "time",    timeparams,  "Must specify time-stepping parameters dt,T!");
  getRequired(result, "field",   fieldparams, "Must specify field parameters freq,amp,ncyc!");
  getRequired(result, "delay",   delay,       "Must specify delay!");
  getRequired(result, "infile",  inpath,      "Must specify input file!");
  getRequired(result, "outpath", outpath,     "Must specify output path!");
  
  // Optional parameter
  getOptional(result, "tol", tol, 1e-6);
  getOptional(result, "verbose", verbose, false);
  
  // pass the raw parameters into the constructor, in keeping with the other (time & field) params
  int ncyc = static_cast<int>(std::round(fieldparams[2]));
  memoryModel mm(timeparams[0], timeparams[1], fieldparams[0], fieldparams[1], ncyc, delay, tol, num_threads, inpath, outpath, verbose);
  Eigen::VectorXcd ic(mm.getdrcCI());
  ic.setZero();
  ic[0] = 1.0;
  mm.tdseProp(ic);
  mm.exact1RDMS();
  mm.filterIndices();
  mm.buildPCC();
  mm.qpropALLV2();
  mm.saveResults();
  return 0;
}

