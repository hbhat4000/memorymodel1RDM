// memory model
// field-on case
// redo the code to cache products
// also, no more incremental SVD

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

__global__ void init_pcc_j1_kernel(
    cuDoubleComplex* d_pcc,
    const cuDoubleComplex* d_fprops,
    int delay, int N, int pccsize, int delay_start)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    int num_J = pccsize - delay_start + 1;
    size_t total = (size_t)num_J * N * N;

    if (idx < total) {
        int J_offset = idx / (N * N);
        int J = delay_start + J_offset;
        int mat_idx = idx % (N * N);

        // d_pcc layout: [J] [delay] [N * N]
        // We write to block j=1 (which is index 0)
        size_t pcc_idx = (size_t)J * delay * N * N + mat_idx;
        size_t f_idx = (size_t)(J - 1) * N * N + mat_idx;

        d_pcc[pcc_idx] = d_fprops[f_idx];
    }
}

__global__ void compute_exact1rdms_kernel(
    cuDoubleComplex* d_true1rdms,
    const cuDoubleComplex* d_coeffs,
    const double* d_Bmat,
    int drcCI,
    int drc2,
    int nsteps)
{
    // Each block handles a single time step 'k'
    int k = blockIdx.x;
    if (k > nsteps) return;

    // Dynamically allocate shared memory for the current step's coefficient vector
    extern __shared__ cuDoubleComplex s_c[];

    // 1. Collaboratively load coeffs.col(k) into Shared Memory
    for (int i = threadIdx.x; i < drcCI; i += blockDim.x) {
        s_c[i] = d_coeffs[k * drcCI + i];
    }
    __syncthreads();

    // 2. Each thread computes one (or more) rows 'r' of the resulting true1rdms column
    for (int r = threadIdx.x; r < drc2; r += blockDim.x) {
        cuDoubleComplex sum = {0.0, 0.0};

        int j = 0; // Flat index representing the reshaped transposed outer product

        for (int col = 0; col < drcCI; ++col) {
            cuDoubleComplex c_col = s_c[col];

            for (int row = 0; row < drcCI; ++row) {
                cuDoubleComplex c_row = s_c[row];

                // Mathematically evaluate: op(row, col) of transpose = c_col * conj(c_row)
                // (x + iy) * (u - iv) = (xu + yv) + i(yu - xv)
                double P_real = c_col.x * c_row.x + c_col.y * c_row.y;
                double P_imag = c_col.y * c_row.x - c_col.x * c_row.y;

                // Bmat is Column-Major: drc2 rows, drcCI2 columns
                double B_val = d_Bmat[j * drc2 + r];

                sum.x += B_val * P_real;
                sum.y += B_val * P_imag;

                j++;
            }
        }

        // Write the reduced sum directly to the final 1RDM matrix
        d_true1rdms[k * drc2 + r] = sum;
    }
}

// Kernel 1: Construct the Hamiltonian matrices for all field-on steps
__global__ void build_H_batched_kernel(
    cuDoubleComplex* d_H, const double* d_H0, const double* d_dimatz,
    double amp, double freq, double h, double pi_val, int m, int batchSize)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total_elements = (size_t)batchSize * m * m;
    
    if (idx < total_elements) {
        int k = idx / (m * m);
        int rem = idx % (m * m);
        int col = rem / m;
        int row = rem % m;

        double t = k * h;
        double field = amp * sin(2.0 * pi_val * freq * t);

        // d_dimatz is read from memory assuming Column-Major (standard Eigen)
        double val_real = field * d_dimatz[col * m + row]; 
        double val_imag = 0.0;

        if (row == col) {
            val_real += d_H0[row];
        }
        d_H[idx] = {val_real, val_imag};
    }
}

// Kernel 2: Reconstruct the unitary propagators U = V * exp(-ihD) * V^H
__global__ void build_props_batched_kernel(
    cuDoubleComplex* d_props, const cuDoubleComplex* d_V, const double* d_W,
    double h, int m, int batchSize)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)batchSize * m * m;
    
    if (idx < total) {
        int k = idx / (m * m);
        int rem = idx % (m * m);
        int col = rem / m;
        int row = rem % m;

        cuDoubleComplex sum = {0.0, 0.0};
        
        for (int j = 0; j < m; ++j) {
            cuDoubleComplex v1 = d_V[(size_t)k * m * m + j * m + row]; // V(row, j)
            cuDoubleComplex v2 = d_V[(size_t)k * m * m + j * m + col]; // V(col, j)

            double phase = -h * d_W[k * m + j];
            double cos_p = cos(phase);
            double sin_p = sin(phase);

            // term = v1 * exp(-i * h * W) * conj(v2)
            double r1 = v1.x * cos_p - v1.y * sin_p;
            double i1 = v1.x * sin_p + v1.y * cos_p;

            sum.x += r1 * v2.x - i1 * (-v2.y);
            sum.y += r1 * (-v2.y) + i1 * v2.x;
        }
        d_props[idx] = sum;
    }
}

// Kernel 3: Generate the diagonal field-off propagators
__global__ void build_props_fieldoff_kernel(
    cuDoubleComplex* d_props, const double* d_H0,
    double h, int m, int offstep, int nsteps)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)(nsteps - offstep) * m * m;
    
    if (idx < total) {
        int k_rel = idx / (m * m);
        int k = k_rel + offstep;
        int rem = idx % (m * m);
        int col = rem / m;
        int row = rem % m;

        if (row == col) {
            double phase = -h * d_H0[row];
            d_props[(size_t)k * m * m + col * m + row] = {cos(phase), sin(phase)};
        } else {
            d_props[(size_t)k * m * m + col * m + row] = {0.0, 0.0};
        }
    }
}

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

    // pcc is a stacked array of contiguous N x N matrices.
    // We want the block for time J, delay block j. (j is 1-indexed).
    // THIS IS SAFE
    size_t pcc_idx = (size_t)J * (size_t)delay * (size_t)N2 + (size_t)(j-1) * (size_t)N2;
    const cuDoubleComplex* pcc_J_j_ptr = d_pcc + pcc_idx;
    // const cuDoubleComplex* pcc_J_j_ptr = d_pcc + J * (delay * N * N) + (j - 1) * N * N;

    for (int i = tid; i < N2; i += blockDim.x) {
        int u = i / N; // row index of R_k
        int v = i % N; // col index of R_k

        cuDoubleComplex sum = {0.0, 0.0};

        for (int x = 0; x < N; ++x) {
            // F(u, x) -> row u, col x
            cuDoubleComplex F_ux = pcc_J_j_ptr[x * N + u];
            cuDoubleComplex L_ux = {F_ux.x, -F_ux.y}; // Conjugate

            for (int y = 0; y < N; ++y) {
                // Bk is RowMajor
                cuDoubleComplex B_xy = Bk_ptr[x * N + y];

                // F^T(y, v) = F(v, y) -> row v, col y
                cuDoubleComplex F_vy = pcc_J_j_ptr[y * N + v];

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

class memoryModel
{
  const double h, T, tol, infreq, amp;
  const int ncyc, nsteps, delay, numthreads;
  const std::string inpath, outpath;
  const bool verbose, savetraj;
  double freq;                             // actual frequency
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
  std::vector<int> goodStatesVec;          // non-trivial indices of coefficient vector "coeffs"
  std::vector<int> goodCols;               // columns in "drcCI**2" space to retain
  Eigen::MatrixXcd kronprop;               // one-step field-free Kronecker product propagator
  // propagator chain cache (pcc)
  // the index for the std::vector is J, a discrete time step
  // each complex matrix in the cache is of size (ell*N) x N
  cuDoubleComplex* d_pcc = nullptr;
 
  public:
    //constructor
    memoryModel(double dt, double T, double infreq, double amp, int ncyc, int delay, double svdtol, int numthreads, std::string infile, std::string outpath, bool savetraj, int g0, int g1, bool verbose);

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

    // print bigmat to screen
    int bigmatPrint(Eigen::MatrixXcd bigmat);

    // propagate 1RDM with memory model for all delay values at once
    int qpropALLV2(void);
    int saveResults(void);
};

memoryModel::memoryModel(double dt, double T, double infreq, double amp, int ncyc, int delay, double svdtol, int numthreads, std::string infile, std::string outpath, bool savetraj, int g0, int g1, bool verbose)
   : h(dt), T(T), infreq(infreq), amp(amp), ncyc(ncyc), 
     delay(delay),
     tol(svdtol), numthreads(numthreads), 
     inpath(std::move(infile)), outpath(std::move(outpath)), 
     nsteps(static_cast<int>(std::ceil(T/h))), savetraj(savetraj), verbose(verbose)
{
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

  // Use negative frequency as a sentinel
  if (infreq < 0)
  {
    freq = (H0(g1) - H0(g0))/(2 * EIGEN_PI);
  }
  else freq = infreq;
  std::cout << "Frequency = " << freq << "\n";
  offstep = static_cast<int>(std::ceil(ncyc / (dt * freq)));
  // Prevent offstep from exceeding the total simulation length
  offstep = std::min(offstep, nsteps);
  std::cout << "Field will be on for " << offstep << " time steps\n";

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
    int m = drcCI;
    int m2 = m * m;

    // Initialize the CPU objects to standard sizes
    props.resize(nsteps, Eigen::MatrixXcd::Identity(m, m));
    for (int j = 0; j < m; ++j) coeffs(j, 0) = ic(j);

    // 1. Initialize GPU Handles
    cusolverDnHandle_t cusolverH = nullptr;
    cusolverDnCreate(&cusolverH);
    syevjInfo_t syevj_params = nullptr;
    cusolverDnCreateSyevjInfo(&syevj_params);
    cublasHandle_t cublasH = nullptr;
    cublasCreate(&cublasH);

    // 2. Allocate Device Memory
    double *d_H0, *d_dimatz;
    cuDoubleComplex *d_props, *d_coeffs;
    CHECK_CUDA(cudaMalloc(&d_H0, m * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_dimatz, m2 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_props, nsteps * m2 * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_coeffs, (nsteps + 1) * m * sizeof(cuDoubleComplex)));

    cuDoubleComplex *d_H = nullptr;
    double *d_W = nullptr;
    int *d_info = nullptr;
    if (offstep > 0) {
        CHECK_CUDA(cudaMalloc(&d_H, offstep * m2 * sizeof(cuDoubleComplex)));
        CHECK_CUDA(cudaMalloc(&d_W, offstep * m * sizeof(double)));
        CHECK_CUDA(cudaMalloc(&d_info, offstep * sizeof(int)));
    }

    // 3. Upload Static Data
    CHECK_CUDA(cudaMemcpy(d_H0, H0.data(), m * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_dimatz, dimatz.data(), m2 * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_coeffs, coeffs.col(0).data(), m * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // 4. Batched Field-On Calculations
    if (offstep > 0) {
        // Build H matrices
        size_t threads_H = 256;
        size_t blocks_H = (offstep * m2 + threads_H - 1) / threads_H;
        build_H_batched_kernel<<<blocks_H, threads_H>>>(d_H, d_H0, d_dimatz, amp, freq, h, EIGEN_PI, m, offstep);

        // Batched Eigendecomposition Workspace
        int lwork = 0;
        cusolverDnZheevjBatched_bufferSize(
            cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
            m, d_H, m, d_W, &lwork, syevj_params, offstep);

        cuDoubleComplex* d_work;
        CHECK_CUDA(cudaMalloc(&d_work, lwork * sizeof(cuDoubleComplex)));

        // Run Batched EigenSolver (d_H is overwritten with Eigenvectors V)
        cusolverDnZheevjBatched(
            cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
            m, d_H, m, d_W, d_work, lwork, d_info, syevj_params, offstep);

        // Form Propagators
        size_t blocks_P = (offstep * m2 + threads_H - 1) / threads_H;
        build_props_batched_kernel<<<blocks_P, threads_H>>>(d_props, d_H, d_W, h, m, offstep);

        cudaFree(d_work);
    }

    // 5. Field-Off Calculations
    if (nsteps > offstep) {
        size_t threads_off = 256;
        size_t blocks_off = ((nsteps - offstep) * m2 + threads_off - 1) / threads_off;
        build_props_fieldoff_kernel<<<blocks_off, threads_off>>>(d_props, d_H0, h, m, offstep, nsteps);
    }

    // 6. Sequential State Propagation
    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta = {0.0, 0.0};
    for (int k = 0; k < nsteps; ++k) {
        cublasZgemv(cublasH, CUBLAS_OP_N, m, m,
                    &alpha, d_props + k * m2, m,
                    d_coeffs + k * m, 1,
                    &beta, d_coeffs + (k + 1) * m, 1);
    }

    // 7. Download Results Directly to CPU Objects
    std::vector<cuDoubleComplex> h_props_flat(nsteps * m2);
    CHECK_CUDA(cudaMemcpy(h_props_flat.data(), d_props, nsteps * m2 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    for (int k = 0; k < nsteps; ++k) {
        memcpy(props[k].data(), &h_props_flat[k * m2], m2 * sizeof(cuDoubleComplex));
    }

    CHECK_CUDA(cudaMemcpy(coeffs.data(), d_coeffs, (nsteps + 1) * m * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    // 8. Cleanup
    cudaFree(d_H0); cudaFree(d_dimatz); cudaFree(d_props); cudaFree(d_coeffs);
    if (offstep > 0) {
        cudaFree(d_H); cudaFree(d_W); cudaFree(d_info);
    }
    cusolverDnDestroySyevjInfo(syevj_params);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);

    havecoeffs = true;

    // 9. Save Trajectory (from original code)
    if (savetraj)
    {
        std::filesystem::path p(inpath);
        std::string stem = p.stem().string();
        std::string filename = stem + "_" + std::to_string(h);
        filename += "_" + std::to_string(freq);
        filename += "_" + std::to_string(amp);
        filename += "_" + std::to_string(ncyc) + "_coeffs.txt";
        std::filesystem::path dir(outpath);
        std::filesystem::path outfile = dir / filename;
        std::ofstream out(outfile);
        for (int k = 0; k <= nsteps; ++k)
        {
            for (int l = 0; l < drcCI; ++l)
                out << coeffs(l, k).real() << "+" << coeffs(l, k).imag() << "j" << (l < (drcCI - 1) ? "," : "");
            out << "\n";
        }
    }
    return 0;
}

int memoryModel::exact1RDMS(void)
{
    if (!havecoeffs) {
        std::cout << "TDCI coefficients have not been computed yet!\n";
        return 1;
    }

    // Matrix dimensions
    int m = drcCI;
    int m2 = m * m;
    int steps = nsteps + 1;

    // 1. Allocate GPU memory
    cuDoubleComplex* d_coeffs;
    cuDoubleComplex* d_true1rdms;
    double* d_Bmat;

    CHECK_CUDA(cudaMalloc(&d_coeffs, m * steps * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_true1rdms, drc2 * steps * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_Bmat, drc2 * m2 * sizeof(double)));

    // 2. Upload Data
    CHECK_CUDA(cudaMemcpy(d_coeffs, coeffs.data(), m * steps * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_Bmat, Bmat.data(), drc2 * m2 * sizeof(double), cudaMemcpyHostToDevice));

    // 3. Launch the Fused Kernel
    // We launch exactly one block per time step.
    int blocks = steps;

    // We want roughly one thread per row of true1rdms (drc2).
    // We round up to the nearest multiple of 32 for warp efficiency.
    int threads = ((drc2 + 31) / 32) * 32;
    if (threads > 1024) threads = 1024; // Safety cap

    // Dynamically allocate shared memory based on drcCI
    size_t shared_mem_size = m * sizeof(cuDoubleComplex);

    compute_exact1rdms_kernel<<<blocks, threads, shared_mem_size>>>(
        d_true1rdms, d_coeffs, d_Bmat, m, drc2, nsteps
    );
    CHECK_CUDA(cudaDeviceSynchronize());

    // 4. Download and Cleanup
    CHECK_CUDA(cudaMemcpy(true1rdms.data(), d_true1rdms, drc2 * steps * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    cudaFree(d_coeffs);
    cudaFree(d_true1rdms);
    cudaFree(d_Bmat);

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

  // 1. Compute all squared row norms simultaneously
  const double thresh_sq = 1e-20;
  Eigen::VectorXd sq_norms = coeffs.rowwise().squaredNorm();

  // 2. Clear and pre-allocate our state vector
  goodStatesVec.clear();
  goodStatesVec.reserve(drcCI);

  N = 0;
  for (int j = 0; j < drcCI; ++j)
  {
    if (sq_norms(j) >= thresh_sq)
    {
      goodStatesVec.push_back(j);
      ++N;
    }
  }

  // 3. Clear, pre-allocate, and generate goodCols
  goodCols.clear();
  goodCols.reserve(N * N);

  for (int r : goodStatesVec)
    for (int c : goodStatesVec)
      goodCols.push_back(r * drcCI + c);

  N2 = N * N;
  if (verbose)
    std::cout << "Retaining " << N2 << " or " << goodCols.size() << " entries\n";

  // reduce the B matrix
  BmatR = Bmat(Eigen::placeholders::all, goodCols);

  // filter all the propagators
  fprops.resize(nsteps, Eigen::MatrixXcd::Identity(N, N));
  #pragma omp parallel for
  for (int k = 0; k < nsteps; ++k)
    fprops[k].noalias() = props[k](goodStatesVec, goodStatesVec);

  // we need this field-free Kronecker propagator below
  std::complex<double> coeff = -1.0i * h;
  Eigen::VectorXcd factor = coeff * H0.array();
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
    if (!filtered) {
        std::cout << "Filtered indices and propagators must be computed first!\n";
        return 1;
    }

    int upper = delay + offstep;
    int pccsize = (nsteps < upper) ? nsteps : upper;

    // 1. Setup handles and allocations
    cublasHandle_t cublasH = nullptr;
    cublasCreate(&cublasH);

    // Flatten fprops for upload
    std::vector<cuDoubleComplex> fprops_flat((nsteps + 1) * N * N, {0.0, 0.0});
    for (int k = 0; k < nsteps; ++k) {
        memcpy(&fprops_flat[k * N * N], fprops[k].data(), N * N * sizeof(cuDoubleComplex));
    }

    size_t pcc_elements_per_J = (size_t)delay * N * N;
    size_t total_pcc_elements = (size_t)(pccsize + 1) * pcc_elements_per_J;

    if (verbose) std::cout << "Allocating " << (total_pcc_elements * 16.0 / 1e9) << " GB for d_pcc on GPU...\n";

    cuDoubleComplex *d_fprops;
    CHECK_CUDA(cudaMalloc(&d_fprops, fprops_flat.size() * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_pcc, total_pcc_elements * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMemset(d_pcc, 0, total_pcc_elements * sizeof(cuDoubleComplex)));

    CHECK_CUDA(cudaMemcpy(d_fprops, fprops_flat.data(), fprops_flat.size() * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    // 2. Initialize j=1 blocks for all J
    int threads = 256;
    size_t blocks = ((pccsize - delay + 1) * N * N + threads - 1) / threads;
    init_pcc_j1_kernel<<<blocks, threads>>>(d_pcc, d_fprops, delay, N, pccsize, delay);
    CHECK_CUDA(cudaDeviceSynchronize());

    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta  = {0.0, 0.0};

    // 3. Sequential Chain for J = delay
    for (int j = 2; j <= delay; ++j) {
        cuDoubleComplex* d_C = d_pcc + delay * pcc_elements_per_J + (j - 1) * N * N;
        cuDoubleComplex* d_A = d_pcc + delay * pcc_elements_per_J + (j - 2) * N * N;
        cuDoubleComplex* d_B = d_fprops + (delay - j) * N * N;

        cublasZgemm(cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
                    N, N, N, &alpha,
                    d_A, N, d_B, N,
                    &beta, d_C, N);
    }
    CHECK_CUDA(cudaDeviceSynchronize());

    // 4. Batched Generation for J > delay
    for (int J = delay + 1; J <= pccsize; ++J) {

        cuDoubleComplex* d_A = d_fprops + (J - 1) * N * N;
        cuDoubleComplex* d_B_array = d_pcc + (J - 1) * pcc_elements_per_J;
        cuDoubleComplex* d_C_array = d_pcc + J * pcc_elements_per_J + N * N; // Starts at j=2

        cublasZgemmStridedBatched(
            cublasH, CUBLAS_OP_N, CUBLAS_OP_N,
            N, N, N,
            &alpha,
            d_A, N, 0,                      // strideA = 0 (Broadcast same fprop)
            d_B_array, N, N * N,            // strideB = N*N (Shift by one block)
            &beta,
            d_C_array, N, N * N,            // strideC = N*N (Shift by one block)
            delay - 1                       // Number of batches
        );
    }
    CHECK_CUDA(cudaDeviceSynchronize());

    // 6. Cleanup
    cudaFree(d_fprops);
    cublasDestroy(cublasH);

    builtpcc = true;
    if (verbose)
        std::cout << "Built propagator chain cache on GPU!\n";
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
    cuDoubleComplex *d_pred1rdms, *d_BmatR, *d_fprops;
    CHECK_CUDA(cudaMalloc(&d_pred1rdms, drc2 * (nsteps + 1) * sizeof(cuDoubleComplex)));
    CHECK_CUDA(cudaMalloc(&d_BmatR, drc2 * N2 * sizeof(cuDoubleComplex)));
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

    cudaFree(d_pred1rdms); cudaFree(d_BmatR); cudaFree(d_fprops);
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
  std::string filename = stem + "_" + std::to_string(h);
  filename += "_" + std::to_string(delay);
  filename += "_" + std::to_string(freq);
  filename += "_" + std::to_string(amp);
  filename += "_" + std::to_string(ncyc) + ".txt";
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
  ("savetraj", "Whether to save TDCI trajectory to disk", cxxopts::value<bool>()->default_value("false"))
  ("g0", "Lower eigenvalue/frequency (if specifying gap)", cxxopts::value<int>()->default_value("0"))
  ("g1", "Higher eigenvalue/frequency (if specifying gap)", cxxopts::value<int>()->default_value("1"))
  ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
  
  auto result = options.parse(argc, argv);
  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0;
  }
  double tol;
  std::vector<double> timeparams, fieldparams;
  int delay, g0, g1;
  std::string inpath, outpath;
  bool verbose, savetraj;
  
  // Required parameters
  getRequired(result, "time",    timeparams,  "Must specify time-stepping parameters dt,T!");
  getRequired(result, "field",   fieldparams, "Must specify field parameters freq,amp,ncyc!");
  getRequired(result, "delay",   delay,       "Must specify delay!");
  getRequired(result, "infile",  inpath,      "Must specify input file!");
  getRequired(result, "outpath", outpath,     "Must specify output path!");
  
  // Optional parameter
  getOptional(result, "tol", tol, 1e-6);
  getOptional(result, "savetraj", savetraj, false);
  getOptional(result, "g0", g0, 0);
  getOptional(result, "g1", g1, 0);
  getOptional(result, "verbose", verbose, false);
  
  // pass the raw parameters into the constructor, in keeping with the other (time & field) params
  int ncyc = static_cast<int>(std::round(fieldparams[2]));
  memoryModel mm(timeparams[0], timeparams[1], fieldparams[0], fieldparams[1], ncyc, delay, tol, num_threads, inpath, outpath, savetraj, g0, g1, verbose);
  Eigen::VectorXcd ic(mm.getdrcCI());
  ic.setZero();
  ic[0] = 1.0;
  std::cout << "Just before tdseProp\n";
  mm.tdseProp(ic);
  std::cout << "Just after tdseProp\n";
  mm.exact1RDMS();
  std::cout << "Just after exact1RDMS\n";
  mm.filterIndices();
  std::cout << "Just after filterIndices\n";
  mm.buildPCC();
  std::cout << "Just after buildPCC\n";
  mm.qpropALLV2();
  mm.saveResults();
  return 0;
}

