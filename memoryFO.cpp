// memory model
// field-on case

#define EIGEN_USE_BLAS 
#define EIGEN_USE_LAPACKE

#include <omp.h>
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

using namespace std::complex_literals;

// Function to compute the Moore–Penrose pseudoinverse
Eigen::MatrixXcd pseudoInverse(const Eigen::MatrixXcd &mat, double tolerance)
{
  // Compute SVD: mat = U * Σ * Vᵀ
  Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd( mat );

  const Eigen::VectorXd &singularValues = svd.singularValues();
  Eigen::VectorXd singularValuesInv(singularValues.size());

  // Invert singular values with tolerance to avoid division by zero
  for (int i = 0; i < singularValues.size(); ++i)
  {
      if (singularValues(i) > tolerance)
          singularValuesInv(i) = 1.0 / singularValues(i);
      else
          singularValuesInv(i) = 0.0;
  }

  // Pseudoinverse formula: V * Σ⁺ * U^\dagger
  return svd.matrixV() * singularValuesInv.asDiagonal() * svd.matrixU().adjoint();
}

class memoryModel
{
  using MatrixXcdRowMajor = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  const double h, T, tol, freq, amp;
  const int ncyc, nsteps, delaystart, delaystep, numdelays, numthreads;
  const std::string outfile;
  int drc, drc2, drcCI, drcCI2;
  int N, N2;                               // N = number of good states
  int offstep;                             // time step at which the field is switched off
  bool havecoeffs = false;                 // have we computed time-dependent coefficients yet?
  bool have1rdms = false;                  // have we computed 1-RDMS yet?
  bool filtered = false;                   // have we filtered out bad indices yet?
  Eigen::VectorXd H0;                      // core Hamiltonian
  Eigen::MatrixXd dimatz;                  // dipole moment matrix in the z-direction
  Eigen::MatrixXd Bmat;                    // B tensor, in matrix form, reshaped to size (drc2 x drcCI2)
  MatrixXcdRowMajor BmatR;                 // reduced B tensor, in matrix form, reshaped to size (drc2 x N2)
  Eigen::MatrixXcd coeffs;                 // time-dependent coefficients
  std::vector<Eigen::MatrixXcd> props;     // all time-dependent propagators "exp(-i H_n dt)"
  std::vector<Eigen::MatrixXcd> fprops;    // all time-dependent propagators "exp(-i H_n dt)"
  Eigen::MatrixXcd true1rdms;              // ground truth 1-RDMs
  std::vector<Eigen::MatrixXcd> pred1rdms; // predicted 1-RDMs (at each of the delay values in delayrange)s
  std::unordered_set<int> goodStates;      // non-trivial indices of coefficient vector "coeffs"
  std::vector<int> goodCols;               // columns in "drcCI**2" space to retain
  MatrixXcdRowMajor bigmat;                // stacked matrix that gets pseudoinverted at each step
  // state variables for lastblock
  int ell;                                 // what delay are we currently working on?
  Eigen::MatrixXcd Amat;                   // running chain of propagators
  // return variables for mySVD
  Eigen::MatrixXcd u1, v1;
  Eigen::VectorXd s1;
  // return variables for updateSVD
  Eigen::MatrixXcd u2, v2;
  Eigen::VectorXd s2;
  // temporary storage
  std::vector<Eigen::MatrixXcd> thread_temps;
  
  public:
    //constructor
    memoryModel(double dt, double T, double freq, double amp, int ncyc, int delaystart, int delaystep, int numdelays,
                double svdtol, int numthreads, std::string infile, std::string ioutfile);

    int getdrc(void) { return drc; }
    int getdrc2(void) { return drc2; }
    int getdrcCI(void) { return drcCI; }
    int getN(void) { return N; }
    int getnsteps(void) { return nsteps; }
    Eigen::MatrixXcd getTrue1rdms(void) { return true1rdms; }
    Eigen::MatrixXcd getPred1rdms(int j) { return pred1rdms[j]; }

    // solve TDSE forward in time starting from a particular initial condition,
    // saving propagators as we go
    int tdseProp(Eigen::VectorXcd& ic);

    // compute and store all exact 1RDMS
    int exact1RDMS(void);

    // figure out which indices (among the drcCI**2 linear indices) we should retain/delete
    int filterIndices(void);

    // form bigmat at particular point in time for particular delay value
    int bigmatBuild(int J, int jell);

    // print bigmat to screen
    int bigmatPrint(void);

    // this function advances the chain of propagators from ell to lastell
    // it forms and returns the next block in bigmat 
    // note that this function alters the state of both ell and Amat
    Eigen::MatrixXcd nextBlock(int J, int lastell);

    // here we incrementally update the SVD to account for the new block
    // this function relies on ell being up to date
    // it assumes that u1,s1,v1 contains an OLD SVD that will be updated and
    // stored in class/member variables u2,s2,v2
    int updateSVD(const Eigen::MatrixXcd& r);

    // a little helper function that computes the SVD of mat and stores it in u1,s1,v1
    int mySVD(const Eigen::MatrixXcd &mat);

    // compute the pseudoinverse multiplied by vec (with given tolerance)
    // the pseudoinverse is computed using u1, s1, v1
    Eigen::VectorXcd pinvvec(const Eigen::VectorXcd &vec, double tolerance);

    // propagate 1RDM with memory model for all delay values at once
    int qpropALL(void);
    int qpropALLV2(void);
};


memoryModel::memoryModel(double dt, double T, double freq, double amp, int ncyc, int delaystart, int delaystep, int numdelays,
                         double svdtol, int numthreads, std::string infile, std::string ioutfile)
   : h(dt), T(T), freq(freq), amp(amp), ncyc(ncyc), delaystart(delaystart), delaystep(delaystep), numdelays(numdelays), 
     tol(svdtol), numthreads(numthreads), outfile(std::move(ioutfile)), nsteps(static_cast<int>(std::ceil(T/h)))
{
  offstep = static_cast<int>(std::ceil(ncyc / (dt * freq)));
  std::cout << "Field will be on for " << offstep << " time steps\n";
  // Load the entire .npz file into a map-like structure
  cnpy::npz_t my_npz = cnpy::npz_load(infile);

  // Load ham
  cnpy::NpyArray arr = my_npz["ham"];
  double* hamdata = arr.data<double>();
  size_t length = arr.shape[0];
  Eigen::Map<Eigen::VectorXd> H0_map(hamdata, length);
  H0 = H0_map;
  drcCI = (int) length;
  drcCI2 = drcCI*drcCI;
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
  std::cout << "drc = " << drc << "\n";
  Eigen::Map<Eigen::VectorXd> Bten(Btendata, length);
  Eigen::Map<const Eigen::MatrixXd> Bmat_map(Bten.data(), drc2, drcCI2);
  Bmat = Bmat_map;

  // initialize various objects
  coeffs.setZero(drcCI, nsteps+1);
  true1rdms.setZero(drc2, nsteps+1);
  pred1rdms.resize(numdelays, Eigen::MatrixXcd::Zero(drc2, nsteps+1));
}


int memoryModel::tdseProp(Eigen::VectorXcd& ic)
{
  double field = 0.0;
  Eigen::MatrixXcd prop;
  Eigen::VectorXcd D, Dexp;
  std::complex<double> scalarfac = -1.0i * h;

  // std::cout << "Received initial condition " << ic << "\n";
  
  // copy initial condition into coeffs
  for (int j=0; j<drcCI; ++j)
    coeffs(j, 0) = ic(j);
    
  // initialize propagators to be the identity
  props.resize(nsteps, Eigen::MatrixXcd::Identity(drcCI, drcCI));

  for (int k=0; k<nsteps; ++k)
  {
    Eigen::MatrixXcd H = H0.asDiagonal();
    if (k < offstep)
    {
      field = amp * std::sin(2 * EIGEN_PI * freq * k * h);
      // std::cout << "Time step j = " << k << "; field strength = " << field << "\n";
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
    // std::cout << "coeffs.col(" << (k+1) << ") = " << coeffs.col(k+1) << "\n";
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
    true1rdms.col(k) = Bmat * op.transpose().reshaped();
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
    // std::cout << "Row " << j << " has norm = " << rownorm << "\n";
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
  std::cout << "Retaining " << N2 << " or " << goodCols.size() << " entries\n";

  // reduce the B matrix
  BmatR = Bmat(Eigen::placeholders::all, goodCols);

  // filter all the propagators
  fprops.resize(nsteps, Eigen::MatrixXcd::Identity(N, N));
  std::vector<int> goodStatesVec(goodStates.begin(), goodStates.end());
  std::sort(goodStatesVec.begin(), goodStatesVec.end());

  #pragma omp parallel for
  for (int k=0; k<nsteps; ++k)
    fprops[k] = props[k](goodStatesVec, goodStatesVec);

  // allocate temporary storage
  thread_temps.resize(numthreads, Eigen::MatrixXcd(N, N));    
  
  filtered = true;
  return 0;
}

// form bigmat at particular point in time for particular delay value
// this function performs double duty as the "starting block" for the incremental approach
int memoryModel::bigmatBuild(int J, int jell)
{
  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return 1;
  }  
  // resize bigmat and store first block
  bigmat.setZero( (jell+1)*drc2, N2 );
  bigmat.block(0, 0, drc2, N2) = BmatR;

  // this will hold our propagator chain
  Amat.setIdentity( N, N );
    
  for (int j=1; j<=jell; ++j)
  {
    Amat = Amat * fprops[J - j];
    MatrixXcdRowMajor Result(drc2, N2);
    #pragma omp parallel for schedule(dynamic)
    for (int k=0; k<drc2; ++k)
    {
      // Get the ID of the thread currently running this iteration
      int tid = omp_get_thread_num();
      Eigen::Map<MatrixXcdRowMajor> BmatRk(BmatR.row(k).data(), N, N);
      Eigen::Map<MatrixXcdRowMajor> Rk(Result.row(k).data(), N, N);
      // 3. Use .noalias() on BOTH steps to force Eigen to write directly 
      // into our pre-allocated memory, preventing any hidden allocations!
      thread_temps[tid].noalias() = Amat.conjugate() * BmatRk;
      Rk.noalias() = thread_temps[tid] * Amat.transpose();
    }
    bigmat.block(j*drc2, 0, drc2, N2) = Result;
  }
  // update class data to reflect that Amat includes a chain up till this point
  ell = jell;
  return 0;
}

int memoryModel::bigmatPrint(void)
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

// this function advances the chain of propagators from ell to lastell
// it forms and returns the next block in bigmat 
// note that this function alters the state of both ell and Amat

// technical note: let us assume that the first block (computed by bigmatBuild) includes ell=0 to ell=delaystart
// then each successive block (computed by this function) will go from the past ell+1 to lastell,
// so that if lastell = ell+delaystep, we get exactly "delaystep" new pieces per block... 
Eigen::MatrixXcd memoryModel::nextBlock(int J, int lastell)
{
  // set up a matrix to hold the next block
  Eigen::MatrixXcd nb( (lastell - ell)*drc2, N2 );

  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return nb;
  }  
  
  for (int j=ell+1; j<=lastell; ++j)
  {
    Amat = Amat * fprops[J - j];
    MatrixXcdRowMajor Result(drc2, N2);    
    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < drc2; ++k)
    {
      int tid = omp_get_thread_num();
      Eigen::Map<MatrixXcdRowMajor> BmatRk(BmatR.row(k).data(), N, N);
      Eigen::Map<MatrixXcdRowMajor> Rk(Result.row(k).data(), N, N);
      thread_temps[tid].noalias() = Amat.conjugate() * BmatRk;
      Rk.noalias() = thread_temps[tid] * Amat.transpose();
    }
    nb.block((j-(ell+1))*drc2, 0, drc2, N2) = Result;
  }
  ell = lastell;
  return nb;
}

// here we incrementally update the SVD to account for the new block
// this function relies on ell being up to date
// note that this function alters u2, s2, v2
int memoryModel::updateSVD(const Eigen::MatrixXcd& r)
{
  int N = v1.rows();  // The maximum Hilbert space dimension
  int k = s1.size();
  int m = r.rows();
  
  Eigen::MatrixXcd v1dag = v1.adjoint();
  Eigen::MatrixXcd r_par = r * v1;
  Eigen::MatrixXcd rorth = r - r_par * v1dag;
  
  // Calculate a dynamic noise threshold based on the incoming block
  double r_norm = r.norm();
  double tol = (r_norm > 1e-14) ? 1e-12 * r_norm : 1e-12;
    
  if (k >= N || rorth.norm() < tol)
  {
    // --- TALL REGIME / NO EXPANSION ---
    Eigen::MatrixXcd kmat(k + m, k);
    kmat.topRows(k) = s1.cast<std::complex<double>>().asDiagonal();
    kmat.bottomRows(m) = r_par;
    Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(kmat);
    
    int rank_to_keep = std::min(static_cast<int>(svd.singularValues().size()), N);
    this->s2 = svd.singularValues().head(rank_to_keep);
    
    Eigen::MatrixXcd utilde = svd.matrixU().leftCols(rank_to_keep);
    Eigen::MatrixXcd v_tilde = svd.matrixV().leftCols(rank_to_keep);

    this->u2.resize(u1.rows() + m, rank_to_keep);
    this->u2.topRows(u1.rows()) = u1 * utilde.topRows(k);
    this->u2.bottomRows(m) = utilde.bottomRows(m);
    
    this->v2 = v1 * v_tilde;
  }
  else
  {
    // --- EXPANSION REGIME ---
    Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd_rorth(rorth);
    
    // Count meaningful new dimensions above the noise floor
    int valid_dims = (svd_rorth.singularValues().array() > tol).count();
    
    // CRITICAL GUARDRAIL: Cap the added dimensions so k + new_dims <= N
    int new_dims = std::min(valid_dims, N - k);
    
    if (new_dims == 0)
    {
      // Everything was noise or subspace is physically full
      Eigen::MatrixXcd kmat(k + m, k);
      kmat.topRows(k) = s1.cast<std::complex<double>>().asDiagonal();
      kmat.bottomRows(m) = r_par;
      
      Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(kmat);
      
      int rank_to_keep = std::min(static_cast<int>(svd.singularValues().size()), N);
      this->s2 = svd.singularValues().head(rank_to_keep);
      
      Eigen::MatrixXcd utilde = svd.matrixU().leftCols(rank_to_keep);
      Eigen::MatrixXcd v_tilde = svd.matrixV().leftCols(rank_to_keep);

      this->u2.resize(u1.rows() + m, rank_to_keep);
      this->u2.topRows(u1.rows()) = u1 * utilde.topRows(k);
      this->u2.bottomRows(m) = utilde.bottomRows(m);
      
      this->v2 = v1 * v_tilde;
    }
    else
    {
      // Extract exactly the safe number of new orthogonal basis vectors.
      // Note: Eigen returns V directly, so no conj().T() is needed like in NumPy.
      Eigen::MatrixXcd q = svd_rorth.matrixV().leftCols(new_dims);
      
      // Re-orthogonalize against v1 to remove floating point leakage
      for (int i = 0; i < 2; ++i)
      {
        q = q - v1 * (v1dag * q);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(q);
        // Extract thin Q
        q = qr.householderQ() * Eigen::MatrixXcd::Identity(q.rows(), q.cols());
      }
      
      Eigen::MatrixXcd z = rorth * q;
      
      // Build the kmat block matrix
      Eigen::MatrixXcd kmat = Eigen::MatrixXcd::Zero(k + m, k + new_dims);
      kmat.topLeftCorner(k, k) = s1.cast<std::complex<double>>().asDiagonal();
      kmat.bottomLeftCorner(m, k) = r_par;
      kmat.bottomRightCorner(m, new_dims) = z;
      
      Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(kmat);
      
      // Strict truncation to N
      int rank_to_keep = std::min(static_cast<int>(svd.singularValues().size()), N);
      this->s2 = svd.singularValues().head(rank_to_keep);
      
      Eigen::MatrixXcd utilde = svd.matrixU().leftCols(rank_to_keep);
      Eigen::MatrixXcd v_tilde = svd.matrixV().leftCols(rank_to_keep);

      this->u2.resize(u1.rows() + m, rank_to_keep);
      this->u2.topRows(u1.rows()) = u1 * utilde.topRows(k);
      this->u2.bottomRows(m) = utilde.bottomRows(m);
      
      Eigen::MatrixXcd v1_q(v1.rows(), v1.cols() + q.cols());
      v1_q << v1, q;
      this->v2 = v1_q * v_tilde;
    }
  }
  // --- FINAL PERIODIC CLEANUP ---
  if (this->ell % 50 == 0)
  {
    // Thin QR of u2
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr_u(this->u2);
    Eigen::MatrixXcd u_clean = qr_u.householderQ() * Eigen::MatrixXcd::Identity(this->u2.rows(), this->u2.cols());
    Eigen::MatrixXcd R_u = qr_u.matrixQR().topRows(this->u2.cols()).triangularView<Eigen::Upper>();
    
    // Thin QR of v2
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr_v(this->v2);
    Eigen::MatrixXcd v_clean = qr_v.householderQ() * Eigen::MatrixXcd::Identity(this->v2.rows(), this->v2.cols());
    Eigen::MatrixXcd R_v = qr_v.matrixQR().topRows(this->v2.cols()).triangularView<Eigen::Upper>();
    
    // Absorb into center core
    Eigen::MatrixXcd core = R_u * this->s2.cast<std::complex<double>>().asDiagonal() * R_v.adjoint();
    
    Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd_core(core);
    
    this->s2 = svd_core.singularValues();
    this->u2 = u_clean * svd_core.matrixU();
    this->v2 = v_clean * svd_core.matrixV();
  }
  return 0;
}

// compute the SVD of mat and store it in u1, s1, v1
int memoryModel::mySVD(const Eigen::MatrixXcd &mat)
{
  // Compute SVD: mat = U * Σ * Vᵀ
  Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd( mat );
  this->u1 = svd.matrixU();
  this->s1 = svd.singularValues();
  this->v1 = svd.matrixV();
  return 0;
}

Eigen::VectorXcd memoryModel::pinvvec(const Eigen::VectorXcd &vec, double tolerance)
{
  Eigen::VectorXd singularValuesInv(s1.size());
  // Invert singular values with tolerance to avoid division by zero
  for (int i = 0; i < s1.size(); ++i)
  {
      if (s1(i) > tolerance)
          singularValuesInv(i) = 1.0 / s1(i);
      else
          singularValuesInv(i) = 0.0;
  }
  // Pseudoinverse formula: V * Σ⁺ * U^\dagger
  Eigen::VectorXcd temp = u1.adjoint() * vec;
  Eigen::VectorXcd temp2 = singularValuesInv.asDiagonal() * temp;
  return v1 * temp2;
}

// propagate 1RDM with memory model for all delay values at once
int memoryModel::qpropALL(void)
{
  // initialize pred1rdms
  for (int i=0; i<numdelays; ++i)
  {
    int jell = delaystart + i*delaystep;
    // copy in some ground truth 1RDMs
    for (int k=0; k<=jell; ++k)
      pred1rdms[i].col(k) = true1rdms.col(k);    
  }  

  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return 1;
  }
  
  for (int J=delaystart; J<nsteps; ++J)
  {
    std::cout << "Propagating 1RDM at time step " << J << "\n";
    for (int i=0; i<numdelays; ++i)
    {
      int jell = delaystart + i*delaystep;
      // catch the case where we cannot propagate the memory model because we don't have enough memory yet!
      if ((J-jell) < 0)
      {
        pred1rdms[i].col(J+1) = true1rdms.col(J+1);
      }
      else
      {
        if (i==0)
        {
          // compute the bigmat once
          // note that this resets Amat and our internal ell
          auto t1 = std::chrono::steady_clock::now();
          bigmatBuild(J, jell);
          auto t2 = std::chrono::steady_clock::now();
          mySVD(bigmat);
          auto t3 = std::chrono::steady_clock::now();
          std::chrono::duration<double> e1 = t2 - t1;
          std::chrono::duration<double> e2 = t3 - t2;
          /*
          std::cout << "bigmatBuild: " << e1.count() << " seconds\n";
          std::cout << "mySVD: " << e2.count() << " seconds\n"; 
          */
        }
        else
        {
          // grab the next block
          // note that this updates Amat and our internal ell
          auto t1 = std::chrono::steady_clock::now();
          Eigen::MatrixXcd nb = nextBlock(J, jell);
          auto t2 = std::chrono::steady_clock::now();
          // computed updated SVD and note it gets placed in u2,s2,v2
          updateSVD(nb);
          auto t3 = std::chrono::steady_clock::now();
          // manually update u1,s1,v1
          std::chrono::duration<double> e1 = t2 - t1;
          std::chrono::duration<double> e2 = t3 - t2;
          /*
          std::cout << "nextBlock: " << e1.count() << " seconds\n";
          std::cout << "updateSVD: " << e2.count() << " seconds\n";  
          */
          u1 = u2;
          s1 = s2;
          v1 = v2;
        }
        Eigen::MatrixXcd qhist = pred1rdms[i].block(0, J-jell, drc2, jell+1).rowwise().reverse();
        Eigen::VectorXcd PreconVec = pinvvec(qhist.reshaped(), this->tol);
        Eigen::MatrixXcd Precon = PreconVec.reshaped(N, N);
        pred1rdms[i].col(J+1) = BmatR * (fprops[J].adjoint() * Precon * fprops[J].transpose()).reshaped();
      }
    }
  }
  return 0;
}

// propagate 1RDM with memory model for all delay values at once
int memoryModel::qpropALLV2(void)
{
  std::vector<Eigen::MatrixXcd> rdmprops;
  rdmprops.resize(numdelays);

  std::vector<int> goodStatesVec(goodStates.begin(), goodStates.end());
  std::sort(goodStatesVec.begin(), goodStatesVec.end());
  std::complex<double> coeff = -1.0i * h;
  Eigen::VectorXcd factor = coeff*H0.array();
  Eigen::VectorXcd diagvals = factor.array().exp();
  diagvals = diagvals(goodStatesVec);
  Eigen::MatrixXcd outer_product = diagvals * diagvals.adjoint();
  Eigen::VectorXcd kron_diag = outer_product.reshaped<Eigen::RowMajor>();
  Eigen::MatrixXcd kronprop = kron_diag.asDiagonal();
  // std::cout << "Constructed kronprop of size " << kronprop.rows() << ", " << kronprop.cols() << "\n";
  
  // remove this as/when we eventually replace qpropALL with qpropALLV2
  pred1rdms.resize(numdelays, Eigen::MatrixXcd::Zero(drc2, nsteps+1));
  
  // initialize pred1rdms
  for (int i=0; i<numdelays; ++i)
  {
    int jell = delaystart + i*delaystep;
    // copy in some ground truth 1RDMs
    for (int k=0; k<=jell; ++k)
      pred1rdms[i].col(k) = true1rdms.col(k);    
  }

  if (!filtered)
  {
    std::cout << "Filtered indices and propagators must be computed first!\n";
    return 1;
  }

  // suppose that J-jell >= offstep
  // then all the propagators we need to go backwards in time are in fact time-independent
  // this means that for J >= jell + offstep we can build the 1-RDM propagator once and for all!
  int firsttimedep = 0;
  for (int J=delaystart; J<nsteps; ++J)
  {
    std::cout << "Propagating 1RDM at time step " << J << "\n";
    for (int i=0; i<numdelays; ++i)
    {
      int jell = delaystart + i*delaystep;
      // catch the case where we cannot propagate the memory model because we don't have enough memory yet!
      if ((J-jell) < 0)
      {
        pred1rdms[i].col(J+1) = true1rdms.col(J+1);
      }
      else
      {
        bool needPropagation=true;
        if ((J-jell) >= offstep)
        {
          if ((J-jell) == offstep)
          {
            // std::cout << "Building rdmprop for delay " << jell << "\n";
            bigmatBuild(J, jell);
            Eigen::MatrixXcd mypinv = pseudoInverse(this->bigmat, this->tol);
            rdmprops[i] = BmatR * kronprop * mypinv;
            firsttimedep++;
          }
          Eigen::MatrixXcd qhist = pred1rdms[i].block(0, J-jell, drc2, jell+1).rowwise().reverse();
          // std::cout << "Using rdmprop for delay " << jell << "\n";
          pred1rdms[i].col(J+1) = rdmprops[i] * qhist.reshaped();
          needPropagation=false;
        }
        if (i==firsttimedep)
        {
          // if (firsttimedep >= 1)
          //     std::cout << "Entering i==firsttimedep block...\n";
          // compute the bigmat once
          // note that this resets Amat and our internal ell
          auto t1 = std::chrono::steady_clock::now();
          bigmatBuild(J, jell);
          auto t2 = std::chrono::steady_clock::now();
          mySVD(bigmat);
          auto t3 = std::chrono::steady_clock::now();
          std::chrono::duration<double> e1 = t2 - t1;
          std::chrono::duration<double> e2 = t3 - t2;
          /*
          std::cout << "bigmatBuild: " << e1.count() << " seconds\n";
          std::cout << "mySVD: " << e2.count() << " seconds\n";  
          */
        }
        else if ((i > firsttimedep) && (firsttimedep < numdelays))
        {
          // if (firsttimedep >= 1)
          //     std::cout << "Entering firsttimedep < numdelays block...\n";
          // grab the next block
          // note that this updates Amat and our internal ell
          auto t1 = std::chrono::steady_clock::now();
          Eigen::MatrixXcd nb = nextBlock(J, jell);
          auto t2 = std::chrono::steady_clock::now();
          // computed updated SVD and note it gets placed in u2,s2,v2
          updateSVD(nb);
          auto t3 = std::chrono::steady_clock::now();
          // manually update u1,s1,v1
          std::chrono::duration<double> e1 = t2 - t1;
          std::chrono::duration<double> e2 = t3 - t2;
          /*
          std::cout << "nextBlock: " << e1.count() << " seconds\n";
          std::cout << "updateSVD: " << e2.count() << " seconds\n";   
          */
          u1 = u2;
          s1 = s2;
          v1 = v2;
        }
        if (needPropagation)
        {
          Eigen::MatrixXcd qhist = pred1rdms[i].block(0, J-jell, drc2, jell+1).rowwise().reverse();
          Eigen::VectorXcd PreconVec = pinvvec(qhist.reshaped(), this->tol);
          Eigen::MatrixXcd Precon = PreconVec.reshaped(N, N);
          pred1rdms[i].col(J+1) = BmatR * (fprops[J].adjoint() * Precon * fprops[J].transpose()).reshaped();
        }
      }
    }
  }
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
  ("delay", "Start,step,numdelays of delay range", cxxopts::value<std::vector<int>>())
  ("infile", "Input file path", cxxopts::value<std::string>())
  ("outfile", "Output file path (for MAEs)", cxxopts::value<std::string>())
  ("tol", "SVD tolerance", cxxopts::value<double>()->default_value("1e-6"));
  
  auto result = options.parse(argc, argv);
  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0;
  }
  double dt, T, tol;
  std::vector<double> timeparams, fieldparams;
  std::vector<int> delayparams;
  std::string inpath, outpath;
  bool verbose, savemae, saveqprop;
  
  // Required parameters
  getRequired(result, "time",    timeparams,  "Must specify time-stepping parameters dt,T!");
  getRequired(result, "field",   fieldparams, "Must specify field parameters freq,amp,ncyc!");
  getRequired(result, "delay",   delayparams, "Must specify delay parameters start,step,numdelays!");
  getRequired(result, "infile",  inpath,      "Must specify input file!");
  getRequired(result, "outfile", outpath,     "Must specify output file!");
  
  // Optional parameter
  getOptional(result, "tol", tol, 1e-6);
  
  // change this to be mindelay, delaystep, numdelays
  // pass the raw parameters into the constructor, in keeping with the other (time & field) params
  int ncyc = static_cast<int>(std::round(fieldparams[2]));
  memoryModel mm(timeparams[0], timeparams[1], fieldparams[0], fieldparams[1], ncyc, delayparams[0], delayparams[1], delayparams[2], tol, num_threads, inpath, outpath);
  Eigen::VectorXcd ic(mm.getdrcCI());
  ic.setZero();
  ic[0] = 1.0;
  mm.tdseProp(ic);
  mm.exact1RDMS();
  mm.filterIndices();
  
  auto start = std::chrono::steady_clock::now();
  mm.qpropALLV2();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  for (int i=0; i<delayparams[2]; ++i)
  {
    int jell = delayparams[0] + i * delayparams[1];
    std::cout << "\n\ndelay = " << jell << "\n";
    double maeTruth2 = (mm.getTrue1rdms() - mm.getPred1rdms(i)).array().abs().mean();
    std::cout << "MAE( truth - batch predictions(i) ) = " << maeTruth2 << "\n";    
  }
  
  start = std::chrono::steady_clock::now();
  mm.qpropALL();
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  for (int i=0; i<delayparams[2]; ++i)
  {
    int jell = delayparams[0] + i * delayparams[1];
    std::cout << "\n\ndelay = " << jell << "\n";
    double maeTruth2 = (mm.getTrue1rdms() - mm.getPred1rdms(i)).array().abs().mean();
    std::cout << "MAE( truth - batch predictions(i) ) = " << maeTruth2 << "\n";    
  }
  return 0;
}

