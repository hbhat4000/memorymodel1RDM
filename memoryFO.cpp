// memory model
// field-on case

#define EIGEN_USE_MKL_ALL

#include <mkl.h>

// This tells the compiler to intentionally crash if MKL is NOT using 32-bit integers (LP64)
static_assert(sizeof(MKL_INT) == 4, "FATAL ERROR: MKL is using the ILP64 interface instead of LP64!");

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include <cxxopts.hpp>
#include "cnpy.h"
#include <utility>       // for std::move
#include <filesystem>
#include <string>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <omp.h>

using namespace std::complex_literals;

// Function to compute the Moore–Penrose pseudoinverse
Eigen::MatrixXcd pseudoInverse(const Eigen::MatrixXcd &mat, double tolerance)
{
    // Compute SVD: mat = U * Σ * Vᵀ
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd( 
        mat, Eigen::ComputeThinU | Eigen::ComputeThinV
    );

    const Eigen::VectorXd &singularValues = svd.singularValues();
    Eigen::VectorXd singularValuesInv(singularValues.size());

    // Invert singular values with tolerance to avoid division by zero
    for (int i = 0; i < singularValues.size(); ++i) {
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
  const int ncyc, nsteps;
  const Eigen::VectorXi delayrange;
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

  public:
    //constructor
    memoryModel(double dt, double T, double freq, double amp, int ncyc, double svdtol, 
                Eigen::VectorXi idelayrange, std::string infile, std::string ioutfile);

    int getdrc(void) { return drc; }
    int getdrcCI(void) { return drcCI; }
    int getN(void) { return N; }
    Eigen::MatrixXcd getTrue1rdms(void) { return true1rdms; }

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
    
    // propagate 1RDM with memory model for particular delay value
    Eigen::MatrixXcd qprop(int jell);
};

memoryModel::memoryModel(double dt, double T, double freq, double amp, int ncyc, double svdtol, 
                         Eigen::VectorXi idelayrange, std::string infile, std::string ioutfile)
   : h(dt), T(T), freq(freq), amp(amp), ncyc(ncyc), tol(svdtol),
     delayrange(std::move(idelayrange)), outfile(std::move(ioutfile)), nsteps(static_cast<int>(std::ceil(T/h)))
{
  offstep = static_cast<int>(std::ceil(ncyc / (dt * freq)));

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
  pred1rdms.resize(delayrange.size(), Eigen::MatrixXcd::Zero(drc2, nsteps+1));
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

  // std::cout << "coeffs.col(0) = " << coeffs.col(0) << "\n";
    
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
  mkl_set_num_threads(1);

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

  mkl_set_num_threads(1);
  #pragma omp parallel for
  for (int k=0; k<nsteps; ++k)
    fprops[k] = props[k](goodStatesVec, goodStatesVec);
  
  filtered = true;
  return 0;
}

// form bigmat at particular point in time for particular delay value
int memoryModel::bigmatBuild(int J, int jell)
{
  if (!filtered)
  {
    std::cout << "Filtered indices and propagated must be computed first!\n";
    return 1;
  }  
  // resize bigmat and store first block
  bigmat.setZero( (jell+1)*drc2, N2 );
  bigmat.block(0, 0, drc2, N2) = BmatR;

  // this will hold our propagator chain
  Eigen::MatrixXcd Amat = Eigen::MatrixXcd::Identity( N, N );  

  mkl_set_num_threads(1);
  
  for (int j=1; j<=jell; ++j)
  {
    Amat = Amat * fprops[J - j];
    MatrixXcdRowMajor Result(drc2, N2);
    
    // Parallelize the batch processing over M^2 matrices
    #pragma omp parallel for
    for (int k=0; k<drc2; ++k)
    {
      // Map the contiguous memory of the k-th row to an N x N matrix
      Eigen::Map<MatrixXcdRowMajor> BmatRk(BmatR.row(k).data(), N, N);
      Eigen::Map<MatrixXcdRowMajor> Rk(Result.row(k).data(), N, N);
      
      // Evaluate in two steps to guarantee two clean BLAS GEMM calls per batch
      Eigen::MatrixXcd temp = Amat.conjugate() * BmatRk;
      Rk.noalias() = temp * Amat.transpose(); 
    }
    bigmat.block(j*drc2, 0, drc2, N2) = Result;
  }
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

// simple function that propagates 1RDM forward in time 
// - take pseudoinverse of freshly constructed bigmat at each step
// - don't worry about parallelization or incremental updates or anything
Eigen::MatrixXcd memoryModel::qprop(int jell)
{
  // declare and initialize predicted 1RDMS
  Eigen::MatrixXcd preds(drc2, nsteps+1);
  preds.setZero();

  // copy in some ground truth 1RDMs
  for (int k=0; k<=jell; ++k)
    preds.col(k) = true1rdms.col(k);
  
  if (!filtered)
  {
    std::cout << "Filtered indices and propagated must be computed first!\n";
    return preds;
  }
  
  // then loop over history and build bigmat and use pseudoinverse to propagate...
  for (int k=jell; k<nsteps; ++k)
  {
    Eigen::MatrixXcd temp = preds.block(0, k-jell, drc2, jell+1).rowwise().reverse();
    bigmatBuild(k, jell);
    Eigen::MatrixXcd bigmatpinv = pseudoInverse(bigmat, tol);
    Eigen::VectorXcd PreconVec = bigmatpinv * temp.reshaped();
    Eigen::MatrixXcd Precon = PreconVec.reshaped(N, N);
    // std::cout << "Precon = " << std::setprecision(15) << Precon << "\n";
    // Eigen::VectorXcd PtrueVec = (coeffs.col(k) * coeffs.col(k).adjoint()).transpose().reshaped();
    // Eigen::MatrixXcd Ptrue = PtrueVec(goodCols).reshaped(N, N);
    // std::cout << "Ptrue = " << std::setprecision(15) << Ptrue << "\n";    
    // std::cout << "|| Precon - Ptrue || = " << std::setprecision(15) << (Precon - Ptrue).norm() << "\n";
    preds.col(k+1) = BmatR * (fprops[k].adjoint() * Precon * fprops[k].transpose()).reshaped();
    // std::cout << "preds.col(k+1) = " << std::setprecision(15) << preds.col(k+1) << "\n";
  }
  return preds;
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
  cxxopts::Options options("memoryFO", "Field-on memory model for 1RDM propagation");
  options.add_options()
  ("h,help", "Print usage")
  ("time", "Time-stepping parameters dt,T", cxxopts::value<std::vector<double>>())
  ("field", "Freq,amp,ncyc of applied electric field in z direction", cxxopts::value<std::vector<double>>())
  ("delay", "Start,end,step of delay range", cxxopts::value<std::vector<int>>())
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
  getRequired(result, "time",    timeparams,  "Must specify time-stepping parameters!");
  getRequired(result, "field",   fieldparams, "Must specify field parameters!");
  getRequired(result, "delay",   delayparams, "Must specify delay parameters!");
  getRequired(result, "infile",  inpath,      "Must specify input file!");
  getRequired(result, "outfile", outpath,     "Must specify output file!");
  
  // Optional parameter
  getOptional(result, "tol", tol, 1e-6);
  
  int n = (delayparams[1] - delayparams[0]) / delayparams[2] + 1;
  Eigen::VectorXi dr = Eigen::ArrayXi::LinSpaced(n, delayparams[0], delayparams[1]);
  int ncyc = static_cast<int>(std::round(fieldparams[2]));
  memoryModel mm(timeparams[0], timeparams[1], fieldparams[0], fieldparams[1], ncyc, tol, dr, inpath, outpath);
  Eigen::VectorXcd ic(mm.getdrcCI());
  ic.setZero();
  ic[0] = 1.0;
  mm.tdseProp(ic);
  mm.exact1RDMS();
  mm.filterIndices();
  // mm.bigmatBuild(100,5);
  // mm.bigmatPrint();
  Eigen::MatrixXcd preds = mm.qprop(100);
  Eigen::MatrixXcd truth = mm.getTrue1rdms();
  double mae = (truth - preds).array().abs().mean();
  std::cout << "Mean Absolute Error: " << mae << "\n";
  return 0;
}


