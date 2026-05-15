// memory model
// field-on case
// redo the code to cache products
// also, no more incremental SVD

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
using MatrixXcdRowMajor = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class memoryModel
{
  const double h, T, tol, freq, amp;
  const int ncyc, nsteps, delaystart, delaystep, numdelays, numthreads;
  const std::string inpath, outpath;
  int drc, drc2, drcCI, drcCI2;
  int ell;                                 // state variable for grab methods
  int ellmax;                              // maximum delay/memory considered
  int N, N2;                               // N = number of good states
  int offstep;                             // time step at which the field is switched off
  bool havecoeffs = false;                 // have we computed time-dependent coefficients yet?
  bool have1rdms = false;                  // have we computed 1-RDMS yet?
  bool filtered = false;                   // have we filtered out bad indices yet?
  bool builtpcc = false;                   // have we built the propagator chain cache yet?
  bool builtbmc = false;                   // have we built the bigmat cache yet?
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
  Eigen::MatrixXcd kronprop;               // one-step field-free Kronecker product propagator
  // propagator chain cache (pcc)
  // the index for the std::vector is J, a discrete time step
  // each complex matrix in the cache is of size (ell*N) x N
  std::vector<Eigen::MatrixXcd> pcc;
  // bigmat cache (bmc)
  // the index for the std::vector is J, a discrete time step
  // after an initial segment of smaller bigmats,
  // each complex matrix in the cache is of size ((ell+1)*drc2) x N2
  std::vector<Eigen::MatrixXcd> bmc;
 
  public:
    //constructor
    memoryModel(double dt, double T, double freq, double amp, int ncyc, int delaystart, int delaystep, int numdelays, double svdtol, int numthreads, std::string infile, std::string outpath);

    int getdrc(void) { return drc; }
    int getdrc2(void) { return drc2; }
    int getdrcCI(void) { return drcCI; }
    int getN(void) { return N; }
    int getnsteps(void) { return nsteps; }
    Eigen::MatrixXcd getTrue1rdms(void) { return true1rdms; }
    Eigen::MatrixXcd getPred1rdms(int j) { return pred1rdms[j]; }

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
    int bigmatFromCache(const int J, const int jell, Eigen::MatrixXcd& bigmat);
    int bigmatBuildLocal(const int J, const int jell, Eigen::MatrixXcd& bigmat);

    // build the bigmat cache
    int buildBMC(void);

    // print bigmat to screen
    int bigmatPrint(Eigen::MatrixXcd bigmat);

    // propagate 1RDM with memory model for all delay values at once
    int qpropALLV2(void);
    int saveResults(void);
};

memoryModel::memoryModel(double dt, double T, double freq, double amp, int ncyc, int delaystart, int delaystep, int numdelays, double svdtol, int numthreads, std::string infile, std::string outpath)
   : h(dt), T(T), freq(freq), amp(amp), ncyc(ncyc), 
     delaystart(delaystart), delaystep(delaystep), numdelays(numdelays), 
     tol(svdtol), numthreads(numthreads), 
     inpath(std::move(infile)), outpath(std::move(outpath)), 
     nsteps(static_cast<int>(std::ceil(T/h)))
{
  offstep = static_cast<int>(std::ceil(ncyc / (dt * freq)));
  std::cout << "Field will be on for " << offstep << " time steps\n";
  // Load the entire .npz file into a map-like structure
  cnpy::npz_t my_npz = cnpy::npz_load(inpath);

  // max delay/memory length that we consider
  ellmax = delaystart + numdelays*delaystep;

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

int memoryModel::tdseProp(const Eigen::VectorXcd& ic)
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
  int upper = ellmax + offstep;
  int pccsize = (nsteps < upper) ? nsteps : upper;
  pcc.resize(pccsize + 1);
  for (int J=ellmax; J<=(ellmax+offstep); ++J)
  {
    if (J==nsteps) break;
    pcc[J].setZero(ellmax*N, N);
    pcc[J].block(0, 0, N, N) = fprops[J-1];
    if (J==ellmax) // construct a propagator chain that goes back ellmax steps
    {
      for (int j=2; j<=ellmax; ++j)
        pcc[J].block((j-1)*N, 0, N, N) = pcc[J].block((j-2)*N, 0, N, N) * fprops[J-j];
    }
    else 
    {
      #pragma omp parallel for schedule(dynamic)
      for (int j=2; j<=ellmax; ++j)
      {
        pcc[J].block((j-1)*N, 0, N, N).noalias() = fprops[J-1] * pcc[J-1].block((j-2)*N, 0, N, N);
      }
    }
  }
  builtpcc = true;
  std::cout << "Built propagator chain cache!\n";
  return 0;
}

int memoryModel::bigmatFromCache(const int J, const int jell, Eigen::MatrixXcd& bigmat)
{
  Eigen::MatrixXcd thisblock;
  if (!builtpcc)
  {
    std::cout << "Propagator chain cache must be computed first!\n";
    return 1;
  }
  bigmat.block(0, 0, drc2, N2) = BmatR;
  for (int j=1; j<=jell; ++j)
  {
    thisblock = pcc[J].block((j-1)*N,0,N,N);
    bigmat.block(j*drc2, 0, drc2, N2).noalias() = BmatMultNP(thisblock.conjugate(), thisblock.transpose());
  }
  return 0;
}

// form bigmat at particular point in time for particular delay value
int memoryModel::bigmatBuildLocal(const int J, const int jell, Eigen::MatrixXcd& bigmat)
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
  for (int j=1; j<=jell; ++j)
  {
    Amat = Amat * fprops[J - j];
    bigmat.block(j*drc2, 0, drc2, N2) = BmatMultNP(Amat.conjugate(), Amat.transpose());
  }
  return 0;
}

int memoryModel::buildBMC(void)
{
  if (!builtpcc)
  {
    std::cout << "Propagator chain cache must be computed first!\n";
    return 1;
  }
  int bmcsize = pcc.size();
  bmc.resize(bmcsize);
  // these first bigmats are weird, because we don't have enough time steps
  // to go all the way back ellmax steps in time!
  // so, each matrix here will be of a different size...
  // ...as we eventually build up to the "actual" bigmats
  #pragma omp parallel for schedule(dynamic)
  for (int J=delaystart; J<ellmax; ++J)
  {
    bmc[J].resize((J+1)*drc2, N2);
    bigmatBuildLocal(J, J, bmc[J]);
  }
  // now we have reached the point where our pcc holds actual stuff
  // so let us use these cached propagator chains to build bigmat
  // note that loop start and end match the main pcc loop
  #pragma omp parallel for schedule(dynamic)
  for (int J=ellmax; J<bmcsize; ++J)
  {
    bmc[J].resize((ellmax+1)*drc2, N2);
    bigmatFromCache(J, ellmax, bmc[J]);
  }
  builtbmc = true;
  std::cout << "Built bigmat cache!\n";
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

// propagate 1RDM with memory model for all delay values at once
int memoryModel::qpropALLV2(void)
{
  if (!builtbmc)
  {
    std::cout << "Bigmat cache must be computed first!\n";
    return 1;
  }
  std::vector<Eigen::MatrixXcd> rdmprops;
  rdmprops.resize(numdelays);

  for (int i=0; i<numdelays; ++i)
  {
    int jell = delaystart + i*delaystep;
    // copy in some ground truth 1RDMs
    for (int k=0; k<=jell; ++k)
      pred1rdms[i].col(k) = true1rdms.col(k);    
  }
  for (int J=delaystart; J<nsteps; ++J)
  {
    // if ((J % 100) == 0)
    std::cout << "Propagating 1RDM at time step " << J << "\n";
    // main loop, try parallelizing over delays
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<numdelays; ++i)
    {
      int jell = delaystart + i*delaystep;
      // catch the case where we cannot propagate the memory model because we don't have enough memory yet!
      if ((J-jell) < 0) pred1rdms[i].col(J+1) = true1rdms.col(J+1);
      else if ((J-jell) >= offstep)
      {
        if ((J-jell) == offstep)
        {
          Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(bmc[J].block(0,0,(jell+1)*drc2,N2));
          // compute pseudoinverse locally
          Eigen::VectorXd s_inv(svd.singularValues().size());
          for (int k = 0; k < svd.singularValues().size(); ++k)
            s_inv(k) = (svd.singularValues()(k) > this->tol) ? 1.0 / svd.singularValues()(k) : 0.0;

          rdmprops[i] = BmatR * kronprop * svd.matrixV() * s_inv.asDiagonal() * svd.matrixU().adjoint();
        }
        Eigen::MatrixXcd qhist = pred1rdms[i].block(0, J-jell, drc2, jell+1).rowwise().reverse();
        pred1rdms[i].col(J+1) = rdmprops[i] * qhist.reshaped();
      }
      else
      {
        // grab bigmat from the cache, take its SVD
        Eigen::BDCSVD<Eigen::MatrixXcd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(bmc[J].block(0,0,(jell+1)*drc2,N2));
        // compute pseudoinverse locally
        Eigen::VectorXd s_inv(svd.singularValues().size());
        for (int k = 0; k < svd.singularValues().size(); ++k)
          s_inv(k) = (svd.singularValues()(k) > this->tol) ? 1.0 / svd.singularValues()(k) : 0.0;
        
        Eigen::MatrixXcd qhist = pred1rdms[i].block(0, J-jell, drc2, jell + 1).rowwise().reverse();
        Eigen::VectorXcd temp = svd.matrixU().adjoint() * qhist.reshaped();
        Eigen::VectorXcd PreconVec = svd.matrixV() * (s_inv.asDiagonal() * temp);
        Eigen::MatrixXcd Precon = PreconVec.reshaped(N, N);
        
        // propagate one step forward in time
        pred1rdms[i].col(J + 1) = BmatR * (fprops[J].adjoint() * Precon * fprops[J].transpose()).reshaped();          
      }
    }
  }
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
  out << "delaystart = " << delaystart << ", delaystep = " << delaystep << ", numdelays = " << numdelays << "\n";
  // save all the MAEs
  for (int i=0; i<numdelays; ++i)
  {
    int iell = delaystart + i * delaystep;
    double maeTruth = (getTrue1rdms() - getPred1rdms(i)).array().abs().mean();
    std::ostringstream oss;
    oss << std::setprecision(17) << maeTruth;
    out << iell << ", " << oss.str() << "\n";
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
  ("outpath", "Output file path (for MAEs)", cxxopts::value<std::string>())
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
  getRequired(result, "outpath", outpath,     "Must specify output path!");
  
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
  mm.buildPCC();
  mm.buildBMC();
  mm.qpropALLV2();
  mm.saveResults();
  return 0;
}

