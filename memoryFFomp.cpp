// memory model
// field-free case

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#include "cnpy.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

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

// 1RDM propagator with memory length n (time steps) and time step h
Eigen::MatrixXcd qprop(int n, double h, const Eigen::VectorXd& hamiltonian, const std::vector<int>& redCols, const Eigen::MatrixXd& BmatT, double tol=1e-9)
{
  int drc2 = BmatT.rows();

  Eigen::MatrixXcd bigmat((n+1)*drc2, redCols.size());
  for (int j=0; j<=n; ++j)
    bigmat.block(j*drc2, 0, drc2, redCols.size()) = BmatT * redprop(-j, h, hamiltonian, redCols);
    
  Eigen::MatrixXcd bigmatpinv = pseudoInverse(bigmat, tol);
  Eigen::MatrixXcd thisredprop = redprop(1, h, hamiltonian, redCols);
  return BmatT * thisredprop * bigmatpinv;
  /*
  // Solve bigmat * x = (BmatT * thisredprop)^T
  // for x, in least squares sense
  auto rhs = (BmatT * thisredprop).transpose();
  auto x = bigmat.householderQr().solve(rhs);
  return x.transpose();
  */
}

int main(void)
{
  // parameters to convert to command line arguments
  double dt = 0.001;
  int delay = 2000;

  // Load the entire .npz file into a map-like structure
  cnpy::npz_t my_npz = cnpy::npz_load("fci_heh+_6-31g.npz");

  // Load ham
  cnpy::NpyArray arr = my_npz["ham"];
  double* hamdata = arr.data<double>();
  size_t length = arr.shape[0];
  Eigen::Map<Eigen::VectorXd> ham(hamdata, length);
  int drcCI = (int) length;
  int drcCI2 = drcCI*drcCI;
  std::cout << "drcCI = " << drcCI << "\n";

  // Load Bten  
  arr = my_npz["Bten"];
  double* Btendata = arr.data<double>();
  length = arr.shape[0]*arr.shape[1]*arr.shape[2]*arr.shape[3];
  int drc = (int) arr.shape[2];
  int drc2 = drc*drc;
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
  std::cout << "number of good cols = " << goodCols.size() << "\n";

  // create new matrix with selected columns
  Eigen::MatrixXd BmatTgood(drc2, goodCols.size());
  for (size_t j=0; j<goodCols.size(); ++j)
    BmatTgood.col(j) = BmatT.col(goodCols[j]);

  std::cout << "BmatTgood # of rows = " << BmatTgood.rows() << "\n";
  std::cout << "BmatTgood # of cols = " << BmatTgood.cols() << "\n";

  Eigen::MatrixXcd thisqprop = qprop(delay, dt, ham, goodCols, BmatTgood, 1e-7);
  std::cout << "thisqprop # of rows = " << thisqprop.rows() << "\n";
  std::cout << "thisqprop # of cols = " << thisqprop.cols() << "\n";

  // save qprop to disk
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"); 
  std::string name = "thisqpropR.csv";
  std::ofstream fileR(name.c_str());
  fileR << thisqprop.real().matrix().format(CSVFormat);
  name = "thisqpropI.csv";
  std::ofstream fileI(name.c_str());
  fileI << thisqprop.imag().matrix().format(CSVFormat);

  // initialize coeff matrix and set initial condition
  double T = 200.0;
  int nsteps = static_cast<int>(std::ceil(T/dt));
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
  std::cout << "About to propagate 1RDMs for " << nsteps << " steps\n";
  for (int k=0; k<=delay; ++k)
    pred1rdms.col(k) = true1rdms.col(k);

  using VecC = Eigen::VectorXcd;
  using MatC = Eigen::MatrixXcd;
  std::vector<MatC> A(delay + 1);
  for (int i = 0; i <= delay; ++i)
    A[i] = thisqprop.block(0, i * drc2, drc2, drc2);

  for (int k = delay; k < nsteps; ++k)
  {
    // IMPORTANT: disable Eigen threading
    Eigen::setNbThreads(1);
    VecC x_next = bvec;
    #pragma omp parallel
    {
      VecC local = VecC::Zero();

      #pragma omp for schedule(static)
      for (int i = 0; i <= delay; ++i) 
      {
        local.noalias() += A[i] * (pred1rdms.col(k - i) - bvec);
      }

      #pragma omp critical
      {
        x_next += local;
      }
    }
    pred1rdms.col(k + 1) = x_next;
  }
  /*
  for (int k = delay; k < nsteps; ++k)
  {
    VecC x_next = bvec;
    for (int i = 0; i <= delay; ++i) 
        x_next.noalias() += A[i] * (pred1rdms.col(k - i) - bvec);

    pred1rdms.col(k + 1) = x_next;
  }
  */
  /*
  for (int k=delay; k<nsteps; ++k)
  {
    Eigen::MatrixXcd temp = pred1rdms.block(0, k-delay, drc2, delay+1).rowwise().reverse();
    Eigen::MatrixXcd temp2 = (temp.colwise() - bvec).reshaped();
    pred1rdms.col(k+1) = bvec + thisqprop * temp2;
  } 
  */
  double mae = (true1rdms - pred1rdms).array().abs().mean();
  std::cout << "Mean Absolute Error: " << mae << "\n";

  return 0;
}


