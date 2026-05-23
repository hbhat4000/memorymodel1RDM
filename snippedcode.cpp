
    // std::cout << "Precon = " << std::setprecision(15) << Precon << "\n";
    // Eigen::VectorXcd PtrueVec = (coeffs.col(k) * coeffs.col(k).adjoint()).transpose().reshaped();
    // Eigen::MatrixXcd Ptrue = PtrueVec(goodCols).reshaped(N, N);
    // std::cout << "Ptrue = " << std::setprecision(15) << Ptrue << "\n";    
    // std::cout << "|| Precon - Ptrue || = " << std::setprecision(15) << (Precon - Ptrue).norm() << "\n";
    // std::cout << "preds.col(k+1) = " << std::setprecision(15) << preds.col(k+1) << "\n";

  // mm.bigmatBuild(100,5);
  // mm.bigmatPrint();

/*
  auto start = std::chrono::steady_clock::now();
  std::vector<Eigen::MatrixXcd> preds;
  preds.resize(delayparams[2], Eigen::MatrixXcd::Zero(mm.getdrc2(), mm.getnsteps()+1));
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i<delayparams[2]; ++i)
  {
    int jell = delayparams[0] + i * delayparams[1];
    std::cout << "delay = " << jell << "\n";
    preds[i] = mm.qprop(jell);
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
*/

  /*
  for (int i=0; i<delayparams[2]; ++i)
  {
    int jell = delayparams[0] + i * delayparams[1];
    std::cout << "\n\ndelay = " << jell << "\n";
    double maePreds = (mm.getPred1rdms(i) - preds[i]).array().abs().mean();
    std::cout << "MAE( batch predictions(i) - single prediction at i ) = " << maePreds << "\n";
    double maeTruth1 = (mm.getTrue1rdms() - preds[i]).array().abs().mean();
    std::cout << "MAE( truth - single prediction at i ) = " << maeTruth1 << "\n";
    double maeTruth2 = (mm.getTrue1rdms() - mm.getPred1rdms(i)).array().abs().mean();
    std::cout << "MAE( truth - batch predictions(i) ) = " << maeTruth2 << "\n";    
  }
  */
