    if (k == delay) 
    {
      // 1. CPU Ground Truth
      Eigen::MatrixXcd temp = pred1rdms.block(0, k-delay, drc2, delay+1).rowwise().reverse();
      Eigen::VectorXcd temp2 = (temp.colwise() - bvec).reshaped();
          
      // 2. GPU Kernel Execution
      cudaMemset(d_temp2, 0, temp2_size * sizeof(cuDoubleComplex));
      prep_temp2_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_pred1rdms, d_bvec, d_temp2, k, delay, drc2);
      cudaDeviceSynchronize(); // Force wait
          
      // 3. Bring GPU result back
      Eigen::VectorXcd ctemp2(temp2_size);
      cudaMemcpy(ctemp2.data(), d_temp2, temp2_size * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
          
      // 4. Compare and exit
      double max_diff = (temp2 - ctemp2).cwiseAbs().maxCoeff();
      std::cout << "Kernel Debug | Max element-wise difference: " << max_diff << "\n";
          
      if (max_diff > 1e-12) {
        std::cout << "KERNEL FAILING: Index mapping is incorrect.\n";
      } else {
        std::cout << "KERNEL PASSING: Math is perfectly identical to Eigen!\n";
      }
        exit(0); 
    }
    cudaMemset(d_temp2, 0, temp2_size * sizeof(cuDoubleComplex));
