#!/bin/bash
export CUDA_VISIBLE_DEVICES="0"

methods=("fci")
mols=("heh+")
# bases=("sto-3g" "3-21g" "6-31g" "6-31gss")
bases=("6-31gss")

for method in "${methods[@]}"; do
  for mol in "${mols[@]}"; do
    for basis in "${bases[@]}"; do
      for delay in $(seq 8500 500 10000); do
        for smax in $(seq 8 2 8); do
          echo "Running: method=${method}, mol=${mol}, basis=${basis}, delay=${delay}, smax=${smax}"
          ./memoryFF \
            --smax ${smax} \
            --dt 0.01 \
            --delay ${delay} \
            --infile "../psi4data/${method}_${mol}_${basis}.npz" \
            --savemae "./newFF/" \
            --verbose
        done
      done
    done
  done
done
