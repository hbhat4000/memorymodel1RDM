#!/bin/bash
export CUDA_VISIBLE_DEVICES="1"

methods=("cis" "fci")
mols=("heh+" "h2")
bases=("sto-3g" "3-21g" "6-31g" "6-31gss")
# bases=("6-31gss")

for method in "${methods[@]}"; do
  for mol in "${mols[@]}"; do
    for basis in "${bases[@]}"; do
      for delay in $(seq 0 1 0); do
        for smax in $(seq 2 1 2); do
          echo "Running: method=${method}, mol=${mol}, basis=${basis}, delay=${delay}, smax=${smax}"
          ./memoryFF \
            --smax ${smax} \
            --dt 0.01 \
            --delay ${delay} \
            --infile "../psi4data/${method}_${mol}_${basis}.npz" \
            --verbose
            # --savemae "./newFF/" \
        done
      done
    done
  done
done
