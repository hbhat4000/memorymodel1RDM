#!/bin/bash
#SBATCH --job-name=sweep
#SBATCH --time=00-24:00
#SBATCH --partition=gpuA100x8
#SBATCH --mem=64G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --account=bdgj-delta-gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hbhat@ucmerced.edu
#SBATCH --output=sweep.job.%j.out
#SBATCH --error=sweep.job.%j.err
# ---------------------------------------------------

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/hbhat/lib:/u/hbhat/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/12.8/lib64

for ((i=1;i<=10;i++)); do
    for ((j=1;j<=10;j++)); do
        for ((k=0;k<=2;k++)); do

            # delay = 100 + 50*k
            delay=$((100 + 50*k))

            # intensity: logspace(1e-4, 1e-2, 10)
            intensity=$(awk -v j="$j" 'BEGIN {
                start = 1e-4;
                end   = 1e-2;
                ratio = exp((log(end/start))/9);
                printf "%.6f", start * (ratio^(j-1));
            }')

            # freq: logspace(0.001, 0.02, 10)
            freq=$(awk -v i="$i" 'BEGIN {
                start = 1e-3;
                end   = 2e-2;
                ratio = exp((log(end/start))/9);
                printf "%.6f", start * (ratio^(i-1));
            }')

            # Construct the expected output filename
            target_file="./sweep/fci_h2_6-31gss_0.050000_${delay}_${freq}_${intensity}_1.txt"

            # Check if the file exists; if it does, skip the rest of this iteration
            if [ -f "$target_file" ]; then
                echo "Detected existing file: $target_file"
                continue
            fi

            echo $delay $intensity $freq

            ./memoryFO \
                --time 0.05,200.0 \
                --field $freq,$intensity,1 \
                --delay $delay \
                --infile ../psi4data/fci_h2_6-31gss.npz \
                --outpath ./sweep/ \
                --verbose --savetraj

        done
    done
done

