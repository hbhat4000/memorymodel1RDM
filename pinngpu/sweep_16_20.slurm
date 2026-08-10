#!/bin/bash
#SBATCH -p dept.appliedmath
#SBATCH --job-name=sweep_16_20
#SBATCH --gres=gpu:1
#SBATCH --time=35:59:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hbhat@ucmerced.edu
#SBATCH --output=sweep_16_20.out
#SBATCH --error=sweep_16_20.err
# ---------------------------------------------------

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/hbhat/lib:/u/hbhat/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/12.8/lib64
export CUDA_VISIBLE_DEVICES="0"

for ((i=1;i<=40;i++)); do
    for ((j=16;j<=20;j++)); do
        for ((k=1;k<=1;k++)); do

            # delay = 100 + 50*k
            delay=$((50*k))

            # intensity: logspace(1e-4, 2e-3, 40)
            intensity=$(awk -v j="$j" 'BEGIN {
                start = 1e-4;
                end   = 2e-3;
                ratio = exp((log(end/start))/39);
                printf "%.6f", start * (ratio^(j-1));
            }')

            # freq: logspace(0.001, 0.05, 40)
            freq=$(awk -v i="$i" 'BEGIN {
                start = 1e-3;
                end   = 5e-2;
                ratio = exp((log(end/start))/39);
                printf "%.6f", start * (ratio^(i-1));
            }')

            # Construct the expected output filename
            target_file="./sweep/fci_h2_6-31g_0.050000_${delay}_${freq}_${intensity}_1.txt"

            # Check if the file exists; if it does, skip the rest of this iteration
            if [ -f "$target_file" ]; then
                echo "Detected existing file: $target_file"
                continue
            fi

            echo $delay $intensity $freq

            ./memoryFO \
                --time 0.05,1000.0 \
                --field $freq,$intensity,1 \
                --delay $delay \
                --infile ../psi4data/fci_h2_6-31g.npz \
                --outpath ./sweep6-31g/ 

        done
    done
done

