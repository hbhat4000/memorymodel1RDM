#!/bin/bash
#SBATCH --job-name=sweep
#SBATCH --time=00:45:00             # Reduced time! Each task only runs ONE iteration now.
#SBATCH --partition=gpuA100x8
#SBATCH --mem=64G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --account=bdgj-delta-gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hbhat@ucmerced.edu
#SBATCH --output=sweep.job.%A_%a.out # %A is the job array ID, %a is the task ID
#SBATCH --error=sweep.job.%A_%a.err
#SBATCH --array=119-299              # 10 * 10 * 3 = 300 total tasks
# ---------------------------------------------------

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/hbhat/lib:/u/hbhat/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/12.8/lib64

# Map the single SLURM_ARRAY_TASK_ID (0-299) back to the 3D loop variables
TASK_ID=$SLURM_ARRAY_TASK_ID

# k originally ranged 0 to 2 (3 elements)
k=$(( TASK_ID % 3 ))

# j originally ranged 1 to 10 (10 elements)
j=$(( (TASK_ID / 3) % 10 + 1 ))

# i originally ranged 1 to 10 (10 elements)
i=$(( (TASK_ID / 30) % 10 + 1 ))

# ---------------------------------------------------
# Variable Setup
# ---------------------------------------------------

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

# Check if the file exists; if it does, exit this specific array task cleanly
if [ -f "$target_file" ]; then
    echo "Detected existing file: $target_file. Skipping execution."
    exit 0
fi

# ---------------------------------------------------
# Execution
# ---------------------------------------------------

echo "Running Array Task ID: $TASK_ID (Mapped to i=$i, j=$j, k=$k)"
echo "Delay: $delay | Intensity: $intensity | Freq: $freq"

./memoryFO \
    --time 0.05,200.0 \
    --field $freq,$intensity,1 \
    --delay $delay \
    --infile ../psi4data/fci_h2_6-31gss.npz \
    --outpath ./sweep/ \
    --verbose --savetraj
