#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --partition dept.appliedmath
#SBATCH --mem=64G
#SBATCH --time=0-01:00:00
#SBATCH --job-name=fci_h2_6-31g
#SBATCH --array=0-1000:100

./memoryFF --dt 0.01 --delay ${SLURM_ARRAY_TASK_ID} --infile ./psi4data/fci_h2_6-31g.npz --savemae ./fci_h2_6-31g

