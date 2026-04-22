#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --partition bigmem
#SBATCH --mem=512G
#SBATCH --time=0-12:00:00
#SBATCH --job-name=cis_lih_6-31g
#SBATCH --array=48-49

delay=$((5100 + 100 * SLURM_ARRAY_TASK_ID))

./memoryFF --dt 0.001 --delay ${delay} --infile ./psi4data/cis_lih_6-31g.npz --savemae ./cis_lih_6-31g

