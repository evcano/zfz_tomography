#!/bin/bash
#SBATCH --job-name=ZAFRAN08
#SBATCH --partition=batch
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCG --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --hint=nomultithread


#OpenMP settings:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#run the application:
srun ./MCTomo
