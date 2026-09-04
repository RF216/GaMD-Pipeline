#!/bin/bash
#SBATCH --job-name=2_equalization
#SBATCH --account=oz241
#SBATCH --partition=milan-gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=eq.out
#SBATCH --error=eq.err

module purge
module load gcc/11.2.0 openmpi/4.1.1
module load amber/20.13-ambertools-22.3-cuda-11.7.0
cd "$SLURM_SUBMIT_DIR"

chmod +x run_equil.sh
./run_equil.sh
