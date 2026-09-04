#!/bin/bash
#SBATCH --job-name=4_gamd_equil
#SBATCH --account=oz241
#SBATCH --partition=milan-gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=gamd_equil.log
#SBATCH --error=gamd_equil.err

cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/11.2.0 openmpi/4.1.1
module load amber/20.13-ambertools-22.3-cuda-11.7.0

# gamd equilibration
chmod +x run_gamd_prep.sh
./run_gamd_prep.sh
