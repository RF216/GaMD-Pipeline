#!/bin/bash
#SBATCH --job-name=3_initial_cmd
#SBATCH --account=oz241
#SBATCH --partition=milan-gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=16G
#SBATCH --time=15:00:00
#SBATCH --output=prod.out
#SBATCH --error=prod.err

cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/11.2.0 openmpi/4.1.1
module load amber/20.13-ambertools-22.3-cuda-11.7.0

# cMD run for 50ns
chmod +x run_cmd.sh
./run_cmd.sh
