#!/bin/bash
# Run the entire GaMD Pipeline for a user specified molecule

# Step 0: Audit the input PDB file for missing atoms and residues
python3 0_audit_struct/audit_pdb.py $1

# Step 1: Prepare the system structure and topology files
chmod +x 1_prepare_struct/run_step_1.sh
./1_prepare_struct/run_step_1.sh $1

# Step 2: Run the minimisation, heating and equilibration steps
sbatch 2_system_prep/run_equil.slurm

# Step 3: Run the initial CMD simulation
sbatch 3_initial_cmd/run_initial_cmd.slurm

# Step 4: Run the GaMD preparation step
sbatch 4_gamd_prep/run_gamd_prep.slurm

# Step 5: Run the GaMD production simulations
sbatch --job-name=1-5_gamd_prod \
       --output="logs/rep1_seg1.out" \
       --error="logs/rep1_seg1.err" \
       --export=REP=1,SEG=1,NSEG=25 5_gamd_prod/run_prod.slurm

sbatch --job-name=2-5_gamd_prod \
       --output="logs/rep1_seg1.out" \
       --error="logs/rep1_seg1.err" \
       --export=REP=1,SEG=1,NSEG=25 5_gamd_prod/run_prod.slurm

sbatch --job-name=3-5_gamd_prod \
       --output="logs/rep1_seg1.out" \
       --error="logs/rep1_seg1.err" \
       --export=REP=1,SEG=1,NSEG=25 5_gamd_prod/run_prod.slurm
       