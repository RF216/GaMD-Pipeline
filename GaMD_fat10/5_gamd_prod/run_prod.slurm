#!/bin/bash
#SBATCH --job-name=FAT10-prod
#SBATCH --account=oz241
#SBATCH --partition=milan-gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=8G
#SBATCH --time=02:30:00
#SBATCH --output=logs/rep_seg_%j.out
#SBATCH --error=logs/rep_seg_%j.err
 
set -u
module purge
module load gcc/11.2.0 openmpi/4.1.1
module load amber/20.13-ambertools-22.3-cuda-11.7.0
cd "$SLURM_SUBMIT_DIR"
 
REP=${REP:-1}; SEG=${SEG:-1}; NSEG=${NSEG:-25}
DIR="rep${REP}"; mkdir -p "$DIR" logs
PRM=../1_prepare_struct/fat10hmr.prmtop
 
INFILE=prod.in
if [ "$SEG" -eq 1 ]; then
    PREV=../4_gamd_prep/fat10_md-1.rst7
else
    PREV="${DIR}/seg$((SEG-1)).rst7"
fi
[ -f "$PREV" ] || { echo "missing $PREV"; exit 1; }
 
if [ ! -f "${DIR}/seg${SEG}.done" ]; then
    pmemd.cuda -O -i "$INFILE" -o "${DIR}/seg${SEG}.out" -p "$PRM" -c "$PREV" \
               -r "${DIR}/seg${SEG}.rst7" -x "${DIR}/seg${SEG}.nc" \
               -inf "${DIR}/seg${SEG}.mdinfo"
    if [ $? -ne 0 ] || [ ! -f "${DIR}/seg${SEG}.rst7" ]; then
        echo "*** rep $REP seg $SEG FAILED - chain stopped"; tail -40 "${DIR}/seg${SEG}.out"; exit 1
    fi
    touch "${DIR}/seg${SEG}.done"
    grep -E "ns/day" "${DIR}/seg${SEG}.out" | tail -1
fi
 
if [ "$SEG" -lt "$NSEG" ]; then
    sbatch --export=REP=${REP},SEG=$((SEG+1)),NSEG=${NSEG} run_prod.slurm
else
    echo "REPLICATE $REP COMPLETE - $((NSEG*20)) ns"
fi
