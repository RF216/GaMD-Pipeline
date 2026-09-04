#!/bin/bash
#SBATCH --job-name=5_gamd_prod
#SBATCH --account=oz241
#SBATCH --partition=milan-gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=8G
#SBATCH --time=06:00:00

set -u
module purge
module load gcc/11.2.0 openmpi/4.1.1
module load amber/20.13-ambertools-22.3-cuda-11.7.0
cd "$SLURM_SUBMIT_DIR"

chmod +x run_prod.sh
REP=${REP:-1}; SEG=${SEG:-1}; NSEG=${NSEG:-25}
DIR="rep${REP}"; mkdir -p "$DIR" logs
PRM="$(realpath ../1_prepare_struct/system.prmtop)"
INFILE="$(realpath prod.in)"

if [ "$SEG" -eq 1 ]; then
    PREV="$(realpath ../4_gamd_prep/system_gamd_prep.rst7)"
    GAMD_RESTART_SRC="$(realpath ../4_gamd_prep/gamd-restart.dat)"
    [ -f "$GAMD_RESTART_SRC" ] || { echo "missing gamd-restart.dat from prep stage"; exit 1; }
    if [ ! -f "${DIR}/gamd-restart.dat" ]; then
        cp "$GAMD_RESTART_SRC" "${DIR}/gamd-restart.dat"
    fi
else
    PREV="$(realpath ${DIR}/seg$((SEG-1)).rst7)"
fi
[ -f "$PREV" ] || { echo "missing $PREV"; exit 1; }

if [ ! -f "${DIR}/seg${SEG}.done" ]; then
    ( cd "$DIR" && pmemd.cuda -O -i "$INFILE" -o "seg${SEG}.out" -p "$PRM" -c "$PREV" \
               -r "seg${SEG}.rst7" -x "seg${SEG}.nc" \
               -inf "seg${SEG}.mdinfo" )
    if [ $? -ne 0 ] || [ ! -f "${DIR}/seg${SEG}.rst7" ]; then
        echo "*** rep $REP seg $SEG FAILED - chain stopped"; tail -40 "${DIR}/seg${SEG}.out"; exit 1
    fi
    mv "${DIR}/gamd.log" "${DIR}/seg${SEG}.gamd.log"
    touch "${DIR}/seg${SEG}.done"
    grep -E "ns/day" "${DIR}/seg${SEG}.out" | tail -1
fi

if [ "$SEG" -lt "$NSEG" ]; then
    NEXTSEG=$((SEG+1))
    sbatch --job-name=${REP}-5_gamd_prod \
           --output="logs/rep${REP}_seg${NEXTSEG}.out" \
           --error="logs/rep${REP}_seg${NEXTSEG}.err" \
           --export=REP=${REP},SEG=${NEXTSEG},NSEG=${NSEG} run_prod.slurm
else
    echo "REPLICATE $REP COMPLETE - $((NSEG*20)) ns"
fi
