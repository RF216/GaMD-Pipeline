#!/bin/bash
# Runs the full 5-stage minimisation/heating/equilibration pipeline.
# Invoked by equalization.pbs (which handles module load + cd $PBS_O_WORKDIR).

PRM=../1_prepare_struct/system.prmtop
INP=../1_prepare_struct/system_raw.inpcrd

# ---------------------------------------------------------------------------
# Write all mdin files inline so the job is fully self-contained and does
# not depend on separate files under md/.
# ---------------------------------------------------------------------------
cat > min1.in << 'EOF'
Stage 1: minimisation, solute heavy atoms restrained
&cntrl
 imin=1, maxcyc=5000, ncyc=2500,
 ntb=1, cut=10.0,
 ntr=1, restraint_wt=10.0, restraintmask='!:WAT,Na+,Cl- & !@H=',
 ntpr=250, ntxo=2,
/
EOF

cat > min2.in << 'EOF'
Stage 2: unrestrained minimisation (isopeptide C-N relaxes 1.466 -> ~1.335 A here)
&cntrl
 imin=1, maxcyc=10000, ncyc=5000,
 ntb=1, cut=10.0, ntr=0,
 ntpr=250, ntxo=2,
/
EOF

cat > heat.in << 'EOF'
Stage 3: heating 0 -> 310 K, NVT, 100 ps, restrained
&cntrl
 imin=0, irest=0, ntx=1,
 nstlim=50000, dt=0.002,
 ntc=2, ntf=2, ntb=1, ntp=0, cut=10.0,
 ntt=3, gamma_ln=2.0, ig=-1, tempi=0.0, temp0=310.0,
 ntr=1, restraint_wt=10.0, restraintmask='!:WAT,Na+,Cl- & !@H=',
 nmropt=1,
 ntpr=1000, ntwx=5000, ntwr=10000, ioutfm=1, ntxo=2, iwrap=1,
/
&wt type='TEMP0', istep1=0,     istep2=40000, value1=0.0,   value2=310.0 /
&wt type='TEMP0', istep1=40001, istep2=50000, value1=310.0, value2=310.0 /
&wt type='END' /
EOF

cat > density.in << 'EOF'
Stage 4: density equilibration, NPT, 100 ps, UNRESTRAINED (see the ntr bug)
&cntrl
 imin=0, irest=1, ntx=5,
 nstlim=25000, dt=0.004,
 ntc=2, ntf=2, ntb=2, ntp=1, barostat=2, pres0=1.0, cut=10.0,
 ntt=3, gamma_ln=2.0, ig=-1, temp0=310.0, ntr=0,
 ntpr=1000, ntwx=5000, ntwr=10000, ioutfm=1, ntxo=2, iwrap=1,
/
EOF

cat > equil.in << 'EOF'
Stage 5: free equilibration, NPT, 5 ns, dt = 4 fs
&cntrl
 imin=0, irest=1, ntx=5,
 nstlim=1250000, dt=0.004,
 ntc=2, ntf=2, ntb=2, ntp=1, barostat=2, pres0=1.0, cut=10.0,
 ntt=3, gamma_ln=2.0, ig=-1, temp0=310.0, ntr=0,
 ntpr=5000, ntwx=25000, ntwr=50000, ioutfm=1, ntxo=2, iwrap=1,
/
EOF

# ---------------------------------------------------------------------------
# Stage runner
# ---------------------------------------------------------------------------
run_stage () {
    local tag=$1 inp=$2 crd=$3 ref=$4
    if [ -f "${tag}.done" ]; then echo ">>> ${tag}: skipping"; return 0; fi
    echo ">>> ${tag}: start $(date +%T)"; local t0=$SECONDS
    if [ -n "$ref" ]; then
        pmemd.cuda -O -i "$inp" -o "${tag}.out" -p "$PRM" -c "$crd" \
                   -r "${tag}.rst" -x "${tag}.nc" -inf "${tag}.mdinfo" -ref "$ref"
    else
        pmemd.cuda -O -i "$inp" -o "${tag}.out" -p "$PRM" -c "$crd" \
                   -r "${tag}.rst" -x "${tag}.nc" -inf "${tag}.mdinfo"
    fi
    if [ $? -ne 0 ] || [ ! -f "${tag}.rst" ]; then
        echo "*** ${tag} FAILED"; tail -40 "${tag}.out"; exit 1
    fi
    touch "${tag}.done"
    echo ">>> ${tag}: done in $((SECONDS-t0)) s"
    grep -E "ns/day" "${tag}.out" | tail -1
}

run_stage min1    min1.in    "$INP"      "$INP"
run_stage min2    min2.in    min1.rst    ""
run_stage heat    heat.in    min2.rst    min2.rst
run_stage density density.in heat.rst    ""
run_stage equil   equil.in   density.rst ""

cp -f equil.rst system.rst7
echo "EQUILIBRATION COMPLETE"
grep -A9 "A V E R A G E S   O V E R" equil.out | tail -11
