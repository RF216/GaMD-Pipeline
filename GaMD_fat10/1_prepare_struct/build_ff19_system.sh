#!/bin/bash
set -u
cd "$(dirname "$0")"

module load Amber

# ---- PASS 1: solvate and count ----
cat > pass1.leap <<'EOF'
source leaprc.protein.ff19SB
source leaprc.water.opc
mol = loadpdb fat10_clean.pdb
charge mol
solvateOct mol OPCBOX 12.0
savepdb mol pass1.pdb
quit
EOF
tleap -f pass1.leap > pass1.out 2>&1
grep -E "^Exiting LEaP" pass1.out
grep -E "^Total unperturbed charge" pass1.out | head -1
grep -E "frcmod.opc|ionslm" pass1.out

python3 > ions.env <<'PYEOF'
nwat = sum(1 for ln in open("pass1.pdb")
           if ln.startswith(("ATOM","HETATM")) and ln[17:20].strip() in ("WAT","HOH")
           and ln[12:16].strip() == "O")
q = None
for ln in open("pass1.out"):
    if ln.startswith("Total unperturbed charge"):
        q = int(round(float(ln.split()[-1]))); break
V = nwat * 30.0 * 1e-27                       # litres of water
npair = int(round(0.15 * V * 6.02214076e23))  # 0.15 M NaCl
print(f"NWAT={nwat}"); print(f"Q={q}"); print(f"NPAIR={npair}")
print(f"NA={npair + (-q if q<0 else 0)}"); print(f"CL={npair + (q if q>0 else 0)}")
PYEOF
source ions.env
echo "waters $NWAT  charge $Q  ->  Na+ $NA / Cl- $CL"

# ---- PASS 2: final build ----
cat > pass2.leap <<EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
mol = loadpdb fat10_clean.pdb
saveamberparm mol dry.prmtop dry.inpcrd
solvateOct mol OPCBOX 12.0
addIonsRand mol Na+ $NA Cl- $CL
charge mol
check mol
savepdb mol system_raw.pdb
saveamberparm mol system_raw.prmtop system_raw.inpcrd
quit
EOF
tleap -f pass2.leap > pass2.out 2>&1
grep -E "^Exiting LEaP" pass2.out
grep -E "^Total unperturbed charge" pass2.out | tail -1
