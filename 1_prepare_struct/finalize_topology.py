#!/usr/bin/env python3
"""
finalize_topology.py <raw.prmtop> <out.prmtop>

Hydrogen mass repartitioning for a 4 fs timestep.
"""
import sys
import parmed as pmd
from parmed.tools import HMassRepartition

RAW, OUT = sys.argv[1], sys.argv[2]

parm = pmd.load_file(RAW, xyz=RAW.replace(".prmtop", ".inpcrd"))
print(f"atoms {len(parm.atoms):,}  net charge {sum(a.charge for a in parm.atoms):+.4f}")
apm = parm.parm_data["ATOMS_PER_MOLECULE"]
print(f"molecules {len(apm):,}, largest {max(apm):,} atoms")

print("\n--- HMR ---")
m0 = sum(a.mass for a in parm.atoms)
sh = [a for a in parm.atoms if a.atomic_number == 1 and a.residue.name not in ("WAT", "HOH")]
wh = [a for a in parm.atoms if a.atomic_number == 1 and a.residue.name in ("WAT", "HOH")]
HMassRepartition(parm, 3.024).execute()
print(f"total mass {m0:,.3f} -> {sum(a.mass for a in parm.atoms):,.3f} Da")
print(f"solute H {sh[0].mass:.4f} Da   water H {wh[0].mass:.4f} Da (must stay 1.008)")

parm.save(OUT, overwrite=True)
chk = pmd.load_file(OUT)
ok = abs(sum(a.charge for a in chk.atoms)) < 1e-3
print(f"\nre-read: {len(chk.atoms):,} atoms, charge {sum(a.charge for a in chk.atoms):+.4f}")
print("OVERALL:", "PASS" if ok else "FAIL")
