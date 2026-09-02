#!/usr/bin/env python3
"""Structural audit of a protein-conjugate PDB prior to AMBER setup."""
import math, sys
from collections import OrderedDict, Counter
 
PDB = sys.argv[1]
STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
 
atoms = []
with open(PDB) as fh:
    for ln in fh:
        if ln.startswith(("ATOM  ", "HETATM")):
            atoms.append(dict(
                rec=ln[0:6].strip(), serial=int(ln[6:11]), name=ln[12:16].strip(),
                altloc=ln[16], resn=ln[17:20].strip(), chain=ln[21],
                resi=int(ln[22:26]), icode=ln[26],
                x=float(ln[30:38]), y=float(ln[38:46]), z=float(ln[46:54]),
                elem=ln[76:78].strip()))
 
def d(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)
 
print(f"atoms {len(atoms)}   HETATM {sum(1 for a in atoms if a['rec']=='HETATM')}")
alt = Counter(a["altloc"] for a in atoms if a["altloc"] != " ")
print(f"altLoc {dict(alt) if alt else 'none'}")
 
res = OrderedDict()
for a in atoms:
    res.setdefault((a["chain"], a["resi"], a["icode"]), {"resn": a["resn"], "atoms": []})["atoms"].append(a)
 
chains = OrderedDict()
for (ch, ri, ic), v in res.items():
    chains.setdefault(ch, []).append((ri, v["resn"], v["atoms"]))
 
print("\n--- chains ---")
for ch, rl in chains.items():
    nums = [r[0] for r in rl]
    print(f"chain {ch}: {len(rl)} residues {min(nums)}..{max(nums)}, "
          f"first {rl[0][1]}{rl[0][0]}, last {rl[-1][1]}{rl[-1][0]}")
    nonstd = sorted({r[1] for r in rl if r[1] not in STD})
    print(f"   non-standard: {nonstd if nonstd else 'none'}")
 
print("\n--- chain breaks (numbering gaps or C(i)-N(i+1) > 1.8 A) ---")
for ch, rl in chains.items():
    breaks = []
    for i in range(len(rl)-1):
        ri, rn, ra = rl[i]; rj, rnj, rb = rl[i+1]
        C = next((a for a in ra if a["name"] == "C"), None)
        N = next((a for a in rb if a["name"] == "N"), None)
        dist = d(C, N) if (C and N) else None
        if rj != ri+1 or dist is None or dist > 1.8:
            breaks.append(f"{rn}{ri}->{rnj}{rj} ({dist:.2f} A)" if dist else f"{rn}{ri}->{rnj}{rj} (missing)")
    print(f"chain {ch}: {breaks if breaks else 'continuous'}")
 
print("\n--- disulfides (SG-SG < 3.0 A) ---")
sgs = [(ch, ri, a) for (ch, ri, ic), v in res.items() if v["resn"] == "CYS"
       for a in v["atoms"] if a["name"] == "SG"]
print(f"CYS count {len(sgs)}: {[f'{c}{r}' for c,r,_ in sgs]}")
found = False
for i in range(len(sgs)):
    for j in range(i+1, len(sgs)):
        dd = d(sgs[i][2], sgs[j][2])
        if dd < 3.0:
            print(f"  SS: {sgs[i][0]}{sgs[i][1]} - {sgs[j][0]}{sgs[j][1]}  {dd:.2f} A"); found = True
if not found:
    print("  none -> all CYS reduced; use CYS not CYX")
 
print("\n--- histidines (tleap default = HIE) ---")
print([f"{ch}{ri}" for (ch, ri, ic), v in res.items() if v["resn"].startswith("HI")])
 
print("\n--- charge estimate ---")
cnt = Counter(v["resn"] for v in res.values())
print(f"ARG {cnt['ARG']}  LYS {cnt['LYS']}  ASP {cnt['ASP']}  GLU {cnt['GLU']}  HIS {cnt['HIS']}")
print(f"side-chain net = {cnt['ARG']+cnt['LYS']-cnt['ASP']-cnt['GLU']:+d}")
 
print("\n--- extent ---")
xs=[a['x'] for a in atoms]; ys=[a['y'] for a in atoms]; zs=[a['z'] for a in atoms]
print(f"{max(xs)-min(xs):.1f} x {max(ys)-min(ys):.1f} x {max(zs)-min(zs):.1f} A")
