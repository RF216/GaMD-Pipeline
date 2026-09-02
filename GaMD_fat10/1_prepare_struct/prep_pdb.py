#!/usr/bin/env python3
"""
Prepare a conjugate PDB for tleap.
 
1. Strip all hydrogens (tleap rebuilds them with correct ff naming).
2. Rename the acceptor lysine -> LYX and the donor glycine -> GLL.
   Because LYX/GLL are absent from tleap's PDB residue map, NEITHER gets an
   automatic N-/C-terminal variant. That is the whole trick: the donor Gly
   keeps its carbonyl and is never given a terminal OXT.
3. Print the sequential residue indices tleap will assign, so the `bond`
   command is derived rather than guessed.
"""
import sys
from collections import OrderedDict
 
SRC = sys.argv[1]
OUT = sys.argv[2]
# (chain, resseq) -> new name.  EDIT THESE TWO LINES FOR YOUR SYSTEM.
RENAME = {("B", 323): "LYX", ("C", 165): "GLL"}
 
records = []
with open(SRC) as fh:
    for ln in fh:
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        elem = ln[76:78].strip(); name = ln[12:16].strip()
        if elem == "H" or (name and (name[0] == "H" or (name[0].isdigit() and "H" in name))):
            continue
        records.append(ln.rstrip("\n"))
 
residues = OrderedDict()
for ln in records:
    residues.setdefault((ln[21], int(ln[22:26]), ln[26]), []).append(ln)
 
out, serial, prev_chain, idx, report = [], 0, None, 0, {}
for (chain, resi, icode), atom_lines in residues.items():
    if prev_chain is not None and chain != prev_chain:
        serial += 1; out.append(f"TER   {serial:5d}")
    idx += 1
    new = RENAME.get((chain, resi))
    for ln in atom_lines:
        serial += 1
        resn = new if new else ln[17:20]
        out.append(f"{ln[0:6]}{serial:5d}{ln[11:17]}{resn:>3s}{ln[20:76]}{ln[76:78]:>2s}")
    if new:
        report[new] = (chain, resi, idx)
    prev_chain = chain
serial += 1
out.append(f"TER   {serial:5d}"); out.append("END")
 
open(OUT, "w").write("\n".join(out) + "\n")
print(f"wrote {OUT}: {sum(1 for l in out if l.startswith('ATOM'))} heavy atoms, {idx} residues\n")
for resn, (ch, ri, i) in report.items():
    print(f"  {resn} <- {ch}{ri} -> tleap residue #{i}")
if "LYX" in report and "GLL" in report:
    print(f"\nBOND COMMAND:\n  bond mol.{report['LYX'][2]}.NZ mol.{report['GLL'][2]}.C")
