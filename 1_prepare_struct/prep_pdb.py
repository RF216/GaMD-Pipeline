#!/usr/bin/env python3
"""
Prepare a conjugate PDB for tleap.
 
Strip all hydrogens and heteroatoms (tleap rebuilds them with correct ff naming).
"""
import sys

SRC = sys.argv[1]
OUT = sys.argv[2]

out, serial = [], 0
with open(SRC) as fh:
    for ln in fh:
        if not ln.startswith("ATOM  "):
            continue
        elem = ln[76:78].strip()
        name = ln[12:16].strip()
        if elem == "H" or (name and (name[0] == "H" or (name[0].isdigit() and "H" in name))):
            continue
        serial += 1
        out.append(f"{ln[0:6]}{serial:5d}{ln[11:]}".rstrip("\n"))

serial += 1
out.append(f"TER   {serial:5d}")
out.append("END")

open(OUT, "w").write("\n".join(out) + "\n")
print(f"wrote {OUT}: {sum(1 for l in out if l.startswith('ATOM'))} heavy atoms")
