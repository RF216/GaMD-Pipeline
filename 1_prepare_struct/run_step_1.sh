#!/bin/bash

python3 prep_pdb.py $1 complex_dry.pdb
chmod +x build_ff19_system.sh
bash build_ff19_system.sh
python3 finalize_topology.py system_raw.prmtop system.prmtop
