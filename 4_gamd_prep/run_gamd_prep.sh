#!/bin/bash
pmemd.cuda -O -i md_prep.in -p ../1_prepare_struct/system.prmtop -c ../3_initial_cmd/system_cmd.rst7 -o system_gamd_prep.out -x system_gamd_prep.nc -r system_gamd_prep.rst7 -gamd system_gamd_prep.log
