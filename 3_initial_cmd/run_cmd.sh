#!/bin/bash
pmemd.cuda -O -i prod.in -p ../1_prepare_struct/system.prmtop -c ../2_system_prep/system.rst7 -o system.out -x system.nc -r system_cmd.rst7 -inf system.mdinfo
