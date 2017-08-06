#!/bin/bash

# Set some environment variables 
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pmpro/from_ppdoo_solvate
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pr_popc
../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pr_from_dpdoo
../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/ppdoo/3H2O
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/ppdoo/Xtal
../../../hbond_network -f param.dat
