#!/bin/bash

# Set some environment variables 
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pm/xtal_solvate
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pm/3h2o
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pm/ppdoo_3h2o
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pm/ppdoo_solvate
../../../hbond_network -f param.dat
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond/gromacs_new/pmpro/from_ppdoo_solvate
../../../hbond_network -f param.dat
