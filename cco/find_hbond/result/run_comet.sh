#!/bin/bash

# Set some environment variables 
../hbond_network -f ../gromacs_new/pmpro/from_ppdoo_solvate/param.dat
../hbond_network -f ../gromacs_new/pmpro/ppdoo_3h2o/param.dat
../hbond_network -f ../gromacs_new/pmpro/3h2o/param.dat
../hbond_network -f ../gromacs_new/pmpro/xtal_solvate/param.dat

../hbond_network -f ../gromacs_new/pm/3h2o/param.dat
../hbond_network -f ../gromacs_new/pm/ppdoo_solvate/param.dat
../hbond_network -f ../gromacs_new/pm/xtal_solvate/param.dat
