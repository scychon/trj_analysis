#!/bin/bash

# Set some environment variables 
#../hbond_network -f ../gromacs_new/pr_popc/param.dat
#../hbond_network -f ../gromacs_new/pr_from_dpdoo/param.dat
#../hbond_network -f ../gromacs_new/ppdoo/3H2O/param.dat
#../hbond_network -f ../gromacs_new/ppdoo/Xtal/param.dat

# gromacs
../hbond_network -f ../gromacs/pr_popc/param.dat
../hbond_network -f ../gromacs/ppdoo/param.dat
../hbond_network -f ../gromacs/ppdoo_from_pr/param.dat
../hbond_network -f ../gromacs_new/dpdoo/param.dat
../hbond_network -f ../gromacs_new/pr_from_ppdoo_PRA/param.dat
