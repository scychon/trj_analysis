#!/bin/bash
#SBATCH --job-name="cond_racf"  
#SBATCH --error="job.%j.err"  
#SBATCH --output="job.%j.out"  
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --export=ALL  
#SBATCH -t 10:00:00  

# Set some environment variables 

CALC=/oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/conductivity

SETUP=/oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/conductivity
BMIMBF4=/oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/conductivity
cd $BMIMBF4
#echo "0 " | trjconv_mpi -f nvt_300K_4.xtc -pbc mol -o dielect_300.xtc -s ../jpcl/300K/simul/jpcl_nvt.tpr
ibrun -v ./calc_cond -f param.dat
