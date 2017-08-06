#!/bin/bash
#SBATCH --job-name="ua_conductivity"  
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH -p development
##SBATCH -p normal
#SBATCH -n 68
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH --no-requeue


# Set some environment variables 
module load intel/17.0.0
module load impi/17.0.0

#export OPM_STACKSIZE=100M
export OMP_NUM_THREADS=64
export OMP_STACKSIZE=500M

CALC=$HOME/program/utility/trj_analysis/openmm_dcd/testcond/knl/
ILROOT=/work/02572/scychon/IonicLiquid/openmm/neat_il_aa/bmimno3/300K/testcond/bmimno3

ncores=68

cd $ILROOT
time -p $CALC/calc_cond_omp_knl -f param_300.dat >& log_calc_300
