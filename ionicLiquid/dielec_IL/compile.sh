export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd
#cp xtc.mod /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/
#cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/dielec_IL
ifort -qopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 calc_corr.f90 main_ILpair.f90 -o corr_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -CB
#ifort /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/xtc-interface-wrap.o variables.f90 topol.f90 calc_corr.f90 main_ILpair.f90 -o corr -lxdrfile -L /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/lib     
