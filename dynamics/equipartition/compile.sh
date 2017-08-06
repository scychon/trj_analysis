export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd

#ifort -openmp $GMXXTC/xtc-interface-wrap.o variables.f90 topol.f90 fvector.f90 filehandle.f90 conductivity.f90 -o calc_cond_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -CB
ifort $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 velanal.f90 -o vel_anal -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -CB
