export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd

#ifort -qopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 velanal.f90 -o vel_anal_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -CB
gfortran -fopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 velanal.f90 -o vel_anal_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -fcheck=bounds -ffree-line-length-256
