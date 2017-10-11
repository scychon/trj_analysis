export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd

cd $OPENMMDCD
ifort -qopenmp -c dcdmod.f90 trajmod.f90 -lxdrfile -L $GMXXTC/lib -I $GMXXTC -CB
cd testcond
ifort -qopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 conductivity.f90 -o calc_cond_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -CB