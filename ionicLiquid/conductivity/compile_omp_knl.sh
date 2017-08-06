export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc/knl
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd
cd knl
ifort -qopenmp -xMIC-AVX512 $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o ../variables.f90 ../topol.f90 ../fvector.f90 ../filehandle.f90 ../conductivity.f90 -o calc_cond_omp_knl -lxdrfile -L $GMXXTC/../lib -I $GMXXTC -I $OPENMMDCD -CB
