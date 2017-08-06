export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
ifort -qopenmp -c dcdmod.f90 trajmod.f90 -lxdrfile -L $GMXXTC/lib -I $GMXXTC -CB -convert little_endian
