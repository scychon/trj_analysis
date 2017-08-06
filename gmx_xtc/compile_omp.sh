module load intel/compiler 
export CC=icc
export FC=ifort
export GMXXTC=$HOME/program/utility/trj_analysis/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trj_analysis/openmm_dcd
cd $GMXXTC/xdrfile-1.1.1
./configure --prefix=$GMXXTC/
make
make install
