export TRAJMOD=$HOME/program/utility/trajmod
export GMXXTC=$TRAJMOD/gmx_xtc
export OPENMMDCD=$TRAJMOD/openmm_dcd

cd equipartition
bash compile_omp.sh
cd ../ionicLiquid/IonPairCorr
bash compile.sh
cd ../conductivity
bash compile_omp.sh
cd ../dielec_IL
bash compile.sh
