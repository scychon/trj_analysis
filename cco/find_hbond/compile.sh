#cd /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/
#ifort -c xtc-interface-wrap.f90 -lxdrfile -L lib/
#cp xtc.mod /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/cco/find_hbond

OPT="-openmp -static -check bounds -check uninit -check format -warn declarations -traceback"
OPT=""

ifort $OPT /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/xtc-interface-wrap.o variables.f90 hbond.f90 hbond_network.f90 -o hbond_network -lxdrfile -L /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/lib     
