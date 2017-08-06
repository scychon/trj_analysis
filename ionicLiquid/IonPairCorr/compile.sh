cd /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/
#ifort -c xtc-interface-wrap.f90 -lxdrfile -L lib/
#cp xtc.mod /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid
ifort /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/xtc-interface-wrap.o variables.f90 topol.f90 calc_corr.f90 main_ILpair.f90 -o corr -lxdrfile -L /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/lib     
