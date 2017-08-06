cd /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/
#ifort -c xtc-interface-wrap.f90 -lxdrfile -L lib/
#cp xtc.mod /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/
cd /oasis/projects/nsf/wis119/scychon/program/utility/trj_analysis/ionicLiquid/conductivity/testidx
ifort /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/xtc-interface-wrap.o variables.f90 topol.f90 fvector.f90 filehandle.f90 conductivity.f90 -o calc_cond -lxdrfile -L /oasis/projects/nsf/wis119/scychon/program/gmx_xtc/lib -CB
