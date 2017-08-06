#!/bin/bash
#SBATCH --job-name="structure"  
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH -p normal
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --no-requeue


# Set some environment variables 
module load intel/17.0.0
module load impi/17.0.0
export PATH=/home1/02572/scychon/program/gmx-4.6.5_knl/bin:$PATH

#export OPM_STACKSIZE=100M
export OMP_NUM_THREADS=8

STR=/home1/02572/scychon/program/utility/trj_analysis/structure_code_FFT/il_ua
ILROOT=/work/02572/scychon/IonicLiquid/ua_il_jpcb/

## declare an array variable
declare -a arril=("emimbf4" "c6mimbf4" "c8mimbf4" "c10mimbf4" "c12mimbf4" "c14mimbf4" "c16mimbf4" "emimpf6" "bmimpf6" "c6mimpf6" "c8mimpf6" "c10mimpf6" "c12mimpf6" "c14mimpf6" "c16mimpf6")
declare -a arrtemp=("300" "400")

# get length of an array
arraylength=${#arril[@]}

# generate the grofiles for each trajectory
idx=0
for (( i=0; i<${arraylength}; i++ ))
do
	for (( j=0; j<2; j++ ))
	do
		ILDIR="$ILROOT${arril[$i]}/${arrtemp[$j]}K/simul"
		#ILFNAME="${arrgro[$i]/temp/${arrtemp[$j]}}"
		ILFNAME="md_nvt"
		IL="$ILDIR/rename.gro"
		ILTRJ="$ILDIR/md_nvt_trj.gro"
		cd $ILDIR
		#echo "$IL"
		#echo "$ILFNAME"
		#rm "${ILFNAME}_trj.gro"
		if [ ! -f "${ILFNAME}_trj.gro" ]; then
			echo "0" | ibrun -n 1 -o $idx trjconv_mpi -f "$ILFNAME.xtc" -o "${ILFNAME}_trj.gro" -s $ILFNAME -dt 40 &
			((idx++))
		fi
	done
  #echo $i " / " ${arraylength} " : " ${array[$i-1]}
done

wait

cd $STR
pwd

idx=0
for (( i=0; i<${arraylength}; i++ ))
do
	#rm -r ${arril[$i]}
	mkdir ${arril[$i]}
	for (( j=0; j<2; j++ ))
	do
		mkdir "${arril[$i]}/${arrtemp[$j]}K"
		ILDIR="$ILROOT${arril[$i]}/${arrtemp[$j]}K/simul"
		IL="$ILDIR/rename.gro"
		ILTRJ="$ILDIR/md_nvt_trj.gro"
		echo "$IL"
		echo "ibrun -n 1 -o $idx ./main_structure_knl $IL $ILTRJ | tee ${arril[$i]}/${arrtemp[$j]}K/sq_40ps"
		#ibrun -n 1 -o $idx ./main_structure_knl $IL $ILTRJ | tee "${arril[$i]}/${arrtemp[$j]}K/sq_40ps" &
		taskset -c ${idx}-$((idx+7)) ./main_structure_knl $IL $ILTRJ | tee "${arril[$i]}/${arrtemp[$j]}K/sq_40ps_2" &
		idx=$((idx+8))
	done
  echo $i " / " ${arraylength} " : " ${array[$i-1]}
done

wait

idx=0
for (( i=0; i<${arraylength}; i++ ))
do
	#mkdir ${arril[$i]}
	for (( j=0; j<2; j++ ))
	do
		#mkdir "${arril[$i]}/${arrtemp[$j]}K"
		ILDIR="$ILROOT${arril[$i]}/${arrtemp[$j]}K/simul"
		IL="$ILDIR/$ILFNAME.gro"
		taskset -c $idx perl average_Sq.pl "${arril[$i]}/${arrtemp[$j]}K/sq_40ps_2" > "${arril[$i]}/${arrtemp[$j]}K/sq_avg_40ps" &
		((idx++))
	done
  echo $i " / " ${arraylength} " : " ${array[$i-1]}
done

wait
