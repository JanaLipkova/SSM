#!/bin/sh



MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "

EpsList=" 0.05 0.03 0.01 "


work_dir=$(pwd)


for method in $MethodsList; do


	for eps in $EpsList; do

		methodFolder=${method}/eps_${eps}/

		cd $methodFolder

		bsub   -n 24   -W 1:00     -o lsf_out 'mpirun ./simulate_all 10000'	
		
		cd $work_dir

	done
done
