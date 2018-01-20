#! /bin/bash -f
	

work_dir=$(pwd)

MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "
#MethodsList="RLeaping"

EpsList=" 0.05 0.03 0.01 "


mkdir pack


for method in $MethodsList; do

	cd ${method}

	for eps in $EpsList; do

		folder=lacy_lacz2_eps_${eps}_Trajectories
	
		echo "${method}  ${eps}"
		ls -1 ${folder} | wc -l 
	
		tar -czf ${folder}.tar.gz ${folder} &

	done

	cd ${work_dir}

done
