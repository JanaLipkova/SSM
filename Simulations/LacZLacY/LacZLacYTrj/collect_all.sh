#! /bin/bash -f
	

work_dir=$(pwd)

#MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "
MethodsList="RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "

EpsList=" 0.05 0.03 0.01 "

for method in $MethodsList; do

	cd ${method}
	cp ../collect_trajectories.sh .

	for eps in $EpsList; do

		trajFile=lacy_lacz2_Output.txt
		./collect_trajectories.sh  eps_${eps}   ${trajFile}    10000  &


	done

	cd $work_dir

done




