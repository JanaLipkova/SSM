#! /bin/bash -f


work_dir=$(pwd)

MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "
#MethodsList="RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "

EpsList=" 0.05 0.03 0.01 "

mkdir pack
cd pack

for method in $MethodsList; do

	mkdir ${method}

	for eps in $EpsList; do

		echo "Copy ${method}  ${eps}"

		folder=lacy_lacz2_eps_${eps}_Trajectories

		cp ../${method}/${folder}.tar.gz ${method}/.

	done

done


cd ${work_dir}

cp  extract.sh  pack/.

tar -czf pack.tar.gz pack
