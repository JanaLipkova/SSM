#! /bin/bash -f


work_dir=$(pwd)

MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "
#MethodsList="RLeaping"

EpsList=" 1.0 0.5 0.1 "



for method in $MethodsList; do

    cd ${method}

    for eps in $EpsList; do

        folder=BSubtilis_eps_${eps}_Trajectories

        echo "${method}  ${eps}"

        tar -xzf ${folder}.tar.gz &

    done

    cd ${work_dir}

done
