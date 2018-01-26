#! /bin/bash -f


work_dir=$(pwd)

MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  SLeaping_v6 "

EpsList=" 0.05 0.03 0.01 "



for method in $MethodsList; do

    cd ${method}

    for eps in $EpsList; do

        folder=lacy_lacz2_${method}_eps_${eps}_Trajectories

        echo "${method}  ${eps}"

        tar -xzf ${folder}.tar.gz &

    done

    cd ${work_dir}

done
