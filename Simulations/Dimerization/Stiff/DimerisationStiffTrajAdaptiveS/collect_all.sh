#! /bin/bash -f

EpsList="0.05 0.03 0.01"
trajFile=Dimerization_RLeapingJana_Output.txt

for eps in $EpsList; do
./collect_trajectories.sh eps_${eps} ${trajFile} 10024
done
