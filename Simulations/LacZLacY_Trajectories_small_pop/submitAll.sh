#!/bin/sh
MethodsList="SLeaping TauLeaping AdaptiveTau"
#MethodsList="AdaptiveS AdaptiveTau SLeaping TauLeaping RLeapingJana"
EpsList="0.05 0.03 0.01"
myBase=${PWD}
echo My base is: ${myBase}

for method in $MethodsList; do
for eps in $EpsList; do

methodFolder=${method}/eps_${eps}/

cp ssm ${methodFolder}/bin
#cp run_lrz.sh ${methodFolder}

cd ${methodFolder}
sbatch runHybridJob_lrz.sh

cd ${myBase}
done 
done
