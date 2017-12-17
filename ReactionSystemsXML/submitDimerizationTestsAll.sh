#!/bin/sh

MethodsList="SLeaping RLeapingJana TauLeaping TauLeaping-NN"
EpsList="0.05 0.03 0.01"
myBase=${PWD}
echo My base is: ${myBase}

for method in $MethodsList; do
for eps in $EpsList; do

methodFolder=${method}/eps_${eps}/

cp ssm ${methodFolder}
cp run_lrz.sh ${methodFolder}

cd ${methodFolder}
sbatch run_lrz.sh
echo submitted method: ${method} with eps ${eps}

cd ${myBase}
done 
done
