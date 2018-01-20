#! /bin/bash -f

#make output folder
inputFolder=$1
inputName=$2
Nsamples=$3

len=${#inputName}
NameBase=${inputName::len-11}
outputFolder=${NameBase}_${inputFolder}_Trajectories
mkdir ${outputFolder}

for ((sample=0; sample<=${Nsamples};sample++));
do
cp ${inputFolder}/tmpdir.${sample}/${inputName} ${outputFolder}/${NameBase}_traj_g_${sample}.txt
done 

