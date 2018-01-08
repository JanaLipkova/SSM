#!/bin/sh

# 1) For each method and eps, creates folder
# 2) write xml file into the method folder
# 3) copy executables to the method folder


MethodsList="AdaptiveS AdaptiveTau SLeaping RLeapingJana TauLeaping SSA"
EpsList="0.05 0.03 0.01"

for method in $MethodsList; do
mkdir ${method}

for eps in $EpsList; do

#1) create folder
methodFolder=${method}/eps_${eps}/
mkdir ${methodFolder}

# 2) write xml file into the bin folder
./writeLacZLacYXMLFile.sh ${eps} ${method} ${methodFolder}

# 3) copy executables into the method folder
cp ssm ${methodFolder}
cp run_lrz.sh ${methodFolder}

done
done
