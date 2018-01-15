#!/bin/sh

# 1) For each method and eps, creates folder
# 2) write xml file into the bin folder
# 3) write run.sh into the bin
# 4) copy executables into the bin folder
# 5) copy executables to the main folder

MethodsList="SSA"
#MethodsList="AdaptiveS AdaptiveTau SLeaping RLeapingJana TauLeaping SSA"
EpsList="0.05 0.03 0.01"

for method in $MethodsList; do
mkdir ${method}

for eps in $EpsList; do

#1) create folder
methodFolder=${method}/eps_${eps}/
binFolder=${methodFolder}/bin
mkdir ${methodFolder}
mkdir ${binFolder}

# 2) write xml file into the bin folder
./writeLacZLacYXMLFile.sh ${eps} ${method} ${binFolder}

# 3) write run.sh file
XMLscriptName=LacZLacY-${method}-${eps}.xml
./writeRunFile.sh ${XMLscriptName} ${binFolder} 

# 4) copy executables into the bin folder
cp ssm ${binFolder}

# 5) copy executables into the main folder
cp runHybridJob_lrz.sh ${methodFolder}
cp simulate_all ${methodFolder}

done
done
