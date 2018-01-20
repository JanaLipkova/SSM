#!/bin/sh


cp $HOME/SSM/ssm .
cp $HOME/SSM/pi4u_lite/tools/simulate_all .



MethodsList=" SSA  TauLeap  AdaptiveTau  RLeaping  SLeaping_v3  SLeaping_v4  SLeaping_v5  AdaptiveS "

EpsList=" 0.05 0.03 0.01 "

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
	#cp runHybridJob_lrz.sh ${methodFolder}
	cp simulate_all ${methodFolder}

	done
done
