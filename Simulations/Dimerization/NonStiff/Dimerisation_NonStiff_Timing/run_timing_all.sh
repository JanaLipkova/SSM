#!/bin/sh
MethodsList="TauLeaping"
#MethodsList="SLeaping RLeapingJana AdaptiveTau AdaptiveS"
EpsList="0.05 0.03 0.01"
myBase=${PWD}
echo My base is: ${myBase}

for method in $MethodsList; do
for eps in $EpsList; do

	methodFolder=${method}/eps_${eps}/
	cp ssm ${methodFolder}
	cd ${methodFolder}
	
	echo "Running ${method} with ${eps} " 
	./ssm *.xml > output

	cd ${myBase}
done
done

method="SSA"
eps="0.05"

methodFolder=${method}/eps_${eps}/
cp ssm ${methodFolder}
cd ${methodFolder}

echo "Running ${method} with ${eps} "
./ssm *.xml > output

 cd ${myBase}

