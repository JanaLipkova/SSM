#!/bin/sh
EpsList="0.05 0.03 0.01"

for eps in $EpsList; do
cd eps_${eps}/bin
./ssm *.xml > output
cd ../../
done
