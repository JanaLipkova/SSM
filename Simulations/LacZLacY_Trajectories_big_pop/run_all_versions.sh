#!/bin/sh
MethodsList="TauLeaping AdaptiveTau RLeapingJana SLeaping_v3 SLeaping_v4 SLeaping_v5"
for method in $MethodsList; do
cd ${method}
./run_all.sh
cd ../
done
