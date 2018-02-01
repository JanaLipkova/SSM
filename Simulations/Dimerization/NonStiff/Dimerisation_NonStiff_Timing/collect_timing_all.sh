#!/bin/sh
MethodsList="SSA SLeaping TauLeaping AdaptiveS AdaptiveTau RLeapingJana"

for method in $MethodsList; do

echo "collecting time for ${method}"

# collect times
t1=$( ./collect_time.sh ${method}/eps_0.05/) 
t2=$( ./collect_time.sh ${method}/eps_0.03/)
t3=$( ./collect_time.sh ${method}/eps_0.01/)

echo $t1
echo $t2
echo $t3
#write times into file
scriptName=Dimerization-NonStiff-${method}_times.txt

cat > ${scriptName} << EOF
0.05   ${t1}
0.03   ${t2}
0.01   ${t3}
EOF

done


