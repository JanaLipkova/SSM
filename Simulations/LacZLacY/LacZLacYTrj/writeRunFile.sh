#! /bin/bash -f

XMLFile=$1
outputFolder=$2
scriptName=run.sh

cat > ${scriptName} << EOF
#! /bin/bash -f
module load gcc boost
./ssm ${XMLFile} \$1 \$2
rm ssm output run.sh LacZLacY-SSA-0.01.xml lacy_lacz2_histogram.txt
EOF

chmod 777 ${scriptName}
mv ${scriptName}  ${outputFolder}

