#! /bin/bash -f

XMLFile=$1
outputFolder=$2
method=$3
scriptName=clean.sh

cat > ${scriptName} << EOF
rm ssm
rm ${XMLFile}
rm Bsubtilis_${method}_histogram.txt
EOF

chmod 777 ${scriptName}
mv ${scriptName}  ${outputFolder}
