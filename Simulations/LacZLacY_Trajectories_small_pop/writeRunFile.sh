#! /bin/bash -f

XMLFile=$1
outputFolder=$2
scriptName=run.sh

cat > ${scriptName} << EOF
#! /bin/bash -f
./ssm ${XMLFile} \$1 \$2
./clean.sh
EOF

chmod 777 ${scriptName}
mv ${scriptName}  ${outputFolder}

