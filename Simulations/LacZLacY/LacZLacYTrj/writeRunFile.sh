#! /bin/bash -f

XMLFile=$1
outputFolder=$2
scriptName=run.sh

cat > ${scriptName} << EOF
#! /bin/bash -f
module load gcc boost
./ssm ${XMLFile} \$1 \$2
rm ssm output run.sh 
files=`ls -1 ./*.xml`
for f in $files; do
	rm $f
done
files=`ls -1 ./*histogram.txt`
for f in $files; do
	rm $f
done
EOF

chmod 777 ${scriptName}
mv ${scriptName}  ${outputFolder}

