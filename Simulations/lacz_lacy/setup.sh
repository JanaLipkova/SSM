ROOT="/Users/garampat/Desktop/ETH-work/projects/2017-sleap/SSM" 

SBML_BASE="lacz_lacy"

METHOD="SSA"

EPS="0.01"


#==================================================================

for method in $METHOD; do
for eps in $EPS; do

	SSM_EXEC="${ROOT}/ssm"
	PROP_EXEC="${ROOT}/pi4u_lite/tools/simulate_all"

	NOW_DIR=`pwd`


	if [ -d "$method" ]; then
		echo "Base directory '${method}' exists. Delete and rerun."
		exit 1
	fi


	mkdir $method
	mkdir "${method}/${eps}"
	SIM_DIR="${method}/${eps}"


	python write_xml.py "${SBML_BASE}.xml" ${METHOD} ${EPS}  


	cd $SIM_DIR

	mkdir bin
	SBML_NEW="${SBML_BASE}_${METHOD}_${EPS}.xml"
	mv "${NOW_DIR}/${SBML_NEW}" bin/. 

	cp ${SSM_EXEC} bin/.

	cp ${PROP_EXEC} .


	cat > run.sh << EOF
#! /bin/bash -f

module load python
./ssm $SBML_NEW \$1 \$2
rm ssm
EOF

	chmod u+x run.sh
	mv run.sh bin/.


	cd $NOW_DIR



done
done

