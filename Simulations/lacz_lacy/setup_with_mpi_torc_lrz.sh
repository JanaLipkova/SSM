# Setup folders for the simulation of ssm.

ROOT="/home/hpc/txh01/di49zin/SSM"
SBML_BASE="lacz_lacy"
METHOD="SSA TAU" # populate with more methos
EPS="0.01 0.02" # populate with more values for epsilon
#==================================================================

for method in $METHOD; do	
	if [ -d "$method" ]; then
		echo "Base directory '$method' exists. Delete and rerun."
		exit 1
	fi
	mkdir $method

for eps in $EPS; do

	SSM_EXEC="${ROOT}/make_file/ssm"
	PROP_EXEC="${ROOT}/pi4u_lite/tools/simulate_all"

	NOW_DIR=`pwd`

	SIM_DIR="${method}/${eps}"


	if [ -d "$SIM_DIR" ]; then
		echo "Base directory '$SIM_DIR' exists. Delete and rerun."
		exit 1
	fi

	mkdir "${method}/${eps}"

	python write_xml.py "${SBML_BASE}.xml" ${method} ${eps}  


	cd $SIM_DIR

	mkdir bin
	SBML_NEW="${SBML_BASE}_${method}_${eps}.xml"
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

