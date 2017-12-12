#===========================
# Running SSM code
#===========================
# ./StochasticSimulationMethods path_to_your_xml_file
# For parallel runs remeber to specify number of threads:

export OMP_NUM_THREADS=48;
bsub -W 08:00 -n 48 -o SSA ./StochasticSimulationMethods ReactionSystems/Dimerization/NonStiff/Dimerization-SSA.xml

