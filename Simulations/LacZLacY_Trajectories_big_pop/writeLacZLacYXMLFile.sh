#! /bin/bash -f

# writes corresponding xml file with given eps $1, method $2 and stores it at given location $3
eps=$1
method=$2
outputFolder=$3

scriptName=LacZLacY-${method}-${eps}.xml

#---------------------- xml content ---------
cat > ${scriptName} << EOF
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <annotation>
    <SimBiology xmlns="http://www.mathworks.com">
      <Version Major="5" Minor="4" Point="0"/>
    </SimBiology>
  </annotation>
  <model id="model1" name="lacy_lacz2_${method}">

      <annotation>
          <stochSim:inputData
              xmlns:stochSim="http://www.stochSim.org/ns"
              stochSim:TimeStart="0.0"
              stochSim:TimeEnd="100.0"
              stochSim:Method="${method}"
              stochSim:NumberOfSamples="1"
              stochSim:Epsilon="${eps}"
              stochSim:Theta="0.0"
              stochSim:Delta="0.05"
              stochSim:StoreInterval="25"
              stochSim:SortInterval="200"
              stochSim:InitialNoise="0.0"
              stochSim:NoiseIncrement="0.0"
              stochSim:NumberOfNoiseLevels="0"
              ></stochSim:inputData>
      </annotation>

    <listOfUnitDefinitions>
      <unitDefinition id="MWBUILTINUNIT_molecule" name="molecule">
        <listOfUnits>
          <unit kind="mole" exponent="1" multiplier="1.66053872801495e-24"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments>
      <compartment id="Comp1" name="unnamed" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="PLac" name="PLac" compartment="Comp1" initialAmount="100" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="RNAP" name="RNAP" compartment="Comp1" initialAmount="35" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="PLacRNAP" name="PLacRNAP" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrLacZ1" name="TrLacZ1" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="RbsLacZ" name="RbsLacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrLacZ2" name="TrLacZ2" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrLacY1" name="TrLacY1" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="RbsLacY" name="RbsLacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrLacY2" name="TrLacY2" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="Ribosome" name="Ribosome" compartment="Comp1" initialAmount="350" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="RbsribosomeLacZ" name="RbsribosomeLacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="RbsribosomeLacY" name="RbsribosomeLacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrRbsLacZ" name="TrRbsLacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="TrRbsLacY" name="TrRbsLacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="LacZ" name="LacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="LacY" name="LacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="dgrLacZ" name="dgrLacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="dgrLacY" name="dgrLacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="dgrRbsLacZ" name="dgrRbsLacZ" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="dgrRbsLacY" name="dgrRbsLacY" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="lactose" name="lactose" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="LacZlactose" name="LacZlactose" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="product" name="product" compartment="Comp1" initialAmount="50" substanceUnits="MWBUILTINUNIT_molecule" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
    </listOfSpecies>

<listOfReactions>
      <reaction id="R01" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="PLac" stoichiometry="1"/>
          <speciesReference species="RNAP" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="PLacRNAP" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k01 </ci>
              <ci> PLac </ci>
              <ci> RNAP </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k01" name="k01" value="0.17" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R02" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="PLacRNAP" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="PLac" stoichiometry="1"/>
          <speciesReference species="RNAP" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k02 </ci>
              <ci> PLacRNAP </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k02" name="k02" value="10" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
<reaction id="R03" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="PLacRNAP" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="TrLacZ1" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k03 </ci>
              <ci> PLacRNAP </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k03" name="k03" value="1" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R04" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrLacZ1" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="RbsLacZ" stoichiometry="1"/>
          <speciesReference species="PLac" stoichiometry="1"/>
          <speciesReference species="TrLacZ2" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k04 </ci>
              <ci> TrLacZ1 </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k04" name="k04" value="1" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R05" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrLacZ2" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="TrLacY1" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k05 </ci>
              <ci> TrLacZ2 </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k05" name="k05" value="0.015" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
 <reaction id="R06" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrLacY1" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="RbsLacY" stoichiometry="1"/>
          <speciesReference species="TrLacY2" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k06 </ci>
              <ci> TrLacY1 </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k06" name="k06" value="1" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R07" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrLacY2" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="RNAP" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k07 </ci>
              <ci> TrLacY2 </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k07" name="k07" value="0.36" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R08" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="Ribosome" stoichiometry="1"/>
          <speciesReference species="RbsLacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="RbsribosomeLacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k08 </ci>
              <ci> Ribosome </ci>
              <ci> RbsLacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k08" name="k08" value="0.17" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R09" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsribosomeLacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="Ribosome" stoichiometry="1"/>
          <speciesReference species="RbsLacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k09 </ci>
              <ci> RbsribosomeLacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k09" name="k09" value="0.45" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
 <reaction id="R10" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="Ribosome" stoichiometry="1"/>
          <speciesReference species="RbsLacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="RbsribosomeLacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k10 </ci>
              <ci> Ribosome </ci>
              <ci> RbsLacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k10" name="k10" value="0.17" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R11" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsribosomeLacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="Ribosome" stoichiometry="1"/>
          <speciesReference species="RbsLacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k11 </ci>
              <ci> RbsribosomeLacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k11" name="k11" value="0.45" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R12" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsribosomeLacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="TrRbsLacZ" stoichiometry="1"/>
          <speciesReference species="RbsLacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k12 </ci>
              <ci> RbsribosomeLacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k12" name="k12" value="0.4" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R13" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsribosomeLacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="TrRbsLacY" stoichiometry="1"/>
          <speciesReference species="RbsLacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k13 </ci>
              <ci> RbsribosomeLacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k13" name="k13" value="0.4" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
  <reaction id="R14" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrRbsLacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="LacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k14 </ci>
              <ci> TrRbsLacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k14" name="k14" value="0.015" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R15" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="TrRbsLacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="LacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k15 </ci>
              <ci> TrRbsLacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k15" name="k15" value="0.036" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R16" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="LacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="dgrLacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k16 </ci>
              <ci> LacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k16" name="k16" value="6.42e-05" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R17" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="LacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="dgrLacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k17 </ci>
              <ci> LacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k17" name="k17" value="6.42e-05" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
  <reaction id="R18" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsLacZ" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="dgrRbsLacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k18 </ci>
              <ci> RbsLacZ </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k18" name="k18" value="0.3" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R19" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="RbsLacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="dgrRbsLacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k19 </ci>
              <ci> RbsLacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k19" name="k19" value="0.3" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R20" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="LacZ" stoichiometry="1"/>
          <speciesReference species="lactose" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="LacZlactose" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k20 </ci>
              <ci> LacZ </ci>
              <ci> lactose </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k20" name="k20" value="9.52e-05" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction id="R21" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="LacZlactose" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="product" stoichiometry="1"/>
          <speciesReference species="LacZ" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k21 </ci>
              <ci> LacZlactose </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k21" name="k21" value="431" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
 <reaction id="R22" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="LacY" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="lactose" stoichiometry="1"/>
          <speciesReference species="LacY" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci> k22 </ci>
              <ci> LacY </ci>
            </apply>
          </math>
          <listOfParameters>
            <parameter id="k22" name="k22" value="14" constant="true"/>
          </listOfParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>


EOF
#---------------------- end of xml content -----
mv ${scriptName}  ${outputFolder}
