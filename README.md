# Stochastic Simulation Methods
* Software for stochastic simulations of well-mixed continuous time Markov chain models for chemical reaction networks.
* The software takes as input .xml file with the list of chemical reactions to be modelled and simulations parameters
* The supported methods are:
    * **SSA** (Gillespie: *Exact stochastic simulation of coupled chemical reaction* 1977)
   * **Tau-leaping** (Cao et. al: *Avoiding negative populations in explicit poisson tau-leaping* 2005)
   * **R-leaping**   (Auger et al: *R-leaping: accelerating the stochastic simulation algorithm by reaction leaps*, 2006)
   * **S-leaping**   (Lipkova et al: *S-Leaping: An adaptive, accelerated stochastic simulation algorithm, bridging τ-leaping and R-leaping*, 2018)
   * **Adaptive Tau-leaping** (Cao et al: *Adaptiveexplicit-implicittau-leapingmethod with automatic tau selection*, 2007)
   * **Adatpive S-leaping** (Lipkova et al: *S-Leaping: An adaptive, accelerated stochastic simulation algorithm, bridging τ-leaping and R-leaping*, 2018)

### Installation
Following prerequisites are require:
* Expat http://expat.sourceforge.net/
* Blitz++ http://sourceforge.net/projects/blitz/
* Boost http://www.boost.org/
* LibSBML http://sourceforge.net/projects/sbml/files/libsbml/
* C++ compiler (version gcc47 or newer is required to support C++11 random number generator)

If you are install above libraries from scratch, we recommend to install them into one folder. In the provided makefile, all libraries are installed in folder called `ssm-libs`.

**Compilation commands for the above libraries in UNIX enviroment:**
*Expat and Blitz++:*
```sh
./configure
sudo make install
```
*Boost*
```sh
sudo ./bootstrap.sh
sudo ./bjam
```
*LibSBML*
```sh
./configure --with-expat --with-python --with-java
sudo make install
```

### Compilation & Execution
**1)  Set up path to your libraries:**
    The folder `SSM/make_file` contains the main makefile called `Makefile` and small enviroment dependent make files called i.e. `make.anslab`. In the `make.anslab` set up path to your installed libraries. If you installed all libraries in `ssm-libs`, then you just need to change the path in `MYBASE` variable.

**2)  Export path to your libraries:**
You can either export path to the libries in your `.bash_profile` (or similar) or set path to your libraies in the setup file `SSM/make_file/setup_anslab.sh` 

**3)  Compile the code:**
```sh
cd SSM/make_file
source setup_anslab.sh
make clean
make
```
It will create executable called `ssm`


**4) Executation:**
Run the code as follows:
```sh
./ssm InputFile.xml
```
where the file `InputFile.xml` contains informaiton about the system of reactions to be simulated and corresponding input parameters

### Anatomy of Input File
Folder `SSM/ReactionSystemsXML` contains examples of XML file for several reactions systems. Each input file consists of four parts:   

1) 	`<AnnotationField>`  defines input parameters of the method:

    | Parameter        | Description     |
    | ------------- |:-----------------:|
    | `model name`  | Name of your simulation, is used to name the output files. Use a name without spaces!|
    |`TimeStart`    | Intial simulation time |
    | `TimeEnd`     | Final simulation |
    |`Method`| Method to be used, i.e. SSA, TauLeaping,..|
    |`NumberOfSamples`| Number of samples |
    |`numberOfNoiseLevels`| not relevant here, keep as it is or remove|
    |`Epsilon`| The accuracy parameter for the leap methods|
    |`Theta`| The control parameter for avoiding negative population in R- and S-leaping |
    |`StoreInterval`| Number of equally spaced time points in which to save the simulation |
    |`SortInteval`| How often to reorder the channels in the R- and S-leaping mehtod|
2) `<listOfCompartments>` degines volume and amount of compartments (for most simulations only one compartment of volume 1 is used)
3) `<listOfSpecies>` defines initial population
4) `<listOfReactions>` defines system of reactions to be modelled

For more reactions systems (already written in the XML format) see http://www.ebi.ac.uk/biomodels-main/

### Output Files
Each simulation produce two output files with extension:
1) *_Output.txt* : stores trajectory of each species avereged over number of samples, where each column correspond to one species and each line to time point. Control the number of time points by parameter `stochSim:StoreInterval` in the `<AnnotationField>`. The species are reporeted in the same order as they are initialised in the `<listOfSpecies>`
2) *_histogram.txt*: * for each sample, reports number of each specie at the final time. Each column corresponds to one species, each line to one sample.

### References:
Pleace cite:
Lipkova et al., *S-Leaping: An adaptive, accelerated stochastic simulation algorithm, bridging tau-leaping and R-leaping*, Bulletin of Mathematical Biology 80 (459) (2018)

### Acknowledgement
J. Lipkova, G. Arampatzis, B. Bayati, P. Chatelain, P. Koumoutsakos
