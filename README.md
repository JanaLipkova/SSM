# Stochastic Simulation Methods 
 * Software for stochastic simulations of well-mixed continuous time Markov chain models for chemical reaction networks.
 * The software takes as input .xml file with the list of chemical reactions to be modelled and simulations parameters (initial locaton, final time,...)
 * The supported methods are:
   * SSA (by Gillespie: "Exact stochastic simulation of coupled chemical reactions" 1977) 
   * Tau-leaping (Cao et. al: "Avoiding negative populations in explicit poisson tau-leaping." 2005)
   * R-leaping   (Auger et al: "R-leaping: accelerating the stochastic simulation algorithm by reaction leaps", 2006)
   * S-leaping   (Lipkova et al: "S-Leaping: An adaptive, accelerated stochastic simulation algorithm, bridging τ-leaping and R-leaping", 2018)
   * Adaptive Tau- and Adatpive S-leaping

* NOTE
  * the solver contains implementation of methods that currently under revision. We will provide more detailed code discription, installation and test example after the the revision process
  * If you wish to use the solver before, just send mail to (jana.lipkova@tum.de) and we will be happy to provide all the details


# Acknowledgement
J. Lipkova, G. Arampatzis, B. Bayati, P. Chatelain, P. Koumoutsakos
