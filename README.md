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

* NOTE
  * the solver contains implementation of new methods, submitted for a publication. We will provide more detailed code discription, installation and test example after the the paper revision process
  * If you wish to use the solver before, just send mail to (jana.lipkova@tum.de) and we will be happy to provide all the details


# Acknowledgement
J. Lipkova, G. Arampatzis, B. Bayati, P. Chatelain, P. Koumoutsakos
