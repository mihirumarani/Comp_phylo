# Comp_phylo
To find if we can detect a signature of interspecific competition on a phylogenetic trait data

The effects of interspecific interactions on trait evolution have only recently included in comparative phylogenetic models. It is still unclear whether such models can provide a unique signature of trait distributions at the tip of the phylogeny which is different from heavily employed trait evolution models such as Ornstein-Uhlenbeck models. 

As a primary analysis, I modeled trait evolution along predetermined topologies of phylogenies (using yule process) to incorporate the effects of drift, stabilzing selection (using Ornstein-Uhlenbeck process) and two distinct model of interspecific competition (a. Trait-matching model (Nuismer and Harmon, 2015) and b. Lotka-Volterra dynamic based QGen model (D'Andrea and Barabas, 2018)). Using different parameter combinations related to three mechanisms (drift, OU process, competition), I created datasets of trait values of tip species. 

To study where a statistical test can effectively estimate the parameters related to different evolutionary mechanisms, I employed a Sequencial Monte-Carlo Approximate Bayesian Model (SMC_ABC), where I used the following simplistic model to simulate trait evolution:
∆z= delta + psi*(theta-z) - S*(mu-z)
where, delta is the variance induced due to drift, theta is a trait optima for a stabilizing selection, mu is a mean trait across the species in a community and psi and S are scaling coefficients. 




