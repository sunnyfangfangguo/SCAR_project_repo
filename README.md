# SCAR_project_repo

SCAR

SCAR is a recombination-aware phylogeographic inference using the Structured Coalescent with Ancestral Recombination, and using a Bayesian MCMC approach to infer the posterior distribution of demographic parameters. This makes SCAR effciently estimates recombination rate and migration rates from an inferred Ancestral Recombination Graph (ARG). 


Set up

We recommend installing a package manager such as [Anaconda](https://anaconda.org/anaconda/python) to satisfy standard Python library dependencies. [msprime](https://tskit.dev/msprime/docs/stable/intro.html) and [tskit](https://tskit.dev/tskit/docs/stable/) are used to simulate and display ARGs in tree sequence format. You will also need to install the [dendropy](https://dendropy.org/) for working with phylogenetic trees in Python. If you not only run simulation but also want to apply SCAR to the real data, you will also need install [ARGweaver](http://mdrasmus.github.io/argweaver/doc/) to reconstruct ARGs. Once these packages are installed you can run the codes to set up or replicate an analysis.

For a basic overview of how to set up and estimate parameters using SCAR, please see this tutorial.


Source code

The core source code contains three major classes:

SCARLikelihood: a class for calculate parameters' likelihood given the ARG under the structured coalescent with ancestral recombination. The compute_like function can calculate the ARG's likelihood under known/unknown ancestral state by setting suitable dt_step. The functions like_profile, like_profile_rho, and like_profile_M are plot the likelihood profiles of parameters.

MCMC: adopt a Bayesian MCMC approach to infer the posterior distribution of demographic parameters, in particular, a Metropolis-Hastings algorithm is used to sample from the joint posterior distribution of parameters given a fixed ARG. 

ARGSimulator: simulate tskit tree sequence format ARGs, including structured/non-structured. 

Simulations

The simulations include respective or joint estimate parameters from simulated ARGs, test the accuracy of ARGweaver.

EstimatePrameters: under each folder, we have test_*_estimates.py, which is a wrapper code to estimate parameter/parameters. And plot_mcmc_*_ests.py is to plot the estimated results with the true values.

AccuracyARGweaver: we have test the accuracy of ARGweaver under 12 mutation rate / recombination rate ratios, here we only show how to run under one ratio. It is easy to change the values of ratios if you want run others.  


Application to Aspergillus flavus

We use ARGweaver to infer ARGs from the SNP data of lineages IB and IC of A. flavus, and use SCAR to estimate recombination rate and migration rates. Because there are four subpopulations, we change the MCMC code to eatimate more than one migration rate. 





