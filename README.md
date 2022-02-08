# SCAR_project_repo

## SCAR

SCAR allows for recombination-aware phylogeographic inference using the Structured Coalescent with Ancestral Recombination. SCAR uses a Bayesian MCMC approach to infer demographic parameters including migration and recombination rates from a reconstructed Ancestral Recombination Graph (ARG). 


## Set up

We recommend installing a package manager such as [Anaconda](https://anaconda.org/anaconda/python) to satisfy standard Python library dependencies. [msprime](https://tskit.dev/msprime/docs/stable/intro.html) and [tskit](https://tskit.dev/tskit/docs/stable/) are used to simulate and display ARGs in tree sequence format. You will also need to install [dendropy](https://dendropy.org/) for working with phylogenetic trees in Python. If you would like to reconstruct ARGs from sequence data, you will also need to install [ARGweaver](http://mdrasmus.github.io/argweaver/doc/).

For a basic overview of how to set up and estimate parameters using SCAR, please see this [tutorial](./SCAR_tutorial.ipynb).


## Source code

The core source code contains three major classes:

SCARLikelihood: a class for computing the likelihood of an ARG under the structured coalescent with ancestral recombination. The compute_like function calculates the ARG's likelihood under known/unknown ancestral state given a set of demographic parameters. If ancestral states are unknown, a suitable dt_step param is also required to calculate lineage state probabilities by numerical integration. The functions like_profile, like_profile_rho, and like_profile_M can be used to plot likelihood profiles for these parameters.

MCMC: uses a Bayesian MCMC approach to infer the posterior distribution of demographic parameters, in particular, a Metropolis-Hastings algorithm is used to sample from the joint posterior distribution of parameters given a fixed ARG. 

ARGSimulator: simulate ARGs in tskit tree sequence format with population structure. 

## Simulations

This code is provided to test the accuracy of ARGweaver in reconstructing simulated ARGs and to test parameter inference from ARGs using SCAR.

EstimateParameters: under each folder, we have test_*_estimates.py, which is a wrapper code to estimate parameters. plot_mcmc_*_ests.py can be used to plot the estimated parameters alongside their true values.

AccuracyARGweaver: tests the accuracy of ARGweaver under 12 mutation rate / recombination rate ratios.  


## Application to Aspergillus flavus

Provides code to replicate our analysis of the <em>Aspergillus flavus</em> data. We use ARGweaver to infer ARGs from SNP data for lineages IB and IC of <em>A. flavus</em>, and use SCAR to estimate recombination rate and migration rates.




