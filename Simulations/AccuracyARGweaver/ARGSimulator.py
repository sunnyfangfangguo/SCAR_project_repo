"""
Created on Mon Feb  1 15:09:01 2021

ARGSimulator simulates (local) phylogenetic trees within ARGs
Trees are simulated under a coalescent model in msprime

@author: david
"""
import msprime
import pyvolve

def sim_unstructured_ARG(params,min_breakpoints=0,plot=False):
    
    "Unpack params from sim"
    haploid_Ne = params['Ne'] / 2 # divide by two because msprime assumes individuals are diploid
    rho = params['rho']
    genome_length = params['genome_length']
    rho_per_site = rho / genome_length
    k = params['sample_sizes']
        
    "Run simulation in msprime"
    breaks = -1
    while breaks < min_breakpoints:
        ts = msprime.simulate(sample_size=k,Ne=haploid_Ne, length=genome_length, recombination_rate=rho_per_site, record_full_arg=True)
        breaks = len(ts.breakpoints(as_array=True)) - 2 # -2 because tskit counts ends as breakpoints
        
    if plot:
        for tree in ts.trees():
            print("-" * 20)
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print(tree.draw(format="unicode"))
        print(ts.tables.nodes)
        print()
        print(ts.tables.edges)
        print()
        print(ts.tables.migrations)
    
    return ts

def sim_ARG(params,min_breakpoints=0,plot=False):
    
    "Unpack params from sim"
    haploid_Ne = params['Ne'] / 2 # divide by two because msprime assumes individuals are diploid
    M = params['M']
    rho = params['rho']
    genome_length = params['genome_length']
    rho_per_site = rho / genome_length
    sample_sizes = params['sample_sizes']
    
    "Set up sample pop configuration"
    pop_config = [msprime.PopulationConfiguration(sample_size=k) for k in sample_sizes]
    
    "Run simulation in msprime"
    breaks = -1
    while breaks < min_breakpoints:
        ts = msprime.simulate(Ne=haploid_Ne, length=genome_length, migration_matrix=M, population_configurations = pop_config, recombination_rate=rho_per_site, record_migrations=True, record_full_arg=True)
        breaks = len(ts.breakpoints(as_array=True)) - 2 # -2 because tskit counts ends as breakpoints
        
    if plot:
        for tree in ts.trees():
            print("-" * 20)
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print(tree.draw(format="unicode"))
        print(ts.tables.nodes)
        print()
        print(ts.tables.edges)
        print()
        print(ts.tables.migrations)
    
    "Just for illustration"
    #tree_0 = ts.at(0)
    #breaks = ts.breakpoints(as_array=True)
    #table_copy = ts.dump_tables()
    #print(breaks)
    
    return ts

def sim_seqs(tree_file,seq_file,mut_rate=1E-08,seq_length=1000):
    
    "User defined params"    
    freqs = [0.25, 0.25, 0.25, 0.25]
    kappa = 2.75
    
    "Read in phylogeny along which Pyvolve should simulate"
    "Scale_tree sets absolute mutation rate"
    my_tree = pyvolve.read_tree(file = tree_file, scale_tree = mut_rate)
    #pyvolve.print_tree(my_tree) # Print the parsed phylogeny
    
    "Or just use an HKY model with kappa"
    nuc_model = pyvolve.Model( "nucleotide", {"kappa":kappa, "state_freqs":freqs})
    
    "Define a Partition object which evolves set # of positions according to my_model"
    my_partition = pyvolve.Partition(models = nuc_model, size = seq_length)
    
    "Define an Evolver instance to evolve a single partition"
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree) 
    
    "Evolve sequences with custom file names"
    #my_evolver(ratefile = "AMR_ratefile.txt", infofile = "AMR_infofile.txt", seqfile = "AMR-seqsim.fasta" )
    my_evolver(seqfile = seq_file)