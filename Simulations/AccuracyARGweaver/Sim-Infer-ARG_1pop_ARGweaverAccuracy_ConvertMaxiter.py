# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:36:36 2020

@author: fguo7
"""



import ARGSimulator
import logging
import msprime
import AddRecombNodestoTStables
import gzip
import pandas as pd
import re
import shutil
import dendropy
from dendropy.calculate import treecompare
import sys
import tempfile
import subprocess
import math
import numpy as np
from Bio import SeqIO
import csv
import tskit

debug = False
plot = False

path = "./sim_trees/"
convert_path = './convert_maxiter/'
smc_path = './infer/'

smc2arg_executable = "smc2arg"
ARGweaver_executable = "arg-sample"



" function 'concat_alignments','change_seq_id' are from David's ARGweaver.py"
def concat_alignments(seq_files,file_out):

    for idx, file in enumerate(seq_files):
        seq_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta")) # one line alternative
        if idx == 0:
            concat_seqs = seq_dict
        else:
            for key in concat_seqs:
                concat_seqs[key] += seq_dict[key]            

    # Convert to list and write seq records to fasta file
    concat_records = [concat_seqs[key] for key in concat_seqs]
    SeqIO.write(concat_records, file_out, "fasta")

    
def change_seq_id(seq_file,file_out):
    
    "Because leafnodes are labelled with their numerical ID+1, which affect the conversion to smc"
    with open(seq_file) as input, open(file_out, 'w') as output:
        records = SeqIO.parse(input,'fasta')
        for record in records:
            if debug: print(record.id)
            new_id = int(record.id)
            new_id -= 1
            record.id = str(new_id)


            record.description = ""
            if debug: print(record.id)
            SeqIO.write(record, output, 'fasta')

"modify the original parse_smc of David to parse one smc file"
def parse_smc(smc_file):
    
    "Parse all trees from a smc file"    
    "Set up lists to hold data"

    start_positions = []
    end_positions = []
    newicks = []
        
    "didn't return smc_breaks, but in case there is no infer_ts breaks. The only different is the first element of smc is 1, and the first element of converted ts breaks is 0"
    smc_breaks = []
    smc_breaks.append(1)

    f = open(smc_file)
    line = f.readline()
    while line:
        if 'NAMES' in line:
            "Parse taxon naming scheme"
            fields = re.split(r'\t+', line.strip('\n'))
            fields.pop(0) # remove first
            taxa_nums = list(range(len(fields)))
            taxa_map = {str(key):value for key,value in zip(taxa_nums,fields)} # dict maps taxa numbers back to original taxa labels
        if 'TREE' in line: #if a tree
            fields = re.split(r'\t+', line)
            start_positions.append(int(fields[1]))
            end_positions.append(int(fields[2]))
            smc_breaks.append(int(fields[2]))
            newicks.append(fields[3])
        line = f.readline()
    f.close()
    
    "Convert to Pandas dataframe"
    d = {'StartPosition':start_positions,'EndPosition':end_positions,'Newicks':newicks}
    df = pd.DataFrame(d)
    #df.to_csv(csv_file, sep='\t', index=False)
    
    return df, taxa_map



"""
functions 'CyclicalARGError', 'ARGweaver_smc_to_ts_txts','ARGweaver_arg_to_ts_txts',
'ts_txts_to_trees', 'arg_cmd', 'run_argweaver' are from tskit ts_ARGweaver 
with/without modification
"""
class CyclicalARGError(Exception):
    """
    Exception raised when ARG Weaver generates a cyclical ARG. This is a bug in
    ARGWeaver, so there's nothing we can do about it other than catch the
    error and abort the conversion.

    See https://github.com/mdrasmus/argweaver/issues/19
    """

def ARGweaver_smc_to_ts_txts(smc2bin_executable, prefix, nodes_fh, edges_fh):
    """
    convert the ARGweaver smc representation to tree sequence text format
    """
    logging.debug(
        "== Converting the ARGweaver smc output file '{}' to .arg format using '{}' ==".format(
            prefix + ".smc.gz", smc2bin_executable))
    subprocess.call([smc2bin_executable, prefix + ".smc.gz", prefix + ".arg"])
    with open(prefix + ".arg", "r+") as arg_fh:
        return ARGweaver_arg_to_ts_txts(arg_fh, nodes_fh, edges_fh)

def ARGweaver_arg_to_ts_txts(ARGweaver_arg_filehandle, nodes_fh, edges_fh):
    """
    convert the ARGweaver arg representation to tree sequence tables

    We need to split ARGweaver records that extend over the whole genome into sections
    that cover just that coalescence point.

    returns the mapping of ARGweaver node names to TS node names
    """
    logging.debug("== Converting .arg output to tree seq ==")
    ARG_nodes={} #cr[X] = child1:[left,right], child2:[left,right],... : serves as intermediate ARG storage
    ARG_node_times={} #node_name => time
    node_names={} #map of ARGweaver names -> numbers
    tips = set()
    root_node = None

    #first row gives start and end
    ARGweaver_arg_filehandle.seek(0)
    firstline = next(ARGweaver_arg_filehandle)
    m = re.match(r'^start=(\d+)\s+end=(\d+)\s*$', firstline)
    if m:
        start=float(m.group(1))
        end=float(m.group(2))
    else:
        raise ValueError("Could not find start and end positions in .arg file")

    for line_num, fields in enumerate(csv.DictReader(ARGweaver_arg_filehandle, delimiter='\t')):
        assert (fields['name'] not in ARG_node_times), \
                "duplicate node names identified: line {}".format(line_num)
        #HACK: make sure that parent nodes are strictly older than children.
        #This assumes that parents always have a higher node number
        ARG_node_times[fields['name']] = float(fields['age'])
        #we save info about nodes when looking at their children, so we
        # should save info into parent nodes
        if fields['parents'] == '':
            assert(root_node == None)
            root_node = fields['name']
            #don't need to record anything here, as we will grab details of the
            # root when looking at children
        else:
            if fields['event']=='recomb':
                #each recombination event has multiple parents
                for second_parent, parent in enumerate(fields['parents'].split(",")):
                    if parent not in ARG_nodes:
                        ARG_nodes[parent]={}
                    ARG_nodes[parent][fields['name']]=[
                        (float(fields['pos']) if second_parent else start),
                        (end if second_parent else float(fields['pos']))]
            else:
                #these should all have one parent
                if fields['parents'] not in ARG_nodes:
                    ARG_nodes[fields['parents']]={}
                ARG_nodes[fields['parents']][fields['name']]=[start,end]

                if fields['event']=='gene':
                    #we should trust the labels from 
                    node_names[fields['name']] = int(fields['name'])
                    tips.add(fields['name'])
    #now relabel the internal nodes
    for key in ARG_nodes:
        node_names[key]=len(node_names)

    #recursive hack to make times strictly decreasing, using depth-first topological
    # sorting algorithm
    def set_child_times(node_name, node_order, temporary_marks=set()):
        if node_name in ARG_nodes:
            if node_name in temporary_marks:
                raise CyclicalARGError(
                    "ARG has a cycle in it, around node {}. This should not be possible."
                    "Aborting this conversion!".format(node_name))
            if node_name not in node_order:
                temporary_marks.add(node_name)
                for child_name in ARG_nodes[node_name]:
                    set_child_times(child_name, node_order, temporary_marks)
                node_order.append(node_name)
                temporary_marks.remove(node_name)

    node_order = [] #contains the internal nodes, such that parent is always after child
    set_child_times(root_node, node_order)

    max_epsilon = len(node_order)
    for epsilon, nm in enumerate(node_order):
        ARG_node_times[nm] += 0.001 * (epsilon+1) / max_epsilon

    print("id\tis_sample\ttime", file=nodes_fh)
    for node_name in sorted(node_names, key=node_names.get): #sort by id
        print("{id}\t{is_sample}\t{time}".format(
            id=node_names[node_name],
            is_sample=int(node_name in tips),
            time=ARG_node_times[node_name]),
            file=nodes_fh)

    print("left\tright\tparent\tchild", file=edges_fh)
    for node_name in sorted(ARG_node_times, key=ARG_node_times.get): #sort by time
        # look at the break points for all the child sequences, and break up
        # into that number of records
        try:
            children = ARG_nodes[node_name]
            assert all([ARG_node_times[child] < ARG_node_times[node_name] for child in children])
            breaks = set()
            for leftright in children.values():
                breaks.update(leftright)
            breaks = sorted(breaks)
            for i in range(1,len(breaks)):
                leftbreak = breaks[i-1]
                rightbreak = breaks[i]
                #The read_text function allows `child` to be a comma-separated list of children
                children_str = ",".join(map(str, sorted([
                    node_names[cnode] for cnode, cspan in children.items()
                        if cspan[0]<rightbreak and cspan[1]>leftbreak])))
                print("{left}\t{right}\t{parent}\t{children}".format(
                    left=leftbreak, right=rightbreak, parent=node_names[node_name],
                    children=children_str), file=edges_fh)
        except KeyError:
            #these should all be the tips
            assert node_name in tips, (
                "The node {} is not a parent of any other node, but is not a tip "
                "either".format(node_name))
    nodes_fh.flush()
    nodes_fh.seek(0)
    edges_fh.flush()
    edges_fh.seek(0)
    return node_names


def ts_txts_to_trees(ts_nodes, ts_edges, trees_outname):
    
    "the script modified from tskit ts_ARGweaver"
    
    trees_outname = trees_outname
    logging.info("== Converting new ts ARG to .trees ===")
    try:
        ts = msprime.load_text(nodes=ts_nodes, edges=ts_edges)
    except:
        logging.warning("Can't load the texts file properly. Saved copied to 'bad.nodes' & 'bad.edges' for inspection")
        shutil.copyfile(ts_nodes.name, "bad.nodes")
        shutil.copyfile(ts_edges.name, "bad.edges")
        raise

    if trees_outname:
        ts.dump(trees_outname)
   
    return(ts)

    
def arg_cmd(cmd, stdout=sys.stdout):
    
    with tempfile.TemporaryFile() as stderr:
        exit_status = subprocess.call(cmd, stderr=stderr, stdout=stdout)
        stderr.seek(0)
        if exit_status != 0:
            raise ValueError(
                "Error running '{}': status={}:stderr{}".format(
                    " ".join(cmd), exit_status, stderr.read()))


def run_argweaver(
            seq_file, Ne, recombination_rate, mutation_rate, path_prefix, seed, 
            MSMC_samples, sample_step, burnin_iterations, maxtime, ntimes, verbose=False):
        """
        this produces a whole load of .smc files labelled <path_prefix>i.0.smc,
        <path_prefix>i.10.smc, etc, which we then convert into sts format

        Returns the iteration numbers ('0', '10', '20' etc), the name of the
        statistics file, the total CPU time, and the max memory usage.
        
        Modified a little from the original function of tskit evaluation.
        """
        

        burn_prefix = None
        exe = [ARGweaver_executable, '--fasta', seq_file.name if hasattr(seq_file, "name") else seq_file,
               '--popsize', str(Ne),
               '--recombrate', str(recombination_rate),
               '--mutrate', str(mutation_rate),
               '--ntimes', str(ntimes),
               '--maxtime', str(maxtime),
               # '--compress-seq', str(10),   # To speed up the algorithm, sequence alignment can be compressed during sampling into blocks of sizes, such as 10bp
               '--overwrite']
        if not verbose:
            exe += ['--quiet']
        if seed is not None:
            exe += ['--randseed', str(int(seed))]
        if burnin_iterations > 0:
            burn_in = str(int(burnin_iterations))
            burn_prefix = path_prefix+"_burn"
            logging.info("== Burning in ARGweaver MCMC using {} steps ==".format(burn_in))
            logging.debug("== ARGweaver burnin command is {} ==".format(" ".join(exe)))
            arg_cmd(exe + ['--iters', burn_in,
                                   '--sample-step', burn_in,
                                   '--output', burn_prefix])

            #if burn_in, read from the burn in arg file, rather than the original .sites
            exe += ['--arg', burn_prefix+'.'+ burn_in +'.smc.gz']
        else:
            exe += ['--fasta', seq_file]

        new_prefix = path_prefix + '_' + str(sim) #we append a '_sim' to mark iteration number
        
        # MSMC_sample is number of smc files output? because it output ever-Xth iterations
        iterations = int(sample_step * (MSMC_samples-1))
        exe += ['--output', new_prefix]
        exe += ['--iters', str(iterations)]
        exe += ['--sample-step', str(int(sample_step))]
        logging.info("== Running ARGweaver for {} steps to collect {} samples ==".format( \
            int(iterations), MSMC_samples))
        logging.debug("== ARGweaver command is {} ==".format(" ".join(exe)))
        arg_cmd(exe)

        new_stats_file_name = path_prefix+".stats"

        #concatenate all the stats together
        with open(new_stats_file_name, "w+") as stats:
            if burn_prefix:
                shutil.copyfileobj(open(burn_prefix + ".stats"), stats)
                print("\n", file = stats)
            shutil.copyfileobj(open(new_prefix + ".stats"), stats)
        #To Do: translate these to treesequence objects, so that e.g. edges can be calculated
        #as https://github.com/mdrasmus/argweaver/issues/20 is now closed


"This function is only change a name of David's local_consensus_tree"
def find_local_tree(df, site, taxa_map):
        
    "Get local consensus tree based on all sampled trees for a particular site"
    
    "NOTE: Sites are indexed from 1 here"
    newicks = df['Newicks'][(df['StartPosition'] <= site) & (df['EndPosition'] >= site)] # why only one row meet condition? sometime the breakpoints will meet two row.
    # if plot: print(newicks)
    local_idx = newicks.index[0]
    # if plot: print(local_idx)

    "May want to write this to a file rather than joining as one big string"
    trees_as_str = "".join(newicks)
   
    local_tree = dendropy.Tree.get(data=trees_as_str, schema="newick", rooting="force-rooted")
    
    "relabel the leaf nodes"
    for lf in local_tree.leaf_node_iter():
            # map back the leaf label
            lf.taxon.label = taxa_map[lf.taxon.label]
            
            # accordant with the tskit tree.newick, which
            # add 1 to the leaf node ID, so we can calculate tree distance
            lf.taxon.label = str(1 + int(lf.taxon.label))
    
    # if plot: local_tree.print_plot()
    
    return local_tree, local_idx


def add_sample_state(ts, pop_lst, pops):  #this function can be integrate into dataframe to tables
    tables = ts.dump_tables()
    tables.populations.clear()

    # tskit provides direct access to the columns of each table as numpy arrays.
    # Since columns cannot be modified directly as properties of the tables, they must be extracted, modified, then replaced.
    tnp = tables.nodes.population
    tnp = pop_lst
    tables.nodes.set_columns(flags=tables.nodes.flags, population=tnp, time=tables.nodes.time)
        
    "the population ID must refer to a row in the Population Table"
    for i in range(pops):
        tables.populations.add_row(metadata=bytes(pops))  # must update the populations table in order to change populations in the nodes table

    
    return tables.tree_sequence()

def norm_RF_dist(tree1, tree2):
    "normalized RF distance-OCTAL:Optimal completion of gene trees in polynomial time(p10)"
    ne1 = len(tree1.internal_edges(exclude_seed_edge=True))
    ne2 = len(tree2.internal_edges(exclude_seed_edge=True))
    fpfn = treecompare.false_positives_and_negatives(tree1, tree2)
    print("false_positives_and_negatives",fpfn)
    norm_rf = float((fpfn[0] + fpfn[1])/(ne1 + ne2))
    
    return norm_rf


"Sim params"
n_sims = 100

# sample number in each population
pops = 1
sample_sizes = 50 # should be 1xpops vector
samples_num = 50

Ne = 100.0
genome_length = 1e04 # 10M


# recombination rate
rho = 0.025 # recombination rate per genome per generation
rho_per_site = rho / genome_length


"Store params in dict"
params = {'Ne': Ne, 'rho': rho, 'genome_length': genome_length, 'sample_sizes': sample_sizes}


"Simulate ARG in msprime"    
min_breakpoints = 2
min_seg_length = 20


sim_rho_vals = [2.5e-06*genome_length]
mut_rate = 5.12e-03 # mutation rate per base per generation

for ir in range(len(sim_rho_vals)):
    
    "Get rho params value"
    rho = sim_rho_vals[ir]
    rho_per_site = rho / genome_length
    params['rho'] = rho   
    
    "Performance metrics"
    sims = []
    true_intervals = []
    smc_intervals = []
    RF_dists_smc = []    # RF dist between true tree and smc local tree
    RF_dists_conv = []   # RF dist between true tree and the last converted tree
    norm_RF_dists_smc = []  # normalized RF distances
    norm_RF_dists_conv = [] # normalized RF distances
    
    sim_params = np.zeros([n_sims, 2])
    num_rec = np.zeros([n_sims, 2])

    
    for sim in range(1,n_sims+1): 
        print("-"*20)
        print("Simulation = " + str(sim))
        print("-"*20)
    
        "save the simulation parameters"
        sim_params[sim-1] = sim,rho
        np.savetxt('input parameters_rho' + str(ir) + '.csv', sim_params, delimiter=',')
    
        "Run tree simulation"
        sim_flag = True
        while sim_flag:
            sim_flag = False    
            ts = ARGSimulator.sim_unstructured_ARG(params,min_breakpoints=2,plot=False)
            
            # save the simulated ts nodes tables for population information
            with open(path + str(ir)+ "-sim_ts" + str(sim) + "_nodes.txt", "w+") as sim_node_file:
                print(ts.tables.nodes, file=sim_node_file)
                
            breaks = ts.breakpoints(as_array=True)
            segments = len(breaks) - 1 # number of non-recombinant segments between breakpoints
            for tr_num, tree in enumerate(ts.trees()):
                seq_length = round(tree.interval[1] - tree.interval[0])
                if seq_length < min_seg_length:
                    sim_flag = True
                    
            "Save the simulated tree sequence"
            sim_ts = path + str(ir) + "-sim_ts" + str(sim) + ".trees"
            ts.dump(sim_ts) # save the simulated tree sequence for later comparison
                    
        
        print("Non-recombinant segments = " + str(segments))
        
        "Write local tree and seq files"
        tree_files = [path + str(ir) + "-argweaver_test" + str(sim) + "_tree" + str(i) + ".tre" for i in range(segments)] # simulated ts local trees list
        seq_files = [path + str(ir) + "-argweaver_test" + str(sim) + "_tree" + str(i) + ".fasta" for i in range(segments)]
        prior_tree_ratio = []
        for tr_num, tree in enumerate(ts.trees()):
            seq_length = round(tree.interval[1] - tree.interval[0])
            prior_tree_ratio.append(seq_length / genome_length)
            with open(tree_files[tr_num], "w") as text_file:
                print(tree.newick(), file=text_file)
            ARGSimulator.sim_seqs(tree_files[tr_num],seq_files[tr_num],mut_rate,seq_length)
            
    
    
        "Infer an ARG with ARGweaver"    
        concat_fasta = path + str(ir) + "-argweaver_test" + str(sim) + "_concatenated.fasta"
        new_concat_fasta = path + str(ir) + "-argweaver_test" + str(sim) + "_concatenated_new.fasta"
        concat_alignments(seq_files,concat_fasta)
        change_seq_id(concat_fasta,new_concat_fasta)
        smc_files = smc_path + str(ir) + '-infered_smc' + str(sim) + '/smc-out'
        
        # ARGweaver.run_arg_sampler(new_concat_fasta,smc_files,verbose=False)
        
        # parameters of run_argweaver(seq_file, Ne, recombination_rate, 
        # mutation_rate, path_prefix, iterations, sample_step, burnin_iterations, ntimes)
        # iterations = int(sample_step * (MSMC_samples-1))
        run_argweaver(seq_file=new_concat_fasta, Ne=Ne, recombination_rate=rho_per_site, 
                       mutation_rate=mut_rate, path_prefix=smc_files, seed=20, MSMC_samples=1001, 
                       sample_step=10, burnin_iterations=1000, maxtime=200, ntimes=20, verbose=False)
        
        new_smc_files = smc_files + '_' + str(sim)
    
    
        
        # =============================================================================
        #  This part is converting the ARGweawer infered smc.gz file to TS      
        # =============================================================================
    
        "Convert the ARG from format *.smc.gz to tree sequence--select the MCMC sample which has the maximum joint likelihood"
    
        # read the stats file in each simulation to a dataframe (skip the second line)
        stats_df = pd.read_csv(new_smc_files + '.stats',sep='\t')
        # print(stats_df.head(n=1005))
        # print(stats_df.dtypes)
        
        # keep the lines which iter%10 is an integer, because the smc output every 10 iterations
        for idx, row in stats_df.iterrows():
            if not (row['iter']/10).is_integer():
                stats_df = stats_df.drop(idx)
        
        # find out the iteration has the maximum joint likelihood    
        mit_idx = stats_df['iter'].idxmax()
        mit_iter = stats_df['iter'][mit_idx]
        
        # find out the correspond smc.gz file
        select_smcgz = new_smc_files + '.' + str(mit_iter) + '.smc.gz'# smc_path is the path where the best interation trees in ARGweaver inference
        
        # unzip this smc.gz file
        smcgz_unzip = smc_path + str(ir) + '-infered_smc' + str(sim) + '/mit_' + str(mit_iter) + '.smc' 
        with gzip.open(select_smcgz, 'rt') as f_in, open(smcgz_unzip, "wt") as f_out:
            smc_unzip = f_in.read()
            f_out.write(smc_unzip)
        
        # parse a smc file's local trees 
        tree_df, taxa_map = parse_smc(smcgz_unzip) # I plan to get the best iteration smc.gz file
        
        "Convert the ARG from format *.smc.gz to tree sequence"
     
        convert_file = convert_path + str(ir) + '-convert-out' + str(sim)
        prefix = new_smc_files + '.' + str(mit_iter) # convert the smc.gz file with maximum joint likelihood to ts arg
        
        # save the converted TSnodes, TSedges, and TS trees
        try:
            with open(convert_file+".TSnodes", "w+") as ts_nodes, \
                    open(convert_file+".TSedges", "w+") as ts_edges:
                nodenames = ARGweaver_smc_to_ts_txts(
                    smc2arg_executable, prefix, ts_nodes, ts_edges)
                
                trees_outname = convert_file + ".trees"
                inferred_ts = ts_txts_to_trees(ts_nodes, ts_edges, trees_outname)
                
                # save the node name file which mapping the ARG nodes name and converted TS nodes name
                nodename_file = open(convert_path + str(ir) + '-nodes_names' + str(sim) + '.txt', 'w')
                print(nodenames, file=nodename_file)
                nodename_file.close()
        
        except CyclicalARGError as e:
                logging.warning("Cyclical ARG Exception when converting {}: {}".format(
                    str(sim) + ".msp", e))
        
        if plot:
            for tree in inferred_ts.trees():
                print("-" * 20)
                print("tree {}: interval = {}".format(tree.index, tree.interval))
                print(tree.draw(format="unicode"))
            print(inferred_ts.tables.nodes)
            print()
            print(inferred_ts.tables.edges)
            print()
            
            
        "save number of recombination events"
        true_rec_num = segments - 1   # the simulated number of recombination events    
        infer_rec_num = len(inferred_ts.breakpoints(as_array=True)) - 2
        
        num_rec[sim-1] = true_rec_num,infer_rec_num    
        np.savetxt("number of recombination events_maxIteration" + str(ir) + ".csv", num_rec, delimiter=',')
            
    
        # =============================================================================
        # This part is adding the recombination information to the converted TS files        
        # =============================================================================
    
        tsnodes = convert_file + '.TSnodes'
        tsedges = convert_file + '.TSedges'
        arg = prefix + '.arg'
        nodename_txt = convert_path + str(ir) + '-nodes_names' + str(sim) + '.txt'
        
        
        "Add recombination information based on the recombinate event number"    
        
        "convert the nodes file, edges file, arg, and the nodes name file to dataframe"
        nodes_df,edges_df,arg_df,nodename_df = AddRecombNodestoTStables.readfiles(tsnodes,tsedges,arg,nodename_txt,genome_length)
                
        "read the samples' population information, here it is in the nodes table of simulation tree sequence. In real, it should be given by user"
        samples_pop_df = pd.read_csv(path + str(ir) + "-sim_ts" + str(sim) + "_nodes.txt", sep='\t')
        
        "save path"
        convert_recomb_tree = convert_path + str(ir) + '-convert-recomb' + str(sim) + '.trees'
    
    
        
        "*************If there is no recombination***************"
        if len(tree_df) == 1:
            
            "if no recombination events, only need to sort nodes table, and update nodes table and edges table"
            edges_df_recomb, nodes_df_recomb = AddRecombNodestoTStables.no_recomb_sort_update(edges_df, nodes_df)
               
            
    
    
            "*************If there are recombinations***************"
        else:
        
            "find out the recombination nodes and correspond break point position in the ARG, and map to TS"
            recomb_pos_df = AddRecombNodestoTStables.recomb_pos(arg_df,nodename_df)
            
            "add recombination nodes to nodes dataframe"
            nodes_df_recomb = AddRecombNodestoTStables.add_reomb_nodes(nodes_df,recomb_pos_df)
                                    
            "add recombination information to edges dataframe"
            edges_df_recomb = AddRecombNodestoTStables.add_recomb_edges(edges_df,recomb_pos_df,nodes_df_recomb,genome_length)
            
        "convert the dataframe to tskit tables"
        max_right = max(edges_df_recomb['right'])
        total_length = min(max_right, genome_length) # Aviod generating empty tree(multiroot tree)
        try:
            tables = AddRecombNodestoTStables.df2TreeTables(edges_df_recomb, nodes_df_recomb, total_length)
        except tskit.LibraryError as e:
            print(e)
            continue
        
        "tables to tree sequences"
        inf_recomb_ts = tables.tree_sequence()
            
        
            
        "add sample population information to nodes table"
        pop_info = path + str(ir) + "-sim_ts" + str(sim) + "_nodes.txt"
        sim_node_df = pd.read_csv(pop_info, sep='\s+')
        
        #assign values to the pop_lst
        pop_lst = np.zeros(inf_recomb_ts.num_nodes, dtype='int32') 
        for idx in range(len(pop_lst)):
            if idx < samples_num:
                pop_lst[idx] = sim_node_df['population'][idx]
            else:
                pop_lst[idx] = -1
                    
        #add population information function        
        inf_recomb_ts = add_sample_state(inf_recomb_ts, pop_lst, pops)
        
        "convert the ts nodes table and edges table to dataframes, respectively"
        nodes_df_recomb1, edges_df_recomb1 = AddRecombNodestoTStables.TreeTables2df(inf_recomb_ts)
        
        nodes_df_check, edges_df_check = AddRecombNodestoTStables.ances_length_check(nodes_df_recomb1, edges_df_recomb1)
        
        "dataframes to tables and tables to tree sequences"
        tables_check = AddRecombNodestoTStables.df2TreeTables(edges_df_check, nodes_df_check, total_length)
        inf_recomb_ts_check = tables_check.tree_sequence()
        
        "add population information function"       
        inf_recomb_ts_check = add_sample_state(inf_recomb_ts_check, pop_lst, pops)
        
        
        "save tree sequences"    
        inf_recomb_ts_check.dump(convert_recomb_tree)
            

        
        "plot tree sequences"
        if plot:
            for tree in inf_recomb_ts_check.trees():
                print("-" * 20)
                print("tree {}: interval = {}".format(tree.index, tree.interval))
                print(tree.draw(format="unicode"))
            
            print(inf_recomb_ts_check.tables.nodes)
            print()
            print(inf_recomb_ts_check.tables.edges)
            print()
    
        inf_rec_breaks = inf_recomb_ts_check.breakpoints(as_array=True)
        inf_rec_segments = len(inf_rec_breaks) - 1 # number of non-recombinant segments between breakpoints
        # smc converted ts tree list
        convert_tree_files = [convert_path + str(ir) + "-argweaver_test" + str(sim) + "_convertedtree" + str(i) + ".tre" for i in range(inf_rec_segments)]
        for tr_num, tree in enumerate(inf_recomb_ts_check.trees()):
            try:
                with open(convert_tree_files[tr_num], "w") as text_file:
                    print(tree.newick(), file=text_file)
            except ValueError as e: # because the segment start and end problem, will cause the last segment have multiroot which does't exist 
                print(e)    
                break
    
        
    
    
        # ==================================================================================
        # This part is testing the performances of ARGweaver inference and conversion to TS        
        # ==================================================================================
    
        
        "Compare local trees difference between tree pairs at each gene segment using metric RF distance"
    
        # need to adjust the breaks in true and infer
        for ix, val in enumerate(breaks):
            breaks[ix] = math.floor(val) + 1
        
        # I don't want to patch like this......but have to
        inf_rec_breaks[0] = 1
        
        
        # get all the breaks in the true trees and inferred trees
        all_breaks = np.concatenate((breaks, inf_rec_breaks))
        all_breaks = np.unique(all_breaks)
        all_segments = len(all_breaks)-1
        all_breaks[all_segments] = min(breaks[segments],inf_rec_breaks[inf_rec_segments])
        all_breaks = np.unique(all_breaks)
        print("all breaks array", all_breaks)
        print("true tree breaks", breaks)
        print("smc tree breaks", inf_rec_breaks)
        
        for ix, site in enumerate(all_breaks):
            if ix >= 1:
                # find the one true tree
                try: # find out the first meet conditions
                    true_break_idx = next(idx for idx,val in enumerate(breaks) if val >= site)
    
                except StopIteration as e:
                    print(e)
                    break
                print("site", site)
                
                true_tree_seg = true_break_idx - 1
                print("true local tree", true_tree_seg)
                
                # set up a TaxonNamespace()
                tns = dendropy.TaxonNamespace()
                
                # find the one infer tree
                local_tree, local_idx = find_local_tree(tree_df, site, taxa_map)
                print("smc local tree", local_idx)
                
    
                smc_local_tree = dendropy.Tree(local_tree, taxon_namespace=tns)
                
             
                # find the one true tree
                true_tree = dendropy.Tree.get(file=open(tree_files[true_tree_seg], 'r'), schema="newick", taxon_namespace=tns, rooting="force-rooted")
                
                "Distance between true tree and converted tree"
                dist_smc_true = treecompare.symmetric_difference(true_tree,smc_local_tree)
                print("distance_smc_true", dist_smc_true)
                
                "Normalized RF distances"
                norm_dist_smc_true = norm_RF_dist(true_tree,smc_local_tree)
                print("normalized distance_smc_true", norm_dist_smc_true)
                
                "Distance between true tree and converted tree"
                "only retain the sample leaf nodes---but why can't prune"
                convert_tree = dendropy.Tree.get(file=open(convert_tree_files[local_idx], 'r'), schema="newick", rooting="force-rooted")
                leaf_labels = [str(x) for x in range(1,samples_num + 1)]  # creat a string array ['1','2','3','4','5','6',.....]
                convert_tree.retain_taxa_with_labels(leaf_labels)
                convert_tree = dendropy.Tree(convert_tree, taxon_namespace=tns)
    
                dist_conv_true = treecompare.symmetric_difference(true_tree,convert_tree)
                print("distance_conv_true", dist_conv_true)            

                "Normalized RF distances"
                norm_dist_conv_true = norm_RF_dist(true_tree,convert_tree)
                print("normalized distance_conv_true", norm_dist_conv_true)
    
                
                "For ARGWeaver trees"
                sims.append(sim)
                true_intervals.append(true_tree_seg)
                smc_intervals.append(local_idx)
                RF_dists_smc.append(dist_smc_true) 
                RF_dists_conv.append(dist_conv_true)    
                norm_RF_dists_smc.append(norm_dist_smc_true)
                norm_RF_dists_conv.append(norm_dist_conv_true)
                
 
            
        
    dist_data = {'Simulation': sims,
        'True Local Tree': true_intervals,
        'Infer Local Tree': smc_intervals,
        'RF Dist_smc_true': RF_dists_smc,
        'RF Dist_convert_true': RF_dists_conv,
        'normalized RF Dist_smc_true':norm_RF_dists_smc,
        'normalized RF Dist_convert_true':norm_RF_dists_conv}
    
    print("distance data",dist_data)
            
             
    dist_df = pd.DataFrame(dist_data)
    dist_df.to_csv("test_ARGweaver_results_rho_maxIteration" +str(ir) + ".csv", index=False)
    
    


