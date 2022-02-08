# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:36:36 2020

CyclicalARGError, ARGweaver_smc_to_ts_txts, ARGweaver_arg_to_ts_txts,
ts_txts_to_trees(a liitle modification) are from tskit :ts_ARGweaver.py.

Modify the original parse_smc from David to parse one smc file

@author: fguo7
"""


import logging
import msprime
import AddRecombNodestoTStables
import gzip
import pandas as pd
import re
import shutil
import subprocess
import numpy as np
import csv





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


def parse_smc(smc_file):
    ""
    
    "Parse all trees from smc file"
    
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







"paramterss"
plot = True
genome_length = 5033688
sample_number = 51
pops = 4

"paths"
smc2arg_executable = "smc2arg"   # the executable path of smc2arg

path = "./aflavus/"
sample_pop_file = path + "sample_name_population.csv"

smc_path = './infer/'
new_smc_files = smc_path + 'smc-out_aflavus'

convert_path = './convert/' 
convert_file = convert_path + 'convert-out_aflavus'



# =============================================================================
#  This part is converting the ARGweawer infered smc.gz file to TS      
# =============================================================================

"Convert the ARG from format *.smc.gz to tree sequence--select the MCMC sample which has the maximum joint likelihood"

"read the stats file to a dataframe"
stats_df = pd.read_csv(new_smc_files + '.stats',sep='\t')
# print(stats_df.head(n=1005))
# print(stats_df.dtypes)

"if like choosing an iteration from first n_iters--in case the iteration need a very long time to run and reach a convergence before the iterations"
# n_iters = 50
# stats_df = stats_df.head(n_iters)

"keep the lines which iter%10 is an integer, because the smc output every 10 iterations"
for idx, row in stats_df.iterrows():
    if not (row['iter']/10).is_integer():
        stats_df = stats_df.drop(idx)

"find out the iteration has the maximum joint likelihood"   
mjlh_idx = stats_df['joint'].idxmax()
mjlh_iter = stats_df['iter'][mjlh_idx]


"find out the correspond smc.gz file"
select_smcgz = new_smc_files + '.' + str(mjlh_iter) + '.smc.gz'# smc_path is the path where the best interation trees in ARGweaver inference

"unzip this smc.gz file"
smcgz_unzip = smc_path + 'mjlh_' + str(mjlh_iter) + '.smc' 
with gzip.open(select_smcgz, 'rt') as f_in, open(smcgz_unzip, "wt") as f_out:
    smc_unzip = f_in.read()
    f_out.write(smc_unzip)

"parse a smc file's local trees"
tree_df, taxa_map = parse_smc(smcgz_unzip) # I plan to get the best iteration smc.gz file

"Convert the ARG from format *.smc.gz to tree sequence"
prefix = new_smc_files + '.' + str(mjlh_iter) # convert the smc.gz file with maximum joint likelihood to ts arg

"save the converted TSnodes, TSedges, and TS trees"
try:
    with open(convert_file+".TSnodes", "w+") as ts_nodes, \
            open(convert_file+".TSedges", "w+") as ts_edges:
        nodenames = ARGweaver_smc_to_ts_txts(
            smc2arg_executable, prefix, ts_nodes, ts_edges)
        
        trees_outname = convert_file + ".trees"
        inferred_ts = ts_txts_to_trees(ts_nodes, ts_edges, trees_outname)
        
        # save the node name file which mapping the ARG nodes name and converted TS nodes name
        nodename_file = open(convert_path + 'ARG_TS_nodes_names_aflavus.txt', 'w')
        print(nodenames, file=nodename_file)
        nodename_file.close()

except CyclicalARGError as e:
        logging.warning("Cyclical ARG Exception when converting {}: {}".format(
            "aflavus.msp", e))

if plot:
    for tree in inferred_ts.trees():
        print("-" * 20)
        print("tree {}: interval = {}".format(tree.index, tree.interval))
        # print(tree.draw(format="unicode"))
    print(inferred_ts.tables.nodes)
    print()
    print(inferred_ts.tables.edges)
    print()
    
    

    

# =============================================================================
# This part is adding the recombination information to the converted TS files        
# =============================================================================

tsnodes = convert_file + '.TSnodes'
tsedges = convert_file + '.TSedges'
arg = prefix + '.arg'
nodename_txt = convert_path + 'ARG_TS_nodes_names_aflavus.txt'


"Add recombination information based on the recombinate event number"

"convert the nodes file, edges file, arg, and the nodes name file to dataframe"
nodes_df,edges_df,arg_df,nodename_df = AddRecombNodestoTStables.readfiles(tsnodes,tsedges,arg,nodename_txt,genome_length)

"read the samples' population information. In real, it should be given by user"
samples_pop_df = pd.read_csv(sample_pop_file)

"save path"
convert_recomb_tree = convert_path + 'convert-recomb_aflavus.trees'



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
tables = AddRecombNodestoTStables.df2TreeTables(edges_df_recomb, nodes_df_recomb, total_length)

"tables to tree sequences"
inf_recomb_ts = tables.tree_sequence()
    

    
"add sample population information to nodes table"
#assign values to the pop_lst
pop_lst = np.zeros(inf_recomb_ts.num_nodes, dtype='int32') 
for idx in range(len(pop_lst)):
    if idx < sample_number:  # only leaf node has population information
        pop_lst[idx] = samples_pop_df['population'][idx]
    else:
        pop_lst[idx] = -1 # internal node don't know the population information

"add population information function"       
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
        # print(tree.draw(format="unicode"))
    
    print(inf_recomb_ts_check.tables.nodes)
    print()
    print(inf_recomb_ts_check.tables.edges)
    print()


print("Running is OK")