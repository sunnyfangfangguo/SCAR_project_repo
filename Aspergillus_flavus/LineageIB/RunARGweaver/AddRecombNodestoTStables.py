# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 17:12:39 2020
add recombination nodes to the converted TSnodes and TSedges files
@author: fguo7
"""
import pandas as pd
import tskit
import numpy as np
import msprime
import ast

debug = True



def readfiles(tsnodes,tsedges,arg,nodename,genome_length):
    "convert the nodes file, edges file, arg, and the nodes name file to dataframe"
    
    "tsnodes is the converted TSnodes file"
    nodes_df = pd.read_csv(tsnodes,sep='\t')
    nodes_df.columns = ['id','flags','time']
    
    "tsedges is the converted TSedges file"
    edges_df = pd.read_csv(tsedges,sep='\t')
    
    "Accord the positions of start and end with tskit tables"
    for edges_key,edges_row in edges_df.iterrows():
        #in some arg infered smc local tree, the right gene location is larger than gene length
        if edges_row['right'] > genome_length:
            edges_df.at[edges_key,'right'] = genome_length
 
    # in some rows there are two children, so explode the column 'child'
    edges_df = edges_df.assign(child=edges_df['child'].str.split(',')).explode('child')
    edges_df = edges_df.reset_index()
    edges_df['id'] = edges_df.index
    
    "arg is the infered arg file"
    arg_df = pd.read_csv(arg,sep='\t',skiprows=1) # read-through the arg file, and skip the first row. 
    
    "nodename is the nodes name file which map the nodes in ARG created by ARGweaver with TSnodes and TSedges"
    # reading and parsig a text file that has the format of a dictionary results ina dictionary data structure in python
    file = open(nodename,"r")
    nodename_contents = file.read()
    nodename_dict = ast.literal_eval(nodename_contents)
    file.close()
    
    # convert a dictionary to dataframe with columns name
    nodename_df = pd.DataFrame.from_dict(nodename_dict,orient='index',columns=['TS_node'])
    nodename_df['ARG_node'] = nodename_df.index
    nodename_df['ARG_node'] = nodename_df['ARG_node'].astype('int64')
    
    return nodes_df,edges_df,arg_df,nodename_df



def recomb_pos(arg_df,nodename_df):
    "find out the recombination nodes and correspond break point position in the ARG, and map to TS"
    # screen out the recombination nodes and correspond break point position 
    recomb_pos_df = arg_df[arg_df['event']=="recomb"]
    
    # covert the arg nodes to ts nodes
    recomb_pos_df = pd.merge(recomb_pos_df,nodename_df,left_on='name',right_on='ARG_node')
    recomb_pos_df = recomb_pos_df[['TS_node','event','pos']]
    
    return recomb_pos_df



def add_reomb_nodes(nodes_df,recomb_pos_df):
    "add recombination nodes to nodes dataframe"
    # add recombination nodes based on the information from recomb_pos_df
    for id_key,id_row in nodes_df.iterrows():
        for recomb_key,recomb_row in recomb_pos_df.iterrows():
            is_equal = id_row['id']-recomb_row['TS_node'] # the indicator variable to judge if the node is a recomb event
            # if debug: print(is_equal)
            if (is_equal==0.0):
                recomb_row1 = {'id':id_row['id']+200000, 'flags':131072,'time':id_row['time']} # the recombination row1
                recomb_row2 = {'id':id_row['id']+500000, 'flags':131072,'time':id_row['time']} # the recombination row2
                
                nodes_df = nodes_df.append(recomb_row1,ignore_index=True)
                nodes_df = nodes_df.append(recomb_row2,ignore_index=True)
                # if debug: print(id_key)
                
                nodes_df = nodes_df.drop(nodes_df[nodes_df['id']==id_row['id']].index) # can't use the id_key as drop index
    
    # sort by "time" and "id"
    nodes_df = nodes_df.sort_values(by=['time','id'],ascending=True).reset_index(drop=True) 
    
    # add a column of new_id, which will be used to update edges table
    nodes_df['new_id'] = nodes_df.index
    nodes_df_recomb = nodes_df
    
    return nodes_df_recomb 



def add_recomb_edges(edges_df,recomb_pos_df,nodes_df_recomb,genome_length):
    "add recombination information to edges dataframe"
    # if parent is a recomb event, then split the parent into two rows
    for edges_key,edges_row in edges_df.iterrows():
        for recomb_key,recomb_row in recomb_pos_df.iterrows():
            is_equal = edges_row['parent'] - recomb_row['TS_node'] # the indicator variable to judge if the edges is a recom event
            # if debug: print(is_equal)
            if (is_equal==0.0):
                # row1 left=left, right=pos, parent = parent + 200000, child keep the same
                # row2 left=pos, right=right, parent = parent + 500000, child keep the same
                recomb_row1 = {'left':edges_row['left'], 'right':recomb_row['pos'], 'parent':edges_row['parent']+200000, 'child':edges_row['child']} # the recombination row1
                recomb_row2 = {'left':recomb_row['pos'], 'right':edges_row['right'], 'parent':edges_row['parent']+500000, 'child':edges_row['child']} # the recombination row2
                
                edges_df = edges_df.append(recomb_row1,ignore_index=True)
                edges_df = edges_df.append(recomb_row2,ignore_index=True)
                # drop the the original row
                edges_df = edges_df.drop(edges_df[edges_df['id']==edges_row['id']].index) # can't use the edges_key as drop index
    
    # if child is a recomb event, then rename the child based on the gene segment positions
    edges_df = edges_df.astype({'child':'int64'})
    for edges_key,edges_row in edges_df.iterrows():
        for recomb_key,recomb_row in recomb_pos_df.iterrows():
            #if debug: print(type(edges_row['child']), type(recomb_row['TS_node']))
            is_equal = edges_row['child'] - recomb_row['TS_node'] # the indicator variable to judge if the edges is a recom event
            #if debug: print(is_equal)
            if (is_equal==0.0):
                child_interval = pd.Interval(edges_row['left'], edges_row['right'])
                recomb1_interval = pd.Interval(0,recomb_row['pos'], closed='left')
                recomb2_interval = pd.Interval(recomb_row['pos'], genome_length, closed='left')
                
                if recomb1_interval.overlaps(child_interval):
                    edges_df.at[edges_key,'child'] = edges_row['child'] + 200000
                elif recomb2_interval.overlaps(child_interval):
                    edges_df.at[edges_key,'child'] = edges_row['child'] + 500000
                    
    # update parent and child in the edges_df based on the new_id in the nodes_df_recomb
    edges_df = pd.merge(edges_df,nodes_df_recomb,left_on='parent',right_on='id')
    edges_df = edges_df[['left','right','parent','child','new_id']]
    edges_df = edges_df.rename(columns={'new_id':'parent_newid'})
    
    edges_df = pd.merge(edges_df,nodes_df_recomb,left_on='child',right_on='id')
    edges_df = edges_df[['left','right','parent','child','parent_newid','new_id']]
    edges_df = edges_df.rename(columns={'new_id':'child_newid'})
    
    # keep the needed columns
    edges_df = edges_df[['left','right','parent_newid','child_newid']]
    edges_df = edges_df.rename(columns={'parent_newid':'parent', 'child_newid':'child'})
    
    edges_df_recomb = edges_df
    
    return edges_df_recomb



def no_recomb_sort_update(edges_df, nodes_df):
    
    "This function used when there is no recombination infered. We only need sort the nodes table and update edges table"
    
    edges_df = edges_df.astype({'child':'int64'})
    # sort by "time" and "id"
    nodes_df = nodes_df.sort_values(by=['time','id'],ascending=True).reset_index(drop=True) 
    
    # add a column of new_id, which will be used to update edges table
    nodes_df['new_id'] = nodes_df.index
    
    # update parent and child in the edges_df based on the new_id in the nodes_df_recomb
    edges_df = pd.merge(edges_df,nodes_df,left_on='parent',right_on='id')
    edges_df = edges_df[['left','right','parent','child','new_id']]
    edges_df = edges_df.rename(columns={'new_id':'parent_newid'})
    
    edges_df = pd.merge(edges_df,nodes_df,left_on='child',right_on='id')
    edges_df = edges_df[['left','right','parent','child','parent_newid','new_id']]
    edges_df = edges_df.rename(columns={'new_id':'child_newid'})
    
    # keep the needed columns
    edges_df = edges_df[['left','right','parent_newid','child_newid']]
    edges_df = edges_df.rename(columns={'parent_newid':'parent', 'child_newid':'child'})
    
    return edges_df, nodes_df



def df2TreeTables(edge_df,node_df,total_length):
    "convert the dataframe to tskit tables"
    
    tables = tskit.TableCollection(total_length)
    for index, row in edge_df.iterrows():
        tables.edges.add_row(row['left'], row['right'], int(row['parent']), int(row['child']))  
    for index, row in node_df.iterrows():
        tables.nodes.add_row(flags=int(row['flags']), time=row['time'])        
    
    tables.sort()
    tables.edges.squash()
    #ts = tables.tree_sequence()
    #print(ts.first().draw(format="unicode"))
    #print()
    
    return tables

def TreeTables2df(ts):
    "nodes table to dataframe"
    # print(ts.tables.nodes)
    col_flags = ts.tables.nodes.flags
    col_pops = ts.tables.nodes.population
    col_time = ts.tables.nodes.time
    
    nodes_dict = {"flags":col_flags,"population":col_pops,"time":col_time}
    nodes_df = pd.DataFrame.from_dict(nodes_dict)
    # nodes_df.to_csv("nodes_table_from_ts.csv")
    
    "edges table to dataframe"
    # print(ts.tables.edges)
    col_left = ts.tables.edges.left
    col_right = ts.tables.edges.right
    col_parent = ts.tables.edges.parent
    col_child = ts.tables.edges.child
    
    edges_dict = {"left":col_left, "right":col_right, "parent":col_parent, "child":col_child}
    edges_df = pd.DataFrame.from_dict(edges_dict)
    
    return nodes_df, edges_df


def ances_length_check(nodes_df, edges_df):
    """
    This function is used to check the ancestral material length of recombination and coalescent nodes.
    Traversal by order of nodes table, from leaf to root
    """
    for id_key,id_row in nodes_df.iterrows():
    
        if id_key < len(nodes_df)-1:
            
            "Recombination node"
            if (id_row['flags']==131072) & (id_row['time']==nodes_df['time'][id_key+1]): # if this is the first node in a recombination event
            
                "child's circumstances"
                child_id1 = edges_df['child'][edges_df['parent']==id_key].values[0]           
                child_id2 = edges_df['child'][edges_df['parent']==(id_key+1)].values[0]
                assert child_id1 == child_id2
                
                "if child is a leaf node, don't need check"      
                "if child is a coalescent or recombination node"
                if nodes_df['flags'][child_id1] != 1:
                    idx_left = edges_df[edges_df['parent']==id_key].index[0]
                    idx_right = edges_df[edges_df['parent']==(id_key+1)].index[0]
                    edges_df.at[idx_left,'left'] = min(edges_df['left'][edges_df['parent']==child_id1].values)
                    edges_df.at[idx_right,'right'] = max(edges_df['right'][edges_df['parent']==child_id1].values)
                    
    
            "Coalescent node"    
            if (id_row['flags']==0):
                
                "it children"
                child_id1 = edges_df['child'][edges_df['parent']==id_key].values[0]
                child_id2 = edges_df['child'][edges_df['parent']==id_key].values[1]
                "if both children are not leaves"
                if (nodes_df['flags'][child_id1] != 1) and (nodes_df['flags'][child_id2] != 1):                             
                    idx_1 = edges_df[edges_df['parent']==id_key].index[0]
                    edges_df.at[idx_1,'left'] = min(edges_df['left'][edges_df['parent']==child_id1].values)
                    edges_df.at[idx_1,'right'] = max(edges_df['right'][edges_df['parent']==child_id1].values)
                    
                    idx_2 = edges_df[edges_df['parent']==id_key].index[1]
                    edges_df.at[idx_2,'left'] = min(edges_df['left'][edges_df['parent']==child_id2].values)
                    edges_df.at[idx_2,'right'] = max(edges_df['right'][edges_df['parent']==child_id2].values)
                
                "if either child is leaf don't need update"    
                if (nodes_df['flags'][child_id1] == 1) and (nodes_df['flags'][child_id2] != 1):             
                    idx_2 = edges_df[edges_df['parent']==id_key].index[1]
                    edges_df.at[idx_2,'left'] = min(edges_df['left'][edges_df['parent']==child_id2].values)
                    edges_df.at[idx_2,'right'] = max(edges_df['right'][edges_df['parent']==child_id2].values)
                    
                if (nodes_df['flags'][child_id1] != 1) and (nodes_df['flags'][child_id2] == 1):                             
                    idx_1 = edges_df[edges_df['parent']==id_key].index[0]
                    edges_df.at[idx_1,'left'] = min(edges_df['left'][edges_df['parent']==child_id1].values)
                    edges_df.at[idx_1,'right'] = max(edges_df['right'][edges_df['parent']==child_id1].values)
    
    
                
        "update to the root"       
        if id_key == len(nodes_df)-1: # only can be coalescent events to root
            child_id1 = edges_df['child'][edges_df['parent']==id_key].values[0]
            child_id2 = edges_df['child'][edges_df['parent']==id_key].values[1]
            
            idx1 = edges_df[edges_df['parent']==id_key].index[0]
            idx2 = edges_df[edges_df['parent']==id_key].index[1]
            
            "if both children are not leaves"
            if (nodes_df['flags'][child_id1] != 1) and (nodes_df['flags'][child_id2] != 1):                    
                edges_df.at[idx1,'left'] = min(edges_df['left'][edges_df['parent']==child_id1].values)        
                edges_df.at[idx1,'right'] = max(edges_df['right'][edges_df['parent']==child_id1].values)
                            
                edges_df.at[idx2,'left'] = min(edges_df['left'][edges_df['parent']==child_id2].values)
                edges_df.at[idx2,'right'] = max(edges_df['right'][edges_df['parent']==child_id2].values)
                
            "if either child is leaf don't need update"    
            if (nodes_df['flags'][child_id1] == 1) and (nodes_df['flags'][child_id2] != 1):                         
                edges_df.at[idx2,'left'] = min(edges_df['left'][edges_df['parent']==child_id2].values)
                edges_df.at[idx2,'right'] = max(edges_df['right'][edges_df['parent']==child_id2].values)
                
            if (nodes_df['flags'][child_id1] != 1) and (nodes_df['flags'][child_id2] == 1):                                         
                edges_df.at[idx1,'left'] = min(edges_df['left'][edges_df['parent']==child_id1].values)
                edges_df.at[idx1,'right'] = max(edges_df['right'][edges_df['parent']==child_id1].values)


    return nodes_df, edges_df



