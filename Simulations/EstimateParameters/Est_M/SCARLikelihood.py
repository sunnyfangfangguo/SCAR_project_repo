"""
Created on Mon Feb 1st 10:03:34 2021

Clean version of struct coalescent w/ anc recombination (SCAR) likelihood function
Ancestral states can be given (known) or marginalized (uknown)

To do:
    - Move simulation routines to ARG Simulation class
    - Move MCMC code to own class (any maybe use PyMC3)
    

@author: david
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import expm

global debug
debug = False
known_ancestral_states = False
dt_step = 0.01 # was 0.1

        
def compute_like(ts,**params):
    
    "Unpack params from dict"
    states = [st.id for st in ts.populations()]
    pops = ts.num_populations
    samples = ts.num_samples
    
    if 'M' in params:
        M = np.transpose(np.array(params['M'])) # transpose to get reverse time matrix
    else:
        print('Need to specify migration rate matrix M')    
    
    if 'Ne' in params:
        Ne = np.array(params['Ne'])
    else:
        print('Need to specify eff pop sizes Ne')
        
    if 'rho' in params:
        rho = params['rho']
    else:
        print('Need to specify recombination rate rho')
        
    
    if 'genome_length' in params:
        genome_length = params['genome_length']
    else:
        print('Need to specify genome length')
    
    
    "Make sure first event (closest to present) in TreeEvents is a sampling event"
    #assert tree_events[0].type == 'sample'
    #assert ts.tables.nodes
    
    #times = np.array(ts.tables.nodes.time,ndmin=1)
    
    "Get transition rate matrix"
    Q = M - np.diag(np.sum(M,axis=1)) # set diagonals to negative row sums
    
    populations = np.array(ts.tables.nodes.population,ndmin=1)
    
    children = np.array(ts.tables.edges.child,ndmin=1)
    parents = np.array(ts.tables.edges.parent,ndmin=1) 
    lefts = np.array(ts.tables.edges.left,ndmin=1)
    rights = np.array(ts.tables.edges.right,ndmin=1)
    
    "Lineage arrays"
    active_lines = []
    line_state_probs = []
    log_like = 0.0 # log likelihood of full tree
    
    "Iterate through each event/node in tree sequence working backwards through time"
    for idx, event in enumerate(ts.tables.nodes):
        
        "Get time of event and time of next event"
        event_time = event.time
        if (idx+1 < len(ts.tables.nodes)): # if not at final event
            next_time = ts.tables.nodes[idx+1].time
        else:
            next_time = event.time
        t_elapsed = next_time - event_time # time elapsed between events
        
        """
            Flags on ts.nodes:
                1 = samples
                0 = coalescent
                131072 = recombination event
                262144 = common ancestor in ARG (but no coalescence in local trees)
        """
        event_type = None
        if event.flags == 1:
            event_type = 'sample'
        if event.flags == 0:
            event_type = 'coalescent'
        if event.flags == 131072:
            event_type = 'recombination'
        if event.flags == 262144:
            event_type = 'hidden_coalescent'
        if event.flags == 524288:
            event_type = 'migration'
        
        "Initialize prob seeing events or no events"
        event_prob = 1.0
        prob_no_coal = 1.0
        prob_no_mig = 1.0
        prob_no_recomb = 1.0
        
        "Update active lineages based on event type: coalescent/sampling/migration events"
        if 'sample' == event_type:
            
            "Add sampled lineage"
            active_lines.append(idx)
            state_probs = np.zeros(pops)
            state_probs[event.population] = 1.0 # set prob to 1.0 for sampled state
            line_state_probs.append(state_probs)            
        
        if 'coalescent' == event_type:
            
            "See if there other events at this time -- appears not to happen"
            #coal_times = times[times == event_time]
            #if len(coal_times) > 1:
                #print("Multiple coalescent events at single time!")
            
            "Get children of parent node at coalescent event"
            coal_children = children[parents == idx] # parent has id == idx in parent column of edges table
            
            "The same parent/child edge may occur more than once in the tree series if not in contiguous local trees"
            coal_children = np.unique(coal_children)
            
            "Make sure coalescent events only occur among two lineages"
            if len(coal_children) > 2:
                print("ERROR: Parent has more than two children at coalescent node")
            assert len(coal_children) == 2
            child1 = coal_children[0]
            child2 = coal_children[1]
            
            child1_idx = active_lines.index(child1)
            child2_idx = active_lines.index(child2)
            
            if debug:
                print(child1, child1_idx)
                print(child2, child2_idx)
                
            "Compute likelihood of coalescent event"
            p1 = line_state_probs[child1_idx]
            p2 = line_state_probs[child2_idx]
            coal_probs = (p1 * p2) / Ne # NOTE: this is element-wise vector multiplication/division
            lambda_sum = sum(coal_probs)
            event_prob = lambda_sum
            
            "Compute new parent state probs"
            if known_ancestral_states:
                parent_probs = np.zeros(pops)
                parent_probs[event.population] = 1.0
            else:
                parent_probs = coal_probs / lambda_sum # renormalize probs
            
            "Update lineage arrays - overwriting child1 with parent"
            active_lines[child1_idx] = idx # name of parent
            line_state_probs[child1_idx] = parent_probs
            del active_lines[child2_idx]
            del line_state_probs[child2_idx]
        
        if 'hidden_coalescent' == event_type:
            
            "Hidden coalescent in ARG not observed in local trees"
            "Need to update active_lines but nothing else"
            
            coal_children = children[parents == idx]
            coal_children = np.unique(coal_children)
            child1 = coal_children[0]
            child2 = coal_children[1]
            child1_idx = active_lines.index(child1)
            child2_idx = active_lines.index(child2)
            
            "Did not have this before for hidden coal events -- compute likelihood of coalescent event"
            p1 = line_state_probs[child1_idx]
            p2 = line_state_probs[child2_idx]
            coal_probs = (p1 * p2) / Ne
            lambda_sum = sum(coal_probs)
            event_prob = lambda_sum
            
            "Compute new parent state probs"
            if known_ancestral_states:
                parent_probs = np.zeros(pops)
                parent_probs[event.population] = 1.0
            else:
                parent_probs = coal_probs / lambda_sum
            
            "Update lineage arrays - overwriting child1 with parent"
            active_lines[child1_idx] = idx # name of parent
            line_state_probs[child1_idx] = parent_probs
            del active_lines[child2_idx]
            del line_state_probs[child2_idx]
        
        if "recombination" == event_type:
            
            """
                At a recombination event: a child node will have two different parents
                We need to find the child shared among these two parents
                Then replace child with left parent and add right parent
            """
            
            "Find child of parent node"
            child = children[parents == idx]
            child = np.unique(child)
            assert len(child) == 1
            
            
            
            "Remember that child may have already been removed"
            if child in active_lines:
            
                "Get indexes of both (left and right) parent of child"
                recomb_parents = parents[children == child]
                
                "Parents edges may occur more than once in the tree series if not in contiguous trees"
                recomb_parents = np.unique(recomb_parents)
                
                "Make sure recombination event results in a child splitting into two parents"
                assert len(recomb_parents) == 2
                
                left_parent = recomb_parents[0]
                right_parent = recomb_parents[1]
    
                child_idx = active_lines.index(child)
                
                "Compute recomination event prob"
                links = max(rights[children==child]) - min(lefts[children==child]) - 1
                event_prob = rho * links / (genome_length - 1)
                
                "Compute new parent state probs"
                if known_ancestral_states:
                    parent_probs = np.zeros(pops)
                    parent_probs[event.population] = 1.0
                else:
                    """
                        Should parent state probs be weighted by prob of recomb event in each pop?
                    """
                    parent_probs = line_state_probs[child_idx]
                
                "Update lineage arrays - overwriting child with left parent"
                active_lines[child_idx] = left_parent # name of parent
                line_state_probs[child_idx] = parent_probs
                
                "Add other parent"
                active_lines.append(right_parent)
                line_state_probs.append(parent_probs)
        
        if 'migration' == event_type:
            
            "Have not yet added migration"
            mig_child = children[parents == idx] # parent has id == idx in parent column of edges table
            "The same parent/child edge may occur more than once in the tree series if not in contiguous local trees"
            mig_child = np.unique(mig_child)
            
            "Migration info from nodes list"
            curr_state = populations[mig_child[0]]
            new_state = populations[idx]
            
            "Migration info from migrations table - should be the same as getting info from populations array above"
            #migration_nodes = np.array(ts.tables.migrations.node,ndmin=1)
            #curr_state_check = ts.tables.migrations.source[migration_nodes == idx][0]
            #new_state_check = ts.tables.migrations.dest[migration_nodes == idx][0]
            #assert curr_state == curr_state_check
            #assert new_state == new_state_check
            
            migrant_idx = active_lines.index(mig_child) #change this for ts index
            
            "Update lineage arrays"
            active_lines[migrant_idx] = idx # name of parent
            
            "Compute event prob"
            if known_ancestral_states:
                new_probs = np.zeros(pops)
                new_probs[new_state] = 1.0 # event. population
                line_state_probs[migrant_idx] = new_probs
                event_prob = M[curr_state][new_state]
            else:
                event_prob = 1.0 # pretend as if we don't see migration events
            
        "Integrate lineage prob equations backwards" 
        
        "Compute prob of no coalescent over time interval"
        if not np.isclose(t_elapsed, 0):
            
            if known_ancestral_states:
                
                A = np.zeros(pops)
                for probs in line_state_probs: A += probs # sum line probs to get total number of lines in each state
                
                "Compute prob of no coalescent over time interval"
                pairs = (A * (A-1)) / 2 # number of pairs in each pop
                lambdas =  pairs * (1/Ne) # coal rate in each pop   
                prob_no_coal = np.exp(-np.sum(lambdas)*t_elapsed)
            
                "Compute prob of no migration over the time interval"
                sam = 0
                for i in range(pops):
                    for z in range(pops):
                        sam += (A[i])*(M[i][z])
                prob_no_mig = np.exp(-sam*t_elapsed)
                
                "Compute prob of no recombination event over the time interval"                
                #calculate links in each population
                pop_links = np.zeros(pops)
                
                for i in range(len(active_lines)):
                    #calculate links of each lineage
                    line = active_lines[i]
                    
                    if len(children[children==line])==1:
                        line_links = (rights[children==line] - lefts[children==line])[0] - 1
                    
                    elif len(children[children==line])>=2:
                        line_links = max(rights[children==line]) - min(lefts[children==line]) - 1 
                    
                    
                    pop_links += line_links * line_state_probs[i]
                
                # sum links in all populations
                sum_links = np.sum(pop_links)
                          
                prob_no_recomb = np.exp(-sum_links * rho / (genome_length - 1) * t_elapsed) # assumes rho / genome_length is constant across pops
                
            else:
                
                "Should move everything below into seperate function because this will change depending on approximation used"
            
                "Integrate lineage prob equations backwards"
                dt_times = list(np.arange(event_time,next_time,dt_step)) # integration steps going backwards in time
                for idx,tx in enumerate(dt_times):
                    
                    "Fix this so index does not go out of bounds"
                    if (idx+1 < len(dt_times)):
                        dt = dt_times[idx+1] - tx # integration time step
                    else:
                        dt = next_time - tx
                    if debug:
                        print("dttimes",tx,dt,event_time,next_time)

                    "Should not need to exponentiate transition matrix if dt is small enough"
                    expQdt = expm(Q*dt) # exponentiate time-scaled transition rate matrix

                    "Update line state probs using Euler integration"
                    for ldx,probs in enumerate(line_state_probs):
                        line_state_probs[ldx] = np.matmul(probs,expQdt)
                    
                    A = np.zeros(pops)
                    for probs in line_state_probs: A += probs # sum line probs to get total number of lines in each state
                    
                    "Compute prob of no coalescent over time interval"
                    pairs = (A * (A-1)) / 2 # number of pairs in each pop
                    pairs = pairs.clip(min=0) # make sure non are negative
                    lambdas = pairs * (1/Ne) # coal rate in each pop
                    prob_no_coal *= np.exp(-np.sum(lambdas)*dt)
                    
                    "Compute prob of no migration over the time interal"
                    prob_no_mig = 1.0
                                        
                    "Compute prob of no recombination event"
                    #calculate links in each population
                    pop_links = np.zeros(pops)
                    
                    for i in range(len(active_lines)):
                        #calculate links of each lineage
                        line = active_lines[i]
                        
                        if len(children[children==line])==1:
                            line_links = (rights[children==line] - lefts[children==line])[0] - 1
                        
                        elif len(children[children==line])>=2:
                            line_links = max(rights[children==line]) - min(lefts[children==line]) - 1 
                        
                        
                        pop_links += line_links * line_state_probs[i]
                    
                    # sum links in all populations
                    sum_links = np.sum(pop_links)
                                        
                                                     
                    prob_no_recomb *= np.exp(-sum_links * rho / (genome_length - 1) * dt)
        
        #print(event_prob)
        #print(prob_no_coal)
        #print(prob_no_mig)
        #print(line_state_probs)
        log_like += np.log(event_prob) + np.log(prob_no_coal) + np.log(prob_no_mig) + np.log(prob_no_recomb)
        
    return log_like


def like_profile(ts,params,true_val):
    
    profile_vals = list(np.arange(50,400,10)) # for true val = 200
    like_vals = []
    for val in profile_vals:
        params['Ne'] = val
        like_vals.append(compute_like(ts,**params))
    sns.set(style="darkgrid")
    fig, axs = plt.subplots(1, 1)
    sns.lineplot(x=profile_vals, y = like_vals, ax=axs)
    axs.plot([true_val,true_val],axs.get_ylim(), 'darkred')
    axs.set_xlabel('Ne', fontsize=14)
    axs.set_ylabel('Likelihood', fontsize=14)
    
def like_profile_rho(ts,params,true_val):
    
    #profile_vals = np.logspace(-2.0, 0.0, base = 10, num=20)
    profile_vals = list(np.arange(0.01,2.0,0.02))
    like_vals = []
    for val in profile_vals:
        params['rho'] = val
        like_vals.append(compute_like(ts,**params))
    sns.set(style="darkgrid")
    fig, axs = plt.subplots(1, 1)
    sns.lineplot(x=profile_vals, y = like_vals, ax=axs)
    axs.plot([true_val,true_val],axs.get_ylim(), 'darkred')
    axs.set_xlabel('Recombination rate', fontsize=14)
    axs.set_ylabel('Likelihood', fontsize=14)
    
def like_profile_M(ts,params,M):
    
    true_val = M[0][1]
    profile_vals = list(np.arange(0.01,0.5,0.01))
    like_vals = []
    for val in profile_vals:
        M = [[0.0,val,val], [val,0.0,val], [val,val,0.0]]
        params['M'] = M
        like_vals.append(compute_like(ts,**params))
    sns.set(style="darkgrid")
    fig, axs = plt.subplots(1, 1)
    sns.lineplot(x=profile_vals, y = like_vals, ax=axs)
    axs.plot([true_val,true_val],axs.get_ylim(), 'darkred')
    axs.set_xlabel('Migration rate', fontsize=14)
    axs.set_ylabel('Likelihood', fontsize=14)
    fig.set_size_inches(6, 6)
    fig.savefig('coalRecombMig_MLE_migRate_test.png', dpi=200)

if __name__ == '__main__':
       
    import ARGSimulator
    
    "Specify sim params"
    samples = 50, 50
    genome_length = 1e4
    rho = 0.5 #1e-1
    Ne = 1.0  # effective pop sizes
    sample_sizes = np.array(samples) # should be 1xpops vector
    
    "Migration and configuration of samples"
    val = 0.25
    M = [[0.0,val], [val,0.0]]  # migration rate matrix
    # M = [[0]]
    
    "Store params in dict"
    params = {'Ne': Ne, 'rho': rho, 'M': M, 'genome_length': genome_length, 'sample_sizes':  sample_sizes}
    
    "Run sim(s)"
    ts = ARGSimulator.sim_ARG(params,min_breakpoints=1,plot=False)
    
    breaks = len(ts.breakpoints(as_array=True)) - 2 # minus 2 b/c tskit counts ends as breakpoints
    print('recombination breakpoints = ' + str(breaks))
    
    #mig_events = ts.num_migrations
    #print('Num migrations: ', str(ts.num_migrations))
    
    "Check likelihood is valid"
    L = compute_like(ts,**params)
    
    #like_profile(ts,params,Ne)
    like_profile_rho(ts,params,rho)
    #like_profile_M(ts,params,M)

