"""
Created on Mon Feb  1 15:17:13 2021

Test inferences of Ne under SCAR model w/ no migration/recombination

@author: david
"""
import MCMC
import ARGSimulator
import numpy as np

def expected_tree_length(sample_size,Ne):
    
    """
        sample_size = sample size
        Ne = haploid pop size
    """
    
    "Expected total tree length"
    tree_length = 0.0
    for k in range(2,sample_size+1):
        tree_length += k * Ne / (k*(k-1)/2)
    return tree_length


"Specify sim params"
sims = 100
samples = 100
genome_length = 1e4
rho = 0.0 #1e-1
Ne = 1.0  # effective pop sizes
mig_rate = 0.25
M = [[0.0,mig_rate], [mig_rate,0.0]]  # migration rate matrix
sample_sizes = [samples]*2 # should be 1xpops vector

T_length = expected_tree_length(samples,Ne)
expected_mig_events = T_length * mig_rate
print('Expected migration events = ',f"{expected_mig_events:.2f}")

"Store params in dict"
params = {'Ne': Ne, 'rho': rho, 'M': M, 'genome_length': genome_length, 'sample_sizes':  sample_sizes}

"Pass estimated params as dict"
est_params = {'Ne': False, 'rho': True, 'M': True}

true_fit_vals = np.zeros([sims,2])
est_fit_vals = np.zeros([sims,6])
#sim_Ne_vals = np.logspace(-1.0, 3.0, base = 10, num=sims) # e.g. from 10e-1 to 1000
sim_rho_vals = np.linspace(0.1,10.0,num=sims)
sim_mig_vals = np.linspace(0.1,2.0,num=sims)

for s in range(sims):
        
    print('Sim: ', str(s))
    mig_rate = sim_mig_vals[s]
    rho = sim_rho_vals[s]
    params['M'] = [[0.0,mig_rate], [mig_rate,0.0]]
    params['rho'] = rho
    true_fit_vals[s] = rho,mig_rate
    
    "Run tree simulation"
    mig_events = 0
    while mig_events < 1:
        ts = ARGSimulator.sim_ARG(params,min_breakpoints=0,plot=False)
        breaks = len(ts.breakpoints(as_array=True))
        mig_events = ts.num_migrations
    #ts_sim_file = 'test_SCAR_mig' + f"{mig_rate:.2f}" + '_sim' + str(s)
    #ts.dump(ts_sim_file, zlib_compression=False)
    
    print('Max root time: ', str(ts.max_root_time))
    print('Num migrations: ', str(ts.num_migrations))
    
    "Run MCMC inference"
    plot_info = {'sim':str(s), 'true_val_rho':rho,'true_val_mig':mig_rate}
    opt_iters = []
    opt_iters = [200,500,1000,2000,5000]
    estimates = MCMC.sample(params,est_params,ts,iterations=10000,opt_iters=opt_iters,burnin=1000,plot=True,**plot_info)
    est_fit_vals[s,:] = estimates[:,0]

    np.savetxt('test_SCAR_true_rho_mig_unknownStates.csv', true_fit_vals, delimiter=',')
    np.savetxt('test_SCAR_mcmc_rho_mig_unknownStates_ests.csv', est_fit_vals, delimiter=',')
            