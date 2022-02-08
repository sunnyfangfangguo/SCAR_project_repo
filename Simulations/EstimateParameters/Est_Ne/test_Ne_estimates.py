"""
Created on Mon Feb  1 15:17:13 2021

Test inferences of Ne under SCAR model w/ no migration/recombination

@author: david
"""
import MCMC
import ARGSimulator
import numpy as np
    
"Specify sim params"
sims = 100
samples = 100
genome_length = 1e4
rho = 0 #1e-1
Ne = 1.0  # effective pop sizes
M = [[0.0]]
sample_sizes = [samples]*1 # should be 1xpops vector

"Store params in dict"
params = {'Ne': Ne, 'rho': rho, 'M': M, 'genome_length': genome_length, 'sample_sizes':  sample_sizes}

"Pass estimated params as dict"
est_params = {'Ne': True, 'rho': False, 'M': False}

true_fit_vals = np.zeros([sims,1])
est_fit_vals = np.zeros([sims,3])
#sim_Ne_vals = np.logspace(-1.0, 3.0, base = 10, num=sims) # e.g. from 10e-1 to 1000
sim_Ne_vals = np.linspace(0.1,1000,num=sims)
for s in range(sims):
        
    print('Sim: ', str(s))
    Ne = sim_Ne_vals[s]
    params['Ne'] = Ne
    true_fit_vals[s] = Ne
    
    "Run tree simulation"
    ts = ARGSimulator.sim_ARG(params,min_breakpoints=0,plot=False)
    breaks = len(ts.breakpoints(as_array=True))
    ts_sim_file = 'test_SCAR_Ne' + f"{Ne:.2f}" + '_sim' + str(s) + '.trees'
    ts.dump(ts_sim_file, zlib_compression=False)
    
    print('Max root time: ', str(ts.max_root_time))
    
    "Run MCMC inference"
    plot_info = {'sim':str(s), 'true_val':Ne}
    opt_iters = []
    opt_iters = [200,500,1000,2000,5000]
    estimates = MCMC.sample(params,est_params,ts_sim_file,iterations=10000,opt_iters=opt_iters,burnin=1000,plot=True,**plot_info)
    est_fit_vals[s,:] = estimates[:,0]

    np.savetxt('test_SCAR_linearScale_true_Ne.csv', true_fit_vals, delimiter=',')
    np.savetxt('test_SCAR_linearScale_mcmc_Ne_ests.csv', est_fit_vals, delimiter=',')
            