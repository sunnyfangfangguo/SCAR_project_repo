"""
Created on Mon Feb  1 15:17:13 2021

Test inferences of Ne under SCAR model w/ no migration/recombination

@author: david
"""
import MCMC_flattenM
import numpy as np
import tskit

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


"Specify params"
samples = 11,3,11,26
genome_length = 5033688
rho = 6.92e-3 #region recombination rate per generation
Ne = 1554  # effective pop sizes

"initial migrate rates"
mig_rate01 = 0.1
mig_rate02 = 0.1
mig_rate03 = 0.1
mig_rate12 = 0.1
mig_rate13 = 0.1
mig_rate23 = 0.1

"migration matrix M"
M = [[0.0,mig_rate01,mig_rate02,mig_rate03], 
     [mig_rate01,0.0,mig_rate12,mig_rate13], 
     [mig_rate02,mig_rate12,0.0,mig_rate23],
     [mig_rate03,mig_rate13,mig_rate23,0.0]]  # migration rate matrix


sample_sizes = np.array(samples) # should be 1xpops vector

T_length = expected_tree_length(sum(samples),Ne)
expected_mig_events = T_length * (mig_rate01 + mig_rate02 + mig_rate03 + mig_rate12 + mig_rate13 + mig_rate23)
print('Expected migration events = ',f"{expected_mig_events:.2f}")

"Store params in dict"
params = {'Ne': Ne, 'rho': rho, 'M': M, 'genome_length': genome_length, 'sample_sizes':  sample_sizes}

"Pass estimated params as dict"
est_params = {'Ne': False, 'rho': True, 'M': True}

"estimation params number"
n_params = 0
if est_params['Ne']: n_params += 1
if est_params['rho']: n_params += 1
if est_params['M']: 
    M_Len = np.array(M).shape[0]
    n_M = int(M_Len * (M_Len - 1) / 2) # number of elements in upper triangle without diagonal
    n_params += n_M

true_fit_vals = np.zeros([1,n_params])
est_fit_vals = np.zeros([1,3*n_params])

true_fit_vals = rho,mig_rate01,mig_rate02,mig_rate03,mig_rate12,mig_rate13,mig_rate23

"load ts file"
ts = tskit.load("../../convert/convert-recomb_aflavus.trees")
breaks = len(ts.breakpoints(as_array=True))
mig_events = ts.num_migrations    
print('Max root time: ', str(ts.max_root_time))
print('Num migrations: ', str(ts.num_migrations))
    
"Run MCMC inference"
plot_info = {'true_val_rho':rho, 'true_val_M0':mig_rate01, 'true_val_M1':mig_rate02, 'true_val_M2':mig_rate03,
             'true_val_M3':mig_rate12, 'true_val_M4':mig_rate13, 'true_val_M5':mig_rate23}
opt_iters = []
opt_iters = [200,500,1000,3000,5000,8000,10000,12000,14000,16000,18000,20000,25000,30000,34000,38000]
estimates,theta_samples,p_samples,prior_arr = MCMC_flattenM.sample(params,est_params,ts,iterations=40000,opt_iters=opt_iters,burnin=1000,plot=True,**plot_info)
est_fit_vals = estimates[:,0]

np.savetxt('Aspergillus_SCAR_true_rho_mig_unknownStates.csv', true_fit_vals, delimiter=',')
np.savetxt('Aspergillus_SCAR_mcmc_rho_mig_unknownStates_ests.csv', est_fit_vals, delimiter=',')
np.savetxt('Aspergillus_SCAR_mcmc_log_theta_samples.csv', theta_samples, delimiter=',')
np.savetxt('Aspergillus_SCAR_mcmc_log_p_samples.csv', p_samples, delimiter=',')
np.savetxt('Aspergillus_SCAR_mcmc_log_priors.csv', prior_arr, delimiter=',')             
