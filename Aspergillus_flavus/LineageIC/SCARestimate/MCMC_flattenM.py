"""
Created on Mon Feb  1 13:20:27 2021

MCMC sampling from posterior. Only Metropolis-Hastings sampling is implemented

To do:
    -Allow non-uniform priors to be given

@author: david
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import SCARLikelihood
from scipy.stats import expon

debug = False
def theta2params(theta,params,param_map):
    M = np.array(params['M'])
    mig_rate = M[np.triu_indices(M.shape[0], k=1)]
    migrateL = len(mig_rate) # size of flaten numpy array
    
    "Update params dict with values in theta"
    for k,v in param_map.items():
        if 'M' in k:
            for i in range(migrateL):
                if k == ('M'+str(i)):
                    mig_rate[i] = theta[v]
        else:
            params[k] = theta[v]
    
    "convert the upper triangel elements without diagnal to matrix"            
    pops = int(M.shape[0]) # calculate the length of matrix M
    M = np.zeros((pops,pops))
    M[np.triu_indices(M.shape[0], k=1)] = mig_rate
    M = M + M.T
    params['M'] = M

    if debug: print('params', params)
        
    return params

def params2theta(params,theta,param_map):
    
    M = np.array(params['M'])
    mig_rate = M[np.triu_indices(M.shape[0], k=1)]
    migrateL = len(mig_rate) # size of flaten numpy array

    "Flatten params dict to theta vec"
    for k,v in param_map.items():
        
        if 'M' in k:
            for i in range(migrateL):
                if k == ('M'+str(i)):
                    theta[v] = mig_rate[i]
                    
        else:
            theta[v] = params[k]
        
    return theta
    
    
def sample(params,est_params,ts,iterations=1000,opt_iters=[],burnin=0,plot=False,**plot_info):
    
    M = np.array(params['M'])
    mig_rate = M[np.triu_indices(M.shape[0], k=1)] # flatten the upper triangel of matrix M and without diagnal elements
    migrateL = len(mig_rate) # size of flaten numpy array
    
    param_map = {}
    param_map_index = 0
    
    "If estimating Ne"
    if est_params['Ne']:
        param_map['Ne'] = param_map_index
        param_map_index += 1
    
    "If estimating rho -- recombination rate"
    if est_params['rho']:
        param_map['rho'] = param_map_index
        param_map_index += 1
     
    "If estimating M -- migration rate"
    if est_params['M']:
        for i in range(migrateL):
            param_map['M'+str(i)] = param_map_index
            param_map_index += 1
    
    
    "Set up parameters and samples of theta"
    theta_now = np.array([0.0]*param_map_index)
    #theta_now = np.array([params[k] for k in param_map]) # current param value
    theta_now = params2theta(params,theta_now,param_map)
    dim = theta_now.size
    theta_samples = np.zeros((dim,iterations))
    theta_samples[:,0] = theta_now
    
    """
        Specify the mean/variance of our lognormal prior
        Here theta has 7 elements [rho, mig_rate01, mig_rate02, mig_rate03, mig_rate12, mig_rate13, mig_rate23]        
    """
    
    "prior probability for rho and migration rate (m)"
    prior_rho_mu = 1.11e-03 # in exponential distribution, mu = 1 / lambda
    prior_mig_mu = 0.1  # in exponential distribution, mu = 1 / lambda
       
    "prior_now is product of all the parameter prior"
    prior_now = 0.0
    for d_i in range(dim):
        if d_i == 0: # for parameter rho
            prior_now += np.log(expon.pdf(theta_now[d_i],loc = 0,scale = prior_rho_mu))

        else: # for parameter migration rates
            prior_now += np.log(expon.pdf(theta_now[d_i],loc = 0,scale = prior_mig_mu))

    "p_now is likelihood times the prior_now before log" 
    p_now = SCARLikelihood.compute_like(ts,**params) + prior_now # current posterior prob
    p_samples = np.zeros(iterations)
    prior_arr = np.zeros(iterations)

    "Proposal denstiy"
    cov = np.zeros((dim, dim), float)
    np.fill_diagonal(cov, 0.01)
    
    accept = 0 # counter for number of accepted proposals
    
    "Main MCMC loop"
    for i in range(1,iterations): # was range(2,iterations)
        
        if i % 100 == 0:
            print("MCMC iterations " + str(i))
            
        if i in opt_iters:
            "Optimize proposal covariance matrix based on current samples"
            if dim > 1:
                cov = np.cov(theta_samples, rowvar=True)
            else:
                cov[0,0] = np.var(theta_samples)
        
        "Propose new params theta_new"
        theta_new = np.random.multivariate_normal(theta_now, cov)
        while any(theta_new[theta_new<0]):
            theta_new = np.random.multivariate_normal(theta_now, cov)
        #theta_new = theta_new.clip(min=0) # don't allow negative value for rates    
        
        #for k,v in param_map.items():
            #params[k] = theta_new[v]
        params = theta2params(theta_new,params,param_map)
        
        "prior_now is product of all the parameter prior before log"
        prior_new = 0.0
        for d_i in range(dim):
            if d_i == 0: # for parameter rho
                prior_new += np.log(expon.pdf(theta_new[d_i],loc = 0,scale = prior_rho_mu))
                
            else: # for parameter migration rates
                prior_new += np.log(expon.pdf(theta_new[d_i],loc = 0,scale = prior_mig_mu))
        
        p_new = SCARLikelihood.compute_like(ts,**params) + prior_new # posterior prob for proposed theta
        
        "Accept or reject proposal"
        a = np.math.exp(p_new - p_now) # exp since we're working with log likelihoods
        z = np.random.uniform(0,1)
        if z < a:
            #print('Proposal accepted')
            theta_now = theta_new
            accept += 1
            p_now = p_new
        
        "Sample params and probs"
        theta_samples[:,i] = theta_now
        p_samples[i] = p_now
        prior_arr[i] = prior_new
        
    "Plot trace plots for params"
    if plot:
        plot_traces(theta_samples,plot_info)
    
    #post_median = np.median(theta_samples,axis=1)
    estimates = np.percentile(theta_samples[:,burnin:],[2.5,50,97.5],axis=1)
    estimates = np.reshape(estimates, (3*dim,1), order='F') # reshape by dim
    return estimates,theta_samples,p_samples,prior_arr

def plot_traces(theta_samples,plot_info):

    dim,iterations = theta_samples.shape 
    sns.set(style="darkgrid")
    fig, axs = plt.subplots(dim, 1)
    for i in range(dim):
        if dim > 1:
            sns.lineplot(x=list(range(iterations)), y=theta_samples[i,:], ax=axs[i])
        else:
            sns.lineplot(x=list(range(iterations)), y=theta_samples[i,:], ax=axs)
        fig.tight_layout()
        true_val = f"{plot_info['true_val_rho']:.2f} {plot_info['true_val_M0']:.2f} {plot_info['true_val_M1']:.2f} {plot_info['true_val_M2']:.2f}"
        # sim = str(plot_info['sim'])
        fig_name = 'mcmc_trace_' + true_val + 'Aspergillus' + '.png'
        fig.set_size_inches(6, 12)
        fig.savefig(fig_name, dpi=200)

def plot_pair_densities(theta_samples):
    
    "Plot pairwise joint densities and scatter plot"
    import pandas as pd
    data = {'M_12':theta_samples[0,:],'M_21':theta_samples[1,:]}
    df = pd.DataFrame(data)
    g = sns.PairGrid(df, diag_sharey=False)
    g.map_upper(sns.scatterplot)
    g.map_lower(sns.kdeplot)
    g.map_diag(sns.kdeplot, lw=2)

if __name__ == '__main__':
    
    import ARGSimulator


    "Specify sim params"
    sims = 100
    samples = 30,30,40
    genome_length = 1e4
    rho = 0.1 #1e-1
    Ne = 1.0  # effective pop sizes
    mig_rate1 = 0.25
    mig_rate2 = 0.1
    mig_rate3 = 0.5
    M = [[0.0,mig_rate1,mig_rate2], [mig_rate1,0.0,mig_rate3], [mig_rate2,mig_rate3,0.0]]  # migration rate matrix
    sample_sizes = np.array(samples) # should be 1xpops vector
    
 
    
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
        
    
    true_fit_vals = np.zeros([sims,n_params])
    est_fit_vals = np.zeros([sims,3*n_params])
    #sim_Ne_vals = np.logspace(-1.0, 3.0, base = 10, num=sims) # e.g. from 10e-1 to 1000
    sim_mig_vals = np.linspace(0.1,2.0,num=sims)
    sim_rho_vals = np.linspace(0.1,10.0,num=sims)
    for s in range(sims):
            
        print('Sim: ', str(s))
        mig_rate1 = sim_mig_vals[s]
        mig_rate2 = 0.8*sim_mig_vals[s]
        mig_rate3 = 0.5*sim_mig_vals[s]
        params['M'] = [[0.0,mig_rate1,mig_rate2], [mig_rate1,0.0,mig_rate3], [mig_rate2,mig_rate3,0.0]]
        rho = sim_rho_vals[s]
        params['rho'] = rho
        true_fit_vals[s] = rho,mig_rate1,mig_rate2,mig_rate3
        
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
        plot_info = {'sim':str(s), 'true_val_rho':rho, 'true_val_M0':mig_rate1, 'true_val_M1':mig_rate2, 'true_val_M2':mig_rate3}
        opt_iters = []
        opt_iters = [200,500,1000,2000,5000]
        estimates = sample(params,est_params,ts,iterations=2000,opt_iters=opt_iters,burnin=1000,plot=True,**plot_info)
        est_fit_vals[s,:] = estimates[:,0]
    
        np.savetxt('test_SCAR_true_mig_unknownStates.csv', true_fit_vals, delimiter=',')
        np.savetxt('test_SCAR_mcmc_mig_unknownStates_ests.csv', est_fit_vals, delimiter=',')
                