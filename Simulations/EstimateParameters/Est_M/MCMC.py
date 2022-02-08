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

def theta2params(theta,params,param_map):
    
    "Update params dict with values in theta"
    for k,v in param_map.items():
        if k == 'M':
            pops = len(params['M'])
            M = np.ones([pops,pops]) * theta[v]
            params['M'] = M - np.diag(np.diag(M))
        else:
            params[k] = theta[v]
            
    return params

def params2theta(params,theta,param_map):
    
    "Flatten params dict to theta vec"
    for k,v in param_map.items():
        if k == 'M':
            mig_rate = params['M'][0][1]
            theta[v] = mig_rate
        else:
            theta[v] = params[k]
            
    return theta
    
    
def sample(params,est_params,ts,iterations=1000,opt_iters=[],burnin=0,plot=False,**plot_info):
    
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
        param_map['M'] = param_map_index
        param_map_index += 1
    
    prior = 1.0 # assume uniform prior
    p_now = SCARLikelihood.compute_like(ts,**params) * prior # current posterior prob
    p_samples = np.zeros(iterations)
    
    "Set up parameters and samples of theta"
    theta_now = np.array([0.0]*param_map_index)
    #theta_now = np.array([params[k] for k in param_map]) # current param value
    theta_now = params2theta(params,theta_now,param_map)
    dim = theta_now.size
    theta_samples = np.zeros((dim,iterations))
    theta_samples[:,0] = theta_now

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
        p_new = SCARLikelihood.compute_like(ts,**params) * prior # posterior prob for proposed theta
        
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
        
    "Plot trace plots for params"
    if plot:
        plot_traces(theta_samples,plot_info)
    
    #post_median = np.median(theta_samples,axis=1)
    estimates = np.percentile(theta_samples[:,burnin:],[2.5,50,97.5],axis=1)
    return estimates

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
        true_val = f"{plot_info['true_val']:.2f}"
        sim = str(plot_info['sim'])
        fig_name = 'mcmc_trace_' + true_val + '_sim' + sim + '.png'
        fig.set_size_inches(6, 6)
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

    samples = 10
    genome_length = 1e4
    rho = 0 #1e-1
    Ne = 200.0  # effective pop sizes
    M = [[0.0]]
    sample_sizes = [samples]*1 # should be 1xpops vector
    
    "Store params in dict"
    params = {'Ne': Ne, 'rho': rho, 'M': M, 'genome_length': genome_length, 'sample_sizes':  sample_sizes}
    
    "Pass estimate params as dict"
    est_params = {'Ne': True, 'rho': False, 'M': False}

    "Run tree simulation"
    ts = ARGSimulator.sim_ARG(params,min_breakpoints=0,plot=False)
    breaks = len(ts.breakpoints(as_array=True))
    #ts_sim_file = 'coal_Ne' + f"{Ne:.2f}" + '_sim' + str(0)
    #ts.dump(ts_sim_file, zlib_compression=False)
        
    "Run MCMC inference"
    plot_info = {'sim':str(0) , 'true_val':Ne}
    opt_iters = [200,500,1000]
    estimates = sample(params,est_params,ts,iterations=2000,opt_iters=opt_iters,burnin=500,plot=True,**plot_info)    

    print(estimates)
