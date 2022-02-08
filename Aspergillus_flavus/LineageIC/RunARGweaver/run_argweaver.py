# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:36:36 2020
functions 'arg_cmd', 'run_argweaver' are from tskit: ts_ARGweaver.py

@author: fguo7
"""


import logging
import shutil
import sys
import tempfile
import subprocess


    
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
            MSMC_samples, sample_step, burnin_iterations, ntimes, maxtime, verbose=False):
        """
        this produces a whole load of .smc files labelled <path_prefix>i.0.smc,
        <path_prefix>i.10.smc, etc, which we then convert into sts format

        Returns the iteration numbers ('0', '10', '20' etc), the name of the
        statistics file, the total CPU time, and the max memory usage.
        
        Modified a little from the original function of tskit evaluation.
        """
        

        burn_prefix = None
        exe = [ARGweaver_executable, '--sites', seq_file.name if hasattr(seq_file, "name") else seq_file,
               '--popsize', str(Ne),
               '--recombrate', str(recombination_rate),
               '--mutrate', str(mutation_rate),
               '--ntimes', str(ntimes),
               '--maxtime', str(maxtime),
               # '--delta', str(0.05),
               '--compress-seq', str(5),   # To speed up the algorithm, sequence alignment can be compressed during sampling into blocks of sizes, such as 10bp
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
            exe += ['--sites', seq_file]
        
        new_prefix = path_prefix + '_' + 'aflavus' #we append a '_sim' to mark iteration number
        
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






"paramterss"
rho_per_site = 2.17e-10  # recombination rate per site per generation
mutation_rate = 2.82e-09    # mutation rate per site per generation
Ne = 9674  # effective population size
maxtime = 19000     # maxtime of ARG, can set as TMRCA time
MSMC_samples = 2001     # iterations = int(sample_step * (MSMC_samples-1))
burnin_iterations = 1000    # burn_in iterations, can set as 0


"paths"
ARGweaver_executable = "arg-sample"     # the executable path of arg-sample

path = "./aflavus/"     # create a folder named aflavus to store sequence file and population table
smc_path = './infer/'   # create a folder named infer to store infered smc files

aflavus_fasta = path + "48IC_4pops_m07SNP_sites.sites"
smc_files_prefix = smc_path + 'smc-out'


"Infer an ARG with ARGweaver"

"""
# parameters of run_argweaver(seq_file, Ne, recombination_rate, 
# mutation_rate, path_prefix, iterations, sample_step, burnin_iterations, ntimes)
# iterations = int(sample_step * (MSMC_samples-1))
"""

run_argweaver(seq_file=aflavus_fasta, Ne=Ne, recombination_rate=rho_per_site, 
               mutation_rate=mutation_rate, path_prefix=smc_files_prefix, seed=70923, MSMC_samples=MSMC_samples, 
               sample_step=10, burnin_iterations=burnin_iterations, ntimes=20, maxtime=maxtime, verbose=False)



print("Running is OK")






    
    


