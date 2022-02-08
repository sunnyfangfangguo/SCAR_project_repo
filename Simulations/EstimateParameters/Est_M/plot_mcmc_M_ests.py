"""
Created on Thu Apr 16 18:19:22 2020

@author: david
"""

"Plot scatter of true vs. ests"
import numpy as np
import matplotlib.pyplot as plt

path = './'
true_vals = np.loadtxt(path+'test_SCAR_true_mig_unknownStates.csv',delimiter=",")
est_vals = np.loadtxt(path+'test_SCAR_mcmc_mig_unknownStates_ests.csv',delimiter=",")

log_transform = False
if log_transform:
    true_vals = np.log10(true_vals)
    est_vals = np.log10(est_vals)

intervals = [est_vals[:,1] - est_vals[:,0],est_vals[:,2] - est_vals[:,1]]
fig, ax = plt.subplots(figsize=(6,5))
ax.errorbar(true_vals, est_vals[:,1], yerr=intervals, fmt='o', color='black',
         ecolor='lightblue', elinewidth=3, capsize=0)

ax.set_xlabel('True M')
ax.set_ylabel('Estimated M')

xlim = ax.get_ylim()
ylim = ax.get_ylim()
ax.plot( [xlim[0],xlim[1]],[ylim[0],ylim[1]], 'k')

ax.set_xlim([-0.1,2.1])
ax.set_ylim([-0.1,6.5])

fig.tight_layout()
plt.show()
fig.savefig('test_SCAR_mcmc_true_vs_estM_unknownStates.png', dpi=200)