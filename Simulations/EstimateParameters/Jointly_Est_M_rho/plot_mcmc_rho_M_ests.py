"""
Created on Thu Apr 16 18:19:22 2020

@author: david
"""

"Plot scatter of true vs. ests"
import numpy as np
import matplotlib.pyplot as plt

path = './'
true_vals = np.loadtxt(path+'test_SCAR_true_rho_mig_unknownStates1.csv',delimiter=",")
est_vals = np.loadtxt(path+'test_SCAR_mcmc_rho_mig_unknownStates_ests1.csv',delimiter=",")

log_transform = False
if log_transform:
    true_vals = np.log10(true_vals)
    est_vals = np.log10(est_vals)

intervals_rho = [est_vals[:,1] - est_vals[:,0],est_vals[:,2] - est_vals[:,1]]
intervals_mig = [est_vals[:,4] - est_vals[:,3],est_vals[:,5] - est_vals[:,4]]
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

"subplot-r"
ax1.errorbar(true_vals[:,0], est_vals[:,1], yerr=intervals_rho, fmt='o', color='black',
         ecolor='lightblue', elinewidth=3, capsize=0)

ax1.set_xlabel('True r')
ax1.set_ylabel('Estimated r')

xlim1 = ax1.get_ylim()
ylim1 = ax1.get_ylim()
ax1.plot( [xlim1[0],xlim1[1]],[ylim1[0],ylim1[1]], 'k')

ax1.set_xlim([0,1.05e-03])
ax1.set_ylim([0,1.3e-03])
ax1.text(0.07, 0.96, 'A', transform=ax1.transAxes,
      fontsize=18, fontweight='medium', va='top', ha='right')

"subplot-migration"
ax2.errorbar(true_vals[:,1], est_vals[:,4], yerr=intervals_mig, fmt='o', color='black',
         ecolor='lightblue', elinewidth=3, capsize=0)

ax2.set_xlabel('True M')
ax2.set_ylabel('Estimated M')

xlim2 = ax2.get_ylim()
ylim2 = ax2.get_ylim()
ax2.plot( [xlim2[0],xlim2[1]],[ylim2[0],ylim2[1]], 'k')

ax2.set_xlim([-0.1,2.1])
ax2.set_ylim([-0.1,5])
ax2.text(0.07, 0.96, 'B', transform=ax2.transAxes,
      fontsize=18, fontweight='medium', va='top', ha='right')




fig.tight_layout()
plt.show()
fig.savefig('test_SCAR_mcmc_true_vs_estRhoM_unknownStates.png', dpi=200)
