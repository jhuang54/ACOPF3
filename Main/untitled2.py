# -*- coding: utf-8 -*-
"""
Created on Mon May 23 11:59:52 2022

@author: jhuang
"""

Total Pgen Adjustment (Mw): 0.04005862512250084
Total Qgen Adjustment (Mvar): 0.6526696553583429
Total Shunt Adjustment (Mvar): 0.0


iplt+=1
plt.figure(iplt) 
dqsg_a=(qsg_t-qsg)*sbase
id_inac=np.where(abs(qsg_t)<1e-3)
dqsg_r=(qsg_t-qsg)/qsg_max*100
dqsg_r[id_inac]=0

dqvg_a=(qg_t-qg)*sbase
id_inac=np.where(abs(qg_t)<1e-4)
dqvg_r=(qg_t-qg)/qll_g_max[busid_g_LL]*100
dqvg_r[id_inac]=0

dqg_a=np.concatenate((dqvg_a,dqsg_a))
dqg_r=np.concatenate((dqvg_r,dqsg_r))
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xlabel('Shunt index')
ax1.set_ylabel('Mvar',color=color)
ax1.plot(range(1,19),dqg_a[-18:],'.',markersize=3,color=color)
ax1.tick_params(axis='y', labelcolor='tab:red')
plt.ylim([0,1])
ax2=ax1.twinx()
color='tab:green'
ax2.set_ylabel('Percent %', color=color)  # we already handled the x-label with ax1
ax2.plot(range(1,19),dqg_r[-18:],'.', color=color)
ax2.tick_params(axis='y', labelcolor=color)
plt.xlim([1,18])
plt.ylim([0,35])
plt.title('Reactive Power Adjustment')
plt.grid(True)
plt.savefig(path_plt+'/Shunt_Adjust.png', dpi=400) 