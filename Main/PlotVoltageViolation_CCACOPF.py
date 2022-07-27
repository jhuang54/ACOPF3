import pandapower as pp
#import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

import matplotlib.pyplot as plt

from scipy.sparse.csgraph import minimum_spanning_tree

import pandas as pd

import numpy as np

path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_plt = os.path.join(path_cur, 'Plot')
path_input=os.path.join(path_cur, 'Input')

# load voltage violation of CCACOPF/
RECORD_V_VIOLATION_EDCOPF=np.load(path_plt+'/CCACOPF/noise/EDCOPFVViolation.npy')
RECORD_V_VIOLATION_ACOPF=np.load(path_plt+'/CCACOPF/noise/ACOPFVViolation.npy')
RECORD_V_VIOLATION_CCACOPF=np.load(path_plt+'/CCACOPF/noise/CCACOPFVViolation.npy')
n_snaps=len(RECORD_V_VIOLATION_ACOPF)

# # load cost
# np_load_old = np.load

# # modify the default parameters of np.load
# np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)

cost_EDCOPF=np.load(path_plt+'/CCACOPF/noise/EDCOPFCost.npy')
cost_ACOPF=np.load(path_plt+'/CCACOPF/noise/ACOPFCost.npy')
cost_CCACOPF=np.load(path_plt+'/CCACOPF/noise/CCACOPFCost.npy')

path_input_ccacopf=os.path.join(path_input,'CCACOPF.txt')
with open(path_input_ccacopf,'r') as f:
    lines = f.readlines()
    for line in lines:
        if line[0]!='#':
           input_t=line.split(',')
           n_sim=int(input_t[0])     
f.close()

xlabels_time=['00:00','06:00','12:00','18:00','23:55']
ylabels_count=['0','50','100']

iplt=1 
plt.figure(iplt)
fig,ax=plt.subplots(1)
plt.figure(iplt)
plt0=ax.plot(range(0,n_snaps),RECORD_V_VIOLATION_EDCOPF/n_sim,'.r',markersize=3)
plt1=ax.plot(range(0,n_snaps),RECORD_V_VIOLATION_ACOPF/n_sim,'.b',markersize=3)
plt2=ax.plot(range(0,n_snaps),RECORD_V_VIOLATION_CCACOPF/n_sim,'.m',markersize=3)
plt3=ax.plot(range(0,n_snaps),0.5*np.ones(n_snaps),'--',color='orange',linewidth=2)
plt4=ax.plot(range(0,n_snaps),0.05*np.ones(n_snaps),'--',color='g',linewidth=2)
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)

#ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))))
#ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))/n_sim))
ax.set_yticks([0,0.5,1])
ax.set_yticklabels(ylabels_count, rotation=0)

ax.set_xlabel('Time')
ax.set_xlim([0,n_snaps])
ax.set_ylim([0,1.02])
ax.set_ylabel('%')
ax.set_title('Voltage Violation rate')
ax.legend((plt0[0],plt1[0],plt2[0],plt3[0],plt4[0]),('DCOPF','ACOPF','CCACOPF','50%','5%'),fontsize = 6,loc='upper left',bbox_to_anchor=(0.05, 0.4, 0.5, 0.5))
ax.grid(True)
fig.savefig(path_plt+'/CCACOPF/noise'+'/VViolation(OPFs).png', dpi=400) 


iplt+=1 
plt.figure(iplt)
fig,ax=plt.subplots(1)
plt0=ax.plot(range(0,n_snaps),cost_EDCOPF,'r*',markersize=7)
plt1=ax.plot(range(0,n_snaps),cost_ACOPF,'o',color='b',markersize=4)
plt2=ax.plot(range(0,n_snaps),cost_CCACOPF,'m.',markersize=2)
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)

# #ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))))
# #ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))/n_sim))
# ax.set_yticks([0,0.5,1])
# ax.set_yticklabels(ylabels_count, rotation=0)

ax.set_xlabel('Time')
ax.set_xlim([0,n_snaps])
#ax.set_ylim([0,1.02])
ax.set_ylabel('$')
ax.set_title('Generation Cost')
ax.legend((plt0[0],plt1[0],plt2[0]),('DCOPF','ACOPF','CCACOPF'),fontsize = 4)
ax.grid(True)
fig.savefig(path_plt+'/CCACOPF/'+'/Generation Cost(OPFs).png', dpi=400) 

dcost_ACOPF=abs(cost_EDCOPF-cost_ACOPF)/abs(cost_EDCOPF)
dcost_CCACOPF=abs(cost_EDCOPF-cost_CCACOPF)/abs(cost_EDCOPF)

iplt+=1 
plt.figure(iplt)
fig,ax=plt.subplots(1)
plt0=ax.plot(range(0,n_snaps),np.zeros(n_snaps),'r.',markersize=7)
plt1=ax.plot(range(0,n_snaps),dcost_ACOPF*100,'b.',markersize=2)
plt2=ax.plot(range(0,n_snaps),dcost_CCACOPF*100,'m.',markersize=2)
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)

# #ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))))
# #ax.set_yticks(list(range(0,n_sim+1,int(n_sim/(len(ylabels_count)-1)))/n_sim))
# ax.set_yticks([0,0.5,1])
# ax.set_yticklabels(ylabels_count, rotation=0)

ax.set_xlabel('Time')
ax.set_xlim([0,n_snaps])
ax.set_ylim([-0.001,max(dcost_CCACOPF)*1.01*100])
#ax.set_ylim([0,1.02])
ax.set_ylabel('%')
ax.set_title('Deviation of Generation Cost from DC OPF')
ax.legend((plt0[0],plt1[0],plt2[0]),('DCOPF','ACOPF','CCACOPF'),fontsize = 4)
ax.grid(True)
fig.savefig(path_plt+'/CCACOPF/noise'+'/Generations deviation from DC OPF.png', dpi=400) 