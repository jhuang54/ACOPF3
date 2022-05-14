# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 21:36:19 2022

@author: jhuang
"""
import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree

import matplotlib.pyplot as plt
# apparent power flow limit model


# branch table [fbus,tbus]
# brh_fbus=np.concatenate((net.line.from_bus.values,net.trafo.hv_bus.values))
# brh_tbus=np.concatenate((net.line.to_bus.values,net.trafo.lv_bus.values))
brh_fbus=np.array([3,4,1,3,3,2,3,3,5,6])
brh_tbus=np.array([1,1,2,2,4,4,5,6,0,0])
nbrh=len(brh_fbus)

brh_y=np.zeros(nbrh,dtype=complex)
# branch impedance
for i in range(nbrh):
    brh_y[i]=1/(0.1+0.2j)
brh_z=1/brh_y

busid_N=[0,1,2,3,4,5,6]
busid_LL=[1,2,3,4,5,6]
busid_slack=0

brh_t=np.array([],dtype='int64').reshape(0)# tree branch

nbus=len(busid_N)
nLL=len(busid_LL)
# branch_tb[:,0]=fbus
# branch_tb[:,1]=tbus


A0=np.zeros((nbus,nbus))
bus_tp=np.array([],dtype='int64').reshape(0)
busid_N=list(range(nbus))
for ibus in range(nbus):
    busid_tp=busid_N[ibus]
    idbrh_f=np.where(brh_fbus==busid_tp)[0]#branch id of branches whose fbus is busid_tp
    bus_tp=brh_tbus[idbrh_f]
    
    idbrh_t=np.where(brh_tbus==busid_tp)[0]
    bus_tp=np.concatenate((bus_tp,brh_fbus[idbrh_t]))
    
    A0[busid_tp,bus_tp]=1
    A0[bus_tp,busid_tp]=1    

tree=minimum_spanning_tree(A0)
out_ind = np.transpose(np.nonzero(tree))

# tree branch id
# tree branch id
id_brh_tb=np.zeros((nbus,nbus),dtype='int64')
for i in range(nbrh):
    fbus_tp=brh_fbus[i]
    tbus_tp=brh_tbus[i]
    id_brh_tb[fbus_tp,tbus_tp]=i
    id_brh_tb[tbus_tp,fbus_tp]=i

#brh_t=np.zeros(1)
brh_t=np.array([],dtype='int64').reshape(0)# tree branch
for i in range(nLL):
    brh_id_tp=id_brh_tb[out_ind[i,0],out_ind[i,1]]
    brh_t=np.concatenate((brh_t,[brh_id_tp]))
    
    
# Afti
A=np.zeros((nLL,nbrh))
for ibus in range(nLL):
    busid_tp=busid_LL[ibus]
    idbrh_f=np.where(brh_fbus==busid_tp)[0]#branch id of branches whose fbus is busid_tp
    idbrh_t=np.where(brh_tbus==busid_tp)[0]
    A[ibus,idbrh_f]=1
    A[ibus,idbrh_t]=-1
    
# select a tree in a meshed network 
bus_pool=[busid_slack]
Nextbus0=[busid_slack]
Nextbus=np.array([],dtype='int64').reshape(0)
brh_l=np.array([],dtype='int64').reshape(0)# lian branch
#brh_t=np.array([],dtype='int64').reshape(0)# tree branch


# while len(brh_t)<nLL:
#     for i_next in range(len(Nextbus0)):
#         busid_tp=Nextbus0[i_next]# last 
        
#         # its branches
#         brh_t_tp=np.where(brh_fbus==busid_tp)[0]
#         brh_f_tp=np.where(brh_tbus==busid_tp)[0]
        
#         # its adjacent buses
#         Tbus_tp=brh_tbus[brh_t_tp]
#         Fbus_tp=brh_fbus[brh_f_tp]
        
#         # next bus   
#         Nextbus_tp=np.concatenate((Fbus_tp,Tbus_tp))
        
#         # aggregate branches
#         brh_tp=np.concatenate((brh_f_tp,brh_t_tp))
        
#         # # collect lian branch and tree branch
#         # mask_lb=np.isin(Nextbus_tp, bus_pool)
#         # brh_l=np.concatenate((brh_l,brh_tp[mask_lb]))
        
#         mask_tb=np.isin(Nextbus_tp,bus_pool,invert=True)
#         brh_t=np.concatenate((brh_t,brh_tp[mask_tb]))
            
#         # update next bus
#         Nextbus=np.concatenate((Nextbus,Nextbus_tp[mask_tb]))
#         bus_pool=np.concatenate((bus_pool,Nextbus_tp[mask_tb]))
#     Nextbus0=Nextbus  
#     # bus_pool=np.concatenate((bus_pool,Nextbus))
#     Nextbus=np.array([],dtype='int64').reshape(0)
    
# loop matrix
At=A[:,brh_t]#tree branch
brh_l=np.setdiff1d(range(nbrh),brh_t)
Al=A[:,brh_l]#lian branch
Btt=np.matmul(-np.linalg.inv(At), Al)
Bt=np.transpose(Btt)

#brh_t=np.concatenate()
nclp=nbrh-nLL
B0=np.concatenate((Bt,np.identity(nclp)), axis=1)  
id_brh=np.concatenate((brh_t,brh_l))

# re-ōrder column id back to [0,1,2,...,nbrh]
B1=np.zeros((nclp,nbrh))
for i in range(nbrh):
    B1[:,id_brh[i]]=B0[:,i]
    
# B
# conjuage of branch impedance matrix Zf
brh_r=np.real(brh_z)
brh_x=np.imag(brh_z)
Bp=np.matmul(B1,np.diag(brh_r))
Apf=np.concatenate((A,Bp), axis=0)
Bpf=np.linalg.inv(Apf)
Rf=Bpf[:,0:nLL]

Bq=np.matmul(B1,np.diag(brh_x))
Aqf=np.concatenate((A,Bq), axis=0)
Bqf=np.linalg.inv(Aqf)
Xf=Bqf[:,0:nLL]


# power flow model
Ybrh=np.zeros((nbus,nbus),dtype=complex)
for i in range(nbrh):
    Ybrh[brh_fbus[i],brh_tbus[i]]=brh_y[i]
    Ybrh[brh_tbus[i],brh_fbus[i]]=brh_y[i]
    
Y=-Ybrh
for i in range(nbus):
    Y[i,i]=-sum(Y[i,:])
busid_LL=[1,2,3,4,5,6]
YLL=Y[np.ix_(busid_LL,busid_LL)]
Z_pu=np.linalg.inv(YLL)
busid_slack=0
YLS=Y[np.ix_(busid_LL,[busid_slack])]

YY=np.matmul(Z_pu,YLS)
vm0=0.93
va0=np.exp(1j*0*np.pi/180)
v0=vm0*va0

vs=YY*v0
vs=np.squeeze(np.asarray(vs))

w=-vs

vll=np.ones(nLL,dtype=complex)
sll_net_t=1e-1*(-1-1j)*np.ones(nLL,dtype=complex)
vll0=vll

itr_pf=0
itr_pf_max=20
err_vll=1
while (err_vll>1e-5 and itr_pf<itr_pf_max):
    #vll=Z_pu*np.conj(sll_net_t/vll)+vs
    ILL=np.conj(sll_net_t/vll)
    dv=np.matmul(Z_pu,ILL.transpose())
    vll=dv+w
    err_vll=max(abs(vll-vll0))
    vll0=vll
    itr_pf=itr_pf+1
           
if err_vll>1e-3:
    raise Exception('power flow diverge!\n')
    
# nonlinear complex power flow
v_t=np.ones(nbus,dtype=complex)
v_t[busid_slack]=v0
v_t[busid_LL]=vll
vf_t=v_t[brh_fbus]
vt_t=v_t[brh_tbus]
If_t=(vf_t-vt_t)*brh_y
Sf_t=vf_t*np.conjugate(If_t)
Pf_t=np.real(Sf_t)
Qf_t=np.imag(Sf_t)
Sf_t=Pf_t+1j*Qf_t
Sfm_t=abs(Sf_t)

# linear model
# Pf==Cp-Dq
# Qf=Dp+Cq
Pi_t=np.real(sll_net_t)
Qi_t=np.imag(sll_net_t)
# Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
# Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
Pf_e=np.matmul(Rf,Pi_t)
Qf_e=np.matmul(Xf,Qi_t)


dpf=Pf_t-Pf_e
dqf=Qf_t-Qf_e

Sfe_t=Pf_e+1j*Qf_e
Sfe_t=Pf_e+1j*Qf_e

# check
# kvl: A1,B1
np.matmul(Bp,Pf_e)+np.matmul(Bq,Qf_e)
# kcl: A
np.matmul(A,Pf_t)-Pi_t
np.matmul(A,Qf_t)-Qi_t

class result:
    def __init__(self, iter_max):
        self.vmax=np.full((iter_max,1),1.05)
        self.vmin=np.full((iter_max,1),0.95)
        self.lambda_u=np.zeros(iter_max)
        self.lambda_d=np.zeros(iter_max)
        self.dSfm_t=np.zeros((iter_max,nbrh))
        self.dSfmp_t=np.zeros((iter_max,nbrh))
        self.dPfm_t=np.zeros((iter_max,nbrh))
        self.dPfmp_t=np.zeros((iter_max,nbrh))
        self.dQfm_t=np.zeros((iter_max,nbrh))
        self.dQfmp_t=np.zeros((iter_max,nbrh))
iter_max=1
result1=result(iter_max) 
iterat=0
result1.dSfm_t[iterat,:]=abs(Sfm_t-abs(Sfe_t))
id_inac=np.where(Sfm_t<1e-4)
result1.dSfmp_t[iterat,:]=result1.dSfm_t[iterat,:]/Sfm_t
result1.dSfmp_t[iterat,id_inac]=0
result1.dPfm_t[iterat,:]=abs(Pf_t-Pf_e)
result1.dPfmp_t[iterat,:]=abs(Pf_t-Pf_e)/abs(Pf_t)
id_inac=np.where(abs(Pf_t)<1e-4)
result1.dPfm_t[iterat,id_inac]=0
result1.dPfmp_t[iterat,id_inac]=0

result1.dQfm_t[iterat,:]=abs(Qf_t-Pf_e)
result1.dQfmp_t[iterat,:]=abs(Qf_t-Pf_e)/abs(Qf_t)
id_inac=np.where(abs(Qf_t)<1e-4)
result1.dQfm_t[iterat,id_inac]=0
result1.dQfmp_t[iterat,id_inac]=0

plt.figure(1)
plt.plot(result1.dPfm_t[iterat,:])
plt.title('Absolute error of Pfm')
plt.grid(True)
plt.figure(2)
plt.plot(result1.dPfmp_t[iterat,:])
plt.title('Relative error of Pfm')
plt.grid(True)

plt.figure(3)
plt.plot(result1.dQfm_t[iterat,:])
plt.title('Absolute error of Qfm')
plt.grid(True)
plt.figure(4)
plt.plot(result1.dQfmp_t[iterat,:])
plt.title('Relative error of Qfm')
plt.grid(True)

