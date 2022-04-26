# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 21:36:19 2022

@author: jhuang
"""
import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
# apparent power flow limit model


# branch table [fbus,tbus]
# brh_fbus=np.concatenate((net.line.from_bus.values,net.trafo.hv_bus.values))
# brh_tbus=np.concatenate((net.line.to_bus.values,net.trafo.lv_bus.values))
brh_fbus=np.array([3,4,1,3,3,2,3,3,5,6])
brh_tbus=np.array([1,1,2,2,4,4,5,6,0,0])
nbrh=len(brh_fbus)

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
B=np.concatenate((Bt,np.identity(nclp)), axis=1) 