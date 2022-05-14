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
brh_fbus=np.array([0,1,0])
brh_tbus=np.array([1,2,2])
nbrh=len(brh_fbus)

brh_y=np.zeros(nbrh,dtype=complex)
# branch impedance
for i in range(nbrh):
    brh_y[i]=1/(1+1j)
brh_z=1/brh_y

busid_N=[0,1,2]
busid_LL=[1,2]
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
B=np.matmul(B1,np.conjugate(np.diag(brh_z)))
Af=np.concatenate((A,B), axis=0)
Bf=np.linalg.inv(Af)
Bfti=Bf[:,0:nLL]
C=np.real(Bfti)
Ct=np.transpose(C)
D=np.imag(Bfti)
Dt=np.transpose(D)

# power flow model
Ybrh=np.zeros((nbus,nbus),dtype=complex)
for i in range(nbrh):
    Ybrh[brh_fbus[i],brh_tbus[i]]=brh_y[i]
    Ybrh[brh_tbus[i],brh_fbus[i]]=brh_y[i]
    
Y=-Ybrh
for i in range(nbus):
    Y[i,i]=-sum(Y[i,:])
busid_LL=[1,2]
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
Wm=np.diag(w)
Wm_inv=np.diag(1/w)
Kyr=np.matmul(np.matmul(np.matmul(abs(Wm),Wm_inv),Z_pu),np.conj(Wm_inv))
#R=np.real(1/vm0*np.conj(v0)*Zt)
R=np.real(Kyr)
#Rt=csc_matrix.transpose(R)
Rt=np.transpose(R)

#X=np.real(1/vm0*np.conj(v0)*(-1j)*Zt)
X=np.real(-1j*Kyr)
#Xt=csc_matrix.transpose(X)
Xt=np.transpose(X)


vll=np.ones(nLL,dtype=complex)
#sll_net_t=1e-2*(-1-1j)*np.ones(nLL,dtype=complex)
sll_net_t0=np.array([1e-2*(-1-0.5j),1e-2*(2+1j)])
sll_net_t=np.array([1e-2*(-1-0.5j),1e-2*(2+1j)])
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

Sf_t0=vf_t*np.conjugate(If_t)
Pf_t0=np.real(Sf_t)
Qf_t0=np.imag(Sf_t)

# linear model
# Pf==Cp-Dq
# Qf=Dp+Cq
Pi_t=np.real(sll_net_t)
Qi_t=np.imag(sll_net_t)
Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)

dpf=Pf_t-Pf_e
dqf=Qf_t-Qf_e
Sfe_t=Pf_e+1j*Qf_e

# initialize controllable loads
epsi_pg=0.05*np.ones((nLL))
epsi_qg=0.05*np.ones((nLL))
epsi_pg[0]=0
epsi_qg[0]=0
epsi_l=0.001

epsi_u=0.005#0.01 for np.min(Smax_brh)>=6, 0.001 for np.min(Smax_brh)>=5.58

v_l=0.90# bounds for opf
v_u=1.1
vplt_min=v_l*0.9# bounds for plot
vplt_max=v_u*1.1


# initial values of load 
# initial values of generator
pg_t=np.real(sll_net_t)
qg_t=np.imag(sll_net_t)

pll_g_t=np.zeros(nLL)
busid_g_LL=[0,1]
pll_g_t[busid_g_LL]=pg_t# load is negative in linearized power flow model
qll_g_t=np.zeros(nLL)
qll_g_t[busid_g_LL]=qg_t

# desired values of gen
pg=np.real(sll_net_t)
qg=np.imag(sll_net_t)

pll_g=np.zeros(nLL)
pll_g[busid_g_LL]=pg# load is negative in linearized power flow model
qll_g=np.zeros(nLL)
qll_g[busid_g_LL]=qg

# capacity of Generator
pll_g_max=np.zeros(nLL)
pll_g_max[busid_g_LL]=abs(np.real(sll_net_t))*2
pll_g_min=np.zeros(nLL)
pll_g_min[busid_g_LL]=-abs(np.real(sll_net_t))*2

qll_g_max=np.zeros(nLL)
qll_g_max[busid_g_LL]=abs(np.imag(sll_net_t))*2
qll_g_min=np.zeros(nLL)
qll_g_min[busid_g_LL]=-abs(np.imag(sll_net_t))*2


# remove slack bus
vmll=abs(vll)
vall=np.angle(vll)

vmll0=abs(vll)

iplt=0
# plt.figure(iplt)
# plot_vm=plt.plot(vmll,'.')
# plot_ub=plt.plot([v_u]*nLL,'--',linewidth=1)
# plot_lb=plt.plot([v_l]*nLL,'--',linewidth=1)
# plt.ylim(vplt_min, vplt_max)
# plt.title('v (p.u.)')
# plt.savefig(path_plt+'/VProfile0.png', dpi=400) 

print('max v:',vmll.max())
print('min v:',vmll.min())



lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)

u_u=np.zeros(nbrh)
u_l=np.zeros(nbrh)

vll0=vll

iter_max=2000
alpha_ld=0.1# f=alpha*(pld-\hat{pld})^2
alpha_g=0.1# f=alpha*(pg-\hat{pg})^2
alpha_sg=0.1
beta_ld=2.5# capacity of DER=rated value*beta

beta_pgmax=1# capacity of Generator=rated value*beta
beta_pgmin=0# capacity of Generator=rated value*beta
beta_qgmax=1.5
beta_qgmin=-1.5

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
result1=result(iter_max)  

# complex voltage of all the buses including slack bus
v_t=np.ones(nbus,dtype=complex)

Qfr=np.zeros(nbrh)
Pfr=np.zeros(nbrh)

result1.u_u=np.zeros(iter_max)
result1.u_l=np.zeros(iter_max)
iplt=0
sbase=1
Smax_brh=np.array([0.01,0.008,0.05])
for iterat in range(iter_max):
    # derivative of voltage constraints with respect to (p, q)
    dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
    dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
    dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
    dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array

    
    # derivative of apparent power flow with respect to (p_g,p_pv,p_battery,p_load)
    Qfrdu=Qfr*(u_u-u_l)
    Pfrdu=Pfr*(u_u-u_l)
    
    dsf_dp=Dt.dot(Qfrdu)+Ct.dot(Pfrdu)
    dsf_dp=np.squeeze(np.array(dsf_dp))
    dsf_dq=Ct.dot(Qfrdu)-Dt.dot(Pfrdu)
    dsf_dq=np.squeeze(np.array(dsf_dq))
    
    
    # each bus has at most 1 load and generator, but may have multiple static generators
    
    # minimize generation from coal generator
    pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)+dvcnstr_dp+dsf_dp)
    qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)+dvcnstr_dq+dsf_dq)
    
    # project
    # pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
    # pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
    # qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
    # qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
    
    pll_g_t=np.maximum(pll_g_t,pll_g_min)
    pll_g_t=np.minimum(pll_g_t,pll_g_max)
    # qll_g_t=np.maximum(qll_g_t,qll_g_min)
    # qll_g_t=np.minimum(qll_g_t,qll_g_max)
    
    
    # update net injections 
    pll_net_t=pll_g_t
    qll_net_t=qll_g_t
    sll_net_t=pll_net_t+1j*qll_net_t
    
    #pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
    

    # fixed point iteration methods as the power flow solver
    vll0=vll
    
    itr_pf=0
    itr_pf_max=20
    err_vll=1
    while (err_vll>1e-5 and itr_pf<itr_pf_max):
        #vll=Z_pu*np.conj(sll_net_t/vll)+vs
        ILL=np.conj(sll_net_t/vll)
        dv=np.matmul(Z_pu,ILL.transpose())
        vll=np.squeeze(np.asarray(dv))+w
        err_vll=max(abs(vll-vll0))
        vll0=vll
        itr_pf=itr_pf+1
               
    # if err_vll>1e-3:
    #     raise Exception('power flow diverge!\n')
    
         
    # # dispatch load 
    # pld_t=pll_ld_t[busid_ld_LL]
    # qld_t=qll_ld_t[busid_ld_LL]
    
    # net_t.load.p_mw[busid_ld_lds]=-pld_t*sbase# load is positive in net.load
    # net_t.load.q_mvar[busid_ld_lds]=-qld_t*sbase
    
    
    # update dual variables (voltage)
    vmll=abs(vll)
    lambda_u=lambda_u+epsi_l*(vmll-v_u)
    lambda_d=lambda_d+epsi_l*(v_l-vmll)
    
    # project dual variables
    lambda_u=np.maximum(lambda_u,0)
    lambda_d=np.maximum(lambda_d,0)
    
    # track voltage and dual
    result1.lambda_u[iterat]=lambda_u.max()
    result1.lambda_d[iterat]=lambda_d.max()
    # result1.u_u[iterat]=u_u.max()
    # result1.u_l[iterat]=u_l.max()
    result1.vmax[iterat]=vmll.max()
    result1.vmin[iterat]=vmll.min()
    
    #update dual variable (apparent power flow)
    # nonlinear complex power flow
    v_t[busid_slack]=v0
    v_t[busid_LL]=vll
    vf_t=v_t[brh_fbus]
    vt_t=v_t[brh_tbus]
    If_t=(vf_t-vt_t)*brh_y
    Sf_t=vf_t*np.conjugate(If_t)
    
    # # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    # # Pi_t=np.real(sll_net_t)
    # # Qi_t=np.imag(sll_net_t)
    # # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    
    # # dpf=np.real(Sf_t)-Pf_e
    # # dqf=np.imag(Sf_t)-Qf_e
    
    
    Sfm_t=abs(Sf_t)
    Pfr=np.real(Sf_t)/Sfm_t
    Qfr=np.imag(Sf_t)/Sfm_t
    
    # inactive small branch flow
    id_sf0=np.where(Sfm_t<1e-7)
    Pfr[id_sf0]=0
    Qfr[id_sf0]=0
    
    # dual variable
    u_u=u_u+epsi_u*(Sfm_t-Smax_brh)
    #u_l=u_l+epsi_u*(Smin_brh-Sfm_t)
    
    # project dual variables
    u_u=np.maximum(u_u,0)
    #u_l=np.maximum(u_l,0)   
    
    result1.u_u[iterat]=u_u.max()
    
    #result1.u_l[iterat]=u_l.max()
    
    # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    # Pi_t=np.real(sll_net_t)
    # Qi_t=np.imag(sll_net_t)
    
    # # Sn=np.matmul(A,Sf_t)
    # # Pi_t=np.real(Sn)
    # # Qi_t=np.imag(Sn)
    
    # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    # Sfe_t=Pf_e+1j*Qf_e
    
    
    # # result1.dSfm_t[iterat,:]=abs(Sfm_t-abs(Sfe_t))
    # id_inac=np.where(Sfm_t<1e-4)
    # result1.dSfmp_t[iterat,:]=result1.dSfm_t[iterat,:]/Sfm_t
    # result1.dSfmp_t[iterat,id_inac]=0
    # result1.dPfm_t[iterat,:]=abs(np.real(Sf_t)-Pf_e)
    # result1.dPfmp_t[iterat,:]=abs(np.real(Sf_t)-Pf_e)/abs(np.real(Sf_t))*100
    # id_inac=np.where(abs(np.real(Sf_t))<1e-4)
    # result1.dPfm_t[iterat,id_inac]=0
    # result1.dPfmp_t[iterat,id_inac]=0
    
    # result1.dQfm_t[iterat,:]=abs(np.imag(Sf_t)-Qf_e)
    # result1.dQfmp_t[iterat,:]=abs(np.imag(Sf_t)-Qf_e)/abs(np.imag(Sf_t))*100
    # id_inac=np.where(abs(np.imag(Sf_t))<1e-4)
    # result1.dQfm_t[iterat,id_inac]=0
    # result1.dQfmp_t[iterat,id_inac]=0
    
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dPfm_t[iterat,:])
    # plt.title('Absolute error of Pfm')
    # plt.ylim([0,0.18])
    # plt.grid(True)
    # plt.savefig(path_plt+'/Aerror_Pfm.png', dpi=400) 
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dPfmp_t[iterat,:])
    # plt.title('Relative error of Pfm')
    # plt.ylim([0,30])
    # plt.grid(True)
    # plt.savefig(path_plt+'/Rerror_Pfm.png', dpi=400)
    
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dQfm_t[iterat,:])
    # plt.title('Absolute error of Qfm')
    # plt.grid(True)
    # plt.ylim([0,0.5])
    # plt.savefig(path_plt+'/Aerror_Qfm.png', dpi=400)
    # #plt.ylim([0,0.18])
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dQfmp_t[iterat,:])
    # plt.title('Relative error of Qfm')
    # #plt.ylim([0,0.3])
    # plt.grid(True) 
    # plt.savefig(path_plt+'/Rerror_Qfm.png', dpi=400)
    
    # Sn=np.matmul(A,Sf_t)
    # Pi_t=np.real(Sn)
    # Qi_t=np.imag(Sn)
   
    # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    # Sfe_t=Pf_e+1j*Qf_e
    
    
    # # result1.dSfm_t[iterat,:]=abs(Sfm_t-abs(Sfe_t))
    # id_inac=np.where(Sfm_t<1e-4)
    # result1.dPfm_t[iterat,:]=abs(np.real(Sf_t)-Pf_e)
    # result1.dPfmp_t[iterat,:]=abs(np.real(Sf_t)-Pf_e)/abs(np.real(Sf_t))*100
    # id_inac=np.where(abs(np.real(Sf_t))<1e-4)
    # result1.dPfm_t[iterat,id_inac]=0
    # result1.dPfmp_t[iterat,id_inac]=0
   
    # result1.dQfm_t[iterat,:]=abs(np.imag(Sf_t)-Qf_e)
    # result1.dQfmp_t[iterat,:]=abs(np.imag(Sf_t)-Qf_e)/abs(np.imag(Sf_t))*100
    # id_inac=np.where(abs(np.imag(Sf_t))<1e-4)
    # result1.dQfm_t[iterat,id_inac]=0
    # result1.dQfmp_t[iterat,id_inac]=0
    
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dPfm_t[iterat,:])
    # plt.title('Absolute error of Pfm')
    # plt.grid(True)
    # plt.ylim([0,0.18])
    # plt.savefig(path_plt+'/Aerror_PfmLs.png', dpi=400) 
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dPfmp_t[iterat,:])
    # plt.title('Relative error of Pfm')
    # plt.ylim([0,30])
    # plt.grid(True)  
    # plt.savefig(path_plt+'/Rerror_PfmLs.png', dpi=400) 
    
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dQfm_t[iterat,:])
    # plt.title('Absolute error of Qfm')
    # plt.grid(True)
    # plt.ylim([0,0.5])
    # plt.savefig(path_plt+'/Aerror_QfmLs.png', dpi=400) 
    # iplt+=1
    # plt.figure(iplt)
    # plt.plot(result1.dQfmp_t[iterat,:])
    # plt.title('Relative error of Qfm')
    # #plt.ylim([0,0.3])
    # plt.grid(True) 
    # plt.savefig(path_plt+'/Rerror_QfmLs.png', dpi=400)  
    
    # # check equations
    # # (A,Sf)=Si
    # dkcl_t=np.matmul(A,Sf_t)-sll_net_t
    # iplt+=1
    # plt.figure(iplt)
    # plt_r=plt.plot(range(0,nLL),sbase*abs(np.real(dkcl_t)))
    # plt.title('Error of real part of kcl_t')
    # plt.savefig(path_plt+'/MisRKcl_t.png', dpi=400)  
    # iplt+=1
    # plt.figure(iplt)
    # plt_x=plt.plot(range(0,nLL),sbase*abs(np.imag(dkcl_t)))
    # plt.title('Error of imaginary part of kcl_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisXKcl_t.png', dpi=400) 
    
    # # (B,Sf)=0
    # dkvl_t=np.matmul(B,Sf_t)
    # iplt+=1
    # plt.figure(iplt)
    # plt_r=plt.plot(range(0,nclp),sbase*abs(np.real(dkvl_t)))
    # plt.title('Error of real part of kvl_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisRKvl_t.png', dpi=400) 
    # iplt+=1
    # plt.figure(iplt)
    # plt_x=plt.plot(range(0,nclp),sbase*abs(np.imag(dkvl_t)))
    # plt.title('Error of imaginary part of kvl_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisXKvl_t.png', dpi=400) 
    
    # # (A,Sf)=Si
    # dkcle_t=np.matmul(A,Sfe_t)-sll_net_t
    # iplt+=1
    # plt.figure(iplt)
    # plt_r=plt.plot(range(0,nLL),abs(np.real(dkcle_t)))
    # plt.title('Error of real part of kcle_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisRKcle_t.png', dpi=400) 
    # iplt+=1
    # plt.figure(iplt)
    # plt_x=plt.plot(range(0,nLL),abs(np.imag(dkcle_t)))
    # plt.title('Error of imaginary part of kcl_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisXKcle_t.png', dpi=400)
    
    # # (B,Sf)=0
    # dkvle_t=np.matmul(B,Sfe_t)
    # iplt+=1
    # plt.figure(iplt)
    # plt_r=plt.plot(range(0,nclp),abs(np.real(dkvle_t)))
    # plt.title('Error of real part of kvle_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisRKvle_t.png', dpi=400)    
    # iplt+=1
    # plt.figure(iplt)
    # plt_x=plt.plot(range(0,nclp),abs(np.imag(dkvle_t)))
    # plt.title('Error of imaginary part of kvle_t')
    # plt.grid(True)
    # plt.savefig(path_plt+'/MisXKvle_t.png', dpi=400)    
iterations=list(range(iter_max))

Pf_t=np.real(Sf_t)
Qf_t=np.imag(Sf_t)

#path_plt = os.path.join(path_cur, 'Plot')
iplt=0
plt.figure(iplt)    
plot_vmax=plt.plot(iterations, result1.vmax,linestyle='--',color='red')
plot_vmin=plt.plot(iterations, result1.vmin,linestyle='--',color='blue')
plt.legend((plot_vmax[0], plot_vmin[0]), ('Max V', 'Min V'))
plt.xlabel('Iteration No.')
plt.ylabel('Voltage magnitude, p.u.')
plt.title('Controlled Voltage')
plt.grid(True)
#plt.savefig(path_plt+'/v.png', dpi=400)    

iplt+=1
plt.figure(iplt)    
plot_lambda_u=plt.plot(iterations, result1.lambda_u,linestyle='--',color='red')
plot_lambda_d=plt.plot(iterations, result1.lambda_d,linestyle='--',color='blue')
plt.legend((plot_lambda_u[0], plot_lambda_d[0]), ('lambda up', 'lambda lower'))
plt.xlabel('Iteration No.')
plt.ylabel('Lambda')
plt.title('Lambda')
plt.grid(True)
#plt.savefig(path_plt+'/lambda.png', dpi=400)    


iplt+=1
plt.figure(iplt)    
plot_u_u=plt.plot(iterations, result1.u_u,linestyle='--',color='red')
#plot_u_d=plt.plot(iterations, result1.u_l,linestyle='--',color='blue')
#plt.legend((plot_u_u[0], plot_u_d[0]), ('u up', 'u lower'))
#plt.legend(plot_u_u[0], 'u up')
plt.xlabel('Iteration No.')
plt.ylabel('u')
plt.title('u')
plt.grid(True)
#plt.savefig(path_plt+'/u.png', dpi=400)   


# DER optimal vs intial (p,q)
#plt.savefig(path_plt+'/Qld.png', dpi=400)  

# plt.figure(5)
# plt.plot(pll_g_t,linewidth=1)
# plt.plot(pll_g_min,'.', markersize=2)
# plt.plot(pll_g_max,'.', markersize=2)
# plt.title('real power generation')
# plt.savefig(path_plt+'/PgenvsCapacity.png', dpi=400)  

# plt.figure(6)
# plt.plot(qll_g_t,linewidth=1)
# plt.plot(qll_g_min,'.', markersize=2)
# plt.plot(qll_g_max,'.', markersize=2)
# plt.title('reactive power generation')
# plt.savefig(path_plt+'/QgenvsCapacity.png', dpi=400) 

# combine gen with sgen
pg_t=pll_g_t[busid_g_LL]
qg_t=qll_g_t[busid_g_LL]

pg_max=pll_g_max[busid_g_LL]
qg_max=qll_g_max[busid_g_LL]


# # gen optimal vs intial (p,v)
# plt.figure(5)
# plot_pg=plt.plot(pg_agr*sbase,'.',markersize=2.5)
# plot_pgt=plt.plot(pgt_agr*sbase,'--',markersize=2.5)
# plt.title('Pgen (mw)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_pgt[0], plot_pg[0]), ('Optimal', 'Initial'))
# plt.savefig(path_plt+'/Pgen.png', dpi=400) 


# plt.figure(6)
# plot_qg=plt.plot(qg_agr*sbase,markersize=2.5)
# plot_qgt=plt.plot(qgt_agr*sbase,'--',markersize=2.5)
# plt.title('Qgen (mvar)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_qgt[0], plot_qg[0]), ('Optimal', 'Initial'))
# plt.savefig(path_plt+'/Qgen.png', dpi=400)  

# gen optimal vs intial (p,v)
iplt+=1
plt.figure(iplt) 
plot_pg=plt.plot(pg,'.')
plot_pgt=plt.plot(pg_t,'.')
plt.title('Pgen (mw)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_pgt[0], plot_pg[0]), ('Optimal', 'Initial'))
plt.grid(True)
#plt.savefig(path_plt+'/Pgen.png', dpi=400)  

# plt.figure(6)
# plot_vg=plt.plot(vg,linewidth=1)
# plot_vgt=plt.plot(vg_t,'--',linewidth=1)
# plt.title('Vg (p.u.)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_vgt[0], plot_vg[0]), ('Optimal', 'Initial'))
# #plt.savefig(path_plt+'/Vg.png', dpi=400)

iplt+=1
plt.figure(iplt) 
plot_qg=plt.plot(qg,'.')
plot_qgt=plt.plot(qg_t,'.')
# plot_ub=plt.plot(qll_g_max[busid_g_LL],linewidth=1)
# plot_lb=plt.plot(qll_g_min[busid_g_LL],linewidth=1)
plt.title('Qgen (mvar)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
#plt.legend((plot_qgt[0], plot_qg[0],plot_lb[0],plot_ub[0]), ('Optimal', 'Initial','Lower Bound','Upper Bound'))
plt.legend((plot_qgt[0], plot_qg[0]), ('Optimal', 'Initial'))
plt.grid(True)
#plt.savefig(path_plt+'/Qgen.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
plot_vm=plt.plot(vmll,marker='o', markersize=0.5)
plot_vm0=plt.plot(vmll0,marker='o',markersize=0.5)
plot_ub=plt.plot([v_u]*nLL,'--',linewidth=2)
plot_lb=plt.plot([v_l]*nLL,'--',linewidth=2)
plt.ylim(vplt_min, vplt_max)
plt.title('v (p.u.)')
plt.legend(['Optimal', 'Initial','upper','lower'])
plt.grid(True)
#plt.savefig(path_plt+'/VProfile.png', dpi=400) 

iplt+=1
plt.figure(iplt) 
plot_Sf=plt.plot(range(0,nbrh),abs(Sf_t),marker='o', markersize=0.5)
plot_Sf0=plt.plot(range(0,nbrh),abs(Sf_t0),marker='o',markersize=0.5)
plot_Sf_max=plt.plot(Smax_brh,'--',color='r')
plt.ylim(0, 1.5*np.max(Smax_brh))
plt.title('Apparent Power')
plt.xlabel(['Bus index'])
plt.ylabel(['S (p.u.)'])
plt.legend(['Optimal', 'Initial','upper'])
plt.grid(True)

iplt+=1
plt.figure(iplt) 
plot_Sf=plt.plot(range(0,nbrh),abs(Pf_t),marker='o', markersize=0.5)
plot_Sf0=plt.plot(abs(Pf_t0),marker='o',markersize=0.5)
plt.ylim(0, 1.5*np.max(Smax_brh))
plt.title('Real Power')
plt.xlabel(['Bus index'])
plt.ylabel(['P (p.u.)'])
plt.legend(['Optimal', 'Initial'])
plt.grid(True)

iplt+=1
plt.figure(iplt) 
plot_Sf=plt.plot(abs(Qf_t),marker='o', markersize=0.5)
plot_Sf0=plt.plot(abs(Qf_t0),marker='o',markersize=0.5)
plt.ylim(0, 1.5*np.max(Smax_brh))
plt.title('Reactive Power')
plt.xlabel(['Bus index'])
plt.ylabel(['Q (p.u.)'])
plt.legend(['Optimal', 'Initial'])
plt.grid(True)
