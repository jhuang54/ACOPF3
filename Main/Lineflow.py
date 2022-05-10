# task: update generator and static generator, dispatch aggregated generators. All the generators on one bus are aggregated to one generator.
# Pgen, Psgen>=0, Qsgen[-18:]>=Qsgen_min[-18], Qsgen[-18:]<=Qsgen_max[-18]; No other constraints
# load is negative in linearized power flow model, but positive in net.load.p_mw, net.load.q_mvar
# generation is positive in linearized power flow model, and positive in net.gen.p_mw, net.gen.q_mvar, net.sgen. p_mw,net.sgen.q_mvar
# net.load.p_mw, net.load.q_mwar, net.gen.p_mw, net.gen.q_mvar, net.sgen. p_mw,net.sgen.q_mvar are real value with unit mw and mwar
# linearized power flow model use p.u., sbase=10mva
# dispatch net.gen (P,V), net.sgen (P,Q), and net.load (P,Q)
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

# parameter setting
iter_max=3000
alpha_ld=0.1# f=alpha*(pld-\hat{pld})^2
alpha_g=0.1# f=alpha*(pg-\hat{pg})^2
alpha_sg=0.1
beta_ld=2.5# capacity of DER=rated value*beta

beta_pgmax=1# capacity of Generator=rated value*beta
beta_pgmin=0# capacity of Generator=rated value*beta
beta_qgmax=1.5
beta_qgmin=-1.5


# convert matpower case file version 2 to a pandapower net.
# path_cur = Path(os.getcwd())
# path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
#icase = 'Maui2022dm_rd_v33.mat'# convert switched shunts into generators, one transformer has phase shift
icase='Maui2022dm_rd_v33_SwitchShuntsNoPhaseshift'# convert switched shunts into generators, transformers don't have phase shift
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition
icase= 'Maui2022dm_rd_AggregateGens.mat'# physical system simulator
net_t=pc.from_mpc(path_output + '\\' + icase, f_hz=60)


sbase=10#10 MVA

nbus = len(net.bus)
nLL=nbus-1

# run power flow to get initial condition
pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)


# Ybus
Ybus = net._ppc['internal']['Ybus']

# # power flow solution
vm = net.res_bus.vm_pu.values
va_dg = net.res_bus.va_degree.values
va = np.exp(1j * va_dg * np.pi / 180)
v = vm * va

vm0=net.ext_grid.vm_pu[0]
va0=np.exp(1j*net.ext_grid.va_degree[0]*np.pi/180)
v0=vm0*va0

# check Ybus is correct or not by calculating power injection with Ybus
mappd2ppc = net._pd2ppc_lookups["bus"]
Ybus = Ybus[mappd2ppc, 0:nbus]
Ybus = Ybus[0:nbus, mappd2ppc]
I = Ybus.dot(v)
S = v * np.conj(I)


# FPL model
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv


busid_slack=net.ext_grid.bus[0]

busid_LL=list(range(0,nbus))
busid_LL.remove(busid_slack)

# generator global bus id and local id in LL
busid_gen=list(net.gen.bus)
#busid_gen.remove(busid_slack)
busid_g_LL=np.array(busid_LL).searchsorted(busid_gen)

# load global bus id and local id in LL
busid_ld=net.load.bus
busid_ld=list(busid_ld)
nlsbus=len(busid_ld)
# map from load buse to load buses+slack bus
busid_s_ld=np.array(busid_ld).searchsorted(busid_slack)
busid_ld_lds=list(range(0,nlsbus))
busid_ld_lds.remove(busid_s_ld)

# map from load bus to index in LL
busid_ld.remove(busid_slack)

busid_ld_LL=np.array(busid_LL).searchsorted(busid_ld)

# index of sgen
busid_sg=net.sgen.bus
busid_sg_LL=np.array(busid_LL).searchsorted(busid_sg)

nsgen=len(busid_sg)
sgen_to_LL=np.zeros((nLL,nsgen), dtype=int)
for isgen in range(0,nsgen):
    sgen_to_LL[busid_sg_LL[isgen],isgen]=1

id_sgen_OutSer=np.where(net.sgen.in_service==False)

YLL=Ybus[np.ix_(busid_LL,busid_LL)]
#YLL=YLL.todense()
Z_pu=inv(YLL)
Z_pu=Z_pu.todense()
Zt=Z_pu/np.conj(v0)

YLS=Ybus[np.ix_(busid_LL,[busid_slack])]
YLS=YLS.todense()

#YY=np.multiply(Z_pu,YLS)
#YY=Z_pu.dot(YLS)
YY=np.matmul(Z_pu,YLS)
vs=YY*v0
#vs=vs.todense()
vs=np.squeeze(np.asarray(vs))

# My=[Zt -1i*Zt]
from scipy.sparse import hstack
# My=hstack((Zt, -1j*Zt))
# Ky=1/vm0*np.real(np.conj(v0)*My)
#const_vm_FPL=np.full((nbus, 1), vm0)
w=-vs
Wm=np.diag(w)
Wm_inv=np.diag(1/w)
Kyr=np.matmul(np.matmul(np.matmul(abs(Wm),Wm_inv),Z_pu),np.conj(Wm_inv))
#R=np.real(1/vm0*np.conj(v0)*Zt)
R=np.real(Kyr)
#Rt=csc_matrix.transpose(R)
Rt=np.transpose(R)
Rt=Rt

Rsgen=R[np.ix_(range(0,nLL),busid_sg_LL)]# update sgen

#X=np.real(1/vm0*np.conj(v0)*(-1j)*Zt)
X=np.real(-1j*Kyr)
#Xt=csc_matrix.transpose(X)
Xt=np.transpose(X)
Xt=Xt

Xsgen=X[np.ix_(range(0,nLL),busid_sg_LL)]# update sgen


# apparent power flow limit model
nbrh=len(net.trafo)+len(net.line)

# branch table [fbus,tbus]
branch_tb=np.zeros((nbrh,2))
brh_fbus=np.concatenate((net.line.from_bus.values,net.trafo.hv_bus.values))
brh_tbus=np.concatenate((net.line.to_bus.values,net.trafo.lv_bus.values))
brh_serv=pd.concat([net.line.in_service, net.trafo.in_service])
id_brh_serv=np.where(brh_serv==True)[0]
nbrh=len(id_brh_serv)
brh_fbus=brh_fbus[id_brh_serv]
brh_tbus=brh_tbus[id_brh_serv]

# bus voltage 
vn=net.bus.vn_kv.to_numpy()
line_fbus=net.line.from_bus.values
id_line_serv=np.where(net.line.in_service==True)[0]
line_fbus=line_fbus[id_line_serv]
vn_line_fbus=vn[line_fbus]

brh_y=np.zeros(nbrh,dtype=complex)
# branch impedance
for i in range(nbrh):
    brh_y[i]=-Ybus[brh_fbus[i],brh_tbus[i]]
brh_z=1/brh_y
# branch_tb[:,0]=fbus
# branch_tb[:,1]=tbus

# Afti
A=np.zeros((nLL,nbrh))
for ibus in range(nLL):
    busid_tp=busid_LL[ibus]
    idbrh_f=np.where(brh_fbus==busid_tp)#branch id of branches whose fbus is busid_tp
    idbrh_t=np.where(brh_tbus==busid_tp)
    A[ibus,idbrh_f]=1
    A[ibus,idbrh_t]=-1

# select a tree in a meshed network 
bus_pool=[busid_slack]
Nextbus0=[busid_slack]
Nextbus=np.array([],dtype='int64').reshape(0)
brh_l=np.array([],dtype='int64').reshape(0)# lian branch
brh_t=np.array([],dtype='int64').reshape(0)# tree branch

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
for i in range(nLL):
    brh_id_tp=id_brh_tb[out_ind[i,0],out_ind[i,1]]
    brh_t=np.concatenate((brh_t,[brh_id_tp]))


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

# # available load
# pld=net.load.p_mw.to_numpy()
# qld=net.load.q_mvar.to_numpy()


# initialize controllable loads
epsi_pld=0
epsi_qld=0
epsi_pg=0.05
epsi_qg=0.05
epsi_psg=0.05
epsi_qsg=0.05
# epsi_pg=0
# epsi_qg=0
# epsi_psg=0
epsi_l=0.1

epsi_u=0.1

v_l=0.95# bounds for opf
v_u=1.05
vplt_min=0.88# bounds for plot
vplt_max=1.07

Sfm_u=10
Sfm_l=0


# initial values of load 
pld_t=-net.load.p_mw.to_numpy()/sbase
pld_t=np.delete(pld_t, busid_s_ld)
qld_t=-net.load.q_mvar.to_numpy()/sbase
qld_t=np.delete(qld_t, busid_s_ld)

pll_ld_t=np.zeros(nLL)
pll_ld_t[busid_ld_LL]=pld_t# load is negative in linearized power flow model
qll_ld_t=np.zeros(nLL)
qll_ld_t[busid_ld_LL]=qld_t

# desired values of load
pld=-net.load.p_mw.to_numpy()/sbase
pld=np.delete(pld,busid_s_ld)
qld=-net.load.q_mvar.to_numpy()/sbase
qld=np.delete(qld,busid_s_ld)

pll_ld=np.zeros(nLL)
pll_ld[busid_ld_LL]=pld

qll_ld=np.zeros(nLL)
qll_ld[busid_ld_LL]=qld

# capacity of DER
pll_ld_max=np.zeros(nLL)
pll_ld_max[busid_ld_LL]=abs(pld)*beta_ld
pll_ld_min=np.zeros(nLL)
pll_ld_min[busid_ld_LL]=-abs(pld)*beta_ld

qll_ld_max=np.zeros(nLL)
qll_ld_max[busid_ld_LL]=abs(qld)*beta_ld
qll_ld_min=np.zeros(nLL)
qll_ld_min[busid_ld_LL]=-abs(qld)*beta_ld

# initial values of generator
pg_t=net.res_gen.p_mw.to_numpy()/sbase
qg_t=net.res_gen.q_mvar.to_numpy()/sbase

pll_g_t=np.zeros(nLL)
pll_g_t[busid_g_LL]=pg_t# load is negative in linearized power flow model
qll_g_t=np.zeros(nLL)
qll_g_t[busid_g_LL]=qg_t

# desired values of gen
pg=net.res_gen.p_mw.to_numpy()/sbase
qg=net.res_gen.q_mvar.to_numpy()/sbase
vg=net.res_gen.vm_pu

pll_g=np.zeros(nLL)
pll_g[busid_g_LL]=pg# load is negative in linearized power flow model
qll_g=np.zeros(nLL)
qll_g[busid_g_LL]=qg

# capacity of Generator
pll_g_max=np.zeros(nLL)
pll_g_max[busid_g_LL]=net.gen.max_p_mw.to_numpy()/sbase
pll_g_min=np.zeros(nLL)
pll_g_min[busid_g_LL]=net.gen.min_p_mw.to_numpy()/sbase

qll_g_max=np.zeros(nLL)
qll_g_max[busid_g_LL]=net.gen.max_q_mvar.to_numpy()/sbase
qll_g_min=np.zeros(nLL)
qll_g_min[busid_g_LL]=net.gen.min_q_mvar.to_numpy()/sbase


# initial values of sgen
psg_t=net.res_sgen.p_mw.to_numpy()/sbase
qsg_t=net.res_sgen.q_mvar.to_numpy()/sbase
nsg=len(qsg_t)
# epsi_qsg=np.zeros(nsg)
# epsi_qsg[-18:]=0.05

# desired values of sgen
psg=net.res_sgen.p_mw.to_numpy()/sbase
qsg=net.res_sgen.q_mvar.to_numpy()/sbase

# capacity of sGenerator
psg_max=net.sgen.max_p_mw.to_numpy()/sbase
psg_max[id_sgen_OutSer]=0#inactive generators which are out of service
psg_min=net.sgen.min_p_mw.to_numpy()/sbase
psg_min[id_sgen_OutSer]=0#inactive generators which are out of service


qsg_max=net.sgen.max_q_mvar.to_numpy()/sbase
qsg_max[id_sgen_OutSer]=0#inactive generators which are out of service
qsg_min=net.sgen.min_q_mvar.to_numpy()/sbase
#qsg_min[-18:]=-net.sgen.max_q_mvar.to_numpy()[-18:]/sbase
qsg_min[id_sgen_OutSer]=0#inactive generators which are out of service


# remove slack bus
vmll=vm[busid_LL]
vall=va[busid_LL]
vll=v[busid_LL]

vmll0=vm[busid_LL]

iplt=1
plt.figure(iplt)
plot_vm=plt.plot(vmll,'.')
plot_ub=plt.plot([v_u]*nLL,'--',linewidth=1)
plot_lb=plt.plot([v_l]*nLL,'--',linewidth=1)
plt.ylim(vplt_min, vplt_max)
plt.title('v (p.u.)')
plt.savefig(path_plt+'/VProfile0.png', dpi=400) 

print('max v:',vmll.max())
print('min v:',vmll.min())


# initial net injections of LL nodes: in panderpower, load is positive
pll_net_t=np.real(S[busid_LL])# DON'T use net.res_bus.p_mw[busid_LL].to_numpy()/sbase, because capacitors are included in Ybus matrix
qll_net_t=np.imag(S[busid_LL])
sll_net_t=S[busid_LL]


lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)

vll0=vll


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
result1=result(iter_max)  

# complex voltage of all the buses including slack bus
v_t=np.ones(nbus,dtype=complex)

Qfr=np.zeros(nbrh)
Pfr=np.zeros(nbrh)
u_u=np.zeros(nbrh)
u_l=np.zeros(nbrh)
for iter in range(iter_max):
    # derivative of voltage constraints with respect to (p, q)
    dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
    dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
    dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
    dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
    
    dvcnstr_dsp=dvcnstr_dp[np.ix_(busid_sg_LL)]
    dvcnstr_dsq=dvcnstr_dq[np.ix_(busid_sg_LL)]
    
    # derivative of apparent power flow with respect to (p_g,p_pv,p_battery,p_load)
    Qfrdu=Qfr*(u_u-u_l)
    Pfrdu=Pfr*(u_u-u_l)
    
    dsf_dp=Dt.dot(Qfrdu)+Ct.dot(Pfrdu)
    dsf_dp=np.squeeze(np.array(dsf_dp))
    dsf_dq=Ct.dot(Qfrdu)-Dt.dot(Pfrdu)
    dsf_dq=np.squeeze(np.array(dsf_dq))
    
    dsf_dsp=dsf_dp[np.ix_(busid_sg_LL)]
    dsf_dsq=dsf_dq[np.ix_(busid_sg_LL)]
    
    # each bus has at most 1 load and generator, but may have multiple static generators
    # minimize deviation from (pll_ld,qll_ld)
    pll_ld_t=pll_ld_t-epsi_pld*(2*alpha_ld*(pll_ld_t-pll_ld)+dvcnstr_dp+dsf_dp)
    qll_ld_t=qll_ld_t-epsi_qld*(2*alpha_ld*(qll_ld_t-qll_ld)+dvcnstr_dq+dsf_dq)
    
    # minimize generation from coal generator
    pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)+dvcnstr_dp+dsf_dp)
    qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)+dvcnstr_dq+dsf_dq)
    
    # minimize sgeneration from coal sgenerator 
    psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp+dsf_dsp)
    qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq+dsf_dsq)
    
    
    # project
    # pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
    # pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
    # qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
    # qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
    
    pll_g_t=np.maximum(pll_g_t,pll_g_min)
    # pll_g_t=np.minimum(pll_g_t,pll_g_max)
    # qll_g_t=np.maximum(qll_g_t,qll_g_min)
    # qll_g_t=np.minimum(qll_g_t,qll_g_max)
    
    psg_t=np.maximum(psg_t,psg_min)
    # psg_t=np.minimum(psg_t,psg_max)
    qsg_t[-18:]=np.maximum(qsg_t[-18:],qsg_min[-18:])
    qsg_t[-18:]=np.minimum(qsg_t[-18:],qsg_max[-18:])
    
    # sgen to sll
    pll_sg_t=np.matmul(sgen_to_LL,psg_t)
    qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
    
    # update net injections 
    pll_net_t=pll_ld_t+pll_g_t+pll_sg_t
    qll_net_t=qll_ld_t+qll_g_t+qll_sg_t 
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
               
    if err_vll>1e-3:
        raise Exception('power flow diverge!\n')
    
         
    # # dispatch load 
    # pld_t=pll_ld_t[busid_ld_LL]
    # qld_t=qll_ld_t[busid_ld_LL]
    
    # net_t.load.p_mw[busid_ld_lds]=-pld_t*sbase# load is positive in net.load
    # net_t.load.q_mvar[busid_ld_lds]=-qld_t*sbase
    
    # # dispatch (p,v) for generator
    # pAg_t=pAg_net_t[busid_g_LL]
    # vg_t=abs(vll[busid_g_LL])
    
    # net_t.gen.p_mw=pAg_t*sbase
    # net_t.gen.vm_pu=vg_t   
    
    # pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)

    # # check whether generator tracks dispatch 
    # # load
    # pld_t_m=-net_t.res_load.p_mw.to_numpy()/sbase
    # pld_t_m=np.delete(pld_t_m, busid_s_ld)
    # qld_t_m=-net_t.res_load.q_mvar.to_numpy()/sbase
    # qld_t_m=np.delete(qld_t_m, busid_s_ld)
    
    # # generator
    # pg_t_m=net_t.res_gen.p_mw.to_numpy()/sbase
    # # qg_t_m=net_t.res_gen.q_mvar/sbase
    # vg_t_m=net_t.res_gen.vm_pu.to_numpy()
    
    # # sgenerator
    # # psg_t_m=net_t.res_sgen.p_mw.to_numpy()/sbase
    # # # qg_t_m=net_t.res_gen.q_mvar/sbase
    # # qsg_t_m=net_t.res_sgen.q_mvar.to_numpy()/sbase
    
    # # mismatch between dispatch and measurements
    # dpld=abs(pld_t_m-pld_t)
    # dqld=abs(qld_t_m-qld_t)
    # dpg=abs(pg_t_m-pAg_t)
    # dvg=abs(vg_t_m-vg_t)
    
    # if dpld.max()>1e-2 or dqld.max()>1e-2 or dpg.max()>1e-3 or dvg.max()>2e-3:
    #     print('iteration:%d: large mismatch' %iter)
        
    # print('maximum mismatch:')
    # print('pld:%.8f' %dpld.max())
    # print('qld:%.8f' %dqld.max())
    # print('pg:%.8f' %dpg.max())
    # print('v:%.8f' %dvg.max())
    
    # print('predicted:')
    # print('max v:',abs(vll).max())
    # print('min v:',abs(vll).min())
    
    # vmll=net_t.res_bus.vm_pu.values[busid_LL]
    
    # update dual variables (voltage)
    vmll=abs(vll)
    lambda_u=lambda_u+epsi_l*(vmll-v_u)
    lambda_d=lambda_d+epsi_l*(v_l-vmll)
    
    # project dual variables
    lambda_u=np.maximum(lambda_u,0)
    lambda_d=np.maximum(lambda_d,0)
    
    # track voltage and dual
    result1.lambda_u[iter]=lambda_u.max()
    result1.lambda_d[iter]=lambda_d.max()
    result1.vmax[iter]=vmll.max()
    result1.vmin[iter]=vmll.min()
    
    #update dual variable (apparent power flow)
    # nonlinear complex power flow
    # v_t[busid_slack]=v0
    # v_t[busid_LL]=vll
    # vf_t=v_t[brh_fbus]
    # vt_t=v_t[brh_tbus]
    # If_t=(vf_t-vt_t)*brh_y
    # Sf_t=vf_t*np.conjugate(If_t)
    SfLine_t=net.res_line.p_from_mw.to_numpy()+1j*net.res_line.q_from_mvar.to_numpy()
    SfTran_t=net.res_trafo.p_hv_mw.to_numpy()+1j*net.res_trafo.q_hv_mvar.to_numpy()
    Sf_t=np.concatenate((SfLine_t,SfTran_t))
    Sf_t=Sf_t[id_brh_serv]/sbase
    
    # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    # Pi_t=np.real(sll_net_t)
    # Qi_t=np.imag(sll_net_t)
    # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    
    # dpf=np.real(Sf_t)-Pf_e
    # dqf=np.imag(Sf_t)-Qf_e
    
    Sfm_t=abs(Sf_t)
    # Pfr=np.real(Sf_t)/Sfm_t
    # Qfr=np.imag(Sf_t)/Sfm_t
    
    # # inactive small branch flow
    # id_sf0=np.where(Sfm_t<1e-7)
    # Pfr[id_sf0]=0
    # Qfr[id_sf0]=0
    
    # # dual variable
    # u_u=u_u+epsi_u*(Sfm_t-Sfm_u)
    # u_l=u_l+epsi_u*(Sfm_l-Sfm_t)
    
    # # project dual variables
    # u_u=np.maximum(u_u,0)
    # u_l=np.maximum(u_l,0)
    
    # accuracy
    
    Pi_t=np.real(sll_net_t)
    Qi_t=np.imag(sll_net_t)
    #Sfe_t=(C.dot(Pi_t)-D.dot(Qi_t))+1j*(D.dot(Pi_t)+C.dot(Qi_t))  
    Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    Sfe_t=Pf_e+1j*Qf_e
    
    result1.dSfm_t[iter,:]=abs(Sfm_t-abs(Sfe_t))
    id_inac=np.where(Sfm_t<1e-4)
    result1.dSfmp_t[iter,:]=result1.dSfm_t[iter,:]/Sfm_t
    result1.dSfmp_t[iter,id_inac]=0
    result1.dPfm_t[iter,:]=abs(np.real(Sf_t)-Pf_e)
    result1.dPfmp_t[iter,:]=abs(np.real(Sf_t)-Pf_e)/abs(np.real(Sf_t))
    id_inac=np.where(abs(np.real(Sf_t))<1e-4)
    result1.dPfm_t[iter,id_inac]=0
    result1.dPfmp_t[iter,id_inac]=0
    plt.figure(1)
    plt.plot(result1.dPfm_t[iter,:])
    plt.title('Absolute error of Pfm')
    plt.grid(True)
    plt.figure(2)
    plt.plot(result1.dPfmp_t[iter,:])
    plt.title('Relative error of Pfm')
    plt.grid(True)
iterations=list(range(iter_max))

#path_plt = os.path.join(path_cur, 'Plot')
iplt+=1
plt.figure(iplt)    
plot_vmax=plt.plot(iterations, result1.vmax,linestyle='--',color='red')
plot_vmin=plt.plot(iterations, result1.vmin,linestyle='--',color='blue')
plt.legend((plot_vmax[0], plot_vmin[0]), ('Max V', 'Min V'))
plt.xlabel('Iteration No.')
plt.ylabel('Voltage magnitude, p.u.')
plt.title('Controlled Voltage')
plt.grid(True)
plt.savefig(path_plt+'/v.png', dpi=400)    

iplt+=1
plt.figure(iplt)    
plot_lambda_u=plt.plot(iterations, result1.lambda_u,linestyle='--',color='red')
plot_lambda_d=plt.plot(iterations, result1.lambda_d,linestyle='--',color='blue')
plt.legend((plot_lambda_u[0], plot_lambda_d[0]), ('lambda up', 'lambda lower'))
plt.xlabel('Iteration No.')
plt.ylabel('Lambda')
plt.title('Lambda')
plt.grid(True)
plt.savefig(path_plt+'/lambda.png', dpi=400)    


# DER optimal vs intial (p,q)
iplt+=1
plt.figure(iplt) 
plot_pldt=plt.plot(-pld_t,'.')
plot_pld=plt.plot(-pld,'.')
plt.title('real load (mw)')
plt.legend((plot_pldt[0], plot_pld[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Pld.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
plot_qldt=plt.plot(-qld_t,'.')
plot_qld=plt.plot(-qld,'*')
plt.title('reactive load (mvar)')
plt.legend((plot_qldt[0], plot_qld[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Qld.png', dpi=400)  

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

pg_agr=np.concatenate((pg, psg))
pgt_agr=np.concatenate((pg_t, psg_t))

qg_agr=np.concatenate((qg, qsg))
qgt_agr=np.concatenate((qg_t, qsg_t))

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
plt.savefig(path_plt+'/Pgen.png', dpi=400)  

# plt.figure(6)
# plot_vg=plt.plot(vg,linewidth=1)
# plot_vgt=plt.plot(vg_t,'--',linewidth=1)
# plt.title('Vg (p.u.)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_vgt[0], plot_vg[0]), ('Optimal', 'Initial'))
# plt.savefig(path_plt+'/Vg.png', dpi=400)

iplt+=1
plt.figure(iplt) 
plot_qg=plt.plot(qg,'.')
plot_qgt=plt.plot(qg_t,'.')
plot_ub=plt.plot(qll_g_max[busid_g_LL],linewidth=1)
plot_lb=plt.plot(qll_g_min[busid_g_LL],linewidth=1)
plt.title('Qgen (mvar)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_qgt[0], plot_qg[0],plot_lb[0],plot_ub[0]), ('Optimal', 'Initial','Lower Bound','Upper Bound'))
plt.savefig(path_plt+'/Qgen.png', dpi=400)  


# sgen optimal vs intial (p,q)
iplt+=1
plt.figure(iplt) 
# plot_psg=plt.plot(psg,linewidth=1,'*')
# plot_psgt=plt.plot(psg_t,linewidth=1,'*')
# plot_psgmin=plt.plot(psg_min,'.', markersize=2)
# plot_psgmax=plt.plot(psg_max,'.', markersize=2)
plot_psg=plt.plot(psg,'.')
plot_psgt=plt.plot(psg_t,'.')
plt.title('Psgen (mw)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_psgt[0], plot_psg[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Psgen.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
plot_qsg=plt.plot(qsg,'.')
plot_qsgt=plt.plot(qsg_t,'.')
plot_ub=plt.plot(qsg_max,linewidth=1)
plot_lb=plt.plot(qsg_min,linewidth=1)
plt.title('Qsgen (mvar)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_qsgt[0], plot_qsg[0],plot_lb[0],plot_ub[0]), ('Optimal', 'Initial','Lower Bound','Upper Bound'))
plt.savefig(path_plt+'/Qsgen.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
plot_vm=plt.plot(vmll,marker='o', markersize=0.5)
plot_vm0=plt.plot(vmll0,marker='o',markersize=0.5)
plot_ub=plt.plot([v_u]*nLL,'--',linewidth=1)
plot_lb=plt.plot([v_l]*nLL,'--',linewidth=1)
plt.ylim(vplt_min, vplt_max)
plt.title('v (p.u.)')
plt.legend(['Optimal', 'Initial','upper','lower'])
plt.savefig(path_plt+'/VProfile.png', dpi=400) 