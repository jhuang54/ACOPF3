# task: update load, generator and static generator, dispatch load and aggregated generators. All the generators on one bus are aggregated to one generator.
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

# parameter setting
iter_max=1100
alpha_ld=1# f=alpha*(pld-\hat{pld})^2
alpha_g=0.1# f=alpha*(pld-\hat{pld})^2
alpha_sg=0.1
beta_ld=2# capacity of DER=rated value*beta
beta_g=1.5# capacity of Generator=rated value*beta
beta_sg=1.5


# convert matpower case file version 2 to a pandapower net.
path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition
#net_t=pc.from_mpc(path_output + '\\' + icase, f_hz=60)
icase= 'Maui2022dm_rd_AggregateGens.mat'# physical system simulator
net_t=pc.from_mpc(path_output + '\\' + icase, f_hz=60)
# icase = 'Maui2022dm_rd_NoPV.mat'
# net_sim=pc.from_mpc(path_output + '\\' + icase, f_hz=60)

sbase=10#10 MVA

nbus = len(net.bus)
nLL=nbus-1

# run power flow to get initial condition
pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)

#pp.runpp(net_sim, algorithm='nr', calculate_voltage_angles=True)

# Ybus
Ybus = net._ppc['internal']['Ybus']

#import math
#import cmath
import numpy as np

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

# # available load
# pld=net.load.p_mw.to_numpy()
# qld=net.load.q_mvar.to_numpy()


# initialize controllable loads
epsi_pld=0.05
epsi_qld=0.05
epsi_pg=0.05
epsi_qg=0.05
epsi_psg=0.05
epsi_qsg=0.05
epsi_l=0.15

v_l=0.95
v_u=1.05

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
pll_g_max[busid_g_LL]=abs(pg)*beta_g
pll_g_min=np.zeros(nLL)

qll_g_max=np.zeros(nLL)
qll_g_max[busid_g_LL]=abs(qg)*beta_g
qll_g_min=np.zeros(nLL)
qll_g_min[busid_g_LL]=-abs(qg)*beta_g


# initial values of sgen
psg_t=net.res_sgen.p_mw.to_numpy()/sbase
qsg_t=net.res_sgen.q_mvar.to_numpy()/sbase

# desired values of sgen
psg=net.res_sgen.p_mw.to_numpy()/sbase
qsg=net.res_sgen.q_mvar.to_numpy()/sbase

# capacity of sGenerator
psg_max=abs(psg)*beta_sg
psg_max[id_sgen_OutSer]=0#inactive generators which are out of service
psg_min=np.zeros(nsgen)
psg_min[id_sgen_OutSer]=0#inactive generators which are out of service


qsg_max=abs(qsg)*beta_sg
qsg_max[id_sgen_OutSer]=0#inactive generators which are out of service
qsg_min=-abs(qsg)*beta_sg
qsg_min[id_sgen_OutSer]=0#inactive generators which are out of service

# power flow solution
# vm = net.res_bus.vm_pu.values
# va_dg = net.res_bus.va_degree.values
# va = np.exp(1j * va_dg * np.pi / 180)
# v = vm * va

# slack bus
# vs=vm[busid_slack]

# remove slack bus
vmll=vm[busid_LL]
vall=va[busid_LL]
vll=v[busid_LL]

print('max v:',vmll.max())
print('min v:',vmll.min())


# initial net injections of LL nodes: in panderpower, load is positive
# pll_net_t=-net.res_bus.p_mw[busid_LL].to_numpy()/sbase
# qll_net_t=-net.res_bus.q_mvar[busid_LL].to_numpy()/sbase

pll_net_t=np.real(S[busid_LL])# DON'T use net.res_bus.p_mw[busid_LL].to_numpy()/sbase, because capacitors are included in Ybus matrix
qll_net_t=np.imag(S[busid_LL])
sll_net_t=S[busid_LL]

# psgen_net_t=pll_net_t-pll_ld_t-pll_g_t
# qsgen_net_t=qll_net_t-qll_ld_t-qll_g_t

lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)

# vll0=vll
# itr_pf=0
# itr_pf_max=20
# err_vll=1
# while (err_vll>1e-5 and itr_pf<itr_pf_max):
#     vll=Z_pu.dot(np.conj(sll_net_t/vll))+vs
#     err_vll=max(abs(vll-vll0))
#     vll0=vll
#     itr_pf=itr_pf+1
# plt.figure(1)
# plt.plot(vmll)
# plt.plot(abs(vll))
# update controllable loads (pld_t,qld_t) and dual variables

vll0=vll

# itr_pf=0
# itr_pf_max=20
# err_vll=1
# while (err_vll>1e-5 and itr_pf<itr_pf_max):
    # #vll=Z_pu*np.conj(sll_net_t/vll)+vs
    # ILL=np.conj(sll_net_t/vll)
    # dv=np.matmul(Z_pu,ILL.transpose())
    # vll=np.squeeze(np.asarray(dv))+w
    # err_vll=max(abs(vll-vll0))
    # vll0=vll
    # itr_pf=itr_pf+1
           
# if err_vll>1e-3:
    # raise Exception('power flow diverge!\n')

class result:
    def __init__(self, iter_max):
        self.vmax=np.full((iter_max,1),1.05)
        self.vmin=np.full((iter_max,1),0.95)
        self.lambda_u=np.zeros(iter_max)
        self.lambda_d=np.zeros(iter_max)
result1=result(iter_max)  
for iter in range(iter_max):
    # dpll=pll_t-pll
    # dqll=qll_t-qll
    # print('max dpll: \n',dpll[busid_ld_LL].max())
    # print('min dpll: \n',dpll[busid_ld_LL].min())
    # print('max dqll: \n',dqll[busid_ld_LL].max())
    # print('min dqll: \n',dqll[busid_ld_LL].min())
    
    # net injection excluding DER and generators
    # pll_net_t1=pll_net_t
    # qll_net_t1=qll_net_t
    # pll_net_t=pll_net_t-pll_ld_t-pll_g_t
    # qll_net_t=qll_net_t-qll_ld_t-qll_g_t
    
    # derivative of voltage constraints with respect to (p, q)
    dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
    dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
    dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
    dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
    
    dvcnstr_dsp=dvcnstr_dp[np.ix_(busid_sg_LL)]
    dvcnstr_dsq=dvcnstr_dq[np.ix_(busid_sg_LL)]
    # dvcnstr_dp=np.matmul(Rt,(lambda_u-lambda_d))
    # dvcnstr_dq=np.matmul(Xt,(lambda_u-lambda_d))
    
    # minimize deviation from (pll_ld,qll_ld)
    pll_ld_t=pll_ld_t-epsi_pld*(2*alpha_ld*(pll_ld_t-pll_ld)+dvcnstr_dp)
    qll_ld_t=qll_ld_t-epsi_qld*(2*alpha_ld*(qll_ld_t-qll_ld)+dvcnstr_dq)
    
    # minimize generation from coal generator
    pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)+dvcnstr_dp)
    qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)+dvcnstr_dq)
    
    # minimize sgeneration from coal sgenerator 
    psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp)
    qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq)
    
    # lambda_u=lambda_u+epsi_l*(vmll-v_u)
    # lambda_d=lambda_d+epsi_l*(v_l-vmll)
    
    # project
    pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
    pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
    qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
    qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
    
    pll_g_t=np.maximum(pll_g_t,pll_g_min)
    pll_g_t=np.minimum(pll_g_t,pll_g_max)
    qll_g_t=np.maximum(qll_g_t,qll_g_min)
    qll_g_t=np.minimum(qll_g_t,qll_g_max)
    
    psg_t=np.maximum(psg_t,psg_min)
    psg_t=np.minimum(psg_t,psg_max)
    qsg_t=np.maximum(qsg_t,qsg_min)
    qsg_t=np.minimum(qsg_t,qsg_max)
    
    # sgen to sll
    pll_sg_t=np.matmul(sgen_to_LL,psg_t)
    qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
    
    # update net injections
    # pll_net_t=psgen_net_t+pll_ld_t+pll_g_t
    # qll_net_t=qsgen_net_t+qll_ld_t+qll_g_t 
    
    pll_net_t=pll_ld_t+pll_g_t+pll_sg_t
    qll_net_t=qll_ld_t+qll_g_t+qll_sg_t 
    sll_net_t=pll_net_t+1j*qll_net_t
    
    pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
    

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

         
    # dispatch load 
    pld_t=pll_ld_t[busid_ld_LL]
    qld_t=qll_ld_t[busid_ld_LL]
    
    net_t.load.p_mw[busid_ld_lds]=-pld_t*sbase# load is positive in net.load
    net_t.load.q_mvar[busid_ld_lds]=-qld_t*sbase
    
    # dispatch (p,v) for generator
    pAg_t=pAg_net_t[busid_g_LL]
    vg_t=abs(vll[busid_g_LL])
    
    net_t.gen.p_mw=pAg_t*sbase
    net_t.gen.vm_pu=vg_t
    
    # dispatch (p,q) for sgen
    # net_t.sgen.p_mw=psg_t*sbase # from p.u to real value
    # net_t.sgen.q_mvar=qsg_t*sbase
    
    pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)

    # check whether generator tracks dispatch 
    # load
    pld_t_m=-net_t.res_load.p_mw.to_numpy()/sbase
    pld_t_m=np.delete(pld_t_m, busid_s_ld)
    qld_t_m=-net_t.res_load.q_mvar.to_numpy()/sbase
    qld_t_m=np.delete(qld_t_m, busid_s_ld)
    
    # generator
    pg_t_m=net_t.res_gen.p_mw.to_numpy()/sbase
    # qg_t_m=net_t.res_gen.q_mvar/sbase
    vg_t_m=net_t.res_gen.vm_pu.to_numpy()
    
    # sgenerator
    # psg_t_m=net_t.res_sgen.p_mw.to_numpy()/sbase
    # # qg_t_m=net_t.res_gen.q_mvar/sbase
    # qsg_t_m=net_t.res_sgen.q_mvar.to_numpy()/sbase
    
    # mismatch between dispatch and measurements
    dpld=abs(pld_t_m-pld_t)
    dqld=abs(qld_t_m-qld_t)
    dpg=abs(pg_t_m-pAg_t)
    dvg=abs(vg_t_m-vg_t)
    # dpsg=abs(psg_t_m-psg_t)
    # dqsg=abs(qsg_t_m-qsg_t)
    
    if dpld.max()>1e-2 or dqld.max()>1e-2 or dpg.max()>1e-3 or dvg.max()>2e-3:
        print('iteration:%d: large mismatch' %iter)
        
    # print('maximum mismatch:')
    # print('pld:%.8f' %dpld.max())
    # print('qld:%.8f' %dqld.max())
    # print('pg:%.8f' %dpg.max())
    # print('v:%.8f' %dvg.max())
    
    # print('predicted:')
    # print('max v:',abs(vll).max())
    # print('min v:',abs(vll).min())
    
    vmll=net_t.res_bus.vm_pu.values[busid_LL]
    
    # update dual variables
    lambda_u=lambda_u+epsi_l*(vmll-v_u)
    lambda_d=lambda_d+epsi_l*(v_l-vmll)
    
    # project dual variables
    lambda_u=np.maximum(lambda_u,0)
    lambda_d=np.maximum(lambda_d,0)
    
    # plt.figure(1)
    # plt.plot(abs(vll))
    # plt.plot(vmll)
    # plt.show()
    
    # print('after dispatch:')
    # print('max v:',vmll.max())
    # print('min v:',vmll.min())
    
    # track voltage and dual
    result1.lambda_u[iter]=lambda_u.max()
    result1.lambda_d[iter]=lambda_d.max()
    result1.vmax[iter]=vmll.max()
    result1.vmin[iter]=vmll.min()
iterations=list(range(iter_max))

path_plt = os.path.join(path_cur, 'Plot')
plt.figure(1)    
plot_vmax=plt.plot(iterations, result1.vmax,linestyle='--',color='red')
plot_vmin=plt.plot(iterations, result1.vmin,linestyle='--',color='blue')
plt.legend((plot_vmax[0], plot_vmin[0]), ('Max V', 'Min V'))
plt.xlabel('Iteration No.')
plt.ylabel('Voltage magnitude, p.u.')
plt.title('Controlled Voltage')
plt.grid(True)
plt.savefig(path_plt+'/v.png', dpi=400)    

plt.figure(2)    
plot_lambda_u=plt.plot(iterations, result1.lambda_u,linestyle='--',color='red')
plot_lambda_d=plt.plot(iterations, result1.lambda_d,linestyle='--',color='blue')
plt.legend((plot_lambda_u[0], plot_lambda_d[0]), ('lambda up', 'lambda lower'))
plt.xlabel('Iteration No.')
plt.ylabel('Lambda')
plt.title('Lambda')
plt.grid(True)
plt.savefig(path_plt+'/lambda.png', dpi=400)    

# # capacity of generators
# plt.figure(3)
# plt.plot(-pll_ld_t,linewidth=1)
# plt.plot(pll_ld_min,'.', markersize=2)
# plt.plot(pll_ld_max,'.', markersize=2)
# plt.title('real load')
# plt.savefig(path_plt+'/PldvsCapacity.png', dpi=400)  

# plt.figure(4)
# plt.plot(-qll_ld_t,linewidth=1)
# plt.plot(qll_ld_min,'.', markersize=2)
# plt.plot(qll_ld_max,'.', markersize=2)
# plt.title('reactive load')
# plt.savefig(path_plt+'/QldvsCapacity.png', dpi=400)  

# DER optimal vs intial (p,q)
plt.figure(3)
plot_pldt=plt.plot(-pld_t*sbase,linewidth=1)
plot_pld=plt.plot(-pld*sbase,linewidth=1)
plt.title('real load (mw)')
plt.legend((plot_pldt[0], plot_pld[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Pld.png', dpi=400)  

plt.figure(4)
plot_qldt=plt.plot(-qld_t*sbase,linewidth=1)
plot_qld=plt.plot(-qld*sbase,linewidth=1)
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

# gen optimal vs intial (p,v)
plt.figure(5)
plot_pg=plt.plot(pg*sbase,'.',markersize=2.5)
plot_pgt=plt.plot(pg_t*sbase,'.',markersize=2.5)
plt.title('Pgen (mw)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_pgt[0], plot_pg[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Pgen.png', dpi=400)  

plt.figure(6)
plot_vg=plt.plot(vg,linewidth=1)
plot_vgt=plt.plot(vg_t,'--',linewidth=1)
plt.title('Vg (p.u.)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_vgt[0], plot_vg[0]), ('Optimal', 'Initial'))
plt.savefig(path_plt+'/Vg.png', dpi=400)

# # sgen optimal vs intial (p,q)
# plt.figure(7)
# # plot_psg=plt.plot(psg,linewidth=1,'*')
# # plot_psgt=plt.plot(psg_t,linewidth=1,'*')
# # plot_psgmin=plt.plot(psg_min,'.', markersize=2)
# # plot_psgmax=plt.plot(psg_max,'.', markersize=2)
# plot_psg=plt.plot(psg*sbase,'.',markersize=2.5)
# plot_psgt=plt.plot(psg_t*sbase,'.',markersize=2.5)
# plt.title('Psgen (mw)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_psgt[0], plot_psg[0]), ('Optimal', 'Initial'))
# plt.savefig(path_plt+'/Psgen.png', dpi=400)  

# plt.figure(8)
# plot_qsg=plt.plot(qsg*sbase,linewidth=1)
# plot_qsgt=plt.plot(qsg_t*sbase,'--',linewidth=1)
# plt.title('Qsgen (mvar)')
# #plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
# plt.legend((plot_qsgt[0], plot_qsg[0]), ('Optimal', 'Initial'))
# plt.savefig(path_plt+'/Qsgen.png', dpi=400)  