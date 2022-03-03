# load is negative in linearized power flow model, but positive in net.load
import pandapower as pp
#import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

import matplotlib.pyplot as plt

# parameter setting
iter_max=200
alpha_ld=1# f=alpha*(pld-\hat{pld})^2
alpha_g=0.5# f=alpha*(pld-\hat{pld})^2
beta=2# capacity of DER=rated value*beta

# convert matpower case file version 2 to a pandapower net.
path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)
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

YLL=Ybus[np.ix_(busid_LL,busid_LL)]
#YLL=YLL.todense()
Z_pu=inv(YLL)
Zt=Z_pu/np.conj(v0)

YLS=Ybus[np.ix_(busid_LL,[busid_slack])]
#YY=np.multiply(Z_pu,YLS)
YY=Z_pu.dot(YLS)
vs=YY*v0
vs=vs.todense()
vs=np.squeeze(np.asarray(vs))

# My=[Zt -1i*Zt]
from scipy.sparse import hstack
My=hstack((Zt, -1j*Zt))
Ky=1/vm0*np.real(np.conj(v0)*My)
const_vm_FPL=np.full((nbus, 1), vm0)
R=np.real(1/vm0*np.conj(v0)*Zt)
Rt=csc_matrix.transpose(R)

X=np.real(1/vm0*np.conj(v0)*(-1j)*Zt)
Xt=csc_matrix.transpose(X)


# # available load
# pld=net.load.p_mw.to_numpy()
# qld=net.load.q_mvar.to_numpy()


# initialize controllable loads
epsi_pld=0.5
epsi_qld=0.5
epsi_pg=0
epsi_qg=0
epsi_l=0.5
v_l=0.95
v_u=1.05

# initial values of load 
pld_t=-net_t.load.p_mw.to_numpy()/sbase
pld_t=np.delete(pld_t, busid_s_ld)
qld_t=-net_t.load.q_mvar.to_numpy()/sbase
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
pll_ld_max[busid_ld_LL]=abs(pld)*beta
pll_ld_min=np.zeros(nLL)
pll_ld_min[busid_ld_LL]=-abs(pld)*beta

qll_ld_max=np.zeros(nLL)
qll_ld_max[busid_ld_LL]=abs(qld)*beta
qll_ld_min=np.zeros(nLL)
qll_ld_min[busid_ld_LL]=-abs(qld)*beta

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

pll_g=np.zeros(nLL)
pll_g[busid_g_LL]=pg# load is negative in linearized power flow model
qll_g=np.zeros(nLL)
qll_g[busid_g_LL]=qg

# capacity of Generator
pll_g_max=np.zeros(nLL)
pll_g_max[busid_g_LL]=abs(pg)
pll_g_min=np.zeros(nLL)
pll_g_min[busid_g_LL]=-abs(pg)

qll_g_max=np.zeros(nLL)
qll_g_max[busid_g_LL]=abs(qg)
qll_g_min=np.zeros(nLL)
qll_g_min[busid_g_LL]=-abs(qg)


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

psgen_net_t=pll_net_t-pll_ld_t-pll_g_t
qsgen_net_t=qll_net_t-qll_ld_t-qll_g_t

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
    dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
    
    # minimize deviation from (pll_ld,qll_ld)
    pll_ld_t=pll_ld_t-epsi_pld*(2*alpha_ld*(pll_ld_t-pll_ld)+dvcnstr_dp)
    qll_ld_t=qll_ld_t-epsi_qld*(2*alpha_ld*(qll_ld_t-qll_ld)+dvcnstr_dq)
    
    # minimize generation from coal generator
    pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)+dvcnstr_dp)
    qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)+dvcnstr_dq)
    
    lambda_u=lambda_u+epsi_l*(vmll-v_u)
    lambda_d=lambda_d+epsi_l*(v_l-vmll)
    
    # project
    pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
    pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
    qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
    qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
    
    pll_g_t=np.maximum(pll_g_t,pll_g_min)
    pll_g_t=np.minimum(pll_g_t,pll_g_max)
    qll_g_t=np.maximum(qll_g_t,qll_g_min)
    qll_g_t=np.minimum(qll_g_t,qll_g_max)
    
    lambda_u=np.maximum(lambda_u,0)
    lambda_d=np.maximum(lambda_d,0)
    
    # update net injections
    pll_net_t=psgen_net_t+pll_ld_t+pll_g_t
    qll_net_t=qsgen_net_t+qll_ld_t+qll_g_t 
    sll_net_t=pll_net_t+1j*qll_net_t
    

    # fixed point iteration methods as the power flow solver
    vll0=vll
    
    itr_pf=0
    itr_pf_max=20
    err_vll=1
    while (err_vll>1e-5 and itr_pf<itr_pf_max):
        vll=Z_pu*np.conj(sll_net_t/vll)+vs
        err_vll=max(abs(vll-vll0));
        vll0=vll
        itr_pf=itr_pf+1;
               
    if err_vll>1e-3:
        raise Exception('power flow diverge!\n')
        
    
    # print('predicted:')
    # print('max v:',abs(vll).max())
    # print('min v:',abs(vll).min())
         
    # dispatch load 
    pld_t=-pll_ld_t[busid_ld_LL]*sbase
    qld_t=-qll_ld_t[busid_ld_LL]*sbase
    
    net_t.load.p_mw[busid_ld_lds]=pld_t# load is positive in net.load
    net_t.load.q_mvar[busid_ld_lds]=qld_t
    
    # dispatch (p,v) for generator
    pg_t=pll_g_t[busid_g_LL]*sbase
    vg_t=abs(vll[busid_g_LL])
    
    net_t.gen.p_mw=pg_t
    net_t.gen.vm_pu=vg_t
    
    pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)

    # check whether generator tracks dispatch 
    # load
    pld_t_m=net_t.res_load.p_mw.to_numpy()
    pld_t_m=np.delete(pld_t_m, busid_s_ld)
    qld_t_m=net_t.res_load.q_mvar.to_numpy()
    qld_t_m=np.delete(qld_t_m, busid_s_ld)
    
    # generator
    pg_t_m=net_t.res_gen.p_mw.to_numpy()
    # qg_t_m=net_t.res_gen.q_mvar/sbase
    vg_t_m=net_t.res_gen.vm_pu.to_numpy()
    
    # mismatch between dispatch and measurements
    dpld=abs(pld_t_m-pld_t)
    dqld=abs(qld_t_m-qld_t)
    dpg=abs(pg_t_m-pg_t)
    dvg=abs(vg_t_m-vg_t)
    
    # print('maximum mismatch:')
    # print('pld:%.8f' %dpld.max())
    # print('qld:%.8f' %dqld.max())
    # print('pg:%.8f' %dpg.max())
    # print('v:%.8f' %dvg.max())
    
    vmll=net_t.res_bus.vm_pu.values[busid_LL]
    
    # print('after dispatch:')
    # print('max v:',vmll.max())
    # print('min v:',vmll.min())
    
    plt.figure(1)
    plt.plot(abs(vll))
    plt.plot(vmll)
    plt.show()
    dv=vmll-abs(vll)
    print('Mismatch:%.8f'%abs(dv).max())
    
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