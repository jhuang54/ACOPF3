import pandapower as pp
#import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

import matplotlib.pyplot as plt

# convert matpower case file version 2 to a pandapower net.
path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)
net_t=pc.from_mpc(path_output + '\\' + icase, f_hz=60)
nbus = len(net.bus)
nLL=nbus-1

# run power flow
pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)

# Ybus
Ybus = net._ppc['internal']['Ybus']

#import math
#import cmath
import numpy as np

# power flow solution
vm = net.res_bus.vm_pu.values
va_dg = net.res_bus.va_degree.values
va = np.exp(1j * va_dg * np.pi / 180)
v = vm * va

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

busid_ld=net.load.bus
busid_ld=list(busid_ld)
nlsbus=len(busid_ld)
#busid_s_ld=np.array(busid_ld).searchsorted(busid_slack)
#busid_ld.remove(busid_s_ld)
# map from load buse to load buses+slack bus
busid_s_ld=np.array(busid_ld).searchsorted(busid_slack)
busid_ld_lds=list(range(0,nlsbus))
busid_ld_lds.remove(busid_s_ld)

# map from load bus to N
busid_ld.remove(busid_slack)
busid_ld_LL=np.array(busid_LL).searchsorted(busid_ld)

YLL=Ybus[np.ix_(busid_LL,busid_LL)]
Z_pu=inv(YLL)
vm0=net.ext_grid.vm_pu[0]
va0=np.exp(1j*net.ext_grid.va_degree[0]*np.pi/180)
v0=vm0*va0


Zt=Z_pu/np.conj(v0)
# My=[Zt -1i*Zt]
from scipy.sparse import hstack
My=hstack((Zt, -1j*Zt))
Ky=1/vm0*np.real(np.conj(v0)*My)
const_vm_FPL=np.full((nbus, 1), vm0)
R=np.real(1/vm0*np.conj(v0)*Zt)
Rt=csc_matrix.transpose(R)

X=np.real(1/vm0*np.conj(v0)*(-1j)*Zt)
Xt=csc_matrix.transpose(X)

#vm_FPL=const_vm_FPL+Ky*[pLL-pLL0;qLL-qLL0]#given operation point, like(0,0)

# # available load
# pld=net.load.p_mw.to_numpy()
# qld=net.load.q_mvar.to_numpy()


# initialize controllable loads
epsi_p=0.5
epsi_q=0.5
epsi_l=0.5
v_l=0.95
v_u=1.05

pld_t=net_t.load.p_mw.to_numpy()
pld_t=np.delete(pld_t, busid_s_ld)
qld_t=net_t.load.q_mvar.to_numpy()
qld_t=np.delete(qld_t, busid_s_ld)

pll_t=np.zeros(nLL)
pll_t[busid_ld_LL]=pld_t
qll_t=np.zeros(nLL)
qll_t[busid_ld_LL]=qld_t


pld=net.load.p_mw.to_numpy()
pld=np.delete(pld,busid_s_ld)
qld=net.load.q_mvar.to_numpy()
qld=np.delete(qld,busid_s_ld)

pll=np.zeros(nLL)
pll[busid_ld_LL]=pld
qll=np.zeros(nLL)
qll[busid_ld_LL]=qld

lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)

v=net.res_bus.vm_pu 
vll=v.drop(busid_slack)

print('max v:',vll.max())
print('min v:',vll.min())
# update controllable loads (pld_t,qld_t) and dual variables
iter_max=2000
alpha=1
beta=2
class result:
    def __init__(self, iter_max):
        self.vmax=np.full((iter_max,1),1.05)
        self.vmin=np.full((iter_max,1),0.95)
        self.lambda_u=np.full((iter_max,1),0)
        self.lambda_d=np.full((iter_max,1),0)
result1=result(iter_max)  
for iter in range(iter_max):
    # dpll=pll_t-pll
    # dqll=qll_t-qll
    # print('max dpll: \n',dpll[busid_ld_LL].max())
    # print('min dpll: \n',dpll[busid_ld_LL].min())
    # print('max dqll: \n',dqll[busid_ld_LL].max())
    # print('min dqll: \n',dqll[busid_ld_LL].min())
    pll_t=pll_t-epsi_p*(2*alpha*(pll_t-pll)-Rt.dot(lambda_u-lambda_d))
    qll_t=qll_t-epsi_q*(2*alpha*(qll_t-qll)-Xt.dot(lambda_u-lambda_d))
    
    lambda_u=lambda_u+epsi_l*(vll-v_u)
    lambda_d=lambda_d+epsi_l*(v_l-vll)
    
    # project
    pll_t=np.maximum(pll_t,-abs(pll)*beta)
    pll_t=np.minimum(pll_t,abs(pll)*beta)
    qll_t=np.maximum(qll_t,-abs(qll)*beta)
    qll_t=np.minimum(qll_t,abs(qll)*beta)
    
    lambda_u=np.maximum(lambda_u,0)
    lambda_d=np.maximum(lambda_d,0)
    
    # update load for power flow calculation
    pld_t=pll_t[busid_ld_LL]
    qld_t=qll_t[busid_ld_LL]
    
    net_t.load.p_mw[busid_ld_lds]=pld_t
    net_t.load.q_mvar[busid_ld_lds]=qld_t
    pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)
    v=net_t.res_bus.vm_pu.values
    vll=v[busid_LL]
    
    # track voltage and dual
    result1.lambda_u[iter]=lambda_u.max()
    result1.lambda_d[iter]=lambda_d.max()
    result1.vmax[iter]=vll.max()
    result1.vmin[iter]=vll.min()
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