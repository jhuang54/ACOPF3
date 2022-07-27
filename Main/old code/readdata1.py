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
iter_max=1000
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
icase = 'Maui2022dm_rd_v33_shunt.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition
# icase= 'Maui2022dm_rd_AggregateGens.mat'# physical system simulator
# net_t=pc.from_mpc(path_output + '\\' + icase, f_hz=60)
icase = 'Maui2022dm_rd_v33_shunt_OnlyLoad.mat'
net_nt = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition

sbase=net.sn_mva#10 MVA

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

# # generator global bus id and local id in LL
# busid_gen=list(net.gen.bus)
# #busid_gen.remove(busid_slack)
# busid_g_LL=np.array(busid_LL).searchsorted(busid_gen)

# # load global bus id and local id in LL
# busid_ld=net.load.bus
# busid_ld=list(busid_ld)
# nlsbus=len(busid_ld)
# # map from load buse to load buses+slack bus
# busid_s_ld=np.array(busid_ld).searchsorted(busid_slack)
# busid_ld_lds=list(range(0,nlsbus))
# busid_ld_lds.remove(busid_s_ld)

# # map from load bus to index in LL
# busid_ld.remove(busid_slack)

# busid_ld_LL=np.array(busid_LL).searchsorted(busid_ld)

# generator information
# generators' bus id in code
path_maui=os.path.join(path_par, 'maui_dauc_rted')
geninfo_name='Gen_Info.csv'

df_geninfo = pd.read_csv (path_maui+'\\'+geninfo_name)

# convert generator bus name to busid in code
GName=df_geninfo.loc[:,'genname'].to_numpy()# gen name
GbusName=df_geninfo.loc[:,'busname'].to_numpy()-1# gen bus name in code
busid_gen=np.array(net.bus.name).searchsorted(GbusName)#busid in Maui code

GbusIDMapFGn=dict(zip(GName, busid_gen))
GbusIDMapFGbn=dict(zip(GbusName, busid_gen))#[bus name]


# generator real power limits
pmax=df_geninfo.loc[:,'Pmax'].to_numpy()/sbase# gen name
pmin=df_geninfo.loc[:,'Pmin'].to_numpy()/sbase# gen name

# dataframe: generator name, busid in Maui, Pmax, Pmin
Geninfo = {"GenName":GName,"bus name":GbusName, "busid": busid_gen, "Pmax": pmax, "Pmin":pmin}
Geninfo = pd.DataFrame(Geninfo)
Geninfo=Geninfo.set_index("GenName")


# import dispatch power of real-time economic dispatch
rted_name='RTEDOutput_Day0.xlsx'
df_rted = pd.read_excel (path_maui+'\\'+rted_name,sheet_name='Dispatch', index_col=[0])
GName_rt=df_rted.columns.values

# Dispatchable generator name; DPV name
GName_g=[]
clnid_g=[]
GName_DPV=[]
clnid_DPV=[]
for count, gname in enumerate(GName_rt):
    if gname[0:3]!='DPV':
        GName_g.append(gname)
        clnid_g.append(count)
    else:
        GName_DPV.append(gname)
        clnid_DPV.append(count)

# dispatched generator, DPV bus id in Maui code
busid_gen=[Geninfo.loc[gname,'busid'] for gname in GName_g]
busid_DPV=[Geninfo.loc[gname,'busid'] for gname in GName_DPV]

# generator p limits
pg_min=[Geninfo.loc[gname,'Pmin'] for gname in GName_g]
pg_max=[Geninfo.loc[gname,'Pmax'] for gname in GName_g]


# # dispatched generator bus id in Maui code exluding slack bus
busid_gen=list(busid_gen)
# slack bus relative id in busid_gen
SinGen=busid_gen.count(busid_slack)
if SinGen>0:
    busid_s_SGen=busid_gen.index(busid_slack)

    # remove slack bus in busid_gen
    busid_gen.remove(busid_slack)
    # remove capacity
    del pg_min[busid_s_SGen]
    del pg_max[busid_s_SGen]
# gen relative bus id in busid_LL
busid_g_LL=np.array(busid_LL).searchsorted(busid_gen)

# slack bus relative id in busid_DPV
busid_DPV=list(busid_DPV)
# slack bus relative id in busid_DPV
SinDPV=busid_DPV.count(busid_slack)
if SinDPV>0:
    # slack bus relative id in {slack, DPV} 
    busid_s_SDPV=busid_DPV.index(busid_slack)

    # remove slack bus in busid_DPV
    busid_DPV.remove(busid_slack)
    # gen relative bus id in busid_LL
busid_DPV_LL=np.array(busid_LL).searchsorted(busid_DPV)

# gen to bus in LL set
ngen=len(busid_gen)
gen_to_LL=np.zeros((nLL,ngen), dtype=int)
for igen in range(0,ngen):
    gen_to_LL[busid_g_LL[igen],igen]=1
    
# DPV to bus in LL set
nDPV=len(busid_DPV)

# DPV to bus in LL set
DPV_to_LL=np.zeros((nLL,nDPV),dtype=int)
for iDPV in range(0,nDPV):
    DPV_to_LL[busid_DPV_LL[iDPV],iDPV]=1
    
# import time series load
ts_fld='Time_Series'
rtLd_name='TIMESERIES_RT_LOAD.csv'
df_rtLd= pd.read_csv (path_maui+'\\'+ts_fld+'\\'+rtLd_name, index_col=[0])
df_rtLd=df_rtLd.set_index('time')
LdbusName=df_rtLd.columns.to_numpy()# load bus name in csv file
LdbusName=[int(LN)-1 for LN in LdbusName]#load bus name in code
busid_ld=np.array(net.bus.name).searchsorted(LdbusName)#load bus id in code

busid_ld=list(busid_ld)
SinLd=busid_ld.count(busid_slack)
if SinLd>0:
    # map from load buse to load buses+slack bus
    busid_s_SLd=busid_ld.index(busid_slack)

    # map from load bus to index in LL
    busid_ld.remove(busid_slack)

busid_ld_LL=np.array(busid_LL).searchsorted(busid_ld)
# power factor
class powerfactor:
    def __init__(self):
        self.Ld=0.9
        self.DPV=0.95
        self.Gen=0.95
pf=powerfactor()

# linear power flow model
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

#X=np.real(1/vm0*np.conj(v0)*(-1j)*Zt)
X=np.real(-1j*Kyr)
#Xt=csc_matrix.transpose(X)
Xt=np.transpose(X)
Xt=Xt

# the following part is in time series loop
id_t=50# time index
# import generator data
# initial values of dispatchable gen
pg_t=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
pg=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
# remove generation at slack bus
if SinGen>0:
    pg_t=np.delete(pg_t,busid_s_SGen)# variable
    pg=np.delete(pg,busid_s_SGen)# prefered value
qg_t=pg_t*np.tan(np.arccos(pf.Gen))
qg=pg*np.tan(np.arccos(pf.Gen)) 

pll_g_t=np.matmul(gen_to_LL,pg_t)
qll_g_t=np.matmul(gen_to_LL,qg_t)

# initial values of DPV
pdpv=df_rted.iloc[id_t,clnid_DPV].to_numpy()/sbase
if SinDPV>0:
    pdpv=np.delete(pdpv,busid_s_SDPV)
 
# sDPV to sLL 
qdpv=np.tan(np.arccos(pf.DPV))*pdpv
pll_dpv=np.matmul(DPV_to_LL,pdpv)
qll_dpv=np.matmul(DPV_to_LL,qdpv)

#time_index=df_rtLd.index[id_t]
pld=df_rtLd.iloc[id_t,:].to_numpy()/sbase
if SinLd>0:
    pld=np.delete(pld, busid_s_SLd)# delete load at slack bus
pll_ld=np.zeros(nLL)
pll_ld[busid_ld_LL]=pld# load is negative in linearized power flow model
qll_ld=np.zeros(nLL)
qll_ld[busid_ld_LL]=np.tan(np.arccos(pf.Ld))*pld# calculate Q based on power factor

# voltage of DCOPF
# net apparent power injection
# update net injections 
pll_net_t=-pll_ld+pll_dpv+pll_g_t# +: in; -: out 
qll_net_t=-qll_ld+qll_dpv+qll_g_t 
sll_net_t=pll_net_t+1j*qll_net_t

# fixed point iteration methods as the power flow solver
vll=np.ones((nLL,1))*v0
vll=vll[0]
vll0=vll

itr_pf=0
itr_pf_max=50
err_vll=1
while (err_vll>1e-5 and itr_pf<itr_pf_max):
    #vll=Z_pu*np.conj(sll_net_t/vll)+vs
    ILL=np.conj(sll_net_t/vll)
    dv=np.matmul(Z_pu,ILL.transpose())
    vll=np.squeeze(np.asarray(dv))+w
    err_vll=max(abs(vll-vll0))
    vll0=vll
    itr_pf=itr_pf+1

# dispatch load data
pp.runpp(net_nt, algorithm='nr', calculate_voltage_angles=True)
net_nt.load.p_mw=-np.real(sll_net_t[:-1])*sbase
net_nt.load.q_mvar=-np.imag(sll_net_t[:-1])*sbase
print('id_t:\n',id_t)
print('net load:\n',sum(sll_net_t))
pp.runpp(net_nt, algorithm='nr', calculate_voltage_angles=True)
