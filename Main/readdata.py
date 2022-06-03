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

# reactive power outpout of generators and DPV
busid_Gen=net.gen.bus
Q_Gen=net.res_gen.q_mvar.to_numpy()

busid_Sgen=net.sgen.bus
Q_Sgen=net.res_sgen.q_mvar.to_numpy()

busid0_GSg=busid_Gen+busid_Sgen
id_GSg=np.argsort(busid_GSg)
Q0_GSg=np.concatenate((Q_Gen,Q_Sgen))

n_GSg=len(busid_Gen)+len(busid_Sgen)
Q_GSg=np.zeros(n_GSg,1)
for i in range(n_GSg):
    Q_GSg[i]=Q0_GSg[id_GSg[i]]

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
        self.Ld=0.95
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

# bound
v_u = float(input('What is the upper bound of voltage (p.u.)? '))
v_l = float(input('What is the lower bound of voltage (p.u.)? '))
vplt_min=v_l*0.9# bounds for plot
vplt_max=v_u*1.1

epsi_pg=0.05
epsi_qg=0.05
epsi_l=0.05
epsi_u=0.001#0.01 for np

# the following part is in time series loop
#id_t=50# time index
id_t=1# time index
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
# print('Baseline voltage:\n')
# print('vmax:\n',max(abs(vll)))
# print('vmax:\n',min(abs(vll)))

vm_net_t=net_nt.res_bus.vm_pu.to_numpy()
vmll_net_t=np.delete(vm_net_t, busid_slack)
dvmll=abs(abs(vll)-vmll_net_t)
print('dvmll:',max(dvmll))

va_net_t=net_nt.res_bus.va_degree.to_numpy()
vall_net_t=np.delete(va_net_t, busid_slack)
vall=np.angle(vll)*180/np.pi
dvall=abs(vall-vall_net_t)
print('dvall:',max(dvall))    

vmll=abs(vll)
vall=np.angle(vll)*180/np.pi

vmll0=vm[busid_LL]

iplt=0
print('max v:',vmll.max())
print('min v:',vmll.min())

# initial net injections of LL nodes: in panderpower, load is positive
# pll_net_t=np.real(S[busid_LL])# DON'T use net.res_bus.p_mw[busid_LL].to_numpy()/sbase, because capacitors are included in Ybus matrix
# qll_net_t=np.imag(S[busid_LL])
# sll_net_t=S[busid_LL]


lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)

# u_u=np.zeros(nbrh)
# u_l=np.zeros(nbrh)

vll0=vll


class result:
    def __init__(self, iter_max):
        self.vmax=np.full((iter_max,1),1.05)
        self.vmin=np.full((iter_max,1),0.95)
        self.lambda_u=np.zeros(iter_max)
        self.lambda_d=np.zeros(iter_max)
result1=result(iter_max)  

# complex voltage of all the buses including slack bus
v_t=np.ones(nbus,dtype=complex)

# Qfr=np.zeros(nbrh)
# Pfr=np.zeros(nbrh)

result1.u_u=np.zeros(iter_max)
result1.u_l=np.zeros(iter_max)
for iter in range(iter_max):
    # derivative of voltage constraints with respect to (p, q)
    dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
    dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
    dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
    dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
    
    dvcnstr_dgp=dvcnstr_dp[np.ix_(busid_g_LL)]
    dvcnstr_dgq=dvcnstr_dq[np.ix_(busid_g_LL)]
    
    
    # dvcnstr_dsp=dvcnstr_dp[np.ix_(busid_sg_LL)]
    # dvcnstr_dsq=dvcnstr_dq[np.ix_(busid_sg_LL)]
    
    # # derivative of apparent power flow with respect to (p_g,p_pv,p_battery,p_load)
    # Qfrdu=Qfr*(u_u-u_l)
    # Pfrdu=Pfr*(u_u-u_l)
    
    # dsf_dp=Dt.dot(Qfrdu)+Ct.dot(Pfrdu)
    # dsf_dp=np.squeeze(np.array(dsf_dp))
    # dsf_dq=Ct.dot(Qfrdu)-Dt.dot(Pfrdu)
    # dsf_dq=np.squeeze(np.array(dsf_dq))
    
    # dsf_dgp=dsf_dp[np.ix_(busid_g_LL)]
    # dsf_dgq=dsf_dq[np.ix_(busid_g_LL)]
    
    # dsf_dsp=dsf_dp[np.ix_(busid_sg_LL)]
    # dsf_dsq=dsf_dq[np.ix_(busid_sg_LL)]
    
    # # each bus has at most 1 load and generator, but may have multiple static generators
    # # minimize deviation from (pll_ld,qll_ld)
    # # pll_ld_t=pll_ld_t-epsi_pld*(2*alpha_ld*(pll_ld_t-pll_ld)+dvcnstr_dp+dsf_dp)
    # # qll_ld_t=qll_ld_t-epsi_qld*(2*alpha_ld*(qll_ld_t-qll_ld)+dvcnstr_dq+dsf_dq)
    
    # # minimize generation from coal generator
    # # pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)*w_pll_g_max+dvcnstr_dp+dsf_dp)
    # # qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)*w_pll_g_max+dvcnstr_dq+dsf_dq)
    # # pll_g_t=pll_g_t-epsi_pg*(2*alpha_g*(pll_g_t-pll_g)+dvcnstr_dp+dsf_dp)
    # # qll_g_t=qll_g_t-epsi_qg*(2*alpha_g*(qll_g_t-qll_g)+dvcnstr_dq+dsf_dq)
    
    # pg_t=pg_t-epsi_pg*(2*alpha_g*(pg_t-pg)+dvcnstr_dgp+dsf_dgp)
    # qg_t=qg_t-epsi_qg*(2*alpha_g*(qg_t-qg)+dvcnstr_dgq+dsf_dgq)
    
    pg_t=pg_t-epsi_pg*(2*alpha_g*(pg_t-pg)+dvcnstr_dgp)
    qg_t=qg_t-epsi_qg*(2*alpha_g*(qg_t-qg)+dvcnstr_dgq)
    
    # minimize sgeneration from coal sgenerator 
    # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)*w_psg_max+dvcnstr_dsp+dsf_dsp)
    # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)*w_psg_max+dvcnstr_dsq+dsf_dsq)
    # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp+dsf_dsp)
    # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq+dsf_dsq) 
    
    # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp)
    # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq) 
    
    
    # project1.095
    # pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
    # pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
    # qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
    # qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
    
    # pll_g_t=np.maximum(pll_g_t,pll_g_min)
    # pll_g_t=np.minimum(pll_g_t,pll_g_max)
    # qll_g_t=np.maximum(qll_g_t,qll_g_min)
    # qll_g_t=np.minimum(qll_g_t,qll_g_max)
    
    # pg_t=np.maximum(pg_t,pg_min)
    # pg_t=np.minimum(pg_t,pg_max)
    # qg_t=np.maximum(qg_t,pg_min)
    # qg_t=np.minimum(qg_t,pg_max)
    
    # psg_t=np.maximum(psg_t,psg_min)
    # psg_t=np.minimum(psg_t,psg_max)
    # qsg_t=np.maximum(qsg_t,qsg_min)# only project switched shunts
    # qsg_t=np.minimum(qsg_t,qsg_max)# only project switched shunts
    # qsg_t[-18:]=np.maximum(qsg_t[-18:],qsg_min[-18:])# only project switched shunts
    # qsg_t[-18:]=np.minimum(qsg_t[-18:],qsg_max[-18:])# only project switched shunts
    
    # gen to sll
    pll_g_t=np.matmul(gen_to_LL,pg_t)
    qll_g_t=np.matmul(gen_to_LL,qg_t)  
    
    # # sgen to sll
    # pll_sg_t=np.matmul(sgen_to_LL,psg_t)
    # qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
    
    # update net injections 
    pll_net_t=-pll_ld+pll_dpv+pll_g_t# +: in; -: out
    qll_net_t=-qll_ld+qll_dpv+qll_g_t
    sll_net_t=pll_net_t+1j*qll_net_t
    
    #pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
    

    # fixed point iteration methods as the power flow solver       
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
               
    # if err_vll>1e-3:
    #     raise Exception('power flow diverge!\n')
    
         
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
    # result1.u_u[iter]=u_u.max()
    # result1.u_l[iter]=u_l.max()
    result1.vmax[iter]=vmll.max()
    result1.vmin[iter]=vmll.min()
    
    # #update dual variable (apparent power flow)
    # # nonlinear complex power flow
    # v_t[busid_slack]=v0
    # v_t[busid_LL]=vll
    # vf_t=v_t[brh_fbus]
    # vt_t=v_t[brh_tbus]
    # If_t=(vf_t-vt_t)*brh_y
    # Sf_t=vf_t*np.conjugate(If_t)
    
    # # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    # # Pi_t=np.real(sll_net_t)
    # # Qi_t=np.imag(sll_net_t)
    # # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
    
    # # dpf=np.real(Sf_t)-Pf_e
    # # dqf=np.imag(Sf_t)-Qf_e
    
    # Sfm_t=abs(Sf_t)
    # Pfr=np.real(Sf_t)/Sfm_t
    # Qfr=np.imag(Sf_t)/Sfm_t
    
    # # inactive small branch flow
    # id_sf0=np.where(Sfm_t<1e-7)
    # Pfr[id_sf0]=0
    # Qfr[id_sf0]=0
    
    # # # dual variable
    # # u_u=u_u+epsi_u*(Sfm_t-Smax_brh)
    # # #u_l=u_l+epsi_u*(Smin_brh-Sfm_t)
    
    # # # project dual variables
    # # u_u=np.maximum(u_u,0)
    # # #u_l=np.maximum(u_l,0)   
    
    # # result1.u_u[iter]=u_u.max()
    # #result1.u_l[iter]=u_l.max()
iterations=list(range(iter_max))

print('Mismatch:',err_vll)
print('max v:',vmll.max())
print('min v:',vmll.min())

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


# gen optimal vs intial (p,v)
iplt+=1
plt.figure(iplt) 
plot_pg=plt.plot(pg,'.')
plot_pgt=plt.plot(pg_t,'.')
plt.title('Pgen (mw)')
#plt.legend((plot_psgt[0], plot_psg[0],plot_psgmin[0], plot_psgmax[0]), ('Optimal', 'Initial','Min', 'Max'))
plt.legend((plot_pgt[0], plot_pg[0]), ('Optimal', 'Initial'))
plt.grid(True)
plt.savefig(path_plt+'/Pgen.png', dpi=400) 

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
plt.savefig(path_plt+'/Qgen.png', dpi=400)  
