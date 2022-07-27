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
iter_max=1500
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
#icase = 'Maui2022dm_rd_v33_shunt.mat'
icase = 'Maui2022dm_rd_v33_shunt_QGen.mat'
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

# reactive power outpout and maximum value of generators and DPV
busid_Gen=net.gen.bus.to_numpy()
Q_Gen=net.res_gen.q_mvar.to_numpy()/sbase
Qmax_Gen=net.gen.max_q_mvar.to_numpy()/sbase

busid_Sgen=net.sgen.bus.to_numpy()
Q_Sgen=net.res_sgen.q_mvar.to_numpy()/sbase
Qmax_Sgen=net.sgen.max_q_mvar.to_numpy()/sbase

busid0_GSg=np.concatenate((busid_Gen,busid_Sgen))
busid_GSG=np.sort(busid0_GSg)

id_GSg=np.argsort(busid0_GSg)

Q0_GSg=np.concatenate((Q_Gen,Q_Sgen))
Qmax0_GSg=np.concatenate((Qmax_Gen,Qmax_Sgen))
n_GSg=len(busid_Gen)+len(busid_Sgen)
Q_GSg=np.zeros(n_GSg)
Qmax_GSg=np.zeros(n_GSg)
for i in range(n_GSg):
    Q_GSg[i]=Q0_GSg[id_GSg[i]]
    Qmax_GSg[i]=Qmax0_GSg[id_GSg[i]]
#Qmax_GSg=np.maximum(abs(Qmax_GSg),abs(Q_GSg)*1.05)
QGen_tab = {"busid": busid_GSG,"q_mvar":Q_GSg,'qmax_mvar':Qmax_GSg}
QGen_tab = pd.DataFrame(QGen_tab)
QGen_tab=QGen_tab.set_index("busid")

busid_ld=list(net.load.bus.to_numpy())
QLd_tab={"busid":busid_ld,"q_mvar":net.load.q_mvar.to_numpy()/sbase}
QLd_tab = pd.DataFrame(QLd_tab)
QLd_tab=QLd_tab.set_index("busid")

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
busid_genrt=[Geninfo.loc[gname,'busid'] for gname in GName_g]
busid_DPVrt=[Geninfo.loc[gname,'busid'] for gname in GName_DPV]

# generator p limits
pg_min=[Geninfo.loc[gname,'Pmin'] for gname in GName_g]
pg_max=[Geninfo.loc[gname,'Pmax'] for gname in GName_g]


# # dispatched generator bus id in Maui code exluding slack bus
busid_genrt=list(busid_genrt)
# slack bus relative id in busid_gen
SinGen=busid_genrt.count(busid_slack)
if SinGen>0:
    busid_s_SGen=busid_genrt.index(busid_slack)
    GenName_s=GName_g[busid_s_SGen]

    # remove slack bus in busid_genrt
    busid_genrt.remove(busid_slack)
    GName_g.remove(GenName_s)
    # remove capacity
    del pg_min[busid_s_SGen]
    del pg_max[busid_s_SGen]
# gen relative bus id in busid_LL
busid_g_LL=np.array(busid_LL).searchsorted(busid_genrt)

# slack bus relative id in busid_DPVrt
busid_DPVrt=list(busid_DPVrt)
# slack bus relative id in busid_DPVrt
SinDPV=busid_DPVrt.count(busid_slack)
if SinDPV>0:
    # slack bus relative id in {slack, DPV} 
    busid_s_SDPV=busid_DPVrt.index(busid_slack)

    # remove slack bus in busid_DPVrt
    busid_DPVrt.remove(busid_slack)
    # gen relative bus id in busid_LL
busid_DPVrt_LL=np.array(busid_LL).searchsorted(busid_DPVrt)

# shunt data
ShuntName=np.array([9,13,42,125,129,134,135,150,216,225,235,405,415,723,803,817,834,917])-1
nshunt=len(ShuntName)
busid_shunt=np.array(net.bus.name).searchsorted(ShuntName)#load bus id in code
qshunt_min=np.zeros(nshunt)/sbase
qshunt_max=np.array([0.600,3.600,0.600,3.600,3.600,3.600,3.600,3.600,1.800,3.600,3.600,3.600,3.600,3.600,3.600,1.200,3.600,3.600])/sbase

# slack bus relative id in busid_shunt
busid_shunt=list(busid_shunt)
SinShunt=busid_shunt.count(busid_slack)
if SinShunt>0:
    busid_s_Shunt=busid_shunt.index(busid_slack)

    # remove slack bus in busid_shunt
    busid_shunt.remove(busid_slack)

    # remove capacity
    del qshunt_min[busid_s_Shunt]
    del qshunt_max[busid_s_Shunt]
# gen relative bus id in busid_LL
busid_shunt_LL=np.array(busid_LL).searchsorted(busid_shunt)

# gen to bus in LL set
ngen=len(busid_genrt)
gen_to_LL=np.zeros((nLL,ngen), dtype=int)
for igen in range(0,ngen):
    gen_to_LL[busid_g_LL[igen],igen]=1
    
# DPV to bus in LL set
nDPV=len(busid_DPVrt)

# DPV to bus in LL set
DPV_to_LL=np.zeros((nLL,nDPV),dtype=int)
for iDPV in range(0,nDPV):
    DPV_to_LL[busid_DPVrt_LL[iDPV],iDPV]=1
    
# shunt to bus in LL set
shunt_to_LL=np.zeros((nLL,nshunt),dtype=int)
for ishunt in range(0,nshunt):
    shunt_to_LL[busid_shunt_LL[ishunt],ishunt]=1    
        
# import time series battery
df_rted_ESSDis = pd.read_excel (path_maui+'\\'+rted_name,sheet_name='ESS Dis', index_col=[0])
df_rted_ESSCh = pd.read_excel (path_maui+'\\'+rted_name,sheet_name='ESS Ch', index_col=[0])
nbat=len(df_rted_ESSDis.columns)
BatName_rt=list(np.array([12033,1203,93101,93102])-1)
busid_Batrt=np.array(net.bus.name).searchsorted(BatName_rt)
busid_Batrt=list(busid_Batrt)    

SinBat=busid_Batrt.count(busid_slack)
if SinBat>0:
    busid_s_Batrt=busid_Batrt.index(busid_slack)

    # remove slack bus in busid_shunt
    busid_s_Batrt.remove(busid_slack)

    # # remove capacity
    # del pbat_min[busid_s_Batrt]
    # del pbat_max[busid_s_Batrt]
# gen relative bus id in busid_LL
busid_bat_LL=np.array(busid_LL).searchsorted(busid_Batrt)

# battery to bus in LL set
bat_to_LL=np.zeros((nLL,nbat),dtype=int)
for ibat in range(0,nbat):
    bat_to_LL[busid_bat_LL[ibat],ibat]=1 


# import time series load
ts_fld='Time_Series'
rtLd_name='TIMESERIES_RT_LOAD.csv'
df_rtLd= pd.read_csv (path_maui+'\\'+ts_fld+'\\'+rtLd_name, index_col=[0])
df_rtLd=df_rtLd.set_index('time')
LdbusName=df_rtLd.columns.to_numpy()# load bus name in csv file
LdbusName=[int(LN)-1 for LN in LdbusName]#load bus name in code
busid_ldrt=np.array(net.bus.name).searchsorted(LdbusName)#load bus id in code

busid_ldrt=list(busid_ldrt)
SinLd=busid_ldrt.count(busid_slack)
if SinLd>0:
    # map from load buse to load buses+slack bus
    busid_s_SLdrt=busid_ldrt.index(busid_slack)

    # map from load bus to index in LL
    busid_ldrt.remove(busid_slack)

busid_ldrt_LL=np.array(busid_LL).searchsorted(busid_ldrt)





##power factor
class powerfactor:
    def __init__(self):
        self.Ld=0.95
        self.DPV=1
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

#epsi_u=0.001#0.01 for np
# stepsize
ACorDCOPF=input('Do you run E-DCOPF or ACOPF? Type 0 for E-DC and 1 for AC.')
if ACorDCOPF=='0':
    epsi_pg=0
    epsi_qg=0
    epsi_qshunt=0.05 
    OPF_type='EDCOPF'
elif ACorDCOPF=='1':
    epsi_qg=0.05
    OPF_type='ACOPF'
    
    # control pgen
    ControlPgen = input('Does ACOPF control real power of generators? Type 0 for No and 1 for yes')
    if ControlPgen=='0':
       epsi_pg=0 
    elif ControlPgen=='1':
       epsi_pg=0.05 
       
    ControlSwitchS = input('Does ACOPF control switched shunts? Type 0 for No and 1 for yes')
    if ControlSwitchS=='0':
       epsi_qshunt=0 
    elif ControlSwitchS=='1':
       epsi_qshunt=0.05 
    
epsi_l=0.35# dual variables associating with voltage constraints

# the following part is in time series loop
#id_t=50# time index
n_snaps=int(24*60/5)
        
class result_interm:
    def __init__(self, nLL, iter_max):
        self.vmax=np.full((iter_max,1),1.05)
        self.vmin=np.full((iter_max,1),0.95)
        self.lambda_u=np.zeros((nLL,iter_max))
        self.lambda_d=np.zeros((nLL,iter_max))
        self.v=np.zeros((nLL,iter_max))
        self.qshunt_t=np.zeros((nshunt,iter_max))      
        
class result:
    def __init__(self,n_snaps,ngen,nshunt):
        self.vmax0=np.zeros(n_snaps)
        self.vmin0=np.zeros(n_snaps)
        self.vmax1=np.zeros(n_snaps)
        self.vmin1=np.zeros(n_snaps)
        self.v0=np.zeros((nLL,n_snaps))# before integer projection
        self.v1=np.zeros((nLL,n_snaps))# after integer projection
        
        self.dTPg=np.zeros(n_snaps)
        self.rdTPg=np.zeros(n_snaps)
        self.dPg=np.zeros((ngen,n_snaps))
        self.rdPg=np.zeros((ngen,n_snaps))
        
        self.dTQg=np.zeros(n_snaps)
        self.rdTQg=np.zeros(n_snaps)
        self.dQg=np.zeros((ngen,n_snaps))
        self.rdQg=np.zeros((ngen,n_snaps))
        
        self.Qs0=np.zeros((nshunt,n_snaps)) 
        self.Qs1=np.zeros((nshunt,n_snaps)) 

        self.TQsopt0=np.zeros(n_snaps)  
        self.TQsopt1=np.zeros(n_snaps)
        self.TQsini1=np.zeros(n_snaps)      
        
        self.dTQs0=np.zeros(n_snaps)
        self.rdTQs0=np.zeros(n_snaps)
        self.dQs0=np.zeros((nshunt,n_snaps))
        self.rdQs0=np.zeros((nshunt,n_snaps))
        self.dQs0=np.zeros((nshunt,n_snaps))
        self.rdQs0=np.zeros((nshunt,n_snaps))
        
        self.dTQs1=np.zeros(n_snaps)
        self.rdTQs1=np.zeros(n_snaps)
        self.dQs1=np.zeros((nshunt,n_snaps))
        self.rdQs1=np.zeros((nshunt,n_snaps))
        self.dQs1=np.zeros((nshunt,n_snaps))
        self.rdQs1=np.zeros((nshunt,n_snaps))        
        
        self.dTsteps_Qs=np.zeros(n_snaps)
        self.adTsteps_Qs=np.zeros(n_snaps)        
        self.dsteps_Qs=np.zeros((nshunt,n_snaps))
        self.adsteps_Qs=np.zeros((nshunt,n_snaps+1))
        
        self.err_vll=np.zeros(n_snaps)
        self.gap_p=np.zeros(n_snaps)
        self.gap_q=np.zeros(n_snaps)
        
        self.rgap_p=np.zeros(n_snaps)
        self.rgap_q=np.zeros(n_snaps)

class Inputdata:
    def __init__(self,n_snaps):
        self.pld_sum=np.zeros(n_snaps)  
        self.qld_sum=np.zeros(n_snaps)
        self.pf=np.zeros(n_snaps)

result_int=result_interm(nLL,iter_max)       
        
#result_interm1=result_interm(iter_max)  
result0=result(n_snaps,ngen,nshunt)
result1=result(n_snaps,ngen,nshunt)
Input0=Inputdata(n_snaps)

vll_b=np.ones((nLL,1))*v0
vll_b=vll_b[0]
vll0_b=vll_b

vll=np.ones((nLL,1))*v0
vll=vll[0]
vll0=vll

qshunt_t=np.zeros(nshunt)
qshunt_t[0]=0.6/sbase
qshunt=np.zeros(nshunt)
qshunt[0]=0.6/sbase
if SinShunt>0:
    qshunt_t=np.delete(qshunt_t,busid_s_Shunt)
    qshunt=np.delete(qshunt,busid_s_Shunt)
qll_shunt_t0=np.matmul(shunt_to_LL,qshunt_t)# for DCOPF
qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)# for EDCOPF
Bi=0.06# admitance increment of each of Ni steps in block i

id_observe=228
for id_t in range(0,n_snaps):
    # import generator data
    # initial values of dispatchable gen
    pg_t=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
    pg=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
    # remove generation at slack bus
    if SinGen>0:
        pg_t=np.delete(pg_t,busid_s_SGen)# variable
        pg=np.delete(pg,busid_s_SGen)# prefered value
    qg_t=np.zeros(ngen)
    qg=np.zeros(ngen)
    qgmax=np.zeros(ngen)
    idG_bus=np.zeros(nbus,dtype=int)
    for i in range(ngen):
        idbus_tp=busid_genrt[i]
        if isinstance(QGen_tab.loc[idbus_tp,'q_mvar'], pd.Series):
            qg_t[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
            qg[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
            qgmax[i]=QGen_tab.loc[idbus_tp,'qmax_mvar'].to_numpy()[idG_bus[idbus_tp]]
            idG_bus[idbus_tp]=idG_bus[idbus_tp]+1
        else:
            qg_t[i]=QGen_tab.loc[idbus_tp,'q_mvar']
            qg[i]=QGen_tab.loc[idbus_tp,'q_mvar']
            qgmax[i]=QGen_tab.loc[idbus_tp,'qmax_mvar']

    pll_g_t=np.matmul(gen_to_LL,pg_t)
    qll_g_t=np.matmul(gen_to_LL,qg_t)

    # initial values of DPV
    pdpv=df_rted.iloc[id_t,clnid_DPV].to_numpy()/sbase
    if SinDPV>0:
        pdpv=np.delete(pdpv,busid_s_SDPV)
     
    qdpv=np.zeros(nDPV)   
    idG_bus=np.zeros(nbus,dtype=int)
    for i in range(nDPV):
        idbus_tp=busid_DPVrt[i]
        if isinstance(QGen_tab.loc[idbus_tp,'q_mvar'], pd.Series):
            qdpv[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
            idG_bus[idbus_tp]=idG_bus[idbus_tp]+1
        else:
            qdpv[i]=QGen_tab.loc[idbus_tp,'q_mvar']
        
    # sDPV to sLL 
    # # test
    # if id_t>=212:
    #     #qdpv=np.tan(np.arccos(pf.DPV))*pdpv
    #     qdpv=qdpv*0.2
    pll_dpv=np.matmul(DPV_to_LL,pdpv)
    qll_dpv=np.matmul(DPV_to_LL,qdpv)
    
    # battery
    pbat_t=(df_rted_ESSDis.iloc[id_t,:].to_numpy()-df_rted_ESSCh.iloc[id_t,:].to_numpy())/sbase
    pbat=(df_rted_ESSDis.iloc[id_t,:].to_numpy()-df_rted_ESSCh.iloc[id_t,:].to_numpy())/sbase

    if SinBat>0:
        pbat_t=np.delete(pbat_t,busid_s_Batrt)
        pbat=np.delete(pbat,busid_s_Batrt)
    pll_bat_t=np.matmul(bat_to_LL,pbat_t)

    #time_index=df_rtLd.index[id_t]
    pld=df_rtLd.iloc[id_t,:].to_numpy()/sbase
    if SinLd>0:
        pld=np.delete(pld, busid_s_SLdrt)# delete load at slack bus
    pll_ld=np.zeros(nLL)
    pll_ld[busid_ldrt_LL]=pld# load is negative in linearized power flow model
    qll_ld=np.zeros(nLL)
    #qll_ld[busid_ld_LL]=np.tan(np.arccos(pf.Ld))*pld# calculate Q based on power factor
    for i in range(len(busid_ldrt)):
        busid_tp=busid_ldrt[i]
        LdrtinLd=busid_ld.count(busid_tp)
        if LdrtinLd>0:# if real-time load is in pandapower load data
           qll_ld[busid_ldrt_LL[i]]=QLd_tab.loc[busid_tp,'q_mvar']# calculate Q based on power factor

    # voltage of DCOPF
    # net apparent power injection
    # update net injections 
    pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out 
    qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t0 
    sll_net_t=pll_net_t+1j*qll_net_t

    # fixed point iteration methods as the power flow solver
    itr_pf=0
    itr_pf_max=50
    err_vll=1
    while (err_vll>1e-5 and itr_pf<itr_pf_max):
        #vll=Z_pu*np.conj(sll_net_t/vll)+vs
        ILL=np.conj(sll_net_t/vll_b)
        dv=np.matmul(Z_pu,ILL.transpose())
        vll_b=np.squeeze(np.asarray(dv))+w
        err_vll=max(abs(vll_b-vll0_b))
        vll0_b=vll_b
        itr_pf=itr_pf+1

    # # dispatch load data
    # pp.runpp(net_nt, algorithm='nr', calculate_voltage_angles=True)
    # net_nt.load.p_mw=-np.real(sll_net_t[:-1])*sbase
    # net_nt.load.q_mvar=-np.imag(sll_net_t[:-1])*sbase
    # print('id_t:\n',id_t)
    # print('net load:\n',sum(sll_net_t))
    # pp.runpp(net_nt, algorithm='nr', calculate_voltage_angles=True)
    # print('Baseline voltage:\n')
    # print('vmax:\n',max(abs(vll)))
    # print('vmin:\n',min(abs(vll)))
    result0.err_vll[id_t]=err_vll
    result0.gap_p[id_t]=sum(pll_net_t)
    result0.gap_q[id_t]=sum(qll_net_t)
    result0.rgap_p[id_t]=sum(pll_net_t)/sum(pll_ld)
    result0.rgap_q[id_t]=sum(qll_net_t)/sum(qll_ld)
    result0.vmin1[id_t]=min(abs(vll_b))
    result0.vmax1[id_t]=max(abs(vll_b))
    
    Input0.pld_sum[id_t]=sum(pll_ld) 
    Input0.qld_sum[id_t]=sum(qll_ld) 
    Input0.pf[id_t]=abs(sum(pll_ld))/np.sqrt(sum(pll_ld)**2+sum(qll_ld)**2)
    
    # vm_net_t=net_nt.res_bus.vm_pu.to_numpy()
    # vmll_net_t=np.delete(vm_net_t, busid_slack)
    # dvmll=abs(abs(vll)-vmll_net_t)
    # print('dvmll:',max(dvmll))

    # va_net_t=net_nt.res_bus.va_degree.to_numpy()
    # vall_net_t=np.delete(va_net_t, busid_slack)
    # vall=np.angle(vll)*180/np.pi
    # dvall=abs(vall-vall_net_t)
    # print('dvall:',max(dvall))    

    # vmll=abs(vll)
    # vall=np.angle(vll)*180/np.pi

    # vmll0=vm[busid_LL]

    # iplt=0
    # print('max v:',vmll.max())
    # print('min v:',vmll.min())

    # # initial net injections of LL nodes: in panderpower, load is positive
    # # pll_net_t=np.real(S[busid_LL])# DON'T use net.res_bus.p_mw[busid_LL].to_numpy()/sbase, because capacitors are included in Ybus matrix
    # # qll_net_t=np.imag(S[busid_LL])
    # # sll_net_t=S[busid_LL]


    lambda_u=np.zeros(nLL)
    lambda_d=np.zeros(nLL)

    # u_u=np.zeros(nbrh)
    # u_l=np.zeros(nbrh)

    #vll0=vll

    # # complex voltage of all the buses including slack bus
    # v_t=np.ones(nbus,dtype=complex)

    # Qfr=np.zeros(nbrh)
    # Pfr=np.zeros(nbrh)

    # result_interm1.u_u=np.zeros(iter_max)
    # result_interm1.u_l=np.zeros(iter_max)
    for iter in range(iter_max):
        # derivative of voltage constraints with respect to (p, q)
        dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
        dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
        dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
        dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
        
        # traditional generator
        dvcnstr_dgp=dvcnstr_dp[np.ix_(busid_g_LL)]
        dvcnstr_dgq=dvcnstr_dq[np.ix_(busid_g_LL)]
        
        # shunt
        dvcnstr_dshunt=dvcnstr_dq[np.ix_(busid_shunt_LL)]
        
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
        
        # generators
        # pg_t=pg_t-epsi_pg*(2*alpha_g*(pg_t-pg)+dvcnstr_dgp)
        # qg_t=qg_t-epsi_qg*(2*alpha_g*(qg_t-qg)+dvcnstr_dgq)
        
        # shunt
        qshunt_t=qshunt_t-epsi_qshunt*(2*alpha_g*(qshunt_t-qshunt)+dvcnstr_dshunt)
        
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
        
        qshunt_t=np.maximum(qshunt_t,qshunt_min)
        qshunt_t=np.minimum(qshunt_t,qshunt_max)
        
        result1.Qs0[:,id_t]=qshunt_t
        
        # gen to sll
        pll_g_t=np.matmul(gen_to_LL,pg_t)
        qll_g_t=np.matmul(gen_to_LL,qg_t)  
        
        # # sgen to sll
        # pll_sg_t=np.matmul(sgen_to_LL,psg_t)
        # qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
        
        # shunt to sll
        qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)  
        
        # update net injections 
        pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out
        qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t
        sll_net_t=pll_net_t+1j*qll_net_t
        
        #pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
        

        # fixed point iteration methods as the power flow solver       
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
        
        if id_t==id_observe:
            result_int.lambda_u[:,iter]=lambda_u
            result_int.lambda_d[:,iter]=lambda_d
            
            result_int.v[:,iter]=vmll
       
            result_int.qshunt_t[:,iter]=qshunt_t

        
        
    #     # track voltage and dual
    #     result_interm1.lambda_u[iter]=lambda_u.max()
    #     result_interm1.lambda_d[iter]=lambda_d.max()
    #     # result_interm1.u_u[iter]=u_u.max()
    #     # result_interm1.u_l[iter]=u_l.max()
    #     result_interm1.vmax[iter]=vmll.max()
    #     result_interm1.vmin[iter]=vmll.min()
        
    #     # #update dual variable (apparent power flow)
    #     # # nonlinear complex power flow
    #     # v_t[busid_slack]=v0
    #     # v_t[busid_LL]=vll
    #     # vf_t=v_t[brh_fbus]
    #     # vt_t=v_t[brh_tbus]
    #     # If_t=(vf_t-vt_t)*brh_y
    #     # Sf_t=vf_t*np.conjugate(If_t)
        
    #     # # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    #     # # Pi_t=np.real(sll_net_t)
    #     # # Qi_t=np.imag(sll_net_t)
    #     # # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    #     # # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
        
    #     # # dpf=np.real(Sf_t)-Pf_e
    #     # # dqf=np.imag(Sf_t)-Qf_e
        
    #     # Sfm_t=abs(Sf_t)
    #     # Pfr=np.real(Sf_t)/Sfm_t
    #     # Qfr=np.imag(Sf_t)/Sfm_t
        
    #     # # inactive small branch flow
    #     # id_sf0=np.where(Sfm_t<1e-7)
    #     # Pfr[id_sf0]=0
    #     # Qfr[id_sf0]=0
        
    #     # # # dual variable
    #     # # u_u=u_u+epsi_u*(Sfm_t-Smax_brh)
    #     # # #u_l=u_l+epsi_u*(Smin_brh-Sfm_t)
        
    #     # # # project dual variables
    #     # # u_u=np.maximum(u_u,0)
    #     # # #u_l=np.maximum(u_l,0)   
        
    #     # # result_interm1.u_u[iter]=u_u.max()
    #     # #result_interm1.u_l[iter]=u_l.max()
    # iterations=list(range(iter_max))

    # print('Mismatch:',err_vll)
    # print('max v:',vmll.max())
    # print('min v:',vmll.min())

    result1.v0[:,id_t]=vmll
    
    result1.vmax0[id_t]=max(vmll)
    result1.vmin0[id_t]=min(vmll)
    
    result1.dQs0[:,id_t]=qshunt_t-qshunt
    result1.rdQs0[:,id_t]=result1.dQs0[:,id_t]/qshunt_max

    result1.dTQs0[id_t]=sum(abs(result1.dQs0[:,id_t]))      
    result1.rdTQs0[id_t]=result1.dTQs0[id_t]/sum(qshunt_max)  

    result1.TQsopt0[id_t]=sum(qshunt_t)
    
    # project switch shunt into integer steps
    qshunt_t=np.round(qshunt_t/Bi)*Bi 
    
    qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)
    
    # update net injections 
    pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out
    qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t
    sll_net_t=pll_net_t+1j*qll_net_t

    # fixed point iteration methods as the power flow solver       
    # vll0=vll_b
    # vll=vll_b
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
        
    vmll=abs(vll)
    result1.vmax1[id_t]=vmll.max()
    result1.vmin1[id_t]=vmll.min()
    result1.v1[:,id_t]=vmll
       
    # pg       
    result1.dPg[:,id_t]=pg_t-pg
    pg_a=pg
    pg_a[np.where(abs(pg)<1e-4)]=1    
    result1.rdPg[:,id_t]=result1.dPg[:,id_t]/pg_a
    
    result1.dTPg[id_t]=sum(abs(result1.dPg[:,id_t]))
    if sum(abs(pg))>0:
        result1.rdTPg[id_t]=result1.dTPg[id_t]/sum(abs(pg))
        
    # qg    
    result1.dQg[:,id_t]=qg_t-qg
    qg_a=qg
    qg_a[np.where(abs(qg)<1e-4)]=1  
    result1.rdQg[:,id_t]=result1.dQg[:,id_t]/pg_a
    
    result1.dTQg[id_t]=sum(abs(result1.dQg[:,id_t]))
    if sum(abs(qg))>0:
        result1.rdTQg[id_t]=result1.dTQg[id_t]/sum(abs(qgmax))
        
    # qshunt
    result1.Qs1[:,id_t]=qshunt_t
    
    result1.dQs1[:,id_t]=qshunt_t-qshunt
    result1.rdQs1[:,id_t]=result1.dQs1[:,id_t]/qshunt_max

    result1.dTQs1[id_t]=sum(abs(result1.dQs1[:,id_t]))      
    result1.rdTQs1[id_t]=result1.dTQs1[id_t]/sum(qshunt_max)
    result1.TQsopt1[id_t]=sum(qshunt_t)
    result1.TQsini1[id_t]=sum(qshunt)
    
    
    
    
    result1.dsteps_Qs[:,id_t]=result1.dQs1[:,id_t]/Bi
    result1.adsteps_Qs[:,id_t]=result1.adsteps_Qs[:,id_t]+result1.dsteps_Qs[:,id_t]
    result1.adsteps_Qs[:,id_t+1]=result1.adsteps_Qs[:,id_t]
    
    result1.dTsteps_Qs[id_t]=result1.dTQs1[id_t]/Bi
    result1.adTsteps_Qs[id_t]=sum(result1.dTsteps_Qs)
    qshunt=qshunt_t

print('Max voltage:',max(result1.vmax1)) 
print('Min voltage:',max(result1.vmin1))     
iplt=1   
xlabels_time=['00:00','06:00','12:00','18:00','23:55']   
# plt.figure(iplt) 
# fig, ax1 = plt.subplots()
# color='tab:red'
# ax1.set_xlabel('Iteration index')
# ax1.set_ylabel('Load (Mw)',color=color)
# ax1.plot(range(1,n_snaps+1),Input0.pld_sum,'.',markersize=2.5,color=color)
# ax1.tick_params(axis='y', labelcolor='tab:red')
# #ax1.set_yticks([0,0.0005,0.001,0.0015])
# #plt.ylim(Padja_ylim)
# ax2=ax1.twinx()
# color='tab:green'
# ax2.set_ylabel('Power factor', color=color)  # we already handled the x-label with ax1
# ax2.plot(range(1,n_snaps+1),Input0.pf,'*',markersize=2.5, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# plt.xlim([1,n_snaps+1])
# #plt.ylim(Padjr_ylim)
# plt.title('Load vs Net injection')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/PldvsPFactor.png', dpi=400)   

    
    
# iplt+=1
# fig,ax = plt.subplots()
# plt.figure(iplt)    
# plot_vmax_DCOPF=ax.plot(range(1,n_snaps+1), result0.vmax,'*',markersize=3,color='red')
# plot_vmin_DCOPF=ax.plot(range(1,n_snaps+1), result0.vmin,'*',markersize=3,color='blue')
# plot_vmax_EDCOPF=ax.plot(range(1,n_snaps+1), result1.vmax,'.',markersize=3,color='green')
# plot_vmin_EDCOPF=ax.plot(range(1,n_snaps+1), result1.vmin,'.',markersize=3,color='orange')
# plot_upper=ax.plot(range(1,n_snaps+1), v_u*np.ones(n_snaps),linestyle='dashed',color='black')
# plot_lower=ax.plot(range(1,n_snaps+1), v_l*np.ones(n_snaps),linestyle='dashed',color='black')
# #ax.legend((plot_vmax_DCOPF[0], plot_vmin_DCOPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
# ax.legend((plot_vmax_DCOPF[0], plot_vmin_DCOPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0]), ('DCOPF Max V', 'DCOPF Min V',OPF_type+' Max V', OPF_type+' Min V'))
# ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
# ax.set_xticklabels(xlabels_time, rotation=0)
# ax.set_xlabel("Time")
# ax.set_ylabel('Voltage magnitude, p.u.')
# ax.set_ylim([0.7,1.15])
# ax.set_title('Voltage')
# ax1.set_xlabel('Time')
# ax.grid(True)
# fig.savefig(path_plt+'/'+OPF_type+'/v.png', dpi=400)    

# before project onto integer
# voltage profile
iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:blue'
ax.set_xlabel('Bus index')
ax.set_ylabel('Voltage magnitude, p.u.')
for i in range(nLL):
    plot_v0=ax.plot(range(1,n_snaps+1),result1.v0[i,:],'.',markersize=1.5,color=color)
plot_upper=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashed',color='red')
plot_lower=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashed',color='black')


ax.legend((plot_v0[0],plot_upper[0],plot_lower[0]), (OPF_type+' V', 'Upper bound','Lower bound'))
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel("Time")
ax.set_ylabel('Voltage magnitude, p.u.')
ax.set_xlim([0,n_snaps+1])
ax.set_ylim([0.7,1.15])
ax.set_title('Voltage Profile (continuous switch shunt control)')
ax.set_xlabel('Time')
ax.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/VProfile (continuous)', dpi=400)  
    
# optimal total qshunt
iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xlabel('Generator index')
ax1.set_ylabel('Total reactive power adjustment, MVar',color=color)
ax1.plot(range(1,n_snaps+1),result1.dTQs0*sbase,'.',markersize=2.5,color=color)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_xlabel("Time")
ax1.set_ylim([0,13])
#plt.ylim(Padja_ylim)
ax2=ax1.twinx()
color='tab:green'
ax2.set_ylabel('Percent of available capacity, %', color=color)  # we already handled the x-label with ax1
ax2.plot(range(1,n_snaps+1),result1.rdTQs0*100,'*',markersize=2.5, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0,20])
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('Switch Shunt Total Reactive Power Adjustment (Continuous Control)')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Qgen_Adjust(Continuous).png', dpi=400)  

# time table
time_table=[]#9 clock is 0900
for i_h in range(24):
    hr=str(int(i_h))
    if i_h<10:
        hr='0'+str(int(i_h))
    for i_m in range(0,56,5):
        ms=str(i_m)
        if i_m<10:
            ms='0'+str(i_m)
        time_table.append(hr+ms)
# # optimal individual qshunt
# for id_t in range(1):
    # iplt+=1
    # plt.figure(iplt) 
    # # figure out time
    # hr=int(id_t/12)
    # ms=id_t%12
    # fig, ax = plt.subplots()
    # color='tab:red'
    # ax.set_xlabel('Shunt index')
    # ax.set_ylabel('Reactive power injection, MVar')
    # plot_qshunt=ax.plot(range(1,nshunt+1),result1.Qs0[:,id_t]*sbase,'o', markerfacecolor='none',markersize=4,color='blue')
    # plot_dqshunt=ax.plot(range(1,nshunt+1),result1.dQs0[:,id_t]*sbase,'r.',markersize=5)
    # plot_max=ax.plot(range(1,nshunt+1), qshunt_max*sbase,linestyle='dashed',color='orange')
    # ax.tick_params(axis='y')
    # plt.xlim([0.8,nshunt+0.2])
    # plt.ylim([-0.1,max(qshunt_max)*sbase*1.5])
    # ax.legend((plot_qshunt[0], plot_dqshunt[0],plot_max[0]), ('Optimal reactive power injection', 'Reactive power injection adjustment','Shunt capacity'))
    # ax.set_xticks(list(range(1,nshunt+1,1)))
    # ax.set_xticklabels(list(range(1,nshunt+1)), rotation=0)
    # #plt.ylim(Padjr_ylim)
    # plt.title('Shunt Reactive Power Injection (Continuous Control)')
    # plt.grid(True)
    # plt.savefig(path_plt+'/'+OPF_type+'/Qshunt/continuous/Qshunt_'+time_table[id_t]+' (Continuous).png', dpi=400)  


# after project onto integer
iplt+=1
fig,ax = plt.subplots()
plt.figure(iplt)   
plot_vmax0=ax.plot(range(0,n_snaps), result1.vmax0,'.',markersize=3,color='green')
plot_vmin0=ax.plot(range(0,n_snaps), result1.vmin0,'.',markersize=3,color='orange') 
plot_vmax1=ax.plot(range(0,n_snaps), result1.vmax1,'.',markersize=3,color='blue')
plot_vmin1=ax.plot(range(0,n_snaps), result1.vmin1,'.',markersize=3,color='red')
plot_upper=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashed',color='black')
plot_lower=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashed',color='black')
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_DCOPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_vmax0[0], plot_vmin0[0],plot_vmax1[0], plot_vmin1[0]), (OPF_type+' Max V (continuous)', OPF_type+' Min V (continuous)',OPF_type+' Max V (discrete)', OPF_type+' Min V (discrete)'))
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel("Time")
ax.set_ylabel('Voltage magnitude, p.u.')
ax.set_xlim([0,n_snaps+1])
ax.set_ylim([0.7,1.15])
ax.set_title('Voltage Range')
ax.set_xlabel('Time')
ax.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/v.png', dpi=400)   


iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:blue'
ax.set_xlabel('Bus index')
ax.set_ylabel('Voltage magnitude, p.u.')
for i in range(nLL):
    plot_v0=ax.plot(range(1,n_snaps+1),result1.v1[i,:],'.',markersize=1.5,color=color)
plot_upper=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashed',color='red')
plot_lower=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashed',color='black')


ax.legend((plot_v0[0],plot_upper[0],plot_lower[0]), (OPF_type+' V', 'Upper bound','Lower bound'))
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel("Time")
ax.set_ylabel('Voltage magnitude, p.u.')
ax.set_xlim([0,n_snaps+1])
ax.set_ylim([0.7,1.15])
ax.set_title('Voltage Profile (Discrete switch shunt control)')
ax.set_xlabel('Time')
ax.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/VProfile (Discrete)', dpi=400)  

# Qshunt
iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xlabel('Generator index')
ax1.set_ylabel('Reactive power, MVar',color=color)
#plt0=ax1.plot(range(1,n_snaps+1),result1.TQsini1*sbase,'.',markersize=2.5,color='orange')
plt1=ax1.plot(range(1,n_snaps+1),result1.TQsopt0*sbase,'.',markersize=2.5,color='b')
plt2=ax1.plot(range(1,n_snaps+1),result1.TQsopt1*sbase,'.',markersize=1,color='r')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_xlabel("Time")
#ax1.set_ylim([min(result1.TQsopt0)*0.9,max(result1.TQsopt0)*1.1])
#plt.legend((plt0[0],plt1[0],plt2[0]),('Initial','Optimal (Continuous)','Optimal (Discrete)'))
plt.legend((plt1[0],plt2[0]),('Optimal (Continuous)','Optimal (Discrete)'))
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('Switch Shunt Total Reactive Power')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Qgen_total(Discrete).png', dpi=400)  


    
# optimal total qshunt
iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xlabel('Generator index')
ax1.set_ylabel('Total reactive power adjustment, MVar',color=color)
ax1.plot(range(1,n_snaps+1),result1.dTQs1*sbase,'.',markersize=2.5,color=color)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_xlabel("Time")
ax1.set_ylim([0,13])
#plt.ylim(Padja_ylim)
ax2=ax1.twinx()
color='tab:green'
ax2.set_ylabel('Percent of available capacity, %', color=color)  # we already handled the x-label with ax1
ax2.plot(range(1,n_snaps+1),result1.rdTQs1*100,'*',markersize=2.5, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0,20])
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('Switch Shunt Total Reactive Power Adjustment (Discrete Control)')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Qgen_Adjust(Discrete).png', dpi=400)  

# # optimal individual qshunt
# for id_t in range(1):
    # iplt+=1
    # plt.figure(iplt) 
    # # figure out time
    # hr=int(id_t/12)
    # ms=id_t%12
    # fig, ax = plt.subplots()
    # color='tab:red'
    # ax.set_xlabel('Shunt index')
    # ax.set_ylabel('Reactive power injection, MVar')
    # plot_qshunt=ax.plot(range(1,nshunt+1),result1.Qs1[:,id_t]*sbase,'o', markerfacecolor='none',markersize=4,color='blue')
    # plot_dqshunt=ax.plot(range(1,nshunt+1),result1.dQs1[:,id_t]*sbase,'r.',markersize=5)
    # plot_max=ax.plot(range(1,nshunt+1), qshunt_max*sbase,linestyle='dashed',color='orange')
    # ax.tick_params(axis='y')
    # plt.xlim([0.8,nshunt+0.2])
    # plt.ylim([-0.1,max(qshunt_max)*sbase*1.5])
    # ax.legend((plot_qshunt[0], plot_dqshunt[0],plot_max[0]), ('Optimal reactive power injection', 'Reactive power injection adjustment','Shunt capacity'))
    # ax.set_xticks(list(range(1,nshunt+1,1)))
    # ax.set_xticklabels(list(range(1,nshunt+1)), rotation=0)
    # #plt.ylim(Padjr_ylim)
    # plt.title('Shunt Reactive Power Injection Qs (Discrete Control)')
    # plt.grid(True)
    # plt.savefig(path_plt+'/'+OPF_type+'/Qshunt/discrete/Qshunt_'+time_table[id_t]+' (Discrete).png', dpi=400)  

# # continuous vs discrete
# for id_t in range(1):
    # iplt+=1
    # plt.figure(iplt) 
    # # figure out time
    # hr=int(id_t/12)
    # ms=id_t%12
    # fig, ax = plt.subplots()
    # color='tab:red'
    # ax.set_xlabel('Shunt index')
    # ax.set_ylabel('Reactive power injection, MVar')
    # plot_qshunt0=ax.plot(range(1,nshunt+1),result1.Qs0[:,id_t]*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
    # plot_qshunt1=ax.plot(range(1,nshunt+1),result1.Qs1[:,id_t]*sbase,'.',markersize=5,color='red')
    # plot_max=ax.plot(range(1,nshunt+1), qshunt_max*sbase,linestyle='dashed',color='orange')
    # ax.tick_params(axis='y')
    # plt.xlim([0.8,nshunt+0.2])
    # plt.ylim([-0.1,max(qshunt_max)*sbase*1.5])
    # ax.legend((plot_qshunt[0], plot_dqshunt[0],plot_max[0]), ('Continuous control', 'Discrete control','Shunt capacity'))
    # ax.set_xticks(list(range(1,nshunt+1,1)))
    # ax.set_xticklabels(list(range(1,nshunt+1)), rotation=0)
    # #plt.ylim(Padjr_ylim)
    # plt.title('Shunt Optimal Reactive Power Injection')
    # plt.grid(True)
    # plt.savefig(path_plt+'/'+OPF_type+'/Qshunt/continusousVSdiscrete/Qshunt_'+time_table[id_t]+'.png', dpi=400)  


# iplt+=1
# plt.figure(iplt) 
# fig, ax1 = plt.subplots()
# color='tab:red'
# ax1.set_xlabel('Generator index')
# ax1.set_ylabel('Mw',color=color)
# ax1.plot(range(1,n_snaps+1),result1.dTPg*sbase,'.',markersize=2.5,color=color)
# ax1.tick_params(axis='y', labelcolor='tab:red')
# ax1.set_yticks([0,0.0005,0.001,0.0015])
# #plt.ylim(Padja_ylim)
# ax2=ax1.twinx()
# color='tab:green'
# ax2.set_ylabel('Percent %', color=color)  # we already handled the x-label with ax1
# ax2.plot(range(1,n_snaps+1),result1.rdTPg*100,'*',markersize=2.5, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# plt.xlim([1,n_snaps+1])
# #plt.ylim(Padjr_ylim)
# plt.title('Real Power Adjustment')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/Pgen_Adjust.png', dpi=400) 

iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:red'

ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel('Time')
ax.set_ylabel('Step change')
ax.plot(range(0,n_snaps),result1.dTsteps_Qs,'.',markersize=4,color=color)
ax.tick_params(axis='y')
#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
plt.xlim([-3,n_snaps+1])
plt.ylim([0,round(1.1*max(result1.dTsteps_Qs))])
#plt.ylim(Padjr_ylim)
plt.title('Switch Shunt Step Changes')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Number_QshuntStepChange.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:red'

ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel('Time')
ax.set_ylabel('Step changes')
ax.plot(range(1,n_snaps+1),result1.adTsteps_Qs,'.',markersize=1.5,color=color)
ax.tick_params(axis='y', labelcolor='tab:red')
#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
plt.xlim([0,n_snaps+1])
plt.ylim([0,round(max(result1.adTsteps_Qs)*1.2)])
#plt.ylim(Padjr_ylim)
plt.title('Accumulated Switch Shunt Step Changes')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/AccumulatedNumber_QshuntStepChange.png', dpi=400)   

# switch shunt step change at one snapshot
id_snap=np.argmin(result1.vmin1)
iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:red'
ax.set_xticks(list(range(1,nshunt+1)))
ax.set_xticklabels(list(range(1,nshunt+1)), rotation=0)
ax.set_xlabel('Switch shunt index')
ax.set_ylabel('Step changes',color=color)
ax.plot(range(1,nshunt+1),result1.dsteps_Qs[:,id_snap],'.',markersize=2.5,color=color)

#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
plt.xlim([0.8,nshunt+0.2])
plt.ylim([-0.2,15])
#plt.ylim(Padjr_ylim)
plt.title('The Switch Shunt Step Changes')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Number_QshuntStepChange_vmin.png', dpi=400)   

####
# switch shunt step change of 288 snapshots
# time table
time_table=[]#9 clock is 0900
for i_h in range(24):
    hr=str(int(i_h))
    if i_h<10:
        hr='0'+str(int(i_h))
    for i_m in range(0,56,5):
        ms=str(i_m)
        if i_m<10:
            ms='0'+str(i_m)
        time_table.append(hr+ms)
        
for id_snap in range(0,n_snaps):
    iplt+=1
    plt.figure(iplt) 
    fig, ax = plt.subplots()
    color='tab:red'

    ax.set_xticks(list(range(1,nshunt+1)))
    ax.set_xticklabels(list(range(1,nshunt+1)), rotation=0)
    ax.set_xlabel('Switch shunt index')
    ax.set_ylabel('Step changes')
    ax.plot(range(1,nshunt+1),result1.dsteps_Qs[:,id_snap],'.',markersize=2.5,color=color)

    #ax1.set_yticks([0,0.0005,0.001,0.0015])
    #plt.ylim(Padja_ylim)
    plt.xlim([0.8,nshunt+0.2])
    plt.ylim([-0.2,15])
    #plt.ylim(Padjr_ylim)
    plt.title('The Switch Shunt Step Changes')
    plt.grid(True)
    plt.savefig(path_plt+'/'+OPF_type+'/Number_QshuntStepChange/'+time_table[id_snap]+'.png', dpi=400) 
####


# convergenc of primal and dual at id_t=228
id_vmin=np.argmin(result_int.v[:,-1])
iplt+=1
plt.figure(iplt) 
#plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
plot_ld=plt.plot(result_int.lambda_d[id_vmin,:],'.',markersize=3)
plt.title('Dual Variable Convergence')
#plt.legend((plot_lu[0], plot_ld[0]), ('$\overline{\u03BB}$', '$\underline{\u03BB}$'))
plt.xlim([0,iter_max])
plt.xlabel('Iteration number')
plt.ylabel('\u03BB',color='black')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Convergence/Convergence_lambda.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
#plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
plot_ld=plt.plot(result_int.v[id_vmin,:],'.',markersize=3)
plt.title('Voltage Convergence')
#plt.legend((plot_lu[0], plot_ld[0]), ('$\overline{\u03BB}$', '$\underline{\u03BB}$'))
plt.xlim([0,iter_max])
plt.xlabel('Iteration number')
plt.ylabel('Voltage magnitude, p.u.')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Convergence/Convergence_v.png', dpi=400)  

id_MaxQs=np.argmax(qshunt_max)
iplt+=1
plt.figure(iplt) 
#plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
plot_ld=plt.plot(result_int.qshunt_t[id_MaxQs,:]*sbase,'.',markersize=3)
plt.title('Reactive Power Injection Convergence')
plt.xlim([0,iter_max])
plt.xlabel('Iteration number')
plt.ylabel('Reactive power injection, MVar')
plt.grid(True)
plt.savefig(path_plt+'/'+OPF_type+'/Convergence/Convergence_qshunt.png', dpi=400)  




    
    

