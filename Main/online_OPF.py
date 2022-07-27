import pandapower as pp
#import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

import matplotlib.pyplot as plt

from scipy.sparse.csgraph import minimum_spanning_tree

import pandas as pd

import numpy as np

from func.util import ACOPF
from func.util import data_rand
from func.util import ACPF

path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_plt = os.path.join(path_cur, 'Plot')
path_input=os.path.join(path_cur, 'Input')

# parameter setting
iter_max=1500
alpha_ld=0.1# f=alpha*(pld-\hat{pld})^2
alpha_pg=0.5# f=alpha*(pg-\hat{pg})^2
alpha_qg=0.1# f=alpha*(pg-\hat{pg})^2
alpha_qs=0.1
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
pg_min=np.array(pg_min)
pg_max=np.array(pg_max)
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

# read OPF input: v_u,v_l,robustness,dc or ac(0 or 1),control pgen, control switch qshunt
path_opf=os.path.join(path_input,'OPF.txt')
with open(path_opf) as f:
     lines=f.readlines()
     for line in lines:
        if line[0]!="#":
           input_t=line.split(',')
           v_u=float(input_t[0])
           v_l=float(input_t[1])
           V_ROBUSTNESS=float(input_t[2])
           ControlQgen=input_t[3]
           ControlPgen=input_t[4]
           ControlSwitchS=input_t[5]
           Bi=float(input_t[6])/sbase#Bi=0.06# admitance increment of each of Ni steps in block i
f.close()
# bound
# v_u = float(input('What is the upper bound of voltage (p.u.)? '))
# v_l = float(input('What is the lower bound of voltage (p.u.)? '))
# V_ROBUSTNESS=float(input('What is the voltage robustness?'))
vplt_min=v_l*0.9# bounds for plot
vplt_max=v_u*1.1

#epsi_u=0.001#0.01 for np
# stepsize
#ACorDCOPF=input('Do you run E-DCOPF or ACOPF? Type 0 for E-DC and 1 for AC.')
if ControlQgen=='0':
   epsi_qg=0
else:
   epsi_qg=0.05
   
   
if ControlPgen=='0':
   epsi_pg=0
else:
   epsi_pg=0.05
   
if ControlSwitchS=='0':
   epsi_qshunt=0   
else:
   epsi_qshunt=0.05 
  
if ControlPgen=='0' and ControlQgen=='0': 
   OPF_type='EDCOPF'
else:
   if abs(V_ROBUSTNESS)<1e-4:
      OPF_type='ACOPF'
   else:      
      OPF_type='CCACOPF'
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
        self.pg_t=np.zeros((ngen,iter_max))
        self.qg_t=np.zeros((ngen,iter_max))
class result:
    def __init__(self,n_snaps,ngen,nshunt):
        self.vmax=np.zeros(n_snaps)
        self.vmin=np.zeros(n_snaps)
        self.TdPg=np.zeros(n_snaps)
        self.rTdPg=np.zeros(n_snaps)
        self.dPg=np.zeros((ngen,n_snaps))
        self.rdPg=np.zeros((ngen,n_snaps))
        
        self.TdQg=np.zeros(n_snaps)
        self.rTdQg=np.zeros(n_snaps)
        self.dQg=np.zeros((ngen,n_snaps))
        self.rdQg=np.zeros((ngen,n_snaps))
        
        self.TdQs=np.zeros(n_snaps)
        self.rTdQs=np.zeros(n_snaps)
        self.dQs=np.zeros((nshunt,n_snaps))
        self.rdQs=np.zeros((nshunt,n_snaps))
        
        self.Tdsteps_Qs=np.zeros(n_snaps)
        self.aTdsteps_Qs=np.zeros(n_snaps)        
        self.dsteps_Qs=np.zeros((nshunt,n_snaps))
        self.adsteps_Qs=np.zeros((nshunt,n_snaps+1))
        
        self.Pgt=np.zeros((ngen,n_snaps))
        self.Qgt=np.zeros((ngen,n_snaps))
        self.Qst=np.zeros((nshunt,n_snaps))
        self.vt=np.zeros((nLL,n_snaps))
        

        self.Pg=np.zeros((ngen,n_snaps))
        self.Qg=np.zeros((ngen,n_snaps))
        self.Qs=np.zeros((nshunt,n_snaps))
        self.v=np.zeros((nLL,n_snaps)) 
 
        
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
        
class alg_par:
    def __init__(self,iter_max,epsi_l,epsi_pg,epsi_qg,epsi_qshunt,itr_pf_max,tol_pf,nLL):
        self.iter_max=iter_max
        self.epsi_l=epsi_l
        self.epsi_pg=epsi_pg
        self.epsi_qg=epsi_qg
        self.epsi_qshunt=epsi_qshunt
        self.itr_pf_max=50
        self.tol_pf=1e-5
        
class mdl_par:
    def __init__(self,ngen,nshunt,nLL):
        self.ngen=ngen
        self.nshunt=nshunt
        self.nLL=nLL
        
        self.pg_min=np.zeros(ngen)
        self.pg_max=np.zeros(ngen)
        self.qg_min=np.zeros(ngen)
        self.qg_max=np.zeros(ngen) 
        self.qshunt_min=np.zeros(nshunt)
        self.qshunt_max=np.zeros(nshunt)  

        self.v_l=0
        self.v_u=0
        self.V_ROBUSTNESS=0


        self.Rt=np.zeros((nLL,nLL))
        self.Xt=np.zeros((nLL,nLL))
        self.Z_pu=np.zeros((nLL,nLL))
        self.w=np.zeros(nLL)
        
        self.busid_g_LL=[]
        self.busid_shunt_LL=[]
        
        self.gen_to_LL=np.zeros((nLL,ngen)) 
        self.shunt_to_LL=np.zeros((nLL,nshunt)) 
        
        self.alpha_pg=0
        self.alpha_qg=0
        self.alpha_qs=0
        
        self.pll_ld=np.zeros(nLL)
        self.pll_dpv=np.zeros(nLL)
        self.pll_bat_t=np.zeros(nLL)
        self.qll_ld=np.zeros(nLL)
        self.qll_dpv=np.zeros(nLL)
        
        self.vll=np.zeros(nLL,dtype=complex)
        self.vll0=np.zeros(nLL,dtype=complex)
        self.sll_net_t=np.zeros(nLL,dtype=complex)
        
class var:
    def __init__(self,ngen,nshunt,nLL):
        var.pg_t=np.zeros(ngen)
        var.qg_t=np.zeros(ngen)
        var.qshunt_t=np.zeros(nshunt)
        var.lambda_u=np.zeros(nLL)
        var.lambda_d=np.zeros(nLL)
        
    
# results        
result_int=result_interm(nLL,iter_max)  
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


# qgen
qg_t=np.zeros(ngen)
qg=np.zeros(ngen)
qg_max=np.zeros(ngen)
idG_bus=np.zeros(nbus,dtype=int)
for i in range(ngen):
    idbus_tp=busid_genrt[i]
    if isinstance(QGen_tab.loc[idbus_tp,'q_mvar'], pd.Series):
        qg_t[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
        qg[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
        qg_max[i]=QGen_tab.loc[idbus_tp,'qmax_mvar'].to_numpy()[idG_bus[idbus_tp]]
        idG_bus[idbus_tp]=idG_bus[idbus_tp]+1
    else:
        qg_t[i]=QGen_tab.loc[idbus_tp,'q_mvar']
        qg[i]=QGen_tab.loc[idbus_tp,'q_mvar']
        qg_max[i]=QGen_tab.loc[idbus_tp,'qmax_mvar']
        
qll_g_t=np.matmul(gen_to_LL,qg_t)   

qdpv=np.zeros(nDPV)   
idG_bus=np.zeros(nbus,dtype=int)
for i in range(nDPV):
    idbus_tp=busid_DPVrt[i]
    if isinstance(QGen_tab.loc[idbus_tp,'q_mvar'], pd.Series):
        qdpv[i]=QGen_tab.loc[idbus_tp,'q_mvar'].to_numpy()[idG_bus[idbus_tp]]
        idG_bus[idbus_tp]=idG_bus[idbus_tp]+1
    else:
        qdpv[i]=QGen_tab.loc[idbus_tp,'q_mvar']         

qll_dpv=np.matmul(DPV_to_LL,qdpv)

lambda_u=np.zeros(nLL)
lambda_d=np.zeros(nLL)
id_observe=156# where the minimum voltage appear


# model parameters
# algorithm parameters: iter_max,epsi_l,epsi_pg,epsi_qg,epsi_qshunt,itr_pf_max,tol_pf
itr_pf_max=50
tol_pf=1e-5
alg_par1=alg_par(iter_max,epsi_l,epsi_pg,epsi_qg,epsi_qshunt,itr_pf_max,tol_pf,nLL) 
   
mdl_par1=mdl_par(ngen,nshunt,nLL)

mdl_par1.pg_min=pg_min
mdl_par1.pg_max=pg_max
if ControlPgen==0:# Do Not control real power output of generators
    mdl_par1.pg_min=-100000*np.ones(ngen)
    mdl_par1.pg_max=100000*np.ones(ngen)

mdl_par1.qg_min=-qg_max
mdl_par1.qg_max=qg_max
if ControlQgen==0:
    mdl_par1.qg_min=-100000*np.ones(ngen)
    mdl_par1.qg_max=100000*np.ones(ngen)

mdl_par1.qshunt_min=qshunt_min
mdl_par1.qshunt_max=qshunt_max
if ControlSwitchS==0:
    mdl_par1.qshunt_min=-100000*np.ones(nshunt)
    mdl_par1.qshunt_max=100000*np.ones(nshunt)


mdl_par1.Rt=Rt
mdl_par1.Xt=Xt

mdl_par1.busid_g_LL=busid_g_LL
mdl_par1.busid_shunt_LL=busid_shunt_LL

mdl_par1.gen_to_LL=gen_to_LL
mdl_par1.shunt_to_LL=shunt_to_LL

mdl_par1.alpha_pg=alpha_pg
mdl_par1.alpha_qg=alpha_qg
mdl_par1.alpha_qs=alpha_qs

mdl_par1.qg=qg        

mdl_par1.v_l=v_l
mdl_par1.v_u=v_u
mdl_par1.V_ROBUSTNESS=V_ROBUSTNESS
mdl_par1.Z_pu=Z_pu
mdl_par1.w=w
var1=var(ngen,nshunt,nLL)  
     
# AC PF
alg_par_pf=alg_par(0,0,0,0,0,itr_pf_max,tol_pf,nLL)     
mdl_par_pf=mdl_par(ngen,nshunt,nLL) 
mdl_par_pf.Z_pu=Z_pu
mdl_par_pf.w=w

RECORD_V_VIOLATION_OPF=np.zeros(n_snaps)   
# read input of chance constraint acopf
#path_input=os.path.join(path_cur, 'Input')
path_input_ccacopf=os.path.join(path_input,'CCACOPF.txt')
with open(path_input_ccacopf,'r') as f:
    lines = f.readlines()
    for line in lines:
        if line[0]!='#':
           input_t=line.split(',')
           n_sim=int(input_t[0])
           mu=float(input_t[1])
           sigma=float(input_t[2])        
f.close()


for id_t in range(0,n_snaps):
#for id_t in range(0,1):
    # import generator data
    # initial values of dispatchable gen
    pg_t=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
    pg=df_rted.iloc[id_t,clnid_g].to_numpy()/sbase
    # remove generation at slack bus
    if SinGen>0:
        pg_t=np.delete(pg_t,busid_s_SGen)# variable
        pg=np.delete(pg,busid_s_SGen)# prefered value

    pll_g_t=np.matmul(gen_to_LL,pg_t)

    # initial values of DPV
    pdpv=df_rted.iloc[id_t,clnid_DPV].to_numpy()/sbase
    if SinDPV>0:
        pdpv=np.delete(pdpv,busid_s_SDPV)
        

    pll_dpv=np.matmul(DPV_to_LL,pdpv)

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

    # # fixed point iteration methods as the power flow solver
    # itr_pf=0
    # itr_pf_max=50
    # err_vll=1
    # while (err_vll>1e-5 and itr_pf<itr_pf_max):
    #     #vll=Z_pu*np.conj(sll_net_t/vll)+vs
    #     ILL=np.conj(sll_net_t/vll_b)
    #     dv=np.matmul(Z_pu,ILL.transpose())
    #     vll_b=np.squeeze(np.asarray(dv))+w
    #     err_vll=max(abs(vll_b-vll0_b))
    #     vll0_b=vll_b
    #     itr_pf=itr_pf+1

    # result0.err_vll[id_t]=err_vll
    # result0.gap_p[id_t]=sum(pll_net_t)
    # result0.gap_q[id_t]=sum(qll_net_t)
    # result0.rgap_p[id_t]=sum(pll_net_t)/sum(pll_ld)
    # result0.rgap_q[id_t]=sum(qll_net_t)/sum(qll_ld)
    # result0.vmin[id_t]=min(abs(vll_b))
    # result0.vmax[id_t]=max(abs(vll_b))  
    
    # result0.Pg[:,id_t]=pg
    
    # Input0.pld_sum[id_t]=sum(pll_ld) 
    # Input0.qld_sum[id_t]=sum(qll_ld) 
    # Input0.pf[id_t]=abs(sum(pll_ld))/np.sqrt(sum(pll_ld)**2+sum(qll_ld)**2)
   
    # initialize state variables
    var1.pg_t=pg_t
    var1.qg_t=qg_t
    var1.qshunt_t=qshunt_t
    var1.lambda_u=lambda_u
    var1.lambda_d=lambda_d    
    
    # load, DPV, battery data
    mdl_par1.pll_ld=pll_ld
    mdl_par1.pll_dpv=pll_dpv
    mdl_par1.pll_bat_t=pll_bat_t
    mdl_par1.qll_ld=qll_ld
    mdl_par1.qll_dpv=qll_dpv
    mdl_par1.pg=pg#DC OPF solution, updated in each snapshot
    mdl_par1.qshunt=qshunt# optimal qshunt at last snapshot,updated in each snapshot
    
    # parameter
    mdl_par1.vll=vll  
    mdl_par1.vll0=vll0  
    
    res2=ACOPF(alg_par1,mdl_par1,var1)
    
    # primal solution
    pg_t=res2.pg_t
    qg_t=res2.qg_t
    qshunt_t=res2.qshunt_t
    qshunt=qshunt_t 
      
    # dual solution
    lambda_u=res2.lambda_u
    lambda_d=res2.lambda_d  
    
    vll0=res2.vll
    vll=res2.vll   
    
    vmll=abs(vll)
    result1.vmax[id_t]=vmll.max()
    result1.vmin[id_t]=vmll.min()
    #result1.vt[:,id_t]=vmll
    
    # simulation tests: ac power flow
    mdl_par_pf.vll=vll
    mdl_par_pf.vll0=vll0
    
    # project pg_t, qg_t, qshunt_t
    pll_g_t=np.matmul(gen_to_LL,pg_t)
    qll_g_t=np.matmul(gen_to_LL,qg_t)  
        
    # shunt to sll
    # project switch shunt into integer steps
    qshunt_t=np.round(qshunt_t/Bi)*Bi 
    qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)  
    
       
    # pg       
    result1.Pgt[:,id_t]=pg_t
    result1.dPg[:,id_t]=pg_t-pg
    result1.rdPg[:,id_t]=result1.dPg[:,id_t]/pg_max
    
    result1.TdPg[id_t]=sum(abs(result1.dPg[:,id_t]))
    result1.rTdPg[id_t]=result1.TdPg[id_t]/sum(abs(pg_max))
        
    # qg 
    result1.Qgt[:,id_t]=qg_t   
    result1.dQg[:,id_t]=qg_t-qg
    qg_max_a=qg_max
    qg_max_a[np.where(qg_max_a<1e-4)]=1
    result1.rdQg[:,id_t]=result1.dQg[:,id_t]/qg_max_a
    
    result1.TdQg[id_t]=sum(abs(result1.dQg[:,id_t]))
    result1.rTdQg[id_t]=result1.TdQg[id_t]/sum(abs(qg_max))
        
    # qshunt
    result1.Qst[:,id_t]=qshunt_t 
    result1.dQs[:,id_t]=qshunt_t-qshunt
    result1.rdQs[:,id_t]=result1.dQs[:,id_t]/qshunt_max

    result1.TdQs[id_t]=sum(abs(result1.dQs[:,id_t]))      
    result1.rTdQs[id_t]=result1.TdQs[id_t]/sum(qshunt_max)
    
    result1.dsteps_Qs[:,id_t]=result1.dQs[:,id_t]/Bi
    result1.adsteps_Qs[:,id_t]=result1.adsteps_Qs[:,id_t]+result1.dsteps_Qs[:,id_t]
    result1.adsteps_Qs[:,id_t+1]=result1.adsteps_Qs[:,id_t]
    
    result1.Tdsteps_Qs[id_t]=result1.TdQs[id_t]/Bi
    result1.aTdsteps_Qs[id_t]=sum(result1.Tdsteps_Qs)
    qshunt=qshunt_t
        
    violation=0
    for i_test in range(n_sim):
    #for i_test in range(0,0):
        # generate load data
        pll_ld_r = data_rand(pll_ld, mu, sigma)
        qll_ld_r = data_rand(qll_ld, mu, sigma)
        
        # run power flow
        # update net injection
        pll_net_t=-pll_ld_r+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out 
        qll_net_t=-qll_ld_r+qll_dpv+qll_g_t+qll_shunt_t 
        sll_net_t=pll_net_t+1j*qll_net_t       

        mdl_par_pf.sll_net_t=sll_net_t

        vll_r=ACPF(alg_par_pf,mdl_par_pf)
        vmll_r=abs(vll_r)
        result1.vt[:,id_t]=vmll_r
        
        if max(vmll_r)>v_u or min(vmll_r)<v_l:
            violation+=1
           
    # RECORD_V_VIOLATION_OPF[id_t]= violation   
    # for iter in range(iter_max):
        # # derivative of voltage constraints with respect to (p, q)
        # dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
        # dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
        # dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
        # dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
        
        # # traditional generator
        # dvcnstr_dgp=dvcnstr_dp[np.ix_(busid_g_LL)]
        # dvcnstr_dgq=dvcnstr_dq[np.ix_(busid_g_LL)]
        
        # # shunt
        # dvcnstr_dshunt=dvcnstr_dq[np.ix_(busid_shunt_LL)]
        
        
        # # generators
        # pg_t=pg_t-epsi_pg*(2*alpha_pg*(pg_t-pg)+dvcnstr_dgp)
        # qg_t=qg_t-epsi_qg*(2*alpha_qg*(qg_t-qg)+dvcnstr_dgq)
        
        # # shunt
        # qshunt_t=qshunt_t-epsi_qshunt*(2*alpha_qs*(qshunt_t-qshunt)+dvcnstr_dshunt)
        
        # # minimize sgeneration from coal sgenerator 
        # # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)*w_psg_max+dvcnstr_dsp+dsf_dsp)
        # # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)*w_psg_max+dvcnstr_dsq+dsf_dsq)
        # # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp+dsf_dsp)
        # # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq+dsf_dsq) 
        
        # # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp)
        # # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq) 
        
        
        # # project1.095
        # # pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
        # # pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
        # # qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
        # # qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
        
        # # pll_g_t=np.maximum(pll_g_t,pll_g_min)
        # # pll_g_t=np.minimum(pll_g_t,pll_g_max)
        # # qll_g_t=np.maximum(qll_g_t,qll_g_min)
        # # qll_g_t=np.minimum(qll_g_t,qll_g_max)
        
        # # pg_t=np.maximum(pg_t,pg_min)
        # # pg_t=np.minimum(pg_t,pg_max)
        # # qg_t=np.maximum(qg_t,pg_min)
        # # qg_t=np.minimum(qg_t,pg_max)
        
        # # psg_t=np.maximum(psg_t,psg_min)
        # # psg_t=np.minimum(psg_t,psg_max)
        # # qsg_t=np.maximum(qsg_t,qsg_min)# only project switched shunts
        # # qsg_t=np.minimum(qsg_t,qsg_max)# only project switched shunts
        # # qsg_t[-18:]=np.maximum(qsg_t[-18:],qsg_min[-18:])# only project switched shunts
        # # qsg_t[-18:]=np.minimum(qsg_t[-18:],qsg_max[-18:])# only project switched shunts
        
        # # project
        # pg_t=np.maximum(pg_t,pg_min)
        # pg_t=np.minimum(pg_t,pg_max)
        # qg_t=np.maximum(qg_t,-qg_max)
        # qg_t=np.minimum(qg_t,qg_max)
        
        # qshunt_t=np.maximum(qshunt_t,qshunt_min)
        # qshunt_t=np.minimum(qshunt_t,qshunt_max)
       
        
        # # gen to sll
        # pll_g_t=np.matmul(gen_to_LL,pg_t)
        # qll_g_t=np.matmul(gen_to_LL,qg_t)  
        
        # # # sgen to sll
        # # pll_sg_t=np.matmul(sgen_to_LL,psg_t)
        # # qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
        
        # # shunt to sll
        # qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)  
        
        # # update net injections 
        # pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out
        # qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t
        # sll_net_t=pll_net_t+1j*qll_net_t
        
        # #pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
        

        # # fixed point iteration methods as the power flow solver       
        
        # itr_pf=0
        # itr_pf_max=50
        # err_vll=1
        # while (err_vll>1e-5 and itr_pf<itr_pf_max):
            # #vll=Z_pu*np.conj(sll_net_t/vll)+vs
            # ILL=np.conj(sll_net_t/vll)
            # dv=np.matmul(Z_pu,ILL.transpose())
            # vll=np.squeeze(np.asarray(dv))+w
            # err_vll=max(abs(vll-vll0))
            # vll0=vll
            # itr_pf=itr_pf+1        
                   
        # # if err_vll>1e-3:
        # #     raise Exception('power flow diverge!\n')
        
             
        # # # dispatch load 
        # # pld_t=pll_ld_t[busid_ld_LL]
        # # qld_t=qll_ld_t[busid_ld_LL]
        
        # # net_t.load.p_mw[busid_ld_lds]=-pld_t*sbase# load is positive in net.load
        # # net_t.load.q_mvar[busid_ld_lds]=-qld_t*sbase
        
        # # # dispatch (p,v) for generator
        # # pAg_t=pAg_net_t[busid_g_LL]
        # # vg_t=abs(vll[busid_g_LL])
        
        # # net_t.gen.p_mw=pAg_t*sbase
        # # net_t.gen.vm_pu=vg_t   
        
        # # pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)

        # # # check whether generator tracks dispatch 
        # # # load
        # # pld_t_m=-net_t.res_load.p_mw.to_numpy()/sbase
        # # pld_t_m=np.delete(pld_t_m, busid_s_ld)
        # # qld_t_m=-net_t.res_load.q_mvar.to_numpy()/sbase
        # # qld_t_m=np.delete(qld_t_m, busid_s_ld)
        
        # # # generator
        # # pg_t_m=net_t.res_gen.p_mw.to_numpy()/sbase
        # # # qg_t_m=net_t.res_gen.q_mvar/sbase
        # # vg_t_m=net_t.res_gen.vm_pu.to_numpy()
        
        # # # sgenerator
        # # # psg_t_m=net_t.res_sgen.p_mw.to_numpy()/sbase
        # # # # qg_t_m=net_t.res_gen.q_mvar/sbase
        # # # qsg_t_m=net_t.res_sgen.q_mvar.to_numpy()/sbase
        
        # # # mismatch between dispatch and measurements
        # # dpld=abs(pld_t_m-pld_t)
        # # dqld=abs(qld_t_m-qld_t)
        # # TdPg=abs(pg_t_m-pAg_t)
        # # dvg=abs(vg_t_m-vg_t)
        
        # # if dpld.max()>1e-2 or dqld.max()>1e-2 or TdPg.max()>1e-3 or dvg.max()>2e-3:
        # #     print('iteration:%d: large mismatch' %iter)
            
        # # print('maximum mismatch:')
        # # print('pld:%.8f' %dpld.max())
        # # print('qld:%.8f' %dqld.max())
        # # print('pg:%.8f' %TdPg.max())
        # # print('v:%.8f' %dvg.max())
        
        # # print('predicted:')
        # # print('max v:',abs(vll).max())
        # # print('min v:',abs(vll).min())
        
        # # vmll=net_t.res_bus.vm_pu.values[busid_LL]
        
        # # update dual variables (voltage)
        # vmll=abs(vll)       
        #result1.vt[:,id_t]=vmll
        
        # lambda_u=lambda_u+epsi_l*(vmll-v_u)
        # lambda_d=lambda_d+epsi_l*(v_l-vmll)
        
        # # project dual variables
        # lambda_u=np.maximum(lambda_u,0)
        # lambda_d=np.maximum(lambda_d,0)
        
        # if id_t==id_observe:
            # result_int.lambda_u[:,iter]=lambda_u
            # result_int.lambda_d[:,iter]=lambda_d
            
            # result_int.v[:,iter]=vmll

            # result_int.pg_t[:,iter]=pg_t        
            # result_int.qg_t[:,iter]=qg_t
                 
    # #     # track voltage and dual
    # #     result_interm1.lambda_u[iter]=lambda_u.max()
    # #     result_interm1.lambda_d[iter]=lambda_d.max()
    # #     # result_interm1.u_u[iter]=u_u.max()
    # #     # result_interm1.u_l[iter]=u_l.max()
    # #     result_interm1.vmax[iter]=vmll.max()
    # #     result_interm1.vmin[iter]=vmll.min()
        
    # #     # #update dual variable (apparent power flow)
    # #     # # nonlinear complex power flow
    # #     # v_t[busid_slack]=v0
    # #     # v_t[busid_LL]=vll
    # #     # vf_t=v_t[brh_fbus]
    # #     # vt_t=v_t[brh_tbus]
    # #     # If_t=(vf_t-vt_t)*brh_y
    # #     # Sf_t=vf_t*np.conjugate(If_t)
        
    # #     # # # linear model: Pf==Cp-Dq, Qf=Dp+Cq
    # #     # # Pi_t=np.real(sll_net_t)
    # #     # # Qi_t=np.imag(sll_net_t)
    # #     # # Pf_e=np.matmul(C,Pi_t)-np.matmul(D,Qi_t)
    # #     # # Qf_e=np.matmul(D,Pi_t)+np.matmul(C,Qi_t)
        
    # #     # # dpf=np.real(Sf_t)-Pf_e
    # #     # # dqf=np.imag(Sf_t)-Qf_e
        
    # #     # Sfm_t=abs(Sf_t)
    # #     # Pfr=np.real(Sf_t)/Sfm_t
    # #     # Qfr=np.imag(Sf_t)/Sfm_t
        
    # #     # # inactive small branch flow
    # #     # id_sf0=np.where(Sfm_t<1e-7)
    # #     # Pfr[id_sf0]=0
    # #     # Qfr[id_sf0]=0
        
    # #     # # # dual variable
    # #     # # u_u=u_u+epsi_u*(Sfm_t-Smax_brh)
    # #     # # #u_l=u_l+epsi_u*(Smin_brh-Sfm_t)
        
    # #     # # # project dual variables
    # #     # # u_u=np.maximum(u_u,0)
    # #     # # #u_l=np.maximum(u_l,0)   
        
    # #     # # result_interm1.u_u[iter]=u_u.max()
    # #     # #result_interm1.u_l[iter]=u_l.max()
    # # iterations=list(range(iter_max))

    # # print('Mismatch:',err_vll)
    # # print('max v:',vmll.max())
    # # print('min v:',vmll.min())
    
    # # # project switch shunt into integer steps
    # # qshunt_t=np.round(qshunt_t/Bi)*Bi 
    
    # # qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)
    
    # # # update net injections 
    # # pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out
    # # qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t
    # # sll_net_t=pll_net_t+1j*qll_net_t

    # # # fixed point iteration methods as the power flow solver       
    # # vll0=vll_b
    # # vll=vll_b
    # # itr_pf=0
    # # itr_pf_max=50
    # # err_vll=1
    # # while (err_vll>1e-5 and itr_pf<itr_pf_max):
        # # #vll=Z_pu*np.conj(sll_net_t/vll)+vs
        # # ILL=np.conj(sll_net_t/vll)
        # # dv=np.matmul(Z_pu,ILL.transpose())
        # # vll=np.squeeze(np.asarray(dv))+w
        # # err_vll=max(abs(vll-vll0))
        # # vll0=vll
        # # itr_pf=itr_pf+1 

    # qshunt=qshunt_t   
    # lambda_u=res2.lambda_u
    # lambda_d=res2.lambda_d
    
    # vmll=abs(vll)
    # result1.vmax[id_t]=vmll.max()
    # result1.vmin[id_t]=vmll.min()
       
    # # pg       
    # result1.Pgt[:,id_t]=pg_t
    # result1.dPg[:,id_t]=pg_t-pg
    # result1.rdPg[:,id_t]=result1.dPg[:,id_t]/pg_max
    
    # result1.TdPg[id_t]=sum(abs(result1.dPg[:,id_t]))
    # result1.rTdPg[id_t]=result1.TdPg[id_t]/sum(abs(pg_max))
        
    # # qg 
    # result1.Qgt[:,id_t]=qg_t   
    # result1.dQg[:,id_t]=qg_t-qg
    # qg_max_a=qg_max
    # qg_max_a[np.where(qg_max_a<1e-4)]=1
    # result1.rdQg[:,id_t]=result1.dQg[:,id_t]/qg_max_a
    
    # result1.TdQg[id_t]=sum(abs(result1.dQg[:,id_t]))
    # result1.rTdQg[id_t]=result1.TdQg[id_t]/sum(abs(qg_max))
        
    # # qshunt
    # result1.Qst[:,id_t]=qshunt_t 
    # result1.dQs[:,id_t]=qshunt_t-qshunt
    # result1.rdQs[:,id_t]=result1.dQs[:,id_t]/qshunt_max

    # result1.TdQs[id_t]=sum(abs(result1.dQs[:,id_t]))      
    # result1.rTdQs[id_t]=result1.TdQs[id_t]/sum(qshunt_max)
    
    # result1.dsteps_Qs[:,id_t]=result1.dQs[:,id_t]/Bi
    # result1.adsteps_Qs[:,id_t]=result1.adsteps_Qs[:,id_t]+result1.dsteps_Qs[:,id_t]
    # result1.adsteps_Qs[:,id_t+1]=result1.adsteps_Qs[:,id_t]
    
    # result1.Tdsteps_Qs[id_t]=result1.TdQs[id_t]/Bi
    # result1.aTdsteps_Qs[id_t]=sum(result1.Tdsteps_Qs)
    #qshunt=qshunt_t    
    
print('Max voltage:',max(result1.vmax)) 
print('Min voltage:',min(result1.vmin))  

xlabels_time=['00:00','06:00','12:00','18:00','23:55']
if V_ROBUSTNESS>0:
    CC='CC'
else:
    CC=''
iplt=1 
fig,ax=plt.subplots(iplt)
plt.figure(iplt)
plt.plot(RECORD_V_VIOLATION_OPF,'.')  
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel('Time')
ax.set_xlim([0,n_snaps])
ax.set_ylim([0,n_sim])
ax.set_ylabel('Voltage violation number')
ax.set_title('Voltage Violation Number')
ax.grid(True)
fig.savefig(path_plt+'/CCACOPF/'+'/'+CC+OPF_type+'VViolation.png', dpi=400) 

np.save(path_plt+'/CCACOPF/'+CC+OPF_type+'VViolation', RECORD_V_VIOLATION_OPF)

iplt+=1      
fig,ax = plt.subplots()
plt.figure(iplt)    
plot_vmax_EDCOPF=ax.plot(range(0,n_snaps), result1.vmax,'.',markersize=3,color='green')
plot_vmin_EDCOPF=ax.plot(range(0,n_snaps), result1.vmin,'.',markersize=3,color='orange')
plot_upper=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashed',color='black')
plot_lower=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashed',color='black')
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_DCOPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0]), (OPF_type+' Max V', OPF_type+' Min V'))
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel("Time")
ax.set_ylabel('Voltage magnitude, p.u.')
ax.set_xlim([0,n_snaps+1])
ax.set_ylim([0.7,1.15])
ax.set_title('Voltage')
ax.set_xlabel('Time')
ax.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/v.png', dpi=400)   


# voltage profile
iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:blue'
ax.set_xlabel('Bus index')
ax.set_ylabel('Voltage magnitude, p.u.')
for i in range(nLL):
    plot_v0=ax.plot(range(1,n_snaps+1),result1.vt[i,:],'.',markersize=1.5,color=color)
    
plot_upper0=ax.plot(range(0,n_snaps), (v_u-V_ROBUSTNESS)*np.ones(n_snaps),linestyle='dashed',color='red',linewidth=2.5)
plot_upper1=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashdot',color='m',linewidth=2.5)
plot_lower0=ax.plot(range(0,n_snaps), (v_l+V_ROBUSTNESS)*np.ones(n_snaps),linestyle='dashed',color='black',linewidth=2.5)
plot_lower1=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashdot',color='orange',linewidth=2.5)

ax.legend((plot_v0[0],plot_upper0[0],plot_upper1[0],plot_lower0[0],plot_lower1[0]), (OPF_type, 'Upper bound, 1.088 p.u.','Upper bound, 1.1 p.u.','Lower bound, 0.912 p.u.','Lower bound, 0.9 p.u.'))
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel("Time")
ax.set_ylabel('Voltage magnitude, p.u.')
ax.set_xlim([0,n_snaps+1])
ax.set_ylim([0.7,1.15])
ax.set_title('Voltage Profile')
ax.set_xlabel('Time')
ax.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/VProfile', dpi=400)

# # voltage profile
# iplt+=1
# plt.figure(iplt) 
# fig, ax = plt.subplots()
# color='tab:blue'
# ax.set_xlabel('Bus index')
# ax.set_ylabel('Voltage magnitude, p.u.')
# for i in range(nLL):
#     plot_v0=ax.plot(range(1,n_snaps+1),result1.vt[i,:],'.',markersize=1.5,color=color)
    
# plot_upper=ax.plot(range(0,n_snaps), v_u*np.ones(n_snaps),linestyle='dashed',color='red')
# plot_lower=ax.plot(range(0,n_snaps), v_l*np.ones(n_snaps),linestyle='dashed',color='black')


# ax.legend((plot_v0[0],plot_upper[0],plot_lower[0]), (OPF_type+' V', 'Upper bound','Lower bound'))
# ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
# ax.set_xticklabels(xlabels_time, rotation=0)
# ax.set_xlabel("Time")
# ax.set_ylabel('Voltage magnitude, p.u.')
# ax.set_xlim([0,n_snaps+1])
# ax.set_ylim([0.7,1.15])
# ax.set_title('Voltage Profile')
# ax.set_xlabel('Time')
# ax.grid(True)
# fig.savefig(path_plt+'/'+OPF_type+'/VProfile', dpi=400)    


iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_ylabel('Pgen, MW')
ax1.set_xlabel('Time')
#plt_1=ax1.plot(range(1,n_snaps+1),result1.TdQg*sbase,'.',markersize=2,color='r')
plt_1=ax1.plot(range(1,n_snaps+1),np.sum(result0.Pg*sbase,axis=0),'o',markerfacecolor='none',markersize=2,color='b')

plt_2=ax1.plot(range(1,n_snaps+1),np.sum(result1.Pgt*sbase,axis=0),'.',markersize=1.5,color='r')

plt_3=ax1.plot(range(1,n_snaps+1),np.ones(n_snaps)*sum(abs(pg_max))*sbase,'--',color='orange')

plt_4=ax1.plot(range(1,n_snaps+1),np.ones(n_snaps)*sum(abs(pg_min))*sbase,'--',color='gray')

plt.xlim([0,n_snaps+1])
plt.xlabel('Time')
plt.title('Generator Real Power (Pgen)')
ax1.legend((plt_1[0], plt_2[0],plt_3[0],plt_4[0]), ('Total Pgen of DC OPF', 'Total Pgen of AC OPF','Total Pgen upper bound','Total Pgen lower bound'))
plt.grid(True)

fig.savefig(path_plt+'/'+OPF_type+'/TotalPgen.png', dpi=400)  


iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_ylabel('Qgen, MVar')
ax1.set_xlabel('Time')
#plt_1=ax1.plot(range(1,n_snaps+1),result1.TdQg*sbase,'.',markersize=2,color='r')
plt_1=ax1.plot(range(1,n_snaps+1),np.ones(n_snaps)*sum(abs(qg))*sbase,'.',markersize=2,color='b')

plt_2=ax1.plot(range(1,n_snaps+1),np.sum(result1.Qgt*sbase,axis=0),'.',markersize=2,color='r')

plt_3=ax1.plot(range(1,n_snaps+1),np.ones(n_snaps)*sum(abs(qg_max))*sbase,linestyle='dashed',color='orange')

plt_4=ax1.plot(range(1,n_snaps+1),-np.ones(n_snaps)*sum(abs(qg_max))*sbase,linestyle='dashed',color='gray')

plt.xlim([0,n_snaps+1])
plt.xlabel('Time')
plt.title('Generator Reactive Power (Qgen)')
ax1.legend((plt_1[0], plt_2[0],plt_3[0],plt_4[0]), ('Total Qgen of default values', 'Total Qgen of AC OPF','Total Qgen upper bound','Total Qgen lower bound'))
plt.grid(True)

fig.savefig(path_plt+'/'+OPF_type+'/TotalQgen.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
fig, ax1 = plt.subplots()
color='tab:red'
ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax1.set_xticklabels(xlabels_time, rotation=0)
ax1.set_ylabel('Mvar',color=color)
ax1.plot(range(1,n_snaps+1),result1.TdQs*sbase,'.',markersize=2.5,color=color)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.set_xlabel('Time')
#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
ax2=ax1.twinx()
color='tab:green'
ax2.set_ylabel('Percent %', color=color)  # we already handled the x-label with ax1
ax2.plot(range(1,n_snaps+1),result1.rTdQs*100,'*',markersize=2.5, color=color)
ax2.tick_params(axis='y', labelcolor=color)
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('Shunt Reactive Power Adjustment')
plt.grid(True)

fig.savefig(path_plt+'/'+OPF_type+'/Qshunt_Adjust.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:red'
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel('Time')
ax.set_ylabel('Step change',color=color)
ax.plot(range(1,n_snaps+1),result1.Tdsteps_Qs,'.',markersize=2.5,color=color)
ax.tick_params(axis='y', labelcolor='tab:red')
#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('Switch shunt step changes')
plt.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/Number_QshuntStepChange.png', dpi=400)  

iplt+=1
plt.figure(iplt) 
fig, ax = plt.subplots()
color='tab:red'
ax.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
ax.set_xticklabels(xlabels_time, rotation=0)
ax.set_xlabel('Time')
ax.set_ylabel('Step changes',color=color)
ax.plot(range(1,n_snaps+1),result1.aTdsteps_Qs,'.',markersize=2.5,color=color)
ax.tick_params(axis='y', labelcolor='tab:red')
#ax1.set_yticks([0,0.0005,0.001,0.0015])
#plt.ylim(Padja_ylim)
plt.xlim([0,n_snaps+1])
#plt.ylim(Padjr_ylim)
plt.title('The accumulated number of switch shunt step changes')
plt.grid(True)
fig.savefig(path_plt+'/'+OPF_type+'/AccumulatedNumber_QshuntStepChange.png', dpi=400)  



# # # plot dPg at maximum TdPg
id_maxTdPg=np.argmax(result1.TdPg) 

iplt+=1
fig,ax = plt.subplots()
plt.figure(iplt)    
plot_upper=ax.plot(range(1,ngen+1), pg_max*sbase,linestyle='dashed',color='orange')
plot_lower=ax.plot(range(1,ngen+1), pg_min*sbase,linestyle='dashed',color='gray')
plot_Pg0=ax.plot(range(1,ngen+1), result0.Pg[:,id_maxTdPg]*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
plot_Pgt=ax.plot(range(1,ngen+1), result1.Pgt[:,id_maxTdPg]*sbase,'.',markersize=5,color='red')
ax.set_xlim([1,31])
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('DCOPF Pg','ACOPF Pg', 'Max Pg',' Min Pg'))
ax.set_xlabel("Generator index")
ax.set_ylabel('Real power injection, MW')
ax.set_title('Optimal Real Power Injection of Generators')
ax.grid(True) 
fig.savefig(path_plt+'/'+OPF_type+'/OptPgen_TdPmax.png', dpi=400)  
 

iplt+=1
fig,ax = plt.subplots()
plt.figure(iplt)    
plot_upper=ax.plot(range(1,ngen+1), qg_max*sbase,linestyle='dashed',color='orange')
plot_lower=ax.plot(range(1,ngen+1), -qg_max*sbase,linestyle='dashed',color='gray')
plot_Pg0=ax.plot(range(1,ngen+1), qg*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
plot_Pgt=ax.plot(range(1,ngen+1), result1.Qgt[:,id_maxTdPg]*sbase,'.',markersize=5,color='red')
ax.set_xlim([1,31])
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('Initial Qg','ACOPF Qg', 'Max Qg',' Min Qg'))
ax.set_xlabel("Generator index")
ax.set_ylabel('Reactive power injection, MVar')
ax.set_title('Optimal Reactive Power Injection of Generators')
ax.grid(True)  
fig.savefig(path_plt+'/'+OPF_type+'/OptQgen_TdPmax.png', dpi=400)  


# # plot dPg at minimum TdPg
id_minTdPg=np.argmin(result1.TdPg)


iplt+=1
fig,ax = plt.subplots()
plt.figure(iplt)    
plot_upper=ax.plot(range(1,ngen+1), pg_max*sbase,linestyle='dashed',color='orange')
plot_lower=ax.plot(range(1,ngen+1), pg_min*sbase,linestyle='dashed',color='gray')
plot_Pg0=ax.plot(range(1,ngen+1), result0.Pg[:,id_minTdPg]*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
plot_Pgt=ax.plot(range(1,ngen+1), result1.Pgt[:,id_minTdPg]*sbase,'.',markersize=5,color='red')
ax.set_xlim([1,31])
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('DCOPF Pg','ACOPF Pg', 'Max Pg',' Min Pg'))
ax.set_xlabel("Generator index")
ax.set_ylabel('Real power injection, MW')
ax.set_title('Optimal Real Power Injection of Generators')
ax.grid(True) 
fig.savefig(path_plt+'/'+OPF_type+'/OptPgen_TdPmin.png', dpi=400)  


iplt+=1
fig,ax = plt.subplots()
plt.figure(iplt)    
plot_upper=ax.plot(range(1,ngen+1), qg_max*sbase,linestyle='dashed',color='orange')
plot_lower=ax.plot(range(1,ngen+1), -qg_max*sbase,linestyle='dashed',color='gray')
plot_Pg0=ax.plot(range(1,ngen+1),qg,'o',markerfacecolor='none',markersize=5,color='blue')
plot_Pgt=ax.plot(range(1,ngen+1), result1.Qgt[:,id_minTdPg]*sbase,'.',markersize=5,color='r')
ax.set_xlim([1,31])
#ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('Initial Qg','ACOPF Qg', 'Max Qg',' Min Qg'))
ax.set_xlabel("Generator index")
ax.set_ylabel('Reactive power injection, MVar')
ax.set_title('Optimal Reactive Power Injection of Generators')
ax.grid(True)  
fig.savefig(path_plt+'/'+OPF_type+'/OptQgen_TdPmin.png', dpi=400)  


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
        
for id_t in range(n_snaps):
    iplt+=1
    fig,ax = plt.subplots()
    plt.figure(iplt)    
    plot_upper=ax.plot(range(1,ngen+1), pg_max*sbase,linestyle='dashed',color='orange')
    plot_lower=ax.plot(range(1,ngen+1), -pg_max*sbase,linestyle='dashed',color='gray')
    plot_Pg0=ax.plot(range(1,ngen+1), result0.Pg[:,id_t]*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
    plot_Pgt=ax.plot(range(1,ngen+1), result1.Pgt[:,id_t]*sbase,'.',markersize=3.5,color='red')
    ax.set_xlim([1,31])
    #ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
    ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('Initial Pg','ACOPF Pg', 'Max Pg',' Min Pg'))
    ax.set_xlabel("Generator index")
    ax.set_ylabel('Real power injection, MW')
    ax.set_title('Optimal Real Power Injection of Generators')
    ax.grid(True)  
    fig.savefig(path_plt+'/'+OPF_type+'/Pgen/OptPgen'+time_table[id_t]+'.png', dpi=400)  
    
# for id_t in range(n_snaps):
    # iplt+=1
    # fig,ax = plt.subplots()
    # plt.figure(iplt)    
    # plot_upper=ax.plot(range(1,ngen+1), qg_max*sbase,linestyle='dashed',color='orange')
    # plot_lower=ax.plot(range(1,ngen+1), -qg_max*sbase,linestyle='dashed',color='gray')
    # plot_Pg0=ax.plot(range(1,ngen+1), qg*sbase,'o', markerfacecolor='none',markersize=5,color='blue')
    # plot_Pgt=ax.plot(range(1,ngen+1), result1.Qgt[:,id_t]*sbase,'.',markersize=3,color='red')
    # ax.set_xlim([1,31])
    # #ax.legend((plot_vmax_DCOPF[0], plot_vmin_D1.1COPF[0],plot_vmax_EDCOPF[0], plot_vmin_EDCOPF[0],plot_upper[0],plot_lower[0]), ('DCOPF Max V', 'DCOPF Min V','EDCOPF Max V', 'EDCOPF Min V','Upper bound','Lower bound'))
    # ax.legend((plot_Pg0[0], plot_Pgt[0], plot_upper[0], plot_lower[0]), ('Initial Qg','ACOPF Qg', 'Max Qg',' Min Qg'))
    # ax.set_xlabel("Generator index")
    # ax.set_ylabel('Reactive power injection, MVar')
    # ax.set_title('Optimal Reactive Power Injection of Generators')
    # ax.grid(True)  
    # fig.savefig(path_plt+'/'+OPF_type+'/Qgen/OptQgen'+time_table[id_t]+'.png', dpi=400) 
    

# iplt+=1
# plt.figure(iplt) 
# fig, ax1 = plt.subplots()
# color='tab:red'
# ax1.set_xticks(list(range(0,n_snaps+1,int(n_snaps/(len(xlabels_time)-1)))))
# ax1.set_xticklabels(xlabels_time, rotation=0)
# ax1.set_ylabel('Mw',color=color)
# ax1.set_xlabel('Time')
# ax1.plot(range(1,n_snaps+1),result1.TdPg*sbase,'.',markersize=2.5,color=color)
# ax1.tick_params(axis='y', labelcolor='tab:red')
# #ax1.set_yticks([0,0.0005,0.001,0.0015])
# #plt.ylim(Padja_ylim)
# ax2=ax1.twinx()
# color='tab:green'
# ax2.set_ylabel('Percent %', color=color)  # we already handled the x-label with ax1
# ax2.plot(range(1,n_snaps+1),result1.rTdPg*100,'*',markersize=2.5, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# plt.xlim([0,n_snaps+1])
# #plt.ylim(Padjr_ylim)
# plt.title('Total Generator Real Power')
# plt.grid(True)

# fig.savefig(path_plt+'/'+OPF_type+'/Pgen_Adjust.png', dpi=400) 

# # convergenc of primal and dual at id_t=228
# id_vmin=np.argmin(result_int.v[:,-1])
# iplt+=1
# plt.figure(iplt) 
# #plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
# plot_ld=plt.plot(result_int.lambda_d[id_vmin,:],'.',markersize=3)
# plt.title('Dual Variable Convergence')
# #plt.legend((plot_lu[0], plot_ld[0]), ('$\overline{\u03BB}$', '$\underline{\u03BB}$'))
# plt.xlim([0,iter_max])
# plt.xlabel('Iteration number')
# plt.ylabel('\u03BB')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/Convergence_lambda.png', dpi=400)  

# iplt+=1
# plt.figure(iplt) 
# #plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
# plot_ld=plt.plot(result_int.v[id_vmin,:],'.',markersize=3)
# plt.title('Voltage Convergence')
# #plt.legend((plot_lu[0], plot_ld[0]), ('$\overline{\u03BB}$', '$\underline{\u03BB}$'))
# plt.xlim([0,iter_max])
# plt.xlabel('Iteration number')
# plt.ylabel('Voltage magnitude, p.u.')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/Convergence_v.png', dpi=400)  

# id_MaxPg=np.argmax(pg_max)
# iplt+=1
# plt.figure(iplt) 
# #plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
# plot_ld=plt.plot(result_int.pg_t[id_MaxPg,:]*sbase,'.',markersize=3)
# plt.title('Real Power Injection Convergence')
# #plt.legend((plot_lu[0], plot_ld[0]), ('$\overline{\u03BB}$', '$\underline{\u03BB}$'))
# plt.xlim([0,iter_max])
# plt.xlabel('Iteration number')
# plt.ylabel('Real power injection, MW')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/Convergence_pg.png', dpi=400) 

# iplt+=1
# plt.figure(iplt) 
# #plot_lu=plt.plot(result_int.lambda_u[123,:],'.',markersize=3)
# plot_ld=plt.plot(result_int.qg_t[id_MaxPg,:]*sbase,'.',markersize=3)
# plt.title('Reactive Power Injection Convergence')
# plt.xlim([0,iter_max])
# plt.xlabel('Iteration number')
# plt.ylabel('Reactive power injection, MVar')
# plt.grid(True)
# plt.savefig(path_plt+'/'+OPF_type+'/Convergence_qg.png', dpi=400) 