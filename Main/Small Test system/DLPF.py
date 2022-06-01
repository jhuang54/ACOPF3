import pandapower.networks as pn

import pandapower as pp

import numpy as np

import pandas as pd

from numpy.linalg import inv

net = pn.case4gs()
nbus = len(net.bus)

pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)

# reference bus
vm0=net.ext_grid.vm_pu.to_numpy()
va0=net.ext_grid.va_degree.to_numpy()*np.pi/180
# base value:
# apparent power: net.sn_mva
# voltage: net.bus.vn_kv
# consider shunt elements in diagonal elements
sbase=net.sn_mva

Ybus = net._ppc['internal']['Ybus']
mappd2ppc = net._pd2ppc_lookups["bus"]
Ybus = Ybus[mappd2ppc, 0:nbus]
Ybus = Ybus[0:nbus, mappd2ppc]
Ybus=Ybus.todense()

# Do Not consider shunt elements in diagonal elements
Ybusp = net._ppc['internal']['Ybus']
mappd2ppc = net._pd2ppc_lookups["bus"]
Ybusp = Ybusp[mappd2ppc, 0:nbus]
Ybusp = Ybusp[0:nbus, mappd2ppc]
Ybusp=Ybusp.todense()

# bus
busid_slack=net.ext_grid.bus[0]

busid_LL=list(range(0,nbus))
busid_LL.remove(busid_slack)

nLL=nbus-1

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

# Aitf
Aitf=np.zeros((nbrh,nLL))
for ibus in range(nLL):
    busid_tp=busid_LL[ibus]
    idbrh_f=np.where(brh_fbus==busid_tp)#branch id of branches whose fbus is busid_tp
    idbrh_t=np.where(brh_tbus==busid_tp)
    Aitf[idbrh_f,ibus]=1
    Aitf[idbrh_t,ibus]=-1

# gl,bl
gl=np.real(brh_y)
gl=np.matmul(np.diag(gl),Aitf)

bl=np.imag(brh_y)
bl=np.matmul(np.diag(bl),Aitf)

# G, B
G=np.real(Ybus)
B=np.imag(Ybus)

GLL=G[busid_LL,:]
GLL=GLL[:,busid_LL]

BLL=B[busid_LL,:]
BLL=BLL[:,busid_LL]

# G', B'
for ibus in range(nbus):
    Ybusp[ibus,ibus]=0
    Ybusp[ibus,ibus]=-np.sum(Ybusp[ibus,:])

Gp=np.real(Ybusp)
Bp=np.imag(Ybusp)

GLLp=Gp[busid_LL,:]
GLLp=GLLp[:,busid_LL]

BLLp=Bp[busid_LL,:]
BLLp=BLLp[:,busid_LL]

# H, N, M, L
H=-BLLp
N=GLL
M=-GLL
L=-BLL

Htl=H-np.matmul(N,np.matmul(inv(L),M))
Ltl=L-np.matmul(M,np.matmul(inv(H),N))

# pl=A*Ptl+B*Qtl
Ritf=-(np.matmul(gl,np.matmul(inv(Ltl),np.matmul(M,inv(H)))))
Xitf=np.matmul(gl,inv(Ltl))+np.matmul(bl,np.matmul(inv(Htl),np.matmul(N,inv(L))))

# pnet
Pnet=-net.res_bus.p_mw.to_numpy()/sbase
Pnet=np.delete(Pnet, busid_slack, 0)

Qnet=-net.res_bus.q_mvar.to_numpy()/sbase
Qnet=np.delete(Qnet, busid_slack, 0)

# Ptl, Qtl
BLRp=Bp[busid_LL,busid_slack]
GLR=G[busid_LL,busid_slack]
BLR=B[busid_LL,busid_slack]
Ptl=Pnet+np.matmul(BLRp,va0)-np.matmul(GLR,vm0)
Ptl=Ptl.reshape((nLL,1))
Qtl=Qnet+np.matmul(GLR,va0)+np.matmul(BLR,vm0)
Qtl=Qtl.reshape((nLL,1))

# Pftl, Qftl
Pftl=np.matmul(Ritf,Ptl)+np.matmul(Xitf,Qtl)

# theta, v
angtl=np.matmul(inv(Htl),Ptl)-np.matmul(inv(Htl),np.matmul(N,np.matmul(inv(L),Qtl)))
vmtl=np.matmul(inv(Ltl),Qtl)-np.matmul(inv(Ltl),np.matmul(M,np.matmul(inv(H),Ptl)))