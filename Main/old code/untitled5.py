# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:51:08 2022

@author: jhuang
"""

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

# reactive power outpout of generators and DPV
busid_Gen=net.gen.bus.to_numpy()
Q_Gen=net.res_gen.q_mvar.to_numpy()

busid_Sgen=net.sgen.bus.to_numpy()
Q_Sgen=net.res_sgen.q_mvar.to_numpy()

busid0_GSg=np.concatenate((busid_Gen,busid_Sgen))
busid_GSg=np.sort(busid0_GSg)

id_GSg=np.argsort(busid0_GSg)

Q0_GSg=np.concatenate((Q_Gen,Q_Sgen))
n_GSg=len(busid_Gen)+len(busid_Sgen)
Q_GSg=np.zeros(n_GSg)
for i in range(n_GSg):
    Q_GSg[i]=Q0_GSg[id_GSg[i]]
Qgen_tab={}