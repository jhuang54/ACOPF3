import pandapower as pp
import pandapower.networks

import pandapower.converter as pc
# convert matpower case file version 2 to a pandapower net.
#'Maui2022dm_rd_v33.mat'#'mpc.mat'#'mpc_maui_21Q3'#'mpc_maui.mat'#'case39.mat#'case240_cost.mat'#case9.mat#case240_21Q3
icase='Maui2022dm_rd_v33.mat'
path_icase='C:\\Users\\jianq\\Dropbox\\Research\\Stability-constrainted AC-OPF\\Matlab files\\output file\\'
#net=pc.from_mpc(path_icase+icase,f_hz=60)
net=pc.from_mpc(path_icase+icase,f_hz=60)
nbus=len(net.bus)

# run power flow
pp.runpp(net,algorithm='nr',calculate_voltage_angles=True)

# Ybus
Ybus=net._ppc['internal']['Ybus']


import math
import cmath
import numpy as np
vm=net.res_bus.vm_pu.values
va_dg=net.res_bus.va_degree.values
va=np.exp(1j*va_dg*np.pi/180)
v=vm*va

mappd2ppc=net._pd2ppc_lookups["bus"]
Ybus=Ybus[mappd2ppc,0:nbus]
Ybus=Ybus[0:nbus,mappd2ppc]
I=Ybus.dot(v)
S=v*np.conj(I)