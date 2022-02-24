import pandapower as pp
import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

# convert matpower case file version 2 to a pandapower net.
path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)
nbus = len(net.bus)

# run power flow
pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)

# Ybus
Ybus = net._ppc['internal']['Ybus']

import math
import cmath
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

print(S)

print(I)
