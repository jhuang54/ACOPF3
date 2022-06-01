# -*- coding: utf-8 -*-
"""
Created on Mon May 30 23:20:49 2022

@author: jhuang
"""
# create a psse raw file does NOT have generator, but have load on every bus

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


# convert matpower case file version 2 to a pandapower net.
# path_cur = Path(os.getcwd())
# path_par = path_cur.parent.absolute()
path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33_shunt.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition

nbus = len(net.bus)
Loadfile=os.path.join(path_par, 'maui_dauc_rted\FullLoads.dss')
#13,'1 ',1,   1,  13,     1.711,     0.437,     0.000,     0.000,     0.000,     0.000, 613,1,0
with open(Loadfile,'w') as f:
    for i in range(nbus-1):
        busname=str(net.bus.name[i]+1)
        n_spc=6-len(busname)
        spc=" "*n_spc
        f.writelines(spc+busname+',\'1 \',1,   1,  13,     0.000,     0.000,     0.000,     0.000,     0.000,     0.000,1\n')
f.close()

    