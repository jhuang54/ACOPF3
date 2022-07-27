# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 07:31:29 2022

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
icase = 'Maui2022dm_rd_v33_shunt_QGen.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition
pp.runpp(net, algorithm='nr', calculate_voltage_angles=True)