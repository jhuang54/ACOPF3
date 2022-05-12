# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 22:24:44 2022

@author: jianq
"""

# import numpy as np
# import pandapower.plotting as pt
# import pandapower as pp
# import pandapower.networks as pn
# import plotly
# import plotly.graph_objs as go
# import plotly.graph_objs as go
# import plotly.io as pio
# #pio.renderers.default = 'browser'

# net = pn.create_kerber_vorstadtnetz_kabel_1()
# fig = pt.simple_plotly(net)

# how_much_buses = 5
# color_buses = np.random.choice(net.bus.index, how_much_buses)

# #color_buses = np.random.choice(range(5), how_much_buses)

# fig.add_trace(go.Scatter(x=net.bus_geodata.loc[color_buses, 'x'],
#                           y=net.bus_geodata.loc[color_buses, 'y'],
#                           mode='markers'))

# fig.show()

# import pandapower.networks
# net = pandapower.networks.mv_oberrhein("generation")
# pandapower.plotting.simple_plot(net)

import pandapower as pp

import pandapower.converter as pc

import pandapower.networks

import os

from pathlib import Path

import pandas as pd

import numpy as np
import pandapower.plotting as pt

from pandapower.plotting import simple_plot, simple_plotly, pf_res_plotly

path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_plt = os.path.join(path_cur, 'Plot')

path_output = os.path.join(path_par, 'Matlab files\output file')
icase = 'Maui2022dm_rd_v33.mat'
net = pc.from_mpc(path_output + '\\' + icase, f_hz=60)# initial condition


import seaborn
colors=seaborn.color_palette()
color_tb=[]
nbus=len(net.bus)
for i in range(nbus):
    color_tb.append(colors[7])

id_shunt=np.array([6,8,24,57,59,61,62,69,88,91,93,110,112,138,141,148,154,167])
id_shunt=list(id_shunt-1)
for i in id_shunt:
    color_tb[i]=colors[9]
 
color_tb[96]=colors[3]    
simple_plot(net, bus_color=color_tb)



# buscolor={'blue'}
# for i in range(nbus-1):
#     buscolor.add('blue')
# pp.plotting.simple_plot(net,bus_color)
#pp.plotting.simple_plot(net,bus_color=buscolor)
# color_buses = np.random.choice(net.bus.index, len(net.bus.index))
# color_trace = pt.plotly.create_bus_trace(net, color_buses, color='red',
#                                           trace_name='special buses')
# pt.plotly.draw_traces(color_trace)

# from pandapower.plotting import cmap_discrete, create_line_trace, draw_traces
# cmap_list = [((20, 50), "green"), ((50, 70), "yellow"), ((70, 100),"red")]
# cmap, norm = cmap_discrete(cmap_list)
# lc = create_line_trace(net, cmap=cmap)
# draw_traces([lc])

# import numpy as np
# import pandapower.plotting as pt
# import pandapower as pp
# import pandapower.networks as pn

# net = pn.create_kerber_vorstadtnetz_kabel_1()
# pt.simple_plotly(net)

# color_buses = np.random.choice(net.bus.index, 33)
# color_trace = pt.plotly.create_bus_trace(net, color_buses, color='red',
#                                           trace_name='special buses')
# pt.plotly.draw_traces(color_trace)