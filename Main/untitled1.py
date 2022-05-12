# -*- coding: utf-8 -*-
"""
Created on Sun May  1 10:22:21 2022

@author: jianq
"""

# import pandapower.networks as nw
# import pandapower.plotting as plot
# import matplotlib.pyplot as plt
# import seaborn
# colors=seaborn.color_palette()
# net=nw.mv_oberrhein()
# bc=plot.create_bus_collection(net,buses=net.bus.index,color=colors[0],size=80,zorder=1)
# lc=plot.create_line_collection(net,lines=net.line.index,color='grey',zorder=2)
# plt.show()

# import pandapower.networks as nw
# from pandapower.plotting import simple_plot, simple_plotly, pf_res_plotly
# net=nw.mv_oberrhein()
# simple_plot(net)
# simple_plotly(net)
# pf_res_plotly(net)

# import pandapower.networks as nw

# from pandapower.plotting import simple_plot, simple_plotly, pf_res_plotly, create_bus_collection, create_line_collection

# %matplotlib inline

# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt

# import seaborn

# colors=seaborn.color_palette()
# net=nw.mv_oberrhein()
# bc=create_bus_collection(net,buses=net.bus.index,color=colors,size=80,zorder=1)
# lc=create_line_collection(net,lines=net.line.index,color='grey',zorder=2)
# plt.show()

# from pandapower.plotting.plotly import pf_res_plotly
# from pandapower.networks import mv_oberrhein
# # net = mv_oberrhein()
# # pf_res_plotly(net)

# net = mv_oberrhein()
# pf_res_plotly(net, on_map=True, projection='epsg:31467', map_style='dark')

import pandapower.networks as nw
from pandapower.plotting import simple_plot, simple_plotly, pf_res_plotly
net=nw.mv_oberrhein()

import seaborn
colors=seaborn.color_palette()
color_tb=[]
nbus=len(net.bus)
for i in range(nbus):
    color_tb.append(colors[7])
color_tb[6]=colors[3]
simple_plot(net, bus_color=color_tb)


