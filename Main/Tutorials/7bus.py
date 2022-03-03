import pandapower as pp
import pandapower.networks

net = pandapower.networks.example_simple()

pp.runpp(net)