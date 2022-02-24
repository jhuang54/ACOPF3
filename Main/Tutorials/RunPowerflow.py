import pandapower as pp
import pandapower.networks

net=pandapower.networks.example_simple()

pp.runpp(net)

net.res_bus

net.res_bus[net.bus.vn_kv==20.].vm_pu.min()

load_or_generation_buses=set(net.load.bus.values) | set(net.sgen.bus.values) | set(net.gen.bus.values)
net.res_bus.vm_pu.loc[load_or_generation_buses].max()
