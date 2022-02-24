# import pandapower as pp
# net=pp.create_empty_network()
# b1=pp.create_bus(net,vn_kv=20.)
# b2=pp.create_bus(net,vn_kv=20.)
# pp.create_line(net,from_bus=b1,to_bus=b2,length_km=2.5,std_type='NAYY 4x50 SE')
# pp.create_ext_grid(net,bus=b1)
# pp.create_load(net,bus=b2,p_mw=1)

# pp.runpp(net)

# print(net.res_bus.vm_pu)
# print(net.res_line.loading_percent)

# import pandapower
# import pandapower.networks
# import pandapower.topology
# import pandapower.plotting
# import pandapower.converter
# import pandapower.estimation

# import pandapower.test 
# pandapower.test.run_all_tests()

import pandapower as pp
#create empty net
net=pp.create_empty_network()

# create buses
b1=pp.create_bus(net,vn_kv=20.,name='Bus 1')
b2=pp.create_bus(net,vn_kv=0.4,name='Bus 2')
b3=pp.create_bus(net,vn_kv=0.4,name='Bus 3')

#create bus elements
pp.create_ext_grid(net,bus=b1,vm_pu=1.02,name='Grid Connection')
pp.create_load(net, bus=b3,p_mw=0.1,q_mvar=0.05,name='load')

# create branch elements
trafo=pp.create_transformer(net,hv_bus=b1,lv_bus=b2,std_type='0.4 MVA 20/0.4 kV',name='Trafo')
line=pp.create_line(net,from_bus=b2,to_bus=b3,length_km=0.1,name='Line', std_type='NAYY 4x50 SE')

pp.runpp(net)

net.trafo.tap_pos.at[trafo]=-1
pp.runpp(net)

net.res_bus

pp.create_switch(net,bus=b3,element=line,et='l',closed=False)

pp.runpp(net)
net.res_bus
net.res_load