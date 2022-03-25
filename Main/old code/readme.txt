PFlow.py:
run power flow with pandapower; obtain bus admittance matrix

DispatchLd.py:
dispatch load

DispatchGenLd.py:
dispatch net_t.load, net_t.gen

DispatchSGenLd.py: 
(1) use primal and dual gradient algorithm to update load, generators and static generators
(2) dispatch net_t.load (real and reactive load), net_t.gen (real power generation, voltage) and net_t.sgen (real and reactive generation);

DispatchAGenLd.py:
(1) use primal and dual gradient algorithm to update load, generators and static generators
(2) dispatch net_t.load (real and reactive load), net_t.gen (real power generation, voltage), where net_t.gen is the aggregated generator which is the sum of one generator and static generators on the same bus. 

C:\Users\jianq\Documents\GitHub\ACOPF3\Matlab files\output file\Maui2022dm_rd_AggregateGens.mat aggregates all the generators on the same bus to one generator by
C:\Users\jianq\Documents\GitHub\ACOPF3\Matlab files\convert matpower to matpower\mpc2mpc_AggregateGens.m

