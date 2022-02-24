clear all;
clc;
rawfile_name='Maui2022dm_rd_v33.raw';
% rawfile_name='Maui2022dm_v4_v33.raw';
verbose=1;
mpc_name='Maui2022dm_rd_v33.mat';
[mpc, warnings] = psse2mpc(rawfile_name, mpc_name, verbose);