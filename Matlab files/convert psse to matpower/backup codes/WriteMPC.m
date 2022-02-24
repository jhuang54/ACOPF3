clear all;
clc;
% % converts a matpower case file version 2 to the matpower case file for pandapower.converter.from_mpc
% % matpower case file must be named as mpc, like mpc.gen, mpc.bus,
% %   mpc.branch

% step 0: address of output file '.mat' and input file 
mydir=pwd;
idcs=strfind(mydir,'\');
newdir=mydir(1:idcs(end)-1);

dir_input=[newdir '\Input file'];

% step 1: create the output file '.mat'
% case 1: input: '.m' file
% mpc = case9;
% save('case9.mat','mpc');

% mpc=case39;
% save([newdir '\case39.mat'],'mpc');

% case 2: input: '.mat' file
% %mpc=importdata('mpc_240_21Q3.mat');
% mpc=importdata('mpc_240_cost.mat');
% save('case240_cost.mat','mpc');
% name_file='mpc_maui';
% name_file='mpc_maui_21Q3';
name_file='Maui2022dm_rd_v33';
mpc=importdata([dir_input '\' name_file '.mat']);
mpc=mpc.mpc;
mpc.version='2';
save([newdir '\' name_file '.mat'],'mpc');

[groups, isolated, islandStatus] = find_islands(mpc);


 


