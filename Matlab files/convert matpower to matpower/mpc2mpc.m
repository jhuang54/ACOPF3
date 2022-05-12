clear all;
clc;
mydir=pwd;
idcs=strfind(mydir,'\');
newdir=mydir(1:idcs(end)-1);

matdir=[newdir '\mat file\'];
%'Maui2022dm_rd_v33.mat'#'mpc.mat'#'mpc_maui_21Q3'#'mpc_maui.mat'#'case39.mat#'case240_cost.mat'#case9.mat#case240_21Q3
% matpower_name='mpc_maui_21Q3.mat';
matpower_name='Maui2022dm_rd_v33.mat';
mpc=importdata([matdir matpower_name]);
mpc.version='2';

outdir=[newdir '\output file'];
save([outdir '\' matpower_name],'mpc');
