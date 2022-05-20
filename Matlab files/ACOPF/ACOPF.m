clear all;
clc;
mydir=pwd;
idcs=strfind(mydir,'\');
pardir=mydir(1:idcs(end)-1);

matdir=[pardir '\output file\'];
%'Maui2022dm_rd_v33.mat'#'mpc.mat'#'mpc_maui_21Q3'#'mpc_maui.mat'#'case39.mat#'case240_cost.mat'#case9.mat#case240_21Q3
matpower_name='mpc_maui_21Q3';
%  matpower_name='Maui2022dm_rd_v33.mat';
%matpower_name='Maui2022dm_rd_v33_SwitchShuntsNoPhaseshift.mat';

mpc=importdata([matdir matpower_name '.mat']);
mpc.version='2';

% % add gencost
% ngen=size(mpc.gen,1);
% nclmn=7; 
% mpc.gencost=zeros(ngen,nclmn);
% mpc.gencost(:,1)=2;
% mpc.gencost(:,4)=3;
% mpc.gencost(:,5)=0.01;
% mpc.gencost(:,6)=0.3;
% mpc.gencost(:,6)=0.2;
% % outdir=[pardir '\output file'];
% % save([outdir '\' matpower_name],'mpc');

% adjust voltage limits
id_VMAX=12;
id_VMIN=13;
mpc.bus(:,id_VMAX)=1.092;
mpc.bus(:,id_VMIN)=0.910;

% adjust power flow limits
id_rateA=6;
id_rateB=7;
mpc.branch(:,id_rateA)=10000;% no line flow limits
mpc.branch(:,id_rateB)=10000;% no line flow limits
[RESULTS, SUCCESS]=runopf(mpc);
vm=RESULTS.bus(:,8);

outdir=[pardir '\output file'];
save([outdir '\MatpowerResult\vm_' matpower_name '_vm.m'], 'vm');