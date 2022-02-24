clear all;
clc;
%'Maui2022dm_rd_v33.mat'#'mpc.mat'#'mpc_maui_21Q3'#'mpc_maui.mat'#'case39.mat#'case240_cost.mat'#case9.mat#case240_21Q3
matpower1_name='mpc_maui_21Q3.mat';
mpc1=importdata([matpower1_name]);
branch1=mpc1.branch;
matpower2_name='mpc_maui.mat';
mpc2=importdata([matpower2_name]);
branch2=mpc2.branch;
dbrh=branch1;
dbrh(:,3:end)=branch1(:,3:end)-branch2(:,3:end);
d3=find(abs(dbrh(:,3))>1e-10);
d4=find(abs(dbrh(:,4))>1e-10);
d5=find(abs(dbrh(:,5))>1e-10);

d6=find(abs(dbrh(:,6))>1e-10);
d7=find(abs(dbrh(:,7))>1e-10);
d8=find(abs(dbrh(:,8))>1e-10);

% [branch1(d3,3) branch2(d3,3) branch1(d4,4) branch2(d4,4) branch1(d6,6) branch2(d6,6) branch1(d7,7) branch2(d7,7)]
d9=find(abs(dbrh(:,9))>1e-10);
d10=find(abs(dbrh(:,10))>1e-10);
d11=find(abs(dbrh(:,11))>1e-10);

d12=find(abs(dbrh(:,12))>1e-10);
d13=find(abs(dbrh(:,13))>1e-10);


matpower3_name='Maui2022dm_rd_v33.mat';
mpc3=importdata([matpower3_name]);