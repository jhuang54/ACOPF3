clear all;
clc;
mydir=pwd;
idcs=strfind(mydir,'\');
newdir=mydir(1:idcs(end)-1);

pssedir=[newdir '\psse raw file\'];

% rawfile_name='Maui2022dm_rd_v33.raw';
%rawfile_name='Maui2022dm_rd_v33_shunt.raw';
rawfile_name='Maui2022dm_rd_v33_shunt_OnlyLoad.raw';
% rawfile_name='Maui2022dm_v4_v33.raw';
verbose=1;
%mpc_name='Maui2022dm_rd_v33_shunt.mat';
mpc_name='Maui2022dm_rd_v33_shunt_OnlyLoad.mat';

[mpc, warnings] = psse2mpc([pssedir rawfile_name], mpc_name, verbose);

mpc.version='2';

% revise r, x, rateA, rateB of branches
clmn_pool=[3,4,6,7];
bnd_pool=[1e-4,1e-4,100,100];
tol=1e-10;
for i=1:length(clmn_pool)
    clmn_tp=clmn_pool(i);
    rw_tp=find(abs(mpc.branch(:,clmn_tp))<tol);
    mpc.branch(rw_tp,clmn_tp)=bnd_pool(i);
end


outdir=[newdir '\output file'];
save([outdir '\' mpc_name],'mpc');

% check isolated buses
[groups, isolated, islandStatus] = find_islands(mpc);