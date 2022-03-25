clear all;
clc;
% aggregate multiple generators on one bus to one generator and output 'Maui2022dm_rd_AggregateGens.mat'
mydir=pwd;
idcs=strfind(mydir,'\');
newdir=mydir(1:idcs(end)-1);

matdir=[newdir '\output file\'];
%'Maui2022dm_rd_v33.mat'#'mpc.mat'#'mpc_maui_21Q3'#'mpc_maui.mat'#'case39.mat#'case240_cost.mat'#case9.mat#case240_21Q3
matpower_name='Maui2022dm_rd_v33.mat';
mpc=importdata([matdir matpower_name]);
% mpc.version='2';
% 
% outdir=[newdir '\output file'];
% save([outdir '\' matpower_name],'mpc');

% aggregate generators on one bus to one generator
% status
status=mpc.gen(:,8);
id_actgen=find(status>0);
n_actgen=length(id_actgen);
actgen=mpc.gen(id_actgen,:);

tab=zeros(n_actgen,size(mpc.gen,2));
ngen_tp=0;
id.vg=6;
id.status=8;
for i=1:n_actgen
    % bus id
    id_bus=actgen(i,1);
    [~,id_rw]=ismember(id_bus,tab(:,1));
    if id_rw>0
        % bus was registered
        vg_tp=tab(id_rw,id.vg);
        tab(id_rw,2:end)=tab(id_rw,2:end)+actgen(i,2:end);
        
        tab(id_rw,id.vg)=vg_tp;
        tab(id_rw,id.status)=1;
    else
        % Not registered yet
        ngen_tp=ngen_tp+1;
        tab(ngen_tp,:)=actgen(i,:);  
    end
end
tab=tab(1:ngen_tp,:);
mpc.gen=tab;
outdir=[newdir '\output file'];
matpower_aggregate_name='Maui2022dm_rd_AggregateGens.mat';
save([outdir '\' matpower_aggregate_name],'mpc');