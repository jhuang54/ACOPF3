import pandapower as pp
#import pandapower.networks

import pandapower.converter as pc

import os

from pathlib import Path

import matplotlib.pyplot as plt

from scipy.sparse.csgraph import minimum_spanning_tree

import pandas as pd

import numpy as np

path_cur = Path(os.getcwd())
path_par = path_cur.parent.absolute()
path_plt = os.path.join(path_cur, 'Plot')


path_maui=os.path.join(path_par, 'maui_dauc_rted')
geninfo_name='Gen_Info.csv'

df_geninfo = pd.read_csv (path_maui+'\\'+geninfo_name)

# convert generator bus name to busid in code
GName=df_geninfo.loc[:,'genname'].to_numpy()# gen name
GbusName=df_geninfo.loc[:,'busname'].to_numpy()# gen bus name in code

# import dispatch power of real-time economic dispatch
rted_name='RTEDOutput_Day0.xlsx'
df_rted = pd.read_excel (path_maui+'\\'+rted_name,sheet_name='Dispatch', index_col=[0])
GName_rt=df_rted.columns.values

Geninfo = {"GenName":GName,"bus name":GbusName}
Geninfo = pd.DataFrame(Geninfo)
Geninfo=Geninfo.set_index("GenName")

# Dispatchable generator name; DPV name
GName_g=[]
clnid_g=[]
GName_DPV=[]
clnid_DPV=[]
for count, gname in enumerate(GName_rt):
    if gname[0:3]!='DPV':
        GName_g.append(gname)
        clnid_g.append(count)
    else:
        GName_DPV.append(gname)
        clnid_DPV.append(count)
        
GenbusName=[Geninfo.loc[gname,'bus name'] for gname in GName_g]
Gens_tab=[]
# create generator data using the PSS\E raw file format
with open('GensDPVs.txt','r') as f:
    lines=f.readlines()
    for i in range(len(lines)-1):
        linetp=lines[i]
        linetp=linetp.split(',')
        busname=int(linetp[0])
        if busname in GenbusName:
            Gens_tab.append(lines[i])
f.close()            
with open('Gens.txt','w') as f:
     for i in Gens_tab:
        f.writelines(i)
f.close()
    
DPVbusName=[Geninfo.loc[gname,'bus name'] for gname in GName_DPV]

# read generators
DPVs_tab=np.ones((1,8),dtype='int')#busname, PG,QG,QT,QB,PT,PB,line id
with open('GensDPVs.txt','r') as f:
    lines=f.readlines()
    for i in range(len(lines)-1):
        linetp=lines[i]
        linetp = linetp.split(',')
        busname=int(linetp[0])
        if busname in DPVbusName:
            PG=float(linetp[2])
            QG=float(linetp[3])
            QT=float(linetp[4])
            QB=float(linetp[5])
            PT=float(linetp[16])
            PB=float(linetp[17])
            temp=np.array([busname,PG,QG,QT,QB,PT,PB,i])#i: line id
            DPVs_tab = np.vstack([DPVs_tab,temp])
        
        # linetp=[linetp for linetp in linetp if linetp!='']
        # if linetp[4][0:3]=='kw=':
            # load_new.append(linetp[0]+' '+linetp[1]+' '+linetp[2]+' '+linetp[3]+' kw=0.05'+' '+linetp[5]+' '+linetp[6]+' '+linetp[7])
        # else:
            # print('wrong')
f.close()
DPVs_tab=np.delete(DPVs_tab,0,0)

# lengths of PG,QG,QT,QB,PT,PB
linetp=lines[0]
linetp = linetp.split(',')
len_PG=len(linetp[2])
len_QG=len(linetp[3])
len_QT=len(linetp[4])
len_QB=len(linetp[5])
len_PT=len(linetp[16])
len_PB=len(linetp[17])

DPV_tab=np.ones((1,7),dtype='int')
n_DPVs=len(GName_DPV)
line_DPVs=[]
for i in range(n_DPVs):
    busname_tp=DPVbusName[i]
    ids_DPV=np.where(DPVs_tab[:,0]==busname_tp)
    line_id=int(DPVs_tab[ids_DPV[0][0],-1])# line id in GensDPVs.txt
    DPVs=DPVs_tab[ids_DPV,:][0]# sum of PG,QG,QT,QB,PT,PB
    DPV_tp=np.sum(DPVs,axis=0)  
    
    PG=round(DPV_tp[1],3)
    QG=round(DPV_tp[2],3)
    QT=round(DPV_tp[3],3)
    QB=round(DPV_tp[4],3)
    PT=round(DPV_tp[5],3)
    PB=round(DPV_tp[6],3)
    
    # write DPV file
    linetp=lines[line_id]
    linetp =linetp.split(',')
    # update PG,QG,QT,QB,,PT,PB
    linetp[2]=' '*(len_PG-len(str(PG)))+str(PG)
    linetp[3]=' '*(len_QG-len(str(QG)))+str(QG)
    linetp[4]=' '*(len_QT-len(str(QT)))+str(QT)
    linetp[5]=' '*(len_QB-len(str(QB)))+str(QB)
    linetp[16]=' '*(len_PT-len(str(PT)))+str(PT)
    linetp[17]=' '*(len_PB-len(str(PB)))+str(PB)
    
    linetp1=linetp[0]
    for j in linetp[1:]:
        linetp1=linetp1+','+j

    line_DPVs.append(linetp1)
    
with open('DPVs.txt','w') as f:
    for i in range(n_DPVs):
        f.writelines(line_DPVs[i])
f.close()