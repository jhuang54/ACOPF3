# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:56:07 2022

@author: jhuang
"""
mpc_maui
plt.figure(1)
plt_qgt=plt.plot(dqvg_a,'.')
plt_qgmax=plt.plot(qll_g_max[busid_g_LL],'*')
plt.legend(['qg', 'qg_max'])
plt.savefig(path_plt+'/Qgen_change_compare.png', dpi=400) 

plt.figure(2)
plt_qg0=plt.plot(qg)
plt_qgt=plt.plot(qg_t)
plt_qgmax=plt.plot(qll_g_max[busid_g_LL],'.')
# plt.legend(['qg0','qgt','qg_max'])
plt.legend(['qg0','qgt','qg_max'])
plt.title('Reactive power output of Qgen')
plt.savefig(path_plt+'/Qgen_compare.png', dpi=400) 

plt.figure(2)
plt_qgt=plt.plot(dqvg_a,'.')
plt.title('Change')

dqvg_a=(qg_t-qg)*sbase

qg_t=qll_g_t[busid_g_LL]
qg=net.res_gen.q_mvar.to_numpy()/sbase

dqvg_a=(qg_t-qg)*sbase

plt.figure(3)
qmax=np.concatenate((qll_g_max[busid_g_LL],qsg_max[:-18]))
plt.plot(qmax)