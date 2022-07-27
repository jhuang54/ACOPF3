import numpy as np
def ACOPF(alg_par,mdl_par,var):
    
    # algorithm parameter
    iter_max=alg_par.iter_max
    epsi_l=alg_par.epsi_l
    epsi_pg=alg_par.epsi_pg
    epsi_qg=alg_par.epsi_qg
    epsi_qshunt=alg_par.epsi_qshunt
    itr_pf_max=alg_par.itr_pf_max
    tol_pf=alg_par.tol_pf   
    
    # model parameters
    vll=mdl_par.vll
    vll0=mdl_par.vll0
    pg_min=mdl_par.pg_min
    pg_max=mdl_par.pg_max
    qg_min=mdl_par.qg_min
    qg_max=mdl_par.qg_max 
    qshunt_min=mdl_par.qshunt_min
    qshunt_max=mdl_par.qshunt_max   
    
    busid_g_LL=mdl_par.busid_g_LL
    busid_shunt_LL=mdl_par.busid_shunt_LL
        
    gen_to_LL=mdl_par.gen_to_LL  
    shunt_to_LL=mdl_par.shunt_to_LL
    
    nLL=mdl_par.nLL
    ngen=mdl_par.ngen
    nshunt=mdl_par.nshunt
    
    Rt=mdl_par.Rt
    Xt=mdl_par.Xt
    Z_pu=mdl_par.Z_pu
    w=mdl_par.w
    
    alpha_pg=mdl_par.alpha_pg
    alpha_qg=mdl_par.alpha_qg
    alpha_qs=mdl_par.alpha_qs
    
    pll_ld=mdl_par.pll_ld
    pll_dpv=mdl_par.pll_dpv
    pll_bat_t=mdl_par.pll_bat_t
    qll_ld=mdl_par.qll_ld
    qll_dpv=mdl_par.qll_dpv
    
    pg=mdl_par.pg
    qg=mdl_par.qg
    qshunt=mdl_par.qshunt
    
    V_ROBUSTNESS=mdl_par.V_ROBUSTNESS
    v_l=mdl_par.v_l+V_ROBUSTNESS
    v_u=mdl_par.v_u-V_ROBUSTNESS
    
    # intial primal and dual variables
    pg_t=var.pg_t
    qg_t=var.qg_t
    qshunt_t=var.qshunt_t
    lambda_u=var.lambda_u
    lambda_d=var.lambda_d
    
    # 
    
    
    class res:
        def __init__(self, nLL,ngen,nshunt):
            self.vmll=np.zeros(nLL)
            self.pg_t=np.zeros(ngen)
            self.qg_t=np.zeros(ngen)
            self.qshunt_t=np.zeros(nshunt)  
            
    res1=res(nLL,ngen,nshunt)
    
    for iter in range(iter_max):
        # derivative of voltage constraints with respect to (p, q)
        dvcnstr_dp=Rt.dot(lambda_u-lambda_d)
        dvcnstr_dp=np.squeeze(np.array(dvcnstr_dp))# convert to array
        dvcnstr_dq=Xt.dot(lambda_u-lambda_d)
        dvcnstr_dq=np.squeeze(np.array(dvcnstr_dq))# convert to array
        
        # traditional generator
        dvcnstr_dgp=dvcnstr_dp[np.ix_(busid_g_LL)]
        dvcnstr_dgq=dvcnstr_dq[np.ix_(busid_g_LL)]
        
        # shunt
        dvcnstr_dshunt=dvcnstr_dq[np.ix_(busid_shunt_LL)]
        
        
        # generators
        pg_t=pg_t-epsi_pg*(2*alpha_pg*(pg_t-pg)+dvcnstr_dgp)
        qg_t=qg_t-epsi_qg*(2*alpha_qg*(qg_t-qg)+dvcnstr_dgq)
        
        # shunt
        qshunt_t=qshunt_t-epsi_qshunt*(2*alpha_qs*(qshunt_t-qshunt)+dvcnstr_dshunt)
        
        # minimize sgeneration from coal sgenerator 
        # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)*w_psg_max+dvcnstr_dsp+dsf_dsp)
        # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)*w_psg_max+dvcnstr_dsq+dsf_dsq)
        # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp+dsf_dsp)
        # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq+dsf_dsq) 
        
        # psg_t=psg_t-epsi_psg*(2*alpha_sg*(psg_t-psg)+dvcnstr_dsp)
        # qsg_t=qsg_t-epsi_qsg*(2*alpha_sg*(qsg_t-qsg)+dvcnstr_dsq) 
        
        
        # project1.095
        # pll_ld_t=np.maximum(pll_ld_t,pll_ld_min)
        # pll_ld_t=np.minimum(pll_ld_t,pll_ld_max)
        # qll_ld_t=np.maximum(qll_ld_t,qll_ld_min)
        # qll_ld_t=np.minimum(qll_ld_t,qll_ld_max)
        
        # pll_g_t=np.maximum(pll_g_t,pll_g_min)
        # pll_g_t=np.minimum(pll_g_t,pll_g_max)
        # qll_g_t=np.maximum(qll_g_t,qll_g_min)
        # qll_g_t=np.minimum(qll_g_t,qll_g_max)
        
        # pg_t=np.maximum(pg_t,pg_min)
        # pg_t=np.minimum(pg_t,pg_max)
        # qg_t=np.maximum(qg_t,pg_min)
        # qg_t=np.minimum(qg_t,pg_max)
        
        # psg_t=np.maximum(psg_t,psg_min)
        # psg_t=np.minimum(psg_t,psg_max)
        # qsg_t=np.maximum(qsg_t,qsg_min)# only project switched shunts
        # qsg_t=np.minimum(qsg_t,qsg_max)# only project switched shunts
        # qsg_t[-18:]=np.maximum(qsg_t[-18:],qsg_min[-18:])# only project switched shunts
        # qsg_t[-18:]=np.minimum(qsg_t[-18:],qsg_max[-18:])# only project switched shunts
        
        # project
        pg_t=np.maximum(pg_t,pg_min)
        pg_t=np.minimum(pg_t,pg_max)
        qg_t=np.maximum(qg_t,-qg_max)
        qg_t=np.minimum(qg_t,qg_max)
        
        qshunt_t=np.maximum(qshunt_t,qshunt_min)
        qshunt_t=np.minimum(qshunt_t,qshunt_max)
       
        
        # gen to sll
        pll_g_t=np.matmul(gen_to_LL,pg_t)
        qll_g_t=np.matmul(gen_to_LL,qg_t)  
        
        # # sgen to sll
        # pll_sg_t=np.matmul(sgen_to_LL,psg_t)
        # qll_sg_t=np.matmul(sgen_to_LL,qsg_t)
        
        # shunt to sll
        qll_shunt_t=np.matmul(shunt_to_LL,qshunt_t)  
        
        # update net injections 
        pll_net_t=-pll_ld+pll_dpv+pll_g_t+pll_bat_t# +: in; -: out
        qll_net_t=-qll_ld+qll_dpv+qll_g_t+qll_shunt_t
        sll_net_t=pll_net_t+1j*qll_net_t
        
        #pAg_net_t=pll_g_t+pll_sg_t# real power output of aggregated Generator
        

        # fixed point iteration methods as the power flow solver       
        
        itr_pf=0
        err_vll=1
        while (err_vll>tol_pf and itr_pf<itr_pf_max):
            #vll=Z_pu*np.conj(sll_net_t/vll)+vs
            ILL=np.conj(sll_net_t/vll)
            dv=np.matmul(Z_pu,ILL.transpose())
            vll=np.squeeze(np.asarray(dv))+w
            err_vll=max(abs(vll-vll0))
            vll0=vll
            itr_pf=itr_pf+1        
                   
        # if err_vll>1e-3:
        #     raise Exception('power flow diverge!\n')
        
             
        # # dispatch load 
        # pld_t=pll_ld_t[busid_ld_LL]
        # qld_t=qll_ld_t[busid_ld_LL]
        
        # net_t.load.p_mw[busid_ld_lds]=-pld_t*sbase# load is positive in net.load
        # net_t.load.q_mvar[busid_ld_lds]=-qld_t*sbase
        
        # # dispatch (p,v) for generator
        # pAg_t=pAg_net_t[busid_g_LL]
        # vg_t=abs(vll[busid_g_LL])
        
        # net_t.gen.p_mw=pAg_t*sbase
        # net_t.gen.vm_pu=vg_t   
        
        # pp.runpp(net_t, algorithm='nr', calculate_voltage_angles=True)

        # # check whether generator tracks dispatch 
        # # load
        # pld_t_m=-net_t.res_load.p_mw.to_numpy()/sbase
        # pld_t_m=np.delete(pld_t_m, busid_s_ld)
        # qld_t_m=-net_t.res_load.q_mvar.to_numpy()/sbase
        # qld_t_m=np.delete(qld_t_m, busid_s_ld)
        
        # # generator
        # pg_t_m=net_t.res_gen.p_mw.to_numpy()/sbase
        # # qg_t_m=net_t.res_gen.q_mvar/sbase
        # vg_t_m=net_t.res_gen.vm_pu.to_numpy()
        
        # # sgenerator
        # # psg_t_m=net_t.res_sgen.p_mw.to_numpy()/sbase
        # # # qg_t_m=net_t.res_gen.q_mvar/sbase
        # # qsg_t_m=net_t.res_sgen.q_mvar.to_numpy()/sbase
        
        # # mismatch between dispatch and measurements
        # dpld=abs(pld_t_m-pld_t)
        # dqld=abs(qld_t_m-qld_t)
        # TdPg=abs(pg_t_m-pAg_t)
        # dvg=abs(vg_t_m-vg_t)
        
        # if dpld.max()>1e-2 or dqld.max()>1e-2 or TdPg.max()>1e-3 or dvg.max()>2e-3:
        #     print('iteration:%d: large mismatch' %iter)
            
        # print('maximum mismatch:')
        # print('pld:%.8f' %dpld.max())
        # print('qld:%.8f' %dqld.max())
        # print('pg:%.8f' %TdPg.max())
        # print('v:%.8f' %dvg.max())
        
        # print('predicted:')
        # print('max v:',abs(vll).max())
        # print('min v:',abs(vll).min())
        
        # vmll=net_t.res_bus.vm_pu.values[busid_LL]
        
        # update dual variables (voltage)
        vmll=abs(vll)       
        #result1.vt[:,id_t]=vmll
        
        lambda_u=lambda_u+epsi_l*(vmll-v_u)
        lambda_d=lambda_d+epsi_l*(v_l-vmll)
        
        # project dual variables
        lambda_u=np.maximum(lambda_u,0)
        lambda_d=np.maximum(lambda_d,0)
        
        # if id_t==id_observe:
            # result_int.lambda_u[:,iter]=lambda_u
            # result_int.lambda_d[:,iter]=lambda_d
            
            # result_int.v[:,iter]=vmll

            # result_int.pg_t[:,iter]=pg_t        
            # result_int.qg_t[:,iter]=qg_t 
            
    res1.pg_t=pg_t
    res1.qg_t=qg_t
    res1.qshunt_t=qshunt_t
    
    res1.vll=vll
    res1.vmll=vmll
    res1.pg_t=pg_t
    res1.qg_t=qg_t
    res1.qshunt_t=qshunt_t
    res1.lambda_u=lambda_u
    res1.lambda_d=lambda_d
    
    return res1
    
def data_rand(da, mu, sigma):
    n_da=len(da)
    rd=np.random.normal(mu, sigma, n_da) 
    da=da*(1+rd)
    return da
    
def ACPF(alg_par,mdl_par):
    tol_pf=alg_par.tol_pf
    itr_pf_max=alg_par.itr_pf_max
    
    sll_net_t=mdl_par.sll_net_t
    vll=mdl_par.vll
    vll0=mdl_par.vll0
    Z_pu=mdl_par.Z_pu
    w=mdl_par.w

    itr_pf=0
    err_vll=1
    while (err_vll>tol_pf and itr_pf<itr_pf_max):
        #vll=Z_pu*np.conj(sll_net_t/vll)+vs
        ILL=np.conj(sll_net_t/vll)
        dv=np.matmul(Z_pu,ILL.transpose())
        vll=np.squeeze(np.asarray(dv))+w
        err_vll=max(abs(vll-vll0))
        vll0=vll
        itr_pf=itr_pf+1   
    return vll