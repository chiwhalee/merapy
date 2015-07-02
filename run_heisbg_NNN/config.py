#!/usr/bin/env python
#coding=utf8

#from vmps.config import CONFIG_HEISENBG_BASIC, copy_config
#from vmps.config import CONFIG_HEISENBG_BASIC, copy_config, Config
from merapy.config import CFG_HEISBG_BASIC, copy_config, Config 
from merapy.top_level import top_level_product_state, top_level_product_state_u1, top_level_eigenstate
from merapy.schedule import schedule_prod_state, copy_schedule

if 1: 
    from merapy.quantum_number_py import init_System_QSp, symmetry_to_Qsp 
    from merapy.tensor_py import iTensor, iTensorFactory 

    def op_heisenberg_pinning(**param):      
        """
            the ground state of the spin-1/2 chain can be constructed exactly by the Bethe ansatz; we therefore know its ground state energy exactly. In the thermodynamic limit  the energy per site is given by
            E0 = 1 / 4 − ln2 = − 0.4431471805599...        
            4*E0 = 1.772588722236
        """
        symmetry = param["symmetry"]
        qspclass =symmetry_to_Qsp(symmetry)
        #QN_idendity, QSp_base = param["QN_idendity"], param["QSp_base"]

        QN_idendity, QSp_base, qsp_null = init_System_QSp(symmetry=symmetry, 
            combine_2site=param['combine_2site'])
     
        Jzz = param["Jzz"]  #Jzz has only implemented for NN interaction
        J_NN = param["J_NN"]
        J_NNN = param["J_NNN"]
        
        only_NNN = param['only_NNN']
        only_NN = param['only_NN']


        h0, Sz, Sp, Sm = {}, {}, {}, {}
        h0[2]=iTensor(4, QSp_base.copy_many(4, reverse=[2, 3]), QN_idendity.copy())
        h0[2].data[:] = 0.0
        h0[3]=iTensor(6, QSp_base.copy_many(6, reverse=[3, 4, 5]), QN_idendity.copy())
        h0[3].data[:] = 0.0

        from math import pi,cos
        momentum = pi/9.0*0.0
        twist = cos(momentum)

        if 1:
            pauli = iTensorFactory.pauli_mat_2site(symmetry=symmetry)
            if symmetry in ["Z2"]:
                ii, szz, sxx, syy= pauli["ii"], pauli["szz"], pauli["sxx"], pauli["syy"]
                sigma_z1, sigma_z2 = pauli["sigma_z1"], pauli["sigma_z2"]
                sigma_x1, sigma_x2 = pauli["sigma_x1"], pauli["sigma_x2"]
                sigma_y1, sigma_y2 = pauli["sigma_y1"], pauli["sigma_y2"]
                
            elif symmetry in ["U1", "Travial"]:
                ii, szz, spm, smp = pauli["ii"], pauli["szz"], pauli["spm"], pauli["smp"]
                sigma_z1, sigma_z2 = pauli["sigma_z1"], pauli["sigma_z2"]
                sigma_p1, sigma_p2 = pauli["sigma_p1"], pauli["sigma_p2"]
                sigma_m1, sigma_m2 = pauli["sigma_m1"], pauli["sigma_m2"]

            #ss =  sz_i ×sz_{i+1} + 2.0*(sp_i × sm_{i+1} + sm_i ×sp_{i+1} )
            #note sz, sp, sm are Pauli mat. not spin ops, so
            #这是标准的能量为E = -0.443的反铁磁海森堡模型ss的4倍 = -1.773...

        if symmetry  in ["Z2"]:
            #s0s1 = Jzz*szz + 2.0*(spm + smp) 
            ss = sxx + (1.0)*syy + Jzz*szz 
            #print 'xxxx', ss, exit()
            ssii = ss.direct_product(ii)
            iiss= ii.direct_product(ss)
            issi = (sigma_x2.direct_product(sigma_x1)  
                    + (-1.0)* sigma_y2.direct_product(sigma_y1)  
                    +  Jzz*sigma_z2.direct_product(sigma_z1))
            if 0: 
                print 'replaced szz  '*100
                ss = szz
                ssii = ss.direct_product(ii)
                iiss = ii.direct_product(ss)
                issi = sigma_z2.direct_product(sigma_z1)
               
          
        elif symmetry in ["U1", "Travial"]:
            ss = Jzz*szz + 2.0*twist*(spm + smp)
            ssii = ss.direct_product(ii)
            iiss = ii.direct_product(ss)
            issi = Jzz*sigma_z2.direct_product(sigma_z1) + 2.0*twist*(
                        sigma_p2.direct_product(sigma_m1) + 
                        sigma_m2.direct_product(sigma_p1))
        
            #print 'may be problem on ii' * 50
       
        s0s1 = 0.5*issi + 0.25*(ssii + iiss)  #final h0 is the average of ssii, issi, iiss with the RIGHT weights(1/4, 1/2, 1/4)  #这一weight的原因是，想象N个spin的链子连接到U1的mera下面，则连接ssii, issi, iiss方式出现的比率为1:2:1, 除以4即次weight
        s0s1 *= J_NN
        h0[2] += s0s1 
        
        
        if 1:  #next neareast neighbour interaction
            if symmetry in ["Z2"]:
                sisi = sigma_z1.direct_product(sigma_z1) + 2.0*(
                        sigma_x1.direct_product(sigma_x1) + 
                        sigma_y1.direct_product(sigma_y1))

                isis = sigma_z2.direct_product(sigma_z2) + 2.0*(
                        sigma_x2.direct_product(sigma_x2) + 
                        sigma_y2.direct_product(sigma_y2))
                

            elif symmetry in ["U1", "Travial"]:
                sisi = sigma_z1.direct_product(sigma_z1) + 2.0*(
                        sigma_p1.direct_product(sigma_m1) + 
                        sigma_m1.direct_product(sigma_p1))

                isis = sigma_z2.direct_product(sigma_z2) + 2.0*(
                        sigma_p2.direct_product(sigma_m2) + 
                        sigma_m2.direct_product(sigma_p2))
                
            s0s2 = 0.5*(sisi + isis)
            s0s2 *=  J_NNN
            h0[2] += s0s2
            
        #this is a temporary workaround for meta stable state for J1J2 model 
        #pinning_term = s0s2.copy()
        #pinning_term.data = -pinning_term.data 
        # todo: pinning_term is better defined in config.py 
        pinning_term = {}
        pinning_term['sisi'] = s0s2.copy()
        pinning_term['issi'] = 0.5*issi 
        pinning_term['ssii'] = 0.25*(ssii + iiss)
                
        if symmetry in ["Z2"]:
            #identity= sigma_0.direct_product(sigma_0)
            identity= ii.direct_product(ii)
        
        elif symmetry in ["U1", "Travial"]:
            identity=ii.direct_product(ii)
            
        h0[2].data -= identity.data
        
        return   pinning_term 


def make_config(J2, algorithm='mera', alg_surfix='', which_top_state='scale_invar_state',  backup_parpath=None, root1='', surfix=''): 
    cfg= copy_config(CFG_HEISBG_BASIC)


    if which_top_state == 'scale_invar_state' : 
        pass
    elif which_top_state == 'prod_state' : 
        cfg['schedule']['schedule'] = copy_schedule(schedule_prod_state)
        cfg['updaters']['rho_top_func'] = top_level_product_state_u1 
    elif which_top_state == 'eigen_state' : 
        cfg['schedule']['schedule'] = copy_schedule(schedule_prod_state)
        cfg['updaters']['rho_top_func'] = top_level_eigenstate 
        pass 
    else: 
        raise 
    
    cfg['model_name'] = 'heisbg_NNN'
    cfg.update(algorithm=algorithm, algorithm_surfix=alg_surfix, combine_2site=1)
    cfg['model_param'].update(J_NN=1.0, J_zz=1.0, J_NNN=J2)
    cfg['do_measure'] = 1
    cfg['measurement_args'] = {
            'exclude_which': [
                'entanglement_brute_force', 
                #'entanglement_brute_force_6',  
                'entanglement_brute_force_9', 
                'entanglement_brute_force_9_aver'
                ], 
            'which': None, 
            } 
    
    
    if 1:
        exact = None 
        if J2 == 0.2411:   
            exact = -1.607784749
        if J2 == 0.5: 
            exact = -1.5 
        cfg['energy_exact'] = exact 
  
    project_name = 'run-heisbg-NNN'
    dir_name = 'J2=%s'%J2 
    
    Config.set_backup_parpath(cfg, project_name, fn=dir_name, backup_parpath=backup_parpath, root1=root1, surfix=surfix)
        
    return cfg
   
if __name__ == '__main__'    : 
    cfg = make_config(J2=1.5,  root1='', which_top_state='scale_invar_state')
    aaa = ['backup_parpath', 'backup_parpath_local']
    for a in aaa: 
        print a,  cfg[a]
   


