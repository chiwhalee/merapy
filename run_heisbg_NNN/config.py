#!/usr/bin/env python
#coding=utf8

#from vmps.config import CONFIG_HEISENBG_BASIC, copy_config
#from vmps.config import CONFIG_HEISENBG_BASIC, copy_config, Config
from merapy.config import CFG_HEISBG_BASIC, copy_config, Config 
from merapy.top_level import top_level_product_state, top_level_product_state_u1, top_level_eigenstate
from merapy.schedule import schedule_prod_state, copy_schedule


def make_config(J2, algorithm='mera', alg_surfix='', which_top_state='scale_invar',  backup_parpath=None, root1='', surfix=''): 
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
   


