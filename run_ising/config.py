#!/usr/bin/env python
#coding=utf8

from merapy.config import (Config, CFG_ISING_BASIC, updaters_u1 as _updaters_u1 , 
        updaters_z2 as _updaters_z2, 
        copy_config)

from merapy.top_level import (top_level_product_state, 
        top_level_product_state_u1, top_level_eigenstate)
from merapy.schedule import schedule_prod_state, copy_schedule

from merapy.finite_site import finite_site_u1
import merapy.schedule as schedule_module


def make_config(h, backup_parpath=None, which_top_state='scale_invar_state',  surfix='', root1=''): 
    cfg = copy_config(CFG_ISING_BASIC)
    cfg['updaters'] = dict(_updaters_u1)
    cfg['updaters']['finite_range_func'] = finite_site_u1
    
    if which_top_state == 'scale_invar_state' : 
        pass
    elif which_top_state == 'prod_state' : 
        cfg['schedule']['schedule'] = copy_schedule(schedule_prod_state)
        #cfg['updaters']  =  list(_updaters_z2)
        cfg['updaters']['rho_top_func'] = top_level_product_state_u1 
        #cfg['updaters']['rho_top_func'] = top_level_product_state
    elif which_top_state == 'eigen_state' : 
        cfg['schedule']['schedule'] = copy_schedule(schedule_prod_state)
        cfg['updaters']['rho_top_func'] = top_level_eigenstate 
        pass 
    else: 
        raise 
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
    
    project_name =  'run-ising'
    dir_name = 'h=%s'%h 
    cfg['model_param'].update(h=h)
    Config.set_backup_parpath(cfg, project_name, fn=dir_name, backup_parpath=backup_parpath, root1=root1, surfix=surfix)
    
    if 0: 
        MERA_BACKUP_ROOT = config['BACKUP_BASE_DIR']
        MERA_BACKUP_ROOT_LOCAL = config['BACKUP_BASE_DIR_LOCAL']
        root = '/'.join([MERA_BACKUP_ROOT, root_name])
        config['backup_root']  = root
        root_local = '/'.join([MERA_BACKUP_ROOT_LOCAL, root_name])
        root1 = 'ternary/z2-symm/scale-invar' if root1 is None else root1 
        if backup_parpath is None: 
            fn = 'h=%s/'%h  
            config['backup_parpath'] = '/'.join([root, root1, fn]) 
            backup_parpath_local = '/'.join([root_local, root1, fn])     
            config['backup_parpath_local'] = backup_parpath_local
        else: 
            config['backup_parpath'] = backup_parpath
        
    return cfg 

if __name__ == '__main__': 
    c=make_config(h=0.2)
    print  c['backup_parpath']
    print  c['backup_parpath_local']
    













