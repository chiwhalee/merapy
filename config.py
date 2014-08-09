import unittest
import os, argparse
import platform
import socket

if 1: 
    from merapy.ascending import ascending_ham
    from merapy.descending import descending_ham
    from merapy.iteration import iterative_optimize_all
    from merapy.all_in_once import update_all_in_once
    from merapy.top_level import top_level_product_state, top_level_product_state_u1, top_level_eigenstate
    from merapy.finite_site import finite_site, finite_site_u1
    #import merapy.finite_site as finite_site_module
    
    from merapy.context_util import LOCAL_IP, LOCAL_USERNAME 
  
import merapy.updaters_binary as upbin

__all__ = ["CFG_ISING_BASIC", "CFG_HEISBG_BASIC", "CFG_POTTS_BASIC", "ising_model_param", "heisbg_model_param", 
        "updaters_z2", "updaters_u1", "updaters_binary", 
        "cfg_graph", "cfg_binary", "cfg_modified_binary", "cfg_quaternary", "cfg_quinary"]

uname=platform.uname()[1]
home = os.path.expanduser('~')
if uname  ==  'QTG-WS1-ubuntu': 
    MERA_BACKUP_DIR = '/home/zhli/Documents/mera_backup_tensor' 
else: 
    MERA_BACKUP_DIR = '/'.join([home, 'mera_backup_tensor'])

BACKUP_BASE_DIR =  '/'.join([home, 'backup_tensor_dir'])
BACKUP_BASE_DIR_LOCAL = '/home/zhli/backup_tensor_dir' #'/'.join([home, 'backup_tensor_dir'])

MERA_BACKUP_DIR_LOCAL = '/home/zhli/Documents/mera_backup_tensor' 


updaters_u1 = {"ascending_func":ascending_ham, "descending_func":descending_ham, "update_mera_func":iterative_optimize_all, 
        "finite_range_func":finite_site_u1, "rho_top_func":top_level_product_state_u1}

updaters_binary = updaters_u1.copy()
updaters_binary.update({"ascending_func":upbin.ascending_ham, "descending_func":upbin.descending_ham, "update_mera_func":upbin.iterative_optimize_all})


updaters_z2= [ascending_ham, None, update_all_in_once,finite_site, top_level_product_state]

#updaters_u1= [ascending_ham, descending_ham, iterative_optimize_all, finite_site_u1, top_level_product_state_u1]



ising_model_param = {"h":-1.0, "J_NN":-1.0, "J_NNN":0.0}
heisbg_model_param = {"J_NN":1.0, "J_NNN":0.241186}

import collections
def recursive_update_cfg(d, u):
    """
        this will be used later 
    """
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

def copy_config(cfg): 
    """
        a simple deep copy of config
    """
    res= dict(cfg)
    #model_param = dict(cfg['model_param'])
    model_param = dict(cfg.get('model_param', {}))
    res['model_param'] = model_param
    return res

#meta cfg for everything  
CFG_BASE = {
    'NUM_OF_THREADS':1, 
    #for dist comput
        'hostname': socket.gethostname(), 
        'pid': os.getpid(), 
        'LOCALHOSTNAME': 'QTG-WS1-ubuntu',  
        'LOCAL_USERNAME': LOCAL_USERNAME, 
        'LOCAL_IP': LOCAL_IP, 
        'BACKUP_BASE_DIR':BACKUP_BASE_DIR,       
        'BACKUP_BASE_DIR_LOCAL':BACKUP_BASE_DIR_LOCAL,       
        'register_job': 0, 

        'use_local_storage': False,  #store and save to local (local means center, my computer)
       
        }

#issue: only_NN and only_NNN are ambiguous, they dont specify which layer.  this may cause problem when dealing with e.g. binary graph, in which
#in bottom layer it should be only_NN, but in higher layers should be only_NNN

CFG_MERA = copy_config(CFG_BASE)
CFG_MERA.update({
        'algorithm': 'mera', 
        'algorithm_surfix': '', 
        
        'USE_CUSTOM_RAND':False,  #always need try different seeds 
        'USE_REFLECTION': False,
        'RAND_SEED':1234,  #this is seed for custom rand
        'use_player':True, 
        'backup_fn': None,
        'filename': None,
        'model_param': {}, 
        'only_NN':True, 
        'only_NNN': False,
        'q_iter': 110,
        'q_lay': 1,
        'message': None,
        'info': 0,
        'tot_layer': 4,
        'trunc_dim': 4,
        
        'run_schedule': False,
        
        'backup_parpath': None, 
        'backup_parpath_local': None,   #used for distributive compute
        'backup_relpath': [], 
        'backup_root': None, 
         
        'do_measure': False, 
        'measurement_args': {
            'exclude_which': [
                'entanglement_brute_force', 
                'entanglement_brute_force_6',  
                'entanglement_brute_force_9', 
                'entanglement_brute_force_9_aver'
                ], 
            'which': None, 
            }, 
        
        #'schedule':{
        #    'schedule': None, 
        #    'mera_shape_min': None, 
        #    'mera_shape_max': None, 
        #    'dim_diff_remap': {}, 
        #    } , 
        #'BACKUP_BASE_DIR': MERA_BACKUP_DIR,         
        #'BACKUP_BASE_DIR_LOCAL': MERA_BACKUP_DIR_LOCAL,         
        'MERA_BACKUP_DIR': MERA_BACKUP_DIR, 
        'MERA_BACKUP_DIR_LOCAL': MERA_BACKUP_DIR_LOCAL, 
 })



#CFG_ISING_BASIC = CFG_MERA.copy()
CFG_ISING_BASIC = copy_config(CFG_MERA)
CFG_ISING_BASIC.update({
     # model specific
     'MODEL': 'Ising',
     'USE_REFLECTION': False,
     'SYMMETRY': "Z2", 
     'model_param':{"h":-1.0, "J_NN":-1.0, "J_NNN":0.0, 'gamma':1.0},
     'unitary_init': 'unit_tensor',
     'updaters':None, 
     'combine_2site':False, 
     'trunc_dim': 4,
     })




CFG_HEISBG_BASIC = copy_config(CFG_MERA)
CFG_HEISBG_BASIC.update({
    'unitary_init': 'random_unit_tensor',
    'updaters':None, 
    #model specific
    'combine_2site':True, 
    'MODEL': 'Heisenberg',
    'SYMMETRY': 'U1',
    'model_param': {'J_NN': 1.0, 'J_NNN': 0.0, 'alpha': 2.0, 'beta': 1.0, 'Jzz': 1.0},
    'energy_exact':-1.73 , 
     })


CFG_POTTS_BASIC = CFG_MERA.copy()
CFG_POTTS_BASIC.update({
    'unitary_init': 'unit_tensor',
    'updaters':None, 
    'trunc_dim': 3,  
    #model specific
    'combine_2site':False, 
    'MODEL': 'Potts',
    'SYMMETRY': 'Z3',
    'model_param':dict(h=1.0, J_NNN=0.0), 
     })


if 1:
    def _cfg_modified_binary():
        from merapy.diagrams.modified_binary import graph_modified_binary
        tensor_defs = {
                "V":{"type":(2, 1)}, 
                "V_dag":{"type":(1, 2)}, 
                }
        mera_kwargs= {"tensor_defs":tensor_defs}
        sys_kwargs= {"graph_module":graph_modified_binary}
        new = {"mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs}
        #current_cfg.update(new)
        return new
    cfg_modified_binary = _cfg_modified_binary()

    def _cfg_binary():
        from merapy.diagrams.V21 import graph_binary
        tensor_defs = {
                "V":{"type":(2, 1)}, 
                "V_dag":{"type":(1, 2)}, 
                }
        mera_kwargs= {"tensor_defs":tensor_defs}
        sys_kwargs= {"graph_module":graph_binary, }
        new = {"updaters":updaters_binary, "mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs, "only_NN":False, "only_NNN":True}
        return new
    cfg_binary = _cfg_binary()

    def _cfg_quaternary():
        from merapy.diagrams.V41 import graph_quaternary 
        tensor_defs = {
                "V":{"type":(4, 1)}, 
                "V_dag":{"type":(1, 4)}, 
                }
        mera_kwargs = {"tensor_defs":tensor_defs}
        sys_kwargs= {"graph_module":graph_quaternary}
        new = {"mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs}
        return new
    cfg_quaternary = _cfg_quaternary()

    def _cfg_quinary():
        from merapy.diagrams.V51 import graph_quinary as graph_module
        tensor_defs = {
                "V":{"type":(5, 1)}, 
                "V_dag":{"type":(1, 5)}, 
                }
        mera_kwargs = {"tensor_defs":tensor_defs}
        sys_kwargs= {"graph_module":graph_module}
        new = {"mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs}
        return new
    cfg_quinary = _cfg_quinary()

    def _cfg_septenary():
        from merapy.diagrams.V71 import graph_septenary as graph_module
        tensor_defs = {
                "V":{"type":(7, 1)}, 
                "V_dag":{"type":(1, 7)}, 
                }
        mera_kwargs = {"tensor_defs":tensor_defs}
        sys_kwargs= {"graph_module":graph_module}
        new = {"mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs}
        return new
    cfg_septenary = _cfg_septenary()
    
def cfg_graph(which):
    dic = {
            'binary':         { 'parpath':'V21',                }, 
            'modified_binary':{ 'parpath':'modified_binary',    }, 
            'ternary':        { 'parpath':'V31',                }, 
            'quaternary':     { 'parpath':'V41',                }, 
            'quinary':        { 'parpath':'V51',                }, 
            'septenary':      { 'parpath':'V71',                }, 
            }
    temp="from merapy.diagrams.%(parpath)s import graph_%(which)s as graph_module"%dict(
            parpath=dic[which]['parpath'], which=which)
    exec(temp)

    G_2_2 = graph_module.G_2_2
    k = G_2_2.keys()[0]
    vnode = G_2_2[k].find_node('V')[0]
    vrank = G_2_2[k].nodes[vnode]
    tensor_defs = {
            "V":{"type":(vrank-1, 1)}, 
            "V_dag":{"type":(1, vrank-1)}, 
            }
    print 'using %(which)s graph, V and V_dag are defined as %(tensor_defs)s'%vars()
    mera_kwargs = {"tensor_defs":tensor_defs}
    sys_kwargs = {"graph_module":graph_module}
    new = {"mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs}

    if which == 'binary':
        new.update(only_NN=False, only_NNN=True)

    return new



if 1: 
    MERA_CMD_LINE_PARSER = argparse.ArgumentParser(description='generic args for mera')
    #MERA_CMD_LINE_PARSER.add_argument('-a', '--alpha', type=float, nargs=1, default=None)
    #MERA_CMD_LINE_PARSER.add_argument('-a', '--alpha', type=float, required=True)
    MERA_CMD_LINE_PARSER.add_argument('-dir', '--dir', default=None)
    MERA_CMD_LINE_PARSER.add_argument('-schedule', '--schedule', default=None)

    #ARGS_LONG = MERA_CMD_LINE_PARSER.parse_args()
    #ARGS_LONG = vars(ARGS_LONG)
   


def check_config(cfg): 
    """
        remind which param need to be set
    """

# eventually,  I would make config into a class
class Config(dict): 
    @staticmethod
    def set_backup_parpath_old(cfg, project_name, fn, backup_parpath=None, root1='', surfix=''): 
        if cfg['algorithm'] == 'mera':  
            ROOT = cfg['MERA_BACKUP_DIR']
            ROOT_LOCAL = cfg['MERA_BACKUP_DIR_LOCAL']
        elif cfg['algorithm'] in ['vmps', 'idmrg']: 
            ROOT = cfg['MPS_BACKUP_DIR']
            ROOT_LOCAL = cfg['MPS_BACKUP_DIR_LOCAL']
        else: 
            raise
           
            #ROOT = cfg['BACKUP_BASE_DIR']
            #ROOT_LOCAL = cfg['BACKUP_BASE_DIR_LOCAL']
        root = '/'.join([ROOT, project_name, root1]) 
        root_local = '/'.join([ROOT_LOCAL, project_name, root1])  
        cfg['root'] = root
        cfg['root_local'] = root_local
        if backup_parpath != 0: #CONVENTION: 0 means not set it
            if backup_parpath is None :  
                #fn =  'alpha=%s'%alpha  + surfix
                backup_parpath = '/'.join([root, fn]).replace('//', '/')
                backup_parpath_local = '/'.join([root_local, fn]).replace('//', '/')
                
            cfg['backup_parpath'] = backup_parpath
            cfg['backup_parpath_local'] = locals().get('backup_parpath_local')
            
            if cfg['hostname'] != cfg['LOCALHOSTNAME'] : 
                cfg['use_local_storage'] = 1
    
    @staticmethod       
    def set_backup_parpath(cfg, project_name, fn, backup_parpath=None, root1='', surfix=''): 
        ROOT = cfg['BACKUP_BASE_DIR']
        ROOT_LOCAL = cfg['BACKUP_BASE_DIR_LOCAL']
        alg = cfg['algorithm']
        alg_sur = cfg['algorithm_surfix']
        if alg_sur: 
            alg = '-'.join([alg, alg_sur])
        root = '/'.join([ROOT, project_name, alg,  root1]) 
        root_local = '/'.join([ROOT_LOCAL, project_name,alg,  root1])  
        cfg['root'] = root
        cfg['root_local'] = root_local
        #surfix  = '-'+surfix if surfix != ''  else ''
        if surfix != '': 
            fn = '-'.join([fn, surfix])
        if backup_parpath != 0: #CONVENTION: 0 means not set it
            if backup_parpath is None :  
                #fn =  'alpha=%s'%alpha  + surfix
                backup_parpath = '/'.join([root, fn ]).replace('//', '/')
                backup_parpath_local = '/'.join([root_local, fn]).replace('//', '/')
                
            cfg['backup_parpath'] = backup_parpath
            cfg['backup_parpath_local'] = locals().get('backup_parpath_local')
            
            if cfg['hostname'] != cfg['LOCALHOSTNAME'] : 
                cfg['use_local_storage'] = 1
           

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    def test_copy_config(self): 
        pass
    def test_temp(self): 
        pass
    def test_set_backup_parpath(self): 
        cfg = copy_config(CFG_ISING_BASIC)
        Config.set_backup_parpath(cfg, 'proj', fn='h=1.0', root1='middle', surfix='surf')
        path  = cfg['backup_parpath_local']
        a = '/home/zhli/backup_tensor_dir/proj/mera/middle/h=1.0-surf'
        print path 
        self.assertTrue(path==a)
        

if __name__ == '__main__': 
    pass

    if 1: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #TestIt('test_temp'), 
           TestIt('test_ising'), 
           TestIt('test_heisbg') 
        ]
        for a in add_list: 
            suite.addTest(a)
        #suite.addTest(TestIt('test_ising'))
        #suite.addTest(TestIt('test_heisbg'))
        unittest.TextTestRunner().run(suite)
       

