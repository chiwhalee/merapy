from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
import unittest
import warnings
import os, argparse
import platform
import socket
import collections

CFG_LOCAL = None
try:
    from .config_local import  config as CFG_LOCAL
except:
    warnings.warn('config_local.py not found')



#LOCAL_IP = '210.45.74.88'
#LOCAL_USERNAME = 'zhli' 
#LOCAL_HOSTNAME = 'qtgc30'

LOCAL_IP = '192.168.50.1'
#LOCAL_IP = '202.200.98.230'

LOCAL_USERNAME = 'ws' 
LOCAL_HOSTNAME = 'ws-Precision-Tower-7910'



__all__ = ["CFG_ISING_BASIC", "CFG_HEISBG_BASIC", "CFG_POTTS_BASIC", "ising_model_param", "heisbg_model_param", 
        'gen_backup_base_dir', 
        "cfg_graph", "cfg_binary", "cfg_modified_binary", "cfg_quaternary", "cfg_quinary"]

def gen_backup_base_dir():
    """
        dynamically generate it,  depend on machines
    """
    HOME = os.path.expanduser('~')
    HOME=HOME.replace('\\', '/') 
    if 'cygwin' in HOME:  #properly deal with cygwin path
        assert 'cygwin64/home' in HOME
        #HOME = HOME.replace('cygwin64/home', 'Users')
        HOME = HOME.replace('cygwin64/home', 'Users')
        HOME  += '/Dropbox' 
    res = '/'.join([HOME, 'backup_tensor_dir'])
    return res

BACKUP_BASE_DIR =  gen_backup_base_dir()

BACKUP_BASE_DIR_LOCAL =  '/'.join(['', 'home', LOCAL_USERNAME, 'backup_tensor_dir'])


ising_model_param = {"h":-1.0, "J_NN":-1.0, "J_NNN":0.0}
heisbg_model_param = {"J_NN":1.0, "J_NNN":0.241186}

def recursive_update_cfg(d, u):
    """
        this will be used later 
    """
    for k, v in u.items():
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
        'LOCAL_HOSTNAME': LOCAL_HOSTNAME ,  
        'LOCAL_USERNAME': LOCAL_USERNAME, 
        'LOCAL_IP': LOCAL_IP, 
        #'BACKUP_BASE_DIR':BACKUP_BASE_DIR,       
        'BACKUP_BASE_DIR_LOCAL':BACKUP_BASE_DIR_LOCAL,       
        'register_job': 0, 

        'use_local_storage': False,  #store and save to local (local means center, my computer)
       
        }

if CFG_LOCAL is not None:
    CFG_BASE.update(CFG_LOCAL)

#issue: only_NN and only_NNN are ambiguous, they dont specify which layer.  this may cause problem when dealing with e.g. binary graph, in which
#in bottom layer it should be only_NN, but in higher layers should be only_NNN


CFG_MERA = copy_config(CFG_BASE)
CFG_MERA.update({
        'algorithm': 'mera', 
        'algorithm_surfix': '', 
        'updaters':{ "ascending_func":None, "descending_func":None, "update_mera_func":None, "finite_range_func":None, "rho_top_func":None}, 
        
        'USE_CUSTOM_RAND':False,  #always need try different seeds 
        'USE_REFLECTION': False,
        'rand_seed':1234,  #this is seed for custom rand
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
        
        #'run_schedule': False,
        
        'backup_parpath': None, 
        'backup_parpath_local': None,   #used for distributive compute
        'parpath_relative':None, 
         
        'do_measure': 1, 
        
        'schedule':{
            #'schedule': None, 
            #'schedule': schedule_module.schedule_scale_invar,  
            'mera_shape_min': (4, 4), 
            'mera_shape_max': (12, 4), 
            'dim_diff_remap': {}, 
            } , 
        
        #'BACKUP_BASE_DIR_LOCAL': MERA_BACKUP_DIR_LOCAL,         
        #'MERA_BACKUP_DIR': MERA_BACKUP_DIR, 
        #'MERA_BACKUP_DIR_LOCAL': MERA_BACKUP_DIR_LOCAL, 
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
     #'updaters':None, 
     'combine_2site':False, 
     'trunc_dim': 4,
     })


CFG_HEISBG_BASIC = copy_config(CFG_MERA)
CFG_HEISBG_BASIC.update({
    'unitary_init': 'random_unit_tensor',
    #'updaters':updaters_u1, 
    #'updaters':'updaters_u1', 
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
    #'updaters':None, 
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
        #sys_kwargs= {"graph_module":graph_binary, }
        sys_kwargs= {"graph_module":'merapy.diagrams.V21.graph_binary', }
        new = {
                #"updaters":updaters_binary, 
                #"updaters":'updaters_binary', 
                "mera_kwargs":mera_kwargs, "sys_kwargs":sys_kwargs, "only_NN":False, "only_NNN":True}
        return new
    cfg_binary = _cfg_binary()

    def _cfg_quaternary():
        from merapy.diagrams.V41 import graph_quaternary 
        #print "set QspU1.MaxQNNum to 30"
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
    k = list(G_2_2.keys())[0]
    vnode = G_2_2[k].find_node('V')[0]
    vrank = G_2_2[k].nodes[vnode]
    tensor_defs = {
            "V":{"type":(vrank-1, 1)}, 
            "V_dag":{"type":(1, vrank-1)}, 
            }
    print('using %(which)s graph, V and V_dag are defined as %(tensor_defs)s'%vars())
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
    def set_backup_parpath(cfg, project_name, fn, root1='', surfix=''): 
        """
            naming rule for path of storage 
        """
        alg = cfg['algorithm']
        alg_sur = cfg['algorithm_surfix']
        if alg_sur: 
            alg = '-'.join([alg, alg_sur])
        
        fn=fn.replace('=-', '=m')
        if surfix != '': 
            fn = '-'.join([fn, surfix])
            
        if 0: 
            root = '/'.join([BACKUP_BASE_DIR,  project_name, alg,  root1]) 
            root_local = '/'.join([cfg['BACKUP_BASE_DIR_LOCAL'], project_name,alg,  root1])  
            backup_parpath = '/'.join([root, fn ]).replace('//', '/')
            backup_parpath_local = '/'.join([root_local, fn]).replace('//', '/')
                
            cfg['backup_parpath'] = backup_parpath
            cfg['backup_parpath_local'] = backup_parpath_local
        else:
            cfg['parpath_relative']=  '/'.join([project_name, alg,  root1, fn]) 
            
     
    @staticmethod 
    def filter(cfg, an=None,  db_class=None, 
            from_energy_rec=True, strict=False, info=0):
        if db_class is None:
            from merapy.measure_and_analysis.result_db import (ResultDB, ResultDB_idmrg, 
                    ResultDB_vmps, ResultDB_tdvp)
            alg = cfg['algorithm']
            db_class = ResultDB.algorithm_name_to_rdb(alg)
        if cfg['parpath_relative'] is None:
             return True
        parpath = '/'.join([ BACKUP_BASE_DIR, cfg['parpath_relative']])
        if platform.system()=='Linux':
            parpath = parpath.replace('backup_tensor_dir', 'resultdb_dir')
        else:
            parpath = parpath.replace('backup_tensor_dir', 'Dropbox/resultdb_dir')
        if an is None:
            db=db_class(parpath)
        else:
            dir_name = os.path.basename(parpath)
            alpha = an.parse_dir_name(dir_name)
            if alpha in an.alpha_parpath_dict:
                db = an[alpha]
                if dir_name !=  os.path.basename(db.parpath):
                    a = os.path.dirname(cfg['parpath_relative'])
                    b = os.path.basename(db.parpath)
                    parpath_relative_new = '/'.join([a, b])
                    print('\tchanged parpath_relative from "{}" to identicle \n\tone "{}"'.format(dir_name, b))
                    cfg['parpath_relative'] = parpath_relative_new
            else:
                db=db_class(parpath)
                
        allow = True 
        N = 0 if alg == 'idmrg' else cfg['N']
        
        msg = [os.path.basename(parpath), 'N=%d'%N,  
                'threads=%d'%cfg['NUM_OF_THREADS']]
        if alg == 'vmps' :
            which_minimize = cfg.get('which_minimize', '1site')
            v = db.fetch_easy('variance', (N, 'max'), default=10.0)
            Dmax = db.get_dim_max_for_N(N)
            schedule = cfg['schedule']
            if cfg.get('auto_resume', True):
                schedule = [s for s in schedule if s['D'] > Dmax]
                cfg['schedule'] = schedule 
            
            if which_minimize =='1site':  #old for 1site algrithm and save file for each D 
                if  (len(schedule)==0 or 
                    abs(v) <= cfg['variance_lim']):
                    allow = False
                if info>0:
                    print('\tvariance',  '%1.2e'%v)
                    
            else:
                Dmax_schedule = max([s['D'] for s in schedule])
                
                tem = db.fetch_easy('run_info', (N, 'max'), 
                        sub_key_list=['trunc_err_max'])
                tem = tem if tem is not None else 1
                #if (tem <= cfg['trunc_err_TOL'] or 
                #        v <= cfg['variance_lim'] or 
                #        Dmax_schedule <= Dmax ):
                #    allow = False
                if info>0:
                    print('\tvariance',  '%1.2e'%v)
                    print('\ttrunc_err',  '%1.2e'%tem) 
                
                conditions = [
                        tem <= cfg['trunc_err_TOL'], 
                        abs(v) <= cfg['variance_lim']
                        ]
                if strict:
                    if all(conditions):
                        allow  = False
                else:
                    if any(conditions):
                        allow = False
                    
        elif alg == 'idmrg':
            schedule = cfg['schedule']
            Dmax = max(schedule)
            sh = (0, Dmax)
            if db.has_shape(sh, from_energy_rec=from_energy_rec):
                allow = False
        elif alg == 'tdvp' :
            trunc_dim, the_time_lim, tlim = cfg['trunc_dim'], cfg['the_time_lim'], cfg['the_time_lim']
            msg.insert(2, 'td=%d, tlim=%s'%(trunc_dim, tlim))
            sh = (N, trunc_dim)
            the_time = db.fetch_easy('the_time', sh, default=0)
            if info>0:
                print(the_time, db.get('the_time'))
            if abs(the_time ) >=  the_time_lim:
                allow = False 
            else:
                allow = True
        else:
            raise  
        
        
        if not cfg.get('auto_resume', True):
            allow = True
        if cfg.get('test_run', False):
            allow = True
            
        if not allow:
            msg = ['FOUND '] + msg[:2]
        else:
            msg = ['ADD '] + msg 
            
        print(' '.join(msg))
           
        return allow 

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
        #path  = cfg['backup_parpath_local']
        path  = cfg['parpath_relative']
        path = '/'.join([cfg['BACKUP_BASE_DIR_LOCAL'], path])
        a = '/'.join(['', 'home', LOCAL_USERNAME, 'backup_tensor_dir/proj/mera/middle/h=1.0-surf']) 
        print('path', path)
        print('a', a)
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
       

