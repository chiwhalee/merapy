#!/usr/bin/env python
#coding=utf8
#PYTHON_ARGCOMPLETE_OK 
import numpy as np 
import unittest
import argparse, argcomplete
import sys
import os
import pprint
import cPickle as pickle
from collections import OrderedDict
import warnings
import dispy 
from  multiprocessing import Process
import tempfile 

from merapy.hamiltonian import System
from merapy.utilities import load 
from merapy.context_util import rpyc_load 
#from merapy.brokest import queue 
#from mera_wigner_crystal import ham_wigner_crystal

#from vmps.mps import MPS
if 0: 
    from merapy.measure_and_analysis.central_charge import * #central_charge, entanglement_entropy, entanglement_spectrum, entanglement_extra
    from merapy.measure_and_analysis.correlation import correlation, correlation_extra
    from merapy.measure_and_analysis.scaling_dimension import calc_scaling_dim_2site as scaling_dim
    from merapy.measure_and_analysis.magnetization import magnetization
#import merapy.measure_and_analysis.all as all_mera
from merapy.utilities import print_vars
import merapy.measure_and_analysis as mmm 
import merapy as kkk 
import merapy.measure_and_analysis.all as all_mera
import vmps.measure_and_analysis.all as all_mps
try: 
    import vmps.measure_and_analysis.measurement_idmrg as all_idmrg 
    import vmps.measure_and_analysis.measurement_idmrg_mcc as all_idmrg_mcc
except ImportError as err: 
    print err 

from merapy.measure_and_analysis.result_db import ResultDB, FIELD_NAME_LIST
from merapy.decorators import timer


#this class should be moved to another file
class Measurement(object):
    """
        the object to be measured has to be a System instance stored in pickle file
    """
    def __init__(self, path):
        #self.parpath = os.path.pardir(path)
        path = os.path.abspath(path)
        self.parpath = os.path.dirname(path)  + '/'
        self.rdb = ResultDB(self.parpath)


def measure_decorator(result_db_class=None): 
    """
        status: 
            incomplete
        complete this in future 
    """
    raise NotImplemented
    def inner(func):
        def wrapper(*args, **kwargs):
           
            res = func(*args, **kwargs)
            

            return retval
        return wrapper
    return inner

def measure_S(S=None, parpath=None, path=None, which=None, exclude_which=None, force=0, 
        fault_tolerant=1, param=None,  field_surfix='', rdb=None, **kwargs):
    """
        todo: 
            parpath should be deprecated, extract it from path instead  
    """
    use_local_storage =  kwargs.get('use_local_storage', False)
    if S is None: 
        assert path is not None
        #S= load(path)
        S = rpyc_load(path, use_local_storage=use_local_storage)
        parpath = os.path.dirname(os.path.abspath(path))
    
    if rdb is None: 
        rdb = ResultDB(parpath, dbname=None, use_local_storage=use_local_storage)
    is_converged = True
    if 1: 
        #if isinstance(S, MPS): 
        if S.__class__.__name__ == 'MPS':                  
            algorithm = 'mps'
        #elif isinstance(S, np.ndarray): 
        elif isinstance(S, dict): 
            if S.has_key('mps'): 
                algorithm = 'mps'
            elif S.has_key('A'): 
                algorithm = 'idmrg'
        else: 
            algorithm = 'mera'
        
        if algorithm == 'mera':  
            field = list(FIELD_NAME_LIST)
            all_func = all_mera
            if S.model == 'wigner_crystal': 
                db_version = 1.0 
                #rdb['version'] = 1.0
            else: 
                db_version = None
           
            if S.symmetry == 'U1': 
                try: 
                    field.remove('entanglement_brute_force_9')
                    field.remove('entanglement_brute_force_9_aver')
                except: 
                    pass
            
            propt = S.key_property()
            dim,  layer = propt['trunc_dim'], propt['num_of_layer']
            iter = S.iter
            fn = System.backup_fn_gen(dim=dim, layer=layer-1, nqn = len(propt['qns']))
            path = parpath + '/' + fn
            
            if db_version is None: 
                k = (dim, layer, iter)
                if not rdb.has_key(k): 
                    rdb[k] = OrderedDict()
            
            shape = (dim, layer)
            corr_param = dict(direction=None, force_update=False, distance_max=10**4)
        
        elif algorithm == 'mps': 
            field = ['correlation']
            all_func = all_mps
            S = S['mps']  #fuck, this is bad, but I like bad!
            mps= S 
            db_version = 1.0
            N, D = mps.N, mps.D
            if mps.symmetry is not None : 
                D = mps.bond_dim_max 
            fn = "N=%d-D=%d.pickle"%(N, D)
            path = '/'.join([parpath, fn])
            dim,  layer = None, None
            shape =  N,  D
            iter = -1
            #field = ['correlation']
            corr_param = {}
            
        elif algorithm == 'idmrg': 
            field = ['energy', 'correlation', 'correlation_length', 'magnetization', 
                    'entanglement', 'entanglement_spectrum']
            if S.get('which_minimize')=='psi_guess': 
                diff = S.get('energy_diff_2')
            elif S.get('which_minimize')=='lam_guess': 
                diff = S.get('energy_diff')
            else: 
                diff = None
            diff_tol = S.get('energy_diff_tol')
            if diff is not None and diff_tol is not None : 
                if not abs(diff) <= abs(diff_tol) and not force: 
                    is_converged = False 
                    #print 'state is not converged. not allowing measure. return None'
                    #return 
            
            if S.has_key('lam_prev_inv'): 
                all_func = all_idmrg_mcc
            else: 
                all_func = all_idmrg
            A = S['A']
            symm = None if isinstance(A, np.ndarray) else A.symmetry 
            db_version = 1.0
            #N, D = 0, A.shape[0]
            N = 0
            if symm is None: 
                D = A.shape[0]
            else: 
                D = A.shape[0].totDim 
            fn = "N=%d-D=%d.pickle"%(N, D)
            path = '/'.join([parpath, fn])
            dim,  layer = None, None
            shape =  N, D 
            #iter = S['N']
            iter = -1
            #field = ['correlation']
            corr_param = {}
   
    which = which if which is not None else field 
    if exclude_which is not None: 
        print  'exlude this from measure_S: ', exclude_which
        which = [i for i  in which if i not in exclude_which]
    
    func = {i: all_func.__getattribute__(i) for i in which}
    
    if param is None:     
        param = {f: {} for f in which}
        param.update(correlation=corr_param)
    else: 
        for f in which: 
            if not param.has_key(f):  
                param[f] = {}
                
    func_name_map = {
            'calc_scaling_dim_2site': 'scaling_dim', 
            }
    
   
    #------------------------  start measureing ----------------------------------------#
    
    print  'meassuring %s at %s %s'%(path, shape, iter)
    if not is_converged: 
        print '\tstate is not converged. not allowing measure. return None'
        return 

    failed_list = []
    for w in which: 
        ff = func[w]
        w = func_name_map.get(ff.__name__, ff.__name__)
        field_name = w if field_surfix  == ''  else w + '_' +  field_surfix  
        print '\t%s...   '%field_name,
        
        if db_version is None: 
            is_found = rdb[k].has_key(field_name)
        else:
            key_list = [field_name] + [shape, iter]
            is_found = rdb.has_entry(key_list)
        if is_found and not force: 
            msg ='%20s'%('found in db')
        else: 
            try:
                res = ff(S, **param[w])
                msg = '%20s'%('done')
                if kwargs.get('show', False): 
                    pprint.pprint(res, indent=1, width=200)
                #if w in ['entanglement_entropy', 'entanglement_spectrum', 'correlation_extra']: 
                if db_version is None: 
                    rdb.put(field_name=field_name, key=k, res=res, sub_key_list=None)
                else: 
                    rdb.insert(field_name=field_name, sh=shape, iter=iter, val=res)
                    
            except Exception as err:
                failed_list.append((w, shape, iter, err))
                msg = '%20s'%('FAILED skip')
                if not fault_tolerant: 
                    raise 
        print msg
    
    if algorithm == 'mera' :
        print '\tenergy...     done'
        if db_version is None: 
            rdb[k].update({'energy': S.energy})
        elif db_version == 1.0:
            if not rdb.has_key('energy'): rdb.add_key_list(['energy', shape])
            rdb['energy'][shape] = S.energy_record
    elif algorithm == 'mps': 
        rdb['version'] = 1.0
        rdb['algorithm'] = 'mps'
    elif algorithm == 'idmrg': 
        rdb['version'] = 1.0
        rdb['algorithm'] = 'idmrg'
        #temp=all_idmrg.config(S)
        #rdb.update(temp)
        
    rdb.commit(info=1)
    if len(failed_list)>0: 
        print 'failed_list: %s'%failed_list


#now I undertand that, this func in fact realized a decorator for individal measure funcs!!
def measure_all_bac(dir_list, mera_shape_list, sh_min=None, sh_max=None,  
        which=None, exclude_which=None, param=None,  force=0, field_surfix='', 
        fault_tolerant=1, recursive=False,  **kwargs): 
    
    use_local_storage = kwargs.get('use_local_storage', False)
    path_not_found = []
    __mera_shape_list = mera_shape_list
    if isinstance(which, str): 
        which = [which]
    if isinstance(dir_list, str): 
        dir_list = [dir_list]
    if recursive: 
        dir_list_1 = list(dir_list)
        dir_list = []
        for dir in dir_list_1: 
            temp=mera_backup_dir_finder(root=dir)
            dir_list.extend(temp)
        #print 'recursively found these dirs %s'%(dir_list, )
    aaa = os.path.abspath(dir_list[0])
    if 'mera' in aaa: 
        algorithm = 'mera'
    elif 'mps' in aaa and 'idmrg' not in aaa: 
        algorithm = 'mps'
    elif 'idmrg' in aaa:
        algorithm = 'idmrg'
    else: 
        raise
    
    for dir in dir_list: 
        #fori d in dim_list: 
        if mera_shape_list  == 'all': 
            db = ResultDB(dir, algorithm=algorithm)
            __mera_shape_list = db.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
        
        for sh in __mera_shape_list: 
            if algorithm  == 'mera':  
                #sh[1] -= 1  #layer->layer-1
                sh = list(sh);  sh[1] -= 1
                fn = System.backup_fn_gen(*sh)
            elif algorithm == 'mps': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            elif algorithm == 'idmrg': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            else: 
                raise 
                
            path = '/'.join([dir, fn])
            try: 
                S= load(path)
            except IOError as err: 
                path_not_found.append(path)
                if fault_tolerant: 
                    continue
                else: 
                    raise
            except Exception as e: 
                path_not_found.append((path, e))
                if fault_tolerant: 
                    continue 
                else: 
                    raise
            
            measure_S(S, dir, which, exclude_which=exclude_which, param=param, force=force, 
                    fault_tolerant=fault_tolerant, field_surfix=field_surfix,  **kwargs)
            
    if len(path_not_found)>0: 
        print 'path_not_found is %s'%path_not_found

def make_measure_many_args(dir_list, sh_list, sh_min=None, sh_max=None,  
        which=None, exclude_which=None, param=None,  force=0, field_surfix='', 
        fault_tolerant=1, recursive=False,  **kwargs): 
    use_local_storage = kwargs.get('use_local_storage', False)
    path_not_found = []
    __sh_list = sh_list
    if isinstance(which, str): 
        which = [which]
    if isinstance(dir_list, str): 
        dir_list = [dir_list]
    if recursive: 
        dir_list_1 = list(dir_list)
        dir_list = []
        for dir in dir_list_1: 
            temp=mera_backup_dir_finder(root=dir)
            dir_list.extend(temp)
        #print 'recursively found these dirs %s'%(dir_list, )
    aaa = os.path.abspath(dir_list[0])
    if 'mera' in aaa: 
        algorithm = 'mera'
    elif 'mps' in aaa and 'idmrg' not in aaa: 
        algorithm = 'mps'
    elif 'idmrg' in aaa:
        algorithm = 'idmrg'
    else: 
        raise
    
    path_list = []
    for dir in dir_list: 
        #fori d in dim_list: 
        if sh_list  == 'all': 
            db = ResultDB(dir, algorithm=algorithm)
            __sh_list = db.get_shape_list(sh_max=sh_max, sh_min=sh_min)
        
        for sh in __sh_list: 
            if algorithm  == 'mera':  
                #sh[1] -= 1  #layer->layer-1
                sh = list(sh);  sh[1] -= 1
                fn = System.backup_fn_gen(*sh)
            elif algorithm == 'mps': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            elif algorithm == 'idmrg': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            else: 
                raise 
                
            path = '/'.join([dir, fn])
            path_list.append(path)
    
    res= []
    for path in path_list:    
        if 0: 
            try: 
                S= load(path)
            except IOError as err: 
                path_not_found.append(path)
                if fault_tolerant: 
                    continue
                else: 
                    raise
            except Exception as e: 
                path_not_found.append((path, e))
                if fault_tolerant: 
                    continue 
                else: 
                    raise
        temp = dict(path=path, which=which, exclude_which=exclude_which, 
            fault_tolerant=fault_tolerant, field_surfix=field_surfix, 
            force=force, **kwargs)  
        res.append(temp) 
    
    return res 

def compute(n):
    import time, socket
    n = 3 
    time.sleep(n)
    host = socket.gethostname()
    return (host, n)

def measure_all_new(arg_group): 
    nodes= ['qtg7501']
    nodes= ['localhost']
    
    #cluster = dispy.JobCluster(measure_S,  nodes=nodes, depends=[])
    #cluster = dispy.JobCluster(measure_it,  nodes=nodes, depends=[])
    import merapy
    import vmps
    cluster = dispy.JobCluster(compute,  nodes=nodes, depends=[merapy, vmps])
    #for i in range(10):
    #    print  'iiiiiiiiiiiiiiiiiiiii', i 
    #    cluster.submit(i)
    job_list = []
    for i, a in enumerate(arg_group): 
        path = a['path']
        #job=cluster.submit(( ), **a)
        #job=cluster.submit((), a)
        job=cluster.submit(i) 
        job.id = i
        job_list.append(job)
    
    print len(job_list)
    for job in job_list:
        #host, n = job() # waits for job to finish and returns results
        #print '%s executed job %s at %s with %s' % (host, job.id, job.start_time, n)
        # other fields of 'job' that may be useful:
        # job.stdout, job.stderr, job.exception, job.ip_addr, job.end_time
        print job.stdout, job.stderr , job.id, job.exception  
    cluster.stats()

def measure_all_by_brokest(arg_group): 
    #for server_port in server_ports:
    #    Process(target=server, args=(server_port,)).start()
    from mypy.brokest.brokest import queue 
    Process( )
    for i in arg_group: 
        print  i['path']
        #print queue(measure_S, args=tuple(), kwargs=i, host='222.195.73.191', port=90900)
        #print queue(measure_S, args=tuple(), kwargs=i, host='210.45.121.30', port=90900)
        print queue(measure_S, args=tuple(), kwargs=i, host='127.0.0.1', port=90900)
        

def measure_all(dir_list, mera_shape_list, sh_min=None, sh_max=None,  
        which=None, exclude_which=None, param=None,  force=0, field_surfix='', 
        fault_tolerant=1, recursive=False,  **kwargs): 
    
    use_local_storage = kwargs.get('use_local_storage', False)
    path_not_found = []
    __mera_shape_list = mera_shape_list
    if isinstance(which, str): 
        which = [which]
    if isinstance(dir_list, str): 
        dir_list = [dir_list]
    if recursive: 
        dir_list_1 = list(dir_list)
        dir_list = []
        for dir in dir_list_1: 
            temp=mera_backup_dir_finder(root=dir)
            dir_list.extend(temp)
        #print 'recursively found these dirs %s'%(dir_list, )
    aaa = os.path.abspath(dir_list[0])
    if 'mera' in aaa: 
        algorithm = 'mera'
    elif 'mps' in aaa and 'idmrg' not in aaa: 
        algorithm = 'mps'
    elif 'idmrg' in aaa:
        algorithm = 'idmrg'
    else: 
        raise
    
    path_list = []
    for dir in dir_list: 
        #fori d in dim_list: 
        if mera_shape_list  == 'all': 
            db = ResultDB(dir, algorithm=algorithm)
            __mera_shape_list = db.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
        
        for sh in __mera_shape_list: 
            if algorithm  == 'mera':  
                #sh[1] -= 1  #layer->layer-1
                sh = list(sh);  sh[1] -= 1
                fn = System.backup_fn_gen(*sh)
            elif algorithm == 'mps': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            elif algorithm == 'idmrg': 
                N, D = sh
                fn = "N=%d-D=%d.pickle"%(N, D)
            else: 
                raise 
                
            path = '/'.join([dir, fn])
            path_list.append(path)
    
    for path in path_list:    
        if 1: 
            try: 
                #S= load(path)
                measure_S(path=path, which=which, exclude_which=exclude_which, param=param, force=force, 
                        fault_tolerant=fault_tolerant, field_surfix=field_surfix,  **kwargs)
            except IOError as err: 
                path_not_found.append(path)
                if fault_tolerant: 
                    continue
                else: 
                    raise
            except Exception as e: 
                path_not_found.append((path, e))
                if fault_tolerant: 
                    continue 
                else: 
                    raise
         
            
    if len(path_not_found)>0: 
        print 'path_not_found is %s'%path_not_found

def mera_backup_dir_finder(root): 
    res=[]
    for dir, subFolders, files in os.walk(root):
        for fn in files:
            if 'pickle' in fn:
                #print dir, files[:5]
                res.append(dir)
                break
    return res

class TestIt(unittest.TestCase): 
    def setUp(self): 
        args= {}
        args['which'] = ['correlation']
        args['force'] = 0
        args['show'] = 1
        if 0: 
            dir = '/home/zhli/Documents/mera_backup_tensor/run-wigner-crystal/a=2.0-h=1.4142135-V=0.0/'
            args['dir_list'] = [dir]
            args['db_version'] = None #1.0
            args['mera_shape_list'] = [(4, 4)]
        else: 
            dir = '/home/zhli/mps_backup_folder/run-ising/h=1.0/' 
            args['dir_list'] = [dir]
            args['mera_shape_list'] = [(40, 40)]
            args['fault_tolerant'] = 0
        self.args= args 
    
    def test_temp(self): 
        pass
        args= self.args
        dir = '/home/zhli/backup_tensor_dir/run-long-better/mera/new_run/alpha=1.9'
        args.update(dir_list=[dir], mera_shape_list=[(4, 4)], which=['entanglement_brute_force_6'])
        
        measure_all( **args )
   
    def test_idmrg_correlation(self): 
        args= self.args
        dir = '/home/zhli/mps_backup_folder/run-heisbg-long-spin1/alpha=2.0/'     
        dir = '/home/zhli/mps_backup_folder/run-ising/idmrg/h=1.0'
        args.update(dir_list=[dir], mera_shape_list=[(0, 8)], which=['correlation'])
        
        measure_all( **args )
   
    def test_idmrg_magnetization(self): 
        args= self.args
        dir = '/home/zhli/mps_backup_folder/run-heisbg-long-spin1/alpha=2.0/'     
        dir = '/home/zhli/mps_backup_folder/run-ising/idmrg/h=1.0'
        args.update(dir_list=[dir], mera_shape_list=[(0, 8)], which=['magnetization'])
        
        measure_all( **args )
    
    def test_vmps_all(self): 
        args = self.args
        dir = '/home/zhli/mps_backup_folder/run-heisbg-long-spin1/alpha=2.0/'     
        #which=['energy'],
        which = None
        args.update(dir_list=[dir],which=which, force=0, mera_shape_list=[(40, 10)])
        #measure_all( **args )
        measure_all( **args )
    
    def test_idmrg_all(self): 
        if 0: 
            args = self.args
            dir = '/home/zhli/mps_backup_folder/run-heisbg-long-spin1/alpha=2.0/'     
            which = None
            args.update(dir_list=[dir], which=which, force=0, mera_shape_list=[(40, 10)])
            #measure_all( **args )
            measure_all( **args )
        dir = tempfile.mkdtemp()
        from vmps.idmrg_mcc import iDMRG_mcc 
        i=iDMRG_mcc.example(model_name='ising', symmetry='Z2', backup_parpath=dir, 
                pinv_tol=1e-15, do_measure=0)
        i.minimize()
        path = '/'.join([dir, 'N=0-D=4.pickle'])
        measure_S(path=path, which=['correlation'], fault_tolerant=0)
         
    def test_make_measure_many_args(self): 
        dir = '/home/zhli/backup_tensor_dir/run-ising/vmps/main/h=0.0/'
        res = make_measure_many_args(dir, 'all', which=['energy'])
        for r in res: 
            measure_S(**r)
    
    def test_measure_all_new(self): 
        dir = '/home/zhli/backup_tensor_dir/run-ising/vmps/main/h=1.0/'
        res = make_measure_many_args(dir, 'all', 
                use_local_storage=1, 
                sh_max=(200, 20), which=['correlation'], force=1)
        if 0: 
            from merapy.utilities import save , load 
            save(res[0], '/tmp/mea_arg')
            save(measure_S, '/tmp/mea_func')
            
            mf = load('/tmp/mea_func')
            ma = load('/tmp/mea_arg')
            mf(**ma)
        
        if 0: 
            for r in res[: 2]: 
                measure_S(**r)
        if 1: 
            measure_all_new(res)
    
    def test_measure_all_by_brokest(self): 
        dir = '/home/zhli/backup_tensor_dir/run-ising/vmps/main/h=1.0/'
        args= make_measure_many_args(dir, 'all', 
                use_local_storage=1, 
                sh_max=(200, 20), which=['correlation'], force=1)
        
        measure_all_by_brokest(args[2: 10])
        #measure_S(**args[2])
        
    def tearDown(self): 
        pass
    


if __name__ == '__main__':
    
    warnings.filterwarnings('once')
    if 0:    # del  
        parser = argparse.ArgumentParser(description='dont know')
        parser.add_argument('dir_list', nargs='*', default=['./'])
        parser.add_argument('-w', '--which', nargs='*', default=None)
        parser.add_argument('-ew', '--exclude_which', nargs='*', default=None)
        parser.add_argument('-d', '--dim_list', nargs='*',  type=int,  default=None)
        parser.add_argument('-l', '--layer_list', nargs='*', type=int, default=None)
        parser.add_argument('-sh', '--mera_shape_list', nargs='*', default='all')
        parser.add_argument('-sh_min', '--sh_min', nargs=2, type=int,  default=None)
        parser.add_argument('-sh_max', '--sh_max', nargs=2, type=int, default=None)
        parser.add_argument('-f', '--force', action='store_true', default=False) 
        parser.add_argument('-ft', '--fault_tolerant', type=int, default=1) 
        parser.add_argument('-r', '--recursive', action='store_true', default=False) 
        #parser.add_argument('-f', '--force', action='store_true')
        args = parser.parse_args()
        args = vars(args)
        if not args['mera_shape_list'] in [None, 'all']: 
            dim_list = args['dim_list']; args.pop('dim_list')
            layer_list = args['layer_list']; args.pop('layer_list')
            mera_shape_list = [(d, l) for d in dim_list for j in layer_list]
            args['mera_shape_list'] = mera_shape_list
   
    if len(sys.argv)>1: 
    #if 0: 
        from cmd_line_args import parser
        parser.add_argument('-w', '--which', nargs='*', default=None)
        #some issue with mps algs, temporarily disable choices
        #parser.add_argument('-w', '--which', nargs='*', choices=FIELD_NAME_LIST,  default=None)
        
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
        args = vars(args)
        if not args['mera_shape_list'] in [None, 'all']: 
            dim_list = args['dim_list']; args.pop('dim_list')
            layer_list = args['layer_list']; args.pop('layer_list')
            mera_shape_list = [(d, l) for d in dim_list for j in layer_list]
            args['mera_shape_list'] = mera_shape_list
    
        measure_all( **args )
        
    else: 
        if 0: #examine
            #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
            #unittest.TextTestRunner(verbosity=0).run(suite)    
            TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
            unittest.main()
        else: 
            suite = unittest.TestSuite()
            add_list = [
                #'test_temp', 
                #'test_vmps_all', 
                'test_idmrg_all', 
                #'test_make_measure_many_args', 
                #'test_measure_all_new', 
                #'test_measure_all_by_brokest'
        
            ]
            for a in add_list: 
                suite.addTest(TestIt(a))
            unittest.TextTestRunner().run(suite)
           
    
    
    



