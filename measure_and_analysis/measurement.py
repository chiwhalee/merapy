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
try: 
    import dispy 
except: 
    warnings.warn('not able to import dispy')
from  multiprocessing import Process
import tempfile 

from merapy.common_util import set_num_of_threads 
#from merapy.hamiltonian import System
from merapy.utilities import load , dict_to_object 
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
#from merapy.config import BACKUP_BASE_DIR_LOCAL   #import error !

try: 
    import vmps.measure_and_analysis.measurement_idmrg as all_idmrg 
    import vmps.measure_and_analysis.measurement_idmrg_mcc as all_idmrg_mcc
    import vmps.measure_and_analysis.measurement_vmps as all_mps 
    #import vmps.measure_and_analysis.all as all_mps
except ImportError as err: 
    print err 

from merapy.measure_and_analysis.result_db import (ResultDB, ResultDB_vmps, ResultDB_mera, 
        ResultDB_idmrg, FIELD_NAME_LIST, BACKUP_STATE_DIR, RESULTDB_DIR, 
        algorithm_name_to_rdb, 
        )
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

def measure_S(S=None, parpath=None, path=None,
        measure_func=None, 
        which=None, 
        exclude_which=None, force=0, 
        fault_tolerant=1, num_of_threads=None,  param=None, check_convergence=True, 
        field_surfix='', rdb=None, **kwargs):
    """
        params: 
            measure_func: arbitrary function can be supported if it is provided 
                it must have one positional arg  S 
    """
    if measure_func is not None: 
        which = measure_func.__name__ 
    if num_of_threads is not None:  
        set_num_of_threads(num_of_threads)
    use_local_storage =  kwargs.get('use_local_storage', False)
    if S is None: 
        assert path is not None
        
        try: 
            S = rpyc_load(path, use_local_storage=use_local_storage)
            parpath = os.path.dirname(os.path.abspath(path))
            #issue: it can happen ~/backup_tensor_dir, ~/dropbox/resultdb_dir/ 
            #then it becomes ~/resultdb_dir so cause mismatch !
            parpath = parpath.replace(BACKUP_STATE_DIR, RESULTDB_DIR)
        except IOError as err: 
            if fault_tolerant: 
                return 
            else: 
                raise 
        except Exception as err:  # 比如 pickle file danmaged, exist but still cant load 
            if fault_tolerant: 
                err_str = str(err)
                warnings.warn('error'  + err_str[: 300]  + '......'  + err_str[-300: ])
                return 
            else: 
                raise 
    
    allow_measure = True
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
            elif S.has_key('mera'): 
                algorithm = 'mera'
                S= dict_to_object(S)
        else: 
            algorithm = 'mera'
        
        if algorithm == 'mera':  
            #field = list(FIELD_NAME_LIST)
            all_func = all_mera
           
            #if S.symmetry == 'U1': 
            #    try: 
            #        field.remove('entanglement_brute_force_9')
            #        field.remove('entanglement_brute_force_9_aver')
            #    except: 
            #        pass
            
            from merapy.hamiltonian import System 
            propt = System.key_property.im_func(S)
            #propt = S.key_property()
            dim,  layer = propt['trunc_dim'], propt['num_of_layer']
            nqn = len(propt['qns'])
            iter = S.iter
            if 1: 
                #fn = System.backup_fn_gen(dim=dim, layer=layer-1, nqn = len(propt['qns']))
                fn = backup_fn_gen(dim=dim, layer=layer-1, nqn = len(propt['qns']))

               
            path = parpath + '/' + fn
            
            #if db_version is None: 
            #    k = (dim, layer, iter)
            #    if not rdb.has_key(k): 
            #        raise NotImplemented   #后来修改了，不支持，需要再修改
            #        rdb[k] = OrderedDict()
            
            if nqn <= 3: 
                shape = (dim, layer)
            else: 
                shape = (dim, layer, nqn)
            
            rdb_class= ResultDB_mera 
            
            
        elif algorithm == 'mps': 
            
            #field = all_mps.DEFAULT_FIELD_LISt 
            all_func = all_mps
            #S = S['mps']  #fuck, this is bad, but I like bad!
            #mps= S 
            mps= S['mps']
            #db_version = 1.0
            N, D = mps.N, mps.D
            #if mps.symmetry is not None : 
            if hasattr(mps, 'symmetry'):  
                D = mps.bond_dim_max 
            fn = "N=%d-D=%d.pickle"%(N, D)
            path = '/'.join([parpath, fn])
            dim,  layer = None, None
            shape =  N,  D
            iter = -1
            
            
            rdb_class= ResultDB_vmps
            
        elif algorithm == 'idmrg': 
            #field = ['energy', 'correlation', 'correlation_length', 'magnetization', 
                    #'entanglement', 'entanglement_spectrum']
            if check_convergence: 
                c1 = S.get('converge_count')
                c2 = S.get('converge_count_lim')
                if c1 is not None  and c2 is not None: 
                    if not c1 >= c2:  
                        allow_measure = False 
                
            if S.has_key('lam_prev_inv'): 
                all_func = all_idmrg_mcc
            else: 
                all_func = all_idmrg
            A = S['A']
            symm = None if isinstance(A, np.ndarray) else A.symmetry 
            #db_version = 1.0
            N = 0
            #D = A.shape[0] if symm is None else A.shape[0].totDim 
            D = S.get('trunc_dim')
            fn = "N=%d-D=%d.pickle"%(N, D)
            path = '/'.join([parpath, fn])
            shape =  N, D 
            iter = -1
            
            rdb_class= ResultDB_idmrg 
            
    #which = which if which is not None else field 
    which = which if which is not None else []
    
    if isinstance(which, str): 
        which = [which]
    if exclude_which:
        print  'exlude this from measure_S: ', exclude_which
        which = [i for i  in which if i not in exclude_which]
    
    if measure_func is None: 
        func = {i: all_func.__getattribute__(i) for i in which}
    else: 
        func = {measure_func.__name__: measure_func}
    
    if param is None: 
        param = {}
    #elif isinstance(param, dict): 
    #    assert len(which)==1, 'when param is provided, only allow measure one field at a time'
    #else: 
    #    raise 
    
    func_name_map = {
            'calc_scaling_dim_2site': 'scaling_dim', 
            }
   
    #------------------------  start measureing ----------------------------------------#
    
    print  'meassuring %s at %s %s'%(path, shape, iter)
    if not allow_measure and not force: 
        print '\tstate is not converged. not allowing measure. return None'
        return 

    if rdb is None: 
        rdb = rdb_class(parpath, use_local_storage=use_local_storage)
    
    #if algorithm == 'mera' :  # backward compatibiliy for ResultDB version = None
    db_version = rdb.get('version')
    if not force: 
        found = []
        for w in which: 
            if db_version is None: 
                is_found = rdb[k].has_key(w)
            else:
                #key_list = [w] + [shape, iter]
                key_list = [w] + [shape]
                is_found = rdb.has_key_list(key_list)
            if is_found: 
                found.append(w)
        if found:
            print 'found these fields {}, omit them.'.format(found)
        which = [w for w in which if w not in found]
     
    failed_list = []
    results= []
    for w in which: 
        ff = func[w]
        w = func_name_map.get(ff.__name__, ff.__name__)
        field_name = w if field_surfix  == ''  else w + '_' +  field_surfix  
        print '\t%s...   '%field_name,
        
        try:
            _param = param.get(w, {})
            res = ff(S, **_param)
            msg = '%20s'%('done')
            if kwargs.get('show', False): 
                pprint.pprint(res, indent=1, width=200)
            
            if db_version is None: 
                results.append(dict(field_name=field_name, key=k, res=res, sub_key_list=None))
            else: 
                #results.append(dict(field_name=field_name, sh=shape, iter=iter, val=res))
                results.append(dict(field_name=field_name, sh=shape, val=res))
                
        except Exception as err:
            failed_list.append((w, shape, iter, err))
            msg = '%20s'%('FAILED skip')
            if not fault_tolerant: 
                raise 
        print msg
    
    if 1: 
        #load rdb again! so that avoid chance of confliction 
        rdb = rdb_class(parpath, use_local_storage=use_local_storage) 
        changed = False
        if len(results)>0:
            for r in results: 
                if db_version is None: 
                    rdb.put(**r)
                else: 
                    rdb.insert(**r)
                
            changed = True
        
        if 1: 
            if not rdb.has_key('algorithm'): 
                rdb['algorithm'] = algorithm
                changed = True
                
            if algorithm in ['mps', 'idmrg']: 
                if not rdb.has_key('dim_max'): 
                    rdb['dim_max'] = {}
                if rdb.has_key_list(['dim_max', shape[0]], info=0): 
                    dmax = rdb.get_dim_max_for_N(shape[0])   #note this need access to local file system
                    if dmax<shape[1]: 
                        rdb['dim_max'][shape[0]] = shape[1]
                        changed = True
                else: 
                    dmax = shape[1]
                    rdb['dim_max'][shape[0]] = dmax 
                    changed = True
            
        if changed: 
            rdb.commit(info=1)
        else: 
            print 'No change to commit.'
    if len(failed_list)>0: 
        print 'failed_list: %s'%failed_list

def backup_fn_gen(dim, layer, nqn=None):
    nqn = "5_" if nqn == 5 else "" 
    layer = "-%dlay"%layer if layer != 3 else ""
    res= "%(nqn)s%(dim)d%(layer)s.pickle"%vars()
    return res

def make_measure_many_args(dir_list, sh_list, sh_min=None, sh_max=None,  
        which=None, exclude_which=None,  force=0, field_surfix='', 
        fault_tolerant=1, algorithm=None, recursive=False,  **kwargs): 
    use_local_storage = kwargs.get('use_local_storage', False)
    path_not_found = []
        
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
    if algorithm is None: 
        if 'mera' in aaa: 
            algorithm = 'mera'
        elif 'mps' in aaa and 'idmrg' not in aaa: 
            algorithm = 'mps'
            rdb_class= ResultDB_vmps
        elif 'idmrg' in aaa:
            algorithm = 'idmrg'
            rdb_class= ResultDB_idmrg 
        else: 
            raise
    rdb_class= algorithm_name_to_rdb(algorithm)
        
   
    path_list = []
    for dir in dir_list: 
        #db = rdb_class(dir, algorithm=algorithm)
        db = rdb_class(dir.replace(BACKUP_STATE_DIR, RESULTDB_DIR), algorithm=algorithm)
        if sh_list  == 'all': 
            ss = db.get_shape_list(sh_max=sh_max, sh_min=sh_min)
        else:
            vv = set(db.get_shape_list())
            ss = list(set(sh_list).intersection(vv))
        for sh in ss: 
            if algorithm  == 'mera':  
                #sh[1] -= 1  #layer->layer-1
                sh = list(sh);  sh[1] -= 1
                #fn = System.backup_fn_gen(*sh)
                fn = backup_fn_gen(*sh)
            elif algorithm in ['mps', 'idmrg']: 
                N, D = sh
                if D == 'max' : 
                    D = db.get_dim_max_for_N(N)
                fn = "N=%d-D=%d.pickle"%(N, D)
            else: 
                raise 
                
            path = '/'.join([dir, fn])
            path_list.append(path)
    
    res= []
    for path in path_list:    
        temp = dict(path=path, which=which, exclude_which=exclude_which, 
            fault_tolerant=fault_tolerant, field_surfix=field_surfix, 
            force=force, **kwargs)  
        res.append(temp) 
    
    return res 

#now I undertand that, this func in fact realized a decorator for individal measure funcs!!
def measure_all(dir_list, mera_shape_list, use_dist_comp=False, submit=1, delay_send=10, servers=None, **kwargs): 
    arg_group = make_measure_many_args(dir_list, mera_shape_list, **kwargs)
    if not use_dist_comp: 
        for i in arg_group: 
            measure_S(**i)
    else: 
        from mypy.brokest.brokest import queue, run_many 
        from mypy.brokest.task_center import submit_one 
        if not submit:     
            tasks = [(measure_S, (), i) for i in arg_group] 
            run_many(tasks, servers, info=kwargs.get('info', 1), 
                    try_period=kwargs.get('try_period', 1))
        else: 
            for c in arg_group: 
                c.update(delay_send=delay_send)
                submit_one(measure_S, kwargs=c)

def compute(n):
    import time, socket
    n = 3 
    time.sleep(n)
    host = socket.gethostname()
    return (host, n)

def measure_all_by_dispy(arg_group): 
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
        args= {
                'which':['correlation'], 
                'force':0, 
                'show':1, 
                'fault_tolerant':0, 
                }
        if 1:
            dir = '/home/zhli/mps_backup_folder/run-ising/h=1.0/' 
            args['dir_list'] = [dir]
            args['mera_shape_list'] = [(40, 40)]
        self.args= args 
    
    def test_measure_all(self): 
        from vmps.minimize import VMPS 
        dir = tempfile.mkdtemp()
        v = VMPS.example(h=0.0, N=16, Jzz=1.0, model_name='heisbg', D=None, 
                symmetry='U1', num_of_sweep_tot=6, backup_parpath = dir, 
                do_measure=1, measure_only=['energy']) 
        
        v.minimize()
        print '\n'*2
        args = {
                'which':['correlation', ],  #['correlation', 'variance', 'energy', 'entanglement_entropy', 'magnetization']
                'param':{'correlation':{'i0':4}}, 
                'force':1, 
                'show':1, 
                'submit':0, 
                'dir_list':[dir], 
                'algorithm':'mps', 
                #'mera_shape_list':'all',
                'mera_shape_list':'all',
                'fault_tolerant':0, 
                }
        measure_all( **args )
   
    def test_vmps_all(self): 
        from vmps.minimize import VMPS 
        v = VMPS.example(N=20, D=5, model_name='ising', backup_parpath='temp')
        v.minimize()
        dir = v.backup_parpath
        args = self.args.copy()
        which = ['variance', 'energy', 'entanglement_entropy', 'magnetization']
        args.update(dir_list=[dir],which=which, force=0, show=0, 
                algorithm='mps', mera_shape_list='all')
        #measure_all( **args )
        measure_all( **args )
        db = ResultDB_vmps(dir)
        print_vars(vars(),  ['db'])
    
    def test_idmrg_all(self): 
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
    
    def test_measure_all_by_dispy(self): 
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
            measure_all_by_dispy(res)
    
        
    def tearDown(self): 
        pass
    
    def test_temp(self): 
        raise
        #from projects_mps.run_long_sandvik.analysis import an_vmps,  an_idmrg_psi 
        #xx = an_vmps.an_main_symm 
        #
        #aa=xx.filter_alpha(g=0.0)
        ##aa = [(0.3,  0.0)]
        #xx.measure_all(aa, [(60, 'max')], which='entanglement_entropy_ln', fault_tolerant=1)

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
                'test_measure_all', 
                #'test_vmps_all', 
                #'test_idmrg_all', 
                #'test_make_measure_many_args', 
                #'test_measure_all_by_dispy', 
        
            ]
            for a in add_list: 
                suite.addTest(TestIt(a))
            unittest.TextTestRunner().run(suite)
           
    
    
    



