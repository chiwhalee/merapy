#!/usr/bin/env python
#coding=UTF8

import warnings 
import time
import datetime 
import pprint
import pickle
import os
import sys
import socket
from collections import OrderedDict
import argparse
import unittest
import random
import multiprocessing, threading
from multiprocessing import Pool
import traceback 
import psutil
import pdb
import itertools
import signal 
import tempfile 
from tempfile import mkdtemp
import platform 

uname=platform.uname()[1]
if uname in ['QTG-WS1-ubuntu']: 
    import pp


from merapy.utilities import print_vars, save 
from merapy.context_util import (make_temp_dir, rpyc_conn_local, 
        rpyc_conn_local_zerodeploy, rpyc_conn, rpyc_conn_zerodeploy)
from merapy.quantum_number import *

from merapy.hamiltonian import System
from merapy.tensor_py import iTensor, iTensorFactory
from merapy.scale_invariant import ScaleInvar
from merapy import common_util
from merapy import crandom
import merapy.schedule
from merapy.schedule import copy_schedule, schedule_scale_invar, schedule_prod_state
from merapy.decorators import * #tensor_player
import merapy.main 


"""
    methodology: 
        WHAT IS THE ESSENSE OF A MAIN FUNC OR CLASS?
            1. main func is the interface that connect the program and the hardware 
                e.g. determing the runtime resorces for parallel runing 
            2. it explicitly indicates where the program starts
            3. parallel computing and distribute computing
    q:
        gamma=1.0 or 0.2
        in hamiltonian.py
            q387,  
        main.py
            q895
        Finisite.f90
            a.call TopLevel_EigenState(M,S)
            b. ascending,  descending
    issue:
        in mera.py:
            following is not complete
            #QSp_V[3] = QSp_U[0].copy()
            QSp_V[3]=QSp_V[0].add(QSp_V[1]).add(QSp_V[2])
            QSp_V[3].update(qsp_max = qsp_max)
        

    generally, main func has 3 effects: 
        1. init all--- pass all params
        2. excute; 
        3. collect and manage result
        Q:
            should I use keyword args or kargs?
        
    
    parallel wraper for running programs
    
    one main problem involved is picklability.  
        commone problems: 
            under pickle under  __main__ module
                it may fail when unpickle in another module, resulting AttributeError 
                    
                the reason is following: 
                Remember that pickle doesn't actually store information about how a
                class/object is constructed, and needs access to the class when
                unpickling. 
                the reason and workaround see 
                    https://skycoders.wordpress.com/2014/04/19/python-pickle-cannot-be-used-properly-in-__main__/
                    http://stackoverflow.com/questions/27732354/unable-to-load-files-using-pickle-and-multipile-modules
    
    
"""


if 1: 
    import copy_reg
    import types
 
    # to make  multiprocessing compatable with pickleing of method funcs
        # This call to copy_reg.pickle allows you to pass methods as the first arg to
        # mp.Pool methods. If you comment out this line, `pool.map(self.foo, ...)` results in
        # PicklingError: Can't pickle <type 'instancemethod'>: attribute lookup
        # __builtin__.instancemethod failed
    
    #the following correctly handling classmethods
    def _pickle_method(method):
        # Author: Steven Bethard
        # http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
        func_name = method.im_func.__name__
        obj = method.im_self
        cls = method.im_class
        if func_name.startswith('__') and not func_name.endswith('__'):
            #deal with mangled names
            cls_name = cls.__name__.lstrip('_')
            func_name = '_%s%s' % (cls_name, func_name)
        return _unpickle_method, (func_name, obj, cls)

    def _unpickle_method(func_name, obj, cls):
        # Author: Steven Bethard
        # http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
        if obj and func_name in obj.__dict__:
            cls, obj = obj, None # if func_name is classmethod
        for cls in cls.__mro__:
            try:
                func = cls.__dict__[func_name]
            except KeyError:
                pass
            else:
                break
        return func.__get__(obj, cls)

    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method) 
    
class Main(object): 
    LOCALHOSTNAME = 'QTG-WS1-ubuntu'  #used in distributed computing 
    def __init__(self, **config): 
        for k, v in config.iteritems(): 
            setattr(self, k, v)
        self.config = config
        self.initialized = False 
    
    def init_alg(self): 
        config = self.config
        #print_vars(vars(),  ['USE_CUSTOM_RAND'], head='\n\nCCCCCCCC\n\n')
        if config.get('USE_CUSTOM_RAND'):
            crandom.mrandom.csrand(config['rand_seed'])
        else:
            random.seed()
            crandom.rand = random.random
        
        args = self.config.copy()
        args['info'] = self.info - 1
        #self.show_config()
        Main.set_num_of_threads(self.NUM_OF_THREADS, config.get('info', 1))
        
        system_class = config.get('system_class', System)
        #self.S = system_class(None, ham, **args )
        args = config.copy()
        self.S = system_class(args['SYMMETRY'], None, args.get('MODEL'), **args )
    
    def run(self, **kwargs): 
        if not self.initialized: 
            self.init_alg()
            self.initialized = True 
        
        #Main.set_num_of_threads(self.NUM_OF_THREADS)
        self.S._minimize_finite_size(**kwargs)
        return self.S

    def run_scale_invar(self, **kwargs): 
        if not self.initialized: 
            self.init_alg()
            self.initialized = True 
        
        #Main.set_num_of_threads(self.NUM_OF_THREADS)
        self.S._minimize_scale_invar(**kwargs)
        return self.S

    def run_schedule(self, schedule=None, mera_shape_min=(4, 4), mera_shape_max=(12, 4), 
            q_iter_max=None, backup_parpath=None,  filename='auto', backup_fn='auto', 
           do_measure=None, skip_list=[],  skip_dim_list=[],  dim_diff_remap=None, info=0):
        """
            schedule is a OrderedDict like
            (ansatz, dim, layer, q_iter, run_time_max, energy_diff_min)
            
            { 4:  ('ScaleInvar', '')}
            
            }
        """
        #super(self.__class__, self).run_schedule(**kwargs)
        if not self.initialized: 
            self.init_alg()
            self.initialized = True 
            self.S.initialize()
        
        if schedule is None: 
            schedule = schedule_scale_invar
        schedule = copy_schedule(schedule, skip_list=skip_list)
        do_measure  = do_measure if do_measure is not None else self.do_measure
        
        #if len(dim_diff_remap)>0: 
        if dim_diff_remap is not None:  
            #print dim_diff_remap; exit()
            print 'remap of dim with energy_diff_min is %s'%dim_diff_remap
            for i in range(len(schedule)): 
                dim = schedule[i]['shape'][0]
                for d in sorted(dim_diff_remap.keys()): 
                    #print 'ddddd', i, d, dim
                    if dim >= d: 
                        diff  = dim_diff_remap[d]
                        schedule[i]['energy_diff_min'] = diff
            
        symmetry = self.S.symmetry
        #symmetry = self.config['SYMMETRY']
        qsp_class=symmetry_to_QspClass(symmetry)
        
        #for key, val in schedule.items():
        for task in schedule: 
            which_minimize = task['ansatz']
            if task['ansatz'] == 'prod_state':  
                run_func = System._minimize_finite_size 
                run_param_name = ('q_iter', 'energy_diff_min')
            elif task['ansatz'] == 'scale_invar': 
                run_func = System._minimize_scale_invar
                run_param_name = ('q_iter', 'run_time_max', 'energy_diff_min')
            elif task['ansatz'] == 'eigenstate': 
                raise NotImplemented
            else: 
                raise
            
            if q_iter_max is not None : 
                if task['q_iter'] > q_iter_max: 
                    task['q_iter']  = q_iter_max
            shape = task['shape']
            
            
            dim, layer = shape[:2]
            nqn = None
            if len(shape)>2: nqn = shape[2]
            
            msg = '\nSTARTING TASK %s'%(task, )
            if info>0: print msg
            
            if dim in skip_dim_list: 
                msg='dim=%d is in skip_dim_list, skip it'%(dim, )
                if info>0: print msg
                continue
            try: 
                qsp_destination = qsp_class.max(dim, nqn)
            except ValueError as err:
                msg = 'error occurred, skip %s'%(shape, )
                if info>0: print err, '\n',  msg
                continue
            
            if 1: 
                if mera_shape_min is not None : 
                    dim_min, layer_min = mera_shape_min[: 2]
                    nqn_min = None
                    if len(mera_shape_min)>2: nqn_min = mera_shape_min[2]
                    qsp_min = qsp_class.max(dim_min, nqn_min)
                    
                    
                    #if qsp_destination < qsp_min or layer <layer_min: 
                    if (not  qsp_destination >= qsp_min)   or layer <layer_min: 
                        #msg = 'shape (%d, %d) smaller than mera_shape_min %s, skip it'%(dim, layer, mera_shape_min, )
                        msg = 'shape %s smaller than mera_shape_min %s, skip it'%(shape, mera_shape_min, )
                        if info>0: print msg
                        continue
                
                if mera_shape_max is not None :    
                    dim_max, layer_max = mera_shape_max[: 2]
                    nqn_max = 3
                    if len(mera_shape_max)>2: nqn_max = mera_shape_max[2]
                    qsp_max = qsp_class.max(dim_max, nqn_max)
                    if 0: 
                        print qsp_destination
                        print qsp_max
                        print  'qsp_destination>qsp_max', qsp_destination>qsp_max
                        exit()
                    
                    if qsp_destination>qsp_max or (nqn> nqn_max and self.SYMMETRY=='U1')  or layer> layer_max: 
                    #if qsp_destination>qsp_max or nqn> nqn_max or layer> layer_max: 
                    #if qsp_destination>qsp_max  or layer> layer_max: 
                        
                        msg = 'shape %s exceeds mera_shape_max %s, skip it'%(shape, mera_shape_max, )
                        if info>0: 
                            print msg
                        continue

                propt = self.S.key_property()
                dim_current = propt['trunc_dim']
                layer_current = propt['num_of_layer']
                nqn_current = propt['nqn']
                qsp_current = qsp_class.max(dim_current, nqn_current)
                shape_current = (dim_current, layer_current, nqn_current)
                
                
                
                if qsp_current < qsp_destination: 
                    self.S.expand_dim(dim, nqn)
                
                if layer_current <layer: 
                    self.S.expand_layer(layer)
                
                #if dim_current>dim or layer_current> layer: 
                if qsp_current> qsp_destination or layer_current> layer: 
                    msg = 'shape %s smaller than current shape = %s, skip it'%(shape, shape_current)
                    print msg
                    continue
                   
            run_param = dict([(i, task.get(i)) for i in run_param_name])
            
            #run_func(self, backup_parpath=backup_parpath, filename=filename, backup_fn=backup_fn, 
            #        do_measure=do_measure, **run_param)
            
            #run_func(self.S, backup_parpath=backup_parpath, filename=filename, backup_fn=backup_fn, do_measure=do_measure, **run_param)
            self.S.minimize(which_minimize, **run_param)

    def show_config(self): 
        return 
        aaa = [ 'model_param']
        aaa = {a: self.__getattribute__(a) for a in aaa}

        msg = ''
        msg += '\n\n' 
        msg += str( datetime.datetime.today()) + '\n'
        msg += str(aaa)
        if self.filename is not None : 
            out = open(self.filename, 'a')
            out.write(msg)
            out.close()
        print msg
            
    @classmethod
    def run_one(cls, config): 
        schedule = config.get('schedule')
        pid = os.getpid();  config.update(pid=pid)
        #schedule = config['run_schedule_args']['schedule']
        main = cls(**config)
        msg = ''
        is_registered = False 
        job_status = ''
        try: 
            status = 'open'
            if config.get('register_job'): 
                job_id = config.get('job_id')
                if job_id is None: 
                    parpath = config['backup_parpath_local']
                    job_id = (str(config.get('job_group_name')), hash(parpath))
                with rpyc_conn('local', 'service', 17012) as conn:                
                    status = conn.root.register_job(job_id, config)
                if status  == 'open' : 
                    is_registered = True
            if status == 'open': 
                if schedule is None: 
                    main.run(q_iter=None, backup_parpath=None, resume=0, auto_resume=0 )
                else: 
                    #print_vars(vars(),  ['schedule']); raise 
                    #if isinstance(schedule[0], dict):  #mera use this 
                    if config['algorithm']== 'mera':  
                        main.run_schedule(**config['schedule'])
                    else:   #mps use this 
                        main.run_schedule(schedule)  
                job_status = 'COMPLETED'
                msg += 'run one COMPLETED.' 
            elif status== 'locked': 
                msg += 'job %s exists'%(job_id, ) 
        except KeyboardInterrupt: 
            #print 'KeyboardInterrupt in child for %s' %(config.get('model_param'), )
            msg += 'KeyboardInterrupt captured.'
            job_status = 'FAILED'
        except Exception, exception:
            print exception
            traceback.print_exc()
            job_status = 'FAILED'
            raise
        finally: 
            msg += ' EXIT chiled for %s.\n final status is %s'%(config.get('model_param'), job_status)
            print msg
            if is_registered:  
                parpath = config['backup_parpath_local']
                with rpyc_conn('local', 'service', 17012) as conn:                
                    print 'iiiiiiiiiiiiiii', job_id 
                    conn.root.update_job(job_id, 'status', job_status)
                    #end_time = str(datetime.datetime.now()).split('.')[0]
                    #conn.root.update_job(job_id, 'end_time', end_time)
               
            main.stop_player()        
        return main 
    
    @classmethod #@staticmethod 
    def run_many(cls, config_group, nproc=None, func=None, parallel=True): 
        """
            How to about multiprocessing: 
                map is used for one unique function, while apply for different functions
                the former is simple,  while the later is more flexible
        """
        res= []
        #init_worker = lambda: signal.signal(signal.SIGINT, signal.SIG_IGN)
        if func is None: 
            func = cls.run_one
        
        if parallel:
            #nproc = nproc if nproc is not None else len(config_group)
            if nproc is None:   #setting default values of nproc of each machine. this assume NUM_OF_THREADS is uniform for each process
                hostname = socket.gethostname()
                nt = config_group[0]['NUM_OF_THREADS']
                ncpu_tot = psutil.NUM_CPUS 
                if 'node' in hostname:
                    ncpu_tot = min(ncpu_tot, 12)
                    nproc = ncpu_tot//nt 
                elif 'qtg' in hostname: 
                    nproc = ncpu_tot/2//nt
                elif hostname == Main.LOCALHOSTNAME: 
                    #nproc = 2
                    nproc = ncpu_tot/2//nt
                else: 
                    warnings.warn('unfamiliar hostname %s'%(hostname, ))
                    nproc = ncpu_tot/2//nt 
                #assert nproc>0, 'nproc=%d'%(nproc, )
                
                nproc = min(nproc, 12)
                print 'nproc is auto set to %d'%(nproc, )
                
            try: 
                pool = Pool(processes=nproc, initializer=None)
                if 0: 
                    pool.map(func, config_group)
                else: 
                    for arg in config_group: 
                        #print 'aaaaa', arg 
                        pool.apply_async(func,  (arg, ))
                    pool.close()   
                    pool.join()   # the above two lines are needed to prevent main thread exit
            except KeyboardInterrupt: 
                print 'KeyboardInterrupt in main exit'
            except Exception:
                print 'EEEEE'*10
                sys.exit(0)
            finally: 
                #pool.terminate()
                #pool.join()
                #sys.exit(0)
                pass
                
        else: 
            for i in config_group: 
                res.append(func(i))
        return res 
    
    @classmethod #@staticmethod 
    def run_many_dist(cls, config_group, 
            submit=True, job_info=None,  
            servers=None, func=None, need_confirm=True,  **kwargs): 
        """
            status: workable now!
                I dont know why cant run run_many_dist under main.py.
                This will cause problem in serialization  
                so Iam not able to unittest it under __main__ 
                
        """
        #if func is None: 
        #    #func = cls.run_one
        #    func = run_one_dist 
        from brokest.brokest import queue, run_many 
        from brokest.task_center import submit_one, submit_many, LOCAL_IP
        if need_confirm:
            i=raw_input('submit jobs? yes(y)\n')
            if i.lower()=='y':
                print 'allow'
            else:
                print  'canceled'
                return 
        
        if not submit: 
            tasks = [(cls.run_one, (c, )) for c in config_group] 
            run_many(tasks, servers, 
                    info=kwargs.get('info', 10), querry=1, 
                    try_period=kwargs.get('try_period', 1)
                    )
        else: 
            #job_info = {}
            #temp = ['delay_send', 'priority', 'job_group_name', 'job_description']
            #for t in temp: 
            #    if kwargs.has_key(t): 
            #        job_info[t] = kwargs[t]
            #raise  #modify job_info !!
            for c in config_group: 
                status=submit_one(cls.run_one, (c, ),  job_info=job_info, querry_timeout=kwargs.get('querry_timeout', 10))
                print status
            
            #if 1:
            #    print 'updationg job order according to their priority... '
            #    tc = TaskCenter(host=LOCAL_IP)
            #    print tc.update_job_order()
            
    @staticmethod 
    def server_discover(servers=None, exclude_patterns=None,  querry_timeout=1, qsize=8, info=0): 
        from brokest.brokest import server_discover 
        return server_discover(servers=servers, exclude_patterns=exclude_patterns, querry_timeout=querry_timeout, qsize=qsize, info=info)
    
    @staticmethod
    def stop_servers(servers): 
        from brokest.brokest import stop_servers 
        stop_servers(servers)
    
    @staticmethod
    def make_grid(*args): 
        grid = itertools.product(*args)
        return list(grid)  #efficiency doesnt matter,  list is better
    grid = make_grid  # another name
    
    @staticmethod
    def grid_merger(*args): 
        return itertools.chain(*args)

    def thread_compare(self, func, args, thread_list):
        """  this idea is stupid and not work!!"""
        raise NotImplemented  
        res= []
        for i in thread_list:
            NUM_OF_THREADS = i
            os.environ["OMP_NUM_THREADS"] = str(NUM_OF_THREADS)
            t1 = time.time()
            ta = time.clock()
            func(**args)
            tb = time.clock()
            t2 = time.time()
            temp = (i, tb-ta, t2-t1)
            
            res.append(temp)
        pprint.pprint(res) 
    
    @staticmethod
    def set_num_of_threads(n, info=1):
        #import common_64_ifort as c64
        common_util.set_num_of_threads(n)
        if info>0: 
            print 'set_num_of_threads to %d'%n
    
    @staticmethod
    def get_num_of_threads(): 
        #import common_64_ifort as c64
        raise NotImplemented('not pass test')
        n=common_util.get_num_of_threads()
        print 'get_num_of_threads %d'%n
        return n
    
    @staticmethod
    def get_cpu_usage(): 
        p = psutil.Process(os.getpid())
        #p.get_cpu_times()
        #cputimes(user=0.06, system=0.03)
        res=p.get_cpu_percent(interval=0.01)
        return res
        
    def transfer_pickle_file(self, fn=None): 
        """
            when complete, transfer file back to QTG-WS1-ubuntu
            result is always supposed transfer back to center after computation
            
            to make it work, 强制要求 远程和local 的 backup tensor目录结构必须完全一致！！！ 这牺牲了一点灵活性，但是值得的
            具体地remote and local 的目录都必须有 mera_backup_tensor 这一层; 在这一层的左侧，目录依赖于系统，
            而在右侧，必须相同
            
        """
        if self.hostname != self.LOCALHOSTNAME: 
            if fn is None: 
                fn = self.S.backup_fn_auto()
            remote_path = '/'.join([self.backup_parpath, fn])
            local_path = '/'.join([self.backup_parpath_local, fn]) 
            temp = self.config['LOCAL_USERNAME']  + '@'  + self.config['LOCAL_IP']  + ":"
            #local_path = 'zhli@210.45.117.30:' + local_path
            local_path = temp + local_path
            cmd = 'scp %s %s'%(remote_path, local_path)
            msg = 'try to transfer pickle file back to local:  '
            msg += cmd  
            status = os.system(cmd)
            if status == 0: 
                msg += '---succeed'
            else: 
                msg += '---failed'
            print msg
            with open(self.filename, 'a') as inn: 
                inn.write(msg)

    def transfer_result(self): 
        host_current = socket.gethostname() 
        if host_current == 'QTG-WS1-ubuntu': 
            return
        else: 
            host1 = socket.gethostname()  
            host2 = 'QTG-WS1-ubuntu'

        client = scp.Client(host=host, user=user, keyfile=keyfile)
        
        #option 2
        client = scp.Client(host=host, user=user)
        client.use_system_keys()
        
        #option 3
        client = scp.Client(host=host, user=user, password=password)

        # and then
        client.transfer('/etc/local/filename', '/etc/remote/filename') 
    
    def stop_player(self):
        set_STATE_end_1(iter=None, record_at=None, stop_at=None, power_on=False)

    def __del__(self, ): 
        
        if self.config.get('info', 1)>-1: 
            print 'exit main, and run main.stop_player'
        self.stop_player()


class TestMain(unittest.TestCase): 
    def setUp(self): 
        if 1: 
            from config import CFG_HEISBG_BASIC, updaters_u1
            temp = dict(USE_CUSTOM_RAND=True, updaters=updaters_u1, trunc_dim=4, tot_layer=4, use_player=True, 
                    SYMMETRY="Travial", 
                    NUM_OF_THREADS=1, do_measure=0)
            cfg_heisbg = CFG_HEISBG_BASIC.copy();   cfg_heisbg.update(temp)
        if 1: 
            from merapy.config import CFG_ISING_BASIC, copy_config, updaters_u1 
            from merapy.finite_site import finite_site_u1
            cfg_ising = copy_config(CFG_ISING_BASIC) 
            cfg_ising['updaters'] = updaters_u1
            
        #model_param={"J_NN":1.0, "J_NNN":0.241186}
        #cfg_heisbg['model_param'].update(model_param)
        #if 1:  #test run_schedule
        try: 
            os.system('rm ./mera_backup_test_folder/*')    
            os.system('rm 4.pickle')
        except: 
            pass 
        
        self.config_heisbg = cfg_heisbg

    def test_grid(self): 
        a = [1, 2]
        b = ['a', 'b']
        g1 = Main.grid(a, b)
        g2 = Main.grid(['x'], [0])
        gg = Main.grid_merger(g1, g2)
        for g in gg: 
            print g 
        
    def test_run(self): 
        config = self.config_heisbg 
        config['SYMMETRY'] = 'Travial'
        config.update(q_iter=10)
        main = Main(**config)
        main.run(q_iter=10, backup_parpath='./mera_backup_test_folder', do_measure=0)
    
    def xtest_run_on_remote(self): 
        #with rpyc_conn('sugon') as conn: 
        with rpyc_conn_zerodeploy('sugon') as conn: 
            conn.namespace['self'] = self
            code = r"""from merapy.main import TestMain; self.test_run()"""
            #Main_remote = conn.modules['merapy.main'].Main
            conn.execute(code)
            print conn.modules.os.getcwd()
    
    def xtest_resume_func(self): 
        config = self.config_heisbg 
        config['SYMMETRY'] = 'Travial'
        with make_temp_dir() as d1, make_temp_dir() as d2: 
            config['backup_parpath'] = d1
            config['backup_parpath_local'] = d1
            main = Main(**config)
            main.run(q_iter=5)
            main.run_scale_invar(q_iter=5)
        
            print '\n'*3
            with rpyc_conn_zerodeploy('sugon') as conn: 
                main.use_local_storage = 1
                conn.namespace['main'] = main
                code = r"""from merapy.main import TestMain; self.test_run()"""
                code = 'main.run_scale_invar(q_iter=5, resume=1)'
                #code = 'main.resume_func()'
                #Main_remote = conn.modules['merapy.main'].Main
                conn.execute(code)
                print conn.modules.os.getcwd()

    def test_resume(self): 
        config = self.config_heisbg.copy()
        config['SYMMETRY'] = 'Travial'
        config.update(USE_CUSTOM_RAND=1)
        dir = tempfile.mkdtemp()
        config['backup_parpath'] = dir 
        main = Main(**config)
        main.run(q_iter=5)
        
        self.assertTrue(os.path.exists(dir + '/4.pickle'))
        
        if 1:  
            main = Main(**config)
            main.run(q_iter=10)
        self.assertAlmostEqual(main.S.energy, -1.15128521352, 10) 

    def xtest_make_dir(self): 
        config = self.config_heisbg 
        print config['backup_parpath']
        print config['backup_parpath_local']
        a ='/tmp/lksjdfoi/oiweffw' 
        b ='/tmp/lksjdfoi/asdfoi' 
        main = Main(**config)
        main.make_dir(a, b)
        os.rmdir(a)
        with rpyc_conn_local() as conn: 
            os_local = conn.modules['os']
            print os_local
            exist = os_local.path.exists(b)
            os_local.rmdir(b)
            
        self.assertTrue(exist)
        
    def test_run_scale_invar(self): 
        config = self.config_heisbg
        #config['SYMMETRY'] = 'U1'
        config['SYMMETRY'] = 'Travial'
        #config['measurement_args']['exclude_which'].append('correlation_extra')  # too slow
        main = Main(**config)
        with make_temp_dir() as dd: 
            main.run_scale_invar(q_iter=5, backup_parpath=dd, do_measure=1)
        
    def test_run_schedule(self): 

        config = self.config_heisbg
        config['SYMMETRY'] = 'Travial'
        config['schedule']['schedule'] = copy_schedule(schedule_prod_state)
        dir = tempfile.mkdtemp()
        config.update(tot_layer=3, backup_parpath=dir )
        main = Main(**config)
        main.run_schedule(schedule=schedule_scale_invar, do_measure=0, #backup_fn='auto', 
                q_iter_max = 5, mera_shape_min = (4, 3), mera_shape_max=(8, 4), skip_dim_list=[14])
        S = main.S
        self.assertEqual(S.iter, 15)
        self.assertEqual(S.iter1, 5)


    def xtest_set_num_of_threads_1(self): 
        config = self.config_heisbg
        #config['NUM_OF_THREADS'] = 4
        main = Main(**config)
        main.run(q_iter=5)
        print Main.get_cpu_usage()
    
    def xtest_run_many(self): 
        """  
            not completed
        """
        from merapy.config import CFG_ISING_BASIC, copy_config, updaters_u1
        NNN = 1
        #config_list = [config_copy(CFG_ISING_BASIC) for i in range(NNN)]
        #num_job = len(copy_config)
        job_server = pp.Server( ncpus=NNN)
        for i in range(NNN): 
            config = copy_config(CFG_ISING_BASIC)
            config['updaters'] = updaters_u1
            main = Main(**config)
            def run(q_iter): 
                main.run(q_iter=q_iter, backup_parpath=None)
            job_server.submit(run, (10, ), modules=('main', )) 
            #run(10)
    
    def test_run_many(self): 
        h = [1.0, 1.1] 
        NNN = len(h)
        from merapy.config import CFG_ISING_BASIC, copy_config, updaters_u1 
        config = copy_config(CFG_ISING_BASIC) 
        from merapy.finite_site import finite_site_u1
        config['updaters'] = updaters_u1
        config_group = [copy_config(config) for i in range(NNN)]
        
        for c in config_group: 
            c['SYMMETRY'] = 'Travial'
            c.update(tot_layer=3)
            c['backup_parpath'] = mkdtemp()
            c['schedule']['schedule'] = copy_schedule(schedule_prod_state)
            c['schedule'].update(q_iter_max=10, mera_shape_min=(3, 3),  mera_shape_max=(4, 3))
        Main.run_many(config_group, parallel=1)

    def xtest_transfer_pickle_file(self): 
        config = self.config_heisbg
        config['SYMMETRY'] = 'Travial'
        print config['backup_parpath']
        print config['backup_parpath_local']
        
        with make_temp_dir() as d1, make_temp_dir() as d2: 
            config['backup_parpath'] = d1
            config['backup_parpath_local'] =   d2
            
            bac = Main.LOCALHOSTNAME
            Main.LOCALHOSTNAME = 'changed it for a while'
            main=Main(**config)
            main.run_schedule(mera_shape_max=(4, 4), q_iter_max=5)
            self.assertTrue(os.path.exists(d2 + '/4.pickle'))
        print 'file exits  ---pass '
        Main.LOCALHOSTNAME = bac
    
    def xtest_run_many_dist(self): 
        ################## NOTE that it will raise if run this test under main.py  #################
        #### see the doc string of the module for details 
        if 1: 
            h = [1.0, 1.1] 
            NNN = len(h)
            from merapy.config import CFG_ISING_BASIC, copy_config, updaters_u1 
            config = copy_config(CFG_ISING_BASIC) 
            from merapy.finite_site import finite_site_u1
            config['updaters'] = updaters_u1
            config_group = [copy_config(config) for i in range(NNN)]
            
            for c in config_group: 
                c['SYMMETRY'] = 'Travial'
                c.update(tot_layer=3)
                c['backup_parpath'] = mkdtemp()
                c['schedule']['schedule'] = copy_schedule(schedule_prod_state)
                c['schedule'].update(q_iter_max=10, mera_shape_min=(3, 3),  mera_shape_max=(4, 3))
            servers=[('localhost', 90900),] 
            Main.run_many_dist(config_group, servers=servers)
        
        if 0:  #submit 
            h = [1.0, 1.1] 
            NNN = len(h)
            from merapy.config import CFG_ISING_BASIC, copy_config, updaters_u1 
            config = copy_config(CFG_ISING_BASIC) 
            from merapy.finite_site import finite_site_u1
            config['updaters'] = updaters_u1
            config_group = [copy_config(config) for i in range(NNN)]
            
            for c in config_group: 
                c['SYMMETRY'] = 'Travial'
                c.update(tot_layer=3)
                c['backup_parpath'] = mkdtemp()
                c['schedule']['schedule'] = copy_schedule(schedule_prod_state)
                c['schedule'].update(q_iter_max=10, mera_shape_min=(3, 3),  mera_shape_max=(4, 3))
            Main.run_many_dist(config_group, submit=1)
        
        
    def test_temp(self): 
        pass
        
        
class TestRpyc(unittest.TestCase): 
    def setUp(self): 
        pass
    def teardown(self): 
        pass

if __name__=="__main__":
    warnings.filterwarnings('ignore')

    
    if 0: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        TestMain.test_temp=unittest.skip("skip test_temp")(TestMain.test_temp) 
        unittest.main()
       
    else:
        suite = unittest.TestSuite()
        add_list = [
        
        #'test_run',
        #'test_run_scale_invar',
        #'test_run_schedule',
        #'test_run_many',
        #'xtest_run_many_dist',
        #'test_resume',
        
        #'xtest_make_dir',
        #'xtest_run_on_remote',
        
        #'xtest_resume_func',
        #'xtest_transfer_pickle_file',
        #'xtest_set_num_of_threads_1',
        
        'test_temp',
                ]
        for a in add_list: 
            suite.addTest(TestMain(a))
        
        unittest.TextTestRunner().run(suite)
       







