#!/usr/bin/env python
#coding=UTF8

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
        
"""
import warnings 
import time
import datetime 
import pprint
import inspect
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
import rpyc
from plumbum import SshMachine
import itertools
import signal 
from tempfile import mkdtemp

uname=os.uname()[1]
if uname in ['QTG-WS1-ubuntu']: 
    import pp

from merapy.utilities import print_vars
from merapy.context_util import (make_temp_dir, rpyc_conn_local, 
        rpyc_conn_local_zerodeploy, rpyc_conn, rpyc_conn_zerodeploy)
from merapy.quantum_number import *
from merapy.mera import Mera
from merapy.hamiltonian import System
from merapy.tensor_py import iTensor, iTensorFactory
from merapy.system_parameter import SystemParam 
from merapy.scale_invariant import ScaleInvar
from merapy import common_util
from merapy import crandom
import merapy.schedule
from merapy.schedule import copy_schedule, schedule_scale_invar
from merapy.decorators import *
from merapy.measure_and_analysis.measurement import measure_S


"""
    generally, main func has 3 effects: 
        1. init all--- pass all params
        2. excute; 
        3. collect and manage result
    Q:
        should I use keyword args or kargs?
"""

"""
sys_param
model_param
mera_param = {"trunc_dim":2, "num_of_layer":6, "symmetry":None}
iter_param
other_param
"""

__all__ = ["Main"]

if 0:  #archive
    def run_many_1(common_cfg, param_name_list, custom_param_list,  
            schedule_list, root_dir, q_iter_max = None, parallel=False, do_measure=True, **kwargs): 
        num_param = len(param_name_list)
        if not parallel: 
            for i in range(len(custom_param_list)): 
                #main = main_list[i]
                cfg = common_cfg.copy()
                param = custom_param_list[i]
                if not hasattr(param, '__iter__'): param = [param]
                dic = dict(zip(param_name_list, param))
                cfg['model_param'].update(dic)
                if 1: 
                    param = map(str, param)
                    temp = zip(param_name_list, param)
                    temp=map(lambda x: '='.join(x), temp)
                    dir = '-'.join(temp)
                    dir = root_dir  + dir
                    if not os.path.exists(dir): 
                        os.system('mkdir %s'%dir)
                schedule = schedule_list[i]
                main = Main(**cfg)
                main.run_schedule(schedule=schedule, q_iter_max=q_iter_max, backup_parpath=dir, do_measure=do_measure, **kwargs)
                iTensor.buff_free()
                main.stop_player()           
                del main.S  # or else tensor buffer would exhaust?
                del main
                iTensor.buff_free()

        else: 
            # Create jobserver
            job_server = pp.Server()
            ncpu = 4
            jobs = []
            start_time = time.time()
            job_server.submit(part_sum, (starti, endi))
            print "Time elapsed: ", time.time() - start_time, "s"

    def run_many_2(common_cfg, param_name_list, custom_param_list,  
            schedule_list, root_dir, ncpu=1, q_iter_max=None, parallel=False, do_measure=True, **kwargs): 
        
        num_param = len(custom_param_list)
        print 'nnnn', num_param
        
        def run_one(i): 
                cfg = common_cfg.copy()
                param = custom_param_list[i]
                if not hasattr(param, '__iter__'): param = [param]
                dic = dict(zip(param_name_list, param))
                cfg['model_param'].update(dic)
                if 1: 
                    param = map(str, param)
                    temp = zip(param_name_list, param)
                    temp=map(lambda x: '='.join(x), temp)
                    dir = '-'.join(temp)
                    dir = root_dir  + dir
                    if not os.path.exists(dir): 
                        os.system('mkdir %s'%dir)
                schedule = schedule_list[i]
                main = Main(**cfg)
                main.run_schedule(schedule=schedule, q_iter_max=q_iter_max, backup_parpath=dir, do_measure=do_measure, **kwargs)
                iTensor.buff_free()
                main.stop_player()           
                del main.S  # or else tensor buffer would exhaust?
                del main
                iTensor.buff_free()
        
        def run_group(group): 
            print 'gggggg', group
            for i in group: 
                #run_one(i)
                print group
        
        #if not parallel: 
        if ncpu == 1:  
            for i in range(len(custom_param_list)): 
                run_one(i)
        else: 
            # Create jobserver
            job_server = pp.Server(ncpus=ncpu)
            print "Starting pp with", job_server.get_ncpus(), "workers"
            
            start_time = time.time()
            num_group = num_param/ncpu if num_param%ncpu == 0 else num_param/ncpu + 1 
            
            print 'num_param = %d, num_group = %d'%(num_param, num_group)
            for g in range(num_group): 
               
                group = range(g*ncpu, (g+1)*ncpu)
                print 'gggg', g, group
                
                #job_server.submit(run_group,  args=(group, ), depfuncs=(run_one, ), modules=())
                job_server.submit(run_group,  args=(group, ), depfuncs=(), modules=())
            print "Time elapsed: ", time.time() - start_time, "s"
            job_server.print_stats()

    def run_mainy_simple(config_list, schedule_list):
        num_job = len(config_list)
        job_server = pp.Server(ncpus=6)
        job_server.submit


    class Main_new(object):
        def __init__(self, updaters, trunc_dim=None, tot_layer=None, q_iter=8,q_lay=3, info=0,
            MODEL="Heisenberg", SYMMETRY="Z2", only_NN=True, only_NNN=False, use_player=True,
            filename=None, J_NNN = None, USE_REFLECTION=None, 
            model_param=None, qsp_max=None):

            if filename == "default":
                import __main__
                filename = __main__.__file__.replace("main", "res").replace("run", "res").replace(".py", ".dat")
                #filename = __main__.__file__.replace(".py", ".dat")
                #filename = "res-" + filename
            else:
                pass
            

            if 0:
                args= ["tot_layer", "trunc_dim", "USE_REFLECTION"]
                loc = locals()
                print loc
                for a in args:
                    if a in loc:
                        if loc[a]:
                            SystemParam.__setattr__(a, loc[a])
            
            if 1:
                #if lang_tensor:
                #    SystemParam.lang_tensor = lang_tensor
                #could be removed
                if tot_layer is not None:
                    SystemParam.tot_layer = tot_layer
                if trunc_dim is not None:
                    SystemParam.trunc_dim = trunc_dim
                
                #can't be removed
                if USE_REFLECTION is not None:
                    SystemParam.USE_REFLECTION = USE_REFLECTION
                    if SystemParam.USE_REFLECTION  and trunc_dim != 2:
                        raise ValueError("USE_REFLECTION is not compatable with trunc_dim")
                
                args=["MODEL", "SYMMETRY", "only_NN"]
                temp=locals()
                for a in args:
                    if temp[a] is not None:
                        setattr(SystemParam, a, temp[a])

                SystemParam.system_cfg()
                if model_param is not None:
                    SystemParam.model_param.update(model_param)

                SystemParam.show(filename)    
            
            LayStart=0; q_one=1 
            layer = tot_layer
            symmetry = SYMMETRY

           
            varables= vars()
            print varables.keys()
            
            sys_param_name = ["model", "symmetry", "only_NN", "only_NNN"]
            self.sys_param = {i:vars()[i] for i in sys_param_name if vars().has_key(i)}
            
            #mera_param_name = ["trunc_dim", "tot_layer", "symmetry"]
            #self.mera_param = {i:varables[i] for i in mera_param_name}
            mera_inspect = inspect.getargspec(Mera.__init__)
            self.mera_param_name = mera_inspect.args
            print self.mera_param_name
            self.param_main_names = inspect.getargspec(Main.__init__).args
            #print Main.__init__.__dict__
            self.param_main_names.remove("self")
            print self.param_main_names
            print updaters
            self.param_main = {i:varables[i] for i in self.param_main_names}
            self.param_main.update(vars())

            print self.param_main #.args

            #self.mera_param = {i:varables[i]}

            iter_param_name = ["q_iter", "q_one", "q_lay"]
            self.iter_param = {i:varables[i] for i in iter_param_name}

            self.init_system()

        def init_system(self):
            """
                set the attribute self.S,  self.M
            """
            SYMMETRY, trunc_dim = self.param_main["SYMMETRY"], self.param_main["trunc_dim"]
            qn_identity, QSp_base, qsp_null=init_System_QSp(symmetry = SYMMETRY)
            #print QSp_base
            qsp_0 = QSp_base.copy()
            reset_System_QSp(symmetry=SYMMETRY)
            qsp_max, qsp_max2=QSp_base.__class__.max(trunc_dim)

            self.S= System(SYMMETRY, **self.sys_param)
            

            SystemParam.ReAssumed= False
            topQN= qn_identity.copy()

            
            nTop = None
            self.param_main.update(vars())
            crandom.mrandom.csrand(1234)

            param_mera_names = inspect.getargspec(Mera.__init__).args
            param_mera_names.remove("self")
            self.param_mera_dic = {i:self.param_main[i] for i in param_mera_names}
            
            print param_mera_names
            if not SystemParam.ReAssumed:
                self.M= Mera(**self.param_mera_dic)
                self.S.init_Hamiltonian(self.M)
            else:
                hamiltonian.Restore_System(SystemParam.BackupFile,M,sys)
            

            SystemParam.ReAssumed=False
            
            #每次更改耦合常数算新的值，要重新初始化
            if self.S.model == "Ising" :
                #self.S.init_Hamiltonian(self.M)
                for iLayer in range(self.M.num_of_layer-1): 
                    #top V shouldn't inited by init_layer
                    self.M.init_layer(iLayer, True)

        def run(self):
            lstart=0
            lstep= 0
            cont = None
            q_one=1 
            
            ascending, descending, update_mera, finite_site, rho_top_func= tuple(self.param_main["updaters"])
            M = self.M
            #S= self.sys
            self.param_main.update(vars())

            self.param_iter_names = inspect.getargspec(finite_site).args        
            self.param_iter_dic= {i:self.param_main[i] for i in self.param_iter_names}
            finite_site(**self.param_iter_dic)

if 1: 
    import copy_reg
    import types
 
    # to make  multiprocessing compatable with pickleing of method funcs
    if 0: 
        def _pickle_method(method):
            # Author: Steven Bethard
            # http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
            func_name = method.im_func.__name__
            obj = method.im_self
            cls = method.im_class
            cls_name = ''
            if func_name.startswith('__') and not func_name.endswith('__'):
                cls_name = cls.__name__.lstrip('_')
            if cls_name:
                func_name = '_' + cls_name + func_name
            return _unpickle_method, (func_name, obj, cls)

        def _unpickle_method(func_name, obj, cls):
            # Author: Steven Bethard
            # http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
            for cls in cls.mro():
                try:
                    func = cls.__dict__[func_name]
                except KeyError:
                    pass
                else:
                    break
            return func.__get__(obj, cls)

        # This call to copy_reg.pickle allows you to pass methods as the first arg to
        # mp.Pool methods. If you comment out this line, `pool.map(self.foo, ...)` results in
        # PicklingError: Can't pickle <type 'instancemethod'>: attribute lookup
        # __builtin__.instancemethod failed

        copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
    
    else: 
        #the following correctly handling classmethods
        def _pickle_method(method):
            func_name = method.im_func.__name__
            obj = method.im_self
            cls = method.im_class
            if func_name.startswith('__') and not func_name.endswith('__'):
                #deal with mangled names
                cls_name = cls.__name__.lstrip('_')
                func_name = '_%s%s' % (cls_name, func_name)
            return _unpickle_method, (func_name, obj, cls)

        def _unpickle_method(func_name, obj, cls):
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
    #used in distributed computing 
    LOCALHOSTNAME = 'QTG-WS1-ubuntu'
    
    def __init__(self, system_class=None, mera_class=None, updaters=None, trunc_dim=None, tot_layer=None, unitary_init=None, q_iter=8, q_lay=3, info=0,
        MODEL="Heisenberg", SYMMETRY="Z2", only_NN=True, only_NNN=False, use_player=True, message=None, 
        filename=None, backup_fn='auto', backup_parpath=None, resume=False, auto_resume=True, J_NNN = None, USE_REFLECTION=None, 
        NUM_OF_THREADS=1, RAND_SEED=1234, USE_CUSTOM_RAND=True, 
        model_param=None, energy_exact=None, combine_2site=True, qsp_base=None, qsp_max=None, 
        mera_kwargs={}, sys_kwargs={}, **kwargs):
        
        self.hostname = socket.gethostname()
        if 1: 
            self.system_class = system_class if system_class is not None else System
            self.mera_class = mera_class if mera_class is not None else Mera
            self.updaters= updaters
            self.trunc_dim = trunc_dim
            self.tot_layer = tot_layer
            self.unitary_init = unitary_init
            self.q_iter = q_iter
            self.q_lay = q_lay
            self.info = info
            self.MODEL= MODEL
            self.SYMMETRY= SYMMETRY
            self.only_NN = only_NN
            self.only_NNN = only_NNN
            self.use_player = use_player
            self.message = message
            self.filename = filename
            self.backup_fn = backup_fn
            self.backup_parpath = backup_parpath
            self.USE_REFLECTION = USE_REFLECTION
            self.NUM_OF_THREADS= NUM_OF_THREADS
            self.RAND_SEED = RAND_SEED
            self.USE_CUSTOM_RAND = USE_CUSTOM_RAND
            self.model_param = model_param
            self.resume = resume
            self.auto_resume = auto_resume
            self.energy_exact = energy_exact
            self.combine_2site = combine_2site
            self.qsp_base = qsp_base
            self.qsp_max = qsp_max
            self.mera_kwargs = mera_kwargs
            self.sys_kwargs= sys_kwargs
            self.do_measure = kwargs.get('do_measure')
            if self.do_measure is None: 
                self.do_measure = False
        
        #above is tedious, below init self in a more auto fashion
        temp = ['backup_parpath_local', 'measurement_args', 'use_local_storage']
        for t in temp: 
            setattr(self, t, kwargs.get(t))

        if self.filename == "default":
            import __main__
            self.filename = __main__.__file__.replace("main", "res").replace("run", "res").replace(".py", ".dat")
        
        #issue in each main, should use independent random state!     
        if USE_CUSTOM_RAND:
            crandom.mrandom.csrand(RAND_SEED)
        else:
            random.seed()
            crandom.rand = random.random

        if 1:
            #if lang_tensor:
            #    SystemParam.lang_tensor = lang_tensor
            #could be removed
            if tot_layer is not None:
                SystemParam.tot_layer = tot_layer
            if trunc_dim is not None:
                SystemParam.trunc_dim = trunc_dim
            
            #can't be removed
            if USE_REFLECTION is not None:
                SystemParam.USE_REFLECTION = USE_REFLECTION
                if SystemParam.USE_REFLECTION  and self.trunc_dim != 2:
                    raise ValueError("USE_REFLECTION is not compatable with trunc_dim")
            
            if 0: #del
                if NUM_OF_THREADS is None:
                    NUM_OF_THREADS= 1
                elif NUM_OF_THREADS== "auto":
                    if trunc_dim<= 4:
                        NUM_OF_THREADS= 1
                    else:
                        NUM_OF_THREADS= 6
                elif isinstance(NUM_OF_THREADS, int):
                    pass
                else:
                    raise ValueError("NUM_OF_THREADS=%s"%NUM_OF_THREADS)
                os.environ["OMP_NUM_THREADS"] = str(NUM_OF_THREADS)
            
            Main.set_num_of_threads(NUM_OF_THREADS)
        if 1:
            args=["MODEL", "SYMMETRY", "only_NN", "only_NNN", "NUM_OF_THREADS"]
            temp=locals()
            for a in args:
                val = self.__getattribute__(a)
                #print a,  val
                if  val is not None:
                    setattr(SystemParam, a, val)
        
        SystemParam.system_cfg()
        
        if model_param is not None:
            SystemParam.model_param.update(model_param)
            
        
        self.init_all()

    def init_all(self):
        qn_identity, qsp_base, qsp_null=init_System_QSp(symmetry = self.SYMMETRY, combine_2site=self.combine_2site)
        
        if self.qsp_base is None: self.qsp_base = qsp_base
        self.qn_identity = qn_identity
        self.qsp_null = qsp_null
        qsp_0 = self.qsp_base.copy()
        if not self.combine_2site and self.SYMMETRY == "U1" :
            qsp_max = qsp_base.__class__.max(self.trunc_dim, combine_2site=self.combine_2site)
        else:
            qsp_max= qsp_base.__class__.max(self.trunc_dim)
        qsp_max2 = qsp_max.copy() 
        
        if self.qsp_max is not None:
            qsp_max = self.qsp_max.copy()
            qsp_max2 = None #self.qsp_max.copy()
            print "trunc_dim is updated to %d"%qsp_max.totDim
            self.trunc_dim = qsp_max.totDim

        topQN = qn_identity.copy()
        
        if 1:
            self.M= self.mera_class(self.tot_layer, SystemParam.nTop, 
                    qsp_0, qsp_max, qsp_max2, topQN, qsp_null, qn_identity, 
                    unitary_init=self.unitary_init, **self.mera_kwargs)
            #print 'sssss', type(self.system_class),  self.sys_kwargs
            
            
            self.S = self.system_class(mera=self.M, model_param=self.model_param,  
                    symmetry=self.SYMMETRY, model=self.MODEL, 
                    only_NN = self.only_NN, only_NNN = self.only_NNN,  
                    energy_exact=self.energy_exact, combine_2site=self.combine_2site, 
                    **self.sys_kwargs)
            
            self.S.init_Hamiltonian()
            

            #每次更改耦合常数算新的值，要重新初始化
            if self.S.model == "Ising":
                self.S.init_Hamiltonian()
                for iLayer in range(self.M.num_of_layer-1): 
                    #top V shouldn't inited by init_layer
                    self.M.init_layer(iLayer, rand_init=True, unitary_init=self.unitary_init)
        
        LayStart=0
        q_one=1 

    def finite_site_u1_del(self, q_iter, q_lay=1, q_one=1, backup_fn=None, 
            lstart=0, lstep=0, resume=False, auto_resume=True, info=0):
        """
            ref:
                参考 Evenbly, G., and G. Vidal. "Algorithms for entanglement renormalization." Physical Review B 79.14 (2009): 144108.
                p.17 A1-A4.  
                与Vidal文中稍有不同的是，see my note n5 at page 17
            todo:
                应该还用Vidal文中的迭代顺序，即先把所有的rho求出来，这样更舒服


            总的来说，要更新h2, rho, mera 这三部分，要合理安排三类算子，以及每类算子的每一曾的更新次序。
            更新h2, 必须要有mera；
            更新mera，必须要有h2和rho；
            更新rho，实际上rho不是独立的，完全由mera决定
            只有h2[0]和mera是实际的矩阵，而其他都是收缩的中间过程, 是虚的
            conduct a complete optimization for the system

            q_iter: wee vidal p.18 qiter, number to sweep  optimizing all u and w in whole mera
            q_lay: see vidal p.18 qlay, the number to repeat optimizing all u and w within that layer
            q_one: see vidal p.16 qone, number to iterate to optimize a single w or u

        """
        if isinstance(self.updaters, list): #only for backward compatable
            ascending, descending, update_mera, finite_site_u1, rho_top_func= tuple(self.updaters)
        elif isinstance(self.updaters, dict):
            ascending = self.updaters["ascending_func"]
            descending = self.updaters["descending_func"]
            update = self.updaters["update_func"]


        use_player = self.use_player
        filename = self.name
        M = self.M
        S= self.S

        #err_tol = 1.E-7
        if resume:    
            auto_resume = False
        iter0 = S.iter
        # iter=0 is init, original value
        #S.iter precisely means the state of S after iter times update, e.g. S.iter = 0 means original state
        
        #note here start is no longer controled by the passed in param start; remove that param in future
        if iter0 == 0:
            if rho_top_func.__name__ in ["top_level_product_state", "top_level_product_state_u1", "top_level_product_state_u1_1site"]:
                rho_top_func(M, S)     #rho_top is never changed when using top_level_product_state
            start = True 
        elif iter0 > 0:
            start = False
        
        
        for iter  in xrange(iter0 + 1, q_iter + iter0 + 1):
            if info>0:
                print str(iter) + "th iteratiion"
            
            if use_player:
                if iter0 == 0:
                    set_STATE_end_1(iter, record_at=iter0 + 2, power_on=True)
                else:
                    if resume:
                        set_STATE_end_1(iter, record_at=iter0 + 1, power_on=True)
                    else:
                        set_STATE_end_1(iter, record_at=-1, power_on=True)  #no more record, start play immediatly


            #第一遍：从下到上把H_2求出来
            #注意对于finite range，在fortran的算法里，最上层的V张量是个摆设，不是实际的最上层!, 它的值没有被更新过
            for ilayer in range(lstart, M.num_of_layer-1):
                #note the update order here, start from H then M then H:  H[0]->M[0]->H[1]->M[1]->....
                #H[0] is inited before finite_site
                if not start:   
                    for iter_lay in range(q_lay):
                        #iterative_optimize_all(M,S,ilayer,j=0, info=info-1)
                        update_mera(M,S,ilayer,j=0, info=info-1)    #update U[ilayer] etc
                ascending(M, S, ilayer, info=info-1)   #calculate H[ilayer+1] for ilayer
            
            #update rho_2[M.num_of_layer-1]        
            if rho_top_func.__name__ == "top_level_eigenstate": 
                rho_top_func(M, S)

            #这里似乎不一定在top layer收缩，可以任意layer
            
            #S.iter = iter
            S.iter += 1 
            S.eng_ham(M.num_of_layer-1)
            S.save_eng(S.iter)
            steps, dt_real= S.display(S.iter, filename)
            if steps*dt_real>86400*10: #(10 days)
                if iter-iter0>10:
                    print "\nENERGY_DIFF IS TOO SMALL, BREAK AT HERE. steps= %(steps)d\n"%vars()
                    #print "energy_diff is too small, break at here. steps= %(steps)d"%vars()
                    break

            if auto_resume and backup_fn is not None:
                try:
                    S0= System.load(backup_fn)
                except IOError:
                    print "try to resume but file %s not found"%backup_fn
                else:
                    if S.energy> S0.energy:
                        print S
                        S= S0
                        M = S0.mera
                        S.iter = S0.iter
                        print "aotu_resumed from %s at a lower enegy=%2.15f"%(backup_fn, S0.energy)
                finally:
                    auto_resume = False


            #这一步更新rho[num_of_layer-2] ,..., [1]
            for ilayer  in range(M.num_of_layer-1, lstart + 1, -1):
                #rho[0]可以求出来，但它在更新U，V时用不到，我们也不用它来收缩算能量，故不求它
                for iter_lay in range(q_lay):
                    #calculate reduced density mat. for each layer
                    descending(M,S,ilayer,  info=info-1)

            start =  False 
            
            if S.iter%50 == 0 and iter-iter0>10:
                if backup_fn is not None:
                    S.save(fn=backup_fn)

    #@profileit(fn="profile.tmp")
    def run(self, resume=True, auto_resume=None, backup_fn=None, filename=None, 
            q_iter=None, do_measure=None, **kwargs):
        LayStart=0
        q_one=1
        backup_fn_local = None
        do_measure  = do_measure if do_measure is not None else self.do_measure
        
        backup_parpath = self.backup_parpath
        backup_parpath_local = self.backup_parpath_local
        if kwargs.get('backup_parpath') is not None : 
            backup_parpath = kwargs.get('backup_parpath')    
        self.make_dir(backup_parpath, backup_parpath_local)
        
        fn = self.S.backup_fn_auto() 
        if backup_parpath is not None:
            
            backup_fn = backup_parpath  + '/' +  fn
        if backup_parpath_local is not None: 
            backup_fn_local = '/'.join([backup_parpath_local, fn])
           
        filename = 'res-'  + self.S.backup_fn_auto()
        filename = filename.split('.')[0] + '.dat'
        if backup_parpath is not None: 
            filename = backup_parpath + '/' + filename
        self.filename = filename 

        if q_iter is None:
            q_iter = self.q_iter
        if resume is not None:
            self.resume  = resume
                
        
        if self.resume and backup_fn is not None:
            if not self.use_local_storage: 
                self.resume_func(backup_fn)
            else: 
                self.resume_func(backup_fn_local, use_local_storage=self.use_local_storage)
            if hasattr(self.S, 'energy_diff_std') and kwargs.get('energy_diff_min') is not None : 
                energy_diff_min = kwargs.get('energy_diff_min')
                if  self.S.energy_diff_std is not None: 
                    if self.S.energy_diff_std <= energy_diff_min: 
                        q_iter = 0
                        do_measure = 0
                        msg = 'self.S.energy_diff_std %1.2e less than energy_diff_min %1.2e, set q_iter=0, do_measure=0'%(
                                self.S.energy_diff_std, energy_diff_min)
                        print msg
           
        if auto_resume is not None:   self.auto_resume = auto_resume
        
        SystemParam.param_extra.update({
                "resume":self.resume, 
                "RAND_SEED":self.RAND_SEED, 
                "USE_CUSTOM_RAND":self.USE_CUSTOM_RAND, 
                "unitary_init":self.M.unitary_init, 
                #"top_level_func":self.updaters[4].func_name, 
                "message":self.message, 
                "energy_exact":self.energy_exact
                })
        SystemParam.show(self.filename)    
        
        if isinstance(self.updaters, list): #only for backward compatable
            
            ascending_func, descending_func, update_mera_func, finite_range_func, rho_top_func= tuple(self.updaters)
        elif isinstance(self.updaters, dict):
            
            finite_range_func = self.updaters["finite_range_func"]
            ascending_func = self.updaters["ascending_func"]
            descending_func = self.updaters["descending_func"]
            update_mera_func = self.updaters["update_mera_func"]
            rho_top_func = self.updaters["rho_top_func"]
        
        
        temp = ['use_local_storage' ]
        args = {t: self.__getattribute__(t) for t in temp}
        args.update(backup_fn_local=backup_fn_local)
        kwargs.update(args)
       
        finite_range_func(self.M, self.S, ascending=ascending_func, descending=descending_func, 
                update_mera=update_mera_func, rho_top_func = rho_top_func, 
                q_one=q_one, q_lay=self.q_lay, q_iter=q_iter, use_player=self.use_player, 
                filename=self.filename, backup_fn=backup_fn, lstart = LayStart, lstep=0, 
                resume=self.resume, auto_resume=self.auto_resume, info=self.info-1, **kwargs)
        
        
        if do_measure:  
            ppp = backup_parpath if not self.use_local_storage else backup_parpath_local
            if ppp is not None : 
                exclude_which = self.measurement_args.get('exclude_which', [])
                exclude_which.append('scaling_dim')
                self.measure_S(self.S, ppp, exclude_which=exclude_which)
       
        return self.M, self.S
   
    def run_scale_invar(self, resume=True, auto_resume=None, backup_fn=None, filename="auto", 
            num_of_SIlayer=3, q_lay=1, q_iter=5, do_measure=None,  **kwargs):
        
        do_measure  = do_measure if do_measure is not None else self.do_measure
        
        q_one = 1
        
        backup_parpath = self.backup_parpath
        backup_parpath_local = self.backup_parpath_local
        if kwargs.get('backup_parpath') is not None : 
            backup_parpath = kwargs.get('backup_parpath')
        if kwargs.get('backup_parpath_local') is not None : 
            backup_parpath_local = kwargs.get('backup_parpath_local')
        
                
        fn = self.S.backup_fn_auto() 
        if backup_parpath is not None : 
           
            backup_fn  = backup_parpath + '/' + fn
        backup_fn_local = None
        if backup_parpath_local is not None: 
            backup_fn_local = '/'.join([backup_parpath_local, fn])
            

        #if filename is not None:   self.filename = filename
        filename = 'res-'  + self.S.backup_fn_auto()
        filename = filename.split('.')[0] + '.dat'
        
        
        if backup_parpath is not None: 
            filename = backup_parpath + '/' + filename
        self.filename = filename 
      

        self.make_dir(backup_parpath, backup_parpath_local)
        temp = backup_fn if not self.use_local_storage else backup_fn_local
        if resume and temp is not None:  
            self.resume_func(temp, use_local_storage=self.use_local_storage)
            if not hasattr(self.S, 'iter1'):
                self.S.iter1 = self.S.iter
            energy_diff_min  = kwargs.get('energy_diff_min')
            if hasattr(self.S, 'energy_diff_std') and kwargs.get('energy_diff_min') is not None : 
                if self.S.energy_diff_std <= energy_diff_min and self.S.energy_diff_std is not None  : 
                    q_iter = 0
                    do_measure = 0
                    msg = 'self.S.energy_diff_std %1.2e less than energy_diff_min %1.2e, set q_iter=0, do_measure=0'%(
                            self.S.energy_diff_std, energy_diff_min)
                    print msg
                
        if auto_resume is not None:   self.auto_resume = auto_resume
        
        temp = ['use_local_storage']
        args = {t: self.__getattribute__(t) for t in temp}
        args['backup_fn_local'] = backup_fn_local
        kwargs.update(args)

        SI = ScaleInvar(self.S, self.updaters, num_of_SIlayer, q_iter, q_one, q_lay=q_lay, 
                filename=filename, backup_fn=backup_fn, 
                resume=resume, auto_resume=auto_resume, 
                use_player=self.use_player, info=self.info, **kwargs)
        
        SI.scale_invariant()
        
        if do_measure: 
            ppp = backup_parpath if not self.use_local_storage else backup_parpath_local
            if ppp is not None : 
                exclude_which = self.measurement_args.get('exclude_which', [])
                self.measure_S(self.S, ppp, exclude_which=exclude_which)
        
        if self.backup_parpath is not None and self.backup_parpath_local is not None and not self.use_local_storage: 
            self.transfer_pickle_file()
        
        return self.S

    def run_schedule(self, schedule=None, mera_shape_min=(4, 4), mera_shape_max=(12, 4), q_iter_max=None, backup_parpath=None,  filename='auto', backup_fn='auto', 
           do_measure=None, skip_list=[],  skip_dim_list=[],  dim_diff_remap={}, info=0):
        """
            schedule is a OrderedDict like
            (ansatz, dim, layer, q_iter, run_time_max, energy_diff_min)
            
            { 4:  ('ScaleInvar', '')}
            
            }
        """
        
        if schedule is None: 
            schedule = schedule_scale_invar
        schedule = copy_schedule(schedule, skip_list=skip_list)
        do_measure  = do_measure if do_measure is not None else self.do_measure
        
        if len(dim_diff_remap)>0: 
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
        qsp_class=symmetry_to_QspClass(symmetry)
        
        #for key, val in schedule.items():
        for task in schedule: 
            if task['ansatz'] == 'prod_state':  
                run_func = self.run
                run_param_name = ('q_iter', 'energy_diff_min')
            elif task['ansatz'] == 'scale_invar': 
                run_func = self.run_scale_invar
                run_param_name = ('q_iter', 'run_time_max', 'energy_diff_min')
            elif task['ansatz'] == 'eigenstate': 
                raise NotImplemented
            else: 
                raise
            
            if q_iter_max is not None : 
                if task['q_iter'] > q_iter_max: 
                    task['q_iter']  = q_iter_max
            shape = task['shape']
            
            msg = '\nSTARTING TASK %s'%(task, )
            if info>0: print msg
            
            dim, layer = shape[:2]
            nqn = None
            if len(shape)>2: nqn = shape[2]
            
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
                    nqn_max = None
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
                #print 'sssss', self.S, type(self.S)
                propt = self.S.key_property()
                dim_current = propt['trunc_dim']
                layer_current = propt['num_of_layer']
                nqn_current = propt['nqn']
                qsp_current = qsp_class.max(dim_current, nqn_current)
                shape_current = (dim_current, layer_current, nqn_current)
                #if dim_current <dim: 
                
                if qsp_current < qsp_destination: 
                    self.expand_dim(dim, nqn)
                
                if layer_current <layer: 
                    self.expand_layer(layer)
                
                
                #if dim_current>dim or layer_current> layer: 
                if qsp_current> qsp_destination or layer_current> layer: 
                    msg = 'shape %s smaller than current shape = %s, skip it'%(shape, shape_current)
                    print msg
                    continue
                   
            run_param = dict([(i, task.get(i)) for i in run_param_name])
            run_func(self, backup_parpath=backup_parpath, filename=filename, backup_fn=backup_fn, 
                    do_measure=do_measure, **run_param)
    
    def make_dir(self, dir, dir_local=None): 
        if dir is not None : 
            if not os.path.exists(dir): 
                os.makedirs(dir)
                print 'dir not found, mkdir %s'%dir
        
        if dir_local is not None: 
            with rpyc_conn_local_zerodeploy() as conn: 
                os_local = conn.modules.os
                if not os_local.path.exists(dir_local): 
                    msg = 'dir %s not found'%(dir_local[-30: ])
                    cmd = 'mkdir %s'%dir_local
                    os_local.makedirs(dir_local)
                    print '\n'.join([msg, cmd])
    
    def measure_S(self, S, parpath, exclude_which=None):
        #state_bac = tensor_player.STATE
        self.stop_player()
        msg = '\nmeasurement after final iter: '
         
        measure_S(S,  parpath, exclude_which=exclude_which, use_local_storage=self.use_local_storage)
        #tensor_player.STATE = state_bac

    def stop_player(self):
        set_STATE_end_1(iter=None, record_at=None, stop_at=None, power_on=False)

    @staticmethod
    def check_hermite(M, S, precision=None):
        print "\ncheck hermite"
        state = tensor_player.STATE
        tensor_player.STATE = "stop"
        for lay in range(M.num_of_layer-0):
        #if 1:
            #lay = M.num_of_layer-1
            print lay
            print "H_2", S.H_2[lay][0].is_hermite(precision=precision, out_more=True)
            if not S.only_NN:
                print "H_3", S.H_3[lay][0].is_hermite(precision=precision, out_more=True)
                if not S.only_NNN:
                    if S.model == "Ising":
                        print "SxA", S.SxA[lay][0].is_hermite(precision)
                        print "SxB", S.SxB[lay][0].is_hermite(precision)
                    else:
                        print "SpA", S.SpA[lay][0].is_hermite(precision)
                        print "SpB", S.SpB[lay][0].is_hermite(precision)
                        print "SmA", S.SmA[lay][0].is_hermite(precision)
                        print "SmB", S.SmB[lay][0].is_hermite(precision)
                        
                        print "SpA adjoint to SmA", S.SpA[lay][0].is_adjoint_to(S.SmA[lay][0], precision=precision, out_more=True)
                        print "SpB adjoint to SmB", S.SpB[lay][0].is_adjoint_to(S.SmB[lay][0], precision=precision, out_more=True)
        tensor_player.STATE = state
    
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
    def set_num_of_threads(n):
        #import common_64_ifort as c64
        common_util.set_num_of_threads(n)
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
    
    @staticmethod
    def get_num_threads_max(): 
        return psutil.NUM_CPUS
        
    @staticmethod
    def init_quaternary_del(init_mera_graph=None):
        if init_mera_graph is None:
            import merapy.init_mera_graph as init_mera_graph
        import merapy.quantum_number_py
        from merapy.diagrams.V41.graph_quaternary import G_2_2, G_3_2, G_22_2, G_22_3

        merapy.quantum_number_py.QspU1.MaxQNNum = 30   #increase this because legs of V increases
        print "set QspU1.MaxQNNum to 30"


        print "replaced init_mera_graph.G_2_2 with quaternary conterpart"
        init_mera_graph.G_2_2 = G_2_2
        init_mera_graph.G_3_2 = G_3_2
        init_mera_graph.G_3_3 = {}
        init_mera_graph.G_22_2 = G_22_2
        init_mera_graph.G_22_3 = G_22_3
        print "using quaternary graph  "*10
    
    def expand_layer(self, lay):
        self.S.expand_layer(lay)
        
    def expand_dim(self, dim, nqn=None):
        self.S.expand_dim(dim, nqn)
    
    @staticmethod
    def load(fn):
        inn = open(fn, "rb")
        res = pickle.load(inn)
        return res
    
    
    def __del__(self, ): 
        print 'exit main, and run main.stop_player'
        self.stop_player()
    
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
                msg += 'job exists' 
                
        except KeyboardInterrupt: 
            #print 'KeyboardInterrupt in child for %s' %(config.get('model_param'), )
            msg += 'KeyboardInterrupt captured.'
            job_status = 'FAILED'
        except Exception, exception:
            print exception
            traceback.print_exc()
            job_status = 'FAILED'
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
                    nproc = 2
                else: 
                    warnings.warn('unfamiliar hostname %s'%(hostname, ))
                    nproc = ncpu_tot/2//nt 
                assert nproc>0, 'nproc=%d'%(nproc, )
                
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
    
    @staticmethod
    def make_grid(*args): 
        grid = itertools.product(*args)
        return list(grid)  #efficiency doesnt matter,  list is better
    grid = make_grid  # another name
    @staticmethod
    def grid_merger(*args): 
        return itertools.chain(*args)
    
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
    
    def resume_func(self, backup_path, use_local_storage=False):
        #if not use_local_storage: 
        if 1: 
            if backup_path is not None : 
                LOCAL = '' if not use_local_storage else 'LOCAL file '
                msg = "\nTRY TO RESUME from %s%s...  "%(LOCAL, backup_path[-30:])
                try:
                    self.S = self.system_class.load(backup_path, use_local_storage=use_local_storage)
                    self.M = self.S.mera
                    iter = self.S.iter1
                    eng = self.S.energy
                    msg  += "RESUMED successfully at iter=%s,  eng=%s"%(iter, eng)
                except IOError as err:
                    if err.errno == 2: 
                        msg += str(err)
                        msg += "discard manually resume"
                    else:
                        raise IOError
                except TypeError as err:
                    msg += str(err)
                    msg += "discard manually resume"
                print(msg)

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
    
    @staticmethod
    def email(): 
        import smtplib
        import email.mime.text
        # my test mail
        #mail_username='chiwhalee@gmail.com'
        mail_username='zhihuali@mail.ustc.edu.cn'
        mail_password='****'
        from_addr = mail_username
        to_addrs=('chiwhalee@gmail.com')

        # HOST & PORT
        #HOST = 'smtp.gmail.com'
        HOST = 'smtp.ustc.edu.cn'
        PORT = 25

        # Create SMTP Object
        smtp = smtplib.SMTP()
        print 'connecting ...'

        # show the debug log
        smtp.set_debuglevel(1)

        # connet
        try:
            print smtp.connect(HOST,PORT)
        except:
            print 'CONNECT ERROR ****'
        # gmail uses ssl
        #smtp.starttls()
        # login with username & password
        try:
            print 'loginning ...'
            smtp.login(mail_username,mail_password)
        except:
            print 'LOGIN ERROR ****'
        # fill content with MIMEText's object 
        msg = email.mime.text.MIMEText('Hi ,I am leehark')
        msg['From'] = from_addr
        msg['To'] = ';'.join(to_addrs)
        msg['Subject']='hello , today is a special day'
        print msg.as_string()
        smtp.sendmail(from_addr,to_addrs,msg.as_string())
        smtp.quit()


class TestMain(unittest.TestCase): 
    def setUp(self): 
        from config import CFG_HEISBG_BASIC, updaters_u1
        temp = dict(USE_CUSTOM_RAND=True, updaters=updaters_u1, trunc_dim=4, tot_layer=4, use_player=True, 
                SYMMETRY="Travial", 
                NUM_OF_THREADS=1, do_measure=0)
        cfg_heisbg = CFG_HEISBG_BASIC.copy();   cfg_heisbg.update(temp)
        #model_param={"J_NN":1.0, "J_NNN":0.241186}
        #cfg_heisbg['model_param'].update(model_param)
        #if 1:  #test run_schedule
        try: 
            os.system('rm ./mera_backup_test_folder/*')    
            os.system('rm 4.pickle')
        except: 
            pass 
        
        self.config = cfg_heisbg
    
    def test_temp(self): 
        from merapy.run_ising.config import make_config 
        import numpy as np 
        hh = np.arange(0, 2.2, 0.1)
        #hh = []
        config_group = []
        for h in hh: 
            c = make_config(h, root1= '')
            c['schedule'].update(mera_shape_min=(4, 4), mera_shape_max=(4, 4), 
                    q_iter_max = 5, do_measure = 0)
            c['NUM_OF_THREADS'] = 1
            config_group.append(c)
            c['backup_parpath_local'] = '/tmp/sdfjl'
            c['register_job'] = 1
            print  c['backup_parpath']

        Main.run_many(config_group, nproc=8)

    def test_grid(self): 
        a = [1, 2]
        b = ['a', 'b']
        g1 = Main.grid(a, b)
        g2 = Main.grid(['x'], [0])
        gg = Main.grid_merger(g1, g2)
        for g in gg: 
            print g 
        
    def test_run(self): 
        config = self.config 
        config['SYMMETRY'] = 'Travial'
        main = Main(**config)
        main.run(q_iter=10, backup_parpath='./mera_backup_test_folder', do_measure=0)
    
    def test_run_on_remote(self): 
        config = self.config
        config['SYMMETRY'] = 'Travial'
        #with rpyc_conn('sugon') as conn: 
        with rpyc_conn_zerodeploy('sugon') as conn: 
            conn.namespace['self'] = self
            code = r"""from merapy.main import TestMain; self.test_run()"""
            #Main_remote = conn.modules['merapy.main'].Main
            conn.execute(code)
            print conn.modules.os.getcwd()
    
    def test_resume_func(self): 
        config = self.config 
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
        
    def test_make_dir(self): 
        config = self.config 
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
        config = self.config
        config['SYMMETRY'] = 'U1'
        config['measurement_args']['exclude_which'].append('correlation_extra')  # too slow
        main = Main(**config)
        with make_temp_dir() as dd: 
            main.run_scale_invar(q_iter=5, backup_parpath=dd, do_measure=1)
        
    def test_run_schedule(self): 
        from merapy import schedule
        print '-------'*30
        #main = self.main
        config = self.config
        config['SYMMETRY'] = 'Travial'
        main = Main(**config)
        main.run_schedule(schedule=schedule.schedule_scale_invar, do_measure=0, backup_fn='auto', 
                backup_parpath='./mera_backup_test_folder', 
                q_iter_max = 10, 
                mera_shape_min = (4, 4), mera_shape_max=(8, 5), skip_dim_list=[14])
        S = main.S
        self.assertEqual(S.iter, 40)
        self.assertEqual(S.iter1, 10)
        print 'pass tests for right num of iter'

    def test_set_num_of_threads_1(self): 
        config = self.config
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
        
        for i in config_group: 
            id = config_group.index(i)
            #i['backup_parpath'] = '/home/zhli/Documents/mera_backup_tensor/mera_backup_test_folder/%d/'%id  
            i['backup_parpath'] = mkdtemp()
        Main.run_many(config_group, parallel=1)

    def test_transfer_pickle_file(self): 
        config = self.config
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
    
    
class TestRpyc(unittest.TestCase): 
    def setUp(self): 
        pass
    def teardown(self): 
        pass

if __name__=="__main__":
    
    if 0: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
       
    else:
        suite = unittest.TestSuite()
        add_list = [
        #'test_temp',
        #'test_run',
        #'test_run_scale_invar',
        #'test_run_schedule',
        #'test_make_dir',
        #'test_run_on_remote',
        #'test_resume_func',
        #'test_transfer_pickle_file',
        #'test_set_num_of_threads_1',
        'test_run_many',
                ]
        for a in add_list: 
            suite.addTest(TestMain(a))
        
        unittest.TextTestRunner().run(suite)
       







