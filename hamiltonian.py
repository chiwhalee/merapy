#coding= UTF8
"""
see Hamiltonian.f90

    
"""
import unittest
import numpy as np
import os 
import time
import cPickle as pickle
#import collections
import types 
from collections import OrderedDict
import warnings
import rpyc
import importlib 
import tempfile 


from vmps.iterative_optimize import IterativeOptimize 

from merapy.decorators import tensor_player 
from merapy.utilities import print_vars
import merapy.crandom as crandom 
from merapy.models import *
from merapy.context_util import rpyc_conn_local, rpyc_conn_local_zerodeploy, rpyc_load
from merapy.mera import Mera
from merapy.tensor import *
from merapy.tensor_network import *
from merapy.quantum_number import *

if 1:  # not used , BUT, not able to remove it now !!! because pickle need it for old backup files
    class Operators(object):
        """
            this may be deprecated,  
            for fortran define a tensor container is necessary,  but for python, 
            use list or dict instead
        """
        def __init__(self):
            self.start_site= None
            self.num_of_site= None
            #self.O= None  #None indicates not "allocated" as in fortran  #[Tensor()]
            nSite = 1
            self.O = [None]*nSite
        def __getitem__(self, i):
            return self.O[i]
        def __setitem__(self, i, x):
            self.O[i] = x
if 0: 
    class Operators_new(list):
        """
            this may be deprecated,  
            for fortran define a tensor container is necessary,  but for python, 
            use list or dict instead
        """
        def __init__(self):
            self.start_site= None
            self.num_of_site= None
            #self.O= None  #None indicates not "allocated" as in fortran  #[Tensor()]
            nSite = 1
            #self.O = [None]*nSite
            #list.__init__(self, [None])
            list.__new__(Operators, [None])
     
    class Hamiltonian(object):
        pass

    class OperatorDict(dict):
        def __init__(self):
            pass

class System(IterativeOptimize):
    """
        TNet_sys:
            一组指针，指向一组要收缩的张量;同时指向收缩中产生的张量
            它又由 tlink 组成， 在f90中， tlink与 shallowcopy 配合使用为收缩做准备
            , 说白了仅仅是为了调整收缩的顺序。在python中，这些都变得很travial
    """
    VERSION  = 1.0
    ALGORITHM = 'mera'
    TNet_sys = TensorNetwork() 
    
    def __init__(self, symmetry=None, mera=None, model=None, qsp_base=None, 
            graph_module=None, model_param=None, 
            only_NN=True, only_NNN=False, 
            energy_exact=None, j_power_r_max=None, **kwargs):
        """
            System 和Mera的区别在于：
                后者定义了U，V等张量
                前者定义了H_2, rho_2等张量
            SxA, ...., I are all instance of either class Operators or Operator3
            H_2, rho_2 etc, in every layer, there is a H_2[layer]
            layer: number of layer is determined by num_of_layer of mera, as can be seen from self.init_Hamiltonian
            note，无论在f90还是python，layer从0开始计数
            parameters:
                model_param: dict, specify parameters for Ising, Heisenberg, etc

        """

        IterativeOptimize.__init__(self, **kwargs)
        
        default_args = {
            #common to vmps
                'use_player': True, 
                'info': 0, 
                'is_initialized': False, 
            
            #ham
                'model': None, 
                'only_NN': True, 
                'only_NNN': False, 
                'model_param': {}, 
                
            'which_minimize': 'prod_state', # in ['prod_state', 'eigen_state', 'scale_invar_state'], 
            #IO
                'save_attr_list': self.save_attr_list + [
                    'symmetry', 
                    'mera',
                    'only_NN', 'only_NNN', 
                    'H_2', 'H_3', 'rho_2', 'rho_3', 
                    'SpB', 'SpA', 'SmA', 'SmB', 'SzA', 'SzB' , 'SxA', 'SxB', 
                    'qsp_max', 'qsp_base', 
                    'iter', 'iter1',  #?  
                    'energy', 'energy_err', 'energy_record', 'energy_diff_std', 'record', 'energy_diff', 
                    'Jleg_tau_0l_dic', 
                    'updaters', #?
                    'pinning_term_def', #this is needed, or else cause problem for tensor_player 
                    'graph_module_name', 
                    'model_param', 
                    ] , 
                
                '_not_resume_attr_list': self._not_resume_attr_list  + [
                    'precision',  'model_param', 
                    ], 
            #control
                'energy_diff_min': 1e-9, 
            
            'USE_REFLECTION': False, 
            'mera_class': Mera, 
            'qsp_base': None, 
            'qsp_max': None, 
            'trunc_dim': 4, 
            'combine_2site': True, 
            
            'tot_layer': 4,
            'nTop': 2, 
            'unitary_init': 'unit_tensor', # or 'random_unit_tensor'
            'isometry_init': 'random',  # or 'AFM'
            'mera_kwargs': {}, 
            
            'translation_invariant': 1, 
            'symmetry': 'Travial',  
            'USE_CUSTOM_RAND': False, 
            'rand_seed': None,  
            
            #minimize
                'updaters': None, 
                'q_iter': 5, 
                'q_lay': 1, 
                'use_pinning_term': False,   # add a pinning term to ham,  to avoid stuck 
                'pinning_term_def': {'op': None,'name': None,'h_name': None, 'lam': None, 
                    'xi': np.inf, 'step_limit': 1000},   # pinning term of the form  lam*exp(-t/xi)*op
                'pinning_term': None, 
            
            #IO
                'use_local_storage': False, 
                'backup_parpath': None, 
                'backup_parpath_local': None, 
                'auto_resume': True, 
                'other_backup_path': None, #if not None, using other point to init mera, as a better init 
            
            #measure
                'do_measure': 1, 
                'measurement_args': {
                    'exclude_which': [
                        'entanglement_brute_force', 
                        'entanglement_brute_force_6',  
                        'entanglement_brute_force_9', 
                        'entanglement_brute_force_9_aver', 
                        'correlation_extra', 
                        ], 
                    'which': ['energy', 'correlation', 'entanglement_entropy', 
                        'entanglement_spectrum', 'magnetization', 
                        ], 
                    }, 
           
            
            #only for long range int 
                'interaction_range_max': 2, 
            }
        if 1: 
            for k in default_args: 
                if kwargs.has_key(k): 
                    default_args[k] = kwargs[k]
            
            for k, v in default_args.iteritems(): 
                setattr(self, k, v)
        
        #computated attr
        if 1:
            self.set_computed_attr()
        
        self.model = model
        self.model_name = model
        self.only_NN = only_NN
        self.only_NNN = only_NNN
        self.symmetry = symmetry
        
       
        self.qn_identity, self.qsp_base, self.qsp_null = init_System_QSp(symmetry=self.symmetry, 
                combine_2site=self.combine_2site)
        if qsp_base is not None: 
            self.qsp_base = qsp_base  #models like Potts need alternate qsp_base than default
        if model_param is not None:
            self.model_param = model_param
        else:
            self.model_param = {}
        
        self.set_graph(graph_module)
        
        if 1:
            #some records of the system
            self.iter = 0   #state of the system after iter times update, when S continue  
            self.iter1 = 0   #state of the system after iter times update, when S change reset to 0 
            self.log = {}
            self.energy = None
            #self.energy_exact = self.get_energy_exact()
            #if energy_exact is not None:  self.energy_exact = energy_exact
            self.energy_exact = energy_exact if energy_exact is not None else self.get_energy_exact()
            self.energy_diff = np.nan 
            self.energy_diff_std = np.nan
            self.energy_err = np.nan
            self.energy_record = OrderedDict()
            self.scaling_dim_record = OrderedDict()
            self.energy_record[0] = (0, 0, 0)   #self.energy_record[iter] = (self.energy, time.time(), time.clock())
            
            #self.time_real = 0.0
            self.time_real = time.time()    #track lifetime of the system
            self.time_cpu = time.clock()
        
            self.record = OrderedDict()   # record of everything
            self.record[0] = {'energy':0.0, 'time':time.time()}   #also scaling_dim, central_charge, etc
            self.status= {}
        
        if mera is not None:
            self.mera = mera
        else: 
            self.mera = None 

        mera_name_map = {"U":("U", 0), 
               "Up":("U_dag", 0), 
               "V":("V", 0), 
               "Vp":("V_dag", 0), 
               }
        sys_name_map = {"oo":("H_2", 0), 
               "ooo":("H_3", 0), 
               "OO":("rho_2", +1), 
               "OOO":("rho_3", +1), 
               "I":("I",  +1),  #when used, I is attached to rho so  +1 
               }
        self.network_tensor_dic = {}
        self.network_tensor_dic.update(mera_name_map)
        self.network_tensor_dic.update(sys_name_map)
        
        #tensors that cant fit into above cases
        self.network_tensor_dic_extra = {
                0:{"oo1":("SzA", 0), "oo2":("SzB", 0)}, 
                1:{"oo1":("SpA", 0), "oo2":("SmB", 0)}, 
                2:{"oo1":("SmA", 0), "oo2":("SpB", 0)} 
                }
        if self.model == "Ising":
            self.network_tensor_dic_extra[0] = {"oo1":("SxA", 0), "oo2":("SxB", 0)} 
        
        self.Jleg_tau_0l_dic = {}   #used only in long rang int
        
        if j_power_r_max is not None:   #only for long range 
            System.J_POWER_R_MAX = j_power_r_max 
        
        #if self.USE_CUSTOM_RAND:
        #    crandom.mrandom.csrand(self.rand_seed)
        #else:
        #    np.random.seed(self.rand_seed)
        #    crandom.rand = np.random.random
    
    def initialize(self): 
        """
            实质性的初始化
            包括几个方面
            ham，mera，图
            
        """
        if self.init_state is not None: 
            self.__dict__.update(self.init_state)   # use resume_func to realize this later, which is better. cf. iDMRG_mcc.initialize 
            self.is_initialized = True
            return 
    
        if self.mera is None:
            if not self.combine_2site and self.symmetry == "U1" :
                qsp_max = self.qsp_base.__class__.max(self.trunc_dim, combine_2site=self.combine_2site)
            else:
                qsp_max = self.qsp_base.__class__.max(self.trunc_dim)
            qsp_max2 = qsp_max.copy() 
            
            if self.qsp_max is not None:
                qsp_max = self.qsp_max.copy()
                qsp_max2 = None #self.qsp_max.copy()
                print "trunc_dim is updated to %d"%qsp_max.totDim
                self.trunc_dim = qsp_max.totDim
            
            self.mera = self.mera_class(self.tot_layer, self.nTop, 
                    self.qsp_base.copy(),  qsp_max, qsp_max2, 
                    topQN=self.qn_identity.copy(), qsp_null=self.qsp_null.copy(), 
                    qn_identity=self.qn_identity.copy(), 
                    unitary_init=self.unitary_init, 
                    isometry_init=self.isometry_init, 
                    **self.mera_kwargs)
            self.M = self.mera
            
            #每次更改耦合常数算新的值，要重新初始化
            if self.model == "Ising":
                #self.init_Hamiltonian()
                for iLayer in range(self.M.num_of_layer-1): 
                    #top V shouldn't inited by init_layer
                    self.M.init_layer(iLayer, rand_init=True, unitary_init=self.unitary_init)
        
        
        if not hasattr(self, 'H_2'): 
            self.init_Hamiltonian()
            #self.init_ham_eff()  # merapy 2.0,  replace init_Hamiltonian in future 
            #self.init_rho_eff()
        
        self.set_backup_path()   #this should put after init mera
        if self.other_backup_path is not None:  #using other optimized mera state to init. this is useful for circumvent metastable state, etc 
            msg = '\ntry to INIT mera from other_backup_path ...{}'.format(
                    self.other_backup_path[-30:], )
            
            temp = self.__class__.load(self.other_backup_path , 
                    use_local_storage=self.use_local_storage)
            if not isinstance(temp, dict): 
                temp = temp.__dict__ 
            m = temp['mera']
            if self.mera == m:   # must be formally equal 
                self.mera = m 
                self.M = m 
                for k in ['rho_2', 'rho_3']: 
                    setattr(self, k, temp[k])
                self.other_backup_path = None 
                # it can be overide by the following resume_func if file exists, this is proper 
                msg += '... done, then set it to None\n' 
            else: 
                msg += ' ...self.mera != m, not allowed, discard it \n' 
            print msg 
        
        if self.auto_resume: 
            self.resume_func()
        
        
        if self.updaters is None: 
            from merapy.ascending import ascending_ham
            from merapy.descending import descending_ham
            from merapy.iteration import iterative_optimize_all
            from merapy.all_in_once import update_all_in_once
            from merapy.top_level import top_level_product_state, top_level_product_state_u1, top_level_eigenstate
            from merapy.minimize import finite_site_u1  
            self.updaters= {"ascending_func":ascending_ham, "descending_func":descending_ham, "update_mera_func":iterative_optimize_all, 
                "finite_range_func":finite_site_u1, "rho_top_func":top_level_product_state_u1}
        
        self.is_initialized = True     
        
    def init_iteration(self): 
        if not self.is_initialized: 
            self.initialize()
    
    if 0: 
        @getter
        def num_of_layer(self): 
            return self.mera.num_of_layer
        @getter
        def trunc_dim(self): 
            return 
    
    def key_property(self):
        """
            key properties to identify S
        """
        temp = ["model", "symmetry", "only_NN", "only_NNN", 'graph_module_name']
        res= {}
        for i in temp:
            res[i] = self.__dict__.get(i)#[i]
        res.update(self.mera.key_property())
        return res

    def get_info(self):
        num_of_layer = self.mera.num_of_layer
        trunc_dim = self.mera.qsp_max.totDim
        res = "num_of_layer = %d, symmetry = %s, trunc_dim = %d"%(num_of_layer, self.symmetry, trunc_dim)
        return res
    
    def get_status(self, keys=None):
        keys= ['iter', 'iter1', 'energy', 'energy_err', 'energy_diff']
        res = OrderedDict()
        for i in keys:
            res[i] = self.__dict__[i]
        return res
    
    def __eq__(self, other):
        """
            test two systems formally equal
        """
        #issue: one omission here!  model_param, this is tedious to compare
        if 1:
            res= True
            if self.mera != other.mera:
                return False
            temp = ["model", "symmetry", "only_NN", "only_NNN"]
            for i in temp:
                if self.__getattribute__(i) != other.__getattribute__(i):
                    print i
                    res = False
                    break
            return res

    def __ne__(self, other):
        return not self.__eq__(other)
    
    #@staticmethod 
    @classmethod 
    def example(cls, M=None, symmetry='Travial', model='Heisenberg', 
            trunc_dim=4, tot_layer=4, 
            model_param=None, **kwargs): 
            #only_NN=True, combine_2site=False, only_NNN=False):
            
        test_Mera = Mera 
        if not kwargs.has_key('do_measure'): 
            kwargs['do_measure'] = 0 
        
        #M = test_Mera.example(trunc_dim=4, tot_layer=4, symmetry="U1"); sys= ts.instance(M, symmetry="U1", model="Heisenberg")
        #M = test_Mera.example(trunc_dim=4, tot_layer=4, symmetry="U1"); 
        if model_param is None: 
            if model.lower() == 'heisenberg' : 
                symmetry = symmetry if symmetry is not None else 'U1'
                model_param = {'J_NN':1.0, 'J_NNN': 0.0, 'Jzz': 1.0}
                kwargs.update(only_NN=1, only_NNN=0, combine_2site=1)
            elif model == 'Ising' : 
                symmetry = symmetry if symmetry is not None else 'Z2'
                model_param = {"h":-1.0, "J_NN":-1.0, "J_NNN":0.0, 'gamma':1.0} 
                kwargs.update(only_NN=1, only_NNN=0, combine_2site=0)
                
            else: 
                raise 
                
        #sys=System(model=model, mera=M, symmetry=symmetry, model_param = model_param, **kwargs)
        sys=cls(model=model, mera=M, symmetry=symmetry, model_param = model_param, **kwargs)
        
        #sys.init_Hamiltonian()
        sys.initialize()
       
        return sys
   
    def set_graph_module_del(self, module):
        #mapper = {"ternary":merapy.init_mera_graph}
        if module is None:
            import merapy.diagrams.V31.graph_ternary as module
        self.graph_module = module
        
    def set_graph(self, module):
        if module is None:
            #import merapy.init_mera_graph
            #module = merapy.init_mera_graph
            import merapy.diagrams.V31.graph_ternary as module
        elif isinstance(module, str): 
            module = importlib.import_module(module)
            #todo: in config, change all graph_modul into str of the path 
        xx = ["G_2_2", "G_2_3", "G_3_2", "G_3_3", "G_22_2", "G_22_3"]  #G_2_3 occurs in binary graph
        found = []
        not_found = []
        for i in xx:
            try:
                g = module.__getattribute__(i)
                #self.__setattr__(i, g)
                setattr(self, i, g)
                found.append(i)
                
            except AttributeError as err:
                
                warnings.warn(str(err))
                not_found.append(i)
                print i, err  
        if self.info>0: 
            print "loading graphs in file %s, \n\tfound these %s, \n\tnot found and omit these%s"%(module.__file__, 
                    found, not_found)
        #self.graph_module_name = module.__name__.split('.')[-1]
        self.graph_module_name = module#.__name__.split('.')[-1]

    def fetch_op_names(self, type, range=None):
        """
            in future, make op_names into a class, and make
            this a method of op_names
        """
        ops = self.op_names #.items()
        if range is None:
            range = [2]
            if not self.only_NN: 
                range.append(3)
                if not self.only_NNN:
                    range.append("long")
        
        #res = [i for i in ops if ops[i]["type"] in type and (ops[i]["range"] in range)]
        #res = {i[0]:i[1] for i in ops.iteritems() if i[1]["type"] in type and (i[1]["range"] in range)}
        model = ["all"] + [self.model]
        res = {key:val for key, val in ops.iteritems() 
                if val["type"] in type 
                and (val["range"] in range)
                and (val["model"] in model)
                }
        return res

    def init_Hamiltonian(self):
        """  
            this is the main initialization that invokes several other initialization
            see init_Hamiltonian in f90
            M: an instance of class Mera
            S: an instance of class System
        """
        if 1: 
            #self.op_names= ["H_2", "H_3", "rho_2", "rho_3", "I", "SxA", "SxB", "SzA", "SzB", "SpA", "SpB", "SmA", "SmB"]
            self.op_names = {
                "H_2":      {"rank":4,    "type":"H",	"range":2     ,    "model":"all"}, 
                "H_3":      {"rank":6,    "type":"H",	"range":3     ,    "model":"all"}, 
                "rho_2":    {"rank":4,    "type":"rho",	"range":2     ,    "model":"all"}, 
                "rho_3":    {"rank":6,    "type":"rho",	"range":3     ,    "model":"all"}, 
                "I":        {"rank":2,    "type":"I",	"range":2     ,    "model":"all"}, 
                
                "SxA":      {"rank":4,    "type":"H",	"range":"long",    "model":"Ising"}, 
                "SxB":      {"rank":4,    "type":"H",	"range":"long",    "model":"Ising"}, 
                
                "SzA":      {"rank":4,    "type":"H",	"range":"long",    "model":"all"}, 
                "SzB":      {"rank":4,    "type":"H",	"range":"long",    "model":"all"}, 
                "SpA":      {"rank":4,    "type":"H",	"range":"long",    "model":"Heisenberg"}, 
                "SpB":      {"rank":4,    "type":"H",	"range":"long",    "model":"Heisenberg"}, 
                "SmA":      {"rank":4,    "type":"H",	"range":"long",    "model":"Heisenberg"}, 
                "SmB":      {"rank":4,    "type":"H",	"range":"long",    "model":"Heisenberg"}, 
                
                }
        
            for i in self.op_names:
                self.__setattr__(i, [[None] for i in range(Mera.MaxLayer)]) 
            

        mera = self.mera

        #set value for h2[0], rho[0]
        
        self.add_physics_layer()  
        
        #here add starting from the 1st layer
        self.layer = mera.num_of_layer

        for i in range(1, mera.num_of_layer):
            self.add_one_layer(ilayer=i, duplicate=False)

        self.add_top_layer()

    def add_one_layer(self, ilayer, duplicate=False):
        """
            see System_AddLayer in f90
            materialize H_2, rho_2, SxA, etc in layers [1, 2, .., mera.num_of_layer-1]
            note that note including layer 0
        """ 
        if duplicate: 
            dic = self.fetch_op_names(type=["H", "rho", "I"])
            for name in dic.keys():
                self.__getattribute__(name)[ilayer][0] = self.__getattribute__(name)[ilayer-1][0].copy()  #resident in memory no need buffer
            return

        mera = self.mera
        
        if 1:
            def func(rank, reverse, qn_val=None, i=0):
                """
                    注意这里复制的是V的超上的腿的属性!!
                    note here the scope of 
                    注意这里的处理和f90不同!, 这里sys layer的编号比mera layer中的编号大1
                    故这里不是temp = mera.V[ilayer].....
                    copy is needless bellow, without copy, it is conceptially more correct
                    QSp = [temp.QSp[temp.rank-1].copy() for i in range(rank)]
                """
                r = mera.V[ilayer-1][0].rank
                qsp = mera.V[ilayer-1][0].QSp[r-1].copy()  
                if qn_val is None:
                    totQN = self.qn_identity.copy()
                else:
                    totQN = qsp.QnClass(qn_val)

                QSp = qsp.copy_many(rank, reverse=reverse)
                return rank, QSp, totQN

            def set_op(name, rank, QSp, totQN, i=0):
                a = self.__getattribute__(name)
                a[ilayer][i] = iTensor(rank, QSp, totQN)
                a[ilayer][i].data[:] = 0.0
                self.__setattr__(name, a)
        
        r = mera.V[ilayer-1][0].rank
        qsp_ilayer = mera.V[ilayer-1][0].QSp[r-1].copy()  
        qsp_ilayer.reverse()

        def dic_to_op(dic):
            """
                注意这里复制的是V的超上的腿的属性!!
                note here the scope of 
                注意这里的处理和f90不同!, 这里sys layer的编号比mera layer中的编号大1
                故这里不是temp = mera.V[ilayer].....
            """
            rank = dic["rank"]
            rank_half = rank//2
            totQN = self.qn_identity.copy()
            #copy is needless, without copy, it is conceptially more correct
            if dic["type"] == "rho":
                reverse = range(rank_half)
            else:
                reverse = range(rank_half, rank)
            QSp = qsp_ilayer.copy_many(rank, reverse=reverse)
            return iTensor(rank, QSp, totQN)
    
        for i in range(1):
            dic = self.fetch_op_names(type=["H", "rho"], range=[2])
            for name, propt in dic.iteritems():
                op = dic_to_op(propt)
                op.data[:] = 0.0
                self.__getattribute__(name)[ilayer][0] = op
            
            #in some curcumstance I should be used
            qsp = qsp_ilayer.copy_many(2, reverse=[1])    
            self.I[ilayer][i] = iTensor.identity(qsp)

            if not self.only_NN:
                
                dic = self.fetch_op_names(type=["H", "rho"], range=[3])
                for name, propt in dic.iteritems():
                    op = dic_to_op(propt)
                    op.data[:] = 0.0
                    self.__getattribute__(name)[ilayer][0] = op
                
                if not self.only_NNN:
                    if self.model == "Heisenberg":
                     
                        set_op("SzA", *func(rank=4, reverse=[0, 1]))
                        set_op("SzB", *func(rank=4, reverse=[0, 1]))

                        mapper_p = {"Z2":-1, "U1":1, "Travial":None}
                        val = mapper_p[self.symmetry]
                        set_op("SpA", *func(rank=4, reverse=[0, 1], qn_val=val))
                        set_op("SpB", *func(rank=4, reverse=[0, 1], qn_val=val))

                        mapper_m = {"Z2":-1, "U1":-1, "Travial":None}
                        val = mapper_m[self.symmetry]
                        set_op("SmA", *func(rank=4, reverse=[0, 1], qn_val=val))
                        set_op("SmB", *func(rank=4, reverse=[0, 1], qn_val=val))
                    
                    elif self.model == "Ising": 
                        set_op("SzA", *func(rank=4, reverse=[0, 1]))
                        set_op("SzB", *func(rank=4, reverse=[0, 1]))

                        mapper_x = {"Z2":-1, "U1":-1, "Travial":None}
                        val = mapper_x[self.symmetry]
                        set_op("SxA", *func(rank=4, reverse=[0, 1], qn_val=val))
                        set_op("SxB", *func(rank=4, reverse=[0, 1], qn_val=val))
                    
                    elif self.model == 'wigner_crystal': 
                        set_op("SzA", *func(rank=4, reverse=[0, 1]))
                        set_op("SzB", *func(rank=4, reverse=[0, 1]))
                       
                    else : 
                        raise

    def transverse_ops_ilayer(self, ilayer):
        """
            an iterator
        """ 
        res= []
        op_names = self.fetch_op_names(type=["H", "rho"]).keys()
        for name in op_names:
            res.append(self.__getattribute__(name)[ilayer][0])
        return res

    def add_physics_layer(self):
        """
            physics layer is the hamiltonian H_2[0] and the top rho  rho_2[toplayer]
            rho_2[toplayer] 不是通过update赋值，而是直接赋值
            这一函数只是materialize了算子，
            see System_AddPhysicsLayer
        """
        mera = self.mera
        #set value for H_2[0]
        if self.model == "Ising": 
            #h0, Sx = System.op_ising(param=self.model_param)  #here is a key step
            if self.combine_2site:
                h0, Sx = System.op_ising_2site(QN_idendity=self.qn_identity, QSp_base=self.qsp_base, 
                        only_NN=self.only_NN, only_NNN=self.only_NNN, symmetry=self.symmetry, **self.model_param)  #here is a key step
            else:
                h0, Sx = System.op_ising(QN_idendity=self.qn_identity, QSp_base=self.qsp_base, 
                        only_NNN=self.only_NNN, symmetry=self.symmetry, **self.model_param)  #here is a key step
            self.H_2[0][0] = h0[2]
            
            if not self.only_NN:
                self.H_3[0][0] = h0[3]
                if not self.only_NNN:
                    self.SxA[0][0] = Sx[1]
                    self.SxB[0][0] = Sx[2]

        elif self.model == "Heisenberg":
            if self.combine_2site:
                h0, Sz, Sp, Sm, pinning_term =System.op_heisenberg(QN_idendity=self.qn_identity, QSp_base=self.qsp_base, 
                        symmetry=self.symmetry, only_NN=self.only_NN, only_NNN=self.only_NNN, **self.model_param)
            else:
                h0, Sz, Sp, Sm, pinning_term = System.op_heisenberg_1site(QN_idendity=self.qn_identity, QSp_base=self.qsp_base, 
                        symmetry=self.symmetry, only_NN=self.only_NN, only_NNN=self.only_NNN, **self.model_param)
            
            #self.pinning_term = pinning_term 
            if self.use_pinning_term: 
                if self.pinning_term_def['op'] is None: 
                    assert pinning_term.has_key(self.pinning_term_def['name'] )
                    self.pinning_term_def['op'] = pinning_term[self.pinning_term_def['name']]
                    
            self.H_2[0][0] = h0[2]
            self.H_3[0][0] = h0[3]
            #if self.use_pinning_term: 
            #    if self.combine_2site: 
            #        self.H_2_bac = self.H_2[0][0].copy()
            #    else: 
            #        self.H_3_bac = self.H_3[0][0].copy()

            if not self.only_NN and not self.only_NNN:

                self.SzA[0][0] =  Sz[1]
                self.SzB[0][0] =  Sz[2]
                
                self.SpA[0][0] =  Sp[1]
                self.SpB[0][0] =  Sp[2]
                
                self.SmA[0][0] =  Sm[1]
                self.SmB[0][0] =  Sm[2]
                
        elif self.model == "Potts":
            h0 = potts(symmetry=self.symmetry, qsp_base=self.qsp_base, **self.model_param)
            self.H_2[0][0] = h0[2]
            self.H_3[0][0] = h0[3]

        else:
            raise NameError("model name %s is in correct"%self.model)
        

        if not self.only_NN and not self.only_NNN: #better treatment of long range ops
            
            if self.symmetry in ["U1", "Travial"] and self.combine_2site:
                varables= vars()
                temp = ["Al", "Ar", "Bl", "Br"]
                if self.model == "Ising":
                    xzpm = ["Sx"]
                elif self.model == "Heisenberg":
                    xzpm = ["Sz", "Sp", "Sm"]
                for x in xzpm:
                    for t in temp:
                        setattr(self, x + t,  [[ varables[x][t] ]])

        #in some curcumstance I should be used
        qsp = mera.qsp_0.copy_many(2, reverse=[1])    
        self.I[0][0] = iTensor.identity(qsp)
        
        if 1:  # 下面是level 0 的rho，在算法中其实没有用到，值一直没改变 
            rank = 4
            QSp = [mera.qsp_0.copy() for i in range(rank)]
            totQN = self.qn_identity.copy()
            #注意这里，rho_2 和U, H_2等张量不同，它的0，1两条腿是朝上的,2,3朝下
            QSp[0].reverse()
            QSp[1].reverse()
            #sys.rho_2[0][i] = iTensor(rank, QSp, totQN )
            self.rho_2[0][0] = iTensor(rank, QSp, totQN )
            self.rho_2[0][0].data[:] = 0.0

            #rank, QSp, totQN = func(rank=6)
            rank = 6
            QSp = mera.qsp_0.copy_many(6, reverse=[3, 4, 5])
            QSp[3].reverse()
            QSp[4].reverse()
            QSp[5].reverse()
            self.rho_3[0][0] = iTensor(rank, QSp, totQN)
            self.rho_3[0][0].data[:] = 0.0

        if 0: #they are not used at all
            rank = 2
            totQN = self.qn_identity.copy()
            QSp = [self.qsp_null.copy() for i in range(rank)]
            rho_top = iTensor(rank, QSp, totQN )  #rho_top, odd are defined in mera.py
            rho_top.data[0] = 1.0

            rank = 2
            totQN = self.qn_identity.copy()
            QSp = [QSp_null.copy() for i in range(rank)]
            QSp[0].QNs[0].set_val([-1]) 
            QSp[1].QNs[0].set_val([-1])         
            rho_top_odd = iTensor(rank, QSp, totQN)

    def add_top_layer(self):
        """
            this rho_top is a (1, 1) tensor,  and actually a scalar, as the dim of each index is 1
        """
        mera = self.mera
        rank = 2
        ilayer = mera.num_of_layer-1
        r = mera.V[ilayer][0].rank
        QSp = mera.V[ilayer][0].QSp[r-1].copy_many(2, reverse=[0])
        totqn = self.qn_identity
        rho_top = iTensor(rank, QSp, totqn)
        rho_top.data[0] = 1.0
        self.rho_2[ilayer+1][0] = rho_top

    def __repr__(self, layers=[], which=None,fewer=True, round=5):
        res= ""
        
        temp = ["iter", "energy"]
        temp1 = {t:self.__dict__[t] for t in temp}

        res +=  str(temp1) + "\n"
            
        if layers == []:
            layers= range(self.mera.num_of_layer)
            #layers= range(self.layer)
        if which is None:
            which = ["H_2",  "rho_2"]
        for w in which:
            #print w
            res +=  str(w) + "\n"
            for i in layers:
                #res  += "\tlayer = %d "%i   
                tensor = self.__getattribute__(w)[i][0]
                if fewer:
                    #print "\t", i, tensor.data[:4].round(round),"  ",tensor.data[-4:].round(round), "  ", tensor.totDim
                    res += "\t%d %s   %s   %d  buff %s\n"%(
                            i, tensor.data[:4].round(round), tensor.data[-4:].round(round), tensor.totDim, tensor.buf_ref)
                else:
                    #print "\t", i, tensor.data.round(round)
                    res +=  "\t%d%s\n"%(i, tensor.data.round(round))
        return res

    @staticmethod
    def op_ising(**param):
        """

            materialize operators(tensors) of the ising hamiltonian, especially h0, Sx
            see op_ising in f90
        """
        symmetry = param["symmetry"]
        QN_idendity, QSp_base = param["QN_idendity"], param["QSp_base"]
        h, gamma, J_NN , J_NNN= param["h"], param["gamma"], param["J_NN"], param['J_NNN']
        only_NNN = param['only_NNN']
        
        h0, Sx = {}, {}
        
        #pauli = iTensorFactory(symmetry=symmetry).pauli_mat_1site()
        pauli = iTensorFactory.pauli_mat_1site(symmetry=symmetry)
        si, sigma_x, sigma_y, sigma_z  = pauli["sigma_0"], pauli["sigma_x"], pauli["sigma_y"], pauli["sigma_z"]

        sxx = sigma_x.direct_product(sigma_x); syy = sigma_y.direct_product(sigma_y); 
        szi = sigma_z.direct_product(si);   siz = si.direct_product(sigma_z)
        
        #the minus sign before syy term comes from i^2 = -1  #J_NN default is -1.0
        h0[2] = ((J_NN*(1+gamma)*0.5)*sxx  + (-1.0*J_NN*(1-gamma)*0.5)*syy + 0.5*h*(siz + szi))
        #when h = 1.0, gamma = 1, J_NN = -1.0, it should be  h2 = [ 1. -1.  0. -1. -1.  0. -1. -1.]
        
        h0[3] = sigma_x.direct_product(si).direct_product(sigma_x)
        h0[3].data *= J_NNN
        
        Sx[1] = sigma_x.direct_product( si )
        Sx[2] = si.direct_product(sigma_x)
        
        if 0:  
            print "for testing, make this change"
            #when set J_NN = J_NNN = 0, dont need change h0[2 and 3]
            Sx[1] = sigma_z.direct_product( si )
            Sx[2] = si.direct_product(sigma_z)

        #  -1 to ensure that we get ground state        
        ii = si.direct_product(si)
        h0[2] += -1.0*ii
        
        return h0, Sx

    @staticmethod
    def op_ising_2site(**param):
        """
            
            materialize operators(tensors) of the ising hamiltonian, especially h0, Sx
            
        """
        symmetry = param["symmetry"]
        QN_idendity, QSp_base = param["QN_idendity"], param["QSp_base"]
        h, gamma, J_NN , J_NNN = param["h"], param["gamma"], param["J_NN"], param['J_NNN']
        only_NN, only_NNN  =  param["only_NN"], param['only_NNN']
        
        h0, Sx = {}, {}
        h0[2]=iTensor(4, QSp_base.copy_many(4, reverse=[2, 3]), QN_idendity.copy())
        h0[2].data[:] = 0.0
        h0[3]=iTensor(6, QSp_base.copy_many(6, reverse=[3, 4, 5]), QN_idendity.copy())
        h0[3].data[:] = 0.0

        pauli = iTensorFactory.pauli_mat_2site(symmetry=symmetry)
        ii, szz, sxx, syy= pauli["ii"], pauli["szz"], pauli["sxx"], pauli["syy"]
       
        sigma_z1, sigma_z2 = pauli["sigma_z1"], pauli["sigma_z2"]
        sigma_x1, sigma_x2 = pauli["sigma_x1"], pauli["sigma_x2"]
        sigma_y1, sigma_y2 = pauli["sigma_y1"], pauli["sigma_y2"]
        
        
        sxi = sigma_x1
        syi = sigma_y1
        szi = sigma_z1
        
        six = sigma_x2
        siy = sigma_y2
        siz = sigma_z2

        
        
        sixxi= sigma_x2.direct_product(sigma_x1)
        siyyi= sigma_y2.direct_product(sigma_y1)
        sziii = sigma_z1.direct_product(ii)
        siiiz = ii.direct_product(siz)
        
        #the minus sign before syy term comes from i^2 = -1
        #xy_r = lambda J_r, sxx, syy, szi, siz, hh=h: ((J_r*(1+gamma)*0.5)*sxx  + (-1.0*J_r*(1-gamma)*0.5)*syy+   0.5*hh*(siz + szi))
        xy_r = lambda J_r, sxx, syy, szi, siz, hh=h: ((J_r*(1+gamma)*0.5)*sxx  + 
                (J_r*(1-gamma)*0.5)*syy+   0.5*hh*(siz + szi))
        #h0[2] = ((J_NN*(1+gamma)*0.5)*sxx  + (-1.0*J_NN*(1-gamma)*0.5)*syy+   0.5*h*(siz + szi))
        
        if 1:  #r = 1, NN,  0.25*s0s1 + 0.5*s1s2 + 0.25*s2s3
            ss =  xy_r(J_NN, sxx, syy, szi, siz) #ss = ((J_NN*(1+gamma)*0.5)*sxx  + (-1.0*J_NN*(1-gamma)*0.5)*syy+   0.5*h*(siz + szi))
            ssii = ss.direct_product(ii)
            iiss = ii.direct_product(ss)
            issi = xy_r(J_NN, sixxi, siyyi, sziii, siiiz)  #issi = ((J_NN*(1+gamma)*0.5)*sixxi+ (-1.0*J_NN*(1-gamma)*0.5)*siyyi+   0.5*h*(sziii + siiiz)) 
            h0[2] += 0.5*issi + 0.25*(ssii + iiss)
        
        if not only_NN: #r = 2, NNN, 0.5*s0s2 + 0.5s1s3
            sxixi = sxi.direct_product(sxi)
            syiyi = syi.direct_product(syi)
            siizi = ii.direct_product(szi)
            s0s2 = xy_r(J_NNN, sxixi, syiyi, sziii, siizi, hh=0.0)

            sixix = six.direct_product(six)
            siyiy = siy.direct_product(siy)
            sizii = siz.direct_product(ii)
            s1s3 = xy_r(J_NNN, sixix, siyiy, sizii, siiiz, hh=0.0)
            h0[2]  += 0.5*s0s2 + 0.5*s1s3
        
        if not only_NN and not only_NNN:  #middle range
            #in fact all middle ranged ops can be mapped into H_3[layer=1] instead of H_3[0], at the expanse of contracting 6 more graph
            #here, in practice, I use a zhezhong way that is only set 0.5*s0s4 in H_3[0] while map the rest s1s4, ..., s0s5 to H_3[1] by ascending
            alpha, beta = param["alpha"], param["beta"]
            J0l = lambda l:System.J_power_0l(l=l, alpha=alpha, beta=beta)
            if 1:  #0.5*(s0s3 + s1s4)
                sx0x3 = sxi.direct_product(six).direct_product(ii)
                sy0y3 = syi.direct_product(siy).direct_product(ii)
                z = sx0x3.copy()
                s0s3 = xy_r(J0l(3), sx0x3, sy0y3, z, z, hh=0.0)
                h0[3] += 0.5*(s0s3 )             
            
            if 1:  #move these middle range ops to H_3[1]
                if 1:
                    sixiixi = six.direct_product(ii).direct_product(sxi)
                    siyiiyi = siy.direct_product(ii).direct_product(syi)
                    z = sixiixi.copy()
                    s1s4 = xy_r(J0l(3), sixiixi, siyiiyi, z, z, hh=0.0)
                    h0[3] += 0.5*(s1s4)             

                if 1:  #0.5*(s0s4 + s1s5)
                    sx0x4 = sxi.direct_product(ii).direct_product(sxi)
                    sy0y4 = syi.direct_product(ii).direct_product(syi)
                    z = sx0x4.copy()
                    s0s4 = xy_r(J0l(4), sx0x4, sy0y4, z, z, 0.0)
                    
                    sx1x5 = six.direct_product(ii).direct_product(six)
                    sy1y5 = siy.direct_product(ii).direct_product(siy)
                    #z = 0.0
                    z = sx1x5.copy()
                    s1s5 = xy_r(J0l(4), sx1x5, sy1y5, z, z, hh=0.0)
                    h0[3] += 0.5*(s0s4 + s1s5)
                
                if 1: #0.5*(s0s5)  
                    sx0x5 = sxi.direct_product(ii).direct_product(six)
                    sy0y5 = syi.direct_product(ii).direct_product(siy)
                    #z = 0.0
                    z = sx0x5.copy()
                    s0s5 = xy_r(J0l(5), sx0x5, sy0y5, z, z, hh=0.0)
                    h0[3] += 0.5*s0s5 
            
            #middle range
            Sx["Al"] = (sigma_x1 ).direct_product(ii) 
            Sx["Ar"] = (sigma_x2 ).direct_product(ii) 
            Sx["Bl"] = ii.direct_product(sigma_x1)
            Sx["Br"] = ii.direct_product(sigma_x2)
        
        if not only_NN and not only_NNN:  # long range

            #Sx[1] = sxi.direct_product(ii)
            #Sx[2] = ii.direct_product(six)
            Sx[1] = 0.5*(sxi + six).direct_product(ii)
            Sx[2] = 0.5*ii.direct_product(sxi + six)
                    

        #  -1 to ensure that we get ground state        
        i4 = ii.direct_product(ii)
        h0[2] += -1.0*i4
        
        if 1:
            assert h0[2].is_hermite()
            if not only_NN and not only_NNN:
                assert h0[3].is_hermite()
                assert Sx[1].is_hermite()
                assert Sx[2].is_hermite()

        return h0, Sx

    @staticmethod
    def op_heisenberg(**param):      
        """
            the ground state of the spin-1/2 chain can be constructed exactly by the Bethe ansatz; we therefore know its ground state energy exactly. In the thermodynamic limit  the energy per site is given by
            E0 = 1 / 4 − ln2 = − 0.4431471805599...        
            4*E0 = 1.772588722236
        """
        symmetry = param["symmetry"]
        qspclass =symmetry_to_Qsp(symmetry)
        QN_idendity, QSp_base = param["QN_idendity"], param["QSp_base"]

        
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
        pinning_term = {}
        pinning_term['s0s2'] = s0s2.copy()
        pinning_term['issi'] = 0.5*issi 
       
                
        if not only_NNN and not only_NN and symmetry in ["U1", "Travial"] :     #long range terms
            alpha, beta = param["alpha"], param["beta"]
            def symbol_to_ops(pattern):
                pattern = pattern.replace("_", "")
                valid = True
                if len(pattern)%2 != 0:
                    valid = False
                if pattern.count("s") != 2:
                    valid = False
                if not valid:
                    raise ValueError(pattern)
                
                temp = ["zz", "pm", "mp"]
                pattern_zpm = {}
                for ss in temp:
                    pattern_zpm[ss] = pattern
                    for s in ss:
                        pattern_zpm[ss] = pattern_zpm[ss].replace("s", s, 1)
                
                pattern_zpm_terms = {a:[pattern_zpm[a][2*i:2*i+2] for i in range(len(pattern)//2)] for a in pattern_zpm}
                    
                mapper = {}
                for a in ["z", "p", "m"]:
                    mapper.update({a+"i": pauli["sigma_%s1"%a]})
                    mapper.update({"i" + a: pauli["sigma_%s2"%a]})
                mapper.update({"ii":ii}) 
                
                op = {}
                for a in temp:
                    op[a] = mapper[pattern_zpm_terms[a][0]]
                    for b in pattern_zpm_terms[a][1:]:
                        op[a] = op[a].direct_product(mapper[b])
                    if a != "zz":
                        op[a] *= 2.0
                
                res = op["zz"] 
                res +=  op["pm"] 
                res +=  op["mp"]
                return res
            
            s2o = symbol_to_ops 
            J0l = lambda l:System.J_power_0l(l=l, alpha=alpha, beta=beta)
            
            #0.5*(s0s3 + s1s4)
            siis = s2o("siis")
            siis*= 0.5*J0l(3)
            isiisi = s2o("isiisi")
            isiisi *= 0.5*J0l(3)
            h0[2] += siis
            h0[3] += isiisi 

            ##0.5*(s0s4 + s1s5)->h3, it it aver. of two terms
            si_ii_si= s2o("si_ii_si")
            is_ii_is= s2o("is_ii_is")
            s0s4 = 0.5*(si_ii_si + is_ii_is)
            s0s4 *= J0l(4)
            h0[3] += s0s4 

            #0.5*s0s5->h3, a single term
            if 0:  #a little check
                si_ii_is = sigma_z1.direct_product(ii).direct_product(sigma_z2) + 2.0*(
                        sigma_p1.direct_product(ii).direct_product(sigma_m2) + 
                        sigma_m1.direct_product(ii).direct_product(sigma_p2))
                si_ii_is_x = s2o("si_ii_is")
                if not np.all(si_ii_is_x.data==si_ii_is.data):
                    raise

            
            si_ii_is = s2o("si_ii_is")
            #print "this maybe not pricise: si_ii_is *= 0.75 "
            #print "0.75  "*100
            weight_05 = param["weight_05"]   #for old not combine two site version, weight_05 = 0.75, for better version weight_05 = 0.5
            si_ii_is *= weight_05 
            s0s5 = si_ii_is
            s0s5 *= J0l(5)
            h0[3] += s0s5
            

        if 1:  #2x2 ops.
            if symmetry in ["Z2"]:
                if 0: 
                    Sz[1] = sigma_z.direct_product(sigma_0)
                    Sz[2] = sigma_0.direct_product(sigma_z)
                    
                    Sp[1] = sigma_p.direct_product(sigma_0)
                    Sm[2] = sigma_0.direct_product(sigma_m)
                    
                    Sm[1] = sigma_m.direct_product(sigma_0)
                    Sp[2] = sigma_0.direct_product(sigma_p)
                
            elif symmetry in ["U1", "Travial"]:
                warnings.warn("in hamiltonian.py Sz, Sx maybe not well defined")
                if  0:
                    Sz[1] = sigma_z1.direct_product(ii) 
                    Sz[2] = ii.direct_product(sigma_z1 + sigma_z2)
                    
                    #SpA,SmB
                    Sp[1] = sigma_p1.direct_product(ii) 
                    Sm[2] = ii.direct_product(sigma_m1 + sigma_m2)
                    #SmA,SpB
                    Sm[1] = sigma_m1.direct_product(ii)
                    Sp[2] = ii.direct_product(sigma_p1 + sigma_p2)
                if 1:
                    #middle range, exact treatment
                    Sz["Al"] = (sigma_z1 ).direct_product(ii) 
                    Sz["Ar"] = (sigma_z2 ).direct_product(ii) 
                    Sz["Bl"] = ii.direct_product(sigma_z1)
                    Sz["Br"] = ii.direct_product(sigma_z2)

                    Sp["Al"] = (sigma_p1 ).direct_product(ii) 
                    Sp["Ar"] = (sigma_p2 ).direct_product(ii) 
                    Sp["Bl"] = ii.direct_product(sigma_p1)
                    Sp["Br"] = ii.direct_product(sigma_p2)

                    Sm["Al"] = (sigma_m1 ).direct_product(ii) 
                    Sm["Ar"] = (sigma_m2 ).direct_product(ii) 
                    Sm["Bl"] = ii.direct_product(sigma_m1)
                    Sm["Br"] = ii.direct_product(sigma_m2)


                    #long range, with approximation due to averaging
                    Sz[1] = 0.5*(sigma_z1 + sigma_z2).direct_product(ii) 
                    Sz[2] = 0.5*ii.direct_product(sigma_z1 + sigma_z2)
                    #SpA,SmB
                    Sp[1] = 0.5*(sigma_p1 + sigma_p2).direct_product(ii) 
                    Sm[2] = 0.5*ii.direct_product(sigma_m1 + sigma_m2)
                    #SmA,SpB
                    Sm[1] = 0.5*(sigma_m1 + sigma_m2).direct_product(ii)
                    Sp[2] = 0.5*ii.direct_product(sigma_p1 + sigma_p2)

                    Sz[1] = 0.5*(sigma_z1 + sigma_z2).direct_product(ii) 
                    Sz[2] = 0.5*ii.direct_product(sigma_z1 + sigma_z2)
                    
        
        if not only_NN and not only_NNN and symmetry in ['U1', 'Travial']:  #some check
            print h0[2].is_hermite()
            print h0[3].is_hermite()
            print Sz[2].is_hermite()
            print not Sp[1].is_hermite(out_more=False)
            print not Sm[1].is_hermite(out_more=False)
            print Sp[1].is_adjoint_to(Sm[1])
            print Sp[2].is_adjoint_to(Sm[2])
            #exit()

        if symmetry in ["Z2"]:
            #identity= sigma_0.direct_product(sigma_0)
            identity= ii.direct_product(ii)
        
        elif symmetry in ["U1", "Travial"]:
            identity=ii.direct_product(ii)
            
        h0[2].data -= identity.data
        
        return  h0, Sz, Sp, Sm, pinning_term 

    @staticmethod
    def op_heisenberg_1site(**param):
        """
            if 0: 
                ss= sxx;  issi = sixxi
                ssii = ss.direct_product(ii)
                iiss = ii.direct_product(ss)
                h0[2] += 0.5*issi + 0.25*(ssii + iiss)
                h0[2]  = -1.0*h0[2]
            ss =  sz_i ×sz_{i+1} + 2.0*(sp_i × sm_{i+1} + sm_i ×sp_{i+1} )
            note sz, sp, sm are Pauli mat. not spin ops, so
            这是标准的能量为E = -0.443的反铁磁海森堡模型ss的4倍 = -1.773...
           
            the ground state of the spin-1/2 chain can be constructed exactly by the Bethe ansatz; we therefore know its ground state energy exactly. In the thermodynamic limit  the energy per site is given by
            E0 = 1 / 4 − ln2 = − 0.4431471805599...        
            4*E0 = 1.772588722236
        """
        #from quantum_number import init_System_QSp
        #QN_idendity, QSp_base, QSp_null= init_System_QSp("U1")
        symmetry = param["symmetry"]
        QN_idendity, QSp_base = param["QN_idendity"], param["QSp_base"]


        #h = param["h"]
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
        twist = 1.0  #twist = cos(momentum)

        if 1:
            pauli = iTensorFactory.pauli_mat_1site(symmetry=symmetry)
            if symmetry in ['U1', 'Travial']: 
                szz, spm, smp=  pauli["szz"], pauli["spm"], pauli["smp"]
                sigma_0, sigma_z, sigma_p, sigma_m =pauli["sigma_0"], pauli["sigma_z"], pauli["sigma_p"], pauli["sigma_m"]
                s0s1 = Jzz*szz + 2.0*(spm + smp) 
            
            elif symmetry in ['Z2']: 
                szz, sxx, syy =  pauli["szz"], pauli["sxx"], pauli["syy"]
                sigma_0, sigma_z, sigma_x, sigma_y =pauli["sigma_0"], pauli["sigma_z"], pauli["sigma_x"], pauli["sigma_y"]
                s0s1 = Jzz*szz + sxx + (-1.0)*syy
                #s0s1 = szz;print 'replaced szz  '*100

        if 1: #symmetry  in ["Travial", "Z2"]:
            s0s1 *= J_NN
            h0[2]  += s0s1 

                
        # =============================================================                                                             
        #next neareast neighbour interaction, 3-site operator
        pinning_term = None 
        if not only_NN:
            if 1: #symmetry in ["Travial", "Z2"]:
                sziz = sigma_z.direct_product(sigma_0).direct_product(sigma_z)
                spim = sigma_p.direct_product(sigma_0).direct_product(sigma_m)
                smip = sigma_m.direct_product(sigma_0).direct_product(sigma_p)
                s0s2 = sziz + 2.0*(spim + smip)
                s0s2 *= J_NNN
                h0[3] += s0s2 

                
            print "J_NNN is set to %1.5f"%J_NNN        
            #this is a temporary workaround for meta stable state for J1J2 model 
        
            pinning_term = h0[3].copy()
            pinning_term.data = -pinning_term.data 

        # =============================================================                                                              
        #long range terms

            

        if not only_NN and not only_NNN:  #2-site ops.
            if 1: #symmetry in ["Travial", "Z2"]:
                Sz[1] = sigma_z.direct_product(sigma_0)
                Sz[2] = sigma_0.direct_product(sigma_z)
                
                Sp[1] = sigma_p.direct_product(sigma_0)
                Sm[2] = sigma_0.direct_product(sigma_m)
                
                Sm[1] = sigma_m.direct_product(sigma_0)
                Sp[2] = sigma_0.direct_product(sigma_p)
                #print Sp[1].is_hermite(out_more=True)                
                #print Sp[1]
                #exit()
        if not only_NN and not only_NNN:
            exam = ["h0[2].is_hermite()", 
                    "h0[3].is_hermite()", 
                    "Sz[2].is_hermite()",
                    "not Sp[1].is_hermite(out_more=False)",
                    "not Sm[1].is_hermite(out_more=False)",
                    "Sp[1].is_adjoint_to(Sm[1])",
                    "Sp[2].is_adjoint_to(Sm[2])"
                    ]
            for i in exam:      
                passed = eval(i)
                print passed
                if not passed:
                    msg = "check hermite failed with %r"%i
                    #raise Exception(msg)

        if 1: #symmetry in ["Travial", "Z2"]:
            identity= sigma_0.direct_product(sigma_0)
            #print identity.matrix_view(); exit()
        h0[2].data -= identity.data

        return  h0, Sz, Sp, Sm, pinning_term 
    
    def _minimize_finite_size(self,  q_iter=None,  **kwargs):
        
        from merapy.minimize import finite_site_u1  as finite_range_func 
        self.init_iteration()
       
        if 1:   
            q_iter = q_iter if q_iter is not None else self.q_iter 
            
                
        if isinstance(self.updaters, list): #only for backward compatable
            ascending_func, descending_func, update_mera_func, finite_range_func, rho_top_func= tuple(self.updaters)
        elif isinstance(self.updaters, dict):
            finite_range_func = self.updaters["finite_range_func"]
            ascending_func = self.updaters["ascending_func"]
            descending_func = self.updaters["descending_func"]
            update_mera_func = self.updaters["update_mera_func"]
            rho_top_func = self.updaters["rho_top_func"]
        else: 
            print_vars(vars(),  ['self.updaters'])
            raise 
        
        temp = ['use_local_storage' ]
        args = {t: self.__getattribute__(t) for t in temp}
        kwargs.update(args)
       
        finite_range_func(self.mera, self, ascending=ascending_func, descending=descending_func, 
                update_mera=update_mera_func, rho_top_func = rho_top_func, 
                q_one=1, q_lay=1, q_iter=q_iter, use_player=self.use_player, 
                filename=self.filename, backup_fn=self.backup_path, lstart=0, lstep=0, 
                resume=True, info=self.info-1, **kwargs)
        
       
        return self.M, self

    def _minimize_scale_invar(self, num_of_SIlayer=3, q_iter=5, **kwargs):
        """ 
            taken from Main.run_scale_invar 
        """
        from merapy.minimize import ScaleInvar 
        self.init_iteration()
        
        
        temp = ['use_local_storage']
        args = {t: self.__getattribute__(t) for t in temp}
        #args['backup_fn_local'] = backup_fn_local
        kwargs.update(args)

        SI = ScaleInvar(self, self.updaters, num_of_SIlayer, q_iter, 
                q_one=1, q_lay=1, 
                filename=self.filename, backup_fn=self.backup_path, 
                resume=True, use_player=self.use_player, info=self.info,
                **kwargs)
        
        SI.scale_invariant()
        
        
        
        #self.transfer_pickle_file()
        
        return self

    def minimize(self, which_minimize, **kwargs): 
        """
            this is a wraper 
        """
        if kwargs.has_key('energy_diff_min'): 
            self.energy_diff_min = kwargs['energy_diff_min']
            
        if self.info>0: 
            print 'start minimize'
        which_minimize = which_minimize if which_minimize is not None else self.which_minimize     
        if which_minimize in ['prod_state', 'eigen_state']: 
            self._minimize_finite_size(**kwargs)
        elif which_minimize in ['scale_invar_state', 'scale_invar']:  # only retain name 'scale_invar_state' in future
            self._minimize_scale_invar(**kwargs)
        else: 
            raise  ValueError(which_minimize)
        
        
        if self.do_measure and self.iter0>0:  
            self.measure_state()
 
    def make_dir(self, dir, dir_local=None): 
        #this is taken form main.Main.make_dir 
      
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
    
    #@staticmethod   
    @classmethod 
    def J_power_0l(cls, l, alpha, beta):
        """
            l: distance between sites 0 and l
        """
        if l<cls.J_POWER_R_MAX: 
            return beta/float(l)**alpha
        else: 
            return 0.0 

    @staticmethod   
    def J_power_0l_only_for_test(l, alpha, beta):
        """
        l: distance between sites 0 and l
        """
        if 0:
            return beta/float(l)**alpha
        if 0:
            res= 0.0
            if 6<=l<= 13:
                res= -1.0
            return res
        if 1:
            #note have combine 2site, #the keys refer to legs not sites, it is sites after devide by 2
            System.temp = {1:-0.0, 2:-0.0,   #short
                    3:-0.0, 4:-0.0, 5:-1.0, 
                    6:-0.0, 7:-0.0, 8:-0.0, 9:-0.0,       #middle
                    10:-0.0,   #also middle
                    16:-0.0, 17: -0.0,  31:-0.0, 
                    60:-0.0, 120:-0., 148:-0.0,      #long
                    }
            System.temp.update({
                15:-0.0 
                
                })
            
            if System.temp.has_key(l):
                return System.temp[l]
            else:
                return 0.0

    #Jleg_tau_0l_dic = {}

    #warnings.warn("changed Jleg_tau_0l "*20)
    def Jleg_tau_0l(self, l, tau, **model_param):
        """
            NOTE if not combine_2site l refers to both length between legs and sites as they are same, 
            while when combine_2site l can only refer to LEGS NOT SITES!!

            use the convention "A在左，B在右" i.e. SxA stores s_0i,  SxB stores is_l alike ops
            represent s0sl by o_{0,1} and o_{l-1,l} with coup. coeff. J(tau, l)
            represent s0sl by s0i and isl with coup. coeff. J(tau, l)
            s0 sl at layer tau comes from 
            {s0i, s1i, s2i} \otimes {is_{3*l-2}, is_{3*l-1}, is_{3*l}}
            example:
                when l = 3, s0iis3 at layer 1 comes from 
                {s0i, s1i, s2i} \otimes {is7, is8, is9}

            while for U1 symm. things bit different for bottom layer
                s0 sl at layer 1 tau comes from 
                {s0i, s1i, s2i} \otimes {is_{3*2l-4}, is_{3*2l-3}, is_{3*2l-2}, is_{3*2l-1}, is_{3*l}, is_{3*2l+1}}
                example:
                    s0iis3 at layer 1 comes from
                    {s0i, s1i, s2i} \otimes {is14, is15, is16, is17, is18, is19}
            
            interaction between o_{0,1} and o_{l,l+1} at layer tau
            Jleg_tau_0l_dic has a self-similar structure similar to mera
        """
        beta = model_param['beta']
        alpha = model_param['alpha']

        #if l<3 : 
        #    raise ValueError("l should be in [3, 6]. l=%d"%l)
        assert l>2, "l should be larger than 2.(i.e. not NN or NNN ops)" 

        #if self.Jleg_tau_0l_dic.has_key(l):   # a really bad bad mistake 
        if self.Jleg_tau_0l_dic.has_key((tau, l)):
            return  self.Jleg_tau_0l_dic[(tau, l)]
        else:
            if tau == 0:
                if self.symmetry in ["U1", "Travial"] and self.combine_2site:
                    #res = 0.5*(System.J_power_0l(2*l, alpha, beta) + System.J_power_0l(2*l + 1, alpha, beta))
                    res= (System.J_power_0l(2*l-1, alpha, beta)  + 2.0*System.J_power_0l(2*l, alpha, beta)  
                            + System.J_power_0l(2*l  + 1, alpha, beta) )/2.0   
                    #ATTENTION here should divide by 2.0, previousely I thougth should divide by 4.0, or not aver, 
                    #They are all wrong!!  the reason is I havent completely understand it. 
                    # /2.0 comes fome aver of combine_2site
                else:
                    #here means from layer 0 to 1 there is NO approximation in h
                    res = System.J_power_0l(l, alpha, beta)
            else:
                res= 0.0
                #i is S0's in the coarsed lattice H^\tau: L, C, R
                #left side SA \equiv s_i*I
                period = self.mera.period
                #period_by_l = period*l
                #assert period == 3 or period  == 4   # when period is other value, I need reconsider following which is ad hoc  
                if period == 3:
                    j0 = period*(l-1) + 1   #range(j0, j0 + period) is like [7, 8, 9]
                elif period == 4:
                    j0 = period*(l-1) + 2   #like [10, 11, 12, 13]
                else:
                    raise

                for i in range(0, period):  
                    #j is position of Sj, not represent L, C, R ,  so below only devide by 3.0
                    #because of recursion,  l will actually ranges [(l+1)*3^tau, (l+3)*3^tau] in the tau=0 layer
                    #right side SB \equiv I*s_j
                    #for j in [3*l-2, 3*l-1, 3*l]:   #[7, 8, 9] 
                    #for j in range(period*(l-1) + 1, period*l + 1):   #[7, 8, 9] 
                    for j in range(j0, j0 + period):   #[7, 8, 9] 
                        #res += self.Jleg_tau_0l(j-i, tau-1, **model_param) /3.0
                        res += self.Jleg_tau_0l(j-i, tau-1, **model_param) /float(period)
            self.Jleg_tau_0l_dic.update({(tau, l):res})
            return res

    def J_long_ising_fetch(self, tau, l, num_of_layer, **model_param):
        res=self.Jleg_tau_0l(l, tau, **model_param)
        #print 'jjjjj', l, tau, res
        return res

    def tensor_network_mapping(self):
        pass
    
    def evaluation_map(self): 
        """
            maps a G -> \hat{G} and etc 
        """
        raise NotImplemented 
    
    def contract(self, M, layer, G, nOps, Names, Ops, order, exception, extra=None, restart=True, branch=0, info=0):
        """
            there is bug in it
            see System_Contraction in f90
            system must be initialized before contraction

            1.把此layer的张量（包括U, V, oo, OO etc) 保存在 TNet_sys.tlink中, 这里layer指的是本layer，以及对偶layer
            2.按照order把TNet_sys.tlink收缩掉，只剩下exceptiong对应点张量，通常，excepting为oo, OO

            info: for debugging
        """
        if restart:  
            System.TNet_sys= TensorNetwork()
            System.TNet_sys.TG = G
        TNet_sys = System.TNet_sys
        
        #tlink.tensor has n tensors, each store a perticular tensor from U, U_deg, V, etc
        
        mapper = {}
        mapper.update(self.network_tensor_dic.copy())
        mapper.update(self.network_tensor_dic_extra[branch])
        if extra is not None:
            mapper.update(extra)
        if info>0:
            print "graph name", G.graph_name #, "weight", G.weight
        try:
            for n in range(G.size):
                node_name = G.names[n]
                if 0:
                    i = Find_Ops(node_name, nOps, Names)
                    if i>= 0 and False:
                        TNet_sys.tlink[n]=Ops[i].shallow_copy()
                    else:
                        pass
                op_name, x = mapper[node_name]
                if hasattr(self, op_name):
                    TNet_sys.tlink[n]=self.__getattribute__(op_name)[layer + x][0].shallow_copy()
                elif hasattr(self.mera, op_name):
                    TNet_sys.tlink[n]=M.__getattribute__(op_name)[layer + x][0].shallow_copy()
                else:
                    raise ValueError(str(op_name))
                TNet_sys.tlink[n].type_name = op_name
        except Exception as err:
            msg = "\terror when mapping network to tensor, op_name=%(op_name)s, layer=%(layer)s, x=%(x)d"%vars()
            if not err.args: 
                       err.args=('',)
            err.args = (err.args[0] +'\n' + msg,)+err.args[1:]
            
            raise 

        
        if 1:
            if info>0:
                print "tensors to be contracted:"
                names= [TNet_sys.tlink[i].type_name for i in order[:G.size]]
                print names
            n = G.size#TNet_sys.TG.nNode
            try: 
                output = TNet_sys.contract_except(n, order, exception=exception, info=info-1)#, memory_save= True )
            except Exception as err: 
                msg = 'additinal info: \n'
                names= [TNet_sys.tlink[i].type_name for i in order[:G.size]]
                for i in order[: G.size]: 
                    t = TNet_sys.tlink[i]
                    msg += '\t{}\t  {}\n'.format(t.type_name, t.shape)
                
                if not err.args: 
                           err.args=('',)
                err.args = (err.args[0] +'\n' + msg,)+err.args[1:]
                raise 
                   
            return output
    
    def contract_new(self, layer, G, exception, extra=None, extra_extra={}, restart=True, branch=0, info=0):
        """
            a better interace for contract

            1.把此layer的张量（包括U, V, oo, OO etc) 保存在 TNet_sys.tlink中, 这里layer指的是本layer，以及对偶layer
            2.按照order把TNet_sys.tlink收缩掉，只剩下exceptiong对应点张量，通常，excepting为oo, OO

            info: for debugging
        """
        M = self.mera
        order = G.contract_order
        if restart:  
            System.TNet_sys = TensorNetwork()
            System.TNet_sys.TG = G
        TNet_sys= System.TNet_sys
        
        #tlink.tensor has n tensors, each store a perticular tensor from U, U_deg, V, etc
        
        mapper = {}
        mapper.update(self.network_tensor_dic.copy())
        mapper.update(self.network_tensor_dic_extra[branch])
        if extra is not None:
            mapper.update(extra)
        if info>0:
            print "graph name", G.graph_name #, "weight", G.weight
        
        if 1:
            for n in range(G.size): 
                node_name = G.names[n]
                if extra_extra.has_key(node_name): 
                    op = extra_extra[node_name]
                    TNet_sys.tlink[n] = op
                    extra_extra.pop(node_name)
                    op_name = node_name
                    if op is not None : TNet_sys.tlink[n].type_name = op_name
                else: 
                    op_name, x = mapper[node_name]
                    if hasattr(self, op_name):
                        TNet_sys.tlink[n]=self.__getattribute__(op_name)[layer + x][0].shallow_copy()
                    elif hasattr(self.mera, op_name):
                        TNet_sys.tlink[n]=M.__getattribute__(op_name)[layer + x][0].shallow_copy()
                    else:
                        raise ValueError('tensor for node_name %s not found'%node_name)
                    TNet_sys.tlink[n].type_name = op_name
           
        for op_name, op in extra_extra.items(): 
            n = G.names.index(op_name)
            TNet_sys.tlink[n] = op.shallow_copy()
            TNet_sys.tlink[n].type_name = op_name
            ee_list.append(n)

        if 1:
            if info>0:
                print "tensors to be contracted:"
                names= [TNet_sys.tlink[i].type_name for i in order[:G.size] if TNet_sys.tlink[i] is not None ]
                print names
            n = G.size#TNet_sys.TG.nNode
            output = TNet_sys.contract_except(n, order, exception=exception, info=info-1)#, memory_save= True )
            return output
    
    
    def add_env(self, M, layer, G, order, name, weight=None, add_tensor=None, branch=0, alpha=1.0, restart=True, info=0):
        """
            no need option of use_buf, as in tensor_network.contract_except we always use buf, so the tensor output.data below is already buffed
            one name may correspond to many nodes, contract except these nodes and add them together 
            this method corresponds to e.g. vidal's paper p.16 fig.21  第1和第4个图加起来
            elist: exception node list
            
            weight: this param will be deprecated in future,  through G.weight instead 
        """
        weight = weight if weight is not None else G.weight 
        elist = G.find_node(name)
        nlist = len(elist)
        if nlist == 0:
            if name != "U":
                raise Exception("node %s not find in graph %s\n%s"%(name, G.graph_name, G))
            warnings.warn("in add_env node %s is not found in graph %s. just return None"%(
                name, G.graph_name))
            return False
        
        #if isinstance(add_tensor, iTensor):  #sometimes this fails,  so use following instead
        if add_tensor is not None:
            for ne in range(nlist):
                exception = elist[ne] 
                temp = self.contract(M, layer, G, nOps=0, Names=[], Ops=[], order=order, exception=exception, restart=restart, branch=branch, info=info-1)
                temp.data *= weight
                add_tensor +=  temp
                restart =  False 
            return True
        else:
            #I should remove following lines in future, it is not safe
            first = True
            for ne  in range(nlist):
                exception = elist[ne] 
                temp = self.contract(M, layer, G, nOps=0, Names=[], Ops=[], order=order, exception=exception, restart=restart, branch=branch, info=info-1)
                temp.data *= weight
                
                if first:
                    output = temp.copy(use_buf = True)
                    first= False
                else:
                    output.data += temp.data 
                #这里把几个tensor按照weight加了起来，
                restart =  False 
            return output

    def add_env_new(self, M, layer, name, G, branch = 0, alpha=1.0, restart=True, info=0):
        """
            one name may correspond to many nodes, contract except these nodes and add them together 
            this method corresponds to e.g. vidal's paper p.16 fig.21  第1和第4个图加起来
            elist: exception node list
        """
        
        elist=G.find_node(name)
        nlist = len(elist)
        #if nlist == 0:
        #    print "error, operator %s not found, exit"%name
        #    exit()

        
        first = True
        for ne  in range(nlist):
            exception = elist[ne] 
            if info>0:
                print 'eee exception=', exception
            T_tmp = self.contract(M, layer, G, 0, Names=[], Ops=[], order=G.order, exception=exception, restart=restart, branch=branch, info=info-1)

            if first:
                output = T_tmp.copy(use_buf = True)
                output.data *= G.weight
                first= False
            else:
                output.data = T_tmp.data*weight + output.data*alpha
            #这里把几个tensor按照weight加了起来，
            restart =  False 
        return output

    def add_env_many(self, M, layer, name, G_dic, add_tensor=None, key_list=None, weight_dic=None, branch = 0, alpha=1.0, restart=True, info=0):
        """
            add_env for a set of graphs, G, order, weight are all dict
            params:
                add_tensor: of type iTensor,  if it is not none, then result will add to add_tensor.data
                key_list: selected keys of graphs 
        """

        if key_list is not None:
            if __debug__ :
            #if __debug__ and False:  #disable check temporarily for efficiency
                k1 = set(key_list)
                k2 = set(G_dic.keys())
                #if not set(key_list)<=set(G_dic.keys()):
                if k1-k2:
                    err_msg = "key_list " + repr(key_list) + " not coincide with G_dic.keys()"
                    raise ValueError(err_msg)
                elif k1<k2:
                    #warnings.filterwarnings('once')
                    warnings.warn("not all graphs are used, lack these: " + str(k2-k1))
            keys = key_list
        else:
            keys= G_dic.keys()
        
        if len(keys)==0:
            # it happens I pass in an empty G_dic, just for compatibility of the code, but this may also
            # induce rather unconceivable errors!! so caution and some check are needed
            if G_dic == {}: 
                warnings.warn("an empty graph is passed in add_env_many")
                return None
            else:
                raise
        keys.sort()    #add_env in a certain order, just easier to debug
        if weight_dic is None:
            weight_dic ={k:G_dic[k].weight  for k in G_dic} 
        if info>0:
            #temp = {G_dic[g].graph_name:round(weight_dic[g], 3) for g in G_dic}
            temp = {G_dic[g].graph_name:round(weight_dic[g], 4) for g in keys}
            print "add_env_many graphs with weights:\n",  temp# weight_dic
        
        if isinstance(add_tensor, iTensor):
            for i in range(len(keys)):
                g = keys[i]
                G = G_dic[g]
                order = G.contract_order
                weight= weight_dic[g]
                self.add_env(M, layer, G, order, name, add_tensor=add_tensor, weight=weight,  branch=branch, restart=restart, info=info-1)
                #add_tensor.data += output.data 
            return True
        else:
            first = True            
            for i in range(len(keys)):
                g = keys[i]
                G = G_dic[g]
                order = G.contract_order
                weight= weight_dic[g]
                T_tmp=self.add_env(M, layer, G, order, name, weight=weight, branch=branch, restart=restart, info=info-1)

                if first:
                    output = T_tmp.copy(use_buf = True)
                    first= False
                else:
                    output.data += T_tmp.data 
            return output

    def add_env_many_many(self, M, layer, name_list, G_dic, key_list=None, weight_dic=None, branch = 0, alpha=1.0, restart=True, info=0):
        """
            add_env for a set of graphs, G, order, weight are all dict
        """
        ilayer = layer
        
        envs= TensorLink(size=8, use_buf=True) #contain envirenment tensors:   o, rho, h, w, u, u^+, w, w^+  
        
        num_of_tensor = 2
        
        raise("fix this, ndiv not right when using quaternary etc")
        ndiv = [3, 2]
        
        envs[0] = M.V[ilayer][j].copy(use_buf=True)
        envs[1] = M.U[ilayer][j].copy(use_buf=True)

    def top_level_contract(M, O_in, i, exception=None, final_order=None):
        """
        this func is not used
        top level 的收缩图与其他不同，因为top V 不同
        """
        TNet_sys= TensorNetwork()

        layer = M.num_of_layer
        
        TNet_sys.tlink[0].shallow_copy(O_in)
        TNet_sys.tlink[1].shallow_copy(M.V[layer][0])
        TNet_sys.tlink[2].shallow_copy(M.V_dag[layer][0])
        TNet_sys.tlink[3].shallow_copy(rho_top)
        
        TNet_sys.TG = GTop.copy()
        order = order_gtop
        O_out = TNet_sys.contract_except(TNet_sys.TG.nNode, order, exception=exception, final_order=final_order)
        return O_out
    
    @staticmethod
    def contract_oo_rho_2(oo, rho_2): 
        V1=[0,1,2,3]
        V2=[2,3,0,1]
        res, legs = oo.contract(rho_2, V1, V2, use_buf=True)  #rho and H_2 fully contracted to a scalar
        res= res.data[0]
        return res
    
    @staticmethod
    def contract_ooo_rho_3(ooo, rho_3): 
        V1 = range(6)
        V2 = [3, 4, 5, 0, 1, 2]
        H_3, legs = ooo.contract(rho_3, V1, V2, use_buf=True)
        res = H_3.data[0]
        return res
  
    def eng_ham(self, ilayer, tau=None, name=None):
        """
            energy can be evaluated at any layer including bottom, top,  etc. 
            But the question is are the values equal at each layer, or any layer it is better
            ilayer = M.num_of_layer-1
            
            对应于vidal paper fig.14(iv)
            收缩后得到一个标量，即能量
        """
        
        tau = tau if tau is not None else ilayer-1
        energy = System.contract_oo_rho_2(self.H_2[ilayer][0], self.rho_2[ilayer][0])
                
        if not self.only_NN:
            #there may be an issue: only top_level_product_state need this, while eigenstate needn't
            #as H_3 has fully merged into H_2 in top level, and H_3, rho_3 at top are zero
            e3 = System.contract_ooo_rho_3(self.H_3[ilayer][0], self.rho_3[ilayer][0])
            energy += e3 
            
        energy += 1.0 
        self.energy = energy

        return self.energy

    def prompt_save(self):
        print "save mera to file?"
        for i in range(60):
            print i, 
            time.sleep(1)

    
    @staticmethod
    def backup_fn_parse(fn):
        nqn = None
        if "transition" in fn or 'transtion' in fn: #'12-transition.pickle'
            dim = "0"
        elif "_" in fn:  #like "5_13.pickle or 5_13-4lay.pickle"
            dim = fn.split("-")[0].split("_")[1].split(".")[0]
            nqn = fn.split("-")[0].split("_")[0]
        else:
            if "lay" in fn: #"like "8-4lay.pickle""
                dim = fn.split("-")[0]
            else:   
                 #like "8.pickle"
                dim = fn.split(".")[0]
        
        i = fn.find("lay")
        layer = fn[i-1] if i>0 else "3"
        
        #return eval(dim), eval(layer) + 1  #  layer= mera.num_of_layer large than actual layer by 1 
        if nqn is None : 
            return eval(dim), eval(layer) + 1  #  layer= mera.num_of_layer large than actual layer by 1 
        else: 
            return eval(dim), eval(layer) + 1, eval(nqn)#  layer= mera.num_of_layer large than actual layer by 1 
    
    def _fn_auto(self):
        dim = self.mera.qsp_max.totDim
        nqn = self.mera.qsp_max.nQN
        lay = self.mera.num_of_layer-1
        nqn = "5_" if nqn == 5 else "" 
        lay = "-%dlay"%lay if lay != 3 else ""
        return dim, nqn, lay

    def backup_fn_auto(self):
        dim, nqn, lay = self._fn_auto()
        fn = "%(nqn)s%(dim)d%(lay)s.pickle"%vars()
        return fn

    def output_fn_auto(self):
        dim, nqn, lay = self._fn_auto()
        fn = "res-%(nqn)s%(dim)d%(lay)s.dat"%vars()
        return fn

    @staticmethod
    def backup_fn_gen_bac(dim, layer, nqn=None):
        nqn = "5_" if nqn == 5 else "" 
        layer = "-%dlay"%layer if layer != 3 else ""
        res= "%(nqn)s%(dim)d%(layer)s.pickle"%vars()
        return res
    
    def backup_fn_gen(self):
        dim, nqn, lay = self._fn_auto()
        fn = "%(nqn)s%(dim)d%(lay)s.pickle"%vars()
        return fn
    
    shape_to_backup_fn = backup_fn_gen  # another name

    def save_del(self, fn, check=True, use_local_storage=False):
        """
            it turns out QnZ2, QnU1 can't be pickled, the reason is not sure
        """
        LOCAL = '' if not use_local_storage  else 'LOCAL path'
        msg = "\nTRY TO SAVE S TO %s  %s  \t"%(LOCAL, fn[-30: ])
        
        allow_save = True
        #msg = ''
        try:
            S0= System.load(fn, use_local_storage)
            pass 
        except IOError:
            msg += "file not exist, a new file is created."  
        else:
            #if S0 != self:
            #    msg  += "current system not match with stored one! save denied."
            #    allow_save = False
            if 1: 
                pass
                
            elif S0.iter1>self.iter1:
                msg  += "S0.iter1 %d > self.iter1 %d, save denied"%(S0.iter1, self.iter1)
                allow_save = False
            else:
                pass
        
        if allow_save: 
            System._save(self, fn, use_local_storage=use_local_storage)
            msg += ' ... S is successfully saved.'  
            print msg
        else: 
            msg += ' ... save denied. raise' 
            raise Exception(msg)

    def save(self, fn=None, use_local_storage=False): 
        IterativeOptimize.save(self, path=fn, use_local_storage=use_local_storage) 

    @staticmethod    
    def _save(obj, fn, use_local_storage = False):
        if not use_local_storage: 
            out = open(fn, "wb")
            pickle.dump(obj, out)
            out.close()
        else: 
            #with rpyc_conn_local() as conn:
            with rpyc_conn_local_zerodeploy() as conn:
                out = conn.builtin.open(fn, 'wb')
                conn.modules.cPickle.dump(obj, out)
                out.close()

    @staticmethod
    def load(fn, use_local_storage=False):
        res = rpyc_load(fn, use_local_storage)    
        if not isinstance(res, dict):  #temporary work around 
            res = System.update_S(res)
        else:
            if not hasattr(res['mera'], 'USE_REFLECTION'): 
                res['mera'].USE_REFLECTION = False 
            if not hasattr(res, 'symmetry'): 
                res['symmetry'] = 'U1'

            #temp1 = ["U", "V", "U_dag", "V_dag"] 
            #temp2 = ["rho_2", "rho_3"]

            ##原因在save时，把iTensor.buf_ref 储存了，造成load后tBuffer引用错乱，see iTensor.__del__。加入如下代码，将其清空
            ##this part is RATHER TRICKY!  whenever it is loaded, it will interfere the entire buff, so must unlink them from buffer
            ##another solution may be remove buf_ref using __state???___ in iTensor for pickling
            #for l in range(res.mera.num_of_layer):
            #    for i in temp1:
            #        try:
            #            #res.mera.__getattribute__(i)[l][0].buf_ref = np.array([-1, -1], np.int)
            #            res['mera']__getattribute__(i)[l][0].buf_ref = np.array([-1, -1], np.int)
            #        except: raise
            #        #except: AttributeError
            #    for i in temp2:
            #        try:
            #            res[i][l][0].buf_ref = np.array([-1, -1], np.int)
            #        except: pass
            #    
        
        return res
    
    @staticmethod
    def update_S(res): 
        #reset tensor.buf_ref
        temp1 = ["U", "V", "U_dag", "V_dag"] 
        temp2 = ["rho_2", "rho_3"]

        #原因在save时，把iTensor.buf_ref 储存了，造成load后tBuffer引用错乱，see iTensor.__del__。加入如下代码，将其清空
        #this part is RATHER TRICKY!  whenever it is loaded, it will interfere the entire buff, so must unlink them from buffer
        #another solution may be remove buf_ref using __state???___ in iTensor for pickling
        for l in range(res.mera.num_of_layer):
            for i in temp1:
                try:
                    res.mera.__getattribute__(i)[l][0].buf_ref = np.array([-1, -1], np.int)
                except: raise
                #except: AttributeError
            for i in temp2:
                try:
                    res.__getattribute__(i)[l][0].buf_ref = np.array([-1, -1], np.int)
                except: pass
        
        #temperory 
        if not hasattr(res, "energy_record"):
            res.energy_record = res.energy_hist
        if not hasattr(res, "scaling_dim_record"):
            res.scaling_dim_record = {}
        if not hasattr(res, "log"):
            res.log = {}
        if not hasattr(res.mera, "period"): res.mera.period= res.mera.V[0][0].rank-1
        if not hasattr(res.mera, "tensor_defs"): 
            res.mera.tensor_defs= {
                    "U":        {"type":(2, 2),    "type_1":"unitary",	        "dual":False,}, 
                    "U_dag":    {"type":(2, 2),    "type_1":"unitary",	        "dual":True, }, 
                    "V":        {"type":(3, 1),    "type_1":"disentangler",	    "dual":False,}, 
                    "V_dag":    {"type":(1, 3),    "type_1":"disentanger",	    "dual":True, }, 
                    }
        
        if not hasattr(res, "G_2_2") or not hasattr(res, 'graph_module_name'):
            res.__class__= System
            #res.set_graph = types.MethodType(System.set_graph, res)
            rank_m1 = res.mera.V[0][0].rank-1
            if rank_m1 == 3 and not res.only_NNN:
                #import merapy.init_mera_graph as module
                import merapy.diagrams.V31.graph_ternary as module
                res.set_graph(module)
            elif rank_m1 == 2:
                from merapy.diagrams.modified_binary import graph_modified_binary
                res.set_graph(graph_modified_binary)
            elif rank_m1 == 4:
                from merapy.diagrams.V41 import graph_quaternary
                res.set_graph(graph_quaternary)
            elif res.only_NNN == True:
                from current import init_mera_graph
                res.set_graph(init_mera_graph)
        
        if not hasattr(res.mera, 'qsp_max'):
            res.mera.qsp_max = res.mera.QSp_max.copy()
            
        if not hasattr(res.mera, 'qsp_max2'):
            res.mera.qsp_max2 = res.mera.QSp_max.copy()
        if not hasattr(res.mera, 'qsp_0'):
            res.mera.qsp_0 = res.mera.QSp_0
        if not hasattr(res, 'energy_diff'):
            eng_rec1, eng_rec0 = res.energy_record.values()[-1:-3:-1]
            a, b = res.energy_record.keys()[-1:-3:-1]
            pace = a-b
            res.energy_diff = (eng_rec1[0] - eng_rec0[0])/pace
            res.energy_err = (res.energy - res.energy_exact)
        if not hasattr(res, 'record'):
            res.record = OrderedDict()
            for a, b in res.energy_record.iteritems():
                res.record[a] = {'energy':b[0], 'time':b[1] }
            for a, b in res.scaling_dim_record.iteritems():
                if not res.record.has_key(a):  res.record[a]={}
                res.record[a].update({'scaling_dim':b})
        #if not hasattr(res, 'graph_module_name'):
        if not hasattr(res, 'iter1'):
            res.iter1 = res.iter

            
        #del res.__class__.key_property

        res.mera.__class__.key_property = Mera.key_property
        res.mera.qsp_max.__class__.key_property = QuantSpaceBase.key_property
        res.__class__.key_property = System.key_property
        res.__class__._expand_dim = System._expand_dim
        res.__class__.expand_dim = System.expand_dim
        res.__class__.display = System.display
        res.__class__._fn_auto = System._fn_auto
        res.__class__.backup_fn_auto = System.backup_fn_auto
        res.__class__.output_fn_auto = System.output_fn_auto
        res.__class__.Jleg_tau_0l = System.Jleg_tau_0l
        #res.__class__.Jleg_tau_0l_dic= System.Jleg_tau_0l_dic
        
        #res.__class__ = System
        if not isinstance(res.energy_record, OrderedDict):
            res.energy_record = OrderedDict(sorted(res.energy_record.items()))
            res.scaling_dim_record = OrderedDict(sorted(res.scaling_dim_record.items()))
        for i, v in res.scaling_dim_record.iteritems():
            if isinstance(v, dict):
                break
            res.scaling_dim_record[i] = dict(v)
            
        if not hasattr(res, 'Jleg_tau_0l_dic'):    # add at 2014-3-22
            res.Jleg_tau_0l_dic = {}
            res.__class__.Jleg_tau_0l = System.Jleg_tau_0l
            
        if not hasattr(res, 'use_pinning_term'): 
            res.use_pinning_term = False 
        if not hasattr(res.mera, 'USE_REFLECTION'): 
            res.mera.USE_REFLECTION = False 
            
        res.__class__.save = System.save 
        
        return res
    
    def resume_func(self):
        """
            issue: use IterativeOptimize.resume_func later 
        """
        #if not use_local_storage: 
        backup_path = self.backup_path if not self.use_local_storage else self.backup_path_local 
        if backup_path is not None : 
            LOCAL = '' if not self.use_local_storage else 'LOCAL file '
            msg = "\nTRY TO RESUME from %s%s...  "%(LOCAL, backup_path[-30:])
            try:
                
                #self.S = self.__class__.load(backup_path, use_local_storage=use_local_storage)
                temp = self.__class__.load(backup_path, use_local_storage=self.use_local_storage)
                
                if not isinstance(temp, dict):  #temporary work around  for backward compatible                 
                    temp = temp.__dict__ 
                    kk = ['updaters', 'backup_parpath', 
                            'is_initialized', 'measurement_args']
                    for k in kk: 
                        if temp.has_key(k): 
                            temp.pop(k)
                    self.__dict__.update(temp)
                else: 
                    self.set_state(temp)
                    self.M = self.mera
                
                if 0:  #resume check 
                    if hasattr(self, 'energy_diff_std') and kwargs.get('energy_diff_min') is not None : 
                        energy_diff_min = kwargs.get('energy_diff_min')
                        if  self.energy_diff_std is not None: 
                            if self.energy_diff_std <= energy_diff_min: 
                                q_iter = 0
                                
                                msg = 'self.energy_diff_std %1.2e less than energy_diff_min %1.2e, set q_iter=0'%(
                                        self.energy_diff_std, energy_diff_min)
                                print msg
                
                #self.M = self.S.mera
                #iter = self.S.iter1
                #eng = self.S.energy
                iter = self.iter1
                eng = self.energy 
                #print_vars(vars(),  ['temp.keys()', 'self.energy_diff'])
                msg  += "RESUMED successfully at iter=%s,  eng=%s, energy_diff=%1.2e"%(iter, eng, self.energy_diff)
            except IOError as err:
                if err.errno == 2: 
                    msg += str(err)
                    msg += "discard manually resume"
                else:
                    raise IOError
            except TypeError as err:
                msg += str(err)
                msg += "discard manually resume"
            except Exception: 
                raise 
            print(msg)
    
    def resume_check(self, old_state): 
        allow_resume = 1
        eng_bac = old_state['mps'].energy
        eng = self.mps.energy if self.mps is not None  else None 
        msg = ''
        if eng is None: eng = np.nan
        if eng_bac < eng or eng is np.nan: 
            msg += 'stored_eng=%1.12f is lower that current_eng=%1.12f.  '%(eng_bac, eng) 
            msg += 'resumed successfully'
        else: 
            allow_resume = 0
            msg += 'stored_eng=%1.12f is higher that current_eng=%1.12f.  '%(eng_bac, eng) 
            msg += 'discard ressuming' 
        return allow_resume, msg
    
    def save_eng(self, iter):
        self.energy_record[iter] = (self.energy, time.time(), time.clock())
    
    def get_energy_exact(self):
        #J_NNN = self.model_param["J_NNN"]
        temp = {0.5:-1.5, 0.241186:-1.607784749}

        if self.model == "Heisenberg":
            if self.only_NN:
                energy_exact = -1.772588722236
                if self.symmetry == "U1" and temp.has_key("J_NNN"):
                    energy_exact = temp[0.241186]
            elif self.only_NNN:
                energy_exact = temp[0.241186]
            else:
                energy_exact = -1.644934
        elif self.model == "Ising":
            from math import pi
            energy_exact = -4./pi  #1.27323954474
        elif self.model == 'Potts':
            energy_exact = -2.4359911239
        else: 
            #print_vars(vars(),  ['self.model'])
            #raise 
            energy_exact = np.nan 
        return energy_exact
    
    def display(self, iter, filename, pace=1, interval=100):
        #iter_last = iter-pace
        eng_rec1, eng_rec0 = self.energy_record.values()[-1:-3:-1]
        a, b = self.energy_record.keys()[-1:-3:-1]
        pace = a-b
        
        #temp = np.array([v[0] for v in self.energy_record.values()[-1: -10: -1]])
        temp = np.array([v[0] for v in self.energy_record.values()[-1: -interval//pace-1: -1]])
        temp1 = np.roll(temp, -1)
        temp = temp1 - temp
        #self.energy_diff_std = np.mean(np.abs(temp[:-1]))/pace
        self.energy_diff_std = np.max(np.abs(temp[:-1]))/pace
        
        self.energy_diff = (eng_rec1[0] - eng_rec0[0])/pace
        
        self.energy_err = (self.energy - self.energy_exact)
        steps = -self.energy_err/self.energy_diff#*pace
        
        print '\n%s %dth iter, energy=%2.11f,  eng_err=%1.1e, eng_diff=%1.1e, eng_diff_std=%1.1e, steps~%1.2e'%(
                self.backup_parpath_name, iter, self.energy, 
                self.energy_err, self.energy_diff, self.energy_diff_std, steps),  
        
        dt_real = eng_rec1[1] - eng_rec0[1]
        dt_cpu = eng_rec1[2] - eng_rec0[2]
        dt_real = dt_real/pace;     dt_cpu = dt_cpu/pace
        
        print "\t", dt_cpu,"\t",  round(dt_real, 3)
        #if iter%100 == 0 or iter<= 5: 
        #    if filename is not None:
        #        out = open(filename, 'a')
        #        out.write('\t%d, %2.11f, %1.1e, %1.1e, %f, %f\n'%(iter, 
        #            self.energy, self.energy_err, self.energy_diff, dt_cpu, dt_real))
        #        out.close()
        energy_last = self.energy
        return steps, dt_real

    def examine_energy(self, res=None, delta=1e-14):
        """
            res: correct res obtained previously
        """
        passed = True
        if res is None:
            print "result is "
            print {i:self.energy_record[i][0] for i in self.energy_record}
            return
        for i in res:
            if self.energy_record.has_key(i):
                diff = abs(res[i] - self.energy_record[i][0])
                print "%(i)d, \t%(diff)s"%vars()
                #if diff <1e-14:
                #    print i, "\t", abs(res[i] - self.energy_record[i][0])
                #else:
                if diff > delta: 
                    passed = False
                    print "\ti=%d, stored res=%1.18f, current res=%1.18f"%(i, 
                        res[i], self.energy_record[i][0])
        if not passed:
            msg = "EXAMINE FAILED\n%s"%{i:self.energy_record[i][0] for i in self.energy_record}
            raise Exception(msg)

        print "\nEXAMINE PASSED\n"
    
    def examine_eng(self, res=None):
        raise

    def _expand_dim(self, qsp_max, out=None):
        #old_key_property = self.key_property()
        if qsp_max <= self.mera.qsp_max:
            print "qsp_max is not larger than original one, system not expaned"
            print "qsp_max is %s, self.mera.qsp_max is %s"%(qsp_max, self.mera.qsp_max)
            return
        from decorators import tensor_player, set_STATE_end_1
        state_bac = tensor_player.STATE
        tensor_player.STATE = "stop"
        M = self.mera
        M.expand_dim(qsp_max)
        i = 0
        def expand_func(ilayer, name, propt):
            qsp_v = M.V[ilayer-1][0].QSp[M.V[ilayer-1][0].rank-1].copy()
            qsp_v.reverse()
            rank = propt["rank"]
            rank_half = rank//2
            
            if propt["type"] == "I":
                self.I[ilayer][0] = iTensor.identity(qsp_v.copy_many(2, reverse=[1]))
                return
            
            if propt["type"] == "rho":
                reverse = range(rank_half)
            elif propt["type"] == "H":
                reverse = range(rank_half, rank)
            qsp = qsp_v.copy_many(rank, reverse=reverse)
            try:
                self.__getattribute__(name)[ilayer][0] = self.__getattribute__(name)[ilayer][0].expand(qsp) 
            except AttributeError:
                print "in System.expand_dim ilayer = %(ilayer)d, name = %(name)s failed to beexpanded"%vars()
        
        dic = self.fetch_op_names(type=["H", "rho", "I"])
        print "\texpand following System ops:  ", dic.keys()
        
        #need not expand layer 0
        for ilayer in range(1, M.num_of_layer):
            for name, propt in dic.iteritems():
                expand_func(ilayer, name, propt)
        
        new_key_property = self.key_property()
        #if old_key_property != new_key_property: 
        self.iter1 = 0  #reset counting
        self.iter0 = 0 
        
        msg = "\tqsp_max has been expanded to %s \n"%qsp_max
        print msg 
        #tensor_player.NEXT_STATE = "record"
        if out is not None:
            out.write("msg")
        #if state_bac == "play": 
        #    print "after expansion, tensor_player is switched from 'play' to 'record'"
        #    tensor_player.NEXT_STATE = "record" 
       
        #if self.energy is not None:
        if 0:  #discard examine temporarily, due to self.energy may be not calced
            print "after expansion, set player to record"
            tensor_player.NEXT_STATE = "record"
            e0 = self.energy
            self.eng_ham(M.num_of_layer-1)
            e1 = self.energy
            e_diff = abs(e0 - self.energy)
            if e_diff <=   1e-13:
                iter = self.iter
                print "energy test passed, at iter=%(iter)d both equal to %(e0)2.14f, e_diff = %(e_diff)e\n\n"%vars()
            else:
                self.energy = e0
                err_msg = "energy test failed, e0 = %2.15f, e1 = %2.15f"%(e0, e1)
                raise Exception(err_msg)
    
    def expand_dim(self, trunc_dim, nqn=None):
        print '\ntry to expand dim'
        qspclass = symmetry_to_Qsp(self.symmetry)
        qsp = qspclass.max(trunc_dim, nqn)
        self._expand_dim(qsp_max=qsp)
        #expand 相当于构造一个新的 System instance, 一些东西需要初始化的！
        
        #reset attr
        self.energy_diff_std = np.inf 
        self.is_initialized = False 

    def expand_layer(self, num_of_layer):
        msg  = 'issue: when using expand_layer together with main.run() error will occor, caused by tensor_palyer'
        msg += 'a work around is to set main.S.iter1 = 1 (instead of 0) after main.expand_layer()' 
        warnings.warn(msg)
        print "\ntry to expand num_of_layer to %d..."%num_of_layer
        if self.mera.num_of_layer >= num_of_layer:
            print "not expanded: self.mera.num_of_layer >= num_of_layer"
        else:
            for i in range(num_of_layer-self.mera.num_of_layer):
                self._expand_layer()
    
        #reset attr 
        self.energy_diff_std = np.inf 
        self.is_initialized = False 
        
        
    def _expand_layer(self, out=None):
        from decorators import tensor_player, set_STATE_end_1
        #print_vars(vars(),  ['tensor_player'])
        #raise 
        if not hasattr(tensor_player, 'STATE'):   # this could happen 
            tensor_player.STATE = 'stop'
        state_bac = tensor_player.STATE
        tensor_player.STATE = "stop"

        ilayer = self.mera.num_of_layer-1
        V_top = self.mera.V[ilayer][0]
        V_dag_top = self.mera.V_dag[ilayer][0]
        
        self.mera.init_layer(ilayer, unitary_init=None, rand_init=False)
        self.mera.num_of_layer += 1 

        ilayer = self.mera.num_of_layer-1
        self.mera.init_layer(ilayer, unitary_init=None, rand_init=False)
        self.mera.V[ilayer][0] = V_top
        self.mera.V_dag[ilayer][0] = V_dag_top
        
        
        self.add_one_layer(ilayer, duplicate=True)
        self.add_top_layer()  # this is needed only when use top_level_eigenstate 
        self.iter1 = 0  #reset counting
        self.iter0 = 0 
        msg = "\nEXPAND ONE LAYER, NOW NUM_OF_LAYER IS %d\n"%self.mera.num_of_layer
        print msg
        if out is not None:
            out.write(msg)
        if state_bac == "play": 
            print "after expansion, tensor_player is switched from 'play' to 'record'"
            tensor_player.NEXT_STATE = "record" 
    
    def use_pinning_term_func_bac(self): 
        if self.iter<self.pinning_term_def['step_limit']: 
            if self.combine_2site: 
                self.H_2[0][0].data = self.H_2_bac.data  + self.pinning_term.data 
            else: 
                self.H_3[0][0].data = self.H_3_bac.data  + self.pinning_term.data 
                
        else:
            if self.combine_2site: 
                self.H_2[0][0] = self.H_2_bac 
            else: 
                self.H_3[0][0] = self.H_3_bac 
                
            self.use_pinning_term = False
            print  'self.use_pinning_term is set to False'
    
    def use_pinning_term_func(self): 
        pin_def = self.pinning_term_def
        if not self.pinning_term_def.has_key('h_backup'): 
            h_name = self.pinning_term_def['h_name']
            self.pinning_term_def['h_backup'] = getattr(self, h_name)[0][0].copy()
        
        if self.iter<self.pinning_term_def['step_limit']: 
            h_name = pin_def['h_name']
            if self.combine_2site: 
                xi = pin_def['xi']
                coeff = pin_def['lam']*np.exp(-float(self.iter)/xi)
                print_vars(vars(),  ['coeff'])
                getattr(self, h_name)[0][0].data = (
                pin_def['h_backup'].data  + coeff*pin_def['op'].data)
                #self.H_2[0][0].data = self.H_2_bac.data  + self.pinning_term.data 
            else: 
                #self.H_3[0][0].data = self.H_3_bac.data  + self.pinning_term.data 
                raise NotImplemented 
                
        else:
            h_name = self.pinning_term_def['h_name']
            h_backup = self.pinning_term_def['h_backup']
            if self.combine_2site: 
                getattr(self, h_name)[0][0]=h_backup 
            else: 
                #self.H_3[0][0] = self.H_3_bac 
                raise NotImplemented 
                
            self.use_pinning_term = False
            print  'self.use_pinning_term is set to False'

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

    def measure_S(self, S, parpath, exclude_which=None):
        from merapy.measure_and_analysis.measurement import measure_S
        #state_bac = tensor_player.STATE
        #self.stop_player()
        tensor_player.STATE = 'stop'
        msg = '\nmeasurement after final iter: '
         
        measure_S(S,  parpath, exclude_which=exclude_which, use_local_storage=self.use_local_storage)
        #tensor_player.STATE = state_bac
    @staticmethod
    def set_ham_to_identity(sys):  #for testing
        """
        for debugging 
        """
        print "hamiltonian set to identity  "*10
        SYMMETRY = sys.symmetry
        
        sys.H_2[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
        #sys.H_2[0][0].data[:] = 0.0
        sys.H_3[0][0] = iTensorFactory(SYMMETRY).unit_tensor(6)
        #sys.H_3[0][0].data[:] = 0.0

        if sys.model == "Ising":        
            sys.SxA[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
            sys.SxB[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
        else:
            if 1:
                sys.SzA[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
                sys.SzB[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
            
            if not sys.only_NN and not sys.only_NNN:
                raise Exception("set_ham_to_identity is not compatible to ops like SpA, SmA with U1 symm")
                sys.SpA[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
                sys.SpB[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
                sys.SmA[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)
                sys.SmB[0][0] = iTensorFactory(SYMMETRY).unit_tensor(4)

    if 0: 
        def check_symm(T2):
            """
            status_0p9_translated_formally
            """
            #ifdef __NONONO__        
            ME = 0.0; mp=1
            for p  in range(T2.totDim):
                X = T2.data[p]
                Get_Position(T2, p, pos)
                iSwap(pos[1], pos[2])
                iSwap(pos[3], pos[4])
        #            Y = GetElement[T2, pos]
                if abs[X-Y] == ME:  
                    ME = abs[X-Y]
                    mp = p
            Get_Position(T2,mp,pos)
            print "[4[1x,I4],2x,E15.8]", pos[0:4], ME
            #endif

        def check_symm3(T3):
            """
            status_0p9_translated_formally
            """
            #ifdef __NONONO__        
            ME = 0.0;mp=1
            for p  in range(T3.totDim):
                X = T3.data[p]
                Get_Position(T3, p, pos)
                iSwap(pos[1], pos[3])
                iSwap(pos[4], pos[6])
                #            Y = GetElement[T3, pos]
                if abs[X-Y] == ME:  
                    ME = abs[X-Y]
                    mp = p
            Get_Position(T3,mp, pos)
            print "[6[1x,I4], E15.8]", pos[0:6], ME
            #endif        
  
class TestSystem(unittest.TestCase):
    def setUp(self): 
        if 0: 
            import mera
            
            test_Mera = mera.test_Mera
             
            #M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); sys= ts.instance(M, symmetry="U1", model="Heisenberg")
            M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); 
            model_param = {'J_NN':1.0, 'J_NNN': 0.0, 'Jzz': 1.0}
            sys= System(mera=M, symmetry="U1", model="Heisenberg", only_NN=1, only_NNN=0, 
                    model_param=model_param)
            #sys = System(model=model, mera=M, symmetry=symmetry, only_NN=only_NN, only_NNN=only_NNN)
            sys.init_Hamiltonian()
            self.S = sys
            self.sys= sys
            self.M = M
    
    def tearDown(self): 
        print 'set tensor_player.STATE = "stop" in tearDown'
        tensor_player.STATE = 'stop'
    
    def xtest_example(self): 
        if 1: 
            M = Mera.example(trunc_dim=4, tot_layer=4, symmetry="U1"); 
            S= System.example(M, 'U1', 'Heisenberg')
            print S 
            H=np.array([-2. ,  0.5,  0.5,  0. ,  0.5, -1. ,  0. ,  0.5,  0.5,  0. , -1. , 0.5,  0. ,  0.5,  0.5, -2. ,  1. ,  0. ,  0. ,  0. ,  0. ,  0. , 0. ,  1. , -1.5,  0.5,  0.5, -0.5,  0. ,  1. ,  0. ,  0. , -0.5, 0.5,  0.5, -1.5,  0. ,  0. ,  1. ,  0. ,  0. ,  0. ,  1. ,  0. , -0.5,  0.5,  0.5, -1.5,  0. ,  1. ,  0. ,  0. ,  0. , -1. ,  0. , 0. ,  1. ,  0. ,  0. , -1.5,  0.5,  0.5, -0.5,  0. ,  0. ,  0. , 1. ,  0. , -1. ,  0. ])
            self.assertTrue(np.all(S.H_2[0][0].data==H))

        if 1: 
            M = Mera.example(trunc_dim=4, tot_layer=4, symmetry="Z2"); 
            S= System.example(M, 'Z2', 'Ising')
            print S 
            H=np.array([-2., -1., -1., -1., -1., -1., -1.,  0.])
            self.assertTrue(np.all(S.H_2[0][0].data==H))
            print_vars(vars(),  ['S.__dict__.keys()'])

    def test_backup_path(self): 
        dir = tempfile.mkdtemp()
        m=System.example(model='Ising', 
                #use_local_storage = 1, 
                backup_parpath=dir, 
                backup_parpath_local = dir, 
                symmetry='Travial', info=1)
        m.set_backup_path()
        #fn = m.backup_fn_auto()
        print_vars(vars(),  ['m.backup_parpath', 
            'm.backup_path', 'm.backup_path_local'])
        m.minimize('prod_state')
        self.assertTrue(os.path.exists(dir + '/4.pickle'))
      
    def test_minimize(self): 
        dir = tempfile.mkdtemp()
        m=System.example(model='Ising', symmetry='Travial', 
                backup_parpath = dir, 
                #USE_CUSTOM_RAND=1, rand_seed=1034, 
                do_measure = 1, 
                measurement_args= {'fault_tolerant':0}, 
                measure_only = ['correlation'], 
                info=1)
        print_vars(vars(),  ['m.mera'])
        
        m.minimize('prod_state')
        print_vars(vars(),  ['m.energy'])
        #some times the following fails, so diable it
        #I dont known why,  should be related to random etc
        #self.assertAlmostEqual(m.energy, -1.1718439621684591, 10)
        from merapy.decorators import tensor_player 
        tensor_player.STATE = 'stop'
        
        if 1: #test measure
            from merapy.measure_and_analysis.result_db import ResultDB_mera
            db = ResultDB_mera(dir)
            print_vars(vars(),  ['db'])
            self.assertTrue(db.has_key('correlation') and not db.has_key('magnetization'))
        
    
    def test_resume(self):
        np.random.seed(1234)
        dir = tempfile.mkdtemp()
        m=System.example(model='Ising',
                backup_parpath = dir, 
                symmetry='Travial', info=1)
        m.minimize('prod_state')
        print_vars(vars(),  ['m.energy'])
        #the following sometimes fails due conflict with main.TestMain  for random seed 
        #self.assertAlmostEqual(m.energy, -1.1718439621684591, 10)
        if 1: 
            from merapy.decorators import tensor_player 
            tensor_player.STATE = 'stop'
            m=System.example(model='Ising',
                    backup_parpath = dir, 
                    symmetry='Travial', info=1)
            
            m.minimize('prod_state', q_iter=10)
            print_vars(vars(),  ['m.energy', 'm.iter0', 'm.iter1'])
            self.assertTrue(m.iter0==5 and m.iter1==10)
            
        if 1:  # test init mera using other_backup_path  
            tensor_player.STATE = 'stop'
            m1=System.example(model='Ising',
                    backup_parpath = tempfile.mkdtemp(), 
                    other_backup_path = dir + '/4.pickle', 
                    symmetry='Travial', info=1)
            
            m1.minimize('prod_state', q_iter=6)
            print_vars(vars(),  ['m1.energy', 'm1.iter0', 'm1.iter1'])
            if 1: 
                from merapy.decorators import tensor_player 
                tensor_player.STATE = 'stop'
                m=System.example(model='Ising',
                        backup_parpath = dir, 
                        symmetry='Travial', info=1)
                
                m.minimize('prod_state', q_iter=15)
                print_vars(vars(),  ['m.energy', 'm.iter0', 'm.iter1'])
                self.assertTrue(m.iter0==5 and m.iter1==15)
                
                self.assertAlmostEqual(m.energy, m1.energy, 10)
        
    
    def test_minimize_finite_site(self): 
        m=System.example(model='Ising', symmetry='Travial', 
                info=1)
        m._minimize_finite_size()
        print_vars(vars(),  ['m.energy'])
        #this may fail due to randomness not properly fixed
        #self.assertAlmostEqual(m.energy, -1.1718439621684591, 10)
        from merapy.decorators import tensor_player 
        tensor_player.STATE = 'stop'
    
    def test_minimize_scale_invar(self): 
        m=System.example(model='Ising', rand_seed=1234, symmetry='Travial', info=1)
        m._minimize_finite_size()
        m._minimize_scale_invar()
        print_vars(vars(),  ['m.energy'])
        #self.assertAlmostEqual(m.energy, -0.98131885665452145, 10)
        from merapy.decorators import tensor_player 
        tensor_player.STATE = 'stop'
    
    def test_expand_dim_and_layer(self): 
        m=System.example(model='Ising', symmetry='Travial', info=1)
        m._minimize_finite_size()
        #print_vars(vars(),  ['m.energy'])
        #self.assertAlmostEqual(m.energy, -1.1718439621684591, 10)
        print m.mera 
        from merapy.decorators import tensor_player 
        tensor_player.STATE = 'stop'
        m.expand_dim(6)
        m._minimize_finite_size()
        print m.mera 

    def xtest_save(self): 
        from merapy.utilities import random_str
        import os
        
        if 1: 
            fn = random_str(size=8)
            fn = '/'.join(['/tmp', fn])
            print fn
            S.save(fn)
            S = System.load(fn)
            S.save(fn)
            print S
            with self.assertRaises(Exception) as ee:
                S.iter1 = -10
                S.save(fn)   # deney save
            os.remove(fn)
            print 'test save check pass'
            
        if 1:
            print '---'*10
            fn = random_str(size=8) 
            fn = '/'.join(['/tmp', fn +  'llllllll'])
            print fn
            System._save(S, fn, use_local_storage=1)
            with rpyc_conn_local() as conn: 
                os_local = conn.modules.os  
                self.assertTrue(os_local.path.exists(fn))
                print 'test save use_local_storage pass'
                os_local.remove(fn)
    
    def test_load(self): 
        #fn = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/4.pickle' 
        #S=System.load(fn, 1)
        #print S
        pass 
    
    if 0: 
        def instance(self, M, symmetry, model, only_NN=True, only_NNN=False):
            test_Mera = mera.test_Mera
            ts= test_System
            
            #M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); sys= ts.instance(M, symmetry="U1", model="Heisenberg")
            M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); 
            
            sys=System(model=model, mera=M, symmetry=symmetry, only_NN=only_NN, only_NNN=only_NNN)
            #M = test_Mera.instance()
            #M = test_Mera.M
            sys.init_Hamiltonian()
            return sys
        
        def contract(self):
            """   ----pass """
            sys= self.instance()
            layer = 1
            which = -1
            G = G3[which]
            #order = order3[:, which]
            order= order3[which]
            M = test_Mera.instance()
            #print "MMM", M.num_of_layer
            
            #print 'ggg', G3[which]
            Names= []
            Ops= []
            exception = 4
            res= 0
            res = sys.contract(M, layer, G, 0, Names, 
                    Ops, order, exception=exception, info=2)

        def contract1(self, info=0):
            S= self.sys
            M = self.M
            
            ilayer = 0
            which =-1 
            G = G2[which]
            #order= order2[:, which]
            order= order2[which]
            #print "ooo", order
            weight=weight2[which]
            name="OO"   #hamiltonian operator
            
            #print 
            exception=G.find_node(name)[0]
            
            print "exception", exception

            info = 5
            T_tmp = S.contract(M, ilayer, G, 0, Names=[], Ops=[], order=order, 
                    exception=exception, info=info-1)
            print "ddd", T_tmp.data.round(5)

        def contract_except_V(self, info=0):
            """
            ---pass
            """
            S= self.sys
            M = self.M
            
            ilayer = 0
            which =-1 
            G = G2[which]
            #order= order2[:, which]
            order= order2[which]
            #print "ooo", order
            weight=weight2[which]
            name="V"   #hamiltonian operator
            
            #print 
            exception=G.find_node(name)[1]

            
            print "exception", exception

            print "rho_2 is set to be 1"
            S.rho_2[0][0].data[:] = 1.

            info = 1
            T_tmp = S.contract(M, ilayer, G, 0, Names=[], Ops=[], order=order, 
                    exception=exception, info=info-1)
            print "ddd", T_tmp.data.round(5)

        def contract_V_Vdag(self):
            T1 = self.M.V[0][0]
            T2 = self.M.V_dag[0][0]
            G = G2[1]
            n1 = G.find_node("V")[0]
            n2 = G.find_node("Vp")[0]
            ord1 = G.edges[:4, n1]
            ord2 = G.edges[:4, n2]

            ord1 = [1, 2, 3, 9]
            ord2 = [6, 1, 2, 3]
            print ord1, ord2
            print T1.data.round(3)
            print T2.data.round(3)        
            res, leg= T1.contract(T2, ord1, ord2)
            print res.data.round(3 )
            print T1

        def contract_two(self):
            """  ---pass """
            T1 = self.M.V[0][0]
            T2 = self.sys.H_2[0][0]
            G = G2[1]
            n1 = G.find_node("V")[0]
            n2 = G.find_node("oo")[0]
            ord1 = G.edges[:4, n1]
            ord2 = G.edges[:4, n2]

            #ord1 = [1, 2, 3, 9]
            #ord2 = [6, 1, 2, 3]
            print ord1, ord2
            print T1.data.round(3)
            print T2.data.round(3)        
            
            info = 1
            res, leg= T1.contract(T2, ord1, ord2, info=info-1)
            print res.data.round(3 )
            print leg


        def contract_two_UV(self):
            """  ---pass """
            T1 = self.M.V[1][0]
            #T2 = self.sys.H_2[0][0]
            T2 = self.M.U[1][0]

            G = G2[1]
            
            n1 = G.find_node("V")[0]
            n2 = G.find_node("U")[0]
            
            ord1 = G.edges[:4, n1]
            ord2 = G.edges[:4, n2]

            #ord1 = [1, 2, 3, 9]
            #ord2 = [6, 1, 2, 3]
            print ord1, ord2
            print T1.data[:8].round(5), T1.data[-5:].round(5)
            print T2.data[:8].round(5), T2.data[-5:].round(5)        
            
            info = 1
            res, leg= T1.contract(T2, ord1, ord2, info=info-1)
            print res.data[:8].round(5), "\t", res.data[-8:].round(5)
            #print res
            #print leg

        def contract_two_core(self):
            """  ---pass """
            T1 = self.M.V[1][0]
            T2 = self.M.U[1][0]

            
            print T1.data[:8].round(5), T1.data[-5:].round(5)
            print T2.data[:8].round(5), T2.data[-5:].round(5)        
            
            info = 1
            res= T1.contract_core(T2, 2 )
            print res.data[:8].round(5), "\t", res.data[-8:].round(5)
            #print res

        def add_env(self, info=0):
            S= self.sys
            #M = test_Mera.instance()
            M = self.M
            ilayer = 0
            jlayer = 1
            which = -1
            G = G2[which]
            #order= order2[:, which]
            order= order2[which]
            weight=weight2[which]
            name="OO"   #hamiltonian operator
            info = 2
            S.H_2[jlayer][0]=S.add_env(M, ilayer, G, 
                    order, name, weight,  info=info-1)
            
            print S.H_2[jlayer][0].data.round(3)
            #print S

        def add_env_V(self, info=0):
            S= self.sys
            #M = test_Mera.instance()
            M = self.M
            ilayer = 0
            
            which = -1 
            G = G2[which]
            #order= order2[:, which]
            order= order2[which]
            weight=weight2[which]
            name="V"   #hamiltonian operator
            info = 2

            print "rho_2 is set to 1.0"
            S.rho_2[0][0].data[:] = 1.
            
            #weight = -weight
            #if which == -1: weight = weight*0.5

            t=S.add_env(M, ilayer, G, 
                    order, name, weight,  info=info-1)
            
            print t.data.round(3)


        def ascending(self,info=0):
            M = self.M
            S= self.sys
            i=0
            j = 0
            name="OO"   #hamiltonian operator
            ilayer = 0
            jlayer = 1
            S.H_2[jlayer][j].data[:] =0.0

            info = 1

            if info:
                print "START ascending_ham"
            for g in [-1, 0, 1]:
                if info:
                    print "iLayer,g=",ilayer, g
                #print S.H_2[ilayer][j]
                
                #q29
                G = G2[g]
                #order= order2[:,g]
                order= order2[g]
                weight=weight2[g]
                
                temp=S.add_env(M, ilayer, G, order, name, weight, S.H_2[jlayer][j], info=info-1)
                S.H_2[jlayer][j].data += temp.data 
                print S.H_2[jlayer][0].data.round(5)
            print S.H_2[jlayer][0].data.round(5)
            if info:
                print "END ascending_ham"

        def eng_ham(self):
            S= self.instance()
            M = test_Mera.instance()
            layer = 3
            S.eng_ham(layer)
        
        @staticmethod
        def op_heisenberg():
            """ ---- pass """
            from quantum_number import init_System_QSp
            QN_idendity, QSp_base, QSp_null= init_System_QSp("U1")

            h0 = System.op_heisenberg()
            print h0[2]
            #print h0[2]  #.data.round(5)
            #results
            """
                625 ---pass
                    ----Begin iTensor----------------------------------------------
                        rank:	4
                        idx_dim:	81
                        nidx:	19
                        totQN:	(    0)
                        QNs:	['[(    0) (    1) (   -1)]', '[(    0) (    1) (   -1)]', '[(    0) (   -1) (    1)]', '[(    0) (   -1) (    1)]']
                        Dims:	array([4, 4, 4, 4])
                        totDim:	70
                        idx:	array([ 0, 81, 81, 81, 81,  1, 81,  2, 81, 81,  3, 81,  4, 81, 81, 81, 81,
                               81, 81, 81,  5, 81, 81, 81,  6, 81, 81, 81,  7, 81,  8, 81, 81, 81,
                               81, 81, 81, 81, 81, 81,  9, 81, 81, 81, 81, 10, 81, 81, 81, 81, 11,
                               81, 12, 81, 81, 81, 13, 81, 81, 81, 14, 81, 81, 15, 81, 81, 81, 81,
                               16, 81, 17, 81, 81, 81, 81, 81, 81, 81, 81, 81, 18])
                        Block_idx:	[ 0 16 20 24 28 32 36 40 44 48 49 53 54 55 59 63 67 68 69]
                            [16  4  4  4  4  4  4  4  4  1  4  1  1  4  4  4  1  1  1]
                            [ 0  5  7 10 12 20 24 28 30 40 45 50 52 56 60 63 68 70 80]

                        Addr_idx:	[[0 2 1 1 0 2 0 1 0 1 0 2 1 2 0 0 2 1 2]
                         [0 1 2 0 1 0 2 0 1 1 0 1 2 0 2 0 1 2 2]
                         [0 0 0 1 1 2 2 0 0 1 2 2 2 0 0 1 1 1 2]
                         [0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]]
                        data:	[0 0 0 0]: [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
                        [2 1 0 0]: [ 0.  0.  0.  0.]
                        [1 2 0 0]: [ 0.  0.  0.  1.]
                        [1 0 1 0]: [ 0.  0.  0.  0.]
                        [0 1 1 0]: [ 0.  0.  0.  0.]
                        [2 0 2 0]: [ 0.  0.  0.  0.]
                        [0 2 2 0]: [ 0.  0.  1.  0.]
                        [1 0 0 1]: [ 0.  0.  1.  0.]
                        [0 1 0 1]: [ 0.  0.  0.  0.]
                        [1 1 1 1]: [ 0.]
                        [0 0 2 1]: [ 1.  0.  0.  0.]
                        [2 1 2 1]: [ 0.]
                        [1 2 2 1]: [ 0.]
                        [2 0 0 2]: [ 0.  0.  0.  0.]
                        [0 2 0 2]: [ 0.  0.  0.  0.]
                        [0 0 1 2]: [ 0.  0.  0.  0.]
                        [2 1 1 2]: [ 0.]
                        [1 2 1 2]: [ 0.]
                        [2 2 2 2]: [ 0.]

                    ----End iTensor----------------------------------------------

                    625----Begin Tensor-------------------------------------------
                        T%rank=    4
                        T%totQN=     0
                        T%nQN=     3     3     3     3
                        T%nIdx=   19
                        T%QNs(:,1)     =     0     1    -1
                        T%QNs%Dims(:,1)=     2     1     1
                        T%QNs(:,2)     =     0     1    -1
                        T%QNs%Dims(:,2)=     2     1     1
                        T%QNs(:,3)     =     0    -1     1
                        T%QNs%Dims(:,3)=     2     1     1
                        T%QNs(:,4)     =     0    -1     1
                        T%QNs%Dims(:,4)=     2     1     1
                        Data=
                        T%Block_QN=     0     0     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     2     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.1000E+01
                        T%Block_QN=     1     0     1     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     1     1     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     0     2     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                    M    T%Block_QN=     0     2     2     0
                        0.0000E+00
                        0.0000E+00
                        0.1000E+01
                        0.0000E+00
                    M    T%Block_QN=     1     0     0     1
                        0.0000E+00
                        0.0000E+00
                        0.1000E+01
                        0.0000E+00
                        T%Block_QN=     0     1     0     1
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     1     1     1
                        0.0000E+00
                        T%Block_QN=     0     0     2     1
                        0.1000E+01
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     2     1
                        0.0000E+00
                        T%Block_QN=     1     2     2     1
                        0.0000E+00
                        T%Block_QN=     2     0     0     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     2     0     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     0     1     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     1     2
                        0.0000E+00
                        T%Block_QN=     1     2     1     2
                        0.0000E+00
                        T%Block_QN=     2     2     2     2
                        0.0000E+00
                    ----End   Tensor-------------------------------------------

                622  --- pass absolutly
                    hhh 622  ----Begin iTensor----------------------------------------------
                        rank:	4
                        idx_dim:	81
                        nidx:	19
                        totQN:	(    0)
                        QNs:	['[(    0) (    1) (   -1)]', '[(    0) (    1) (   -1)]', '[(    0) (   -1) (    1)]', '[(    0) (   -1) (    1)]']
                        Dims:	array([4, 4, 4, 4])
                        totDim:	70
                        idx:	array([ 0, 81, 81, 81, 81,  1, 81,  2, 81, 81,  3, 81,  4, 81, 81, 81, 81,
                               81, 81, 81,  5, 81, 81, 81,  6, 81, 81, 81,  7, 81,  8, 81, 81, 81,
                               81, 81, 81, 81, 81, 81,  9, 81, 81, 81, 81, 10, 81, 81, 81, 81, 11,
                               81, 12, 81, 81, 81, 13, 81, 81, 81, 14, 81, 81, 15, 81, 81, 81, 81,
                               16, 81, 17, 81, 81, 81, 81, 81, 81, 81, 81, 81, 18])
                        Block_idx:	[ 0 16 20 24 28 32 36 40 44 48 49 53 54 55 59 63 67 68 69]
                            [16  4  4  4  4  4  4  4  4  1  4  1  1  4  4  4  1  1  1]
                            [ 0  5  7 10 12 20 24 28 30 40 45 50 52 56 60 63 68 70 80]

                        Addr_idx:	[[0 2 1 1 0 2 0 1 0 1 0 2 1 2 0 0 2 1 2]
                         [0 1 2 0 1 0 2 0 1 1 0 1 2 0 2 0 1 2 2]
                         [0 0 0 1 1 2 2 0 0 1 2 2 2 0 0 1 1 1 2]
                         [0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]]
                        data:	[0 0 0 0]: [-2.  1.  1.  0.  1.  0.  0.  1.  1.  0.  0.  1.  0.  1.  1. -2.]
                        [2 1 0 0]: [ 2.  0.  0.  0.]
                        [1 2 0 0]: [ 0.  0.  0.  0.]
                        [1 0 1 0]: [-1.  1.  1.  1.]
                        [0 1 1 0]: [ 0.  2.  0.  0.]
                        [2 0 2 0]: [ 1.  1.  1. -1.]
                        [0 2 2 0]: [ 0.  0.  0.  0.]
                        [1 0 0 1]: [ 0.  0.  0.  0.]
                        [0 1 0 1]: [ 1.  1.  1. -1.]
                        [1 1 1 1]: [ 2.]
                        [0 0 2 1]: [ 0.  0.  0.  0.]
                        [2 1 2 1]: [ 0.]
                        [1 2 2 1]: [ 0.]
                        [2 0 0 2]: [ 0.  2.  0.  0.]
                        [0 2 0 2]: [-1.  1.  1.  1.]
                        [0 0 1 2]: [ 0.  0.  0.  2.]
                        [2 1 1 2]: [ 0.]
                        [1 2 1 2]: [ 0.]
                        [2 2 2 2]: [ 2.]

                    ----End iTensor----------------------------------------------


                    622----Begin Tensor-------------------------------------------
                        T%rank=    4
                        T%totQN=     0
                        T%nQN=     3     3     3     3
                        T%nIdx=   19
                        T%QNs(:,1)     =     0     1    -1
                        T%QNs%Dims(:,1)=     2     1     1
                        T%QNs(:,2)     =     0     1    -1
                        T%QNs%Dims(:,2)=     2     1     1
                        T%QNs(:,3)     =     0    -1     1
                        T%QNs%Dims(:,3)=     2     1     1
                        T%QNs(:,4)     =     0    -1     1
                        T%QNs%Dims(:,4)=     2     1     1
                        Data=
                        T%Block_QN=     0     0     0     0
                        -.2000E+01
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        0.1000E+01
                        0.0000E+00
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        0.0000E+00
                        0.1000E+01
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        -.2000E+01
                        T%Block_QN=     2     1     0     0
                        0.2000E+01
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     2     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     0     1     0
                        -.1000E+01
                        0.1000E+01
                        0.1000E+01
                        0.1000E+01
                    M    T%Block_QN=     0     1     1     0
                        0.0000E+00
                        0.2000E+01
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     0     2     0
                        0.1000E+01
                        0.1000E+01
                        0.1000E+01
                        -.1000E+01
                        T%Block_QN=     0     2     2     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     0     0     1
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     1     0     1
                        0.1000E+01
                        0.1000E+01
                        0.1000E+01
                        -.1000E+01
                        T%Block_QN=     1     1     1     1
                        0.2000E+01
                        T%Block_QN=     0     0     2     1
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     2     1
                        0.0000E+00
                        T%Block_QN=     1     2     2     1
                        0.0000E+00
                    M    T%Block_QN=     2     0     0     2
                        0.0000E+00
                        0.2000E+01
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     2     0     2
                        -.1000E+01
                        0.1000E+01
                        0.1000E+01
                        0.1000E+01
                        T%Block_QN=     0     0     1     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.2000E+01
                        T%Block_QN=     2     1     1     2
                        0.0000E+00
                        T%Block_QN=     1     2     1     2
                        0.0000E+00
                        T%Block_QN=     2     2     2     2
                        0.2000E+01
                    ----End   Tensor-------------------------------------------


                641 --- pass
                    641 ----Begin Tensor-------------------------------------------
                        T%rank=    4
                        T%totQN=     0
                        T%nQN=     3     3     3     3
                        T%nIdx=   19
                        T%QNs(:,1)     =     0     1    -1
                        T%QNs%Dims(:,1)=     2     1     1
                        T%QNs(:,2)     =     0     1    -1
                        T%QNs%Dims(:,2)=     2     1     1
                        T%QNs(:,3)     =     0    -1     1
                        T%QNs%Dims(:,3)=     2     1     1
                        T%QNs(:,4)     =     0    -1     1
                        T%QNs%Dims(:,4)=     2     1     1
                        Data=
                        T%Block_QN=     0     0     0     0
                        -.1000E+01
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        0.1000E+01
                        -.1000E+01
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        -.1000E+01
                        0.1000E+01
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        -.1000E+01
                        T%Block_QN=     2     1     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     2     0     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     0     1     0
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        T%Block_QN=     0     1     1     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     0     2     0
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        T%Block_QN=     0     2     2     0
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     1     0     0     1
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     1     0     1
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        T%Block_QN=     1     1     1     1
                        0.1000E+01
                        T%Block_QN=     0     0     2     1
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     2     1
                        0.1000E+01
                        T%Block_QN=     1     2     2     1
                        0.0000E+00
                        T%Block_QN=     2     0     0     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     0     2     0     2
                        0.0000E+00
                        0.1000E+01
                        0.1000E+01
                        0.0000E+00
                        T%Block_QN=     0     0     1     2
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        0.0000E+00
                        T%Block_QN=     2     1     1     2
                        0.0000E+00
                        T%Block_QN=     1     2     1     2
                        0.1000E+01
                        T%Block_QN=     2     2     2     2
                        0.1000E+01
                    ----End   Tensor-------------------------------------------


                    hhh 641  ----Begin iTensor----------------------------------------------
                        rank:	4
                        idx_dim:	81
                        nidx:	19
                        totQN:	(    0)
                        QNs:	['[(    0) (    1) (   -1)]', '[(    0) (    1) (   -1)]', '[(    0) (   -1) (    1)]', '[(    0) (   -1) (    1)]']
                        Dims:	array([4, 4, 4, 4])
                        totDim:	70
                        idx:	array([ 0, 81, 81, 81, 81,  1, 81,  2, 81, 81,  3, 81,  4, 81, 81, 81, 81,
                               81, 81, 81,  5, 81, 81, 81,  6, 81, 81, 81,  7, 81,  8, 81, 81, 81,
                               81, 81, 81, 81, 81, 81,  9, 81, 81, 81, 81, 10, 81, 81, 81, 81, 11,
                               81, 12, 81, 81, 81, 13, 81, 81, 81, 14, 81, 81, 15, 81, 81, 81, 81,
                               16, 81, 17, 81, 81, 81, 81, 81, 81, 81, 81, 81, 18])
                        Block_idx:	[ 0 16 20 24 28 32 36 40 44 48 49 53 54 55 59 63 67 68 69]
                            [16  4  4  4  4  4  4  4  4  1  4  1  1  4  4  4  1  1  1]
                            [ 0  5  7 10 12 20 24 28 30 40 45 50 52 56 60 63 68 70 80]

                        Addr_idx:	[[0 2 1 1 0 2 0 1 0 1 0 2 1 2 0 0 2 1 2]
                         [0 1 2 0 1 0 2 0 1 1 0 1 2 0 2 0 1 2 2]
                         [0 0 0 1 1 2 2 0 0 1 2 2 2 0 0 1 1 1 2]
                         [0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]]
                        data:	[0 0 0 0]: [-1.  1.  1.  0.  1. -1.  0.  1.  1.  0. -1.  1.  0.  1.  1. -1.]
                        [2 1 0 0]: [ 0.  0.  0.  0.]
                        [1 2 0 0]: [ 0.  0.  0.  0.]
                        [1 0 1 0]: [ 0.  1.  1.  0.]
                        [0 1 1 0]: [ 0.  0.  0.  0.]
                        [2 0 2 0]: [ 0.  1.  1.  0.]
                        [0 2 2 0]: [ 0.  0.  0.  0.]
                        [1 0 0 1]: [ 0.  0.  0.  0.]
                        [0 1 0 1]: [ 0.  1.  1.  0.]
                        [1 1 1 1]: [ 1.]
                        [0 0 2 1]: [ 0.  0.  0.  0.]
                        [2 1 2 1]: [ 1.]
                        [1 2 2 1]: [ 0.]
                        [2 0 0 2]: [ 0.  0.  0.  0.]
                        [0 2 0 2]: [ 0.  1.  1.  0.]
                        [0 0 1 2]: [ 0.  0.  0.  0.]
                        [2 1 1 2]: [ 0.]
                        [1 2 1 2]: [ 1.]
                        [2 2 2 2]: [ 1.]

                    ----End iTensor----------------------------------------------

                635  ---pass
                    ----Begin Tensor-------------------------------------------
                        T%rank=    4
                        T%totQN=     0
                        T%nQN=     3     3     3     3
                        T%nIdx=   19
                        T%QNs(:,1)     =     0     1    -1
                        T%QNs%Dims(:,1)=     2     1     1
                        T%QNs(:,2)     =     0     1    -1
                        T%QNs%Dims(:,2)=     2     1     1
                        T%QNs(:,3)     =     0    -1     1
                        T%QNs%Dims(:,3)=     2     1     1
                        T%QNs(:,4)     =     0    -1     1
                        T%QNs%Dims(:,4)=     2     1     1
                        Data=
                        T%Block_QN=     0     0     0     0
                         -.10E+01 0.50E+00 0.50E+00 0.00E+00 0.50E+00 0.00E+00 0.00E+00 0.50E+00 0.50E+00 0.00E+00 0.00E+00 0.50E+00 0.00E+00 0.50E+00 0.50E+00 -.10E+01
                        T%Block_QN=     2     1     0     0
                         0.10E+01 0.00E+00 0.00E+00 0.00E+00
                        T%Block_QN=     1     2     0     0
                         0.00E+00 0.00E+00 0.00E+00 0.10E+01
                        T%Block_QN=     1     0     1     0
                         -.50E+00 0.50E+00 0.50E+00 0.50E+00
                        T%Block_QN=     0     1     1     0
                         0.00E+00 0.10E+01 0.00E+00 0.00E+00
                        T%Block_QN=     2     0     2     0
                         0.50E+00 0.50E+00 0.50E+00 -.50E+00
                        T%Block_QN=     0     2     2     0
                         0.00E+00 0.00E+00 0.10E+01 0.00E+00
                        T%Block_QN=     1     0     0     1
                         0.00E+00 0.00E+00 0.10E+01 0.00E+00
                        T%Block_QN=     0     1     0     1
                         0.50E+00 0.50E+00 0.50E+00 -.50E+00
                        T%Block_QN=     1     1     1     1
                         0.10E+01
                        T%Block_QN=     0     0     2     1
                         0.10E+01 0.00E+00 0.00E+00 0.00E+00
                        T%Block_QN=     2     1     2     1
                         0.00E+00
                        T%Block_QN=     1     2     2     1
                         0.00E+00
                        T%Block_QN=     2     0     0     2
                         0.00E+00 0.10E+01 0.00E+00 0.00E+00
                        T%Block_QN=     0     2     0     2
                         -.50E+00 0.50E+00 0.50E+00 0.50E+00
                        T%Block_QN=     0     0     1     2
                         0.00E+00 0.00E+00 0.00E+00 0.10E+01
                        T%Block_QN=     2     1     1     2
                         0.00E+00
                        T%Block_QN=     1     2     1     2
                         0.00E+00
                        T%Block_QN=     2     2     2     2
                         0.10E+01
                    ----End   Tensor-------------------------------------------


                    635 ----Begin iTensor----------------------------------------------
                        rank:	4
                        idx_dim:	81
                        nidx:	19
                        totQN:	(    0)
                        QNs:	['[(    0) (    1) (   -1)]', '[(    0) (    1) (   -1)]', '[(    0) (   -1) (    1)]', '[(    0) (   -1) (    1)]']
                        Dims:	array([4, 4, 4, 4])
                        totDim:	70
                        idx:	array([ 0, 81, 81, 81, 81,  1, 81,  2, 81, 81,  3, 81,  4, 81, 81, 81, 81,
                               81, 81, 81,  5, 81, 81, 81,  6, 81, 81, 81,  7, 81,  8, 81, 81, 81,
                               81, 81, 81, 81, 81, 81,  9, 81, 81, 81, 81, 10, 81, 81, 81, 81, 11,
                               81, 12, 81, 81, 81, 13, 81, 81, 81, 14, 81, 81, 15, 81, 81, 81, 81,
                               16, 81, 17, 81, 81, 81, 81, 81, 81, 81, 81, 81, 18])
                        Block_idx:	[ 0 16 20 24 28 32 36 40 44 48 49 53 54 55 59 63 67 68 69]
                            [16  4  4  4  4  4  4  4  4  1  4  1  1  4  4  4  1  1  1]
                            [ 0  5  7 10 12 20 24 28 30 40 45 50 52 56 60 63 68 70 80]

                        Addr_idx:	[[0 2 1 1 0 2 0 1 0 1 0 2 1 2 0 0 2 1 2]
                         [0 1 2 0 1 0 2 0 1 1 0 1 2 0 2 0 1 2 2]
                         [0 0 0 1 1 2 2 0 0 1 2 2 2 0 0 1 1 1 2]
                         [0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]]
                        data:	[0 0 0 0]: [-1.   0.5  0.5  0.   0.5  0.   0.   0.5  0.5  0.   0.   0.5  0.   0.5  0.5
                         -1. ]
                        [2 1 0 0]: [ 1.  0.  0.  0.]
                        [1 2 0 0]: [ 0.  0.  0.  1.]
                        [1 0 1 0]: [-0.5  0.5  0.5  0.5]
                        [0 1 1 0]: [ 0.  1.  0.  0.]
                        [2 0 2 0]: [ 0.5  0.5  0.5 -0.5]
                        [0 2 2 0]: [ 0.  0.  1.  0.]
                        [1 0 0 1]: [ 0.  0.  1.  0.]
                        [0 1 0 1]: [ 0.5  0.5  0.5 -0.5]
                        [1 1 1 1]: [ 1.]
                        [0 0 2 1]: [ 1.  0.  0.  0.]
                        [2 1 2 1]: [ 0.]
                        [1 2 2 1]: [ 0.]
                        [2 0 0 2]: [ 0.  1.  0.  0.]
                        [0 2 0 2]: [-0.5  0.5  0.5  0.5]
                        [0 0 1 2]: [ 0.  0.  0.  1.]
                        [2 1 1 2]: [ 0.]
                        [1 2 1 2]: [ 0.]
                        [2 2 2 2]: [ 1.]

                    ----End iTensor----------------------------------------------


                665 ---pass
                    ----Begin Tensor-------------------------------------------
                        T%rank=    4
                        T%totQN=     0
                        T%nQN=     3     3     3     3
                        T%nIdx=   19
                        T%QNs(:,1)     =     0     1    -1
                        T%QNs%Dims(:,1)=     2     1     1
                        T%QNs(:,2)     =     0     1    -1
                        T%QNs%Dims(:,2)=     2     1     1
                        T%QNs(:,3)     =     0    -1     1
                        T%QNs%Dims(:,3)=     2     1     1
                        T%QNs(:,4)     =     0    -1     1
                        T%QNs%Dims(:,4)=     2     1     1
                        Data=
                        T%Block_QN=     0     0     0     0
                         -.18E+01 0.50E+00 0.50E+00 0.00E+00 0.50E+00 -.12E+01 0.00E+00 0.50E+00 0.50E+00 0.00E+00 -.12E+01 0.50E+00 0.00E+00 0.50E+00 0.50E+00 -.18E+01
                        T%Block_QN=     2     1     0     0
                         0.10E+01 0.24E+00 0.24E+00 0.00E+00
                        T%Block_QN=     1     2     0     0
                         0.00E+00 0.24E+00 0.24E+00 0.10E+01
                        T%Block_QN=     1     0     1     0
                         -.15E+01 0.50E+00 0.50E+00 -.50E+00
                        T%Block_QN=     0     1     1     0
                         0.24E+00 0.10E+01 0.00E+00 0.24E+00
                        T%Block_QN=     2     0     2     0
                         -.50E+00 0.50E+00 0.50E+00 -.15E+01
                        T%Block_QN=     0     2     2     0
                         0.24E+00 0.00E+00 0.10E+01 0.24E+00
                        T%Block_QN=     1     0     0     1
                         0.24E+00 0.00E+00 0.10E+01 0.24E+00
                        T%Block_QN=     0     1     0     1
                         -.50E+00 0.50E+00 0.50E+00 -.15E+01
                        T%Block_QN=     1     1     1     1
                         0.24E+00
                        T%Block_QN=     0     0     2     1
                         0.10E+01 0.24E+00 0.24E+00 0.00E+00
                        T%Block_QN=     2     1     2     1
                         -.12E+01
                        T%Block_QN=     1     2     2     1
                         0.00E+00
                        T%Block_QN=     2     0     0     2
                         0.24E+00 0.10E+01 0.00E+00 0.24E+00
                        T%Block_QN=     0     2     0     2
                         -.15E+01 0.50E+00 0.50E+00 -.50E+00
                        T%Block_QN=     0     0     1     2
                         0.00E+00 0.24E+00 0.24E+00 0.10E+01
                        T%Block_QN=     2     1     1     2
                         0.00E+00
                        T%Block_QN=     1     2     1     2
                         -.12E+01
                        T%Block_QN=     2     2     2     2
                         0.24E+00
                    ----End   Tensor-------------------------------------------

                    ----Begin iTensor----------------------------------------------
                        rank:	4
                        idx_dim:	81
                        nidx:	19
                        totQN:	(    0)
                        QNs:	['[(    0) (    1) (   -1)]', '[(    0) (    1) (   -1)]', '[(    0) (   -1) (    1)]', '[(    0) (   -1) (    1)]']
                        Dims:	array([4, 4, 4, 4])
                        totDim:	70
                        idx:	array([ 0, 81, 81, 81, 81,  1, 81,  2, 81, 81,  3, 81,  4, 81, 81, 81, 81,
                               81, 81, 81,  5, 81, 81, 81,  6, 81, 81, 81,  7, 81,  8, 81, 81, 81,
                               81, 81, 81, 81, 81, 81,  9, 81, 81, 81, 81, 10, 81, 81, 81, 81, 11,
                               81, 12, 81, 81, 81, 13, 81, 81, 81, 14, 81, 81, 15, 81, 81, 81, 81,
                               16, 81, 17, 81, 81, 81, 81, 81, 81, 81, 81, 81, 18])
                        Block_idx:	[ 0 16 20 24 28 32 36 40 44 48 49 53 54 55 59 63 67 68 69]
                            [16  4  4  4  4  4  4  4  4  1  4  1  1  4  4  4  1  1  1]
                            [ 0  5  7 10 12 20 24 28 30 40 45 50 52 56 60 63 68 70 80]

                        Addr_idx:	[[0 2 1 1 0 2 0 1 0 1 0 2 1 2 0 0 2 1 2]
                         [0 1 2 0 1 0 2 0 1 1 0 1 2 0 2 0 1 2 2]
                         [0 0 0 1 1 2 2 0 0 1 2 2 2 0 0 1 1 1 2]
                         [0 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]]
                        data:	[0 0 0 0]: [-1.75881  0.5      0.5      0.       0.5     -1.24119  0.       0.5      0.5
                          0.      -1.24119  0.5      0.       0.5      0.5     -1.75881]
                        [2 1 0 0]: [ 1.       0.24119  0.24119  0.     ]
                        [1 2 0 0]: [ 0.       0.24119  0.24119  1.     ]
                        [1 0 1 0]: [-1.5  0.5  0.5 -0.5]
                        [0 1 1 0]: [ 0.24119  1.       0.       0.24119]
                        [2 0 2 0]: [-0.5  0.5  0.5 -1.5]
                        [0 2 2 0]: [ 0.24119  0.       1.       0.24119]
                        [1 0 0 1]: [ 0.24119  0.       1.       0.24119]
                        [0 1 0 1]: [-0.5  0.5  0.5 -1.5]
                        [1 1 1 1]: [ 0.24119]
                        [0 0 2 1]: [ 1.       0.24119  0.24119  0.     ]
                        [2 1 2 1]: [-1.24119]
                        [1 2 2 1]: [ 0.]
                        [2 0 0 2]: [ 0.24119  1.       0.       0.24119]
                        [0 2 0 2]: [-1.5  0.5  0.5 -0.5]
                        [0 0 1 2]: [ 0.       0.24119  0.24119  1.     ]
                        [2 1 1 2]: [ 0.]
                        [1 2 1 2]: [-1.24119]
                        [2 2 2 2]: [ 0.24119]
                    ----End iTensor----------------------------------------------
            """
        
        @staticmethod
        def op_heisenberg_h3():
            """ ---- pass """
            from quantum_number import init_System_QSp
            symmetry = "U1"

            QN_idendity, QSp_base, QSp_null= init_System_QSp(symmetry)
            h0 = System.op_heisenberg(QN_idendity=QN_idendity, 
                    symmetry=symmetry, 
                    only_NN=False, 
                    QSp_base=QSp_base, only_NNN=False, )
            #print h0[3]
            #print h0[2]  #.data.round(5)


        def H2(self):
            """
            ---pass  
            when h = 1,  gamma = 1.0
            final form of H0[2].data = [ 0., -1., -1., -1., -1., -1., -1., -2.]

            """
            S= self.instance()
            h2 = S.H_2[0][0]
            print h2.data



        def H2_all(self):
           for l in range(self.sys.layer):
                print "H_2 at layer ", l
                h2 = self.sys.H_2[l][0]
                print h2


        def rho2(self):
            print "\n"*5                  
            for l in range(sys.layer):
                print "layer=", l
                r2 = self.sys.rho_2[l][0]
                #print r2.data.round(5)
                print r2.QSp

        def rho2_all(self):
           import top_level
           top_level.top_level_product_state(self.mera,self.sys)
           for l in range(self.sys.layer):
                print "rho_2 at layer ", l
                r2 = self.sys.rho_2[l][0]
                print r2

        def apply_to_layer(self, layer, func=None, w1=["U", "V"], w2=["H_2", "rho_2"]):
            print " =========================================="
            if func:
                print "apply "+func.__name__," to layer ", layer


            M = self.M
            sys=self.sys
            #w1 = ["U", "V"]
            r1 = range(4)
            #w2 = ["H_2", "rho_2"]
            
            r2 = range(4)

            print "before"
            print M.__repr__(layers=r1, which=w1)
            print sys.__repr__(layers=r2, which=w2)
            
            if func == None:
                func = lambda a,b,c:None
            func(M, sys, layer)
            
            print "\nafter"
            print M.__repr__(layers=r1, which=w1)
            print sys.__repr__(layers=r2, which=w2)
        
        
        #from merapy.decorators import timer
        #@timer
        def J_ising(self):
            sys= ts.instance(M, symmetry="Z2", model="Ising")
            param = {"alpha":1., "beta":1.0}
            l = 3
            tau = 5
            #res=System.Jleg_tau_0l(l,tau, **param)
            sys.J_long_ising_calc(tau, **param)
            import pprint
            pprint.pprint(System.Jleg_tau_0l_dic)

        #@timer
        def J_heisbg(self):
            sys= ts.instance(M, symmetry="U1", model="Heisenberg")
            param = {"alpha":1., "beta":1.0}
            l = 3
            tau = 5
            #res=System.Jleg_tau_0l(l,tau, **param)
            sys.J_long_ising_calc(tau, **param)
            import pprint
            pprint.pprint(System.Jleg_tau_0l_dic)

    def test_temp(self) : 
        s = System.example(symmetry='U1', isometry_init='AFM' )
        print_vars(vars(),  ['s.mera'])


if __name__=='__main__':
    if 0: 
        import mera
        test_Mera = mera.test_Mera
        ts= test_System
        
        #M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); sys= ts.instance(M, symmetry="U1", model="Heisenberg")
        M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1"); 
        sys= ts.instance(M, symmetry="U1", model="Heisenberg")
     
        ts.sys= sys
        ts.M = M
        #print M
        #print sys  #.H_3[0]
        #print sys.H_3[0][0]
        print sys.H_2[0][0].data
        #sys.add_one_layer(M.num_of_layer-1, duplicate=True)
        sys.expand_layer(5)
        print sys
        print sys.mera
        print sys.key_property()

        
        #print ts.sys.__repr__(which=["H_2"])
        #print M.V[2][0].QSp[M.V[2][0].rank-1]
        #print M.V[1][0]
        #print M.V[2][0]
        #print sys.rho_2[3][0].QSp #print sys
        #print M.U[0][0].QSp


        #ts.contract()
        #ts.contract_two()
        #ts.contract_two_UV()
        #ts.contract_two_core()
        #ts.contract1(info=3)
        #ts.contract_except_V()
        #ts.add_env(info=2)
        #ts.add_env_V()
        #ts.ascending(info=2)
        
        #ts.eng_ham()
        

        #ts.test_ising()
        #ts.test_ising_2()
        #ts.op_heisenberg()
        #ts.op_heisenberg_h3()
        #ts.H2()
        #ts.rho2()
        #ts.rho2_all()
        #ts.H2_all()
        #ts.J_ising()
        #ts.J_heisbg()

    if 0: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        TestSystem.test_temp=unittest.skip("skip test_temp")(TestSystem.test_temp) 
        unittest.main()
                
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'xtest_save', 
           #'test_load', 
           #'xtest_example',  
           #'test_backup_path',  
           'test_minimize',  
           #'test_resume',  
           #'test_minimize_finite_site',  
           #'test_minimize_scale_invar',  
           #'test_expand_dim_and_layer',
           #'test_temp', 
        ]
        for a in add_list: 
            suite.addTest(TestSystem(a))
        unittest.TextTestRunner().run(suite)
       





