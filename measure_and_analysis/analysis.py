#coding=utf8
#!/usr/bin/env python
"""
    various analysis utullities.
    while obtaining reults is only half of the work, If you dont carefully analyze the relations among them 
    and gain insight you miss the other half.
    analysis should finaly give a report

    the art of data analysys:
        1. data mongling is invenitable
"""
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import input
from builtins import filter
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import *
from past.utils import old_div
import unittest
import argparse, argcomplete
import os,  sys
import copy 
import pprint
import pickle as pickle
import numpy as np
#import pandas as pd
from collections import (OrderedDict, namedtuple)
import platform
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import inspect
import warnings
import importlib
import datetime
import itertools 
import socket 
import types
from numpy import inf

from merapy.utilities import print_vars, OrderedSet 
from merapy.hamiltonian import System
#from .tabulate import tabulate
from merapy.measure_and_analysis.measurement import mera_backup_dir_finder,  measure_all as measure_all_orig
from merapy.context_util import rpyc_conn 
import merapy.measure_and_analysis.result_db as result_db_module 
from merapy.measure_and_analysis.result_db import (ResultDB, 
        MATPLOTLIBRC, 
        ResultDB_mera, ResultDB_idmrg, ResultDB_vmps, html_border, 
        ResultDB_tdvp, 
        ResultDB_proj_qmc, ResultDB_bethe_ansatz, 
        BACKUP_STATE_DIR, RESULTDB_DIR, RESULTDB_ROOT, LOCAL_HOSTNAME)
from  merapy.measure_and_analysis.result_db import (  MARKER_CYCLE, 
        AnalysisTools, AnalyticFormular
        )


import matplotlib as mpl


mpl.rcParams.update(MATPLOTLIBRC)


IS_PY3 = sys.version_info.major>2

__all__ = ['Analysis_mera', 'Analysis_vmps', 'Analysis_idmrg', 'Analysis_proj_qmc', 
         'MARKER_CYCLE', ]

class OrderedDictLazy(OrderedDict): 
    """
        a walk around for efficiency  
        lazy loading of each ResultDB, 
    """
    def __init__(self, alpha_parpath_dict, result_db_class=ResultDB, result_db_args={}, param_list=None, Alpha=None): 
        OrderedDict.__init__(self)
        self.alpha_parpath_dict = alpha_parpath_dict
        self.result_db_class= result_db_class
        self.result_db_args= result_db_args
        self.param_list = param_list if param_list is not None else []
        self.Alpha = Alpha 
    
    def __getitem__(self, a, fault_tol=True): 
        if isinstance(a, float): 
            a = (a, )
        if a[-1] == '' : 
            a = a[: -1]
        if len(a)<len(self.param_list) and self.Alpha is not None:
            a = self.Alpha(*a)
        #if hasattr(a, '_fields'):   #if isinstance(a,  namedtuple):
        #    a = tuple(a)
            
        if a not in self: 
            try: 
                p = self.alpha_parpath_dict[a]
                db =self.result_db_class(parpath=p,  
                        create_empty_db=False, **self.result_db_args) 
                self[a] = db
            except KeyError:  
                msg = '%s not in self.alpha_parpath_dict'%(a, )  #both parpathdir and db file not exists 
                warnings.warn(msg)
                db = self.result_db_class(parpath='/tmp',  dbname='empty.db.pickle', **self.result_db_args)
            except Exception as err:   #this happen when db file exists but file broken 
                warnings.warn(str(err)  + '\n a=%s'%(a, ))
                db = self.result_db_class(parpath='/tmp', dbname='empty.db.pickle', **self.result_db_args)
            return db
        else: 
            return OrderedDict.__getitem__(self, a)
    
    def copy(self): 
        res=OrderedDict.copy(self)
        res.alpha_parpath_dict = self.alpha_parpath_dict.copy()
        res.result_db_class = self.result_db_class 
        res.result_db_args = self.result_db_args
        return res 
    
    def update(self, other): 
        self.alpha_parpath_dict.update(other.alpha_parpath_dict)
        OrderedDict.update(self, other)

class AlphaList(list): 
    def __init__(alpha_list, param_list=None): 
        pass 
    
    def identify_param(self, aa): 
        if len(self.param_list)==1: 
            return self.param_list[0], 0 
            
        temp = list(zip(*aa))
        if len(aa) == 1: #then cant identify, just return empty str
            return 'NaN', 0
            
        # 这么做需要保证aa只有一个轴是不同的
        i = None   # some times aa is empty 
        for i in range(len(temp)): 
            if len(set(temp[i]))>1: 
                break 
        if i is not None : 
            param_name = self.param_list[i]
            param_name_id = i 
        else: 
            param_name = None 
            param_name_id = None 
        return param_name, param_name_id 

class Analysis(AnalysisTools, AnalyticFormular):
    """
    """
    #REMOTE_HOST = 'zhihuali@211.86.151.102:'
    ALL_FIELD_NAMES = ['energy', 'entanglement_entropy', 'magnetization', 
            'correlation', 'concurrence', 'concurrence_2', 'concurrence_neighbour']
    
    #DISTANCE_LIST = [3**i for i in range(1, 9)]
    def __init__(self,  local_root=None, param_list=None,
            param_name_list=None, 
            param_value_list=None, 
            default_resolution=None, 
            Alpha=None, 
            result_db_class=None, result_db_args=None, algorithm='all'):
        """
        
        
        """
        self.local_root = local_root
        if local_root is not None and not os.path.exists(self.local_root): 
            msg="local_root '%s' not exists, create one? yes(y) "%(self.local_root, )
            i=input(msg)
            if i.lower()=='y':
                os.makedirs(self.local_root)
                msg = 'path created'
            else:
                raise  
        self.Alpha = Alpha 
        self.param_list = Alpha._fields if Alpha is not None else param_list 
        
        self.default_resolution = default_resolution 
        self.result_db_class = result_db_class if result_db_class is not None else ResultDB
        self.result_db_args = {'analysis_class':self}
        if result_db_args is not None:
            self.result_db_args.update(result_db_args)
       
        self.algorithm = algorithm
        self.name = ''  # a name to describe self 
        
        self.alpha_parpath_dict = OrderedDict()   #this attr should always be updated, while not replaced, or else, it affects the next line
        self.alpha_rdb_dict = OrderedDictLazy( self.alpha_parpath_dict, result_db_class=self.result_db_class, 
                result_db_args=self.result_db_args, param_list=self.param_list, Alpha=self.Alpha)
                
        self.rdb_dict = self.alpha_rdb_dict
        
        self.sh_list_max = None 
        self.sub_analysis= []
        if local_root is not None: 
            self._last_modify_time = os.path.getmtime(self.local_root)
            self.set_parpath_dict()   #scan disc
        if 1: #plotting configs
            self.invert_xaxis= 0
            self.invert_alpha_order = 0
    
    def __getitem__(self, k): 
        #if isinstance(k, floa0w): 
        #    k = (k, )
        return self.alpha_rdb_dict[k]
    
    def __setitem__(self, k, v):
        self.alpha_rdb_dict[k] = v
    
    def __add__(self, other): 
        """
            a rough implementation    
            some alpha may overwritten, warn this in future  
        """
        res= self.copy()
        #res.alpha_list.extend(other.alpha_list)
        aa = self.alpha_list  +  other.alpha_list
        res.alpha_list = list(sorted(set(aa)))
        res.alpha_parpath_dict.update(other.alpha_parpath_dict)
        #re.alpha_rdb_dict.alpha_parpath_dict.update(other.alpha_parpath_dict)
        res.alpha_rdb_dict.update(other.alpha_rdb_dict)
        return res 
    
    def __repr__(self): 
        res= 'instance of %s\n'%(self.__class__.__name__, )
        res += 'local_root = %s\n'%(self.local_root, )  
        return res 

    @property
    def last_modify_time(self):
        return os.path.getmtime(self.local_root)
    
    def copy(self): 
        #self.reload_class(info=0)   #reload 一下，避免下面调用__class__时容易报错
        #res= self.__class__()
        #super(Analysis_idmrg, self).__init__(**kwargs)
        #calling above usually contradicts relaod of ipython,  so copy instead 
        res= copy.copy(self)
        res.param_list = self.param_list 
        res.alpha_list = list(self.alpha_list)
        res.alpha_parpath_dict = self.alpha_parpath_dict.copy()
        res.alpha_rdb_dict = self.alpha_rdb_dict.copy()
        return res
    
    def get_sub_analysis(self, name, fuzzy=True, recursive=False):
        res = None 
        if not recursive: 
            for i in self.sub_analysis: 
                #if name in i: 
                if (name in i) & fuzzy or (name == i)&(not fuzzy): 
                    res= getattr(self, i)
                    break 
        if res is None: 
            print('%s not found'%(name, ))
        return res 
    
    def get_param_range(self, param, aa=None, param_list=None): 
        aa = self.alpha_list if aa is None else aa 
        param_list = param_list if param_list is not None else self.param_list 
        aa = [a if hasattr(a, '__iter__') else (a, )  for a in aa]
        k = param_list.index(param)
        x = [a[k] for a in aa]
        x = list(set(x))
        x.sort()
        #x=np.array(x)
        return x

    def get_run_failure_bac(self, sh, aa_all=None): 
        aa_all=self.alpha_list if aa_all is None else aa_all 
        #aa_fail=OrderedDict()
        #DD=[16, 20, 40, 80]
        #for D in DD:
        #    temp= self.filter_alpha(sh=(0,D))
        #    aa=list(set(aa_all)- set(temp))
        #    aa.sort()
        #    aa_fail[D] =  aa
        temp= self.filter_alpha(aa=aa_all, sh=sh)
        aa_fail=list(set(aa_all)- set(temp))
        aa_fail.sort()
        return aa_fail 
    
    def get_run_failure(self, sh, aa=None, resolution=None): 
        resolution = resolution if resolution is not None  else self.default_resolution 
        aa = aa if aa is not None else self.alpha_list 
        for p in self.param_list: 
            temp = self.get_param_range(p, aa) 
            print_vars(vars(),  ['p, temp'])
    
    def filter_alpha(self, aa=None, sh=None, field=None, no_field=None, param_list=None, 
            surfix='', resolution='default', auto_refrese=1,  no_bad=False, field_range=None,   **kwargs): 
        """
            from wigner_crystal 
            parmas: 
                no_bad: remove points if db['status'] == 'bad' 
                field_range: is used specifically for select (mainly for bad) points within the range
        """
        if auto_refrese:
            if self.last_modify_time >  self._last_modify_time:
                self._last_modify_time = self.last_modify_time
                print('set_parpath_dict called in filter_alpha')
                self.set_parpath_dict()
        
        aa = aa if aa is not None else list(self.alpha_list) 
        param_list = param_list if param_list is not None else self.param_list #['a', 'h', 'V']
        if len(aa)<1:
            return []
        num_param = len(param_list)
        #make each a into a tuple  if it is not 
        aa = [a if hasattr(a, '__iter__') else (a, ) for a in aa]
        #if no_surfix: 
        #    #exclude such as (1.0, 'sufix') from aa 
        #    aa = filter(lambda x: not isinstance(x[-1], str),  aa)
        
        if surfix is not None : 
            #surfix = None means all surfix; surfix = '' means empty surfix 
            if surfix == '' :
                aa = [x for x in aa if not isinstance(x[-1], str)]
            elif '<' in surfix or '>' in surfix:
                func = lambda x: eval(surfix.replace('surfix', x[-1]))
                aa = list(filter(func, aa))
                
            else:
                aa = [x for x in aa if x[-1]==surfix]
            
            if len(aa)==0: 
                msg='in filter_alpha arg surfix "{}" may not correct,  return aa=[]'.format(surfix)
                warnings.warn(msg)
                return []

        
        for k, v in list(kwargs.items()): 
            if v is None:  # None means no constraint 
                continue 
            if k not in param_list:
                warnings.warn('param name {} not in param_list'.format(k))
                continue
            i = param_list.index(k)
            
            if isinstance(v, float): 
                #if hasattr(aa[0], '__iter__'): 
                #    func = lambda x: x[i]== v
                #else:
                #    func = lambda x: (x, )[i]== v
                v = '%s==%f'%(k, v)
            #elif isinstance(v, str): 
            if 1: 
                if hasattr(aa[0], '__iter__'): 
                    func = lambda x: eval(v.replace(k, 'x[%d]'%i))
                else: 
                    func = lambda x: eval(v.replace(k, '(x, )[%d]'%i))
            aa = list(filter(func, aa))
            if len(aa)<1:
                return []
        if resolution == 'default' : 
            resolution = getattr(self, 'default_resolution', None)
        
        if resolution is not None:
            if not hasattr(resolution, '__iter__'): 
                resolution = float(resolution)
                #use a[: num_param] is because there may be surfix a = (1.0, 'a') then 'a' is excluded
                func=lambda a:   all([round(old_div(x,resolution), 8).is_integer() for x in a[: num_param]])
            else: 
                ran = range(num_param)
                #round is needed, because for example 0.3/0.1 = 02999999999996 
                func=lambda a:   all([round((old_div(a[i],resolution[i])), 8).is_integer() for i in ran])
            
            aa = list(filter(func, aa))
            
        #aa = map(lambda a:  a if len(a)>1 else a[0], aa)
        
        if sh is not None : 
            if 0:  #改成不需要检查file，直接从db['energy']判断
                def func(a): 
                    dir = self.alpha_parpath_dict.get(a, 'nowhere')
                    fn = self.result_db_class.shape_to_backup_fn(sh)
                    path  = '/'.join([dir, fn])
                    print('ppp', path) 
                    return os.path.exists(path)
            else: 
                if sh[1] != 'max': 
                    func = lambda a:  self[a].has_key_list(['energy', sh])
                else:
                    func = lambda a:  self[a].has_key_list(['dim_max', sh[0]])
            aa = list(filter(func, aa))

        if field is not None: 
            if sh is None: 
                func = lambda a: field in self[a]
            else: 
                func = lambda a: self[a].fetch_easy(field, sh)
            aa = list(filter(func, aa)) 
            
        
        if no_bad: 
            aa = [a for a in aa if self[a]['status'] == 'good']
        
        if no_field is not None: 
            assert field is None, 'field and no_field cant be given together'
            temp = self.filter_alpha(aa, field=no_field, resolution=resolution,  
                    sh=sh, surfix=surfix,  param_list=param_list, **kwargs)
            aa = list(set(aa)-set(temp))
            aa.sort()
        
        if field_range is not None: 
            assert sh is not None , "'sh' is requird when select by field_range"
            aa = [a for a in aa if  
                    field_range[0] <= self[a].fetch_easier(field, sh) <= field_range[1]]
        if auto_refrese:
            temp = [a for a in aa if self[a].last_modify_time>self[a]._last_modify_time]
            if temp:
                print('set_rdb_dict called in filter_alpha')
                self.set_rdb_dict(temp)
            
        return aa
    
    alpha_filter = filter_alpha  
    
    def alpha_list_to_grid(self, aa=None, k1=0, k2=1):
        """
            used in e.g. 3D plot 
        """
        x1 = [a[k1] for a in aa]
        x2 = [a[k2] for a in aa]
        x1 = list(set(x1)); x1.sort(); x1=np.array(x1)
        x2 = list(set(x2)); x2.sort(); x2=np.array(x2)
        return x1, x2
    
    def scan_backup_dir(self, root='./'): 
        """     
            automatically scan and identify folders that store mera tensors
        """
        res = mera_backup_dir_finder(root)
        return res
    
    def scan_alpha(self, root=None, surfix='', signiture=None, info=0): 
        """
            a valid backup_parpath name is like
                    a=1.0-b=2.0 
                    a=1.0-b=m2.0
                or with a surfix in the end 
                    a=1.0-b=2.0-1
                    a=1.0-b=2.0-a
        """
        if root is None: 
            if self.local_root is not None : 
                root = self.local_root
            else: 
                root = './'
        dir_list = os.listdir(root)
        if info>1: 
            print('dir_list is', dir_list)
        signiture = signiture if signiture is not None else '='
        dir_list = [d for d in dir_list if signiture in d]
        alpha_parpath_dict = OrderedDict()
        
        if self.Alpha is not None:  
            def parse(term): 
                if '=' in term: 
                    key, val = term.split('=')
                    if 'm' == val[0]:  #change 'm' to minus sign
                        val = '-'  + val[1: ]
                    try: 
                        val = float(val)  #val=(float(val), )   #todo:  make alpha always a tuple 
                    except: 
                        val = val
                    return (key, val)
                else: 
                    return ('surfix', term)# term is surfix 
            #print_vars(vars(),  ['dir_list'])
            for name in dir_list: 
                
                nn = name.split('-')
                temp = [parse(a) for a in nn]
                if temp[-1][0] != 'surfix' :  #no surfix
                    alpha = self.Alpha(** dict(temp))
                    alpha = tuple(alpha)
                else:
                    surfix = temp[-1][1]
                    alpha = self.Alpha(** dict(temp[:-1]))
                    alpha += (surfix, ) 
                
                alpha_parpath_dict[alpha] =  '/'.join([root, name])
        
        else:
            def parse(xx): 
                if '=' in xx: 
                    xx = xx.split('=')[1]
                    #print('xxxxxxxxxxxxxxxxx', xx)
                    if 'm' == xx[0]:  #change 'm' to minus sign
                        xx = '-'  + xx[1: ]
                else: 
                    return xx   # xx is surfix 
                try: 
                    res=float(xx)
                    #res=(float(xx), )   #todo:  make alpha always a tuple 
                except: 
                    res= xx
                return res
                
            for name in dir_list: 
                aa = name.split('-')
                alpha = tuple([parse(a) for a in aa])   #no mater one or many param,  use tuple as key uniformly 
                if surfix != '': 
                    alpha = str(alpha) + '-' +  surfix
                alpha_parpath_dict[alpha] =  '/'.join([root, name])
        
        return alpha_parpath_dict
    
    def scan_sub_analysis_top(self): 
        temp = os.listdir(self.local_root)
        temp = [x for x in temp if '=' not in x]
        temp = [x for x in temp if os.path.isdir('/'.join([self.local_root, x]))]
        #print_vars(vars(),  ['temp']) 
        temp = [{'sub_dir': i} for i in temp]
        return temp 
    
    def create_db_parpath(self, param_val_tuple):
        param_name = self.param_list
        assert param_name is not None 
        assert len(param_val_tuple)==len(param_name)
        param_val_tuple = list(map(str, param_val_tuple))
        fn = list(zip(param_name, param_val_tuple))
        fn = ['='.join(x) for x in fn]
        fn = '-'.join(fn)
        db_parpath = '/'.join([self.local_root, fn])
        if not os.path.exists(db_parpath): 
            #cmd = 'mkdir %s'%db_parpath
            #os.system(cmd)
            os.makedirs(db_parpath)
            msg = 'create parpath ... done'
        else:
            msg = 'parpath exists, denied'
        print(msg)
        return db_parpath
            
    def get_db(self, param):
        p = self.create_db_parpath(param)
        db =self.result_db_class(parpath=p,  
                create_empty_db=False, **self.result_db_args) 
        
        return db
    
    def add_sub_analysis(self,  hook_list='auto'): 
        """
            pass
            example:
                an_main.add_sub_analysis(hook_list=[{'sub_dir':dir_name}])
        
        """
        if hook_list is None: 
            hook_list = []
        elif hook_list == 'auto': 
            hook_list = self.scan_sub_analysis_top()

        for h in hook_list: 
            #sub_dir = h if isinstance(h, str) else h['sub_dir']   #h may be a dir name or in the old format {'sub_dir': dir_name_str}
            sub_dir = h['sub_dir']
            Alpha = h.get('Alpha', self.Alpha)
            
            sub_root = '/'.join([self.local_root, sub_dir])
            if 'name' in h: 
                name = 'an_' +  h['name']
            else: 
                name = 'an_%s'%(sub_dir.replace('/', '_'), )
            if name in self.sub_analysis: 
                continue 
            antemp = self.__class__(local_root=sub_root, 
                    result_db_class= self.result_db_class, 
                    default_resolution=self.default_resolution, 
                    Alpha=Alpha, 
                    param_list=self.param_list)
            antemp.name = name
            setattr(self, name, antemp)
            self.sub_analysis.append(name)
        
        # below is recursive. may be not safe or very slow!!    
        if hasattr(self, 'an_test'): 
            warnings.warn('auto added sub analisys of an_test')
            self.an_test.add_sub_analysis('auto')
        if hasattr(self, 'an_temp'): 
            warnings.warn('auto added sub analisys of an_temp')
            self.an_temp.add_sub_analysis('auto')
    
    def search_sub_analysis(self, s, recursive=False): 
        for sub in self.sub_analysis: 
            if s in sub: 
                return getattr(self, sub)
        print('%s not found in self.sub_analysis %s'%(s, self.sub_analysis))
    
    def set_rdb_dict(self, alpha_list=None, create_empty_db=1, 
            db_version=None, algorithm=None, reload_db_class=0, info=0):        
        if reload_db_class: 
            self.reload_db_class()
        
        part_update = 0
        if alpha_list is not None: 
            part_update = 1
            
        alpha_list = alpha_list if alpha_list is not None else  sorted(self.alpha_parpath_dict)
        alpha_list = [a if isinstance(a, tuple) else (a, )  for a in alpha_list]
        error = []
        for a in alpha_list: 
            try: 
                p = self.alpha_parpath_dict[a]
                db = self.result_db_class(parpath=p, 
                        create_empty_db=create_empty_db, version=db_version, 
                        algorithm=algorithm)
                #rdb[a] = db
                self.alpha_rdb_dict[a] = db
            except Exception as err: 
                error.append((a, err))
        msg = 'alpha_rdb_dict is reset'
        
        if part_update: 
            msg += '  for %s'%(alpha_list, )
        else: 
            msg += '  for all alpha' 
        if info>0: 
            print(msg)
        
        if len(error)>1: 
            msg = 'error in set_rdb_dict\n'
            msg +=  str(error)
            warnings.warn(msg)
    
    def set_parpath_dict(self, root=None, signiture=None, info=0): 
        if root is None: 
            root = self.local_root
        assert root is not None , 'param root reqired'
        dic = self.scan_alpha(root=root, signiture=signiture)
        self.alpha_parpath_dict.update(dic)
        
        if info>1: 
            msg = 'parpath found are %s'%(dic, )
            print(msg)
        if info>0: 
            msg = 'alpha_parpath_dict has been reset'
            print(msg)
        
        alpha_list = list(dic.keys())
        try:
            alpha_list.sort()  #some times it can not be sorted,  due to not different format of alpha in the list
        except TypeError as err:
            warnings.warn(str(err))
        self.alpha_list = list(alpha_list)
    
    def reset(self, aa=None, root=None):
        """
            just a convenient function to put the two in one func
        """
        self.set_parpath_dict(root=root)
        self.set_rdb_dict(aa)
    
    def search_alpha(self, alpha_str, alpha_list=None):
        if not isinstance(alpha_str, str): 
            alpha_str = str(alpha_str)
        if alpha_list is None: 
            all_key = list(self.alpha_parpath_dict.keys())
        else: 
            all_key = alpha_list
        
        res= []
        for a in all_key: 
            if alpha_str in str(a): 
                res.append(a)
        res.sort()
        return res
    
    def shape_to_backup_fn(self, sh): 
        return System.shape_to_backup_fn(*sh)
    
    def search_field(self, field_name): 
        res = []
        for a in list(self.alpha_parpath_dict.keys()): 
            if field_name in self[a]: 
                res.append(a)
        return res
            
    def upgrade_result_db(self): 
        res = mera_backup_dir_finder(root)
    
    def measure_all(self, aa=None, sh_list=None, sh_min=None, sh_max=None,  which=None, 
            exclude_which=None, force=0, fault_tolerant=1, use_dist_comp=0, recursive=False,  **kwargs): 
        sh_list = sh_list if sh_list is not None else 'all'
        if aa is None: 
            aa = self.alpha_list
        else: 
            aa = [a if not isinstance(a, float) else (a, ) for a in aa]
        aa = list(aa)
        dir_list = [self.alpha_parpath_dict[a] for a in aa]
        if 'C:' in dir_list[0] and use_dist_comp:
            func = ResultDB.parpath_map.__func__
            dir_list = [func(a).replace('dropbox', '') for a in dir_list]
        dir_list = [d.replace(RESULTDB_DIR, BACKUP_STATE_DIR) for d in dir_list]
        if use_dist_comp: 
            kwargs['use_local_storage'] = 1 
            
        measure_all_orig(dir_list, mera_shape_list=sh_list, sh_min=sh_min, sh_max=sh_max, which=which, 
                exclude_which = exclude_which, use_dist_comp=use_dist_comp,  
                fault_tolerant = fault_tolerant, force=force, **kwargs)
        self.set_rdb_dict(aa)

    
    def measure_all_new(self, aa, sh_list=None, sh_min=None, sh_max=None,  
            which=None, num_of_threads=2, force=0, fault_tolerant=1,
            **kwargs): 
        """
            take use of db.measure and use distributive compute 
        """
        if 1:
            aa = [a if not isinstance(a, float) else (a, ) for a in aa]
            aa = list(aa)
            if isinstance(sh_list, tuple):
                sh_list = [sh_list]
            if isinstance(which, str):
                which = [which]
            
        temp = []
        msg = ''
        for a in aa:
            msg += '\n' 
            db = self[a]
            msg += '\n{}\n\t'.format(a)
            
            for sh in sh_list:
                add_it = True
                msg += '  {} '.format(sh)
                if not force:
                    if sh[1] == 'max' or sh[1] == 0:
                        _sh = sh[0], db.get_dim_max_for_N(sh[0])
                    else:
                        _sh = sh
                    has_sh = True if db.has_shape(_sh, from_energy_rec=1) else False
                    if not has_sh:
                        msg += 'not exist'  
                        
                        add_it = 0
                    else:    
                        ww = list(which)
                        for w in ww:
                            if db.has_key_list([w, _sh]):
                                ww.remove(w)
                        if not ww:
                            msg += 'found' 
                            add_it = 0
                else:
                    ww = which
                if add_it:
                    msg += 'added' 
                    
                    temp.append((a, sh, ww))    
        #print msg
        for t in temp:
            print(t)
        i=input('submit %s jobs? yes(y)\n'%(len(temp)))
        if i.lower()!='y':
            print('canceled')
            return 
        for a, sh, w in temp:
            db = self[a]
            db.measure(sh, which=which, force=force, num_of_threads=num_of_threads, **kwargs)

    def measure_two(self, alpha_paier_list): 
        """
            status: 
                to implement 
            measure e.g. overlap
        """
        pass 
    
    def load(self, path): 
        res = None
        with open(path, 'rb') as inn: 
            res=pickle.load(inn)
        return res
    
    def _outline(self, which, alpha_list=None, mera_shape_list=None, 
            show_eng_diff=False,  grid_view=False, info=0): 
        if info>-1: 
            dir_name = '' if self.local_root is None else '...'+self.local_root[-30:]
            print('outline %s at %s'%(dir_name, str(datetime.datetime.now()), ))
        data=[]
        if mera_shape_list is None: 
            mera_shape_list=[(4, 4)] + [(8, i) for i in [4,5]] + [(12, i) for i in [4,5]]
        if alpha_list is None: 
            alpha_list = self.alpha_list
        
        res = []
        for a in alpha_list:
            if which == 'file':  
                sh_list = self.alpha_rdb_dict[a].get_shape_list(from_energy_rec=0)
                
                diff_list = []
            
            elif which == 'db': 
                db=self.alpha_rdb_dict[a]
                if db.get('version')==1.0:
                   
                    keys= db.get('energy')
                    if keys is None: 
                        pass
                    keys= list(keys.keys())
                    keys.sort()
                    sh_list = keys
                    
                else: 
                    keys= list(db.keys())
                    keys= [k[: 2] for k in keys]
                    keys= list(set(keys))
                    keys.sort()
                    sh_list = keys
                    
            if not grid_view: 
                print(a)
                print('\t', sh_list, '\n')
                
            d= [str(a)]
            for sh in mera_shape_list:
                if sh in sh_list:
                    d.append('   *   ')
                else:
                    d.append('       ')
            data.append(d)
        
        if grid_view:     
            print(tabulate(data, headers= ['a']+list(map(str,mera_shape_list)), tablefmt='grid', stralign='center')) 
    
    def outline_db(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0):
        self._outline('db', alpha_list=alpha_list, mera_shape_list=mera_shape_list, info=info)
    
    def outline_file(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0): 
        self._outline('file', alpha_list=alpha_list, mera_shape_list=mera_shape_list, show_eng_diff=show_eng_diff, info=info)
    
    def outline(self, alpha_list=None, indent=0, recursive=1): 
        """
            outline_db, outline_file will be deprecated in future 
        """
        #self.set_parpath_dict()
        alpha_list = alpha_list if alpha_list is not None else self.alpha_list
        tab = '\t'*indent 
        sh_list_max = set()
        for a in alpha_list: 
            db = self[a]
            sh_list = db.get_shape_list()
            sh_list_max = sh_list_max.union(sh_list)
        
        sh_list_max = list(sh_list_max)
        sh_list_max.sort()
        self.sh_list_max = sh_list_max
        msg = ''
        #msg += ''.join([tab, str(alpha_list), '\n']) 
        for p in self.param_list: 
            msg += ''.join([tab, p, ': ']) 
            temp = self.get_param_range(p, alpha_list) 
            msg += ''.join(['\n',tab, '\t', str(temp), '\n']) 
        msg += ''.join([tab, 'sh_list: \n']) 
        msg += ''.join([tab, '\t', str(sh_list_max), '\n'])
        
        if 1: #extract dir surfix 
            surfix_list = self.get_surfix_list(alpha_list, force=1)
            msg += ''.join([tab, 'surfix_list: \n']) 
            msg += ''.join([tab, '\t', str(surfix_list), '\n'])
        
        if 1: 
            msg += ''.join([tab, 'sub an: \n']) 
            msg += ''.join([tab, '\t', str(self.sub_analysis), '\n'])
            
        print(msg) 
        if recursive: 
            for a in self.sub_analysis: 
                an = getattr(self, a)
                print(''.join([tab, a]))
                an.outline(recursive=1, indent=indent+1)
    
    def preprocess_alpha_list(self, alpha_list, alpha_remap=None): 
        parser = lambda a: a if not isinstance(a, str) else  float(a.split('-')[0]) 
        if 1: 
            aa = OrderedDict()
            aa_bad = []
            for a in sorted(alpha_list):
                if a in self.alpha_parpath_dict:
                    aa[parser(a)] = self.alpha_rdb_dict[a]
                else: 
                    aa_bad.append(a)
        
        if len(aa_bad)>0: 
            print('keys %s are not in self.alpha_parpath_dict'%aa_bad)
            #warnings.warn('keys %s are not in self.alpha_parpath_dict'%aa_bad)
        return aa
            
    def fetch_many(self, field_name, aa, sh, key_list=None): 
        aa_dic = self.preprocess_alpha_list(aa)
        res= OrderedDict()
        for a in aa_dic:
            db = aa_dic[a]
            rec=db.fetch_easy(field_name, sh, key_list)
            res[a] = rec
        return res
             
    def merge_remote_db(self, aa, remote_root=None): 
        raise  #deprecated 
        remote_root = remote_root if remote_root is not None else self.remote_root
        for a in aa:
            db=self.alpha_rdb_dict[a]
            dir = self.get_backup_dir_name(a)
            dir = remote_root + '/' + dir
            db.merge_remote(dir)        
    
    def correaltion_exact(self, r=None, which='shastry'): 
        """
            the exact value is: 
                [(3, -0.044424436468189332), (9, -0.014200837573249929), (27, -0.0046643666206658257), (81, -0.0015470704943921049), (243, -0.00051483226149443202)]
        """
        #if isinstance(r, iterab)
        r = r if r is not None else np.array([3**i for i in range(1, 6)])
        pi = np.pi
        if which == 'shastry': 
            res = 1./(1.0*pi*r)*sici(r*pi)[0]*(-1)**r
        #return zip(r, res)
        return res

    def get_backup_dir_name(self, a): 
        dir = self.alpha_parpath_dict[a]
        #dir = os.path.dirname(dir)
        dir = dir.split('/')[-1]
        return dir
    
    def transfer_pickle_file(self, alpha_list, sh, overwrite=1): # remote_root, local_root): 
        raise  #deprecated 
        #remote_root = 'zhihuali@211.86.151.102:mera_backup_tensor/run-long-better/'
        #locale_root = '~/Documents/mera_backup_tensor/run-long-better/'
        
        #A_root = 'zhihuali@211.86.151.102:' + remote_root
        alpha_list = list(alpha_list)
        A_root = self.REMOTE_HOST + self.remote_root
        #B_root = '~/Documents/mera_backup_tensor/run-long-better/'
        B_root = self.local_root
        dim, layer  = sh[:2]
        if len(sh)>2:  
            nqn=sh[2]
        else: 
            nqn = None
        for a in alpha_list: 
            
            
            #alpha_parpath = p1 if p1 is not None else 'alpha=%s/'%a
            alpha_parpath = self.get_backup_dir_name(a)
            A_dir = '/'.join([A_root, alpha_parpath]) 
            B_dir = '/'.join([B_root, alpha_parpath])
            
            if layer == 4:  
                fn = '%s.pickle'%dim
            else: 
                fn = '%s-%slay.pickle'%(dim, layer-1)
            if nqn == 5: 
                fn = '5_' + fn
            
            A_path = A_dir + '/'  +  fn
            B_path = B_dir + '/'  +  fn
            if not os.path.exists(B_dir): 
                os.mkdir(B_dir)
                print('local dir not found, mkdir %s'%B_dir)
            
            #if not only_measure: 
            if 1: 
                cmd = 'scp %s %s'%(A_path, B_dir)                     
                msg = 'scp %s ---> %s'%(A_path[-300:], B_path[-300: ])                     
               
                if not overwrite and os.path.exists(B_path): 
                    msg += '\t --- already exists, skip it' 
                    continue 
                    
                status=os.system('scp %s %s'%(A_path, B_dir))
                if status  ==  0: 
                    msg+= '   --- done'
                else: 
                    msg+= '     --- FAILED status= %d'%(status, )
                    #print 'A_path, %s\nB_path, %s'%(A_path, B_path)
                print(msg)
            
           
            #if measure: 
            if 0: 
                path = B_path
                mm = Measurement(path=path)
                S= System.load(path)
                mm.measure_all(S, path, which=None, show=False, exclude_which=['scaling_dim'])
    
    def get_shape_list_all(self, N=None, D=None, aa=None, sh_min=None, 
            sh_max=None, from_energy_rec=True, only_return_N=False, force=0):
        #if self.sh_list_max is not None and not force : 
        if self.sh_list_max  and not force : 
            sh_list_max = list(self.sh_list_max)
        else:
            sh_list_max = set()
            aa = aa if aa is not None else list(self.alpha_list)
            rdb_cls = self.result_db_class 
            #if not from_energy_rec:
            if socket.gethostname()==LOCAL_HOSTNAME and not from_energy_rec: 
                for a in aa: 
                    #using db meth is slow 
                    #db = self[a]             
                    #sh_list = db.get_shape_list()  
                    dir = self.alpha_parpath_dict[a]
                    fn = rdb_cls.get_fn_list.__func__(self, dir)
                    sh_list = [rdb_cls.parse_fn.__func__(self, f) for f in fn]
                    sh_list_max = sh_list_max.union(sh_list)
            else:
                warnings.warn('get_shape_list_all using energy rec')
                for a in aa: 
                    db = self[a]
                    sh_list = list(db.get('energy', {}).keys())
                    sh_list_max = sh_list_max.union(sh_list)
            
            sh_list_max = list(sh_list_max)
            self.sh_list_max = sh_list_max 
        if N is not None:
            sh_min = (N, 0) if sh_min is None else (N, sh_min[1])
            sh_max = (N, np.inf) if sh_max is None else (N, sh_max[1])
        if D is not None: 
            sh_min = (0, D) if sh_min is None else (sh_min[0], D)
            sh_max = (np.inf, D) if sh_max is None else (sh_max[0], D)
            
        if sh_max is not None : 
            sh_list_max=[x for x in sh_list_max if x[0]<=sh_max[0] and x[1] <= sh_max[1]]
        if sh_min is not None : 
            sh_list_max=[x for x in sh_list_max if x[0]>= sh_min[0] and x[1]>= sh_min[1]]
        if only_return_N:
            sh_list_max = [i[0] for i in sh_list_max]
            sh_list_max = list(set(sh_list_max))
        sh_list_max.sort()
        return list(sh_list_max)
    
    def get_surfix_list(self, aa=None, force=0):
        if hasattr(self, 'surfix_list') and not force: 
            return self.surfix_list 
        else: 
            aa = aa if aa is not None else self.alpha_list 
            temp = [x[-1] if isinstance(x[-1], str) else '' for x in aa]
            surfix_list = list(set(temp))
            #if '' in surfix_list: 
            #    surfix_list.remove('')
            surfix_list.sort()
            self.surfix_list = surfix_list
            return surfix_list 
        
    
    def get_fn_list_del(self):
        nlist = os.listdir(self.parpath)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" in i]
        print("all pickle files under specified path are %s"%pickle_files)
        
        def parse_fn(fn):
            if "transition" in fn: #'12-transition.pickle'
                dim = "0"
            elif "_" in fn:  #like "5_13.pickle or 5_13-4lay.pickle"
                dim = fn.split("-")[0].split("_")[1].split(".")[0]
            else:
                if "lay" in fn: #"like "8-4lay.pickle""
                    dim = fn.split("-")[0]
                else:   
                     #like "8.pickle"
                    dim = fn.split(".")[0]
            
            i = fn.find("lay")
            layer = fn[i-1] if i>0 else "3"
            #print "ddd", [dim, layer]
            return eval(dim), eval(layer)

        temp = []
        for i in pickle_files:
            print(i, parse_fn(i))
            try:
                dim, layer = parse_fn(i)
            except:
                raise Exception(i)
            temp.append((layer, dim, i))
        temp.sort(key=lambda x:x[:2], reverse=True)
        return temp
    
    def delete_db(self, aa): 
        for a in aa: 
            try: 
                #print a
                #self[a].delete_db()
                db_parpath = self.alpha_parpath_dict[a]
                path = '/'.join([db_parpath, 'RESULT.pickle.db'])
                print('removing db %s and its parpath'%(path[-30: ]), end=' ') 
                os.remove(path) 
                os.rmdir(db_parpath)
                print('  ... done')
            except Exception as err:
                print('  ..failed')
                print(err) 
    
    def delete_attr_in_file(self, name, aa,  sh): 
        #aa=self.filter_alpha(sh=(1000, 20))
        for a in aa[2:]:
            db=self[a]
            dir=db.parpath
            try:
                name= db.backup_fn_gen(sh)
                path='/'.join([dir,name])
                print(path)
                s=vmps.load(path)
                s.pop(name)
                vmps.save(s,path)
            except Exception as err:
                print(err)

    def delete_file(self, aa, info=1): 
        """
            
        """
        for a in aa: 
            db = self[a]
            db.delete_file(sh='all', delete_rec=0, need_comfirm=0, info=1)
    
    def delete_dir(self, aa, ):
        """
            delete the data points , i.e delete *both* tensor file and db file
        """
        for a in aa: 
            db = self[a]
            msg = 'a=%s'%(a, )
            try: 
                msg += '\n\trm dir %s ...'%(db.state_parpath, )
                db.delete_file('all', delete_rec=0, need_comfirm=0, info=0)
                os.rmdir(db.state_parpath)
                msg += '\n\trm dir %s ...'%(db.parpath, )
                
                db.delete_db(info=0)
                os.rmdir(db.parpath)
                msg += '\n\tdone' 
            except Exception as err: 
                msg += str(err)
            print(msg) 
            
    
    def rename_folder_surfix(self, aa, new_sur, dry_run=1): 
        #aa = self.filter_alpha(aa=aa, surfix=old_sur)
        for a in aa: 
            db = self[a]
            db.rename_folder_surfix(new_sur, dry_run=dry_run)
            print('\n')
        self.set_parpath_dict(info=1)
    
    def move_folder(self, aa, local_root): 
        for a in aa: 
            db = self[a]
            db.move_folder(local_root)
        self.set_parpath_dict(info=1)  # note the local_root parpath_dict need also reset 
    
    def get_fn_list(self):
        nlist = os.listdir(self.parpath)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" == i.split('.')[-1]]
        print("all pickle files under specified path are %s"%pickle_files)
        return pickle_files
    
    def data_mangle(self, force=False):
        if not 'record.pickle.db' in os.listdir(self.parpath) or force:
            res = OrderedDict()
            if not hasattr(self, 'fn_list'):
                self.fn_list = self.get_fn_list()
            #print self.fn_list
            #exit()
            for i in self.fn_list:
                res[i] = {}
                path = self.parpath + '/'  + i
                print(i, path)
                if 0:
                    inn = open(path, 'rb')
                    S= pickle.load(inn)
                    inn.close()
                S= System.load(path)  #System.load is not exactly pickle.load
                res[i]['energy'] = S.energy
                res[i]['energy_record'] = S.energy_record
            #print res
            out = open(self.parpath +'/' +  'record.pickle.db', 'wb')
            print('file record.pickle.db is generated')
            pickle.dump(res, out)
            out.close()
    
    def dim_vs_energy_err(self):
        pass

    def find_max_dim_file_bac(self, path):
        nlist = os.listdir(path)
        
        pickle_files = [i for i in nlist if "pickle" in i]
        print("all pickle files under specified path are %s"%pickle_files)
        
        def parse_fn(fn):
            if "transition" in fn: #'12-transition.pickle'
                dim = "0"
            elif "_" in fn:  #like "5_13.pickle or 5_13-4lay.pickle"
                dim = fn.split("-")[0].split("_")[1].split(".")[0]
            else:
                if "lay" in fn: #"like "8-4lay.pickle""
                    dim = fn.split("-")[0]
                else:   
                     #like "8.pickle"
                    dim = fn.split(".")[0]
            
            i = fn.find("lay")
            layer = fn[i-1] if i>0 else "3"
            #print "ddd", [dim, layer]
            return eval(dim), eval(layer)

        temp = []
        for i in pickle_files:
            #print i, parse_fn(i)
            try:
                dim, layer = parse_fn(i)
            except:
                raise Exception(i)
            temp.append((layer, dim, i))
        temp.sort(key=lambda x:x[:2], reverse=True)

        res= temp[0][2]
        return res
    
    @staticmethod    
    def find_max_dim_file(path):
        temp = Analysis.get_fn_list(path)
        res= temp[0][2]
        print('head is %s'%res)
        return res
    
    def monitor(self): 
        """
            moniting the excution state of jobs
            what to watch: 
                eng, eng_diff, time spend, etc
        """
        pass
        
    def _plot(self, x, y, **kwargs):
        if IS_PY3:
            fig=ResultDB._plot(None, x, y, **kwargs)  
        else:
            fig=ResultDB._plot.__func__(None, x, y, **kwargs)  
        return fig
    
    def _plot3d(self, xx, yy, data, which_plot=None, zfunc=None,   **kwargs):
        """
            pass
        """
        which_plot = 'imshow' if which_plot is None else which_plot 
        figsize = kwargs.get('figsize')
        figsize = (4, 4) if figsize is None else figsize
        fig = kwargs.get('fig', None)
        ax = kwargs.get('ax', None)
        info = kwargs.get('info', 0)
        projection = '3d' if which_plot in ['surface'] else  None 
        if fig is None and ax is None: 
            fig=plt.figure(figsize=figsize)
        if ax is None: 
            if len(fig.axes)==0: 
                ax=fig.add_subplot(111, projection=projection)
            else: 
                ax_found = False
                #ax = fig.axes[-1]   # improved by below
                for aa in fig.axes: 
                    if len(aa.lines)==0: 
                        ax = aa
                        ax_found = True
                        break 
                if not ax_found: 
                    raise Exception('ax is not defined; examine the layout of fig')
        fig = ax.figure  
        XX,YY=np.meshgrid(xx,yy)
        
        style_dic = {
                'cmap': None, 
                'interpolation': 'nearest', 
                'alpha': 0.3 if which_plot in ['surface'] else 0.9 ,  
                'extent': None, 
                'aspect':'auto',
                }
        for t in style_dic: 
            if t in kwargs: 
                style_dic[t]=kwargs[t]
     
        not_found = []
        data = data.T
       
        if kwargs.get('make_up'):  
            make_up = kwargs['make_up']
            make_up_list = []
            m, n = data.shape 
            #print np.where(np.isnan(data))
            aaa, bbb = np.where(np.isnan(data))
            ab = list(zip(aaa.tolist(), bbb.tolist()))
            if make_up == 'v' : 
                for a, b in ab: 
                    if a>0 and a <m-1: 
                        data[a, b] = (data[a-1, b]  +  data[a+1, b])/2. 
                        make_up_list.append((a, b))
                    elif a == 0: 
                        data[a, b] = data[a+1, b]
                    elif a == m-1: 
                        data[a, b] = data[a-1, b]
                    else: 
                        raise 
            elif make_up == 'h' : 
                for a, b in ab: 
                    if b>0 and b <n-1:
                        #print 'bbbb', b, m, n  
                        data[a, b] = (data[1, b-1]  +  data[a, b + 1])/2. 
                        make_up_list.append((a, b))
                    elif b == 0: 
                        data[a, b] = data[a, b+1]
                    elif b == n-1: 
                        data[a, b] = data[a, b-1]
                    else: 
                        raise 
            else: 
                raise ValueError('make_up={make_up} should be in ["h", "v"]'.format(**vars()))

        X,Y,Z=XX,YY,data
        if zfunc is not None : 
            Z = zfunc(Z)
        if kwargs.get('zmin'): 
            zmin = kwargs['zmin']
            Z[Z<zmin] = np.nan 
        if kwargs.get('zmax'): 
            zmax= kwargs['zmax']
            Z[Z>zmax] = np.nan 
        cb = None 
        if which_plot == 'contour': 
            pass
            cb = ax.contour(X, Y, Z, 
                   colors=kwargs.get('colors', None), 
                   levels=kwargs.get('levels', None), 
                   linestyles=kwargs.get('linestyles', None), 
                    )
            ax.clabel(cb, inline=kwargs.get('inline', True), fontsize=10)
            
        elif which_plot == 'contourf': 
            #levels control intersection at where 
            #cb = ax.contourf(Y, X, Z, zdir='z', offset=0, cmap=plt.cm.coolwarm, alpha=0.9, aspect=1 )
            cb = ax.contourf(X, Y, Z, zdir='z', offset=0, 
                    levels=kwargs.get('levels', None), 
                    cmap=plt.cm.coolwarm, aspect=1 )
            fig.colorbar(cb,   ax=ax)
            clim = kwargs.get('clim')
            if clim is not None : 
                cb.set_clim(clim)
        elif which_plot == 'imshow' : 
            if style_dic.get('extent') is None and (len(xx)>1) and (len(yy)>1): #extent = [min(xx)-0.05, max(xx) + 0.05, min(yy), max(yy)]
                dx, dy = old_div((max(xx)-min(xx))*1.0,(len(xx)-1)), old_div((max(yy)-min(yy))*1.0,(len(yy)-1))
                #dx2, dy2 = dx/2.,  dy/2.
                #style_dic['extent'] = [min(xx)-dx2, max(xx)+dx2, min(yy)-dy2, max(yy)+dy2]
                style_dic['extent'] = [min(xx), max(xx)+dx, min(yy), max(yy)+dy]
            im = ax.imshow(Z,  origin='lower', **style_dic) 
                    
            #cb=fig.colorbar(im,  orientation="horizontal", ax=ax)
            cb=fig.colorbar(im,   ax=ax)
            clim = kwargs.get('clim')
            if clim is not None : 
                cb.set_clim(clim)
        elif which_plot == 'surface': 
            temp = ['interpolation', 'extent', 'aspect']  #these are not accepted by plot_surface 
            for t in temp: 
                style_dic.pop(t)
            ax.plot_surface(XX, YY, data, rstride=1, cstride=1, **style_dic)
        else: 
            raise ValueError("alowed which_plot=%s"%(['contourf', 'imshow', 'surface']))
       
        if 'title' in kwargs and ax.is_first_row(): 
            ax.set_title(kwargs.get('title')) 
        
        temp = ['xlabel', 'ylabel', 'xlim', 'ylim', 'xscale', 'yscale']
        for t in temp: 
            tt = kwargs.get(t)
            if tt is not None :
                if isinstance(tt, tuple): 
                    ax.__getattribute__('set_' + t)(*tt)
                else: 
                    ax.__getattribute__('set_' + t)(tt)
      
        ax.grid(1)
        #ax.invert_yaxis()
        #res={'data':data.T, 'alpha':XX,'g':YY, 'cb':cb}
        res= {'fig': fig, 'ax': ax, 'cb': cb}
        return res 
    

    def search_field_name(self, name, only_first=False, only_one=False, assert_found=False): 
        """
           
        """
        all_field_names =  self.ALL_FIELD_NAMES 
        res = self.search_str(name, all_field_names, only_first=only_first, 
                only_one = only_one, assert_found = assert_found,) 
        return res 

    def get_rec_1d(self, field_name, aa, sh, sub_key_list=None, rec_getter=None, rec_getter_args=None, **kwargs): 
        info = kwargs.get('info', 0)
        
        aa = [a if hasattr(a, '__iter__') else (a, ) for a in aa]
        param_name, param_name_id = AlphaList.identify_param.__func__(self, aa)    
        data=[]
        not_found = []
        for a in aa:
            db=self[a]
            if rec_getter is None: 
                rec=db.fetch_easy(field_name, sh, sub_key_list=sub_key_list)
            else: 
                rec_getter_args= rec_getter_args if rec_getter_args is not None else {}
                rec = rec_getter(db, sh, **rec_getter_args)
            
            a1 = a[param_name_id]
            if rec is not None:
                data.append((a1, rec))
            else:
                data.append((a1, np.nan))
                not_found.append(a)
                
        if not_found and info>0: 
            print('not_found is ', not_found)
        try: 
            x,y=list(zip(*data))        
        except ValueError as err: 
            print(err)
            x, y = np.nan, np.nan 
        x = np.asarray(x); y=np.asarray(y)
        return x,  y

    def plot_field_vs_alpha(self, field_str, aa=None, sh_list=None,
            sub_key_list=None, alpha_db_map=None, empty_to_nan=True,  
            rec_getter=None,  rec_getter_args=None, 
            param_name=None,  fault_tol=1,  **kwargs): 
        
        #field_name = self.result_db_class.search_field_name.im_func(
        #        self.result_db_class, field_str, only_first=False, only_one=True,  
        #        assert_found=True)[0]
        
        field_name = self.search_field_name( field_str, only_first=False, only_one=True,  
                assert_found=False)
        if len(field_name)==1:
            field_name = field_name[0]
        else:
            field_name = field_str
        aa = aa if aa is not None else list(self.alpha_list)
        if isinstance(sh_list, tuple):  
            sh_list = [sh_list]
        sh_min = kwargs.get('sh_min')
        sh_max = kwargs.get('sh_max')
        sh_list = sh_list if sh_list is not None else self.get_shape_list_all(aa=aa, sh_min=sh_min, sh_max=sh_max)
        #sh_list = list(reversed(sorted(sh_list)))
        #sh_list.reverse()
        
        alpha_db_map = alpha_db_map if alpha_db_map is not None else {}
        rec_getter_args= rec_getter_args if rec_getter_args is not None else {}
        info = kwargs.get('info', 0)
        algorithm = self.algorithm 
        if algorithm == 'vmps':
            if field_name == 'entanglement': 
                field_name += '_entropy' 
        
        if param_name is None: 
            #param_name, param_name_id = AlphaList.identify_param.__func__(self, aa)
            param_name, param_name_id = AlphaList.identify_param(self, aa)
        else: 
            param_name_id = self.param_list.index(param_name)
        if param_name in ['alpha', 'beta']: 
            param_name_tex = '\\' + param_name
        else: 
            param_name_tex = param_name 
        
        kwargs_orig = kwargs.copy()
        not_found=[]
        for sh in sh_list: 
            
            if algorithm == 'vmps' and field_name == 'entanglement_entropy' : 
                skl = sub_key_list if sub_key_list is not None else ['EE', sh[0]//2, 1]
            else: 
                skl = sub_key_list 
                
            data=[]
            for a in aa:
                if a in alpha_db_map: 
                    db = alpha_db_map[a]
                else: 
                    db=self[a]
                if rec_getter is None: 
                    rec=db.fetch_easy(field_name, sh, sub_key_list=skl, fault_tolerant=fault_tol)
                else: 
                    if type(rec_getter)==types.MethodType:
                        rec_getter = rec_getter.__func__
                    rec = rec_getter(db, sh, **rec_getter_args)
                
                a1 = a[param_name_id]
                if rec is not None:
                    data.append((a1, rec))
                else:
                    if empty_to_nan: 
                        data.append((a1, np.nan))
                    not_found.append(a)
            if not_found and info>0: 
                print('not_found is ', not_found)
            if not data:
                warnings.warn('empty data for sh=%s'%(sh, ))
                continue
            try: 
                x,y=list(zip(*data))        
            except ValueError as err: 
                warnings.warn(str(err))
                x, y = np.nan, np.nan 
            x = np.asarray(x); y=np.asarray(y)
            
            label = kwargs_orig.get('label', str(sh))
            #marker = kwargs_orig.get('marker', MARKER_CYCLE.next())
            ylabel = field_name if rec_getter is None else (
                    rec_getter.__name__.replace('_get', ''))
            kwargs.update(ylabel=ylabel, return_ax=1, label=label, 
                     xlabel='$%s$'%param_name_tex)
            
            try: 
                fig, ax = self._plot(x, y, **kwargs) 
                
                kwargs['fig'] = fig 
                kwargs['ax'] = ax
            except Exception as err: 
                if fault_tol:
                    warnings.warn(str(err))
                else:
                    raise  
            
        if kwargs.get('return_fig'): 
            return fig 
    
    def plot_field_vs_alpha_3d(self, field_name, aa, sh, xparam, yparam, 
            rec_getter=None, rec_getter_args= None,  fault_tol=True, 
            zmin=None, zmax=None, 
           make_up=False, alpha_db_map=None,  sub_key_list=None, **kwargs): 
        info = kwargs.get('info', 0)
        aa = aa if aa is not None else self.alpha_list
        alpha_db_map = alpha_db_map if alpha_db_map is not None else {}
        rec_getter_args= rec_getter_args if rec_getter_args is not None else {}
        xx = self.get_param_range(xparam, aa=aa)
        yy = self.get_param_range(yparam, aa=aa)
         
        XX,YY=np.meshgrid(xx,yy)
        data=np.ndarray((len(xx),len(yy)))
        data[: ] = np.nan 
        alg = self.algorithm
        not_found = []
        p1=self.param_list.index(xparam) 
        p2=self.param_list.index(yparam) 
        x_to_i = {x:xx.index(x) for x in xx}
        y_to_j = {y:yy.index(y) for y in yy}
        if info>0: 
            grid = itertools.product(xx, yy)
            grid1 = [(a[p1], a[p2]) for a in aa]
            temp=set(grid)-set(grid1)
            temp = list(temp)
            temp.sort()
            if temp: 
                print('points not in alpha_parpath_dict is', temp)
        
        for a in aa: 
            x = a[p1]; y = a[p2]
            i = x_to_i[x]; j = y_to_j[y]
            if a in alpha_db_map: 
                db = alpha_db_map[a]
            else: 
                db=self[a]
            
            if rec_getter is None: 
                rec=db.fetch_easy(field_name, sh, sub_key_list=sub_key_list, fault_tolerant=fault_tol)
            else: 
                rec = rec_getter(db, sh, **rec_getter_args)
            #rec = db.fetch_easy(field_name, sh, sub_key_list, info=0)
            
            if rec is not None : 
                data[i, j] = rec 
            else:
                data[i, j] = np.nan 
                not_found.append((x, y))
        
        if make_up and 0:  #and len(not_found)>0: 
            make_up_list = []
            m, n = data.shape 
            #print np.where(np.isnan(data))
            aaa, bbb = np.where(np.isnan(data))
            ab = list(zip(aaa.tolist(), bbb.tolist()))
            for a, b in ab: 
                if a>0 and a <m: 
                    data[a, b] = (data[a-1, b]  +  data[a+1, b])/2. 
                    make_up_list.append((a, b))
                elif a == 0: 
                    data[a, b] = data[a+1, b]
                elif a == m-1: 
                    data[a, b] = data[a-1, b]
        kwargs['make_up'] = make_up 
        kwargs['zmin'] = zmin
        kwargs['zmax'] = zmax
        if not_found and info>0: 
            print('not_found is ', not_found) 
        kwargs.update(xlabel=xparam, ylabel=yparam)
        
        #if len(xx)>0 and len(yy)>0: 
        if 1: 
            dic = self._plot3d(xx, yy, data, **kwargs)
            fig, ax, cb = dic['fig'], dic['ax'], dic['cb']
        else: 
            pass 
        
        #if kwargs.get('return_fig'): 
        #    return fig 
        return dic
        
    def get_rec_2d(self, field_name, aa, sh, xparam, yparam, 
            rec_getter=None, rec_getter_args= None,  fault_tol=True, 
           make_up=False, alpha_db_map=None,  sub_key_list=None, **kwargs): 
        info = kwargs.get('info', 0)
        aa = aa if aa is not None else self.alpha_list
        alpha_db_map = alpha_db_map if alpha_db_map is not None else {}
        rec_getter_args= rec_getter_args if rec_getter_args is not None else {}
        xx = self.get_param_range(xparam, aa=aa)
        yy = self.get_param_range(yparam, aa=aa)
         
        XX,YY=np.meshgrid(xx,yy)
        data=np.ndarray((len(xx),len(yy)))
        data[: ] = np.nan 
        alg = self.algorithm
        not_found = []
        p1=self.param_list.index(xparam) 
        p2=self.param_list.index(yparam) 
        x_to_i = {x:xx.index(x) for x in xx}
        y_to_j = {y:yy.index(y) for y in yy}
        if info>0: 
            grid = itertools.product(xx, yy)
            grid1 = [(a[p1], a[p2]) for a in aa]
            temp=set(grid)-set(grid1)
            temp = list(temp)
            temp.sort()
            if temp: 
                print('points not in alpha_parpath_dict is', temp)
        
        for a in aa: 
            x = a[p1]; y = a[p2]
            i = x_to_i[x]; j = y_to_j[y]
            if a in alpha_db_map: 
                db = alpha_db_map[a]
            else: 
                db=self[a]
            
            if rec_getter is None: 
                rec=db.fetch_easy(field_name, sh, sub_key_list=sub_key_list, fault_tolerant=fault_tol)
            else: 
                rec = rec_getter(db, sh, **rec_getter_args)
            #rec = db.fetch_easy(field_name, sh, sub_key_list, info=0)
            
            if rec is not None : 
                data[i, j] = rec 
            else:
                data[i, j] = np.nan 
                not_found.append((x, y))
        
        if not_found and info>0: 
            print('not_found is ', not_found) 
        
        return xx, yy, data

    def show_fig(self): 
        plt.show()
    
    def _scaling_dim_exact(self, model_name):
        exact = {}

        exact["xx"] = {
                0:(0.0, ) + (1.0, )*4   +  (2.0, )*8, 
                1:(0.25,) + (1.25, )*4  +  (2.25, )*8 , 
                2:          (1.0, )*4  + (2.0, )*8}
        
        exact["Ising"] = {
                1:(0, ) + (1, )   +  (2, )*4, 
                -1:(1/8.,) + (1 + 1/8., )*2  +  (2+1/8., )*3 
                }
        exact["heisbg_long"] = {
                0:(0, ) + (0.5, )*2  + (1.0, )*2  + (1.5, )*2, 
                1:(0.5, ) + (1, )*2  + (1.5, )*2, 
                2:(2.0, )*3, 
                }
        exact["heisbg"] = exact['heisbg_long']

        #ee = self.obj.energy_exact
        #mapper = {-1.2732:"xx", -1.6077:"heisbg_NNN", -1.6449:"heisbg_long", -1.7725:'heisbg'}
        #for i in mapper:
        #    if abs(i-ee)<1e-2:
        #        model = mapper[i]
        #exact_val = exact[model]
        #if self.obj.model == "Ising": 
        #    exact_val = exact["Ising"]

        return exact[model_name]
    
    get_scaling_dim_exact = _scaling_dim_exact
    
    def plot_group(self, an_meth, args_group, figsize=None):
        n = len(self.alpha_group)
        fs1 = figsize if figsize is not None else (4, 3)
        fs= fs1[0]*n, fs1[1]
        fig = plt.figure(figsize=fs)
        #show_title = args.get('show_title')
        i = 0
        for aa in self.alpha_group: 
            args= args_group[i]
            i += 1 
            #db = self.alpha_rdb_dict[a]
            #aa = self.alpha_group[g]
            fig.add_subplot(1, n, i)
            args['fig'] = fig
            an_meth.__call__(aa, **args)
        
            #if show_title: 
            if 0: 
                ax = fig.axes[-1]
                ax.set_title(str(a))
        fig.tight_layout()
 
    def plot_collumn(self, aa, db_meth, args):
        n = len(aa)
        for a in aa: 
            db = self.alpha_rdb_dict[a]
            #fig = plt.figure(figsize=(4*n, 3))
            #args['kwargs'] = {'fig': fig}
            #args['fig'] = fig
            #print 'aaa', args.keys()
            #db_meth.im_func.__call__(self=db, **args, kwargs={'fig': fig})
            db_meth.__func__.__call__(db, **args)
       
    def plot_row(self, aa, db_meth, args):
        n = len(aa)
        fig = plt.figure(figsize=(4*n, 3))
        show_title = args.get('show_title')
        i = 0
        for a in aa: 
            i += 1 
            db = self.alpha_rdb_dict[a]
            fig.add_subplot(1, n, i)
            args['fig'] = fig
            db_meth.__func__.__call__(db, **args)
        
            if show_title: 
                ax = fig.axes[-1]
                ax.set_title(str(a))
        fig.tight_layout()
        
    def plot_square(self, meth, args):
        pass
   
    def plot_scaling_dim(self, qns=None, num_of_fields=None, calc_err=False):
        """
            one rec is like this:
            {1025:
                [(0, [0.0, 1.0, 1.16326228, 1.19447544, 1.19447544, 1.23058869, 1.6277161, 1.6277161, 1.82188553, 1.82188553]),
                 (1, [0.87008947, 1.08398164, 1.12192732, 1.50196967, 1.63179961, 1.74482795, 1.88535985, 2.0160677, 2.13870443, 2.13870443]),
                 (2, [1.0722275, 1.78260076, 2.37668506, 2.60884736, 2.60884736, 2.800528, 2.800528, 2.85390176, 2.93886244, 2.93886244])
                 ]
             } or
             {1025:
                [(1, [-0.0, 1.00009277, 1.99996219, 1.99997006, 2.00000392, 2.00289768, 2.72124425, 2.72124425, 2.72266868, 2.72266868]), 
                [(-1, [0.1250008, 1.12503901, 1.12508933, 2.12244868, 2.12717037, 2.12830034, 2.35310428, 2.35310428, 2.35354695, 2.35354695])}
            
        """
        num_of_fields = num_of_fields if num_of_fields is not None else 2
        try:
            rec = self.obj.scaling_dim_record
        except AttributeError:
            rec = self.obj.conformal_dim_record
        
        if calc_err:
            exact_val = self._scaling_dim_exact()
            print(exact_val)

        x = list(rec.keys())
        x.sort()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax1 = fig.add_subplot(111)
        
        xy = self.get_energy_err()
        ax1.plot(*xy, linestyle="dashed")

        xy = self.get_energy_diff()
        ax1.plot(*xy, linestyle="dashed")
        #ax1.xlim(xmin=x[0])
        
        if qns is None:
            if self.obj.symmetry == "U1":
                qns = list(range(3))
            elif self.obj.symmetry == "Z2":
                qns= [1, -1]
            else:
                qns= [1]
        fontP = FontProperties()
        fontP.set_size('small')
   
        legend = ["eng_err"]
        color_map = {0:"red", 1:"blue", 2:"green"}
        line_sty_map = {0:"solid", 1:"dashed", 2:"dotted"}

        for qn in qns:
            for i in range(num_of_fields):
                qn_ = qns.index(qn)
                if qn_ == 0 and num_of_fields== 1: i = 1
                linewidth = 3.5 if (qn_, i) in [(0, 1), (1, 0), (2, 0)] else 1.0
                y1 = np.array([(rec[ii][qn_][1][i]) for ii in x])
                if calc_err:
                    try:
                        y1 = y1-exact_val[qn][i]
                        y1 = np.log10(np.abs(y1))
                    except IndexError:
                        print("index qn=%d, i=%d not found in exact_val, omit it "%(qn, i))
                        continue
                ax.plot(x, y1, "o", linewidth=linewidth,  
                        color=color_map[qn_], linestyle="solid", marker=None)
                
                #if (qn_, i) in [(0, 1), (1, 0), (2, 0)]:
                legend.append("%d, %d"%(qn, i))
        box = ax.get_position()
        box1 = ax1.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax1.set_position([box.x0, box1.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(legend, loc='center left', bbox_to_anchor=(1, 0.5), prop=fontP)
        plt.xlim(xmin=x[0])
        #plt.ylim(ymin=-0.1)
        #plt.tick_params()
        #plt.minor_tickson()
        #plt.legend(legend, loc=3)
        plt.grid(True)
        
        print("ssss\n", rec[x[-1]])
        if hasattr(self.obj, "log"):
            print("llll\n", self.obj.log)
        x = "_err" if calc_err else ""
        y = self.fn_only.split(".")[0]
        fn=self.path + "scaling_dim-%s--%s.png"%(x, y)
        plt.savefig(fn ,bbox_inches='tight') # save as png
        os.popen("display " + fn)

    def plot_energy(self, aa=None, sh=None, param_name=None, **kwargs): 
        if aa is None: 
            aa = list(self.alpha_list)
        aa_dic = self.preprocess_alpha_list(aa)
        
        for sh in [sh]:
            x = []
            y=[]
            not_found = []
            #for a in aa:
            #print 'aaa', aa_dic
            for a , db in aa_dic.items(): 
                
                #e= an.alpha_rdb_dict[a].fetch('energy',dim=12)
                #db = self.alpha_rdb_dict[a] 
                e= db.fetch_easy('energy',sh)
                if param_name is not None : 
                    a = a[self.param_list.index(param_name)]
                if e is not None :
                    
                    x.append(a)
                    y.append(e)
                else: 
                    not_found.append((a, sh))
            if len(not_found)>0: 
                print('not_found is', not_found)
            #print a,  x, y
            #ax.plot(x, y, '.-')
            fig=self._plot(x, y, **kwargs)
        if kwargs.get('return_fig'): 
            return fig
        
    def get_energy_err(self):
        """
            one record is like this:
            {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        """
        recs = self.obj.energy_record
        ee = self.obj.energy_exact
        x = list(recs.keys())
        x.sort()
        y = [recs[i][0] for i in x]
        y = np.array(y)
        y = np.log10(y-ee)
        #print recs[1]
        #self._plot(x, y,"o", color='red',linestyle='dashed')
        return x, y
    
    def plot_energy_diff_vs_time(self):
        """
        one record is like this:
        {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        
        """
        #recs= self.obj.energy_record
        res = System.load(self.path)
        recs = res.energy_record
        if 1:
            t_eng = [(t, v[0])  for t, v in recs.items()]
            #Et_diff = [(i[0], i[]) for i in t_eng]

            t, eng = list(zip(*t_eng))
            t_diff = [(t_eng[i+1][0], old_div((t_eng[i+1][1]-t_eng[i][1]),(t_eng[i+1][0]-t_eng[i][0]))) 
                    for i in range(len(t_eng)-1)]
            t, diff = list(zip(*t_diff))
            diff = np.array(diff)
            diff = np.log10(abs(diff))
        if 0:
            plt.plot(t, diff, "o", color='red',linestyle='dashed')
            plt.show()
        if 1:
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            
            xx = 1e-3
            e = np.mean(eng[1000:])
            ax1.set_ylim(ymin=e-xx, ymax=e+xx)

            ax1.plot(t, eng[1:])
            
            ax2.plot(t, diff, "o", color='red',linestyle='dashed')
            plt.show()

    def plot_energy_diff_vs_chi(self, aa, sh_list, sh_eng_min, eng_exact=None, log_scale=1, return_fig=0): 
        fig=plt.figure()
        ax=fig.add_subplot(111)
        aa=list(aa)
        aa.sort(reverse=1)
        sh_list.sort()
        
        for a in aa:
            db=self.alpha_rdb_dict[a]
            eng_min = db.fetch_easy('energy', sh_eng_min)
            if eng_exact is not None : 
                eng_min = eng_exact
            if eng_min is None: 
                continue
            y=[]
            x = []
            for sh in sh_list:
                #db=self.alpha_rdb_dict[a]
                e= db.fetch_easy('energy', sh)
                if e is None: 
                    continue
                x.append(sh[0])
                y.append(e-eng_min)
            
            label = '%s'%(a, )
            if log_scale: 
                ax.set_yscale('log')
                ax.set_xscale('log')
            ax.plot(x, y, '.-', label=label)
        #ax.invert_xaxis()
        ax.set_ylabel('energy diff',fontsize=18)
        ax.set_xlabel('$dim$',fontsize=18)
        ax.grid(1)
        ax.legend()
        if return_fig: 
            return fig

    def plot_energy_diff(self, aa, sh_list, sh_eng_min, log_scale=1, return_fig=0): 
        fig=plt.figure()
        ax=fig.add_subplot(111)
        aa=list(aa)
        aa.sort(reverse=1)
     
        for sh in sh_list:
            
            y=[]
            x = []
            for a in aa:
                db=self.alpha_rdb_dict[a]
                eng_min = db.fetch_easy('energy', sh_eng_min)
                if eng_min is None: 
                    continue
                #db=self.alpha_rdb_dict[a]
                e= db.fetch_easy('energy', sh)
                if e is None: 
                    continue
                x.append(a)

                y.append(e-eng_min)
            
            label = '%s'%(sh, )
            if log_scale: 
                ax.set_yscale('log')
            ax.plot(x, y, '.-', label=label)
        ax.invert_xaxis()
        ax.set_ylabel('energy diff',fontsize=18)
        ax.set_xlabel('$\\alpha$',fontsize=18)
        ax.grid(1)
        ax.legend()
        if return_fig: 
            return fig
        
    def get_energy_diff(self):
        """
        one record is like this:
        {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        
        """
        recs= self.obj.energy_record
        x = list(recs.keys())
        x.sort()
        dy = np.array([recs[i+1][0] - recs[i][0] for i in x[:-2]])
        dy = np.log10(-dy)
        #print recs[1]
        #self._plot(x[:-2], dy,"o", color='red',linestyle='dashed')
        return x[:-2], dy

    def plotax_data(self, ax, data, label_list=None): 
        pass

    def plotax_dist_del(self, ax, data, x=None, transpose=1,  label_list=None, info=0): 
        if x is None: 
            x = [i for i in self.DISTANCE_LIST]
        if info>0: print('X=', x)

        data = np.array(data)
        if transpose: 
            data = data.T
        for i in range(data.shape[1]):  #differ with plotax_dist_old at here and next line
            y = data[:, i]  
            ax.loglog(x, y, '.-')
        
        if label_list is not None: 
            ax.legend(label_list, loc='lower right') 
        ax.set_xlabel('$r$', fontsize=18)
        ax.set_xlim(0, 1e5)
        ax.grid(1)
   
    def plotax_corr_extra(self, ax, a, mera_shape, direct): 
        db=self.alpha_rdb_dict[a]
        rec=db.fetch_easy('correlation_extra', mera_shape,[direct])
        if rec is None: 
            print('rec not found for a=%s'%a)
            return
        x,y=list(zip(*rec[0: -1: 2]))        
        y = np.array(y) #y=abs(y)
        if direct == 'pm':  y *= 2 
        ax.loglog(x, y, '.-', label=str(a))
    
    def plot_central_charge(self, aa, sh, alpha_remap={}, alpha_sh_map={}, 
            layer='top', **kwargs):
        
        #aa.sort(reverse=1)
        if layer=='top':
            ll= [sh[1]-2 ] 
        elif layer=='all':
            ll=  list(range(sh[1]-1))
        else: 
            ll = layer
      
        param_name, param_name_id = AlphaList.identify_param.__func__(self, aa)
        #if fig is None: 
        #    fig=plt.figure(figsize=(5, 3))
        #ax=fig.add_subplot(111)
        kwargs_orig = dict(kwargs)
        
        #ll= [sh[1]-2 ] if layer is None else range(sh[1]-1)
        not_found = []
        for layer in ll:
            cc=[]
            for a in aa:
                a_orig = a
                a1 = alpha_remap.get(a)
                a = a1 if a1 is not None else a_orig
               
                db=self.alpha_rdb_dict[a]
                sh1 = alpha_sh_map.get(a_orig)
                
                shape = sh1 if sh1 is not None else sh
                rec = db.fetch_easy('entanglement_entropy', shape)
                if rec is not None : 
                    rec = rec.get(layer)
                if rec is None: 
                    not_found.append((a, layer))
                    continue
                try: 
                    e1, e2= rec[1], rec[2]
                except: 
                    print(db.path)
                    raise
                c= (e2-e1)*3
                #cc.append((a_orig,c))
                cc.append((a_orig[param_name_id],c))
            if len(cc)>0:     
                x,y=list(zip(*cc))    
                x = np.array(x)
                y = np.array(y)
                label = str(layer) if 'label' not in kwargs_orig else kwargs_orig['label']
                marker = '.' if 'marker' not in kwargs_orig else kwargs_orig['marker']
                kwargs.update(marker=marker, label=label, return_ax=1)
                fig, ax=self._plot(x, y, **kwargs)
                #fig, ax=self._plot(x, y, **kwargs)
                kwargs['ax'] = ax
                kwargs['fig'] = fig
        if 'ax' in list(locals().keys()): 
            if ax.is_first_col(): 
                ax.set_ylabel('central change')
            if ax.is_last_row(): 
                ax.set_xlabel('$\\alpha$')
            #ax.invert_xaxis()
        if len(not_found)>0: 
            msg = 'rec for central_charge not_found is :%s'%not_found
            warnings.warn(msg)
        if kwargs.get('return_ax', False): 
            return ax 
        if kwargs.get('return_fig', False): 
            return fig
    
    def plot_central_charge_vs_layer(self, aa, sh, alpha_remap={}, alpha_sh_map={}, **kwargs):
        for a in aa: 
            
            try: 
                db = self[a]
                kwargs.update(return_ax=1, label=str(a) )#, fig=fig)
                fig, ax = db.plot_central_charge_vs_layer(sh=sh, **kwargs)
                if ax is not None :  #plot succecss
                    kwargs['ax'] = ax
                if fig is not None : 
                    kwargs['fig'] = fig
            except Exception as e: 
                msg = 'error: %s %s '%( a, e )
                print(msg, end=' ') 
                #raise(e)
            
        if kwargs.get('return_fig', False): 
            return fig
    
    def plot_correlation_bac(self, alpha_list, mera_shape_list, alpha_remap={}, alpha_shape_dic={},  
            direct='pm', which='correlation_extra', fig=None,  log_scale=1,  return_fig=0, draw=0, **kwargs): 
        assert which in ['correlation_extra', 'correlation', 'correlation_mixed']
        alpha_list.sort(reverse=1)
        if fig is None : 
            fig=plt.figure(figsize=(6, 4))
        if len(fig.axes)==0: 
            ax=fig.add_subplot(111)
        else: 
            ax = fig.axes[-1]
        not_found=[]
        
        if len(mera_shape_list)==2 and isinstance(mera_shape_list, tuple): 
            mera_shape_list = [mera_shape_list]
            
        sub_key_list = [direct] if which == 'correlation_extra' else [direct, 'val'] 
        for a in alpha_list:
            a1 = alpha_remap.get(a)
            a = a1 if a1 is not None else a
            
            db=self.alpha_rdb_dict[a]
            for sh in mera_shape_list: 
                sh_temp = alpha_shape_dic.get(a)
                if sh_temp is not None :  sh = sh_temp
                rec=db.fetch_easy(which, mera_shape=sh, sub_key_list=sub_key_list)
                if rec is None: 
                    not_found.append((a, sh))
                    continue
                
                if which  == 'correlation_extra' and log_scale:  
                    x,y=list(zip(*rec[0:-1:2]))
                else: 
                    x, y = list(zip(*rec[: ]))
                if direct == 'pm':   
                    y=np.array(y);  y=y*2;  
                
                label = ''
                if len(mera_shape_list)<2: 
                    label += ' %s  '%(a, )
                if len(alpha_list)<2: 
                    label += ' %s  '%(sh, )
                if log_scale: 
                    y=np.abs(y)
                    ax.loglog(x,y, '--', label=label)
                else: 
                    #ax.plot(x,y, '--', label=label)
                    ax.plot(x,y,  label=label)
            
        l=ax.legend(title='$\\alpha$')
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transparent
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
           
        ax.grid(1)
        if draw: 
            plt.show()
        if len(not_found)>0:
            print('not_found is ', not_found)
        if return_fig:
            return fig

    def plot_correlation(self, alpha_list, sh, direct='pm', which='correlation_extra', 
            r_min=None, r_max=None, alpha_remap={}, alpha_shape_dic={}, 
            fig=None,  log_scale=1,  return_fig=0, draw=0, **kwargs): 
        assert which in ['correlation_extra', 'correlation', 'correlation_mixed']
        #kwargs['log_scale'] = log_scale
        mera_shape_list = sh
        param_name = kwargs.get('param_name')
        #if self.invert_xaxis: 
        alpha_list = list(alpha_list)
        if self.invert_alpha_order: 
            alpha_list.sort(reverse=1)
        if fig is None : 
            fig=plt.figure(figsize=(6, 4))
        if len(fig.axes)==0: 
            ax=fig.add_subplot(111)
        else: 
            ax = fig.axes[-1]
        not_found = []
        
        if len(mera_shape_list)==2 and isinstance(mera_shape_list, tuple): 
            mera_shape_list = [mera_shape_list]
        mark_list=['d', '^', '>', '<', 'o', 'D', 'H', '.', '*']
        mark_list += mark_list            
        count = -1
        for a in alpha_list:
            count += 1 
            a1 = alpha_remap.get(a)
            a = a1 if a1 is not None else a
            sh  =  mera_shape_list
            sh_temp = alpha_shape_dic.get(a)
            db=self.alpha_rdb_dict[a]
            if sh_temp is not None :  
                sh = sh_temp
            if param_name is not None : 
                param = a[self.param_list.index(param_name)]
                label = str(param)
            else: 
                label = str(a)
            if 'marker' not in kwargs:
                kwargs['marker'] = mark_list[count]
            db.plot_correlation(sh, direct=direct, r_min=r_min, r_max=r_max, which=which, fig=fig, log_scale=log_scale, 
                    label=label, **kwargs)    
        
        #l=ax.legend(title='$\\alpha$')
        #l=ax.legend(alpha_list, title='$\\alpha$')
        l=ax.legend(title='$\\alpha$')
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transparent
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
           
        ax.grid(1)
        if draw: 
            plt.show()
        if len(not_found)>0:
            print('not_found is ', not_found)
        if return_fig:
            return fig

    def plot_entanglement_vs_alpha(self, aa, sh, alpha_remap={}, sh_remap=None, 
            site_list=None, 
            layer=0, switch='entropy', **kwargs):
        """
            switch in [entropy, extra, brute_force, brute_force_6, brute_force_9]
        """
        #-----------------------------------------------------------------------
        #print 'ploting switch=%s'%(switch, )
        if 1: 
            if alpha_remap is None: 
                alpha_remap = {}
            if sh_remap is None: 
                sh_remap = {}
            fault_tol = kwargs.get('fault_tol', True)
            info = kwargs.get('info', 0)
            aa=list(aa)
            aa.sort(reverse=1)
            
            fig = kwargs.get('fig')
            param_name = kwargs.get('param_name')
        
        if 1:  #-----------------------------------------------------------------------
            algorithm = self[aa[0]].get('algorithm')
            if algorithm not in  ['mps', 'idmrg']:  
                if switch == 'entropy': 
                    site_list = [1, 2] if site_list is None else site_list
                    sub_key_list = [layer]
                elif switch == 'extra': 
                    sub_key_list = ['EE']
                    site_list = [0] if site_list is None else site_list
                elif switch in ['brute_force', 'brute_force_6', 'brute_force_9']: 
                    sub_key_list = ['EE']
                    site_list = [10, 12] if site_list is None else site_list
            elif algorithm == 'mps': 
                    site_list = [sh[0]//2] if site_list is None else site_list
            elif algorithm == 'idmrg' : 
                site_list = [0]
            
        for i in site_list: 
            not_found = []
            x = []
            y = []
            for a in aa:
                a_orig = a
                a1 = alpha_remap.get(a)
                a = a1 if a1 is not None else a_orig
                db=self.alpha_rdb_dict[a]
                if algorithm  == 'mps' : 
                    pass
                    name = 'entanglement_entropy'
                    sub_key_list_1 = ['EE', i, 1]
                elif algorithm == 'idmrg': 
                    name = 'entanglement'
                    sub_key_list_1 = []
                    
                else: #mera 
                    if switch in ['entropy']  + ['brute_force', 'brute_force_6', 'brute_force_9']: 
                        sub_key_list_1 = sub_key_list +  [i] 
                    else: 
                        sub_key_list_1 = sub_key_list
                    #e=db.fetch_easy('entanglement_' + switch, sh, sub_key_list=sub_key_list_1, fault_tolerant=1)
                    name = 'entanglement_' + switch
               
                
                shape = sh_remap[a] if a in sh_remap else sh
                e=db.fetch_easy(name, shape, sub_key_list=sub_key_list_1, 
                        info=info, fault_tolerant=fault_tol)
                
                if param_name is not None : 
                    a_orig = a[self.param_list.index(param_name)]
              
                if e is not None : 
                    x.append(a_orig)
                    y.append(e)
                else: 
                    not_found.append(a)
            
            #label = '%d site'%i if kwargs.get('label') is None else kwargs.get('label')
            #args = dict(marker='.', return_ax=1, label=str(i), title='at layer=%d'%layer)
            args = dict(marker='.', return_ax=1) 
            args.update(kwargs)
            #print 'xxxx', x 
            fig, ax=self._plot(x, y, **args)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
            
        if len(not_found)>0 and info>0: 
            print('not_found is %s'%(not_found, ))
        
        if self.invert_xaxis: 
            ax.invert_xaxis()
        if ax.is_first_col(): 
            ax.set_ylabel('entanglement',fontsize=18)
        if ax.is_last_row(): 
            xlable = '$\\alpha$'
            ax.set_xlabel('$\\alpha$')
            if param_name is not None : 
                xlabel = '$%s$'%param_name 
            if 'xlabel' in kwargs: 
                xlabel = kwargs['xlabel']
                ax.set_xlabel(xlabel, fontsize=18)
        #ax.set_title('at layer=%d'%layer)
        
        ax.grid(1)
        ax.legend(loc='lower left')
        
        if kwargs.get('return_fig'): 
            return fig

    def plot_entanglement_vs_layer(self, aa, sh): 
        aa=[a for a in aa]
        aa.sort(reverse=1)
        nsite=1
        fig=plt.figure()
        ax=fig.add_subplot(111)
        for a in aa:
            data=[]
            db=self.alpha_rdb_dict[a]
            rec=db.fetch_easy('entanglement_entropy', mera_shape=sh, info=0)
            if rec is None: continue
            data = [i[1] for i in list(rec.values())]
            ll=list(rec.keys())
                
            ax.plot(ll, data, '.-', label=str(a))
        ax.legend()    
        ax.grid(1)
        ax.set_xlim(0, 7)
        ax.set_ylim(0, 2.0)
        
    def plot_entanglement_scaling(self, aa, sh, log2=1, aver=0, num_site=None,  
            site_start=0, site_period=1, site_end=None, 
            linear_fit=True, show_res=0, **kwargs): 
        aa=list(aa)
        aa.sort(reverse=1)
        not_found = []
        args= kwargs.copy()
        for a in aa:
            db=self.alpha_rdb_dict[a]
            if 0: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
            label = args['label'] if 'label' in kwargs else str(a) 
            args.update(return_ax=1, label=label)
            res = db.plot_entanglement_scaling(alpha=a, sh_list=[sh], aver=aver, 
                    log2=log2, site_start=site_start, site_period=site_period, 
                    site_end=site_end, 
                    num_site=num_site, linear_fit=linear_fit, 
                    show_res=show_res, **args)
            fig, ax = res 
            if ax is not None :  #plot succecss
                args['ax'] = ax
            if fig is not None : 
                args['fig'] = fig
            
        if kwargs.get('return_fig'): 
            return fig
    
    plot_entanglement_vs_length = plot_entanglement_scaling
    
    def plot_entanglement_vs_chi(self, aa=None, **kwargs): 
        if aa is None: 
            aa = list(self.alpha_list) 
        
        args= kwargs.copy()
        for a in aa:
            db=self[a]
            label = args['label'] if 'label' in kwargs else str(a)
            args.update(return_ax=1, label=label)
            fig, ax =db.plot_entanglement_vs_chi(  **args)
            if ax is not None :  #plot succecss
                args['ax'] = ax
            if fig is not None: 
                args['fig'] = fig
        
        if kwargs.get('return_fig'): 
            return fig 
   
    def plot_entanglement_spectrum(self, aa, sh, sub_key_list=None,  which_level=None, 
            num_level=None, param_name=None, **kwargs): 
        num_level= 20  if num_level is None else num_level
        if sub_key_list is None: 
            pass 
        if self.algorithm == 'vmps': 
            which = 'entanglement_entropy'
        elif self.algorithm == 'idmrg' : 
            which = 'entanglement_spectrum'
        else: 
            raise ValueError(self.algorithm)
        
        if which_level is not None : 
            num_level = len(which_level)
        ES=np.ndarray((len(aa), num_level))
        for j,a in enumerate(aa):
            db=self[a]
            es=db.fetch_easy(which, sh, sub_key_list)
            if es is not None:
                if which_level is not None : 
                    ES[j]=es[which_level]
                else: 
                    #len(es) may smaller than num_level 
                    ES[j,:len(es)]=es[:num_level]
            else:
                ES[j,:]=np.nan
        
        param_name, param_id = AlphaList.identify_param.__func__(self, aa) 
        x=[_[param_id] for _ in aa]
        kwargs.update(return_ax=1)
        fig, ax=self._plot(x, -np.log10(ES), **kwargs)
        ax.set_xlabel(param_name)
        ax.set_ylabel('ES')
    
    def plot_magnetization(self, aa=None, sh=None, direct='z', pos=0, **kwargs): 
        if aa is None: 
            aa = list(self.alpha_list)
        if self.algorithm == 'mera': 
            key_list = [0]
        elif self.algorithm == 'idmrg': 
            #key_list = None   #crosswhite version should use this,  below for mcculloch 
            key_list = [direct, pos]
        else: 
            raise NotImplemented
        
        x=[]; y=[]
        not_found = []
        for a in aa:
            db=self.alpha_rdb_dict[a]
            rec=db.fetch_easy('magnetization', sh, sub_key_list=key_list)
            if rec is None: 
                not_found.append(a)
                continue 
            x.append(a)
            y.append(rec)
            #print 'rrrrr', db.keys(), db['magnetization']
            #raise 
        if not_found: 
            print('not_found is ', not_found)
        
        fig = self._plot(x, y, marker='.', **kwargs)
        if kwargs.get('return_fig') : 
            return fig
            
              
    def show_scaling_dim(self, aa, sh_list, with_header=1,  alpha_remap={}, alpha_shape_dic={}): 
        data=[]
        aa = list(aa)
        
        for a in aa:
            for sh in sh_list: 
                a_orig = a
                
                temp = alpha_remap.get(a)
                sh_temp = alpha_shape_dic.get(a)
                if temp is not None: 
                    a = temp
                db=self.alpha_rdb_dict[a]
                rec= db.fetch_easy('scaling_dim', sh)
                sdim= rec
                #print a, '\t', 0, sdim[0][:3], 1, sdim[1][:3], 2, sdim[2][0:1]
                #print a, '\t', 0, sdim[0][:3], 1, sdim[1][:3], 2, sdim[2][0:1]

                if rec is not None:
                    #data.append((a,)+sdim[0][:3]+ sdim[1][:3]+ sdim[2][:1])
                    data.append((str(a)+str(sh), ) + rec[0][:5] + rec[1][:3] +rec[2][:3] )
            
                    if with_header: 
                        headers = [getattr(self, 'name', '')[-12: ]]+[str(i) for i in [0,0,0,0,0,1,1,1,2]]
                    else: 
                        headers= None
        print(tabulate(data, headers=headers, tablefmt='simple'))

    def plot_eta_old(self, aa, eta, fig=None, return_fig=0):    
        y=eta
        #print aa, y
        if fig is None: 
            fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(aa, y, '.-')
        if self.invert_xaxis: 
            ax.invert_xaxis()
        #ax.set_xticklabels(aa)
        ax.grid(1)
        #no=ax.set_xticks(aa)    
        if return_fig: 
            return fig
    
    def _plot_eta(self, aa, eta, fig=None, return_fig=0):    
        y=eta
        #print aa, y
        if fig is None: 
            fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(aa, y, '.-')
        if self.invert_xaxis: 
            ax.invert_xaxis()
        ax.grid(1)
        #ax.set_xticklabels(aa)
        #no=ax.set_xticks(aa)    
        if return_fig: 
            return fig

    def _plot_eng_fluc_count(self, aa, sh, fault_tol=1, **kwargs):
        """
            not for general use, only for trouble shoting,  so add an under core
        """
        bins=[-1, 0,  1e-8, 1e-7, 1e-6, 1e-5, 1]
        aa=list(aa)
        temp={}
        for a in aa:
            db=self[a]
            #diff=db.get_energy( sh=sh, which='diff', iter_max=None)
            #diff_cut=np.histogram(diff, bins=bins)[0]
            #print a , diff_cut[1:]
            #temp[a]=diff_cut[1:]
            try: 
                count, bin = db.get_energy_fluc_count(sh=sh)
            except Exception as err: 
                if fault_tol: 
                    warnings.warn(str(err))
                    continue 
                else: 
                    raise 
            temp[a] = count[1: ]
        df=pd.DataFrame(temp, index=['<1e-8', '<e-7', '<1e-6', '<1e5', '<1'])
        #df=pd.DataFrame(temp, index=[str(i) for i in bins[1: ]])
        #df.unstack()
        #df.keys()
        df_trans = df.transpose()
        df_trans.plot( ax=kwargs.get('ax'))
        return df_trans

    #del this
    def plot_eta(self, aa, sh, direct='pm', func_type='power', **kwargs):    
        if 0:  
            fig = plt.figure()
            ax = fig.add_subplot(111)
        aa = list(aa)
        eta=[]
        for a in aa:
            db=self.alpha_rdb_dict[a]
            xx=db.fit_correlation(sh, direct=direct, func_type=func_type, plot=0)
            if xx is not None : 
                eta.append((a[2], xx[1]))
        x, y = list(zip(*eta))   
        self._plot(x, y, marker='.', **kwargs)
        #ax.plot(x, y, '.-')  
    
    def reload_class(self, info=1): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        importlib.reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        if info>0: 
            print('class %s is reloaded'%(name, ))
    
    def reload_db_class(self, info=1):
        """
            this method is not perfect, it is due to limitation of python
            itself. After reload, old objects are still referred to.
            https://www.dup2.org/node/950 
            https://docs.python.org/3/library/imp.html 
            When reload(module) is executed:
                Python modules’ code is recompiled and the module-level code reexecuted, defining a new set of objects which are bound to names in the module’s dictionary. The init function of extension modules is not called a second time.
                As with all other objects in Python the old objects are only reclaimed after their reference counts drop to zero.
                The names in the module namespace are updated to point to any new or changed objects.
                Other references to the old objects (such as names external to the module) are not rebound to refer to the new objects and must be updated in each namespace where they occur if that is desired.            
                
                
            #sys.getrefcount(an_vmps.an_main_2site.result_db_class)
        
        """
        #print 'old 11', id(result_db_module)
        module_name = self.result_db_class.__module__
        module = importlib.import_module(module_name)
        #reload(module)
        importlib.reload(module)
        class_name = self.result_db_class.__name__
        new_class = module.__getattribute__(class_name)
        self.result_db_class = new_class 
        if 0:  # below has no use 
            #setattr(result_db_module, class_name, new_class)
            print('old 2', id(result_db_module))
            print('new', id(module))
        if info>0: 
            print('%s is reloaded'%(self.result_db_class.__name__, ))
        return module 
        
    
    def plot_fidelity_2d(self, aa, data,  **kwargs): 
        dic = self._plot3d(aa, aa, data, **kwargs)
        fig, ax, cb = dic['fig'], dic['ax'], dic['cb']
        #fig.colorbar(cb, ax=ax) 
        #cb.set_clim(np.min(data), 1.0)
        #cb.set_clim(0, 1.0)
        
        ax.set_aspect(1)

    def plot_fidelity(self, aa, sh_list=None, delta=0.1, fault_tol=1, use_dist_comp=0, **kwargs):  
        sh_list = sh_list if sh_list is not None  else self.get_shape_list_all()
        
        kwargs_orig = kwargs.copy()
        param_name, param_id = AlphaList.identify_param.__func__(self, aa) 
        for sh in sh_list:
            data=self.calc_fidelity(aa[: -1], sh, 
                    delta=delta, fault_tol=fault_tol, force=kwargs_orig.get('force', False), 
                    use_dist_comp=use_dist_comp) #aa[: -1] means not calc fidelity for the last point 
            if data:
                label = kwargs_orig.get('label', str(sh))
                kwargs.update(return_ax=1, label=label, 
                         xlabel='$%s$'%param_name)
                
                x,y=list(zip(*data));  y=np.abs(np.asarray(y))
                _=self._plot(x,y, **kwargs) 
            

    def disregister_job(self, aa): 
        with rpyc_conn('local', 'service', 17012) as conn:     
            for a in aa:
                db=self[a]
                job_id=hash(os.path.normpath(db.parpath)  )
                try:
                    conn.root.disregister_job(job_id)
                except Exception as err:
                    print(err, 'param=', a) 

    def calc_fidelity(self, aa, sh, delta=0.1, force=0, use_dist_comp=False, 
            fault_tol=1, **kwargs):
        xx = self 
        if self.algorithm == 'vmps' : 
            from vmps.measure_and_analysis.measurement_vmps import fidelity
        elif self.algorithm == 'idmrg':
            from vmps.measure_and_analysis.measurement_idmrg_mcc import fidelity_direct  as fidelity
        else: 
            raise 
        
        dbs={}
        ss={}
        param_name, param_name_id = AlphaList.identify_param.__func__(xx, aa)
        data=[]
        
        which = 'fidelity'
        if delta is not None:
            db_change_list = []
            for a1 in aa:
                
                a2 = list(a1)
                ai = a1[param_name_id]
                a2[param_name_id] = round(ai + delta, 5)
                a2 = tuple(a2)
                
                db1 = xx[a1]
                rec = db1.fetch_easy(which, sh, sub_key_list=[a2], fault_tolerant=1) 
                if rec is not None and not force:  
                    data.append((ai, rec))
                else: 
                    #print_vars(vars(),  ['a1, a2'])
                    try:
                        db2 = xx[a2]
                        s1 = db1.load_S(sh)
                        s2 = db2.load_S(sh)
                        if s1 is None or s2 is None:
                            continue
                        #print_vars(vars(),  ['a2'])
                        if not use_dist_comp: 
                            ff=fidelity(s1, s2)
                        else: 
                            from mypy.brokest.brokest import queue
                            #print s1.keys(), s1['mpo'].additional_info.keys()
                            #raise 
                            s1['mpo'].additional_info.pop('J_power_fit')
                            s2['mpo'].additional_info.pop('J_power_fit')
                            ff = queue(fidelity, args=(s1, s2), kwargs={},  querry=0, host='210.45.121.30', port=90903)
                            #ff = queue(fidelity, args=(s1, s2), querry=0, port=90900)
                            print_vars(vars(),  ['a2', 'ff'])
                            ff = ff[1]
                        
                        if 1:
                            N, D = sh
                            D1 = db1.get_dim_max_for_N(N=N) if D == 'max' else D 
                            D2 = db1.get_dim_max_for_N(N=N) if  D == 'max' else D 
                            D = min(D1, D2)
                            db1.insert(which, (N, D),  sub_key_list=[a2], val=ff)
                            db2.insert(which, (N, D),  sub_key_list=[a1], val=ff)
                            db_change_list.extend([a1, a2])
                        data.append((ai, ff))
                    except Exception as err:
                        if kwargs.get('info', 0)>0: 
                            print(a1, a2,  err) 
                        if not fault_tol: 
                            raise 
            if db_change_list: 
                db_change_list = list(set(db_change_list))
                db_change_list.sort()
                print('fidelity updated for %s'%(db_change_list, ))
                for a in db_change_list: 
                    db = xx[a]
                    #print_vars(vars(),  ['db["fidelity"]'])
                    db.commit(info=0)
                self.set_rdb_dict(db_change_list, reload_db_class=0)                      
        else:
            for i in range(len(aa)):
                try:
                    a1=aa[i]
                    a2=aa[i+1]
                    
                    s1=xx[a1].load_S(sh)
                    s2=xx[a2].load_S(sh)
                    if not use_dist_comp: 
                        f=fidelity(s1, s2)    
                    else: 
                        pass 
                    #print_vars(vars(),  ['f'])
                    data.append((a1[param_name_id], f))
                except Exception as err:
                    pass
        return data

    def calc_fidelity_2d(self,aa, sh, which='fidelity_direct', eigs_tol=None, 
            save_to=None, force=0, fault_tolerant=1, info=0):
        if save_to is not None and not force: 
            res= self.load(save_to)
            return res 
        if self.algorithm == 'idrmg' : 
            from vmps.measure_and_analysis.measurement_idmrg_mcc import fidelity, fidelity_direct
        else: 
            raise 
        
        if which=='fidelity':
            fidelity_func=fidelity
        elif which=='fidelity_direct':
            fidelity_func=fidelity_direct
        else:
            raise 
        
        data=np.zeros((len(aa), len(aa)))
        faild_list = []
        db_changed = False
        for i, a1 in enumerate(aa):
            db1=self[a1]
            for j, a2 in enumerate(aa):
                db2=self[a2] 
                if a1 <= a2:
                    
                    #f12 = db1.fetch_easy_new(which, sh, sub_key_list=[a2], fault_tolerant=1)
                    f12 = db1.fetch_easy(which, sh, sub_key_list=[a2], fault_tolerant=1)
                    
                    #if f12 or f21: 
                    if f12 is not None and not force:  
                        ff = f12
                    else: 
                        try:
                            s1 = db1.load_S(sh);  s2=db2.load_S(sh)
                            if info>0: 
                                print('aaa', a1, a2)
                            ff = fidelity_func(s1, s2, eigs_tol=eigs_tol, info=info)
                            db1.insert(which, sh, iter=-1, sub_key_list=[a2], val=ff)
                            db2.insert(which, sh, iter=-1, sub_key_list=[a1], val=ff)
                            db2.commit()                
                            db_changed = True 
                        except Exception as err:
                        
                            if not fault_tolerant: 
                                print('error occured for', a1, a2)
                                raise 
                            else: 
                                faild_list.append((a1, a2, repr(err)))
                                ff=np.nan
                    data[i,j]=data[j,i]=ff
                
            if db_changed:                  
                db1.commit()   
            
        if db_changed: 
           self.set_rdb_dict(aa, reload_db_class=0)  
        
        if faild_list: 
            print('faild_list', faild_list) 
        res={'aa':aa,'data':data, 'faild_list': faild_list}      
        if save_to is not None : 
            from vmps import save
            save(res, save_to)
        return res     

    def submit_job(self, N=None):
        pass

class Analysis_mera(Analysis): 
    def __init__(self, **kwargs): 
        kwargs.update(algorithm='mera', result_db_class=ResultDB_mera)
        Analysis.__init__(self, **kwargs)

class Analysis_vmps(Analysis): 
     def __init__(self, **kwargs): 
        kwargs.update(algorithm='vmps', result_db_class=ResultDB_vmps)
        Analysis.__init__(self, **kwargs)

class Analysis_idmrg(Analysis): 
    
    def __init__(self, **kwargs): 
        """
        """        
        kwargs.update(algorithm='idmrg', result_db_class=ResultDB_idmrg)
        #sometimes in ipython, an error: unbound method __init__ ...; this is caused by relaod of ipython
        Analysis.__init__(self, **kwargs)
        #super(Analysis_idmrg, self).__init__(**kwargs)
    
    
    def calc_EE_del(self, S):
        #lam = S['Lambda']
        lam = S['lam']
        U, s, V = np.linalg.svd(lam)
        spect= s
        #print np.sum(s**2)
        spect2=spect**2
        EE = -np.sum(spect2*np.log(spect2))
        res = EE
        #import vmps 
        #vmps.print_vars(vars(), ['lam', 's', 'spect', 'spect2', 'EE'])
        return res 

    def get_EE(self, aa, D): 
        res=self.measure(self.calc_EE, aa, D=D)
        return res

    
    def update_state_files(self, aa, sh, updater): 
        """
            this should be do with caution 
        """
        from vmps import save
        for a in aa: 
            db = self[a]
            try: 
                state = db.load_S(sh)
            except IOError as err:
                print(err) 
                continue 
            new_state = updater(state)
            fn = db.__class__.shape_to_backup_fn(sh)
            dir = db.parpath 
            path = '/'.join([dir, fn])
            save(new_state, path)
            print('update file ... %s succeeded'%(path[-30: ]))
  
class Analysis_proj_qmc(Analysis): 
     def __init__(self, **kwargs): 
        kwargs.update(algorithm='proj_qmc', result_db_class=ResultDB_proj_qmc)
        Analysis.__init__(self, **kwargs)

class Analysis_tdvp(Analysis): 
     def __init__(self, **kwargs): 
        kwargs.update(algorithm='tdvp', result_db_class=ResultDB_tdvp)
        Analysis.__init__(self, **kwargs)


class TestAnalsysis(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_temp(self): 
        Alpha = namedtuple('Alpha', ['t', 'alpha', 'V', 'W'], 
                defaults = [None, None, None, None])
        
        an_bose_long_tdvp = Analysis_tdvp(
                local_root='/'.join([RESULTDB_ROOT, 'run-bose-long', 'tdvp', ]), 
                Alpha=Alpha, 
                param_list=Alpha._fields, 
                )
        an_bose_long_tdvp.add_sub_analysis()
        xx = an_bose_long_tdvp.an_disorder
        
        #xx.set_parpath_dict()
        #print(xx.alpha_list)
        #db = xx[0.5, 1.0, 1.5, 4.0, '0']
        #print(db)
        
        #xx.show_fig()
    
    def test_preprocess_alpha_list(self): 
        an = self.an
        aa = [0.5, 0.3, 1.0, ]
        print(list(an.alpha_parpath_dict.keys()))
        print(list(an.preprocess_alpha_list(aa).keys()))
    
    def test_alpha_list_operations(self): 
        an = self.an 
        print(an.result_db_class) 
        print(an.alpha_filter(param_list=['a'], a='a>1.0', sh=(3, 4)))
        aa = [1, 2, 3, 1.1, 1.01]
        res=an.alpha_filter(aa, param_list=['a'], resolution=0.1) 
        self.assertEqual(res, [1, 2, 3, 1.1])
        
    def test_plot_energy(self): 
        an = self.an
        an.plot_energy(an.alpha_list, (4, 4))
    
    def test_plot_magnetization(self): 
        pass 





if __name__ == '__main__':
    
    if len(sys.argv)>1: 
        parser = argparse.ArgumentParser(description='dont know')
        parser.add_argument('-dir', '--dir', default='./')
        #parser.add_argument('-w', '--which', nargs='*', default=None)
        parser.add_argument('-w', '--which',  default=None)
        parser.add_argument('-ew', '--exclude_which', nargs='*', default=None)
        parser.add_argument('-d', '--dim_list', nargs='*',  type=int,  default=None)
        parser.add_argument('-l', '--layer_list', nargs='*', type=int, default=None)
        parser.add_argument('-sh', '--mera_shape_list', nargs='*', default='all')
        parser.add_argument('-f', '--force', action='store_true', default=False) 
        parser.add_argument('-show', '--show', action='store_true', default=False) 
        args = parser.parse_args()
        args = vars(args)
        
        an = Analysis()
        aa = [1.0]
        an.alpha_parpath_dict[1.0] = args['dir']
        #an.set_rdb_dict()
        db = an.alpha_rdb_dict[1.0]
        #an.outline_db()
        print(args)
        which = args['which']
        if which in ['correlation_extra', 'correlation']: 
            sh = db.get_mera_shape_list()
            print(sh)
            an.plot_correlation(aa, sh, which=which, draw=1)
        #an = Analysis(path=None, parpath=args['dir'])

        if args['show']: 
            if which == 'scaling_dim' : 
                db.show_scaling_dim()
        #if 
   
    else: 
        
        if 0: 
            from .result_db import ResultDB 
            for a in aaa: 
                db = ResultDB(parpath=a)
                db.upgrade()
            
        if 0:
            unittest.main()
            
        else: 
            suite = unittest.TestSuite()
            add_list = [
               'test_temp', 
               #'test_preprocess_alpha_list', 
               #'test_plot_energy', 
               #'test_plot_magnetization', 
               #'test_alpha_list_operations', 
               
            ]
            for a in add_list: 
                suite.addTest(TestAnalsysis(a))
            unittest.TextTestRunner().run(suite)
           
        

