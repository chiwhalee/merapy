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
import unittest
import argparse, argcomplete
import os,  sys
import copy 
import pprint
import cPickle as pickle
import numpy as np
import pandas as pd
from collections import OrderedDict
import platform
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.font_manager import FontProperties
import inspect
import warnings
import importlib
import datetime
import itertools 

from merapy.hamiltonian import System
from tabulate import tabulate
from merapy.measure_and_analysis.measurement import mera_backup_dir_finder,  measure_all as measure_all_orig
from merapy.context_util import rpyc_conn 
from merapy.utilities import print_vars

from merapy.measure_and_analysis.result_db import (ResultDB, 
        ResultDB_mera, ResultDB_idmrg, ResultDB_vmps)
from result_db import (MATPLOTLIBRC, MARKER_LIST, MARKER_CYCLE, 
        AnalysisTools, 
        )

import matplotlib as mpl
mpl.rcParams.update(MATPLOTLIBRC)

__all__ = ['Analysis_mera', 'Analysis_vmps', 'Analysis_idmrg',
        'MARKER_LIST', 'MARKER_CYCLE', ]

class OrderedDictLazy(OrderedDict): 
    """
        a walk around for efficiency  
        lazy loading of each ResultDB, 
    """
    def __init__(self, alpha_parpath_dict, result_db_class=ResultDB, result_db_args={}): 
        OrderedDict.__init__(self)
        self.alpha_parpath_dict = alpha_parpath_dict
        self.result_db_class= result_db_class
        self.result_db_args= result_db_args
    
    def __getitem__(self, a): 
        if not self.has_key(a): 
            if self.alpha_parpath_dict.has_key(a): 
                p = self.alpha_parpath_dict[a]
                db =self.result_db_class(parpath=p, **self.result_db_args) 
                self[a] = db
            else: 
                msg = '%s not in self.alpha_parpath_dict'%(a, )
                warnings.warn(msg)
                #db = self.result_db_class(parpath='/tmp/empty.db.pickle', **self.result_db_args)
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
        temp = zip(*aa)
        # 这么做需要保证aa只有一个轴是不同的
        i = None   # some times aa is empty 
        for i in xrange(len(temp)): 
            if len(set(temp[i]))>1: 
                break 
        if i is not None : 
            param_name = self.param_list[i]
            param_name_id = i 
        else: 
            param_name = None 
            param_name_id = None 
        return param_name, param_name_id 

class Analysis(AnalysisTools):
    """
    """
    REMOTE_HOST = 'zhihuali@211.86.151.102:'
    
    #DISTANCE_LIST = [3**i for i in range(1, 9)]
    def __init__(self, fn=None, path=None, parpath=None, backup_tensor_root='./', 
            remote_root = None, local_root=None, param_list=None, 
            result_db_class=None, result_db_args= {}, algorithm=None):
        self.remote_root = remote_root
        self.local_root = local_root
        self.param_list= param_list
        self.result_db_class= result_db_class if result_db_class is not None else ResultDB
        self.result_db_args= result_db_args if result_db_args is not None else ResultDB
        self.algorithm = algorithm
        
        if 0:  #remove these later
            if fn is not None:
                path = os.path.abspath(fn); path = fn[-1::-1];  n = path.find("/");   
                path = path[n:]; path = path[-1::-1]
                self.fn_only = path[:n]
            
            elif path is not None:
                self.path = path
            else:
                pass
            self.parpath = parpath
        
        self.alpha_parpath_dict = OrderedDict()   #this attr should always be updated, while not replaced, or else, it affects the next line
        self.alpha_rdb_dict = OrderedDictLazy( self.alpha_parpath_dict, 
                result_db_class=self.result_db_class, result_db_args= self.result_db_args)
        self.rdb_dict = self.alpha_rdb_dict
        
        self.sh_list_max = None 
        self.sub_analysis= []
        if 1: #plotting configs
            self.invert_xaxis= 0
            self.invert_alpha_order = 0
    
    def __getitem__(self, k): 
        
        #if isinstance(k, float): 
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

    def get_run_failure(self, sh, aa_all=None): 
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

    def alpha_filter_bac(self, a=None, h=None, V=None):  #, **kwargs): 
        """
            for wigner_crystal 
        """
        aa = list(self.alpha_list)
        if h == 1.4: 
            h = 1.4142135
        param_list = ['a', 'h', 'V']
        #print vars()
        for i, x in enumerate(param_list): 
            #i = param_list.index(x)
            x_val=vars()[x]
            if x_val is not None : 
                aa = filter(lambda x: x[i]== x_val, aa)
        return aa
    
    def filter_alpha(self, aa=None, sh=None, field=None, param_list=None, 
            surfix='', resolution=None, **kwargs): 
        """
            from wigner_crystal 
        """
        aa = aa if aa is not None else list(self.alpha_list) 
        param_list = param_list if param_list is not None else self.param_list #['a', 'h', 'V']
        
        num_param = len(param_list)
        #make each a into a tuple  if it is not 
        aa = map(lambda a:  a if hasattr(a, '__iter__') else (a, ), aa)
        #if no_surfix: 
        #    #exclude such as (1.0, 'sufix') from aa 
        #    aa = filter(lambda x: not isinstance(x[-1], str),  aa)
        
        if surfix is not None : 
            #exclude such as (1.0, 'sufix') from aa 
            if surfix == '' : 
                aa = filter(lambda x: not isinstance(x[-1], str),  aa)
            else:
                aa = filter(lambda x: x[-1]==surfix,  aa)
        
        for k, v in kwargs.items(): 
            if v is None:  # None means no constraint 
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
            aa = filter(func, aa)
        
        if resolution is not None:
            resolution = float(resolution)
            #use a[: num_param] is because there may be surfix a = (1.0, 'a') then 'a' is excluded
            func=lambda a:   all(map(lambda x:  (x/resolution).is_integer(), a[: num_param]))
            aa = filter(func, aa)
            
        aa = map(lambda a:  a if len(a)>1 else a[0], aa)
        if sh is not None : 
            def func(a): 
                dir = self.alpha_parpath_dict[a]
                fn = self.result_db_class.shape_to_backup_fn(sh)
                path  = '/'.join([dir, fn])
                #print 'ppp', path 
                return os.path.exists(path)
            aa = filter(func, aa)
        if field is not None: 
            if sh is None: 
                func = lambda a: self[a].has_key(field)
            else: 
                func = lambda a: self[a].fetch_easy(field, sh)
            aa = filter(func, aa) 
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
    
    def scan_backup_dir_1_del(self, root='./'): 
        nlist = os.listdir(root)
        alpha_parpath_dic = {} 
        root = '/home/zhli/Documents/mera_backup_tensor/run-wigner-crystal/'    
        for name in nlist: 
            if 'a=' in name: 
                aaa = name.split('-')
                bbb = [a.split('=')[1] for a in aaa]
                alpha = [float(b) for b in bbb]
                #alpha_list.append(tuple(alpha))
                alpha_parpath_dic[tuple(alpha)]=root + '/%s/'%(name, )
       
        return alpha_parpath_dic
    
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
            print 'dir_list is', dir_list
        signiture = signiture if signiture is not None else '='
        
        if 1: 
            def parse(xx): 
                if '=' in xx: 
                    xx = xx.split('=')[1]
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
                
            alpha_parpath_dict = OrderedDict()
            for name in dir_list: 
                if signiture in name: 
                    
                    split = name.split('-')
                    if len(split)>1: #more than one param 
                        alpha = tuple([parse(a) for a in split])
                    else: 
                        #alpha = str(parse(split[0]))
                        alpha = parse(split[0])
                    
                    if surfix != '': 
                        alpha = str(alpha) + '-' +  surfix
                    alpha_parpath_dict[alpha] =  '/'.join([root, name])
            
            res= alpha_parpath_dict
        
        return res
    
    def scan_alpha_1_del(self, root='./'): 
        """
            support for multi param
        """
        nlist = os.listdir(root)
        alpha_list = []
        for name in nlist: 
            if 'a=' in name: 
                aaa = name.split('-')
                bbb = [a.split('=')[1] for a in aaa]
                print bbb
                alpha = [float(b) for b in bbb]
                alpha_list.append(tuple(alpha))
        alpha_list.sort()
        return alpha_list
    
    def scan_sub_analysis_top(self): 
        temp = os.listdir(self.local_root)
        temp = filter(lambda x: '=' not in x, temp)
        temp = filter(lambda x: os.path.isdir('/'.join([self.local_root, x])), temp)
        #print_vars(vars(),  ['temp']) 
        return temp 
    
    def add_sub_analysis(self,  hook_list=None): 
        if hook_list is None: 
            hook_list = []
        if hook_list == 'auto': 
            hook_list = self.scan_sub_analysis_top()
        for h in hook_list: 
            sub_dir = h if isinstance(h, str) else h['sub_dir']   #h may be a dir name or in the old format {'sub_dir': dir_name_str}
            sub_root = '/'.join([self.local_root, sub_dir])
            name = 'an_%s'%(sub_dir.replace('/', '_'), )
            if name in self.sub_analysis: 
                continue 
            antemp = self.__class__(local_root=sub_root, param_list=self.param_list)
            antemp.set_parpath_dict()
            setattr(self, name, antemp)
            self.sub_analysis.append(name)
    
    def set_rdb_dict(self, alpha_list=None, create_empty_db=1, 
            db_version=None, algorithm=None, info=0):        
        if 1: 
            module_name = self.result_db_class.__module__
            module = importlib.import_module(module_name)
            reload(module)
            self.result_db_class= module.__getattribute__(self.result_db_class.__name__)
            if info>0: 
                print '%s is reloaded'%(self.result_db_class.__name__, )
        
        part_update = 0
        if alpha_list is not None: 
            part_update = 1
            
        alpha_list = alpha_list if alpha_list is not None else  sorted(self.alpha_parpath_dict)
        error = []
        for a in alpha_list: 
            try: 
                p = self.alpha_parpath_dict[a]
                db = self.result_db_class(parpath=p, create_empty_db=create_empty_db, version=db_version, algorithm=algorithm)
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
            print msg
        
        if len(error)>1: 
            print error
    
    def set_parpath_dict(self, root=None, signiture=None, info=0): 
        if root is None: 
            root = self.local_root
        assert root is not None , 'param root reqired'
        dic = self.scan_alpha(root=root, signiture=signiture)
        self.alpha_parpath_dict.update(dic)
        if info>0: 
            msg = 'parpath found are %s'%(dic, )
            print(msg)
        
        self.alpha_list_all = dic.keys()
        self.alpha_list_all.sort()
        self.alpha_list = list(self.alpha_list_all)
    
    def search_alpha(self, alpha_str, alpha_list=None):
        if not isinstance(alpha_str, str): 
            alpha_str = str(alpha_str)
        if alpha_list is None: 
            all_key = self.alpha_parpath_dict.keys()
        else: 
            all_key = alpha_list
        
        res= []
        for a in all_key: 
            if alpha_str in str(a): 
                res.append(a)
        return res
    
    def shape_to_backup_fn(self, sh): 
        return System.shape_to_backup_fn(*sh)
    
    def search_field(self, field_name): 
        res = []
        for a in self.alpha_parpath_dict.keys(): 
            if self[a].has_key(field_name): 
                res.append(a)
        return res
            
    def upgrade_result_db(self): 
        res = mera_backup_dir_finder(root)
    
    def measure_all(self, aa=None, sh_list='all', sh_min=None, sh_max=None,  which=None, 
            exclude_which=None, force=0, fault_tolerant=1, recursive=False,  **kwargs): 
        if aa is None: 
            aa = self.alpha_list
        aa = list(aa)
        dir_list = [self.alpha_parpath_dict[a] for a in aa]
        measure_all_orig(dir_list, mera_shape_list=sh_list, sh_min=sh_min, sh_max=sh_max, which=which, 
                exclude_which = exclude_which, fault_tolerant = fault_tolerant, force=force, **kwargs)
        self.set_rdb_dict(aa)
    
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
            print 'outline %s at %s'%(dir_name, str(datetime.datetime.now()), )
        data=[]
        if mera_shape_list is None: 
            mera_shape_list=[(4, 4)] + [(8, i) for i in [4,5]] + [(12, i) for i in [4,5]]
        if alpha_list is None: 
            alpha_list = self.alpha_list
        
        res = []
        for a in alpha_list:
            if which == 'file':  
                sh_list = self.alpha_rdb_dict[a].get_mera_shape_list()
                
                diff_list = []
                if show_eng_diff: 
                    for sh in sh_list: 
                        dir = self.alpha_parpath_dict[a]
                        fn = self.shape_to_backup_fn((sh[0], sh[1]-1))
                        path = '/'.join([dir, fn]) 
                        S=self.load(path)
                        diff = S.energy_diff
                        diff_list.append('%1.2e'%diff)
                        
                    sh_list = zip(sh_list, diff_list)
            
            elif which == 'db': 
                db=self.alpha_rdb_dict[a]
                if db.get('version')==1.0:
                   
                    keys= db.get('energy')
                    if keys is None: 
                        pass
                    keys= keys.keys()
                    keys.sort()
                    sh_list = keys
                    
                else: 
                    keys= db.keys()
                    keys= [k[: 2] for k in keys]
                    keys= list(set(keys))
                    keys.sort()
                    sh_list = keys
                    
            if not grid_view: 
                print a
                print '\t', sh_list, '\n'
                
            d= [str(a)]
            for sh in mera_shape_list:
                if sh in sh_list:
                    d.append('   *   ')
                else:
                    d.append('       ')
            data.append(d)
        
        if grid_view:     
            print tabulate(data, headers= ['a']+map(str,mera_shape_list), tablefmt='grid', stralign='center') 
    
    def outline_db(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0):
        self._outline('db', alpha_list=alpha_list, mera_shape_list=mera_shape_list, info=info)
    
    def outline_file(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0): 
        self._outline('file', alpha_list=alpha_list, mera_shape_list=mera_shape_list, show_eng_diff=show_eng_diff, info=info)
    
    def outline(self, alpha_list=None, indent=0, recursive=1): 
        """
            outline_db, outline_file will be deprecated in future 
        """
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
            temp = map(lambda x: x[-1] if isinstance(x[-1], str) else '', alpha_list)
            surfix_list = list(set(temp))
            if '' in surfix_list: 
                surfix_list.remove('')
            surfix_list.sort()
            msg += ''.join([tab, 'surfix_list: \n']) 
            msg += ''.join([tab, '\t', str(surfix_list), '\n'])
            self.surfix_list = surfix_list
            
                
        print msg 
        if recursive: 
            for a in self.sub_analysis: 
                an = getattr(self, a)
                print ''.join([tab, a])
                an.outline(recursive=1, indent=indent+1)
    
    def preprocess_alpha_list(self, alpha_list, alpha_remap=None): 
        parser = lambda a: a if not isinstance(a, str) else  float(a.split('-')[0]) 
        if 1: 
            aa = OrderedDict()
            aa_bad = []
            for a in sorted(alpha_list):
                if self.alpha_parpath_dict.has_key(a):
                    aa[parser(a)] = self.alpha_rdb_dict[a]
                else: 
                    aa_bad.append(a)
        
        if len(aa_bad)>0: 
            print 'keys %s are not in self.alpha_parpath_dict'%aa_bad
            #warnings.warn('keys %s are not in self.alpha_parpath_dict'%aa_bad)
        return aa
            
    def fetch_many(self, field_name, aa, sh): 
        aa_dic = self.preprocess_alpha_list(aa)
        res= OrderedDict()
        for a in aa_dic:
            db = aa_dic[a]
            rec=db.fetch_easy(field_name, sh)
            res[a] = rec
        return res
             
    def merge_remote_db(self, aa, remote_root=None): 
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
        #remote_root = 'zhihuali@211.86.151.102:mera_backup_tensor/run-long-better/'
        #locale_root = '/home/zhli/Documents/mera_backup_tensor/run-long-better/'
        
        #A_root = 'zhihuali@211.86.151.102:' + remote_root
        alpha_list = list(alpha_list)
        A_root = self.REMOTE_HOST + self.remote_root
        #B_root = '/home/zhli/Documents/mera_backup_tensor/run-long-better/'
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
                print 'local dir not found, mkdir %s'%B_dir
            
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
                print msg
            
           
            #if measure: 
            if 0: 
                path = B_path
                mm = Measurement(path=path)
                S= System.load(path)
                mm.measure_all(S, path, which=None, show=False, exclude_which=['scaling_dim'])
    
    def get_shape_list_all(self, N=None, D=None, aa=None, sh_min=None, sh_max=None):
        if self.sh_list_max is not None : 
            sh_list_max = list(self.sh_list_max)
        else: 
            sh_list_max = set()
            aa = aa if aa is not None else list(self.alpha_list)
            for a in aa: 
                db = self[a]
                sh_list = db.get_shape_list()
                sh_list_max = sh_list_max.union(sh_list)
            sh_list_max = list(sh_list_max)
        if N is not None:
            sh_min = (N, 0)
            sh_max = (N, np.inf)
        if D is not None: 
            sh_min = (sh_min[0], D)
            sh_max = (sh_max[0], D)
        if sh_max is not None : 
            sh_list_max=filter(lambda x: x[0]<=sh_max[0] and x[1] <= sh_max[1] , sh_list_max)
        if sh_min is not None : 
            sh_list_max=filter(lambda x: x[0]>= sh_min[0] and x[1]>= sh_min[1] , sh_list_max)
        sh_list_max.sort()
                
        return sh_list_max  
        
    def get_fn_list_bac(self):
        nlist = os.listdir(self.parpath)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" in i]
        print "all pickle files under specified path are %s"%pickle_files
        
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
            print i, parse_fn(i)
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
                print a
                self[a].delete_db()
            except Exception as err: 
                print err 
    
    def delete_attr_in_file(self, name, aa,  sh): 
        #aa=self.filter_alpha(sh=(1000, 20))
        for a in aa[2:]:
            db=self[a]
            dir=db.parpath
            try:
                name= db.backup_fn_gen(sh)
                path='/'.join([dir,name])
                print path
                s=vmps.load(path)
                s.pop(name)
                vmps.save(s,path)
            except Exception as err:
                print err
    
    def get_fn_list(self):
        nlist = os.listdir(self.parpath)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" == i.split('.')[-1]]
        print "all pickle files under specified path are %s"%pickle_files
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
                print i, path
                if 0:
                    inn = open(path, 'rb')
                    S= pickle.load(inn)
                    inn.close()
                S= System.load(path)  #System.load is not exactly pickle.load
                res[i]['energy'] = S.energy
                res[i]['energy_record'] = S.energy_record
            #print res
            out = open(self.parpath +'/' +  'record.pickle.db', 'wb')
            print 'file record.pickle.db is generated'
            pickle.dump(res, out)
            out.close()
    
    def dim_vs_energy_err(self):
        pass

    def find_max_dim_file_bac(self, path):
        nlist = os.listdir(path)
        
        pickle_files = [i for i in nlist if "pickle" in i]
        print "all pickle files under specified path are %s"%pickle_files
        
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
        print 'head is %s'%res
        return res
    
    def monitor(self): 
        """
            moniting the excution state of jobs
            what to watch: 
                eng, eng_diff, time spend, etc
        """
        pass
        
    def _plot(self, x, y, **kwargs):
        fig=ResultDB._plot.im_func(None, x, y, **kwargs)  
        return fig
    
    def _plot3d(self, xx, yy, data, which_plot=None, zfunc=None,   **kwargs):
        which_plot = 'imshow' if which_plot is None else which_plot 
        figsize = kwargs.get('figsize')
        figsize = (4, 4) if figsize is None else figsize
        fig = kwargs.get('fig', None)
        ax = kwargs.get('ax', None)
        info = kwargs.get('info', 0)
        #
            #if ax is None:
            #    fig=plt.figure(figsize=(5,5))
            #    ax=fig.add_subplot(111,)
            #    #ax=fig.add_subplot(111, projection='3d')
            #else:
            #    fig=ax.figure
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
            #if kwargs.get(t):  
            if kwargs.has_key(t): 
                style_dic[t]=kwargs[t]
     
        not_found = []
        data = data.T
       
        if kwargs.get('make_up'):  
            make_up = kwargs['make_up']
            make_up_list = []
            m, n = data.shape 
            #print np.where(np.isnan(data))
            aaa, bbb = np.where(np.isnan(data))
            ab = zip(aaa.tolist(), bbb.tolist())
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
        cb = None 
        if which_plot == 'contourf': 
            #cb = ax.contourf(Y, X, Z, zdir='z', offset=0, cmap=plt.cm.coolwarm, alpha=0.9, aspect=1 )
            cb = ax.contourf(X, Y, Z, zdir='z', offset=0, cmap=plt.cm.coolwarm, aspect=1 )
            #cb = ax.imshow(X, Y, Z, zdir='z', offset=0, cmap=plt.cm.coolwarm, aspect=1 )
            fig.colorbar(cb, ax=ax)
            clim = kwargs.get('clim')
            if clim is not None : 
                cb.set_clim(clim)
        elif which_plot == 'imshow' : 
            #extent = [min(xx)-0.05, max(xx) + 0.05, min(yy), max(yy)]
            if style_dic.get('extent') is None:
                dx = (max(xx)-min(xx))*1.0/(len(xx)-1)
                dy = (max(yy)-min(yy))*1.0/(len(yy)-1)
                #print 'dddd', dx, dy
                style_dic['extent'] = [min(xx), max(xx)+dx, min(yy), max(yy)+dy]
            im = ax.imshow(Z,  #extent=extent, 
                    origin='lower', 
                     **style_dic) 
                    #aspect=len(xx)*0.5/len(yy), **style_dic) ,aspect=1 )
            cb=fig.colorbar(im, ax=ax)
            clim = kwargs.get('clim')
            if clim is not None : 
                cb.set_clim(clim)
        elif which_plot == 'surface': 
            temp = ['interpolation']
            for t in temp: 
                style_dic.pop(t)
            ax.plot_surface(XX, YY, data, rstride=1, cstride=1, **style_dic)
        else: 
            raise ValueError("alowed which_plot=%s"%(['contourf', 'imshow', 'surface']))
       
        if kwargs.has_key('title') and ax.is_first_row(): 
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
        
    def plot_field_vs_alpha(self, field_name, aa=None, sh_list=None,
            sub_key_list=None, alpha_db_map=None, rec_getter=None,  rec_getter_args=None, 
            param_name=None,  **kwargs): 
        
        aa = aa if aa is not None else list(self.alpha_list)
        if isinstance(sh_list, tuple):  
            sh_list = [sh_list]
        sh_min = kwargs.get('sh_min')
        sh_max = kwargs.get('sh_max')
        sh_list = sh_list if sh_list is not None else self.get_shape_list_all(aa=aa, sh_min=sh_min, sh_max=sh_max)
        sh_list = reversed(sorted(sh_list))
        #sub_key_list = sub_key_list if sub_key_list is not None else []
        alpha_db_map = alpha_db_map if alpha_db_map is not None else {}
        info = kwargs.get('info', 0)
        algorithm = self.algorithm 
        if algorithm == 'vmps':
            if field_name == 'entanglement': 
                field_name += '_entropy' 
        
        if param_name is None: 
            if 0: 
                # 这么做需要保证aa只有一个轴是不同的
                temp = zip(*aa)
                i = None   # some times aa is empty 
                for i in xrange(len(temp)): 
                    if len(set(temp[i]))>1: 
                        break 
                if i is not None : 
                    param_name = self.param_list[i]
                    param_name_id = i 
                else: 
                    param_name = None 
                    param_name_id = None 
            else: 
                param_name, param_id = AlphaList.identify_param.im_func(self, aa)
        else: 
            param_name_id = self.param_list.index(param_name)
            
        kwargs_orig = kwargs.copy()
        not_found=[]
        for sh in sh_list: 
            #print 'sss', sh 
            data=[]
            for a in aa:
                if alpha_db_map.has_key(a): 
                    db = alpha_db_map[a]
                else: 
                    db=self[a]
                if rec_getter is None: 
                    rec=db.fetch_easy(field_name, sh, sub_key_list=sub_key_list)
                else: 
                    rec_getter_args= rec_getter_args if rec_getter_args is not None else {}
                    rec = rec_getter(db, sh, **rec_getter_args)
                
                #if param_name is not None : 
                #    a1 = a[self.param_list.index(param_name)]
                a1 = a[param_name_id]
                if rec is not None:
                    data.append((a1, rec))
                else:
                    data.append((a1, np.nan))
                    not_found.append(a)
                    
            if not_found and info>0: 
                print 'not_found is ', not_found
            
            try: 
                x,y=zip(*data)        
            except ValueError as err: 
                print err
                x, y = np.nan, np.nan 
            x = np.asarray(x); y=np.asarray(y)
            label = kwargs_orig.get('label', str(sh))
            marker = kwargs_orig.get('marker', MARKER_CYCLE.next())
            kwargs.update(ylabel=field_name, return_ax=1, label=label, 
                    maker=marker, xlabel='$%s$'%param_name)
            fig, ax = self._plot(x, y, **kwargs) 
            
            kwargs['fig'] = fig 
            kwargs['ax'] = ax
            
        if kwargs.get('return_fig'): 
            return fig 
    
    def plot_field_vs_alpha_3d(self, field_name, aa, sh, xparam, yparam, 
           make_up=False, sub_key_list=None, **kwargs): 
        info = kwargs.get('info', 0)
        aa = aa if aa is not None else self.alpha_list
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
            grid1 = map(lambda a:(a[p1], a[p2]),   aa)
            temp=set(grid)-set(grid1)
            temp = list(temp)
            temp.sort()
            if temp: 
                print 'points not in alpha_parpath_dict is', temp
        
        for a in aa: 
            x = a[p1]; y = a[p2]
            i = x_to_i[x]; j = y_to_j[y]
            db = self[a]
            rec = db.fetch_easy(field_name, sh, sub_key_list, info=0)
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
            ab = zip(aaa.tolist(), bbb.tolist())
            for a, b in ab: 
                if a>0 and a <m: 
                    data[a, b] = (data[a-1, b]  +  data[a+1, b])/2. 
                    make_up_list.append((a, b))
                elif a == 0: 
                    data[a, b] = data[a+1, b]
                elif a == m-1: 
                    data[a, b] = data[a-1, b]
        kwargs['make_up'] = make_up 
        if not_found and info>0: 
            print  'not_found is ', not_found 
        kwargs.update(xlabel=xparam, ylabel=yparam)
        dic = self._plot3d(xx, yy, data, **kwargs)
        fig, ax, cb = dic['fig'], dic['ax'], dic['cb']
        #ax.set_aspect(1)
        
        if kwargs.get('return_fig'): 
            return fig 
    
    def get_rec_2d(self, field_name, aa, sh, xparam, yparam, 
            sub_key_list=None, rec_getter=None,  **kwargs): 
        info = kwargs.get('info', 0)
        aa = aa if aa is not None else self.alpha_list
        xx = self.get_param_range(xparam, aa=aa)
        yy = self.get_param_range(yparam, aa=aa)
         
        XX,YY=np.meshgrid(xx,yy)
        data=np.ndarray((len(xx),len(yy)))
        #data[:, :] = np.nan 
        alg = self.algorithm
        
        sub_key = ['EE'] if self.algorithm == 'vmps' else None 
        #sub_key = None 
        not_found = []
         
        p1=self.param_list.index(xparam) 
        p2=self.param_list.index(yparam) 
        x_to_i = {x:xx.index(x) for x in xx}
        y_to_j = {y:yy.index(y) for y in yy}
        for a in aa: 
            x = a[p1]; y = a[p2]
            i = x_to_i[x]; j = y_to_j[y]
            db = self[a]
            if rec_getter is None: 
                rec = db.fetch_easy(field_name, sh, sub_key_list, info=0)
            else: 
                rec = rec_getter(db, sh)
            if rec is not None : 
                data[i, j] = rec 
            else:
                data[i, j] = np.nan 
                not_found.append((x, y))
                
        return xx, yy, data  
        #print XX.shape, xx.shape, YY.shape, yy.shape, data.shape
    
    def show_fig(self): 
        plt.show()
    
    def fig_layout(self, ncol=1, nrow=1,  size=(3.5, 2.5), dim=2): 
        return ResultDB.fig_layout.im_func(None, ncol, nrow, size, dim=dim)
   
    def _scaling_dim_exact(self):
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

        ee = self.obj.energy_exact
        mapper = {-1.2732:"xx", -1.6077:"heisbg_NNN", -1.6449:"heisbg_long", -1.7725:'heisbg'}
        for i in mapper:
            if abs(i-ee)<1e-2:
                model = mapper[i]
        exact_val = exact[model]
        if self.obj.model == "Ising": 
            exact_val = exact["Ising"]

        return exact_val
    
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
            db_meth.im_func.__call__(db, **args)
       
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
            db_meth.im_func.__call__(db, **args)
        
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
            print exact_val

        x = rec.keys()
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
                qns = range(3)
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
                        print "index qn=%d, i=%d not found in exact_val, omit it "%(qn, i)
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
        
        print "ssss\n", rec[x[-1]]
        if hasattr(self.obj, "log"):
            print "llll\n", self.obj.log
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
            for a , db in aa_dic.iteritems(): 
                
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
                print 'not_found is', not_found
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
        x = recs.keys()
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
            t_eng = [(t, v[0])  for t, v in recs.iteritems()]
            #Et_diff = [(i[0], i[]) for i in t_eng]

            t, eng = zip(*t_eng)
            t_diff = [(t_eng[i+1][0], (t_eng[i+1][1]-t_eng[i][1])/(t_eng[i+1][0]-t_eng[i][0])) 
                    for i in range(len(t_eng)-1)]
            t, diff = zip(*t_diff)
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
        
    def get_correlation(self, dir='xx', alpha=None, mera_shape=None, info=0): 
        """
        mera_shape: a tuple (dim, layer)
        
        """
        #if dir  == 'xx': dir = 'pm'
        
        db = self.alpha_rdb_dict[alpha]
        dim = mera_shape[0]
        layer = mera_shape[1]
        record = db.fetch(field_name='correlation', dim=dim, num_of_layer=layer)
        items = record[dir]['val']
        x, y = zip(*items)
        x = np.array(x)
        y = np.array(y) 
        if info>0: 
            print 'X=', x
        if dir=='pm': 
            y=y*2 
        return y

    def get_energy_diff(self):
        """
        one record is like this:
        {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        
        """
        recs= self.obj.energy_record
        x = recs.keys()
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
        if info>0: print 'X=', x

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
            print 'rec not found for a=%s'%a
            return
        x,y=zip(*rec[0: -1: 2])        
        y = np.array(y) #y=abs(y)
        if direct == 'pm':  y *= 2 
        ax.loglog(x, y, '.-', label=str(a))
    
    def plot_central_charge(self, aa, sh, alpha_remap={}, alpha_sh_map={}, 
            layer='top', **kwargs):
        
        aa.sort(reverse=1)
        if layer=='top':
            ll= [sh[1]-2 ] 
        elif layer=='all':
            ll=  range(sh[1]-1)
        else: 
            ll = layer
        
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
                    print  db.path
                    raise
                c= (e2-e1)*3
                cc.append((a_orig,c))
            if len(cc)>0:     
                x,y=zip(*cc)    
                label = str(layer) if not kwargs_orig.has_key('label') else kwargs_orig['label']
                marker = '.' if not kwargs_orig.has_key('marker') else kwargs_orig['marker']
                kwargs.update(marker=marker, label=label, return_ax=1)
                fig, ax=self._plot(x, y, **kwargs)
                #fig, ax=self._plot(x, y, **kwargs)
                kwargs['ax'] = ax
                kwargs['fig'] = fig
        if 'ax' in locals().keys(): 
            if ax.is_first_col(): 
                ax.set_ylabel('central change')
            if ax.is_last_row(): 
                ax.set_xlabel('$\\alpha$')
            ax.invert_xaxis()
        if len(not_found)>0: 
            print 'rec for central_charge not_found is :%s'%not_found
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
                print  msg, 
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
                    x,y=zip(*rec[0:-1:2])
                else: 
                    x, y = zip(*rec[: ])
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
            print 'not_found is ', not_found
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
            if not kwargs.has_key('marker'):
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
            print 'not_found is ', not_found
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
            fault_tolerant = kwargs.get('fault_tolerant', True)
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
               
                
                shape = sh_remap[a] if sh_remap.has_key(a) else sh
                e=db.fetch_easy(name, shape, sub_key_list=sub_key_list_1, 
                        info=info, fault_tolerant=fault_tolerant)
                
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
            print 'not_found is %s'%(not_found, )
        
        if self.invert_xaxis: 
            ax.invert_xaxis()
        if ax.is_first_col(): 
            ax.set_ylabel('entanglement',fontsize=18)
        if ax.is_last_row(): 
            xlable = '$\\alpha$'
            ax.set_xlabel('$\\alpha$')
            if param_name is not None : 
                xlabel = '$%s$'%param_name 
            if kwargs.get('xlabel' ): 
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
            data = [i[1] for i in rec.values()]
            ll=rec.keys()
                
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
            label = args['label'] if kwargs.has_key('label') else str(a) 
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
            label = args['label'] if kwargs.has_key('label') else str(a)
            args.update(return_ax=1, label=label)
            fig, ax =db.plot_entanglement_vs_chi(  **args)
            if ax is not None :  #plot succecss
                args['ax'] = ax
            if fig is not None: 
                args['fig'] = fig
        
        if kwargs.get('return_fig'): 
            return fig 
   
    def plot_entanglement_spectrum(self, aa, sh, which_level=None, num_level=None, param_name=None, **kwargs): 
        num_level= 20  if num_level is None else num_level
        #sh=(0,80)
        #ax=fig.axes[i]
        if which_level is not None : 
            num_level = len(which_level)
        ES=np.ndarray((len(aa), num_level))
        for j,a in enumerate(aa):
            db=self[a]
            es=db.fetch_easy('entanglement_spectrum', sh)
            if es is not None:
                if which_level is not None : 
                    ES[j]=es[which_level]
                else: 
                    ES[j,:]=es[:num_level]
            else:
                ES[j,:]=np.nan
        if param_name is not None : 
            ii=self.param_list.index(param_name)
            x=[_[ii] for _ in aa]
        else: 
            x = aa 
        #ax.plot(x, -np.log10(ES), '-o', )
        self._plot(x, -np.log10(ES), **kwargs)
        
        

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
            print  'not_found is ', not_found
        
        fig = self._plot(x, y, marker='.', **kwargs)
        if kwargs.get('return_fig') : 
            return fig
            
              
    def show_scaling_dim(self, aa, sh, alpha_remap={}, alpha_shape_dic={}): 
        data=[]
        aa = list(aa)
        aa.sort(reverse=1)
        sh_orig = sh
        for a in aa:
            a_orig = a
            
            temp = alpha_remap.get(a)
            sh_temp = alpha_shape_dic.get(a)
            
            if sh_temp is not None :  
                sh = sh_temp
            else: 
                sh = sh_orig

            if temp is not None: 
                a = temp
            db=self.alpha_rdb_dict[a]
            rec= db.fetch_easy('scaling_dim', sh)
            sdim= rec
            #print a, '\t', 0, sdim[0][:3], 1, sdim[1][:3], 2, sdim[2][0:1]
            #print a, '\t', 0, sdim[0][:3], 1, sdim[1][:3], 2, sdim[2][0:1]

            if rec is not None:
                #data.append((a,)+sdim[0][:3]+ sdim[1][:3]+ sdim[2][:1])
                data.append((a,) + rec[0][:5] + rec[1][:3] +rec[2][:3] )
        #print tabulate(data, headers=['']+[str(i) for i in [0,0,0,1,1,1,2]], tablefmt='simple') 
        print tabulate(data, headers=['']+[str(i) for i in [0,0,0,0,0,1,1,1,2]], tablefmt='simple')

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

    def _plot_eng_fluc_count(self, aa, sh):
        """
            not for general use, only for trouble shoting,  so add an under core
        """
        bins=[-1, 0,  1e-8, 1e-7, 1e-6, 1e-5, 1]
        aa=list(aa)
        temp={}
        for a in aa:
            db=self.alpha_rdb_dict[a]
            #diff=db.get_energy( sh=sh, which='diff', iter_max=None)
            #diff_cut=np.histogram(diff, bins=bins)[0]
            #print a , diff_cut[1:]
            #temp[a]=diff_cut[1:]
            count, bin = db.get_energy_fluc_count(sh=sh)
            temp[a] = count[1: ]
        df=pd.DataFrame(temp, index=['<1e-8', '<e-7', '<1e-6', '<1e5', '<1'])
        #df=pd.DataFrame(temp, index=[str(i) for i in bins[1: ]])
        #df.unstack()
        #df.keys()
        df_trans = df.transpose()
        df_trans.plot(figsize=(6,3))
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
        x, y = zip(*eta)   
        self._plot(x, y, marker='.', **kwargs)
        #ax.plot(x, y, '.-')  
    
    #del this
    def _fit_correlation_plot(self, func, aa, sh, func_type='power'): 
        for a in aa:
            fig=plt.figure(figsize=(12,4))
            ax=fig.add_subplot(1,2,1)
            ax2=fig.add_subplot(1,2,2)
            data=[]
           
            db=self.alpha_rdb_dict[a]
            rec=db.fetch_easy('correlation_extra', mera_shape=sh, sub_key_list=['pm'])
            if rec is None: 
                continue
                
            if func_type == 'power': 
                x,y=zip(*rec[0:-1:2]); 
                x=np.array(x); y=np.array(y);  y=y*2;  y=np.abs(y)
            elif func_type == 'exp': 
                n = 200
                x,y=zip(*rec[20:n:2])  #len(rec)=77
                x=np.array(x); y=np.array(y);  y=y*2;  y=np.abs(y)
       
            param, cov = curve_fit(func, x, y)
            const, k = param[:2]
            data.append( (a,) + tuple(param) + tuple(cov.diagonal()))
            label = str(a)
            label += '$ \\eta=$%1.2f'%k 
            #print a,  param
            
            if 1:
                ax.loglog(x, y, '--', label=label)
                ax.loglog(x, func(x, *param))
            if 1:
                ax2.plot(x,np.log10(y), '--', label=label)
                y=func(x, *param)
                ax2.plot(x, np.log10(y))
            if 0:
                ax.plot(x, y, '--', label=label)
                y=func(x, *param)
                ax.plot(x, y)
                        
            ax.grid(1)
            ax.legend()
        
    #del this
    def _fit_correlation_noplot(self, func, aa, sh, func_type='power'):
        data=[]
        error = []
        for a in aa:
            db=self.alpha_rdb_dict[a]
            rec=db.fetch_easy('correlation_extra', mera_shape=sh, sub_key_list=['pm'])
            if rec is None: 
                continue
            
            if func_type == 'power': 
                x,y=zip(*rec[0:-1:2]); 
                x=np.array(x); y=np.array(y);  y=y*2;  y=np.abs(y)
            elif func_type == 'exp': 
                n = 200
                x,y=zip(*rec[20:n:2])  #len(rec)=77
                x=np.array(x); y=np.array(y);  y=y*2;  y=np.abs(y)
            try:
                param,cov = curve_fit(func, x, y)
                if isinstance(cov, float): 
                    #cov = np.ndarray((1, 1)); cov
                    cov = np.array([[cov]])
                data.append( (a,) + tuple(param) + tuple(cov.diagonal()))
            except Exception as err: 
                error.append((a, err))
        #print tabulate(data )
        if len(error)>0: 
            print error
        return data 
    
    #del this
    def fit_correlation(self, func, aa, sh, func_type='power', plot=True): 
        aa.sort(reverse=1)
        if plot: 
            res=self._fit_correlation_plot(func, aa, sh, func_type=func_type)
        else: 
            res=self._fit_correlation_noplot(func, aa, sh, func_type=func_type)
        return res

    def fit_lines(ax, func=None, fit_range=None):
        for l in ax.lines:
            l.set_linewidth(0)
        lines=list(ax.lines)
        for l in lines:
            pass
        l.set_linewidth
    
    
    @staticmethod 
    def fit_curve(x, y, func, fit_range=None): 
        param, cov = curve_fit(func, x, y)    
        
    def reload_class(self, info=1): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        if info>0: 
            print('class %s is reloaded'%(name, ))
    
    def plot_fidelity(self, aa, data,  **kwargs): 
        dic = self._plot3d(aa, aa, data, **kwargs)
        fig, ax, cb = dic['fig'], dic['ax'], dic['cb']
        #fig.colorbar(cb, ax=ax) 
        #cb.set_clim(np.min(data), 1.0)
        #cb.set_clim(0, 1.0)
        
        ax.set_aspect(1)

    def disregister_job(self, aa): 
        with rpyc_conn('local', 'service', 17012) as conn:     
            for a in aa:
                db=self[a]
                job_id=hash(os.path.normpath(db.parpath)  )
                try:
                    conn.root.disregister_job(job_id)
                except Exception as err:
                    print err, 'param=', a 
        

class Analysis_mera(Analysis): 
    def __init__(self, **kwargs): 
        kwargs.update(algorithm='mera', result_db_class=ResultDB_mera)
        Analysis.__init__(self, **kwargs)

class Analysis_vmps(Analysis): 
     def __init__(self, **kwargs): 
        #kwargs['result_db_class'] =  ResultDB_vmps
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
    
    
    def calc_EE(self, S):
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

    def measure(self, func, aa, D):
        #aa=self.alpha_list
        res=OrderedDict()
        for a in aa:
            db=self[a]
            try:
                #path='/'.join([db.parpath, 'N=0-D=%d.pickle'%D])
                #S=iDMRG.load( path)
                #S= db.load_S(path)
                S= db.load_S(sh=(0, D))
                val= func.im_func(None, S)
                res[a]=val

            except:
                #print a, path
                raise
        return res
    
    def get_EE(self, aa, D): 
        res=self.measure(self.calc_EE, aa, D=D)
        return res

    def calc_fidelity(self,aa, sh, which='fidelity_direct', eigs_tol=None, 
            save_to=None, force=0, fault_tolerant=1, info=0):
        if save_to is not None and not force: 
            res= self.load(save_to)
            return res 
        from vmps.measure_and_analysis.measurement_idmrg_mcc import fidelity, fidelity_direct
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
                    
                    f12 = db1.fetch_easy_new(which, sh, sub_key_list=[a2], fault_tolerant=1)
                    
                    #if f12 or f21: 
                    if f12 is not None and not force:  
                        ff = f12
                    else: 
                        try:
                            s1 = db1.load_S(sh);  s2=db2.load_S(sh)
                            if info>0: 
                                print 'aaa', a1, a2
                            ff = fidelity_func(s1, s2, eigs_tol=eigs_tol, info=info)
                            db1.insert(which, sh, iter=-1, sub_key_list=[a2], val=ff)
                            db2.insert(which, sh, iter=-1, sub_key_list=[a1], val=ff)
                            db2.commit()                
                            db_changed = True 
                            #self.set_rdb_dict([a2])
                        except Exception as err:
                        
                            if not fault_tolerant: 
                                print 'error occured for', a1, a2
                                raise 
                            else: 
                                faild_list.append((a1, a2, repr(err)))
                                ff=np.nan
                    data[i,j]=data[j,i]=ff
                
            if db_changed:                  
                db1.commit()   
            
        if db_changed: 
           self.set_rdb_dict(aa)  
        
        if faild_list: 
            print 'faild_list', faild_list 
        res={'aa':aa,'data':data, 'faild_list': faild_list}      
        if save_to is not None : 
            from vmps import save
            save(res, save_to)
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
                print err 
                continue 
            new_state = updater(state)
            fn = db.__class__.shape_to_backup_fn(sh)
            dir = db.parpath 
            path = '/'.join([dir, fn])
            save(new_state, path)
            print 'update file ... %s succeeded'%(path[-30: ])
   

class TestAnalsysi(unittest.TestCase): 
    def setUp(self): 
       
        parpath = '/home/zhli/backup_tensor_dir/run-ising/mera/ternary/z2-symm/scale-invar/'
        path = '/home/zhli/backup_tensor_dir/run-long-better/mera/alpha=1.5/8.pickle'
        
        an = Analysis_mera()
        aa_dic = an.scan_alpha(root=parpath)
        #print aa_dic
        an.alpha_parpath_dict.update(aa_dic)
        alpha_list = aa_dic.keys()
        an.alpha_list = sorted(alpha_list)
        
        self.an = an
    
    def test_temp(self): 
        from projects_mps.run_long_sandvik.analysis import an_idmrg_lam 
        an_idmrg_lam.an_test.scan_sub_analysis_top( )
        an_idmrg_lam.sub_analysis= []
        print_vars(vars(),  ['an_idmrg_lam.sub_analysis']) 
        an_idmrg_lam.add_sub_analysis('auto')
        print_vars(vars(),  ['an_idmrg_lam.sub_analysis']) 
        
    def test_preprocess_alpha_list(self): 
        an = self.an
        aa = [0.5, 0.3, 1.0, ]
        print an.alpha_parpath_dict.keys()
        print an.preprocess_alpha_list(aa).keys()
    
    def test_alpha_list_operations(self): 
        an = self.an 
        print  an.result_db_class 
        print an.alpha_filter(param_list=['a'], a='a>1.0', sh=(3, 4))
        aa = [1, 2, 3, 1.1, 1.01]
        res=an.alpha_filter(aa, param_list=['a'], resolution=0.1) 
        self.assertEqual(res, [1, 2, 3, 1.1])
        
    def test_plot_energy(self): 
        an = self.an
        an.plot_energy(an.alpha_list, (4, 4))
    
    def test_plot_magnetization(self): 
        pass 
        #dir = '/home/zhli/backup_tensor_dir/run-ising/idmrg-psi_guess/'
        

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
        print args
        which = args['which']
        if which in ['correlation_extra', 'correlation']: 
            sh = db.get_mera_shape_list()
            print sh
            an.plot_correlation(aa, sh, which=which, draw=1)
        #an = Analysis(path=None, parpath=args['dir'])

        if args['show']: 
            if which == 'scaling_dim' : 
                db.show_scaling_dim()
        #if 
   
    else: 
        
        if 0: 
            from result_db import ResultDB 
            for a in aaa: 
                db = ResultDB(parpath=a)
                db.upgrade()
            
       
        
        if 0:
            unittest.main()
            
        else: 
            suite = unittest.TestSuite()
            add_list = [
               TestAnalsysi('test_temp'), 
               #TestAnalsysi('test_preprocess_alpha_list'), 
               #TestAnalsysi('test_plot_energy'), 
               #TestAnalsysi('test_plot_magnetization'), 
               #TestAnalsysi('test_alpha_list_operations'), 
               
            ]
            for a in add_list: 
                suite.addTest(a)
            unittest.TextTestRunner().run(suite)
           
        

