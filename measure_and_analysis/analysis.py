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

from merapy.hamiltonian import System
from result_db import ResultDB, ResultDB_idmrg, ResultDB_vmps
from tabulate import tabulate
from measurement import mera_backup_dir_finder,  measure_all as measure_all_orig

from result_db import MATPLOTLIBRC

import matplotlib as mpl
mpl.rcParams.update(MATPLOTLIBRC)

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
        
class Analysis(object):
    """
    """
    REMOTE_HOST = 'zhihuali@211.86.151.102:'
    
    #DISTANCE_LIST = [3**i for i in range(1, 9)]
    def __init__(self, fn=None, path=None, parpath=None, backup_tensor_root='./', 
            remote_root = None, local_root=None, param_name_list=None, 
            result_db_class=None, result_db_args= {}, algorithm=None):
        self.remote_root = remote_root
        self.local_root = local_root
        self.param_name_list= param_name_list
        self.result_db_class= result_db_class if result_db_class is not None else ResultDB
        self.result_db_args= result_db_args if result_db_args is not None else ResultDB
        self.algorithm = algorithm
        
        if 1:  #remove these later
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
        
        if 1: #plotting configs
            self.invert_xaxis= 0
            self.invert_alpha_order = 0
    
    def __getitem__(self, k): 
        #try: 
        #    return self.alpha_rdb_dict[a]
        #except: 
        #    print 'all keys are %s'%(self.alpha_rdb_dict.keys())
        #    raise
        return self.alpha_rdb_dict[k]
    
    def __setitem__(self, k, v):
        self.alpha_rdb_dict[k] = v
    
    def copy(self): 
         
        res= self.__class__()
        res.alpha_parpath_dict = dict(self.alpha_parpath_dict)
        res.alpha_rdb_dict = dict(self.alpha_rdb_dict)
        return res
    
    def alpha_filter(self, a=None, h=None, V=None):  #, **kwargs): 
        """
            from wigner_crystal 
        """
        aa = list(self.alpha_list)
        if h == 1.4: 
            h = 1.4142135
        param_name_list = ['a', 'h', 'V']
        #print vars()
        for x in param_name_list: 
            i = param_name_list.index(x)
            x_val=vars()[x]
            if x_val is not None : 
                aa = filter(lambda x: x[i]== x_val, aa)
        return aa

    def scan_backup_dir(self, root='./'): 
        """     
            automatically scan and identify folders that store mera tensors
        """
        res = mera_backup_dir_finder(root)
        return res
    
    def scan_backup_dir_1(self, root='./'): 
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
    
    def scan_alpha(self, root='./', surfix='', signiture=None, info=0): 
        """
            a valid backup_parpath name is like
                    a=1.0-b=2.0 
                or with a surfix in the end 
                    a=1.0-b=2.0-1
                    a=1.0-b=2.0-a
                
        """
        dir_list = os.listdir(root)
        if info>1: 
            print 'dir_list is', dir_list
        signiture = signiture if signiture is not None else '='
        
        if 0: #bac before modify
            alpha_parpath_dict = OrderedDict()
            for name in dir_list: 
                if signiture in name: 
                    
                    alpha = name.split('=')[1]
                    try: 
                        alpha=float(alpha)
                    except: 
                        alpha = alpha
                    if surfix != '': 
                        alpha = str(alpha) + '-' +  surfix
                    #alpha_list.append(alpha)
                    alpha_parpath_dict[alpha] =  '/'.join([root, name])
                    #alpha_parpath_dict[alpha] = os.path.abspath( name)
            
            res= alpha_parpath_dict
        
        if 1: 
            def parse(xx): 
                if '=' in xx: 
                    xx = xx.split('=')[1]
                try: 
                    res=float(xx)
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
                    #alpha_list.append(alpha)
                    alpha_parpath_dict[alpha] =  '/'.join([root, name])
                    #alpha_parpath_dict[alpha] = os.path.abspath( name)
            
            res= alpha_parpath_dict
        
        return res
    
    def scan_alpha_1(self, root='./'): 
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

    def set_rdb_dict(self, alpha_list=None, create_empty_db=1, db_version=None, algorithm=None):        
        if 0: 
            if 0: 
                from result_db import ResultDB  #this is to guranteen update the def ResultDB in ipython notebook
            else: 
                import result_db
                reload(result_db)
                ResultDB = result_db.ResultDB
        else:
            module_name = self.result_db_class.__module__
            module = importlib.import_module(module_name)
            reload(module)
            self.result_db_class= module.__getattribute__(self.result_db_class.__name__)
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
    
    def load(self, path): 
        res = None
        with open(path, 'rb') as inn: 
            res=pickle.load(inn)
        return res
    
    def _outline(self, which, alpha_list=None, mera_shape_list=None, show_eng_diff=False,  grid_view=False, info=0): 
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

    def outline_file(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0): 
        self._outline('file', alpha_list=alpha_list, mera_shape_list=mera_shape_list, show_eng_diff=show_eng_diff, info=info)
    
    def outline_db(self, alpha_list=None, mera_shape_list=None, show_eng_diff=False, info=0):
        self._outline('db', alpha_list=alpha_list, mera_shape_list=mera_shape_list, info=info)

    def preprocess_alpha_list(self, alpha_list): 
        temp = self.alpha_parpath_dict.keys()
        parser =lambda a: a if not isinstance(a, str) else  float(a.split('-')[0]) 
        if 0: 
            aa = filter(lambda x: x in temp, alpha_list)
            aa_bad = filter(lambda x: x not in temp, alpha_list)
        else: 
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
            print  'aaaa', alpha_parpath
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
            self[a].delete_db()
    
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
    
    
    def fig_layout(self, ncol=1, nrow=1,  size=(3.5, 2.5)): 
        return ResultDB.fig_layout.im_func(None, ncol, nrow, size)
   
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
                    a = a[self.param_name_list.index(param_name)]
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
            r_min = None, r_max = None, 
            alpha_remap={}, alpha_shape_dic={}, fig=None,  log_scale=1,  return_fig=0, draw=0, **kwargs): 
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
                param = a[self.param_name_list.index(param_name)]
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

    def plot_entanglement_vs_alpha(self, aa, sh, alpha_remap={}, sh_remap=None, site_list=None, 
            layer=0, switch='entropy', **kwargs):
        """
            switch in [entropy, extra, brute_force, brute_force_6, brute_force_9]
        """
        #-----------------------------------------------------------------------
        print 'ploting switch=%s'%(switch, )
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
            if algorithm != 'mps':  
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
            
        for i in site_list: 
            not_found = []
            x = []
            y = []
            for a in aa:
                a_orig = a
                a1 = alpha_remap.get(a)
                a = a1 if a1 is not None else a_orig
                db=self.alpha_rdb_dict[a]
                if algorithm != 'mps':  
                    if switch in ['entropy']  + ['brute_force', 'brute_force_6', 'brute_force_9']: 
                        sub_key_list_1 = sub_key_list +  [i] 
                    else: 
                        sub_key_list_1 = sub_key_list
                    #e=db.fetch_easy('entanglement_' + switch, sh, sub_key_list=sub_key_list_1, fault_tolerant=1)
                    name = 'entanglement_' + switch
                else: 
                    pass
                    name = 'entanglement_entropy'
                    sub_key_list_1 = ['EE', i, 1]
                
                
                shape = sh_remap[a] if sh_remap.has_key(a) else sh
                e=db.fetch_easy(name, shape, sub_key_list=sub_key_list_1, info=info, fault_tolerant=fault_tolerant)
                
                if param_name is not None : 
                    a_orig = a[self.param_name_list.index(param_name)]
              
                if e is not None : 
                    x.append(a_orig)
                    y.append(e)
                else: 
                    not_found.append(a)
            
            #label = '%d site'%i if kwargs.get('label') is None else kwargs.get('label')
            args = dict(marker='.', return_ax=1, label=str(i), title='at layer=%d'%layer)
            args.update(kwargs)
            fig, ax=self._plot(x, y, **args)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
            
        if len(not_found)>0: 
            print 'not_found is %s'%(not_found, )
        
        if self.invert_xaxis: 
            ax.invert_xaxis()
        if ax.is_first_col(): 
            ax.set_ylabel('entanglement',fontsize=18)
        if ax.is_last_row(): 
            ax.set_xlabel('$\\alpha$')
        if param_name is not None : 
            ax.set_xlabel('$%s$'%param_name,fontsize=18)
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
        
    def plot_entanglement_scaling(self, aa, sh, log2=1, aver=0, num_site=None,  linear_fit=True, 
            fig=None, show_res=0, **kwargs): 
        aa=list(aa)
        aa.sort(reverse=1)
        #sh=(8,4)
        if fig is None : 
            fig=plt.figure(figsize=(5,4))
            ax=fig.add_subplot(111)
        else: 
            ax = fig.axes[0]
        not_found= []
        for a in aa:
            db=self.alpha_rdb_dict[a]
            if 0: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
            db.plot_entanglement_scaling(alpha=a, sh_list=[sh], aver=aver, 
                    log2=log2, fig=fig, num_site=num_site, linear_fit=linear_fit, show_res= show_res)
        if kwargs.get('return_fig'): 
            return fig
    
    def plot_magnetization(self, aa=None, sh=None, **kwargs): 
        if aa is None: 
            aa = list(self.alpha_list)
        if self.algorithm == 'mera': 
            key_list = [0]
        elif self.algorithm == 'idmrg': 
            key_list = None
        else: 
            raise NotImplemented
        
        x=[]; y=[]
        for a in aa:
            db=self.alpha_rdb_dict[a]
            rec=db.fetch_easy('magnetization', sh, sub_key_list=key_list)
            x.append(a)
            y.append(rec)
            
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

    def reload_class(self): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        print('class %s is reloaded'%(name, ))

class Analysis_mera(Analysis): 
    def __init__(self, **kwargs): 
        kwargs['algorithm'] = 'mera'
        Analysis.__init__(self, **kwargs)

class Analysis_vmps(Analysis): 
     def __init__(self, **kwargs): 
        #kwargs['result_db_class'] =  ResultDB_vmps
        kwargs.update(algorithm='mera', result_db_class=ResultDB_vmps)
        Analysis.__init__(self, **kwargs)

class Analysis_idmrg(Analysis): 
    def __init__(self, **kwargs): 
        """
            def __init__(self, fn=None, path=None, parpath=None, backup_tensor_root='./', 
                remote_root = None, local_root=None, param_name_list=None, 
                result_db_class=None, result_db_args= {}):
        """        
        kwargs.update(algorithm='idmrg', result_db_class=ResultDB_idmrg)
        Analysis.__init__(self, **kwargs)
        #super(Analysis, self).__init__(**kwargs)
    
    def get_energy(self, aa, D):
        res=OrderedDict()
        for a in aa:
            db=self[a]
            try:
                path='/'.join([db.parpath, 'N=0-D=%d.pickle'%D])
                S=self.load( path)
                eng=S['energy']
                res[a]=eng
            except Exception as err:
                print a, err
                
        return res
    
    def calc_EE(self, S):
        res=None
        #lam = S['Lambda']
        lam = S['lam']
        U, s, V = np.linalg.svd(lam)
        spect= s
        #print np.sum(s**2)
        spect2=spect**2
        EE = -np.sum(spect2*np.log(spect2))
        res = EE
        return res 

    def measure(self, func, aa, D):
        #aa=self.alpha_list
        res=OrderedDict()
        for a in aa:
            db=self[a]
            try:
                path='/'.join([db.parpath, 'N=0-D=%d.pickle'%D])
                #S=iDMRG.load( path)
                S= db.load_S(path)
                val= func.im_func(None, S)
                res[a]=val

            except:
                #print a, path
                raise
        return res
    
    def get_EE(self, aa, D): 
        res=self.measure(self.calc_EE, aa, D=D)
        return res



class TestAnalsysi(unittest.TestCase): 
    def setUp(self): 
       
        parpath = '/home/zhli/Documents/mera_backup_tensor/run-ising/ternary/z2-symm/scale-invar/'
        path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=1.5/8.pickle'
        
        an = Analysis()
        aa_dic = an.scan_alpha(root=parpath)
        #print aa_dic
        an.alpha_parpath_dict.update(aa_dic)
        alpha_list = aa_dic.keys()
        an.alpha_list = sorted(alpha_list)
        
        self.an = an
    
    def test_preprocess_alpha_list(self): 
        an = self.an
        aa = [0.5, 0.3, 1.0, ]
        print an.alpha_parpath_dict.keys()
        print an.preprocess_alpha_list(aa)
    
    def test_plot_energy(self): 
        an = self.an
        an.plot_energy(an.alpha_list, (4, 4))


if __name__ == '__main__':
    #if 0: 

    
    #else: 
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
               #TestAnalsysi('test_preprocess_alpha_list'), 
               TestAnalsysi('test_plot_energy'), 
            ]
            for a in add_list: 
                suite.addTest(a)
            unittest.TextTestRunner().run(suite)
           
        

