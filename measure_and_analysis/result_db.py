#!/usr/bin/env python
#coding=utf8
#PYTHON_ARGCOMPLETE_OK 

import argparse, argcomplete
import os, sys
import pprint
from tabulate import tabulate
import cPickle as pickle
from collections import OrderedDict
from operator import itemgetter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import pandas
from scipy.optimize import curve_fit
import unittest

from merapy.hamiltonian import System
from merapy.context_util import rpyc_load, rpyc_save

if 1: 
    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = (4, 3)
    mpl.rcParams['axes.labelsize'] = 16
    #mpl.rcParams['legend.frameon'] = False  # whether or not to draw a frame around legend   
    MATPLOTLIBRC = {
            'legend.alpha': 0,    # this will make the box totally transparent 
            'legend.edege_color': 'white',   # this will make the edges of the border white to match the background instead 
            'legend.frameon': 0,# whether or not to draw a frame around legend    
            
            }
    mpl.rcParams.update(MATPLOTLIBRC)


FIELD_NAME_LIST = [
    'central_charge', 'scaling_dim', 'correlation', 'correlation_extra',  
    'entanglement_entropy', 'entanglement_spectrum', 'magnetization', 
    'entanglement_brute_force', 'entanglement_brute_force_6', 
    # these are too slow and not useful,  disable them
    #'entanglement_brute_force_9',   
    #'entanglement_brute_force_aver', 
    #'entanglement_brute_force_6_aver',
    #'entanglement_brute_force_9_aver' , 
    ]


#class ResultDB_Utills

class ResultDB(OrderedDict): 
    """
        a DB to store measurement results 
        
        the structure of the db,  that gives the instruction of the key chains to fetch recs
        mera version = 1.0
            STRUCT = {
                    'entanglement': [shape, iter, layer, num_sites], 
                    }
        MPS, version = 1.0
            STRUCT = {
                    'entanglement': ['shape', 'iter', num_sites], 
                    }
        
    """
    DBNAME = 'RESULT.pickle.db'   #fix the name of db
    def __init__(self, parpath, dbname=None, version=None, algorithm=None, 
            use_local_storage=False, create_empty_db=False,  upgrade=0):
        """
            params: 
                algorithm is necessary, or else it is difficult to determine which alg
        """
        self.use_local_storage = use_local_storage
        self.parpath = parpath
        dbname = dbname if dbname is not None else self.DBNAME
        self.dbname = dbname
        self.path = '/'.join([parpath , dbname])
        
        #self.inited = False  
        self.version = version
        self.upgrade = upgrade
       
        if 0: 
            if not os.path.exists(parpath):  # these lines are useful, when merge remote db
                self.create_parpath()
       
        if not os.path.exists(self.path):  
            if create_empty_db:     
                ResultDB.create_empty_db(self.path, use_local_storage) 
            
            OrderedDict.__init__(self, {})
        else: 
            db = rpyc_load(self.path, use_local_storage)
            OrderedDict.__init__(self, db)
        
        if algorithm is not None : 
            self['algorithm'] = algorithm
        version = self.version if self.version is not None else self.get('version')
        self.version = version
        self['version'] = version
        if version  == 1.0: 
            self.fetch_easy = self.fetch_easy_new
        if self.upgrade: 
            if self.get('version') != 1.0 and dbname != self.dbname + '.new': 
                self.upgrade()
    
    def __repr__(self): 
        return self.keys()
    
    def __str__(self): 
        return str(self.keys())
    
    if 0: 
        def __getitem__(self, k): 
            """
                
            """
            print 'iiiii'
            if not self.inited: 
                self.init_db()
                print 'iiiii'
            else: 
                #return self[k]
                return OrderedDict.__getitem__(self, k)

        def __getitem__(self, k): 
            return self._db[k]
        
        def __setitem__(self, k, v): 
            self._db[k] = v
            pass
    
    def load_S(self, path): 
        with open(path, 'rb') as inn: 
            res=pickle.load(inn)
            return res
    
    def reload_class(self): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        print('class %s is reloaded'%(name, ))
    
    def create_parpath(self): 
        cmd = 'mkdir %s'%self.parpath
        os.system(cmd)
        msg = 'parpath not exist, so %s'%cmd
   
    def has_key_list(self, key_list, info=1): 
        temp = self
        for k in key_list: 
            #print 'kkk', k, key_list
            if not hasattr(temp, 'has_key'):   # val of rec
                if info>0: 
                    print 'not has this key: %s'%(k, )
                return False
            if not temp.has_key(k): 
                if info>0: 
                    print 'not has this key: %s'%(k, )
                return False
            temp = temp[k]
        
        return True
   
    def add_key_list(self, key_list, val=None, verbose=0): 
        temp = self
        for k in key_list[: -1]: 
            if not temp.has_key(k): 
                temp[k] = OrderedDict()
            temp = temp[k]
        last_key = key_list[-1]
        if val is None: 
            temp[last_key] = OrderedDict()
        else: 
            temp[last_key] = val
        #print self.has_key_list(key_list)
        #exit()
        if verbose: 
            print 'add_key_list %s'%(key_list, )
                
    @classmethod
    def create_empty_db(cls, path, use_local_storage=False): 
        if 0: 
            out = open(path,  'wb')  #raise #useless , just use 'ab+' mode
            db = OrderedDict()
            pickle.dump(db, out)
            out.close()
        db = OrderedDict()
        rpyc_save(path, db, use_local_storage=use_local_storage)
        print 'ResultDB %s has been created'%path
    
    def get_fn_list(self):
        nlist = os.listdir(self.parpath)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" == i.split('.')[-1]]
        #print "all pickle files under specified path are %s"%pickle_files
        if not self.get('algorithm')=='mps': 
            pickle_files= [i for i in pickle_files if 'trans' not in i]
        
        return pickle_files
    
    @staticmethod
    def parse_fn_del(fn):
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

    def parse_fn(self, fn):
        if self.get('algorithm') in ['mps', 'idmrg']: 
            ND=fn.split('.')[0]
            N, D = ND.split('-')
            N = N.split('=')[1]; N=eval(N)
            D = D.split('=')[1]; D=eval(D)
            res = N, D
        else:     
            res = System.backup_fn_parse(fn)
        return res

    def get_mera_shape_list(self, sh_max=None, sh_min=None): 
        fn = self.get_fn_list()
        #print 'backup tensor files are %s'%fn
        sh = [self.parse_fn(f) for f in fn]
        #print sh; exit()
        if sh_max is not None : 
            sh=filter(lambda x: x[0]<=sh_max[0] and x[1] <= sh_max[1] , sh)
        if sh_min is not None : 
            sh=filter(lambda x: x[0]>= sh_min[0] and x[1]>= sh_min[1] , sh)
        sh.sort()
        return sh
    
    get_shape_list = get_mera_shape_list  #a shorter name
    
    def get_energy_min_bac(self, sh=None, **kwargs): 
        version = self.get('version')
        
        if version is None: 
            if 1: 
                eng_list = []
                for k, v in self.iteritems(): 
                    #print k
                    eng=v.get('energy')
                    if eng is not None : 
                        eng_list.append((k, eng)) 
                res= min(eng_list, key=lambda x: x[1])
            else: 
                pass
            
        elif version == 1.0: 
            #self['energy'][sh]
            raise NotImplemented 
        
        else: 
            raise
        
        if kwargs.get('show_val'): 
            print 'energy min is: ', res  
        return res
    
    def get_energy_min(self, sh=None, **kwargs): 
        version = self.get('version')
        
        if version is None: 
            if 0: 
                eng_list = []
                for k, v in self.iteritems(): 
                    #print k
                    eng=v.get('energy')
                    if eng is not None : 
                        eng_list.append((k, eng)) 
                res= min(eng_list, key=lambda x: x[1])
            if 1: 
                recs= self.get_energy(sh).itervalues()
                res= min(recs, key=lambda x: x[0])
                res= res[0]
                    
                pass
            
        elif version == 1.0: 
            #self['energy'][sh]
            raise NotImplemented 
        
        else: 
            raise
        
        if kwargs.get('show_val'): 
            print 'energy min is: ', res  
        return res
    
    def get_energy_fluc_count(self, sh, bins=None): 
        if bins is None: 
            bins = [-1, 0,  1e-8, 1e-7, 1e-6, 1e-5, 1] 
     
        diff=self.get_energy( sh=sh, which='diff', iter_max=None)
        count, bin = np.histogram(diff, bins=bins)
        return count, bins
    
    def commit(self, info=0):
        #issue: should change to 'ab+?'
        if 0: 
            out = open(self.path, 'wb')
            #out = open(self.path, 'wb')
            #pickle.dump(self, out)  #this is rather critical, or else got error
            pickle.dump(OrderedDict(self), out)
            out.close()
        rpyc_save(self.path, OrderedDict(self), use_local_storage=self.use_local_storage)
        if info>0: 
            temp = str(self.path)
            print '\tmodification in ...%s has been commited'%(temp[-40: ])
    
    def _fetch(self, dim, num_of_layer=None): 
        keys= self.iterkeys()
        #keys= filter(lambda x: x[0]==dim and x[1]==num_of_layer, keys)
        keys= filter(lambda x: x[0]==dim, keys)
        if num_of_layer is not None: 
            keys= filter(lambda x: x[1]==num_of_layer, keys)
        if len(keys)<1: 
            return None
        else: 
            #print 'kkk', keys
            #print self.keys()
            try: 
                res = max(keys, key=itemgetter(2))  #records with max S.iter
            except: 
                print 'parpath %s'%(self.parpath, )
                print 'keys are, %s'%keys
                raise
            return res
    
    def fetch(self, field_name, dim, num_of_layer=4, step=None, info=1): 
        key = self._fetch(dim, num_of_layer)
        if key is None: 
            if info>0: 
                msg = 'key (%s) not found in ResultDB  at %s\n keys are %s'%(
                    (dim, num_of_layer), self.parpath,  self.keys())
                print msg
            return None
        else: 
            if self[key].has_key(field_name): 
                return  self[key][field_name]
            else: 
                if info>0: 
                    msg = "rec not found for field_name '%s' at %s"%(field_name, key)
                    print msg
                return None

    def fetch_easy(self, field_name, mera_shape,  sub_key_list=None, mera_shape_optional=None, 
             auto_meassure=False, info=0, fault_tolerant=1): 
        """
            interface to make fetch easier
        """
        
        dim, layer = mera_shape[:2]
        if 1: 
            mapper = {
                    'EE':  ['central_charge', 'EE'  ], 
                    'EE1': ['central_charge', 'EE', 1], 
                    'EE2': ['central_charge', 'EE', 2], 
                    }
            if mapper.has_key(field_name): 
                temp = mapper[field_name]
                field = temp[0]
                key_list = temp[1: ]
            else: 
                field = field_name
                key_list = None
            if sub_key_list is not None: 
                if key_list is None: 
                    key_list = sub_key_list
                else: 
                    key_list += sub_key_list
        
        
        rec = self.fetch(field, dim, num_of_layer=layer, info=info)
    
        if key_list is not None and rec is not None: 
            try: 
                for key in key_list: 
                    rec = rec[key]
            except: 
                if fault_tolerant: 
                    print 'error in fetch_easy'
                    print 'keys are', rec.keys()
                else: 
                    raise
        return rec

    def fetch_easy_new(self, field_name, mera_shape,  sub_key_list=None, mera_shape_optional=None, 
             auto_meassure=False, info=0, fault_tolerant=1): 
        """
            fetch_easy for version 1.0
        """
        try: 
            rec = self[field_name][mera_shape]
            key = next(reversed(rec))  #last iter step
            rec = rec[key]
        except KeyError as err: 
            if fault_tolerant: 
                if info>0: 
                    print 'key not found:', err
                #print 'error in fetch_easy'
                #print 'keys are', rec.keys()
                return None
            else:
                raise
        except: 
            raise
        
        if 1: 
            mapper = {
                    'EE':  ['central_charge', 'EE'  ], 
                    'EE1': ['central_charge', 'EE', 1], 
                    'EE2': ['central_charge', 'EE', 2], 
                    }
            if mapper.has_key(field_name): 
                temp = mapper[field_name]
                field = temp[0]
                key_list = temp[1: ]
            else: 
                field = field_name
                key_list = None
            if sub_key_list is not None: 
                if key_list is None: 
                    key_list = sub_key_list
                else: 
                    key_list += sub_key_list
        
        if key_list is not None and rec is not None: 
            try: 
                for key in key_list: 
                    rec = rec[key]
            except: 
                if fault_tolerant: 
                    print 'error in fetch_easy'
                    print 'keys are', rec.keys()
                else: 
                    raise
        return rec
    
    def insert(self, field_name, sh, iter, val): 
        if not self.has_key(field_name): 
            self[field_name] = OrderedDict()
        if not self[field_name].has_key(sh): 
            self[field_name][sh] = OrderedDict()
        self[field_name][sh][iter] = val

    def put(self, field_name, key, res, sub_key_list=None): 
        self[key].update({field_name:res})
    
    #def has_entry(self, field_name, mera_shape, iter): 
    def has_entry(self, key_list): 
        #key_list = [field_name, mera_shape, iter]
        rec = self
        for t in key_list: 
            if not rec.has_key(t): 
                return False
            else: 
                rec = rec[t]
        return True

    def fit_correlation(self, sh, func_type='power', which='correlation_extra', 
            direct='pm', r_min=None, r_max=None, plot=0, plot_orig=1, show_label=1, **kwargs): 
        power_func= lambda x, C, eta: -C*x**eta*(-1)**x
        #power_si_func= lambda x,  C,eta, m: -C*x**eta*sici(pi*m*x)[0]*(-1)**x
        power_si_func= lambda x,  C,eta, m: -C*x**eta*sici(m*x)[0]*(-1)**x
        power_cos_func= lambda x,  C,eta, m, b, k: -C*x**eta*(-1)**x * (b+ np.sin(1.*x/k))   #not usefull, only test
        def exp_func_1(t, a, b, c):                                                         
            return a * np.exp(b * t) + c   

        def exp_func_2(t, A, B):                                                         
            #return A * np.exp(-B * t) 
            return A * np.exp(- t/B) 

        if func_type == 'power': 
            func = power_func; 
        elif func_type == 'exp' : 
            func = exp_func_2
        
        if 1: 
            xy = self.get_correlation( sh, which, direct=direct, r_min=r_min, r_max=r_max,  period=None)
            if xy is None: 
                return None
            else: 
                x, y = xy
        
        if kwargs.get('info'): 
            print 'xxx', x[0], x[-1], len(x)
        x=np.array(x); y=np.array(y);  
        #if direct == 'pm':  y=y*2;  
        y=np.abs(y)
            #print x, y; exit()
        #x, y = self.get_correlation(which=which, direct=direct)
        fit_succeed = 1
        try:
            param,cov = curve_fit(func, x, y)
            if isinstance(cov, float): 
                #cov = np.ndarray((1, 1)); cov
                cov = np.array([[cov]])
            #data.append( (a,) + tuple(param) + tuple(cov.diagonal()))
            
        except Exception as err: 
            #raise
            print 'fitting faild', err 
            #error.append((a, err))
            fit_succeed = 0
            #return None
        if fit_succeed: 
            if plot: 
          
                label = '' if not show_label else param[1]
                fig=self._plot(x, func(x, *param), xscale='log',yscale='log', label=label, return_fig=1, **kwargs)
                if plot_orig: 
                    if not kwargs.has_key('fig'): kwargs['fig'] = fig
                    fig=self._plot(x, y, xscale='log', yscale='log',   **kwargs)
            #if return_fig
            return param
       
        else: 
            return None
            
        
    def get_correlation(self, sh, which='correlation_extra',direct='pm', 
            r_min=None,  r_max=None, period=None, fault_tolerant=1,  info=0): 
        if which == 'correlation_mixed' and self.get('algorithm') != 'mps': 
            rec1=self.fetch_easy('correlation_extra', mera_shape=sh, sub_key_list=[direct])
            rec2=self.fetch_easy('correlation', mera_shape=sh, sub_key_list=[direct, 'val'])
            if rec1 is None  or rec2 is None: 
                #not_found.append(sh)
                return None
            if 0:    
                if log_scale: 
                    x1, y1=zip(*rec1[0:-1:2])
                else: 
                    x1, y1 = zip(*rec1);  
            period  = period if period is not None else 2
            x1, y1=zip(*rec1[0:-1:period])
            x2, y2 = zip(*rec2)
            x1 = np.array(x1);    x2 = np.array(x2);  
            y1 = np.array(y1);    y2 = np.array(y2);  
            x1_max = np.max(x1);  temp = np.where(x2>x1_max) 
            x = np.append(x1, x2[temp])
            y = np.append(y1, y2[temp])
        else: 
            if which == 'correlation': 
                key_list = [direct, 'val'] 
                period  = period if period is not None else 1                
            elif which == 'correlation_extra': 
                key_list = [direct] 
                period  = period if period is not None else 2            
            if self.get('algorithm') == 'mps': 
                which = 'correlation'  #not extra
                key_list = [direct]
                period = period if period is not None else 2            
            rec = self.fetch_easy(which, mera_shape=sh, sub_key_list=key_list, 
                    fault_tolerant=fault_tolerant, info=info-1) 
            if rec is None: 
                return None
            
            if self.get('algorithm') == 'mps': 
                rec = map(lambda x: (x[1]-x[0], x[2]), rec) 
            
            x,y=zip(*rec[0:-1:period])
        
        if direct == 'pm':   
            y=np.array(y);  y=y*2;  
        
        x = np.array(x)
        a = 0
        b = -1
        if r_min is not None : 
            temp = np.where(x>=r_min);  a=temp[0][0]
        if r_max is not None : 
            temp = np.where(x<=r_max);  b=temp[0][-1] + 1
            x = x[a: b]
            y = y[a: b]
           
                   
        return x, y    
    
    def delete_rec(self, key): 
        self.pop(key)
        self.commit()
    
    def delete_db(self): 
        os.remove(self.path)
        msg = 'db is removed'
        print(msg)
    
    def dump(self, file=None): 
        print '---'*30
        print 'all records in %s'%self.dbname
        #print self.items()
        pprint.pprint(self, stream=sys.stdout)

    def merge(self, other, info=0): 
        """
            merge other into self
        """
        if self.get('verion') is not None: 
            raise # need modify 
        self.update(other)
        self.commit()
        if info>0: 
            msg = 'other rdb is merged into self'
            print msg
        
    def merge_remote(self, parpath, remote_host=None): 
        remote_host = 'zhihuali@211.86.151.102:'
        #remote_root = 'current/run-long-better/'
        remote_path = remote_host + parpath   + '/RESULT.pickle.db'
        locale_parpath = self.parpath
        locale_path = locale_parpath + 'temp.pickle'
        cmd_str = 'scp %s %s'%(remote_path, locale_path)
        print cmd_str
        status=os.system(cmd_str)
        if status != 0: 
            msg = 'scp %s ---> %s'%(remote_path[-20:], locale_path[-20: ])                     
            print 'scp failed'
        else: 
            other = ResultDB(self.parpath, 'temp.pickle')
            self.merge(other)
            os.system('rm %s'%locale_path)
            #os.system('ls %s'%self.parpath)
            print 'merge succeed'
    
    
    def fig_layout(self, ncol=1, nrow=1,  size=(3.5, 2.5)): 
        size = ncol*size[0], nrow*size[1]
        fig = plt.figure(figsize=size)
        iii=0
        for i in range(nrow): 
            for j in range(ncol): 
                iii+=1
                fig.add_subplot(nrow, ncol, iii)
        fig.tight_layout()
        return fig

    def _plot(self, x, y, **kwargs): 
        figsize = kwargs.get('figsize')
        figsize = (3.5, 2.5) if figsize is None else figsize
        fig = kwargs.get('fig', None)
        ax = kwargs.get('ax', None)
        
        if fig is None and ax is None: 
            fig=plt.figure(figsize=figsize)
        if ax is None: 
            if len(fig.axes)==0: 
                ax=fig.add_subplot(111)
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
                    
        style_default = {
                'marker':'.', 
                'label': '', 
                'color': None, 
                }    
        
        temp = ['marker', 'label', 'color']
        dic = {i:kwargs.get(i) for i in temp if kwargs.has_key(i)  }
        
        if not dic.has_key('label'): 
            dic['label'] = ''
        if not dic.has_key('marker'): 
            dic['marker'] = '.'
        
        ax.plot(x, y, **dic)
        
        temp = ['xlabel', 'ylabel', 'xlim', 'ylim', 'xscale', 'yscale']
        for t in temp: 
            tt = kwargs.get(t)
            if tt is not None :
                if isinstance(tt, tuple): 
                    ax.__getattribute__('set_' + t)(*tt)
                else: 
                    ax.__getattribute__('set_' + t)(tt)
        if kwargs.has_key('title') and ax.is_first_row(): 
            ax.set_title(kwargs.get('title')) 
        
        ax.legend(title=kwargs.get('legend_title'))
        ax.grid(1)
        if kwargs.get('show_fig') : 
           plt.show() 
        if kwargs.get('return_ax', False): 
            return fig, ax
        else: 
            #return fig
            return ax.figure


    def plot_central_charge_vs_layer(self, sh, alpha_remap={}, alpha_sh_map={}, layer='all',  **kwargs):
        if layer=='top':
            ll= [sh[1]-2 ] 
        elif layer=='all':
            ll=  range(sh[1]-1)
        else: 
            ll = layer
        
        #ll= [sh[1]-2 ] if layer is None else range(sh[1]-1)
        not_found = []
        cc=[]
        for layer in ll:
            if 1: 
                #rec = db.fetch_easy('entanglement_entropy', shape, [layer])
                rec = self.fetch_easy('entanglement_entropy', sh)
                if rec is not None : 
                    rec = rec.get(layer)
                if rec is None: 
                    #not_found.append((a, layer))
                    not_found.append( layer)
                    continue
                #print 'rrrr', rec
                e1, e2 = rec[1], rec[2]
                #print a, e2, e1
                c= (e2-e1)*3
                cc.append((layer,c))
        #if len(cc)==0: 
        #    fig, ax = None, None
        #else: 
        if 1: 
            x,y=zip(*cc)    
            kwargs.update(marker='.', return_ax=1)
            fig, ax = self._plot(x, y, **kwargs)
            #ax.plot(x,y, '.-')    
            ax.legend()
            ax.set_xticks(ll)
            #ax.invert_xaxis()
            ax.grid(1)
        if len(not_found)>0: 
            print 'rec for central_charge not_found is :%s'%not_found
        if kwargs.get('return_ax'): 
            return fig, ax
        if kwargs.get('return_fig'): 
            return fig
        if kwargs.get('show_fig'): 
            plt.show()
   
    
    def plot_entanglement_scaling_old(self, alpha=None, sh_list=None, sh_min=None, sh_max=None, 
            log2=False, num_site=None, aver=False, linear_fit=True,  fit_points= None, 
            **kwargs): 
        #fig=None, show_res=0, return_fig=0, info=0): 
        fig = kwargs.get('fig')
        ax = kwargs.get('ax')
        show_res= kwargs.get('show_res')
        return_fig = kwargs.get('return_fig')
        info = kwargs.get('info', 0)
        if sh_list is None: 
            sh_list = self.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
            if len(sh_list)==0: 
                raise
        if fig is None : 
            fig=plt.figure(figsize=(5,3))
            ax=fig.add_subplot(111)
         
        if len(fig.axes)==0: 
            ax=fig.add_subplot(111)
        else:
            for i in fig.axes: 
                if len(i.lines)==0:
                    ax = i
                    break
            #if ax is None: 
            #    ax = fig.axes[-1]
                
            
        not_found= []
        algorithm =  self.get('algorithm') 
        
        for sh in sh_list: 
            if self.get('algorithm')=='mps': 
                rec = self.fetch_easy('entanglement_entropy', sh, ['EE'], info=info-1)
                if rec is None: 
                    not_found.append(sh)
                    continue
                rec = rec[: len(rec)/2-1]
                x,y=(zip(*rec))
                
            else: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
                if aver: name += '_aver' 
                rec=self.fetch_easy(name, sh, ['EE'],)
                if rec is None: 
                    not_found.append(sh)
                    continue
                x,y=(zip(*rec.items()))
            
            if show_res: 
                print alpha, sh, rec
            x=np.array(x); 
            y=np.array(y); 
            label = ''
            if alpha is not None : label += '%s'%(alpha, )
            if len(sh_list)>1: label += ' %s'%(sh, )
            
            if fit_points is None: 
                xx = x;  yy = y
            else: 
                xx = x[fit_points]
                yy = y[fit_points]
            
            if linear_fit: 
                c, b = np.polyfit(np.log2(xx.real[1: -1]), yy.real[1: -1], deg=1)
                c = c*3
                #print 'central_charge: ', label, c*3
            else: 
                func = lambda x, c, const: c/3.*np.log2(x) + const
                param,cov = curve_fit(func, xx.real[1: -1], yy.real[1: -1])
                #print alpha, sh, param.round(3),  cov.diagonal()
                c = param[0]
            
            label += ' %1.2f'%c 
            ax.plot(x,y, '.-', label=label)
            
        l=ax.legend(loc='lower right')    
        if rec is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transparent
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
        
        if log2 and rec is not None: 
            #ax.set_xticklabels(2**x)
            #tic = x.tolist()
            ax.set_xscale('log')
            if algorithm != 'mps':  
                tic = x[1: ]   # x[0] = 0, wont show in log scale
                ax.set_xticks(tic)
                ax.set_xticklabels(tic)
        ax.grid(1)
        
        if return_fig: 
            return fig
           
    def plot_entanglement_scaling(self, alpha=None, sh_list=None, sh_min=None, sh_max=None, 
            log2=False, num_site=None, aver=False, linear_fit=True,  fit_points= None, 
            **kwargs): 
       
        
        show_res= kwargs.get('show_res')
        return_fig = kwargs.get('return_fig')
        info = kwargs.get('info', 0)
        if sh_list is None: 
            sh_list = self.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
            if len(sh_list)==0: 
                raise
        if 0: 
            if fig is None : 
                fig=plt.figure(figsize=(5,3))
                ax=fig.add_subplot(111)
             
            if len(fig.axes)==0: 
                ax=fig.add_subplot(111)
            else:
                for i in fig.axes: 
                    if len(i.lines)==0:
                        ax = i
                        break
                #if ax is None: 
                #    ax = fig.axes[-1]
                
            
        not_found= []
        algorithm =  self.get('algorithm') 
        
        for sh in sh_list: 
            if self.get('algorithm')=='mps': 
                rec = self.fetch_easy('entanglement_entropy', sh, ['EE'], info=info-1)
                if rec is None: 
                    not_found.append(sh)
                    continue
                rec = rec[: len(rec)/2-1]
                x,y=(zip(*rec))
                
            else: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
                if aver: name += '_aver' 
                rec=self.fetch_easy(name, sh, ['EE'],)
                if rec is None: 
                    not_found.append(sh)
                    continue
                x,y=(zip(*rec.items()))
            
            if show_res: 
                print alpha, sh, rec
            x=np.array(x); 
            y=np.array(y); 
            label = ''
            if alpha is not None : label += '%s'%(alpha, )
            if len(sh_list)>1: label += ' %s'%(sh, )
            
            if fit_points is None: 
                xx = x;  yy = y
            else: 
                xx = x[fit_points]
                yy = y[fit_points]
            
            if linear_fit: 
                c, b = np.polyfit(np.log2(xx.real[1: -1]), yy.real[1: -1], deg=1)
                c = c*3
                #print 'central_charge: ', label, c*3
            else: 
                func = lambda x, c, const: c/3.*np.log2(x) + const
                param,cov = curve_fit(func, xx.real[1: -1], yy.real[1: -1])
                #print alpha, sh, param.round(3),  cov.diagonal()
                c = param[0]
            
            label += ' %1.2f'%c 
            #ax.plot(x,y, '.-', label=label)
            kwargs['label'] = label
            kwargs.update(label=label, marker='.', return_ax=1)
            fig, ax=self._plot(x,y, **kwargs)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
        
        if 'ax' in locals():          
            l=ax.legend(loc='lower right')    
            if rec is not None : 
                l.get_frame().set_alpha(0) # this will make the box totally transparent
                l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
            
            if log2 and rec is not None: 
                #ax.set_xticklabels(2**x)
                #tic = x.tolist()
                ax.set_xscale('log')
                if algorithm != 'mps':  
                    tic = x[1: ]   # x[0] = 0, wont show in log scale
                    ax.set_xticks(tic)
                    ax.set_xticklabels(tic)
            ax.grid(1)

        if return_fig and 'fig' in locals(): 
            return fig
           
       
    def plot_entanglement_scaling_compare(self, num_site_max=6): 
        """
            only for mera
        """
        for n, a in [(None, False),(6, False), (6, True), (9, False), (9,True)]:
            if n>num_site_max: 
                continue
            fig=self.plot_entanglement_scaling(num_site=n, sh_max=(12,6), aver=a, log2=1, return_fig=1)
            fig.axes[0].set_title('%s-%s'%(n,a))
    
    def plot_all(self, filed_name_list, sh, **kwargs): 
        fig = kwargs.get('fig')
        name_list= filed_name_list
        if fig is None: 
            fig=plt.figure(figsize=(4*len(name_list),3))
        i=0
        for n in name_list:
            i+=1
            ax=fig.add_subplot(1,len(name_list), i)
            #self.plot_common(n, sh=(4,4), fig=fig,direct=direct)
            self.plot_common(n, sh=(4,4), fig=fig, **kwargs)
            #if aa.index(a)==0: ax.set_title(n, fontsize=16)
            #if i==1: ax.set_ylabel('V=%s'%(a[2]), fontsize=16)
        fig.tight_layout()
        if kwargs.get('return_fig'): 
            return fig

    def plot_common(self, field_name, sh, **kwargs): 
        methname = 'plot_' + field_name
        if hasattr(self, methname): 
            self.__getattribute__(methname)(sh, **kwargs)
    
    def plot_magnetization(self, sh, shift_y=0, plot_mean=0, **kwargs): 
        rec=self.fetch_easy('magnetization', sh)
        if rec is None: 
            return None
        x, y = zip(*rec.items())
        if shift_y: 
            y = np.asarray(y)
            y = (y + 1.0)/2.0
        fig=self._plot(x, y, **kwargs)
        if plot_mean: 
            mean = np.ndarray(len(y))
            mean[: ] = np.mean(y)
            fig = self._plot(x, mean, fig=fig)
        
        if kwargs.get('return_fig'): 
            return fig
    
    def get_magnetization_mean(self, sh, **kwargs): 
        rec=self.fetch_easy('magnetization', sh)
        if rec is None: 
            return np.nan
        x, y = zip(*rec.items())
        y = np.asarray(y)
        #y = (y + 1.0)/2.0
        return y.mean()  
    
    def plot_correlation(self, sh,  direct='pm', r_min=None, r_max=None, which='correlation_mixed', 
            log_scale=1, fit=None, plot_fit=False,   return_fig=0, draw=0, **kwargs): 
        
        if 1: 
            assert which in ['correlation_extra', 'correlation', 'correlation_mixed']
            sh_list = sh
            fig = kwargs.get('fig')
            ax = kwargs.get('ax')
            if ax is None: 
                if fig is None :
                    fig=plt.figure(figsize=(6, 4))
                if len(fig.axes)==0: 
                    ax=fig.add_subplot(111)
                else: 
                    ax = fig.axes[-1]
        not_found=[]
        is_many_sh = 0
        if len(sh_list)==2 : # and isinstance(sh_list[0], int): 
            sh_list = [sh_list]
            #is_many_sh = 1
        #direct_orig = direct
        #if direct == 'pm': direct = 'xx' 
        sub_key_list = [direct] if which == 'correlation_extra' else [direct, 'val'] 
        for sh in sh_list: 
            
            period = 2 if log_scale else 1
                
            rec = self.get_correlation(sh, which, direct, r_min=r_min, r_max=r_max, period=period )
            if rec is not None: 
                x, y = rec
            else: 
                continue
            
            if direct == 'pm':   
                y=np.array(y);  y=y*2;  

            label = kwargs.get('label')
            
            label = label if label is not None else ''
            #label = ''  #if len(sh_list)<2: label += ' %s  '%(a, )   #if len(alpha_list)<2: label += ' %s  '%(sh, )
            if is_many_sh: label += ' %s  '%(sh, )
            marker = kwargs.get('marker')
            marker = marker if marker is not None else '--'
            if fit is None and len(marker)==1: marker  += '-'
            
            if log_scale: 
                y=np.abs(y)
                #ax.loglog(x, y, marker, label=label)
                ax.loglog(x, y, marker)
            else: 
                ax.plot(x, y)
            
            line = ax.lines[-1]
            color = line.get_color()
            line.set_markersize(5)
            line.set_mfc('w')
            line.set_markeredgewidth(1)
            line.set_mec(color)
            if fit is not None : 
                    #line.set_mfc('w')
                    #line.set_markersize(5)
                res=self.fit_correlation(sh, func_type=fit,  which=which, direct=direct, 
                        r_min=r_min, r_max=r_max, plot=plot_fit, plot_orig=0, show_label=0, fig=fig, 
                        color=color)
                if res is not None : 
                    #if fit == 'power': 
                    A, eta  = res
                    if fit == 'power':  
                        name = '\\eta' 
                    else: 
                        name = '\\xi'
                    label += ' $%s=%1.2f$'%(name, eta)
                       
            line.set_label(label)
            if fit: ax.set_xlim(x[0], x[-1]*10)
        l=ax.legend()   #l=ax.legend(title='$\\alpha$')
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transparent
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
           
        ax.grid(1)
        if kwargs.get('show_fig'): plt.show()
        if draw: plt.show()
        if len(not_found)>0: print 'not_found is ', not_found
        if return_fig: return fig

    def plot_structure_factor_bac(self, sh, direct=None, dist_max=100, show_fig=0, save_fig=1, **kwargs): 
        fft = np.fft.fft
        algorithm = self.get('algorithm')
        if algorithm == 'mps': 
            name = 'correlation' 
            direct = 'xx' if direct is None else direct
        else: 
            name = 'correlation_extra'
            direct = 'zz' if direct is None else direct
        rec = self.fetch_easy(name, sh, sub_key_list=[direct])
        if rec is None: 
            print 'rec is None for sh %s, exit'%(sh, )
            return
        #print 'rrrrr', rec
        temp = [i[1] for i in rec]
        temp = np.array(temp)[: dist_max]
        sf = fft(temp)
        
        fig = kwargs.get('fig')
        if fig is None: 
            fig = plt.figure()
        if len(fig.axes)==0: 
            ax = fig.add_subplot(111)
        else:
            for aa in fig.axes: 
                if len(aa.lines)==0: 
                    ax = aa
                    break 
            #ax = fig.axes[-1]
        ax.plot(sf)
        if kwargs.get('return_fig'): 
            return fig
        if save_fig: 
            fn = 'structure_factor-sh=%s.png'%(tuple(sh), )
            path = self.parpath +'/' +   fn
            fig.savefig(path)
        if show_fig: 
            plt.show()

    def plot_structure_factor(self, sh, direct=None, dist_max=100, show_fig=0, save_fig=1, **kwargs): 
        fft = np.fft.fft
        algorithm = self.get('algorithm')
        if algorithm == 'mps': 
            name = 'correlation' 
            direct = 'xx' if direct is None else direct
            ind = 2
        else: 
            name = 'correlation_extra'
            direct = 'zz' if direct is None else direct
            ind = 1
        rec = self.fetch_easy(name, sh, sub_key_list=[direct])
        if rec is None: 
            print 'rec is None for sh %s, exit'%(sh, )
            return
        
        temp = [i[ind] for i in rec]
        temp = np.array(temp)[: dist_max]
        sf = fft(temp)
        
        fig = self._plot(range(len(sf)), sf, **kwargs)
        if kwargs.get('return_fig'): 
            return fig
        if save_fig: 
            fn = 'structure_factor-sh=%s.png'%(tuple(sh), )
            path = self.parpath +'/' +   fn
            fig.savefig(path)
        if show_fig: 
            plt.show()

    def plot_energy_diff_vs_time_old(self, sh, system_class=None, **kwargs):
        if self.get('version') is None: 
            if system_class is None: 
                from merapy.hamiltonian import System as system_class
            sh = sh[0], sh[1]-1
            path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
            #print path; exit()
            res = system_class.load(path)
            recs = res.energy_record
        else: 
            sh = sh[0], sh[1]
            recs= self['energy'][sh]
        
        fig = kwargs.get('fig')
        if fig is None:
            fig = plt.figure()
        if len(fig.axes)==0: 
            ax = fig.add_subplot(111)
        else: 
            ax = fig.axes[-1]
        if 1:
            t_eng = [(t, v[0])  for t, v in recs.iteritems()]
            t, eng = zip(*t_eng)
            t_diff = [(t_eng[i+1][0], (t_eng[i+1][1]-t_eng[i][1])/(t_eng[i+1][0]-t_eng[i][0])) 
                    for i in range(len(t_eng)-1)]
            t, diff = zip(*t_diff)
            diff = np.array(diff)
            diff_log = np.log10(abs(diff))  # 必须在log下画，否则什么都看不清楚
        ax.plot(t, np.sign(diff), label='$sign(E_{t+1}-E_t)$')           
        ax.plot(t, diff_log, label='$log(E_{t+1}-E_t)$')
        ax.set_xlabel('t',fontsize=16)
        l=ax.legend(prop={'size':14})  #fontsize
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transparent
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
       
        
        if kwargs.get('show_fig'): 
            plt.show()
    
    def plot_energy_vs_time(self, sh, switch=None, system_class=None, **kwargs):
        """
            kwargs: xlim, 
        """
        #print 'kkk', kwargs
        self._plot_energy_vs_time(sh, switch=switch, **kwargs)

    def update_energy_temp(self, sh): 
        """
            a func,  used temporaryly
        """
        aa = ['energy', tuple(sh)]
        if not self.has_key_list(aa): 
            if system_class is None: 
                from merapy.hamiltonian import System as system_class
            if len(sh)==2: 
                sh = sh[0], sh[1]-1
            else: 
                sh = sh[0], sh[1]-1, sh[2]
            path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
            #print path; exit()
            res = system_class.load(path)
            recs = res.energy_record
            self.add_key_list(aa, val=recs, verbose=1)
            self.commit(info=1)
        
    def get_energy(self, sh, which='eng', iter_min=0, iter_max=1e6, system_class=None, **kwargs):
        aa = ['energy', tuple(sh)]
        iter_min = 0 if iter_min is None else iter_min
        iter_max = 1e6 if iter_max is None else iter_max
        if not self.has_key_list(aa): 
            if system_class is None: 
                from merapy.hamiltonian import System as system_class
            if len(sh)==2: 
                sh = sh[0], sh[1]-1
            else: 
                sh = sh[0], sh[1]-1, sh[2]
            #path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
            path = '/'.join([self.parpath, system_class.shape_to_backup_fn(*sh)])
            #print path; exit()
            res = system_class.load(path)
            recs = res.energy_record
            self.add_key_list(aa, val=recs, verbose=1)
            self.commit(info=1)
        else: 
            sh = sh[0], sh[1]
            recs = self['energy'][sh]
        if which in ['energy', 'eng']: 
            return recs
            
        elif which == 'diff': 
            pass
            temp = np.asarray( [(k, v[0]) for k, v in recs.iteritems() if iter_min< k < iter_max])
            iter = temp[: , 0].astype(int)
            eng = temp[: , 1]
            diff = eng[1: ] -eng[0: -1]
            return diff
        else: 
            raise
        
    
    def count_eng_fluctuation(self, sh, iter_min=0, iter_max=1e10, info=0, **kwargs): 
        try: 
            recs = self.get_energy(sh)
        except IOError as err: 
            if info>0: 
                print err
            return None
        temp = np.asarray( [(k, v[0]) for k, v in recs.iteritems() if iter_min< k < iter_max])
        iter = temp[: , 0].astype(int)
        eng = temp[: , 1]
        
        
        diff = eng[1: ] -eng[0: -1]
        
        ind_inc = np.where(diff>0)[0]
        
        num_pos = np.sign(diff) + 1 
        num_pos =np.count_nonzero(num_pos)
        if 0: 
            #print num_pos, len(ind_inc), type(ind_inc), ind_inc, len(ind_dec)
            print ind_inc
            print iter[ind_inc]
            print diff[ind_inc]
        mean_fluc = np.mean(diff[ind_inc[5:]])
        
        return num_pos, mean_fluc
    
    def _plot_energy_vs_time(self, sh, switch='diff',system_class=None, **kwargs):
        if 0: 
            aa = ['energy', tuple(sh)]
            if not self.has_key_list(aa): 
                if system_class is None: 
                    from merapy.hamiltonian import System as system_class
                if len(sh)==2: 
                    sh = sh[0], sh[1]-1
                else: 
                    sh = sh[0], sh[1]-1, sh[2]
                path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
                #print path; exit()
                res = system_class.load(path)
                recs = res.energy_record
                self.add_key_list(aa, val=recs, verbose=1)
                self.commit(info=1)
            else: 
                sh = sh[0], sh[1]
                recs = self['energy'][sh]
                
        recs = self.get_energy(sh)
        if 1:
            t_eng = [(t, v[0])  for t, v in recs.iteritems()]
            t, eng = zip(*t_eng)
        if 1:   
            if switch == 'diff':  
                t_diff = [(t_eng[i+1][0], (t_eng[i+1][1]-t_eng[i][1])/(t_eng[i+1][0]-t_eng[i][0])) 
                        for i in range(len(t_eng)-1)]
                t, diff = zip(*t_diff)
                diff = np.array(diff)
                diff_log = np.log10(abs(diff))  # 必须在log下画，否则什么都看不清楚
                show_fig = kwargs.get('show_fig')
                kwargs['show_fig'] = 0
                fig=self._plot(t, np.sign(diff), label='$sign(E_{t+1}-E_t)$', **kwargs)           
                kwargs['show_fig'] = show_fig
                kwargs['fig'] = fig
                self._plot(t, diff_log, label='$log(E_{t+1}-E_t)$', xlabel='t', **kwargs)
            elif switch  == 'err': 
                #eng_min = db.get_energy_min()
                eng = np.array(eng)
                eng_min = np.min(eng)
                err = eng - eng_min
                #ax.plot(t, np.log10(err) , label='$E-E_{min}$')           
                self._plot(t, np.log10(err), xlabel='t', **kwargs)
  
    def show_scaling_dim(self, sh, alpha_remap={}, alpha_shape_dic={}, **kwargs): 
        data=[]
        if 1: 
            rec= self.fetch_easy('scaling_dim', sh)
            #print 'rrrr', rec
            if rec is not None:
                data.append(('a',) + rec[0][:5] + rec[1][:3] +rec[2][:3] )
                #data.append( rec[0][:5] + rec[1][:3] +rec[2][:3] )
            else: 
                print 'rec not found for sh = %s'%(sh, )
        print tabulate(data, headers=['']+[str(i) for i in [0,0,0,0,0,1,1,1,2]], tablefmt='simple')

    def upgrade(self): 
        #field_list = ['central_charge', 'scaling_dim', 'correlation', 'correlation_extra', 'entanglement_entropy', 'entanglement_spectrum', 'magnetization', 'entanglement_brute_force_6', 'energy']
        if 0: 
            print self.__version__; exit()
            if not hasattr(self, '__version__'): 
                self.__version__ = 0.9
        if self.get('version')==1.0: 
            print 'already newest version'
            return 
        self.backup()
        new = ResultDB(parpath=self.parpath, dbname=self.dbname + '.new')
        
        if 1: 
            for sh_iter, V in self.iteritems(): 
                sh = sh_iter[: -1]
                it = sh_iter[-1]
                for field_name, val in V.iteritems(): 
                    if not new.has_key(field_name): 
                        new[field_name] = OrderedDict()
                    if not new[field_name].has_key(sh): 
                        new[field_name][sh] = OrderedDict()
                    new[field_name][sh][it] = val
            
            #pprint.pprint(new['entanglement_entropy'][12, 5])
            #print new.fetch_easy_new( 'correlation_extra', (12, 5), sub_key_list=['pm', ])
            new['version'] = 1.0
            new.commit(info=1)
            cmd = 'mv  %s %s'%(new.path, self.path)
            os.system(cmd)
            print 'db is upgraded'
    
    def time_cost(self): 
        pass
        
    def backup(self): 
        #bac = ResultDB(self.parpath, dbname=self.dbname +'.bac')
        #bac.update(self)
        #bac.commit(info=1)
        
        cmd = 'cp %s %s'%(self.path, self.path + '.bac')
        os.system(cmd)
        print 'db is backuped: %s'%cmd

class ResultDB_vmps(ResultDB): 
     def __init__(self, parpath,  **kwargs): 
        kwargs['version'] = 1.0
        kwargs.update(algorithm='mps')
        ResultDB.__init__(self, parpath,  **kwargs)
     
     def get_energy(self, sh): 
        
        try:
            path='/'.join([self.parpath, 'N=%d-D=%d.pickle'%(sh[0], sh[1])])
            S=self.load_S( path)
            #eng=S['energy']
            eng = S['mps'].energy

        except:
            raise 
            #print a, path
        return eng
  

class ResultDB_idmrg(ResultDB): 
    def __init__(self, parpath,  **kwargs): 
        kwargs['version'] = 1.0
        kwargs.update(algorithm='idmrg')
        ResultDB.__init__(self, parpath,  **kwargs)
    
    def calc_EE(self, S):
        res=None
        lam = S['Lambda']
        U, s, V = np.linalg.svd(lam)
        spect = s
        #print np.sum(s**2)
        spect2 = spect**2
        EE = -np.sum(spect2*np.log2(spect2))
        res = EE
        return res 
    
    def get_energy(self, dim): 
        
        try:
            path='/'.join([self.parpath, 'N=0-D=%d.pickle'%dim])
            S=self.load_S( path)
            eng=S['energy']

        except:
            raise 
            #print a, path
        return eng
    
    def get_correlation(self, sh, which='correlation_extra', direct='xx', 
            r_min=None,  r_max=None, period=None, fault_tolerant=1,  info=0): 
        

        which = 'correlation'  #not extra
        key_list = [direct]
        period = period if period is not None else 2            
        rec = self.fetch_easy(which, mera_shape=sh, sub_key_list=key_list, 
                fault_tolerant=fault_tolerant, info=info-1) 
        
        if rec is None: 
            return None
       
        rec = rec.items()
        #rec = map(lambda x: (x[1]-x[0], x[2]), rec) 
        
        x,y=zip(*rec[0:-1:period])
    
        if direct == 'pm':   
            y=np.array(y);  y=y*2;  
        
        x = np.array(x)
        a = 0
        b = -1
        if r_min is not None : 
            temp = np.where(x>=r_min);  a=temp[0][0]
        if r_max is not None : 
            temp = np.where(x<=r_max);  b=temp[0][-1] + 1
            x = x[a: b]
            y = y[a: b]
           
                   
        return x, y    
    
   
    def get_EE(self, dim): 
        res=self._measure(self.calc_EE, dim=dim)
        return res
   
    def _measure(self, func, dim):
        res = OrderedDict()
        try:
            path='/'.join([self.parpath, 'N=0-D=%d.pickle'%dim])
            #S=iDMRG.load( path)
            S= self.load_S(path)
            res = func.im_func(None, S)
            
        except:
            #print a, path
            raise
            
        return res


class TestResultDB(unittest.TestCase): 
    def setUp(self):
        self.seq = range(10)

        args= {}
        #parpath = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0'
        parpath = '/tmp/'
        sh = (12, 5);  
        args['dir'] = parpath
        args['mera_shape_list'] = sh
        args['draw'] = 1
        which='fit_correlation'
        args['which'] = which
        parpath = args['dir'] 
        self.db = ResultDB(parpath)
        
        parpath_idmrg = '/home/zhli/mps_backup_folder/run-heisbg-long/alpha=2.0'
    def test_temp(self): 
        pass 
        #parpath =  '/home/zhli/mps_backup_folder/run-heisbg-long/idmrg/alpha=2.0'
        parpath =  '/home/zhli/mps_backup_folder/run-ising/idmrg/h=1.0'
        
        db = ResultDB_idmrg(parpath)
        
        print db.get_correlation((0, 8))
    
    def test_idmrg_get_EE(self): 
        parpath =  '/home/zhli/mps_backup_folder/run-heisbg-long/idmrg/alpha=2.0'
        db = ResultDB_idmrg(parpath)
        
        print db.get_EE(2)
   
    def test_add_key_list(self): 
        xx = [1, 2, 3, 4]
        val = 'xxx'
        self.db.add_key_list(xx, val='xxx')
        self.assertTrue(self.db.has_key_list(xx ))
        self.assertFalse(self.db.has_key_list(xx + [5]))
        

if __name__ == '__main__': 

    if 1: 
        which_list = FIELD_NAME_LIST + ['structure_factor', 'fit_correlation']
        temp  =  ResultDB.__dict__.keys()
        which_list = filter(lambda x: 'plot' in x or 'show' in x or 'get' in x,  temp)
        which_list += ['fit_correlation'] 
        #import types 
        #temp = filter(lambda x: isinstance(temp[x], types.FunctionType) , temp.keys())
        
        parser = argparse.ArgumentParser(description='dont know')
        #parser.add_argument('-dir', '--dir', default='./')
        parser.add_argument('dir', nargs='?', default='./')
        parser.add_argument('-direct', '--direct', default='zz', choices=['pm', 'zz', 'xx', 'yy'])
        parser.add_argument('-w', '--which', choices=which_list,  default=None)
        parser.add_argument('-sw', '--switch', default=None)
        #parser.add_argument('-sh', '--mera_shape_list', nargs='*', default='all')
        parser.add_argument('-sh', '--mera_shape_list', nargs=2, type=int,  default=None)

        #parser.add_argument('-show', '--show', action='store_true', default=False) 
        parser.add_argument('-show', '--show', default=False) 
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
        args = vars(args)
        if args['mera_shape_list'] is not None : 
            args['mera_shape_list'] = tuple(args['mera_shape_list'])

        
        sh = args['mera_shape_list']
        which = args['which']
        if which is not None: 
            
            args.pop('which')
            res=db.__getattribute__(which)(sh=sh, show_fig=1, show_val=1, **args)
            print res
            #db.upgrade()
            #db.backup()

    #if len(sys.argv)<=1 :   
    if 0: 
        args= {}
        parpath = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=0.3'
        #parpath = '/tmp/'
        sh = (4, 4);  
        args['dir'] = parpath
        args['mera_shape_list'] = sh
        args['draw'] = 1
        #args['show_fig'] = 1
        args['which'] = 'plot_structure_factor' 
        args['switch'] = 'diff'
        args['fit'] = 'power'
        args['plot'] = 1
        args['plot_fit'] = 1
        #args['xlim'] = (0, 10000)
        
        parpath = args['dir'] 
        db = ResultDB(parpath)
    
        if 0: 
            sh = args['mera_shape_list']
            #sh[0] = eval(sh[0]); sh[1]=eval(sh[1])
            which = args['which']
            
            #if which in which_list: 
            if 1: 
                args.pop('which')
                res=db.__getattribute__(which)(sh=sh, show_fig=1, show_val=1, **args)
                print res
            else: 
                pass
                #db.upgrade()
                #db.backup()

    
    if len(sys.argv)<=1 :   
      
        if 0:
    
            TestResultDB.test_temp=unittest.skip("skip test_temp")(TestResultDB.test_temp) 
            unittest.main()
        else: 
            suite = unittest.TestSuite()
            add_list = [
               TestResultDB('test_temp'), 
               #TestResultDB('test_idmrg_get_EE'), 
  
            ]
            for a in add_list: 
                suite.addTest(a)
            #suite.addTest(TestResultDB('test_ising'))
            #suite.addTest(TestResultDB('test_heisbg'))
            unittest.TextTestRunner().run(suite)
           
          
                
        
        
    
