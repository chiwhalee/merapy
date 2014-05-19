#!/usr/bin/env python
#coding=utf8

"""
    issue:
        1. as put by Pfeifer 2012 in his thesis, correlator should be calced entirely through scainvar layer tensors, 
        but is he right?
        2. should scaling dim calced through corration and diagonalization of super-scaling operator yielding 
        highly coincide results?
        3. the way calc correlation has not average translation of 3 site, this not good?

"""

import os, sys
import math
import numpy as np
import pprint
import warnings
import matplotlib.pyplot as plt
from collections import OrderedDict
import argparse
from scipy.optimize import curve_fit

from result_db import ResultDB

from merapy.hamiltonian import System
from merapy.tensor_py import iTensorFactory
from merapy.decorators import timer
from merapy.graphics import calc_ascending_site
from merapy.diagrams.V31 import graph_ternary
from merapy.descending import descending_ham
  
from magnetization import __magnetization  


def ascending_middle_site(S, o1, ilayer):
    M = S.mera
    V = M.V[ilayer][0]
    V_dag = M.V_dag[ilayer][0]
    o1.type_name = 'o'
    lab1 = [0, 1, 2, 3]
    lab2 = [6, 1]
    lab3 = [4, 0, 6, 2]
    #print V
    #exit()
    V_o1, lab12 = V.contract(o1, lab1, lab2)
    #V_o1_Vdag, lab123 = V_o1.contract(V_dag, lab12, lab3)
    res, lab123 = V_o1.contract(V_dag, lab12, lab3)
    res=res.permutation([1, 0])
    #print lab123 
    #exit()
    return res


def fetch_correlation(S, direction, info=0): 
    if 1: 
        is_found = False
        record = S.record
        try:
            #print 'try to fetch correlation from S.record...  ', 
            msg = 'try to fetch correlation_%s from S.record...  '%direction
            res = record[S.iter]['correlation'][direction]#['val']
            val = res['val']
        except KeyError as err:
            msg += 'record not found:  ' + str(err.message)
        else:
            msg  += 'record found'  
            is_found = True
        
        if info>0: print msg
    if is_found: 
        res= {direction: res}
    else: 
        res= {}
    return res


def get_spin_ops_2(S, direct, r_mod2):
    """
    s1, s2
    params: 
        r_mod2: used when combine_2site,  to distinguish s2 = is or si
    
    """
    direction = direct
    if 1: 
        combine_2site = S.combine_2site
        #if S.symmetry in ['Z2', 'Travial']:  combine_2site = False
        if S.symmetry in [ 'Travial']:  combine_2site = False
        d1, d2 = direction[0:2]
        if not combine_2site:
            pauli = iTensorFactory.pauli_mat_1site(S.symmetry)
            name1 = "sigma_" + direction[0];  name2 = "sigma_" + direction[1]   
            sigma_0 = pauli['sigma_0']  #used for calc magnetization
        else:
            pauli = iTensorFactory.pauli_mat_2site(S.symmetry)
            #issue there may need fine adjust.  +1 or  +2  which is more precise?
            #name1 = "sigma_" + d1 + '1';  name2 = "sigma_" + d2  + '1'
            name1 = 'sigma_%s1'%d1; 
            #print 'rrrr', r_mod2
            name2 = 'sigma_%s%d'%(d2, r_mod2 + 1) #e.g. as for s0sr,  when r=6,  o1=si, o2=si; when r=7, o1=si, o2=is
            sigma_0 = pauli['sigma_01'] #used for calc magnetization
        o1 = pauli[name1].copy();  o2 = pauli[name2].copy()
    return o1, o2, sigma_0


#@timer
def __correlation(distance_max, direction, S=None, path=None, force_update=False, adjust_z=True, plot=False, verbose=True):
    """
        lay is the numbering of layers of mera. under this convension, lay =0  refers to the first layer
        refs:
            Evenbly and Vidal 2007, 
    """
        
    is_found = False
    name = 'correlation_%s'%direction
    SIlayer = S.mera.num_of_layer-2
    
    o1, o2, sigma_0 = get_spin_ops_2(S, direction, r_mod2=0)  #r_mod2 i.e. assume always r even
    
    if 1: 
        val = []
        rescale_factor = S.mera.period
        lay_max = int(math.log(distance_max, rescale_factor))
        
        assert rescale_factor == 3 
        lab1 = [0, 1, 2, 3]
        warnings.warn("note this maybe wrong 01, 23 which are upward?")
        #issue: why they are all right?
        #lab2 = [2, 0];  lab3 = [3, 1]
        lab2 = [0, 2];  lab3 = [1, 3]
        
        if (direction == 'zz') and adjust_z: #calc <sigma_z>
            sz = o1.copy()
            sz = ascending_middle_site(S, sz, ilayer=0)
            
            s0 = sigma_0.copy()
            s0 = ascending_middle_site(S, s0, ilayer=0)
            
            rho2_o1, lab12 = S.rho_2[1][0].contract(sz, lab1, lab2, info=0)
            rho2_o1_o2, lab123 = rho2_o1.contract(s0, lab12, lab3)
            magnetization = rho2_o1_o2.data[0]     
            #print 'magnetization <sigma_z> = %1.10f'%magnetization #, magnetization-2/np.pi
            #magnetization = 0.0
            
            #magnetization = calc_magnetization(S, 'z'); print 'magnetization_test  = ', magnetization
        
        else: 
            magnetization = 0.0
       
        for lay in range(lay_max):
            distance = rescale_factor**(lay + 1)
            if S.combine_2site: distance *= 2
            lay = lay if lay <= SIlayer else SIlayer
            lay_p1 = lay + 1
            o1=ascending_middle_site(S, o1, lay)
            o2=ascending_middle_site(S, o2, lay)
        
            which = lay_p1 if lay<= SIlayer else SIlayer + 1
            rho2_o1, lab12 = S.rho_2[which][0].contract(o1, lab1, lab2, info=0)
            rho2_o1_o2, lab123 = rho2_o1.contract(o2, lab12, lab3)
            
            val.append((distance, rho2_o1_o2.data[0] - magnetization**2))
        #print "may be a problem\n odd even"*10
        warnings.warn("issue: this may be a problem. in some cases <sz_i> neq <sz_j> , even there sign are different!")

    if 1:
        x = np.array([i[0] for i in val])
        x = np.log10(x)
        y = np.array([i[1] for i in val])
        y = np.log10(y)
        k, b = np.polyfit(x, y, 1)
        nn = SIlayer
        k_first, b_first = np.polyfit(x[:nn], y[:nn], 1)   #only use transitiaon-layer tensors
        k_last, b_last = np.polyfit(x[nn+1:], y[nn+1:], 1) #use more than one SI tensor
    
    if 1: 
        res = {}
        res[direction] = {}
        k = round(k, 6) 
        res[direction]['exponent'] = k
        res[direction]['exponent_first'] = k_first
        res[direction]['exponent_last'] = k_last
        res[direction]['val'] = val

    return res

def eta_fit(x, y, a, b): 
    eta, c= np.polyfit(np.log10(x[a:b]), np.log10(y[a:b]), 1) 
    return eta


#def correlation_bac(S=None, path=None, distance_max=10000, direction=None, dim_list=None, 
    #        force_update=False, adjust_z=True,
    #        plot=False, verbose=True, **kwargs):
    #    """
    #        lay is the numbering of layers of mera. under this convension, lay =0  refers to the first layer
    #        refs:
    #            Evenbly and Vidal 2007, 
    #    """
    #    is_found = False
    #    name = 'correlation_%s'%direction
    #    
    #    path = os.path.abspath(path)
    #    parpath = os.path.dirname(path) 
    #    #print parpath; exit()
    #    
    #    if S is None: 
    #        if path is None: raise Exception('S, path shouldnt be none at the same time')
    #        S = System.load(path)
    #        is_loaded = True
    #        
    #    if direction is None: 
    #        if S.symmetry in ["Z2", 'Travial']: 
    #            direction = ["xx", "yy", "zz"]
    #        if S.symmetry == "U1": 
    #            direction = ["zz", "pm"]
    #    else: 
    #        direction = [direction]
    #    
    #    res= {}
    #    is_update = False
    #    for d in direction: 
    #        temp = fetch_correlation(S, direction=d)
    #        #print temp; exit()
    #        if temp == {} or force_update: 
    #            is_update = True
    #            temp = __correlation(distance_max=distance_max, direction=d, S=S, force_update=force_update, adjust_z=adjust_z)
    #        res.update(temp)
    #    
    #    rdb = ResultDB(path)
    #    propt = S.key_property()
    #    num_of_layer = propt['num_of_layer']
    #    trunc_dim = propt['trunc_dim']
    #
    #    if is_update  and (path is not None): 
    #        if not S.record.has_key(S.iter):  
    #            S.record[S.iter]={}
    #        if S.record[S.iter].has_key('correlation'): 
    #            S.record[S.iter]['correlation'].update(res)
    #        else: 
    #            S.record[S.iter]['correlation'] = res
    #           
    #        if verbose: print 'correlation at step %d is saved'%S.iter
    #        S._save(path)
    #        
    #        if not rdb.has_key((trunc_dim, num_of_layer, S.iter)): 
    #            rdb[(trunc_dim, num_of_layer, S.iter)]={}
    #        if rdb[(trunc_dim, num_of_layer, S.iter)].has_key('correlation'): 
    #            rdb[(trunc_dim, num_of_layer, S.iter)]['correlation'].update(res)
    #        else: 
    #            rdb[trunc_dim, num_of_layer, S.iter]['correlation'] = res
    #        if verbose: print 'correlation is saved in %s'%rdb.DBNAME
    #        rdb.commit()
    #   
    #    if plot:
    #        if dim_list is None: 
    #            dim_list = [trunc_dim]
    #        
    #        #record = []
    #        record = {}
    #        dim_list_not_exist = []
    #        for i in dim_list: 
    #            
    #            #record.append(rdb.fetch(field_name='correlation', dim=i, num_of_layer=None))
    #            record[i]=rdb.fetch(field_name='correlation', dim=i, num_of_layer=None)
    #            if record[i] is None: 
    #                dim_list_not_exist.append(i)
    #        if len(dim_list_not_exist)>0: 
    #            msg = 'these dims not found: %s  , probably pickle files for these dims not exsit'%dim_list_not_exist
    #            warnings.warn(msg)
    #            dim_list = [d for d in dim_list if  not d in dim_list_not_exist]
    #            
    #            
    #        if rdb.has_key((0, 0, 0)): 
    #            print rdb[0, 0, 0]
    #            record['exact'] = rdb[0, 0, 0]['correlation']
    #        else: 
    #            record['exact'] = None
    #        plot_correlation(record, direct_list=direction, dim_list=dim_list, save_parpath=parpath)
    #    
    #    #res= S.record[S.iter]['correlation']
    #    #print 'rrr', res
    #    return res
    #

def correlation(S, distance_max=10000, direction=None, dim_list=None, 
        force_update=False, adjust_z=True,
        plot=False, verbose=True, **kwargs):
    """
        lay is the numbering of layers of mera. under this convension, lay =0  refers to the first layer
        refs:
            Evenbly and Vidal 2007, 
    """
    is_found = False
    name = 'correlation_%s'%direction
    
    if S is None: 
        if path is None: raise Exception('S, path shouldnt be none at the same time')
        S = System.load(path)
        is_loaded = True
        
    if direction is None: 
        if S.symmetry in ["Z2", 'Travial']: 
            direction = ["xx", "yy", "zz"]
        if S.symmetry == "U1": 
            direction = ["zz", "pm"]
    else: 
        direction = [direction]
    
    res= {}
    is_update = False
    for d in direction: 
        temp = __correlation(distance_max=distance_max, direction=d, S=S, force_update=force_update, adjust_z=adjust_z)
        res.update(temp)
    
    if plot and 0:
        if dim_list is None: 
            dim_list = [trunc_dim]
        
        #record = []
        record = {}
        dim_list_not_exist = []
        for i in dim_list: 
            
            #record.append(rdb.fetch(field_name='correlation', dim=i, num_of_layer=None))
            record[i]=rdb.fetch(field_name='correlation', dim=i, num_of_layer=None)
            if record[i] is None: 
                dim_list_not_exist.append(i)
        if len(dim_list_not_exist)>0: 
            msg = 'these dims not found: %s  , probably pickle files for these dims not exsit'%dim_list_not_exist
            warnings.warn(msg)
            dim_list = [d for d in dim_list if  not d in dim_list_not_exist]
            
            
        if rdb.has_key((0, 0, 0)): 
            print rdb[0, 0, 0]
            record['exact'] = rdb[0, 0, 0]['correlation']
        else: 
            record['exact'] = None
        plot_correlation(record, direct_list=direction, dim_list=dim_list, save_parpath=parpath)
    
    return res

def __correlation_extra(S, r, direct, info=0):
    """
        two point correlation
        direct: either be 'pm', 'xx', 'x0', etc
        correlation for more s0sr cases
    """
    
    #print 'dont forget magnetization'
    #n_min = 2 if S.only_NN else 3
    if S.combine_2site and S.symmetry in ['U1', 'Z2'] : 
        r_mod2 = r%2
        r = r//2
    else: 
        r_mod2 = None
    site_list = calc_ascending_site(r, n_min=1)
    if info>0: 
        print 'o1-o2 at each coarsed layer %s'%(site_list, )
    temp = [(s[0]%3, s[1]%3) for s in site_list]
    #print temp
    
    _o1, _o2, sigma_0 = get_spin_ops_2(S, direct, r_mod2)
    o1 = _o1.copy(); o2= _o2.copy()
    o2 = sigma_0.direct_product(o2)
        
        
    numl = len(site_list)  #number of coarse grained o1-o2 entil dist(o1, o2)<=n_min
    num_of_layer = S.mera.num_of_layer
    r_max = 3**(num_of_layer-1)
    
    #assert numl <= num_of_layer, 'r is too large, num_of_layer=%d, r_max=%d'%(num_of_layer, r_max)
    assert 1 <= r <= r_max, 'r exceed range r=%d,  r_max=%d'%(r, r_max)
    
    #o1, o2分别独立ascending to higher layer finally becomes either (0, (2, 3)) or (0, (1, 2))
    for site in site_list[:-2]: 
        l = site_list.index(site)
        if info>0 : print 'oo at layer=%d,  site=%s'%(l, site_list[l], )
        
        site = site_list[l][1]
        site_mod3 = site%3
        graph_key = {0: (2, 3), 1: (0, 1), 2: (1, 2)}[site_mod3]
        
        G2 = graph_ternary.G_2_2[graph_key]; ee = G2.find_node('OO')[0]
        o2 = S.contract_new(layer=l, G=G2, exception=ee, extra_extra={'oo':o2 }, info=info-1)
        G1 = graph_ternary.G_1_1[0]; ee = G1.find_node('O')[0]
        o1 = S.contract_new(layer=l, G=G1, exception=ee, extra_extra={'o':o1, 'O': None}, info=info-1)
    
    #o1, o2 ascengding 并且合to top layer (0, (2, 3)) or (0, (1, 2))->(0, 1)
    for site in site_list[-2: -1]: 
        l = site_list.index(site)
        #print 'ssss', site, site_list
        if info>0:  print 'site is %s'%(site, )
        if  site == (1, 2): 
            G = graph_ternary.G_3_2[(0, 1, 2)]; ee=G.find_node('OO')[0]
            ooo = o1.direct_product(o2)
            oo_final = S.contract_new(layer=l, G=G, exception=ee, extra_extra={'ooo': ooo})
        elif site  == (2, 3): 
            G = graph_ternary.G_22_2[(0, 3)]
            id = iTensorFactory.identity(qsp=[q.copy() for q in  o1.QSp[: 2]])
            o1 = o1.direct_product(id); ee=G.find_node('OO')[0]
            oo_final = S.contract_new(layer=l, G=G, exception=ee, extra_extra={'oo1':o1, 'oo2':o2})
        else: 
            raise
    
    #at top layer o1o2和rho收缩
    l = numl-1
    if l == 0: #i.e r == 1 
        oo_final  = _o1.direct_product(_o2)
        descending_ham(S.mera, S, ilayer=1)
    rho = S.rho_2[l][0]  
    res = System.contract_oo_rho_2(oo_final, rho)
    #if direct == 'zz': 
    if direct in ['zz', 'xx', 'yy']: 
        mag_o1 = __magnetization(S, direct[0], r)
        mag_o2 = __magnetization(S, direct[1], r)
        res= res-mag_o1*mag_o2
    return res

def correlation_extra(S, r_list=None, direct_list=None, fail_skip=1, **kwargs): 
    num_of_layer = S.mera.num_of_layer
    if num_of_layer >= 7: 
        #return None
        raise Exception('too many layers, abort it ')
    
    r_min = 1
    r_max = 3**(num_of_layer-1)
    if S.combine_2site and S.symmetry in ['U1', 'Z2']: 
        r_min *= 2
        r_max *=  2
        
    if r_list is None: 
        r_list = range(r_min, r_max)
        
    if direct_list is None: 
        if S.symmetry in ["Z2", 'Travial']: 
            direct_list = ["xx", "yy", "zz"]
        if S.symmetry == "U1": 
            direct_list = ["zz", "pm"]
    res= {} 
    
    for dir in direct_list: 
        val = []
        r_failed = []
        for r in r_list: 
            try: 
                corr = __correlation_extra(S=S, r=r, direct=dir)
                #print 'vvv', r, corr
                val.append((r, corr))
            except Exception as err: 
                
                if not fail_skip: 
                    #raise(Exception)
                    raise
                r_failed.append((r, err))
        if len(r_failed)>0: 
            print dir, 'r_failed is %s'%(r_failed, ) 
                
        res[dir] = val
    return res

def plot_correlation(record, direct_list=None, dim_list=None, save_parpath=None): 
    pass
    #rec = S.rec[S.iter]['correlation']
    rec = record[dim_list[0]]
    direct_list = direct_list if direct_list is not None else rec.keys()
    #print direct_list; exit()
    
    
    #for i in range(len(dim_list)): 
    
    fig = plt.figure( figsize=(5*len(direct_list), 4))
    ax = {}
    for direct in direct_list: 
        dir = direct
        ax[dir] = fig.add_subplot('1%d%d'%(len(direct_list), direct_list.index(dir) + 1))
        legend = []
        if record['exact'] is not None: 
            pass
            val = record['exact'][direct]
            k = record['exact'][direct]['exponent']
            x, y = zip(*val['val'])
            #plt.loglog(x, y, '--')
            ax[dir].loglog(x, y, '--')
            _direct = direct if direct != 'pm' else 'xx' 
            #legend.append('exact %s $\eta$=%s'%(direct, k))
            legend.append('exact $\eta$=%s'%(k))
        for i in dim_list: 
            #chi = dim_list[i]
            chi = i
            rec = record[i]
            #direct_list = direct_list if direct_list is not None else rec.keys()
        
            val = rec[direct]
            x, y = zip(*val['val'])
            x = np.array(x);  y=np.array(y)
            if direct == 'pm': 
                #direct = 'xx'
                msg = 'direction pm is changed to xx, and Cxx = 2.0*Cpm'
                warnings.warn(msg)
                y = y*2.
            
            k, k_first, k_last = val['exponent'], val['exponent_first'], val['exponent_last']
            k_middle, b= np.polyfit(np.log10(x[2:5]), np.log10(y[2:5]), 1) 
            
            #plt.loglog(x, y, '.-')
            ax[dir].loglog(x, y, '.-')
        
            #legend.append('%s $\chi=%d$ $\eta$=%1.2f (%1.2f, %1.2f)'%(direct, chi, k, k_first, k_last))
            #legend.append('%s $\chi=%d$ $\eta$=%1.2f (%1.2f)'%(direct, chi, k, k_first))
            #if direct == 'pm':  _direct = 'xx' 
            _direct = direct if direct != 'pm' else 'xx' 
            #legend.append('%s $\chi=%d$ $\eta$=%1.2f'%(_direct, chi, k_middle))
            legend.append('$\chi=%d$ $\eta$=%1.2f'%(chi, k_middle))
            ax[dir].legend(legend, loc='lower left', prop={'size': 11}, title=_direct)
            ax[dir].set_xlabel('$r$', fontsize=18)
        ax[dir].grid(1) 
    #plt.legend(legend)
    #plt.legend(legend, loc='lower left', bbox_to_anchor=(1.1, 1.1))
    if 0: 
        plt.legend(legend, loc='lower left')
        plt.xlabel('$r$', fontsize=18)
        plt.ylabel('correlation', fontsize=18)
    
    #trunc_dim = S.key_property()['trunc_dim']
    dim = str(dim_list[0])
    for i in dim_list[1: ]: 
        dim += '_' + str(i)
    dir = str(direct_list[0])
    for i in direct_list[1: ]:
        dir += '_' + i 
    fn = 'correlation-dim=%s-direct=%s.png'%(dim, dir)
    if save_parpath is not None: 
        fn = save_parpath + '/'  + fn
    plt.savefig(fn, bbox_inches='tight') 
    #plt.set_title('$chi=$trunc_dim')
    #txt = "$k=%s$"%k
    #plt.text(.1, .1, txt )
    #plt.show()
    os.popen("display " + fn) 


if __name__ == '__main__':
    warnings.filterwarnings('once')
    
    #if 1: 
    if  len(sys.argv)>1: 
        parser = argparse.ArgumentParser(description='dont know')
        parser.add_argument('path', nargs='*', default=None)
        parser.add_argument('-d', '--direction', default=None)
        parser.add_argument('-f', '--force_update', action='store_true')
        parser.add_argument('-sr', '--show_result', default=False)
        parser.add_argument('-az', '--adjust_z', type=int, default=1)
        parser.add_argument('-dim', '--dim_list', type=int, nargs='*', default=None)
        args = parser.parse_args()
        #print args
        #print vars(args)
        args= vars(args)
        
        if len(args['path'])>0:  
            args['path']=args['path'][0]
            S= System.load(args['path'])
            
            res=correlation(S=S, distance=10000,  plot=True, **args)
            print res
    else: 
        args= {}
        
        p1 = '/home/zhli/Documents/mera_backup_tensor/run-ising/ternary/z2-symm/scale-invar/h=1.0/4.pickle'
        p1 = '/home/zhli/Documents/mera_backup_tensor/run-wigner-crystal/a=2.0-h=1.4142135-V=0.0/4.pickle'
        args['path']=p1
        args['force_update'] = 1
        S = System.load(args['path'])
        #res=correlation(S=S, distance=10000, plot=True, **args); print res
        #S = System.load(p1)
        res=correlation_extra(S, r_list=[27], fail_skip=0 ); print res
        
            


    
    
