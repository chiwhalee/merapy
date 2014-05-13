#!/usr/bin/env python
#coding=utf8

"""
    issue:
        1. as put by Pfeifer 2012 in his thesis, correlator should be calced entirely through scainvar layer tensors, 
        but is he right?
        2. should scaling dim calced through corration and diagonalization of super-scaling operator yielding 
        highly coincide results?

"""

import os
import math
import numpy as np
import sys
import warnings
import matplotlib.pyplot as plt
from collections import OrderedDict
import argparse

from merapy.hamiltonian import System
from merapy.tensor_py import iTensorFactory
from merapy.decorators import timer
from merapy.diagrams.V31 import graph_ternary

def magnetization_bac(S, direction, info=0, tau=None): 
    #G_2_2 = init_mera_graph.G_2_2
    M = S.mera
    
    G_2_2 = S.G_2_2
    
    
    if 1: 
        combine_2site = S.combine_2site
        if S.symmetry == 'Z2':  combine_2site = False
        if not combine_2site:
            pauli = iTensorFactory.pauli_mat_1site(S.symmetry)
            #name = "sigma_" + direction
            o = pauli['sigma_%s'%direction]
            i = pauli['sigma_0']  #used for calc magnetization
            oi = o.direct_product(i)  #either o*i or i*o is ok
        
        else:
            pauli = iTensorFactory.pauli_mat_2site(S.symmetry)
            #issue there may need fine adjust.  +1 or  +2  which is more precise?
            name1 = "sigma_" + d1 + '1';  name2 = "sigma_" + d2  + '1'
            sigma_0 = pauli['sigma_01'] #used for calc magnetization
            o1 = pauli[name1].copy();  o2 = pauli[name2].copy()
   
    i=0
    j = 0
    
    ilayer = 0
    ilayer_p1 = ilayer + 1
   
    h2_bac = S.H_2[ilayer][0]
    H2_bac = S.H_2[ilayer_p1][0]
    
    S.H_2[ilayer][0] = oi  #replace
    name="OO"   
    S.H_2[ilayer_p1].O[j].data[:] =0.0
    
    #S.add_env_many(M, ilayer,name, G_2_2, add_tensor=S.H_2[ilayer_p1][0], branch=0, info=info-1)

    for g in sorted(G_2_2):
        G = G_2_2[g]
        order = G.contract_order
        weight = G.weight
        S.add_env(M, ilayer, G, order, name, weight, add_tensor=S.H_2[ilayer_p1][0], info=info-1)
 


    #if 1: 
        V1=[0,1,2,3]
        V2=[2,3,0,1]

        H_2, legs = S.H_2[ilayer_p1][0].contract(S.rho_2[ilayer_p1][0], V1, V2, use_buf=True)  #rho and H_2 fully contracted to a scalar
        res = H_2.data[0]
        print res

    S.H_2[ilayer][0] = h2_bac
    S.H_2[ilayer_p1][0] = H2_bac
    return res    


def get_spin_ops(S, direct, r_mod2):
    """
    params: 
        r_mod2: used when combine_2site,  to distinguish s2 = is or si
    
    """
    dir = direct
    if 1: 
        
        combine_2site = S.combine_2site
        if S.model == 'Ising': combine_2site = False 
        
        #if S.symmetry == 'Z2':  combine_2site = False
        if not combine_2site:
            pauli = iTensorFactory.pauli_mat_1site(S.symmetry)
            name = 'sigma_%s'%dir
        else:
            pauli = iTensorFactory.pauli_mat_2site(S.symmetry)
            name = 'sigma_%s%d'%(dir, r_mod2 + 1) #e.g. as for s0sr,  when r=6,  o1=si, o2=si; when r=7, o1=si, o2=is
            
        res= pauli[name].copy()
    return res

def magnetization_bac(S, direction='z', r_mod2=0,  info=0, tau=None): 
    """
        r_mod2:  is only relevent for case of combine_2site
    """
    M = S.mera
    G_2_2 = graph_ternary.G_2_2    
    if S.symmetry == 'Z2': 
        if direction in ['x', 'y']: 
            return 0.0
    
    rank = S.mera.V[0][0].rank
    assert  rank - 1  == 3 , 'only suits ternary graph'
   
    #G_2_2 = S.G_2_2
    if 0: 
        combine_2site = S.combine_2site
        if S.symmetry == 'Z2':  combine_2site = False
        if not combine_2site:
            pauli = iTensorFactory.pauli_mat_1site(S.symmetry)
            #name = "sigma_" + direction
            o = pauli['sigma_%s'%direction]
            i = pauli['sigma_0']  #used for calc magnetization
            oi = o.direct_product(i)  #either o*i or i*o is ok, they are exactly same as reflect symm is preserved
        
        else:
            pauli = iTensorFactory.pauli_mat_2site(S.symmetry)
            sigma_z1 = pauli['sigma_z1']; sigma_z2 = pauli['sigma_z2']
            #o = 0.5*(sigma_z1 + sigma_z2)  # this is incorrect !! no average 
            o = sigma_z1 
            #o = sigma_z2 
            i = pauli['sigma_01']  #used for calc magnetization
            oi = o.direct_product(i)
    
    O = get_spin_ops(S, direction, r_mod2=r_mod2)
    I = get_spin_ops(S, '0', r_mod2=r_mod2)       
    oi = O.direct_product(I)
    
    i=0
    j = 0
    
    ilayer = 0
    ilayer_p1 = ilayer + 1
   
    h2_bac = S.H_2[ilayer][0]
    H2_bac = S.H_2[ilayer_p1][0]
    
    S.H_2[ilayer][0] = oi  #replace
    name="OO"   
    S.H_2[ilayer_p1].O[j].data[:] =0.0
    
    #S.add_env_many(M, ilayer,name, G_2_2, add_tensor=S.H_2[ilayer_p1][0], branch=0, info=info-1)
    
    res= 0.0
    xx = np.ndarray(3)
    i = 0
    res= {}
    for g in sorted(G_2_2):
        G = G_2_2[g]
        order = G.contract_order
        weight = G.weight
        branch = 0
        temp = S.contract(M, ilayer, G, nOps=0, Names=[], Ops=[], order=order, exception=None, 
                branch=branch, info=info-1)
        #temp.data *=  weight
        #xx.append(temp.data[0])
        #res += temp.data[0]
        #xx[i] = temp.data[0]
        res[g] = temp.data[0]
    if 0:     
        fluc = xx - xx.mean()
        fluc=fluc.tolist()
        print 'fluctuation of <sigma> at different site %1.5e %1.5e %1.5e'%(tuple(fluc))
        res = np.sum(xx)
        
    S.H_2[ilayer][0] = h2_bac
    S.H_2[ilayer_p1][0] = H2_bac
    return res    

def __magnetization(S, direction, r,   info=0, tau=None): 
    """
        r_mod2:  is only relevent for case of combine_2site
    """
    
    G_2_2 = graph_ternary.G_2_2    
    if S.symmetry == 'Z2': 
        if direction in ['x', 'y']: 
            return 0.0
    
    rank = S.mera.V[0][0].rank
    assert  rank - 1  == 3 , 'only suits ternary graph'
   
    r_mod2 = r%2
    O = get_spin_ops(S, direction, r_mod2=r_mod2)
    I = get_spin_ops(S, '0', r_mod2=r_mod2)       
    oi = O.direct_product(I)
    
    r_mod3 = r%3
    #graph_key = {0: (0, 1), 1: (1, 2), 2: ()}
    graph_key = (r_mod3, r_mod3 + 1)
    G = G_2_2[graph_key]
    res = S.contract_new(layer=0, G=G, exception=None, extra_extra={'oo':oi}, info=info-1)
    res= res.data[0]
    
    return res    

def magnetization(S, direction='z', r_list=None, info=0, **kwargs): 
    if 0: 
        if S.symmetry == 'Z2': 
            direction = 'x'
        else: 
            direction = 'z'
    
    res= OrderedDict()
    if r_list is None: 
        #if S.symmetry == 'U1' : 
        if  S.combine_2site  and S.model != 'Ising' :   #issue:  这里有个问题，算ising模型时，combine_2site = True, 这样，但实际上为False !
            #print 'sss', S.combine_2site, S.symmetry
            r_list = range(6)
        else: 
            r_list = range(3)
        
    for r in r_list: 
        res[r] = __magnetization(S, direction, r)
    return res


if __name__ == '__main__':
    """
        0.212066211149
        0.424117854237
        0.636171429421
        0.636171429421
    """
    
    
    #path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=1.5/8.pickle'
    if 1: 
        parser = argparse.ArgumentParser(description='dont know')
        parser.add_argument('path', nargs='*', default=None)
        parser.add_argument('-d', '--direction', default=None)
        parser.add_argument('-dim', '--dim_list', type=int, nargs='*', default=None)
        args = parser.parse_args()
        #print args
        #print vars(args)
        args= vars(args)
        if len(args['path'])>0:  
            args['path']=args['path'][0]
            S= System.load(args['path'])
            print magnetization(S, 'z')
           
        else: 
            #p1 = '/home/zhli/Documents/mera_backup_tensor/run-ising/ternary/z2-symm/scale-invar/h=1.0/4.pickle'
            p1 = '/home/zhli/Documents/mera_backup_tensor/run-wigner-crystal/a=2.0-h=1.4142135-V=0.0/4.pickle'
               
            args['path']=p1
            S= System.load(args['path'])
            res=magnetization(S, 'z')
            print res
            assert res[0] == 0.63641726101318663 

   

    
    
