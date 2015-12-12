#!/usr/bin/env python
#coding=utf8

#import numpy.linalg as linalg
#import scipy
import sys
import numpy as np
import cPickle as pickle
from collections import OrderedDict
import matplotlib.pyplot as plt
import unittest
import tempfile 

import merapy 
from merapy.utilities import print_vars, dict_to_object 
from merapy import load, iTensorFactory
from merapy.hamiltonian import System
from merapy.tensor_svd import Tensor_svd
from merapy.quantum_number_py import *
from merapy.decorators import timer
from merapy.descending import descending_ham
from merapy.tensor_network import TensorNetwork
from merapy.diagrams.V31.graph_ternary import G_calc_EE_gen
from merapy.graphics import Graphics

__all__ = ['entanglement_entropy', 'entanglement_spectrum', 'central_charge', 
    'entanglement_extra', 
    'entanglement_brute_force', 'entanglement_brute_force_aver',
    'entanglement_brute_force_9', 'entanglement_brute_force_6', 
    'entanglement_brute_force_9_aver', 'entanglement_brute_force_6_aver' , 
    'entanglement_special', 'concurrence'
    ]

"""

    todo: 
        根本的测试，应该用mera构造一个比如GHZ态，来准确测试！ 
        应该从单个(2, 2) 型张量的稀疏对角化开始来测试
        
    issue: 
        when calc ent entropy,  only retain k = 40 eig values of rho,  this may be not enough!!!
        num of quantum num is not enough when eig_sparse
        Z2 symm is not properly treated
"""
def initialize_S(S): 
    from merapy.hamiltonian import System 
    if not isinstance(S, dict) and not isinstance(S, System): 
        S= S.__dict__ 
    if isinstance(S, dict): 
        S = System(symmetry='U1', init_state=S)
        S.info = 2
        System.set_graph.im_func(S, getattr(S, 'graph_module_name', None))
        S.initialize()
    return S

def rho_eig(rho, k=None, info=0): 
    """
        k: number of eigen vals during sparse diag 
    """
    symmetry = rho.QSp[0].QnClass.SYMMETRY
    qn_list = symmetry_to_Qn(symmetry).QNS
    if symmetry == "U1":
        #issue: is this list large enough?
        qn_list = [0, -1, 1, -2, 2, -3, 3, -4, 4]
    qn_instance_list = [qn_factory(symmetry, ii) for ii in qn_list]
    
    v = {}
    
    if k is None: 
        k = 80
    #v = Tensor_svd.eig_sparse(rho, qn_instance_list, k=k)
    v = Tensor_svd.eig_sparse(rho, qn_instance_list, k=k, is_hermite=1)
        

    #print 'eigen val ', v
    for q in v:
        if qn_list.index(q)== 0:
            val = v[q]
        else:
            
            val = np.append(val, v[q])
    sum_val = np.sum(val)           
    if 1.0-sum_val>1e-8: 
        msg = 'error val should sum to 1.0, sum_val = %f'%(sum_val, )
        #raise Exception(msg)  # 1e-8 is to strict for some cases, so disable check for now
    val = val[val != 0.0]    

    if 1: 
        val_neg = val[val < 0.0]

        val.sort()
        val = val[::-1]
        if info>0:
            print val[:10], '\n', val[-5:]
            print val_neg
        
        EE = -np.sum(val*np.log2(val))
        
        try:
            EE = -np.sum(val*np.log2(val))
        except:
            raise ValueError(" neg %s"%(val2[val2<0.0]))
   
    return EE, val

def partial_trace_rho(rho, n, average=0):
    """
        n: num of sites to be traced out 
    """
    if n == 0: 
        return rho
    
    #res= rho.partial_trace(0) + rho.partial_trace(1)
    site_total = range(rho.rank//2)
    #site_left = range(n)
    #site_right = range(site_total-n, site_total)
    site_left = site_total[: n]
    site_right = site_total[-n: ]
    if average: 
        res = rho.partial_trace(site_left) + rho.partial_trace(site_right)  #average between left and right
        res.data *= 0.5
    else: 
        res = rho.partial_trace(site_right)
    return res

def concurrence(S, Jzz=None): 
    """
        below only suitble for xxz model  
        the exact value for the xxx point is 2ln(2)-1 almost equal 0.386
        ref: 
            Gu sj 2003 entanglement, pra
    """
    
    #rho at 0th layer is not updated during optimization, so update it here
    from merapy.hamiltonian import System 
    if not isinstance(S, dict): 
        S= S.__dict__ 
    if isinstance(S, dict): 
        S = System(symmetry='U1', init_state=S)
        S.info = 2
        System.set_graph.im_func(S, getattr(S, 'graph_module_name', None))
        S.initialize()
    
    descending_ham(S.mera, S, ilayer=1)
    rho_01 = S.rho_2[0][0].copy()
    rho_0 = rho_01.partial_trace([1])
    pauli = iTensorFactory.pauli_mat_2site(S.symmetry)
    szz = pauli['szz']
    #the correaltion function, C^zz_{01}
    C01, _= szz.contract(rho_0, ['d', 'u'], ['u', 'd'])
    C01 = C01.data[0]
    
    eng = S.energy 
    if Jzz is None: 
        Jzz = S.model_param['Jzz']
    
    res= 0.5*max(0.0, abs(eng-Jzz*C01)-C01-1.0)
    return res 
    

def __entanglement_entropy(S, layer, k=None, spectrum=0, info=0): 
    """
        1site, 2site EE at effective layer 
        parms: 
            k: num of eigenvalues of rho to keep
    """
    if 1: 
        if layer == 0:   #rho at 0th layer is not updated during optimization, so update it here
            descending_ham(S.mera, S, ilayer=1)
        rho_2 = S.rho_2[layer][0].copy()
        if 0:         
            rho_1 = rho_2.partial_trace(0) + rho_2.partial_trace(1)
            rho_1.data *= 0.5
        rho_1 = partial_trace_rho(rho_2, n=1)
        #if info>0: print 'trace rho_2=%s, trace rho_1=%s'%(rho_2.trace(), rho_1.trace())
        
        rho = {}
        rho[2] = rho_2
        rho[1] = rho_1
        
        val = {}
        val_neg = {}
        EE = {}

        for i in [1, 2]:
            rho_ = rho[i]
            EE[i], val[i] = rho_eig(rho_, k)
        
        #cc = 3*(EE[2]-EE[1]); cc = round(cc, digit)
        digit = 6
        EE[1] = round(EE[1], digit)
        EE[2] = round(EE[2], digit)
               
        if spectrum: 
            return EE, val
        else: 
            return EE

xxxx = __entanglement_entropy 


def entanglement_spectrum(S, k=None, info=0, path=None, verbose=False): 
    numl = S.mera.num_of_layer
    res = OrderedDict()
    for l in range(numl-1): 
       nothing, res[l] = __entanglement_entropy(S, layer=l, spectrum=1)
    return res

def entanglement_entropy(S, spectrum=0, info=0, path=None, verbose=False):
    numl = S.mera.num_of_layer
    res = OrderedDict()
    for l in range(numl-1): 
       res[l]=__entanglement_entropy(S, layer=l, spectrum=spectrum)
    return res

def entanglement_special(S, ): 
    """
        this one is used in U1 symm when combine 2site, 
        rho = rho_{i1 i2 i3 i4}
        this func trace out [i1, i4]
        for calc the valence bond effect on entanglement in AFM 
    """
    from merapy.hamiltonian import System 
    if not isinstance(S, dict): 
        S= S.__dict__ 
    if isinstance(S, dict): 
        S = System(symmetry='U1', init_state=S)
        S.info = 2
        System.set_graph.im_func(S, getattr(S, 'graph_module_name', None))
        S.initialize()
  
    layer = 0 
    if layer == 0 :   #rho at 0th layer is not updated during optimization, so update it here
        descending_ham(S.mera, S, ilayer=1)
    H2 = S.H_2[0][0]
    rho_2 = S.rho_2[layer][0].copy()
    d = rho_2.qsp_class.easy_init([1, -1], [1, 1]) 
    
    def re_index(t): 
        d = t.qsp_class.easy_init([1, -1], [1, 1]) 
        d2 = d**2 
        qsp = d2.copy_many(4, reverse=[2, 3])
        res = t.__class__(QSp=qsp)
        
        for i in range(t.nidx): 
            qn_id_tuple = t.Addr_idx[:, i]
            temp = (qn_id_tuple + 1)%3 
            #print_vars(vars(),  ['i ','qn_id_tuple', 'temp'], head='', sep=', ')
            data = t[qn_id_tuple][0]
            j=res.ravel_qn_id_tuple(temp)
            j = res.idx[j] 
            res.set_block(j, data)
        return res 
    
    rho4 = re_index(rho_2)
    dc = d.conj()
    rho8 = rho4.split_qsp(0, (d, d), 1, (d, d), (2), (dc, dc), 3, (dc, dc))
    res= {}
    if 1: 
        rho6 = rho8.merge_qsp((1, 2), (5, 6))
        rho_red = rho6.partial_trace([0, 2])
        ee, _ = rho_eig(rho_red)
        res[12] = ee
    if 1: 
        rho6 = rho8.merge_qsp((2, 3), (6, 7))
        rho_red = rho6.partial_trace([0, 1])
        ee, _ = rho_eig(rho_red)
        res[23] = ee
    if 1: 
        #rho6 = rho8.merge_qsp((0, 1), (4, 5))
        #rho_red = rho6.partial_trace([1, 2])
        rho6 = rho8.merge_qsp((0, 1), (2, 3), (4, 5), (6, 7))
        rho_red = rho6.partial_trace([1])
        ee, _ = rho_eig(rho_red)
        #res[01] = ee
        res[1] = ee
    
    
    #print_vars(vars(),  ['EE, lam'])
    #res= OrderedDict()    
    #EE = {}
    #EE[1] =  ee  
    #res[0] = EE
    return res 


def central_charge(S, force_update=False, verbose=True, info=0):
    """
        refs:
            Pfeifer and Vidall 2008
    """
    is_found = False
    is_loaded = False
    if S is None:
        if path is None: raise Exception('S, path shouldnt be none at the same time')
        S = System.load(path)
        is_loaded = True
    
    if 0:
        record = S.record
        try:
            #print 'try to fetch central_charge from S.record...  ', 
            msg = 'try to fetch central_charge from S.record...  '
            res = record[S.iter]['central_charge']
        except KeyError as err:
            #print 'record not found'
            msg +=  'record not found:  ' + str(err.message) 
            #print err
            is_loaded = True
        else:
            #print 'record found'
            msg += 'record found'
            is_found = True
        if verbose: print msg
    
    if force_update:  is_found = False            
    
    #if not is_found: 
    if 1: 
        EE = {}
        SIlayer = S.mera.num_of_layer-2
        #EE[1], EE[2] = __entanglement_entropy(S, layer=SIlayer)
        EE = __entanglement_entropy(S, layer=SIlayer)
        #print "EE %s"%EE       
        cc = 3*(EE[2]-EE[1])
        digit = 6
        cc = round(cc, digit)
        EE[1] = round(EE[1], digit)
        EE[2] = round(EE[2], digit)
        res= {}
        res['c'] = cc
        res['EE'] = EE
        
    
    #if not is_found and is_loaded:
    if 0: 
        if not S.record.has_key(S.iter):  
            S.record[S.iter]={}
        S.record[S.iter].update({'central_charge':res})
        if verbose:  print 'central_charge at step %d is saved'%S.iter
        S._save(path)

    return res

def G_calc_EE(nsite): 
    if nsite == 6: 
        
        G_dic = {
                'V_0_0': {1: (None, None),  2: (None, None), 3: ('U_0_0', 3), 4: (None, None) }, 
                'U_0_0': {1: (None, None),  2: (None, None),                  4: ('V_0_1', 1)}, 
                'V_0_1': {                  2: (None, None), 3: (None, None), 4: (None, None) }, 
                }
        contract_order = range(3)
        final_order = [('V_0_0', 1), ('V_0_0', 2), ('U_0_0', 1), ('U_0_0', 2), ('V_0_1', 2), ('V_0_1', 3), 
                ('V_0_0', 4), ('V_0_1', 4)]
        ind_offset = 1
    
    if nsite == 9: 
        G_dic = {
                
                #'rho_1_0': {2: }
                #'V_1_0':  {0: ('V_0_0', 3), 1: ('V_0_1', 3), 2: ('V_0_2', 3)}, #3: ('rho_1_0', 1)
                'V_1_0':  {0: ('V_0_0', 3), 1: ('V_0_1', 3), 2: ('V_0_2', 3), 3: (None, None)},  
                'U_0_0':  {0: (None, None), 1: (None, None), 2: ('V_0_0', 2), 3: ('V_0_1', 0)}, 
                'U_0_1':  {0: (None, None), 1: (None, None), 2: ('V_0_1', 2), 3: ('V_0_2', 0)}, 
                'V_0_0':  {0: (None, None), 1: (None, None)                                  }, 
                'V_0_1':  {                 1: (None, None),                                 }, 
                'V_0_2':  {                 1: (None, None), 2: (None, None)}
            }
        
        
        ind_offset = 0
        if 1: 
            final_order =  [
            ('V_0_0' , 0), ('V_0_0' , 1), ('U_0_0' , 0), ('U_0_0' , 1), ('V_0_1' , 1), 
            ('U_0_1' , 0), ('U_0_1', 1),  ('V_0_2', 1), ('V_0_2', 2), 
            ('V_1_0', 3)
                ]
       
    res = Graphics.init_from_dict(G_dic, ind_offset)
    contract_order = range(res.size)
    res.set_contract_order(contract_order)
    res.set_final_order(final_order, ind_offset)

    return  res
    
    
def entanglement_extra(S, k=80, info=0): 
    G = G_calc_EE_gen()
    tnet_sys = TensorNetwork()
    tnet_sys.TG = G
    
    for node_id in range(G.size): 
        node_name = G.names[node_id]
        name, l, p = node_name.split('_')
        if name == 'rho': 
            op = S.rho_2[2][0].copy()
            #print 'rrrr'; op.reverse_qsp()
           
        if name == 'V': 
            op = S.mera.V[int(l)][0].copy()
        if name == 'Vp': 
            op = S.mera.V_dag[int(l)][0].copy()
        op.type_name = node_name
        tnet_sys.tlink[node_id] = op
    
    
    if 1: #finally ordered according to position of lattice sites
        
        temp =  [
        ('V_0_0' , 2),('V_0_0' , 3), ('V_1_0' , 3), ('V_1_1' , 1), ('V_0_1' , 1), ('V_0_1' , 2), 
        ('Vp_0_0', 3),('Vp_0_0', 4), ('Vp_1_0', 4), ('Vp_1_1', 2), ('Vp_0_1', 2), ('Vp_0_1', 3) 
            ]
        G.set_final_order(temp, ind_offset=1)
        
   
    #contract to reduced density mat
    rho_red = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info)
    EE, val=rho_eig(rho_red, k)
    res= {'EE': EE, 'spectrum': val}
        
    return res
   
def entanglement_brute_force(S, k=80, site_num_max=6, plot=False, average=False, info=0): 
    multiple = 1
    if S.symmetry == 'U1': 
        multiple = 2
    
    EE = {}
    spectrum = {}
    #mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: (9, 8, 7)}
    mapper = {
            2: (2, 1, 0), 
            3: (3, 2, 1, 0),  
            6: (6, 5, 4, 3), 
            9: range(9, -1, -1)}
    def set_val(rho_red, site_num, k, info=info): 
        nn = mapper[site_num]
        for n in nn: 
            n1 = n*multiple
            if n == 0: 
                data = rho_red.data[: ]
                #assert  abs(data-1.0)<1e-10
                EE[n1] = 0.0
                spectrum[n1] = data
                continue
            ee, sp = rho_eig(rho_red, k)
            EE[n1] = ee
            spectrum[n1] = sp
            if info>0: 
                print 'sum of spectrum for n=%d  is  %1.15f'%(n, np.sum(sp))
            #rho_red = partial_trace_rho(rho_red, site_num-n)
            rho_red = partial_trace_rho(rho_red, 1, average=average) 
        
    if 1: 
        site_num = 2
        descending_ham(S.mera, S, ilayer=1)
        rho_red = S.rho_2[0][0].copy()
        set_val(rho_red, site_num, k=100)
   
    if 1: 
        site_num = 3
        if not S.only_NN: 
            descending_ham(S.mera, S, ilayer=1)
            rho_red = S.rho_3[0][0].copy()
            set_val(rho_red, site_num, k=100)
      
    if 1:
        site_num = 6
        G = G_calc_EE(site_num)
        tnet_sys = TensorNetwork()
        tnet_sys.TG = G
        temp = TensorNetwork.make_tensor_link(S, G, trans_invar=1)
        tnet_sys.tlink[: G.size] = temp
        A = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info-1)
       
        
        A_dag = A.conjugate(A.rank-2)
        layer = 1
        rho2 = S.rho_2[layer][0].copy()
        
        temp, leg = A.contract(rho2, range(A.rank), [ 1000, 1001, A.rank-2,  A.rank-1])
        v1 = range(temp.rank)
        v2 = [temp.rank-2, temp.rank-1] + [1000 + i for i in range(A_dag.rank-2)]
        rho_red, leg = temp.contract(A_dag, v1, v2)
      
        set_val(rho_red, site_num, k=200)
    
    site_num = 9  
    if site_num <= site_num_max: 
        
        #if site_num > site_num_max: 
            
        G = G_calc_EE(site_num)
        tnet_sys = TensorNetwork()
        tnet_sys.TG = G
        temp = TensorNetwork.make_tensor_link(S, G, trans_invar=1)
        tnet_sys.tlink[: G.size] = temp
        A = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info-1)
       
        layer = 2
        A_dag = A.conjugate(A.rank-1)
        rho2 = S.rho_2[layer][0].copy()
        rho1 = partial_trace_rho(rho2, n=1, average=average)
        
        temp, leg = A.contract(rho1, range(A.rank), [ 1000,A.rank-1])
        
        v1 = range(temp.rank)
        v2 = [temp.rank-1] + [1000 + i for i in range(A_dag.rank-1)]
        rho_red, leg = temp.contract(A_dag, v1, v2 )
        set_val(rho_red, 9, k=80, info=info)
        
    res = {'EE': EE, 'spectrum': spectrum}
    if plot: 
        x, y = zip(*EE.items());
        x = np.array(x)
        y = np.array(y)
        plt.plot(x, y, '.-')
        plt.show()
        plt.plot(np.log2(x), y, '.-')
        plt.show()
        start = 0
        end = 2
        x1 = np.log2(x[start:end]);
        y1 = y[start: end]
        c = np.polyfit(x1, y1, deg=1)
        print 'central_charge is %s'%c
    
    return res

def entanglement_brute_force_6(S, k=80, site_num_max=6, plot=False, info=0): 
    multiple = 1
    #if S.symmetry == 'U1': multiple = 2
    if S.combine_2site or S.symmetry == 'U1': 
        multiple = 2
    
    EE = {}
    spectrum = {}
    #mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: (9, 8, 7)}
    mapper = {6: (6, 5, 4, 3, 2, 1, 0), 3: (3, 2, 1, 0),  2: (2, 1, 0), 9: range(9, -1, -1)}
    def set_val(rho_red, site_num, k, info=info): 
        nn = mapper[site_num]
        for n in nn: 
            n1 = n*multiple
            #print 'nnnnnnnnnnnnn', n 
            if n == 0: 
                data = rho_red.data[: ]
                #assert  abs(data-1.0)<1e-10
                EE[n1] = 0.0
                spectrum[n1] = data
                continue
            ee, sp = rho_eig(rho_red, k)
            EE[n1] = ee
            spectrum[n1] = sp
            if info>0: 
                print 'sum of spectrum for n=%d  is  %1.15f'%(n, np.sum(sp))
            #rho_red = partial_trace_rho(rho_red, site_num-n)
            rho_red = partial_trace_rho(rho_red, 1) 
        
    if 0: 
        site_num = 2
        descending_ham(S.mera, S, ilayer=1)
        rho_red = S.rho_2[0][0].copy()
        set_val(rho_red, site_num, k=100)
   
    if 0: 
        site_num = 3
        if 1: 
            descending_ham(S.mera, S, ilayer=1)
            rho_red = S.rho_3[0][0].copy()
            #print rho_red 
            set_val(rho_red, site_num, k=100)
        #V = S.mera.V[0][0].copy()
      
    if 1:
        site_num = 6
        G = G_calc_EE(site_num)
        tnet_sys = TensorNetwork()
        tnet_sys.TG = G
        temp = TensorNetwork.make_tensor_link(S, G, trans_invar=1)
        tnet_sys.tlink[: G.size] = temp
        A = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info-1)
       
        
        A_dag = A.conjugate(A.rank-2)
        layer = 1
        rho2 = S.rho_2[layer][0].copy()
        
        temp, leg = A.contract(rho2, range(A.rank), [ 1000, 1001, A.rank-2,  A.rank-1])
        v1 = range(temp.rank)
        v2 = [temp.rank-2, temp.rank-1] + [1000 + i for i in range(A_dag.rank-2)]
        rho_red, leg = temp.contract(A_dag, v1, v2)
      
        set_val(rho_red, site_num, k=200)
    
    res = {'EE': EE, 'spectrum': spectrum}
    if plot: 
        x, y = zip(*EE.items());
        x = np.array(x)
        y = np.array(y)
        plt.plot(x, y, '.-')
        plt.show()
        plt.plot(np.log2(x), y, '.-')
        plt.show()
        start = 0
        end = 2
        x1 = np.log2(x[start:end]);
        y1 = y[start: end]
        c = np.polyfit(x1, y1, deg=1)
        print 'central_charge is %s'%c
    
    return res

def entanglement_brute_force_9(S, k=80, site_num_max=9, plot=False, info=0): 
    multiple = 1
    if S.symmetry == 'U1': 
        multiple = 2
        return None
        
    
    EE = {}
    spectrum = {}
    #mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: (9, 8, 7)}
    mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: range(9, -1, -1)}
    def set_val(rho_red, site_num, k, info=info): 
        nn = mapper[site_num]
        for n in nn: 
            n1 = n*multiple
            if n == 0: 
                data = rho_red.data[: ]
                #assert  abs(data-1.0)<1e-10
                EE[n1] = 0.0
                spectrum[n1] = data
                continue
            ee, sp = rho_eig(rho_red, k)
            EE[n1] = ee
            spectrum[n1] = sp
            if info>0: 
                print 'sum of spectrum for n=%d  is  %1.15f'%(n, np.sum(sp))
            #rho_red = partial_trace_rho(rho_red, site_num-n)
            rho_red = partial_trace_rho(rho_red, 1) 
        
    site_num = 9  
    if site_num <= site_num_max: 
        
        #if site_num > site_num_max: 
            
        G = G_calc_EE(site_num)
        tnet_sys = TensorNetwork()
        tnet_sys.TG = G
        temp = TensorNetwork.make_tensor_link(S, G, trans_invar=1)
        tnet_sys.tlink[: G.size] = temp
        A = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info-1)
       
        layer = 2
        A_dag = A.conjugate(A.rank-1)
        rho2 = S.rho_2[layer][0].copy()
        rho1 = partial_trace_rho(rho2, n=1)
        
        temp, leg = A.contract(rho1, range(A.rank), [ 1000,A.rank-1])
        
        v1 = range(temp.rank)
        v2 = [temp.rank-1] + [1000 + i for i in range(A_dag.rank-1)]
        rho_red, leg = temp.contract(A_dag, v1, v2 )
        set_val(rho_red, 9, k=160, info=info)
        
    res = {'EE': EE, 'spectrum': spectrum}
   
    return res

def entanglement_brute_force_9_test(S, k=80, site_num_max=9, plot=False, info=0): 
    multiple = 1
    if S.symmetry == 'U1': 
        multiple = 2
    print 'hhhhh'
    
    EE = {}
    spectrum = {}
    #mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: (9, 8, 7)}
    mapper = {6: (6, 5, 4, 3), 2: (2, 1, 0), 9: range(9, -1, -1)}
    def set_val(rho_red, site_num, k, info=info): 
        nn = mapper[site_num]
        for n in nn: 
            n1 = n*multiple
            if n == 0: 
                data = rho_red.data[: ]
                #assert  abs(data-1.0)<1e-10
                EE[n1] = 0.0
                spectrum[n1] = data
                continue
            ee, sp = rho_eig(rho_red, k)
            EE[n1] = ee
            spectrum[n1] = sp
            if info>0: 
                print 'sum of spectrum for n=%d  is  %1.15f'%(n, np.sum(sp))
            #rho_red = partial_trace_rho(rho_red, site_num-n)
            rho_red = partial_trace_rho(rho_red, 1) 
        
    site_num = 9  
    if site_num <= site_num_max: 
        
        #if site_num > site_num_max: 
            
        G = G_calc_EE(site_num)
        tnet_sys = TensorNetwork()
        tnet_sys.TG = G
        temp = TensorNetwork.make_tensor_link(S, G, trans_invar=1)
        tnet_sys.tlink[: G.size] = temp
        A = tnet_sys.contract_except(G.size, G.contract_order,final_order=G.final_order, info=info-1)
        return A
       
        layer = 2
        A_dag = A.conjugate(A.rank-1)
        rho2 = S.rho_2[layer][0].copy()
        rho1 = partial_trace_rho(rho2, n=1)
        
        temp, leg = A.contract(rho1, range(A.rank), [ 1000,A.rank-1])
        
        v1 = range(temp.rank)
        v2 = [temp.rank-1] + [1000 + i for i in range(A_dag.rank-1)]
        rho_red, leg = temp.contract(A_dag, v1, v2 )
        set_val(rho_red, 9, k=80, info=info)
        
    res = {'EE': EE, 'spectrum': spectrum}
   
    return res

def test(): 
    raise #
    #from merapy.all import *  # this is not allowed
    qsp=QspU1.easy_init([0,1,-1], [2,1,1])
    qsp=qsp.copy_many(18, reverse=range(9,18))
    qn=QnU1.qn_id()

#t=iTensor(18, qsp,qn)
    t.data.size
    9075135300
    t.data.size

    from merapy.hamiltonian import System
    S=System.load('./alpha=2.0/12-4lay.pickle')
    from merapy.measure_and_analysis.central_charge import entanglement_brute_force_9_test
    A = entanglement_brute_force_9_test(S)

    def partial_trace_rho(rho, n):
        """
            n: num of sites to be traced out 
        """
        if n == 0: 
            return rho
        
        #res= rho.partial_trace(0) + rho.partial_trace(1)
        site_total = range(rho.rank//2)
        #site_left = range(n)
        #site_right = range(site_total-n, site_total)
        site_left = site_total[: n]
        site_right = site_total[-n: ]
        if 0: 
            res = rho.partial_trace(site_left) + rho.partial_trace(site_right)  #average between left and right
            res.data *= 0.5
        else: 
            res = rho.partial_trace(site_right)
        return res

    layer = 2
    A_dag = A.conjugate(A.rank-1)

    rho2 = S.rho_2[layer][0].copy()
    rho1 = partial_trace_rho(rho2, n=1)
    temp, leg = A.contract(rho1, range(A.rank), [ 1000,A.rank-1])
    v1 = range(temp.rank)
    v2 = [temp.rank-1] + [1000 + i for i in range(A_dag.rank-1)]
    rho_red, leg = temp.contract(A_dag, v1, v2 )

#print rho_red.data.size
#9075135300
    set_val(rho_red, 9, k=80, info=info)

    f=open('./rho_red.pickle', 'wb')
    pickle.dump(rho_red, f)    

if 1: 
    def entanglement_brute_force_aver(S, k=80, site_num_max=None, plot=False, info=0): 
        res=entanglement_brute_force(S, k, site_num_max, average=1)
        return res

    def entanglement_brute_force_6_aver(S, k=80, site_num_max=6, plot=False, info=0): 
        res=entanglement_brute_force(S, k, site_num_max, average=1)
        return res
        
    def entanglement_brute_force_9_aver(S, k=150, site_num_max=9, plot=False, info=0): 
        if S.symmetry == 'U1': 
            return None
        res=entanglement_brute_force(S, k, site_num_max, average=1)
        return res

def test_wave_func(): 
    """
        构造一些波函数，来精确测试纠缠测量
    """
    pass
    

def negativity(): 
    pass

class TestIt(unittest.TestCase): 
    INITED = False 
    def setUp(self): 
        path_dic = {}
        path_z2 = '/home/zhli/Documents/mera_backup_tensor/run-ising/ternary/z2-symm/scale-invar/h=1.0/4.pickle'           
        #path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/12-4lay.pickle'
        #path = '/home/zhli/windows/C/Users/zhli/Dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python/merapy/mera_backup_test_folder/2.pickle'
        #path = '/home/zhli/windows/C/Users/zhli/Dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python/merapy/run-ising/ternary-symm=triv/2.pickle'     
        path_dic['z2'] = path_z2
        self.path_dic = path_dic

    def initialize(self): 
        if not TestIt.INITED:  #and 0: 
            from merapy.config import CFG_HEISBG_BASIC 
            from merapy.main import Main 
            parpath = tempfile.mkdtemp()
            cfg = CFG_HEISBG_BASIC.copy();   
            
            cfg.update(USE_CUSTOM_RAND=True, 
                    trunc_dim=4, tot_layer=4, 
                    use_player=True, SYMMETRY="U1", NUM_OF_THREADS=1, 
                    backup_parpath = parpath, 
                    )
            cfg['model_param'].update({"J_NN":1.0, "J_NNN":0.0})
            main = Main(**cfg)
            main.run(q_iter=10, do_measure=0)
            #res= {0: 0, 1: -0.15447477974680113, 2: -0.27316843978691985, 3: -0.38435891607150841, 4: -0.49142890463353184, 5: -0.59574502419909914, 6: -0.69595138685370195, 7: -0.78927036090526448, 8: -0.87373436401172677, 9: -0.9490632254128426, 10: -1.015790781309708}
            #main.S.examine_energy(res)
            self.S= main.S 
            TestIt.INITED = True 
       
    def test_central_charge(self): 
        #path = './states_for_tests/dim=4-symm=U1-model=heisbg.pickle'
        path = './states_for_tests/dim=12-lay=4-symm=U1-model=heisbg.pickle'
        S = merapy.load(path)
         
        res = central_charge(S=S, force_update=True, info=0)
        res_old = {'EE': {1: 1.511427, 2: 1.688231}, 'c': 0.530412}
        print res  
    
    def test_entanglement_extra(self):    
        pass 
        #entanglement_extra(S, k=80, info=2)
        #path = '/home/zhli/Documents/mera_backup_tensor/run-xx/ternary/traiv-symm/prod-state/4.pickle' 
        
        #res=entanglement_brute_force_9(S=S, k=80, plot=1, site_num_max=9, info=2);  print res['EE']
  
    def test_entanglement_brute_force_6(self): 
        path = './states_for_tests/dim=4-symm=U1-model=heisbg.pickle'
        #path = './states_for_tests/dim=12-lay=4-symm=U1-model=heisbg.pickle'
        #path = self.path_dic['z2']

        S = merapy.load(path)
        res=entanglement_entropy(S, spectrum=0, info=2, path=None, verbose=False) 
        print  res 
        print central_charge(S=S, force_update=True, info=0)
        #entanglement_extra(S, k=80, info=2)
        
        res = entanglement_brute_force_6(S=S, k=80, plot=0, info=2);  
        print res['EE']

    def test_concurrence(self): 
        #path= '/home/zhli/backup_tensor_dir/run-long-heisbg/mera/test/rmax12_4812125/alpha=0.88-init1p44/8.pickle'
        #path= '/home/zhli/backup_tensor_dir/run-heisbg-NNN/mera/test/prod_state/J2=0.02-pin2_lamm1_xi2000/4.pickle'
        #from merapy.run_heisbg.analysis import an_mera, an_vmps 
        #xx = an_mera.an_main
        #db = xx[1.0, 'a']
        #sh = (20, 4)
        #S= db.load_S(sh)
        #S = load(path)
        #ee=entanglement_special(S)
        self.initialize()
        S= self.S 
        cc=concurrence(S)
        print_vars(vars(),  ['cc'])
        self.assertAlmostEqual(cc, 0.0236695468515, 10)
        
            


    def test_temp(self): 
   
        #path= '/home/zhli/backup_tensor_dir/run-long-heisbg/mera/test/rmax12_4812125/alpha=0.88-init1p44/8.pickle'
        #path= '/home/zhli/backup_tensor_dir/run-heisbg-NNN/mera/test/prod_state/J2=0.02-pin2_lamm1_xi2000/4.pickle'
        from merapy.run_heisbg.analysis import an_mera, an_vmps 
        xx = an_mera.an_main
        db = xx[1.0, 'a']
        sh = (20, 4)
        S= db.load_S(sh)
        #S = load(path)
        #ee=entanglement_special(S)
        ee=concurrence(S, 1.0)
        print_vars(vars(),  ['ee'])
            


if __name__ == '__main__':

    if 0: 
        if len(sys.argv)>1:
            path = sys.argv[1]   
        else:
            if 1:         
                path = '/home/zhli/Documents/mera_backup_tensor/run-ising/ternary/z2-symm/scale-invar/h=1.0/4.pickle'           
                path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/12-4lay.pickle'
                #path = '/home/zhli/windows/C/Users/zhli/Dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python/merapy/mera_backup_test_folder/2.pickle'
                #path = '/home/zhli/windows/C/Users/zhli/Dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python/merapy/run-ising/ternary-symm=triv/2.pickle'     
                #path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/new_run/alpha=1.6/12-4lay.pickle'
                S = System.load(path)
                 
                print central_charge(S=S, force_update=True, info=0)
                print "res should be {'EE': {1: 0.842492, 2: 1.010302}, 'c': 0.50343}" 
                
                #entanglement_extra(S, k=80, info=2)
                #path = '/home/zhli/Documents/mera_backup_tensor/run-xx/ternary/traiv-symm/prod-state/4.pickle' 
                
                #res=entanglement_brute_force(S=S, k=80, plot=0, site_num_max=6, info=2);  print res['EE']
                res=entanglement_brute_force_6(S=S, k=80, plot=0, site_num_max=6, info=2);  print res['EE']
                #res=entanglement_brute_force_9(S=S, k=80, plot=1, site_num_max=9, info=2);  print res['EE']
    
    if 0: #examine
        TestIt.test_temp = unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
    if 1: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_entanglement_brute_force_6', 
           'test_concurrence', 
           #'test_temp', 
       
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)
       




