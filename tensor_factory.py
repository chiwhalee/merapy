#!/usr/bin/env python
#coding=utf8   

from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
from builtins import object
import unittest
import nose
import warnings
import pickle
import pprint
import numpy as np
import itertools 
import math 
from collections import OrderedDict

from merapy.utilities import print_vars
from merapy.tensor_py import iTensor
from merapy.quantum_number_py import  (QspU1, QspZ2, 
        QspTravial, symmetry_to_Qsp, symmetry_to_QspClass, make_qsp, QnU1, QnZ2)

   
def pauli_mat(): 
    """
        Note these are pauli mats,  not their half--- spin operators
        s_alpha means sigma_alpha
    """
    s0 = np.identity(2)
    sx = np.array([[0, 1],[1, 0]] )
    sy = np.array([[0, -1j],[1j, 0]])
    sz = np.array([[1, 0],[0, -1]])   
    #sp = sx + 1j*sy  #a BAD mistake !!! should be 1/2.0(sx + 1j*sy)
    #sm = sx - 1j*sy
    sp = np.array([[0, 1],[0, 0]])   
    sm = np.array([[0, 0],[1, 0]])   
    
    res = vars()
    res[0] = np.zeros((2, 2)); 
    res['X'] = res['sx'];  res['Z'] = res['sz']; res['Y']=res['sy']; res['I'] = res['s0']
    res['op_type'] = 'spin'
    return res

class iTensorFactory(object):
    """
        define some instances of iTensor,  that are frequently used
    """
    def __init__(self, symmetry):
        pass
        self.symmetry = symmetry
        self.qn_identity, self.qsp_base, self.qsp_null = init_System_QSp(symmetry, dim=2)

    @staticmethod
    def simple(rank, dim, symmetry, nqn=None, reverse=None):
        qnclass = symmetry_to_Qn(symmetry)
        qspclass = symmetry_to_Qsp(symmetry)
        qn = qnclass.qn_id()
        qsp = [qspclass.max(d, nqn) for  d in dim]
        if reverse is not None:
            for i in reverse:
                qsp[i].reverse()
        return iTensor(rank, qsp, qn)
    
    @staticmethod
    def random(qsp, totqn=None, dtype=float): 
        """
            create an iTensor whose data is randomized 
        """
        if isinstance(totqn, int):
            totqn = qsp[0].QnClass(totqn)
        res = iTensor(QSp=qsp, totQN=totqn, dtype=dtype)
        #issue:  maybe here need devide total size?
        n = res.data.size
        if dtype == float: 
            res.data[: ] = np.random.random(n) - 0.5 
        elif dtype == complex:
            res.data[: ] = (np.random.random(n) - 0.5  + 
                    1j*(np.random.random(n) - 0.5))
        else:
            raise  
        return res 
    
    @staticmethod
    def isometry():
        pass
    
    @staticmethod
    def diagonal_tensor_rank2(qsp): 
        """
            a tensor t satisfis:  tt* = t*t  =  I 
            not any qsp can make diagonal_tensor_rank2
            the condision is for each qn and -qn in the QNs, their dim
            must equal, otherwise, it is actually conceptially an isometry tensor 
            e.g.  q = QspU1.easy_init([0, 1, -1], [2, 3, 3]) will do
            while  q = QspU1.easy_init([0, 1, -1], [2, 3, 5]) will raise 
        """
        res= iTensor(QSp=qsp.copy_many(2))
        for i in range(res.nidx): 
            sh=res.get_block_shape(i) 
            d = sh[0]
            temp = np.identity(d, dtype=res.dtype).ravel()            
            res.set_block(i, temp)
        return res 
    
    @staticmethod
    def pauli_mat_2site(symmetry, which=None):
        """
            todo: 
                it is better to implement 2site through 1site and merge_qsp
                in this way, not only code easier, but also the totqn of the ops 
                are more well defiend 
            
        """
        qspclass=symmetry_to_Qsp(symmetry)
        qn_identity, qsp_base, qsp_null = qspclass.set_base(dim=4)  #combine 2site
        temp = lambda: (2, qsp_base.copy_many(2, reverse=[1]),  qn_identity.copy())
       
        rank, qsp, totqn = temp()
        
        sigma_z1=iTensor(*temp());  sigma_z1.data[:]=0.0    #sigma_z1 = sz×I
        sigma_z2=iTensor(*temp());  sigma_z2.data[:]=0.0    #sigma_z2 = I×sz
        ii = iTensor(*temp());  ii.data[:]=0.0  #ii = I×I
        
        qDims= np.ndarray(2, "int")
        iDims= np.ndarray(2, "int")

        if symmetry == "Travial":
            pauli_1 = iTensorFactory.pauli_mat_1site(symmetry)
            p1 = pauli_1
            sigma_z1.data[:] = p1["sigma_z"].direct_product(p1["sigma_0"]).data[:]
            sigma_z2.data[:] = p1["sigma_0"].direct_product(p1["sigma_z"]).data[:]
            ii.data[:] = p1["sigma_0"].direct_product(p1["sigma_0"]).data[:]
            
            sigma_p1 = iTensor(*temp())
            sigma_p2 = iTensor(*temp())
            sigma_p1.data[:] = p1["sigma_p"].direct_product(p1["sigma_0"]).data[:]
            sigma_p2.data[:] = p1["sigma_0"].direct_product(p1["sigma_p"]).data[:]

            sigma_m1 = iTensor(*temp())
            sigma_m2 = iTensor(*temp())
            sigma_m1.data[:] = p1["sigma_m"].direct_product(p1["sigma_0"]).data[:]
            sigma_m2.data[:] = p1["sigma_0"].direct_product(p1["sigma_m"]).data[:]

            sigma_x1 = iTensor(*temp())
            sigma_x2 = iTensor(*temp())
            sigma_x1.data[:] = p1["sigma_x"].direct_product(p1["sigma_0"]).data[:]
            sigma_x2.data[:] = p1["sigma_0"].direct_product(p1["sigma_x"]).data[:]

            sigma_y1 = iTensor(*temp())
            sigma_y2 = iTensor(*temp())
            sigma_y1.data[:] = p1["sigma_y"].direct_product(p1["sigma_0"]).data[:]
            sigma_y2.data[:] = p1["sigma_0"].direct_product(p1["sigma_y"]).data[:]
        elif symmetry == 'Z2':  
            
            if 1 : 
                if 1: 
                    """
                        HowTo: combine two site for Z2 symm 
                        note CONVENTION: 
                            qn 1, -1 are spin up, down resp. 
                            |a><b|  其中 <b|是张量朝下的腿
                        since combine 2 sites, each qn has two dim
                        qsp_base:  1: (uu, dd); -1: (ud, du)
                        
                        merge qsp (index):  
                            |uu>,|dd> => e(ven) => 0,  
                            |uu> ->0, |dd> ->1
                            <uu| ->0, <dd| ->1
                            (uu, uu) => (0, 0), (uu, dd) => (0, 1), (dd, uu) =>(1, 0), (dd, dd)=>(1, 1)
                            
                            |ud>, |du> => o(dd) => 1 
                            <du| ->0, <ud| ->1
                            |ud> ->0, |du> ->1
                            (ud, ud) => (0, 0), (ud, du) => (0, 1), (du, ud) =>(1, 0), (du, du)=>(1, 1)
                    """
                    
                    """
                        sigma_z(1)%data(1:8) = (/1.d0,0.d0,0.d0,-1.d0,  0.d0,1.d0,1.d0, 0.d0/)
                        sigma_z(2)%data(1:8) = (/1.d0, 0.d0, 0.d0,-1.d0, 0.d0,-1.d0,-1.d0, 0.d0/)      

                        sz  = 1*|u><u|+ (-1)*|d><d|;  id = 1*|u><u|+ 1*|d><d|
                        sz \otimes id = 1*|uu><uu| + 1*|ud><du| + (-1)*|du><ud| + (-1)*|dd><dd| 
                       
                    """
                    sigma_z1.set_element(qDims=[0, 0], iDims=[0, 0], X=1.0 ) 
                    sigma_z1.set_element(qDims=[1, 1], iDims=[0, 0], X=1.0 ) 
                    sigma_z1.set_element(qDims=[1, 1], iDims=[1, 1], X=-1.0 ) 
                    sigma_z1.set_element(qDims=[0, 0], iDims=[1, 1], X=-1.0 ) 
                    """
                        id = 1*|u><u|+ 1*|d><d|; sz  = 1*|u><u|+ (-1)*|d><d|;  
                        id \otimes sz = 1*|uu><uu| + (-1)*|ud><du| + 1*|du><ud|  + (-1)*|dd><dd|
                    """
                    sigma_z2.set_element(qDims=[0, 0], iDims=[0, 0], X=1.0 ) 
                    sigma_z2.set_element(qDims=[1, 1], iDims=[0, 0], X=-1.0 ) 
                    sigma_z2.set_element(qDims=[1, 1], iDims=[1, 1], X=1.0 ) 
                    sigma_z2.set_element(qDims=[0, 0], iDims=[1, 1], X=-1.0 ) 
                    
                    #sigma_ii = sigma_z1.copy()
                    #sigma_ii.data[: ] = 0.0
                    sigma_00 = iTensorFactory.identity(qsp_base.copy_many(2, reverse=[1]) )
                    ii  =  sigma_00 
                    
                if 1:
                    rank, qsp, totqn = temp()
                    totqn.set_val(-1)  
                    sigma_x1 = iTensor(rank, qsp, totqn); sigma_x1.data[:] = 0.0
                    
                    sigma_x2 = sigma_x1.copy(); sigma_x2.data[: ]=0.0
                    sigma_y1 = sigma_x1.copy(); sigma_y1.data[: ]=0.0
                    sigma_y2 = sigma_x1.copy(); sigma_y2.data[: ]=0.0
                   
                    """
                          sigma_x(1)%data(1:8) = [1.0,-1.0,1.0,1.0, 1.0,1.0,-1.0,1.0]
                            sigma_x(2)%data(1:8) = [1.0,1.0,1.0,-1.0, 1.0,1.0,1.0,-1.0]

                        sx = 1*|u><d| + 1*|d><u|; id = 1*|u><u|+ 1*|d><d|
                        sx \otimes id = 1*|uu><ud| + 1*|ud><dd| + 1*|du><uu| + 1*|dd><du|
                        4 terms依次对应于 : 
                    """
                       
                    if 1: 
                        sigma_x1.set_element(qDims=[0, 1], iDims=[1,0], X=1.0)   #注意 iDims 第1位对应于 <.|, 第二位对应于 |.>
                        sigma_x1.set_element(qDims=[1, 0], iDims=[1,0], X=1.0) 
                        sigma_x1.set_element(qDims=[1, 0], iDims=[0,1], X=1.0) 
                        sigma_x1.set_element(qDims=[0, 1], iDims=[0,1], X=1.0) 
                    
                    """
                        sy = (-1)*|u><d| + 1*|d><u|; id = 1*|u><u|+ 1*|d><d|
                        sy \otimes id = (-1)*|uu><ud| + (-1)*|ud><dd| + 1*|du><uu| + 1*|dd><du|
                        4 terms依次对应于 : 
                    """
                    #note sigma_y1 与通常的差个因子 i
                    sigma_y1.set_element(qDims=[0, 1], iDims=[1,0], X=-1.0)   #注意 iDims 第1位对应于 <.|, 第二位对应于 |.>
                    sigma_y1.set_element(qDims=[1, 0], iDims=[1,0], X=-1.0) 
                    sigma_y1.set_element(qDims=[1, 0], iDims=[0,1], X=1.0) 
                    sigma_y1.set_element(qDims=[0, 1], iDims=[0,1], X=1.0) 
                   
                    """
                        id = 1*|u><u|+ 1*|d><d|; sx = 1*|u><d| + 1*|d><u| 
                        id \otimes sx = 1*|uu><du| + 1*|ud><uu| + 1*|du><dd| + 1*|dd><ud|
                        4 terms依次对应于 : 
                    """
                    sigma_x2.set_element(qDims=[1, 0], iDims=[0, 0], X= 1.0)  #注意 iDims 第1位对应于 <.|, 第二位对应于 |.>
                    sigma_x2.set_element(qDims=[0, 1], iDims=[0, 0], X= 1.0) 
                    sigma_x2.set_element(qDims=[1, 0], iDims=[1, 1], X= 1.0) 
                    sigma_x2.set_element(qDims=[0, 1], iDims=[1, 1], X= 1.0) 
                    
                    """
                        id = 1*|u><u|+ 1*|d><d|; sy = (-1)*|u><d| + 1*|d><u| 
                        id \otimes sx = (-1)*|uu><du| + 1*|ud><uu| + (-1)*|du><dd| + 1*|dd><ud|
                        4 terms依次对应于 : 
                    """
                   
                    #note sigma_y2 与通常的差个因子 i
                    sigma_y2.set_element(qDims=[1, 0], iDims=[0, 0], X=-1.0) #注意 iDims 第1位对应于 <.|, 第二位对应于 |.>
                    sigma_y2.set_element(qDims=[0, 1], iDims=[0, 0], X= 1.0) 
                    sigma_y2.set_element(qDims=[0, 1], iDims=[1, 1], X=-1.0) 
                    sigma_y2.set_element(qDims=[1, 0], iDims=[1, 1], X= 1.0) 
            
            else :   #both work
                    print('in tensor_py.py changeed '*30)
                    pau1 = iTensorFactory.pauli_mat_1site(symmetry=symmetry)
                    s0, sx, sy, sz = pau1['s0'], pau1['sx'], pau1['sy'], pau1['sz']
                    
                    sigma_00 = (s0.direct_product(s0)).index_merge_simple()
                    
                    sigma_z1 = (sz.direct_product(s0)).index_merge_simple()
                    sigma_x1 = (sx.direct_product(s0)).index_merge_simple() 
                    sigma_y1 = (sy.direct_product(s0)).index_merge_simple()
                                                    
                    sigma_z2 = (s0.direct_product(sz)).index_merge_simple()
                    sigma_x2 = (s0.direct_product(sx)).index_merge_simple()
                    sigma_y2 = (s0.direct_product(sy)).index_merge_simple()
                    
            szi = sigma_z1; sxi = sigma_x1; syi = sigma_y1 
            siz = sigma_z2; six = sigma_x2; siy = sigma_y2 
            
            sigma_zz = sigma_z1.contract_core(sigma_z2,1)
            sigma_xx = sigma_x1.contract_core(sigma_x2,1)
            sigma_yy = -1.0 * sigma_y1.contract_core(sigma_y2,1)   #-1.0 comes from i^2
            
            s00 = sigma_00
            szz = sigma_zz 
            sxx = sigma_xx
            syy = sigma_yy
            sigma_01 = ii.copy()
            sigma_02 = ii.copy()

            #print 'sss', ss; exit()

        elif symmetry  == "U1":
            if 1:   #only meant for folding and editing
                # [-1,-1] -> [(-1)-1,(-1)-1]
                qDims[0:2] = [2,2]
                iDims[0:2] = [0,0]
                sigma_z1.set_element(qDims, iDims, -1.0)
                sigma_z2.set_element(qDims, iDims, -1.0)
                ii.set_element(qDims, iDims, 1.0)
                
                # [1,1] -> [1,1]
                qDims[0:2] = [1,1]
                iDims[0:2] = [0,0]
                sigma_z1.set_element(qDims, iDims, 1.0)
                sigma_z2.set_element(qDims, iDims, 1.0)
                ii.set_element(qDims, iDims, 1.0)
                
                # [-1,1] -> [-1(-1),1]
                qDims[0:2] = [0,0]
                iDims[0:2] = [0,0]
                sigma_z1.set_element(qDims, iDims, -1.0)
                sigma_z2.set_element(qDims, iDims, 1.0)
                ii.set_element(qDims, iDims, 1.0)
                
                # [1,-1] -> [1,-1(-1)]
                qDims[0:2] = [0,0]
                iDims[0:2] = [1,1]
                sigma_z1.set_element(qDims, iDims, 1.0)
                sigma_z2.set_element(qDims, iDims, -1.0)
                ii.set_element(qDims, iDims, 1.0)
                
                sigma_01 = ii.copy()
                sigma_02 = ii.copy()
                
            if 1: 
                rank, qsp, totqn = temp()
                #here val = 1 is absolutely right.
                #_val = 1 just means spin = 1,  not 1/2 !!
                #issue:  here totqn for sp should be -1, see note1 in doc of .spin_one_mat, 
                #but the value wont affect the final results, so not change it at present 
                totqn.set_val(1)  #note this differs with sigma_z
                sigma_p1 = iTensor(rank, qsp, totqn)
                #sigma_p1 := sp×I
                sigma_p1.data[:] = 0.0
                sigma_p2 = sigma_p1.copy()
                #print_vars(vars(),  ['sigma_p1'])
                
                if 0:  # old def of iTensor with totQN reversed 
                    # [-1,-1]->[1,-1]
                    qDims[0:2] = [2,0]   #|ud><dd|
                    iDims[0:2] = [0,1]
                    sigma_p1.set_element(qDims, iDims, 1.0)
                    
                    # [-1,1]->[1,1]
                    qDims[0:2] = [0,1]  #|uu><du|
                    iDims[0:2] = [0,0]
                    sigma_p1.set_element(qDims, iDims, 1.0)
                else: 
                    # [-1,-1]->[1,-1]
                    qDims[0:2] = [1,0]   #|ud><dd|
                    iDims[0:2] = [0,0]
                    sigma_p1.set_element(qDims, iDims, 1.0)
                    
                    # [-1,1]->[1,1]
                    qDims[0:2] = [0,2]  #|uu><du|
                    iDims[0:2] = [1,0]
                    sigma_p1.set_element(qDims, iDims, 1.0)
                
                if 0: 
                    # [-1,-1]->[-1,1]
                    qDims[0:2] = [2,0]  
                    iDims[0:2] = [0,0]
                    sigma_p2.set_element(qDims, iDims, 1.0)
                    # [1,-1]->[1,1]        
                    qDims[0:2] = [0,1]  
                    iDims[0:2] = [1,0]
                    sigma_p2.set_element(qDims, iDims, 1.0)
                else: 
                    # [-1,-1]->[-1,1]
                    qDims[0:2] = [1,0]  
                    iDims[0:2] = [0,1]
                    sigma_p2.set_element(qDims, iDims, 1.0)
                    # [1,-1]->[1,1]        
                    qDims[0:2] = [0,2]  
                    iDims[0:2] = [0,0]
                    sigma_p2.set_element(qDims, iDims, 1.0)
                    
                
                sigma_m1=sigma_p1.conjugate(1)
                sigma_m2=sigma_p2.conjugate(1)
                
                if 0:  #following are wrong, as sigma_p and sigma_m cant be added as U1 covariant tensors
                       # or to say, sigma_x, sigma_y are not U1 covariant
                    sigma_x1 = (sigma_p1 + sigma_m1)
                    sigma_x2 = (sigma_p2 + sigma_m2)
                    
                    sigma_y1 = (sigma_p1 - sigma_m1)  #these are not exactly sigma_y, miss a factor of -I
                    sigma_y2 = (sigma_p2 - sigma_m2)

        #nearest neighbour interaction, 2-site operator 
        #if symmetry in ["Z2", "Travial"]:
        if symmetry in ["Travial"]:
            sxx = sigma_x1.contract_core(sigma_x2,1)
            syy = sigma_y1.contract_core(sigma_y2,1)
            szi = sigma_z1
            siz = sigma_z2
        
        if symmetry in ["U1", "Travial"]:  
            szz = sigma_z1.contract_core(sigma_z2,1) #szz is a rank 2 tensor  # sz_i ×sz_{i+1} 
            spm=sigma_p1.contract_core(sigma_m2, 1)  # sp_i × sm_{i+1}
            smp = sigma_m1.contract_core(sigma_p2, 1)  # sm_i ×sp_{i+1} 
            
        if 1:
            assert sigma_z1.is_hermite()
            assert sigma_z2.is_hermite()
            assert szz.is_hermite()
            
        return vars()
    
    @staticmethod
    def pauli_mat_1site(symmetry):
        qspclass=symmetry_to_Qsp(symmetry)
        qn_identity, qsp_base, qsp_null = qspclass.set_base()
        if symmetry  == "Z2": 
            temp = lambda: (2, qsp_base.copy_many(2, reverse=[1]),  qn_identity.copy())

            rank, qsp, totqn = temp()
            sigma_z=iTensor(rank, qsp, totqn)
            sigma_z.data[0:2] = [1.0,-1.0]
            
            rank, qsp, totqn = temp()
            sigma_0=iTensor(rank, qsp, totqn)
            sigma_0.data[0:2] = [1.0,1.0]
            
            rank, qsp, totqn = temp()
            totqn.set_val(-1)
            sigma_x=iTensor(rank, qsp, totqn)
            sigma_x.data[0:2] = [1.0,1.0]
            
            rank, qsp, totqn = temp()
            totqn.set_val(-1)
            sigma_y=iTensor(rank, qsp, totqn)
            #note this is not standard sigma_y, rather 
            sigma_y.data[0:2] = [-1.0,1.0]
            
            I, sx, sy, sz = sigma_0.copy(), sigma_x.copy(), sigma_y.copy(), sigma_z.copy()    
            I, X, Y, Z = sigma_0.copy(), sigma_x.copy(), sigma_y.copy(), sigma_z.copy()    

            if 1:
            #following are not fully verified
                if 0:  # cant define sigma^\pm with z2 symm !!! 
                    rank, qsp, totqn = temp()
                    totqn.set_val(-1)
                    sigma_p=iTensor(rank, qsp, totqn)
                    sigma_p.data[0:2] = [1.0, 0.0]

                    rank, qsp, totqn = temp()
                    totqn.set_val(-1)
                    sigma_m=iTensor(rank, qsp, totqn)
                    sigma_m.data[0:2] = [0.0, 1.0]

                    spm = sigma_p.tensor_prod(sigma_m)
                    smp = sigma_m.tensor_prod(sigma_p)

                sxx = sigma_x.tensor_prod(sigma_x)
                syy = sigma_y.tensor_prod(sigma_y)
                szz = sigma_z.tensor_prod(sigma_z)

        elif symmetry ==  "Travial":
            temp = lambda: (2, qsp_base.copy_many(2, reverse=[1]),  qn_identity.copy())

            rank, qsp, totqn = temp()
            sigma_z=iTensor(rank, qsp, totqn)
            sigma_z.data[:] = [1.0, 0, 0, -1.0]
            
            rank, qsp, totqn = temp()
            sigma_0=iTensor(rank, qsp, totqn)
            sigma_0.data[:] = [1.0, 0.0, 0.0, 1.0]
            
            rank, qsp, totqn = temp()
            sigma_x=iTensor(rank, qsp, totqn)
            sigma_x.data[:] = [0.0, 1.0, 1.0, 0.0]
            
            rank, qsp, totqn = temp()
            sigma_y=iTensor(rank, qsp, totqn)
            #note this is not standard sigma_y, rather sigma_y/i
            sigma_y.data[:] = [0.0, -1.0, 1.0, 0.0]

            rank, qsp, totqn = temp()
            sigma_p=iTensor(rank, qsp, totqn)
            sigma_p.data[:] = [0.0, 1.0, 0.0, 0.0]

            rank, qsp, totqn = temp()
            sigma_m=iTensor(rank, qsp, totqn)
            sigma_m.data[:] = [0.0, 0.0, 1.0, 0.0]
            i, z, x, y, p, m = sigma_0, sigma_z, sigma_x, sigma_y, sigma_p, sigma_m
            I, Z, X, Y, P, M = sigma_0, sigma_z, sigma_x, sigma_y, sigma_p, sigma_m
            I, sp, sm, sz = sigma_0, sigma_p, sigma_m, sigma_z    
            
            if 1:
                #due to Fortran order,  modified to this
                sigma_p.data[:] = [0.0, 0.0, 1.0, 0.0]
                sigma_m.data[:] = [0.0, 1.0, 0.0, 0.0]

            order = "F"  # this order has no effect
            spm = sigma_p.tensor_prod(sigma_m, order)
            smp = sigma_m.tensor_prod(sigma_p, order)

            sxx = sigma_x.tensor_prod(sigma_x, order)
            syy = sigma_y.tensor_prod(sigma_y, order)
            szz = sigma_z.tensor_prod(sigma_z, order)
            
            ii = i.direct_product(i)
            sigma_z1 = z.direct_product(i)
            sigma_z2 = i.direct_product(z)
            sigma_p1 = p.direct_product(i)
            sigma_p2 = i.direct_product(p)
            sigma_m1 = m.direct_product(i)
            sigma_m2 = i.direct_product(m)
        
        elif symmetry  == "U1":
            #qns1 = [1, -1]
            qsp_base = QspU1.easy_init(qns=(1, -1), dims=(1, 1))
            qn_identity = QnU1.qn_id()
        
            def temp():
                rank = 2
                qsp = qsp_base.copy_many(2, reverse=[1])
                totqn = qn_identity.copy()
                return rank, qsp, totqn
            
            #rank, qsp, totqn = temp()
            sigma_z = iTensor(*temp())
            #sigma_z = iTensor(rank,qsp,totqn)
            sigma_z.data[:] = [1.0, -1.0]
            
            #rank, qsp, totqn = temp()
            sigma_0 = iTensor(*temp())
            #sigma_0 = iTensor(rank, qsp, totqn)
            sigma_0.data[:] = [1.0, 1.0]
            
            rank, qsp, totqn = temp()
            #issue:  here totqn for sp should be -2, see note1 in doc of .spin_one_mat, 
            #but the value wont affect the final results, so not change it at present 
            totqn.set_val(2)  
            
            sigma_p = iTensor(rank, qsp, totqn)
            sigma_p.data[0] = 1.0
            sigma_m = sigma_p.conjugate(1)
            
            szz = sigma_z.tensor_prod(sigma_z)
            spm = sigma_p.tensor_prod(sigma_m)
            smp = sigma_m.tensor_prod(sigma_p)
            
            I, sp, sm, sz = sigma_0, sigma_p, sigma_m, sigma_z    
            #I, X, Y, Z = sigma_0, sigma_x, sigma_y, sigma_z    
            
            I, Z = sigma_0, sigma_z 
            
            #spin_z = 0.5*sigma_z
        F = -sigma_z   # this used in jordan-wigner trans
        temp = vars()
        res= {k: v for k, v in temp.items() if isinstance(v, iTensor)}
        for i in res: 
            res[i].type_name = i 
        zero = res['I'].copy()
        zero.data[:] = 0.0
        res['zero'] = zero
        
        res['op_type'] = 'spin'
            
        return res 
    
    pauli_mat = pauli_mat_1site 
    
    @staticmethod
    def spin_one_mat( symmetry='U1'): 
        """
            irreducible representation matrices of su2 lie algebra labeled by [spin]. 
            these are spin mat (not corresponding to pauli mat)
            note1: 
                if spin = 1
                    o = sum_{m', m in [1, 0, -1]} { a_{m', m} |m'><m| } 
                    for sz,  requirs m' = m
                    for sp, in the textbook of QM, it requirs m' = m + 1. It
                        should be that m' = -(m + 1), so the tot_qn= m' + m = -1,
                        instead of 1. 其实不矛盾，其要点是，base (1>, 0>, -1>)^t conj ->
                        (-1>, 0>, 1>), 定义i (i')为m (m')在基矢中的位置，则 i'= i +
                        1.  换句话说，在QM中，总是用量子数的位置标记矩阵元a_{i',
                        i}，这里直接用量子数标记
                    for sm, m' = -(m-1) => tot_qn=1
        """
        pass 
        qspclass = symmetry_to_Qsp(symmetry)
        qnclass = qspclass.QnClass 
        #qn_identity, qsp_base, qsp_null = qspclass.set_base()
        qn_identity = qspclass.QnClass.qn_id()
        qsp_null = qspclass.null()
        
        if symmetry == "Travial":
            
            qsp_base = qspclass.easy_init([1], [3])
            
            temp = lambda: (2, qsp_base.copy_many(2, reverse=[1]),  qn_identity.copy())

            rank, qsp, totqn = temp()
            sz=iTensor(rank, qsp, totqn)
            sz.data[:] = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dtype=np.float64).ravel(order='F')
            
            rank, qsp, totqn = temp()
            s0=iTensor(rank, qsp, totqn)
            s0.data[:] = np.identity(3, dtype=np.float64).ravel(order='F')
            
            rank, qsp, totqn = temp()
            sx=iTensor(rank, qsp, totqn)
            sx.data[:] = 1/np.sqrt(2)*np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.float64).ravel(order='F')
            
            rank, qsp, totqn = temp()
            sy=iTensor(rank, qsp, totqn)
            #note this is not standard sy, rather sy/i
            sy_complex = 1/np.sqrt(2)*np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]], dtype=np.complex128).ravel(order='F')
            sy.data[:] = np.asarray(sy_complex/1j, dtype=float)

            rank, qsp, totqn = temp()
            sp=iTensor(rank, qsp, totqn)
            sp.data[:] = np.asarray(sx.data + 1j*sy_complex, dtype=float)
            sm = sp.conjugate(1)
            i, z, x, y, p, m = s0, sz, sx, sy, sp, sm
            I, Z, X, Y, P, M = s0, sz, sx, sy, sp, sm
            

            order = "F"  # this order has no effect
            spm = sp.tensor_prod(sm, order)
            smp = sm.tensor_prod(sp, order)

            sxx = sx.tensor_prod(sx, order)
            syy = sy.tensor_prod(sy, order)
            szz = sz.tensor_prod(sz, order)
            
            ii = i.direct_product(i)
            sz1 = z.direct_product(i)
            sz2 = i.direct_product(z)
            sp1 = p.direct_product(i)
            sp2 = i.direct_product(p)
            sm1 = m.direct_product(i)
            sm2 = i.direct_product(m)
        
        elif symmetry == "U1":
            qsp_base = qspclass.easy_init([1, 0, -1], [1, 1, 1])
            qn_identity = QnU1.qn_id()
        
            def temp():
                rank = 2
                qsp = qsp_base.copy_many(2, reverse=[1])
                totqn = qn_identity.copy()
                return rank, qsp, totqn
            
            #rank, qsp, totqn = temp()
            sz = iTensor(*temp())
            #sz = iTensor(rank,qsp,totqn)
            sz.data[:] = [1.0, 0.0, -1.0]

            s0 = iTensor(*temp())
            #s0 = iTensor(rank, qsp, totqn)
            s0.data[:] = [1.0, 1.0, 1.0]
            
            rank, qsp, totqn = temp()
            totqn.set_val(-1)  
            sp = iTensor(rank, qsp, totqn)
            sp.data[:] = math.sqrt(2.)
            sm = sp.conjugate(1)

            if 0:
                print(sz) #.matrix_view()
                print(sp.matrix_view())
                print(sm.matrix_view())
            
            szz = sz.tensor_prod(sz)
            spm = sp.tensor_prod(sm)
            smp = sm.tensor_prod(sp)
            
            #I, X, Y, Z = s0, sx, sy, sz    
            
            I, Z = s0, sz 
        else: 
            raise 
        
        temp = vars()
        res= {k: v for k, v in temp.items() if isinstance(v, iTensor)}
        for i in res: 
            res[i].type_name = i 
            
        return res 
    
    def check_pauli_mat(self):
        if self.symmetry  ==  "U1":
            pauli = self.pauli_mat()
            ii, szz, spm, smp= pauli["ii"], pauli["szz"], pauli["spm"], pauli["smp"]
            sz1, sz2= pauli["sigma_z1"], pauli["sigma_z2"]
            sp1, sp2= pauli["sigma_p1"], pauli["sigma_p2"]
            sm1, sm2= pauli["sigma_m1"], pauli["sigma_m2"]
            sigma_z1, sigma_z2, sigma_p1, sigma_p2, sigma_m1, sigma_m2 = sz1, sz2, sp1, sp2, sm1, sm2
            ii = ii
            zi, iz, pi, ip, mi, im = sz1, sz2, sp1, sp2, sm1, sm2

            
            np.set_printoptions(threshold=np.nan)
        if self.symmetry  == "Travial":
            pauli = self.pauli_mat()
            I_2, szz, spm, smp= pauli["sigma_0"], pauli["szz"], pauli["spm"], pauli["smp"]
            pm, mp = spm, smp
            sigma_z, sigma_p, sigma_m = pauli["sigma_z"], pauli["sigma_p"], pauli["sigma_m"]
            x, y=  pauli["sigma_x"], pauli["sigma_y"]
            zz, xx, yy=  pauli["szz"], pauli["sxx"], pauli["syy"]

            z, p, m = sigma_z, sigma_p, sigma_m
            i = I_2
            print("check using two methods to define h2 and h3")
            h2a = xx  + (-1.0)*yy
            h2b = 2.0*(pm + mp)
            print(np.all(h2a.data==h2b.data))

            h3a = (-1.0)*y.direct_product(i).direct_product(y) + x.direct_product(i).direct_product(x)
            h3b = 2.0*(p.direct_product(i).direct_product(m) + m.direct_product(i).direct_product(p))
            print(np.all(h3a.data==h3b.data))
        if self.symmetry == "Z2":
            pauli = self.pauli_mat()
            I_2, szz, spm, smp= pauli["sigma_0"], pauli["szz"], pauli["spm"], pauli["smp"]
            pm, mp = spm, smp
            sigma_z, sigma_p, sigma_m = pauli["sigma_z"], pauli["sigma_p"], pauli["sigma_m"]
            x, y=  pauli["sigma_x"], pauli["sigma_y"]
            zz, xx, yy=  pauli["szz"], pauli["sxx"], pauli["syy"]

            z, p, m = sigma_z, sigma_p, sigma_m
            i = I_2
            
            print("check using two methods to define h2 and h3")
            h2a = xx  + (-1.0)*yy
            h2b = 2.0*(pm + mp)
            print(h2a.data)
            print(h2b.data)

            h3a = (-1.0)*y.direct_product(i).direct_product(y) + x.direct_product(i).direct_product(x)
            h3b = 2.0*(p.direct_product(i).direct_product(m) + m.direct_product(i).direct_product(p))
            print(np.all(h3a.data==h3b.data))

    @staticmethod
    def fermion_op(symmetry, shift_qn=True):
        """
            ref http://hedrock.ps.uci.edu/docs.cgi?page=tutorials/fermions
            the basis
                {|0>, |u>, |d>, |ud>}
                note a convention here
                    |ud> = cdag_up cdag_dn|0>
                    
            params:
                shift_qn: if true, this assumes add a postive electron at each site
                    this has benifits that at half filling MPS has total qn of 0
                    
        """
        vec = iTensorFactory.base_state(which='fermion', symmetry=symmetry, 
                shift_qn=shift_qn)
        def make_op(xx):
            op = 0.0
            for a, b, c in xx:
                #print_vars(vars(),  ['a', 'b'])
                a = vec[a].insert_1d_qsp(1)  #|a>
                b = vec[b].insert_1d_qsp(1) #|b>
                a_tc = a.T.conj()
                #a_tc.totQN.reverse()
                cba = c*b.dot(a_tc)  #  + c|b><a|
                op += cba
            return op
        
        cdag_up = [('0', 'u', 1.0), ('d', 'ud', 1.0)]  # maps |0> to |u>,  |d> to |ud>
        cdag_up = make_op(cdag_up)
        c_up = cdag_up.T.conj()  #[('u', '0', 1.0), ('ud', 'd', 1.0)]
        
        cdag_dn = [('0', 'd', 1.0), ('u', 'ud', -1.0)]
        cdag_dn = make_op(cdag_dn)
        c_dn = cdag_dn.T.conj()  #[('d', '0', 1.0), ('ud', 'u', -1.0)]
        
        if 1:
            #F = (-1)^n_i  used in the jordan-wigner str 
            F = [('0', '0', 1.0), ('u', 'u', -1.0), ('d', 'd', -1.0), ('ud', 'ud', 1.0)]
            F = make_op(F)
            
            adag_up = [('0', 'u', 1.0), ('d', 'ud', 1.0)]  
            adag_up = make_op(adag_up)
            a_up = adag_up.T.conj()  
            
            adag_dn = [('0', 'd', 1.0), ('u', 'ud', 1.0)]
            adag_dn = make_op(adag_dn)
            a_dn = adag_dn.T.conj()  
        
        if 1:
            n_up = cdag_up.dot(c_up)
            n_dn = cdag_dn.dot(c_dn)
            n_i = n_up + n_dn
            n_up_prod_n_dn = n_up.dot(n_dn)
            
            I = iTensor.identity(n_up.QSp[0].copy_many(2, reverse=[1]) )
        
        temp = ['I', 
                'cdag_up', 'cdag_dn', 'c_up', 'c_dn', 
                'adag_up', 'adag_dn', 'a_up', 'a_dn', 
                'F', 
                'n_i', 'n_up', 'n_dn', 'n_up_prod_n_dn'] 
        
        dic = locals()
        res= {a:dic[a] for a in temp}
        for i in res:
            res[i].type_name = i
        
       
            
        return res 
    
    @staticmethod
    def boson_op(symmetry, nmax, shift_qn=1):
        """
            nmax:
                max num of boson allowed on a site. 
                nmax = 2 is the hard core limit i.e U = infty.
                usually nmax=5 is enough due to impulsion from U
                
                After truncation of nmax [b, b^+1] no longer strictly equals 1 
            
        """
        vec = iTensorFactory.base_state(which='boson', nmax=nmax,  
                symmetry=symmetry, shift_qn=shift_qn)
        def make_op(xx):
            op = 0.0
            for a, b, c in xx:
                a = vec[a].insert_1d_qsp(1)  #|a>
                b = vec[b].insert_1d_qsp(1) #|b>
                a_tc = a.T.conj()
                #a_tc.totQN.reverse()
                cba = c*b.dot(a_tc)  #  + c|b><a|
                op += cba
            return op
        sqrt = math.sqrt
        bdag = [('%d'%n, '%d'%(n+1), sqrt(n+1)) for n in range(nmax) ]  
        
        bdag = make_op(bdag)
        b = bdag.T.conj()  #[('u', '0', 1.0), ('ud', 'd', 1.0)]
        I = [('%d'%i, '%d'%i, 1.0) for i in range(nmax + 1) ]  
        I = make_op(I)
        n_i = bdag.dot(b)
        temp = ['I', 
                'b', 'bdag'  , 
                'F', 
                'n_i', 'n_up', 'n_dn', 'n_up_prod_n_dn'] 
        
        dic = locals()
        res= {a:dic[a] for a in temp if a in dic}
        for i in res:
            res[i].type_name = i
        
        return res 
    
    @staticmethod
    def lattice_op(op_type, symmetry, spin=None, nmax=None):
        """
            wrapper of spin, fermion, boson operators into one func
            params:
                nmax: max num of bosons. required by boson op
        """
        if op_type == 'spin':
            if spin == 'one_half': 
                mapping = iTensorFactory.pauli_mat(symmetry) if symmetry is not None else pauli_mat()
            elif spin == 'one':
                if symmetry is None: 
                    mapping = TensorFactory.spin_mat(1)
                else: 
                    raise  NotImplemented 
        elif op_type == 'fermion':
            mapping = iTensorFactory.fermion_op(symmetry)
        elif op_type == 'boson' :
            mapping = iTensorFactory.boson_op(symmetry, nmax)
        else:
            raise  
        
        mapping['op_type'] = op_type
        
        return mapping 
    
    @classmethod
    def common_op(cls, op_name, which=None, symmetry=None, nmax=None, force=0):
        """
            for convenience 
            params:
                nmax: only for bosons 
                force: only needed for bosons, when changing nmax 
            
        """
        #if not force and hasattr(cls, op_name):
        #    return getattr(cls, op_name)
        
        if which == 'boson' or 'b' in op_name:
            temp = cls.lattice_op('boson', symmetry=symmetry, nmax=nmax)
        elif 'sigma' in op_name:
            temp = cls.lattice_op('spin', symmetry=symmetry, spin='one_half')
        else:
            raise NotImplemented  
        res = temp[op_name]
        
        return res 
    
    @staticmethod
    def base_state(which, symmetry, nmax=None, shift_qn=True, **kwargs):
        """
            params:
                shift_qn: for fermions and bosons
            
        """
        vec = OrderedDict()
        if which == 'spin':
            raise NotImplemented
        elif which == 'fermion':
            names = ['0', 'u', 'd', 'ud']
            if symmetry ==  "Travial":
                q = QspTravial.easy_init([1], [4])
                #e = QspTravial.easy_init([1], [1])
                for i, n in enumerate(names) :
                    #vec[n] = iTensor(QSp=[q.copy(), e.copy()])
                    vec[n] = iTensor(QSp=[q.copy()])
                    vec[n].data[:] = 0.0
                    vec[n].data[i] = 1.0
                    
            elif symmetry == 'U1':  #charge symmetry
                if not shift_qn:
                    q = make_qsp(symmetry, [0, 1, 2], [1, 2, 1])
                    totqn_dic = {'0':0, 'u':1, 'd':1, 'ud':2}
                else:  #I have tested with hubbard model, the following is correct in all cases
                    q = make_qsp(symmetry, [-1, 0, 1], [1, 2, 1])
                    totqn_dic = {'0':-1, 'u':0, 'd':0, 'ud':1}
                    
                for i, n in enumerate(names) :
                    totqn = QnU1(totqn_dic[n])
                    #e = make_qsp(symmetry, [0], [1])
                    #vec[n] = iTensor(QSp=[q.copy(), e.copy()], totQN=totqn)
                    vec[n] = iTensor(QSp=[q.copy()], totQN=totqn)
                    vec[n].data[:] = 0.0
                
                
                vec['0'].data[:] = [1.0]
                vec['u'].data[:] = [1.0, 0.0]
                vec['d'].data[:] = [0.0, 1.0]
                vec['ud'].data[:] = [1.0]
            else:
                raise ValueError
        
        elif which == 'boson':
            assert nmax is not None 
            #the dimension of the basis = nmax+1
            names = ['%d'%i for i in range(nmax + 1)]
            if symmetry == 'Travial':
                q = QspTravial.easy_init([1], [nmax + 1])                
                for i, n in enumerate(names) :
                    vec[n] = iTensor(QSp=[q.copy()])
                    vec[n].type_name = '|%s>'%i 
                    vec[n].data[:] = 0.0
                    vec[n].data[i] = 1.0
            elif symmetry == 'U1':
                for n, name in enumerate(names):
                    if not shift_qn:
                        totqn =  QnU1(n) 
                        q = make_qsp(symmetry, 
                            list(range(nmax + 1)), 
                            [1]*(nmax + 1), 
                            )
                    else:
                        assert nmax%2 == 0 
                        nmax_half = (nmax )//2
                        totqn =  QnU1(n-nmax_half) 
                        q = make_qsp(symmetry, 
                            list(range(-nmax_half, nmax_half + 1)), 
                            [1]*(nmax + 1), 
                            )
                        
                    t = iTensor(QSp=[q.copy()], totQN=totqn)   # a rank-1 tensor 
                    t.data[:] = [1.0]
                    vec[name] = t
            else:
                raise  
        else:
            raise 
        return vec
    
    @staticmethod
    def symbol_tensor_prod(symbol, mapper):
        o = mapper[symbol[0]]
        for a in symbol[1:]:
            o = o.direct_product(mapper[a])
        return o
            
    @staticmethod
    def identity(qsp):
        return iTensor.unit_tensor(2, qsp)
    
    @staticmethod
    def swap(qsp): 
        """
            used especially in fermionic mera
        """
        raise NotImplemented

    @staticmethod
    def travial_tensor(rank, symmetry, dtype=float): 
        """
            scalar 1 with dummy indices 
        """
        qsp_class= symmetry_to_QspClass(symmetry)
        qsp = [qsp_class.null() for i in range(rank)]
        res= iTensor(QSp=qsp, dtype=dtype)
        if dtype == float:  
            res.data[0] = 1.0
        else:  #complex type,  if not doing this, numpy would warn
            res.data[0] = 1.0 + 0j 
        return res 
    
class TestIt(unittest.TestCase): 
    def test_diagonal_tensor_rank2(self): 
        qsp = QspU1.easy_init([0, 1, -1], [4, 2, 2])
        t = iTensorFactory.diagonal_tensor_rank2(qsp )
        t.show_data()
        print(t.to_ndarray())
            
    def test_spin_one_mat(self): 
        #t = iTensorFactory.spin_one_mat('U1')
        t = iTensorFactory.spin_one_mat('Travial')
        #t = iTensorFactory.pauli_mat_1site('U1')
        sp = t['sp']
        print_vars(vars(),  ['t.keys()', 'sp.to_ndarray()'])

    def test_fermion_op(self):
        for symmetry in ['U1', 'Travial']:
            res = iTensorFactory.fermion_op(symmetry)
            cdag_up, cdag_dn, c_up, c_dn = res['cdag_up'], res['cdag_dn'], res['c_up'], res['c_dn']
            n_up_prod_n_dn = res['n_up_prod_n_dn']
            n_up, n_dn = res['n_up'], res['n_dn']
            #test {c, c^+} = 1
            temp=[c_up.commutator_plus(cdag_up).is_close_to(1), 
            c_dn.commutator_plus(cdag_dn).is_close_to(1)
                ]
            self.assertTrue(all(temp))
        
            ni = res['n_i']
            self.assertTrue(np.all(ni.matrix_view().diagonal()==[0, 1, 1, 2]))

    def test_boson_op(self):
        for symm in ['Travial', 'U1']:
            nmax = 4
            res = iTensorFactory.boson_op(symm, nmax, 
                    shift_qn=1)
            bdag, b = res['bdag'], res['b']
            I, n_i  =  res['I'], res['n_i']
            
            # test [b, b^+]  is almost 1 
            c = b.commutator(bdag)  #After truncation of nmax [b, b^+1] no longer strictly equals 1 
            c = c.to_ndarray().diagonal().round(5)
            
            
            self.assertTrue(np.all(c[:-1]==[1]*nmax))
            print_vars(vars(),  ['b.dot(bdag)'])
            
            # test n_i 
            n_i = n_i.to_ndarray().diagonal().round(5)
            print_vars(vars(),  ['n_i'])
            self.assertTrue(np.all(n_i==np.arange(0, nmax + 1, 1.0)))
        
    def test_temp(self):
        
        t = iTensorFactory.common_op('sigma_x', symmetry='Z2', nmax=4)
        print_vars(vars(),  ['t'])
            
        


if __name__ == "__main__":
    
    if 0:
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main(verbosity=-10)
       
    else: 
        suite = unittest.TestSuite()

    add_list_iTF = [
        #'test_diagonal_tensor_rank2', 
        #'test_spin_one_mat', 
        #'test_fermion_op', 
        #'test_boson_op', 
        'test_temp', 
        
            ]
  
    for a in add_list_iTF: 
        suite.addTest(TestIt(a))
    unittest.TextTestRunner().run(suite)

