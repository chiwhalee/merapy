#coding=utf8
#!/usr/bin/env python

"""  
    tensor.f90
    todo: Improvements:
        1. 区分 leg direction
        2. 统一用ndarray 存储
        3. add label for legs

    Questions: q:
        9. 所有leg quantum number都相同的话， 可以作为iTensor的标签，而不必要跟leg绑定，
        10. IMPORTANT  iTensor 稀疏存储在 numpy有现成的没有，或者部分用numpy实现，而非
            at least  data 可以用 ndarray存储
        12  how to deal with rank 0 tensor?
            rank 0 tensor can be taken as arbitrary rank tensor, only that the dim of each index
            is 1. then every thing is ok
        13 什么时候用copy，什么时候不用  this is a QUESTION!!

    2012-07-13 
    it turns out type nTensor may be completely replaced by class np.ndarray :
        tensor product  --> np.multiply.outer
        tensor contract  --> np.tensordot
        exchange index  --> np.swqpaxes  only pairwise
        get_position  --> ? seems no direct method, but simply use np.arange(totDim,Dims), just yeilds what I need
        get_position_rev --> np.unravel_index(linear_ind, shape)

"""
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()


from builtins import str
from builtins import map
from builtins import range
from builtins import *
from builtins import object
import unittest
import warnings
import pickle
import pprint
import time 
import numpy as np
import itertools 
import math 
from collections import OrderedDict
from multiprocessing import Pool 


from scipy.linalg.blas import dgemm, sgemm, zgemm
from numpy.lib.stride_tricks import as_strided
import numpy.core.numeric as _nx
try:        
    from py3nj import (clebsch_gordan, wigner3j, wigner6j, wigner9j)
except:
    warnings.warn('py3nj is not installed')


import scipy

try:
    import cupy as cp
except:
    pass

#from quantum_number import *  #QuantSpace, QN_idendity, QSp_null, QSp_base
from merapy.utilities import (print_vars, save, load)
from merapy.ntensor import TensorBase, nTensor 
#from merapy.quantum_number import *  

#from merapy.quantum_number_py  import (QspU1, QspZ2, QspTravial, qsp_any, symmetry_to_Qsp, make_qsp)
#from merapy.quantum_number_py  import (QspU1, QspZ2, QspTravial, qsp_any, make_qsp)
#from quantum_number_py import (QspU1, QspZ2, QspTravial, qsp_any, symmetry_to_Qsp, make_qsp)
from merapy.quantum_number import *

#from merapy import quantum_number_py

from merapy import common_util


#from merapy import array_permutation


from merapy.utilities import get_local

from merapy.decorators import (tensor_player, decorate_methods, set_player_state_manual, set_player_state_auto)
#from merapy.tensor_player_multiple import decorate_methods,  tensor_player, set_player_state_auto

#from merapy.decorators import (tensor_player, decorate_methods, 
#        reset_tensor_player, set_player_state_auto)
#

#import scipy.weave as weave
#import numexpr


__all__ = ['TensorBase', 'nTensor', 'iTensor', ]


class DataBuffer_bac(object):
    """
    the point of tBUffer is :当tensor用过自动garbage collect后，只要把tensor使用
    的buff标记in_use改为-1，而buff并未必删除，可以继续使用
    """
    def __init__(self, size, dtype, use_gpu=0):
        """
            size: number of tensors
        """
        #self.T=[Tensor() for i in xrange(size)]
        #make self.T simpliy a place holder
        if use_gpu != 1:
            self.T=[np.ndarray(0, dtype=dtype) for i in range(size)]
        else:
            self.T=[cp.ndarray(0, dtype=dtype) for i in range(size)]
            
        #for i in xrange(size): self.T[i].nullify()
        self.in_use=np.ndarray(1024,"bool")
        self.in_use[:] = False
        self.size=size   
        
    def delete(size):
        for i  in range(tB.size):
            tB.T[i].delete()
        self.in_use[:] = False
        sef.size = 0
        self.T = []

class DataBuffer(object):
    """
    the point of tBUffer is :当tensor用过自动garbage collect后，只要把tensor使用
    的buff标记in_use改为-1，而buff并未必删除，可以继续使用
    """
    def __init__(self, size, use_gpu=0):
        """
            size: number of tensors
            note1:
                np.uint8 is 1 Byte (8 bit); Then float32 is 4*unint8,  float64
                is 8*unint8,  ....  This makes the DataBuffer not bind to
                certain dtypes. 
        """
        #self.T=[Tensor() for i in xrange(size)]
        #make self.T simpliy a place holder
        if use_gpu != 1:
            self.T=[np.ndarray(0, dtype=np.uint8) for i in range(size)]  #note1
        else:
            self.T=[cp.ndarray(0, dtype=np.uint8) for i in range(size)]
            
        #for i in xrange(size): self.T[i].nullify()
        self.in_use=np.ndarray(1024, dtype=bool)
        self.in_use[:] = False
        self.size=size   
        
    def delete(size):
        for i  in range(tB.size):
            tB.T[i].delete()
        self.in_use[:] = False
        sef.size = 0
        self.T = []



class ndindex:
    """
    An N-dimensional iterator object to index arrays.

    Given the shape of an array, an `ndindex` instance iterates over
    the N-dimensional index of the array. At each iteration a tuple
    of indices is returned, the last dimension is iterated over first.

    Parameters
    ----------
    `*args` : ints
      The size of each dimension of the array.

    See Also
    --------
    ndenumerate, flatiter

    Examples
    --------
    >>> for index in np.ndindex(3, 2, 1):
    ...     print(index)
    (0, 0, 0)
    (0, 1, 0)
    (1, 0, 0)
    (1, 1, 0)
    (2, 0, 0)
    (2, 1, 0)

    """

    def __init__(self, *shape, order='C'):
        if len(shape) == 1 and isinstance(shape[0], tuple):
            shape = shape[0]
        x = as_strided(_nx.zeros(1), shape=shape,
                       strides=_nx.zeros_like(shape))
        self._it = _nx.nditer(x, flags=['multi_index', 'zerosize_ok'],
                              order=order)

    def __iter__(self):
        return self

    def ndincr(self):
        """
        Increment the multi-dimensional index by one.

        This method is for backward compatibility only: do not use.
        """
        next(self)

    def __next__(self):
        """
        Standard iterator method, updates the index and returns the index
        tuple.

        Returns
        -------
        val : tuple of ints
            Returns a tuple containing the indices of the current
            iteration.

        """
        next(self._it)
        return self._it.multi_index


#meth_names= ["__init__", "set_data_entrance", "contract_core", "permutation"]
#meth_names.pop(0)

#print( 'DISABLED DECORATOR '*10)
@decorate_methods(decorator=tensor_player, meth_names=None)
class iTensor(TensorBase):
    """
        todo: 
            1. 考虑：无论是block index 还是 data index都应该统一改成 C-order 
            2. 考虑：qsp应该弄成 immutable type,  成为一种 const,  避免copy，提高效率
            3. consider:  use qn itsef as identification of block directly,  rather than id of qn 
    
    """
    num_of_instance = 0
    #issue:  these buffer may not compatible with complex type 
    DATA_BUFFER = [DataBuffer(size=20) for  i in range(2)]
    DATA_BUFFER_GPU = [DataBuffer(size=20) for  i in range(2)]
    
    #USE_GPU_FOR_BLOCK = False
    USE_GPU_MUL_LIM = 100**3   # used when use_gpu = 2
    USE_GPU_EIG_LIM = None 
    
    def __init__(self, rank=None, QSp=None, totQN=None, order='F', dtype=np.float64, 
            buffer=None, use_buf=False, index_data=True, 
            has_data=True, init_data=None, use_gpu=0):
        """
            params: 
                totQN : 
                    决定了张量的变换性质 totQN = id for invariant, other for covariant
                    default is qn_id 

                attention_here here the default ordering of quantum number combination is changed from fortran to C 
                
                self.Dims[i]:  
                    在leg i 对应的空间维数
                buf_ref: 
                    buf_ref[0]---which iTensor.DATA_BUFFER, buf_ref[1]---which tensor in the iTensor.DATA_BUFFER
                    default -1 means not using iTensor.DATA_BUFFER
                use_gpu:
                    can take value in [0, 1, 2]. 
                        0 means not use gpu  
                        1 means use gpu  
                        2 means use gpu for certain large blocks of data.
                            When to use gpu depend on USE_GPU_MUL_LIM, USE_GPU_EIG_LIM
                    
            issue:
                in future, QSp only searve as a way to construct iTensor, forbit it being modified, because modifying qsp alone
                is meanless and dangeous if not copyed.
                this raise a question: should qsp immutable, permanante or mutalbe?
                In this way, one only need to construct very fewer distinct QSps,  only reference, no copying
            QUESTION, Q:
                QSp, totQN应该在init外面copy，还是在里面copy？ copy or not this is a question
                -a.  这只是个惯例保持一致即可。但从概念上看，似乎应该在外面更好，因为初始化一个张量是一个具体的张量，那么传
                进来的QSp也应该是某"个"具体的QSp，而非某一"类"QSp。
                并且，这样做的好处是，方便在init外控制 copy，or reference
            theory: 
                the guide line of manipulating symmetric tensor:
                    is that always do operation first at the block level then in the inner-block level 
                about rank-0 tensor. 
                    this is conceptially important. 
                    怎样理解它？
                        it is a scalar 
                        它既是0阶，又是任意阶 
                            其实这样并不稀奇，任意的张量，总可以安上 dummy ind, 也是任意阶 
                            所以，rank 应该有个 strick-rank 这一说，把dummy ind 全去掉
                    how to set its attr?
                        rank = 0
                        QSp = []
                        data = np.ndarray((1,))
                        问题是 idx 应不应该设？—— 需要. 指标的变化范围就是1 
                
        """
        #TensorBase.__init__(self, rank, None)  #comment this only for a little faster
        rank = len(QSp)
        self.rank = rank
        self.QSp = QSp  #NO COPYING CONVENTION 
        self.ind_labels = None 
        self.totQN = totQN if totQN is not None else QSp[0].QnClass.qn_id()
          
        self.ndiv = None   #ndiv 实际是把协变反变腿分开 一个(ndiv, rank-ndiv) tensor
        self.use_buf=use_buf
        self.buf_ref = np.array([-1, -1], int)
        self.type_name = ''
        self.use_gpu = use_gpu
        
        if rank == 0:   #Dims 指的是对应的**dense** tensor 的维数
            self.Dims = np.array((1, ), int) 
        else: 
            self.Dims = [QSp[i].totDim for i in range(rank)]

        if index_data:  #this sets idx, nidx, idx_dim, totDim, Addr_idx, Block_idx 
            if rank == 0 or QSp[0].IS_Abelian:
                self.set_data_entrance(order=order )
            else:
                self.set_data_entrance_non_Abelian(order=order )
        
        
        if has_data:
            if buffer is None and use_buf:   #use internal DATA_BUFFER; else use external buffer or no buffer
                if self.use_gpu != 1:
                    #buffer = self.buffer_assign(data_size=self.totDim if dtype==float else self.totDim*2) 
                    buffer = self.buffer_assign(data_size=self.totDim, item_size=np.dtype(dtype).itemsize) 
                else:
                    #buffer = self.buffer_assign_gpu(data_size=self.totDim 
                    #                                if dtype==float else self.totDim*2)  
                    buffer = self.buffer_assign_gpu(data_size=self.totDim, item_size=np.dtype(dtype).itemsize) 
                   
            if self.use_gpu != 1:
                self.data = np.ndarray(self.totDim, buffer=buffer, dtype=dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered
            else:
                if isinstance(buffer, cp.ndarray):
                    buffer = buffer.data 
                else:
                    assert buffer is None
                self.data = cp.ndarray(self.totDim, memptr=buffer, dtype=dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered
    
    def __setstate__(self, d): 
        self.__dict__.update(d)
        self.buf_ref = np.array([-1, -1], int)   # this is important 
        if not hasattr(self, 'use_gpu'):   # add at 2023/ 08/13
            self.use_gpu = 0
    
    def __eq__(self, other):
        """
            still need wanshan 
        """
        if not isinstance(other, iTensor):
            return False
        if not self.QSp == other.QSp:
            return False
        if not np.allclose(self.data, other.data, atol=1e-14):
            return False
        return True
    
    def __matmul__(self, other):
        """
            the @ operator 
        """
        return self.contract(other, return_v3=False)
    
    def buffer_assign(self, data_size, item_size,  n=None):
        """
            let self point to a iTensor.DATA_BUFFER
            n: which buff to use
                n = 0: for small rank tensor
                n = 1: for lager rank tensor
            #issue-unfixed: one need buffer for rank 4, 6, 8, 10 etc
            as for Ternary graph with oo ops, there could be rank 4, 6 tensors
            with ooo ops, there could be rank 8 tensors
        """
        if n is None:  
            if self.rank <= 4:  
                n = 0
            else:
                n = 1
        for i in range(iTensor.DATA_BUFFER[n].size):
            if not iTensor.DATA_BUFFER[n].in_use[i]:  
                self.buf_ref[0]=n
                self.buf_ref[1]=i                
                iTensor.DATA_BUFFER[n].in_use[i] = True 
                num_of_bytes = data_size * item_size
                if iTensor.DATA_BUFFER[n].T[i].size<num_of_bytes:
                    #iTensor.DATA_BUFFER[n].T[i] = np.empty(data_size, dtype=np.float64)
                    iTensor.DATA_BUFFER[n].T[i] = np.empty(num_of_bytes, dtype=np.uint8)
                
                return iTensor.DATA_BUFFER[n].T[i].data
        raise Exception('Error, All buffer elements are in use, stop %s\n '%(str(iTensor.DATA_BUFFER[0].in_use[:100], )))
    
    def buffer_assign_gpu(self, data_size, item_size,  n=None):
        """
        """
        if n is None:  
            if self.rank <= 4:  
                n = 0
            elif self.rank <= 6:
                n = 1
            else:
                n = 2
        for i in range(iTensor.DATA_BUFFER_GPU[n].size):
            if not iTensor.DATA_BUFFER_GPU[n].in_use[i]:  
                self.buf_ref[0]=n
                self.buf_ref[1]=i                
                iTensor.DATA_BUFFER_GPU[n].in_use[i] = True 
                num_of_bytes = data_size * item_size
                if iTensor.DATA_BUFFER_GPU[n].T[i].size<data_size:
                    iTensor.DATA_BUFFER_GPU[n].T[i] = cp.empty(num_of_bytes, dtype=np.uint8)
                
                return iTensor.DATA_BUFFER_GPU[n].T[i].data
        raise Exception('Error, All buffer elements are in use, stop %s\n '%(str(iTensor.DATA_BUFFER_GPU[0].in_use[:100], )))
    
    
   
    
    def set_data_entrance(self, order="F"):
        """
            实现1维稀疏存储, 对应于下式的第二个->号
            T->(D, S)->(dat, Block_ind, Addr_ind)
            this step is in effect 把Qsp 中的信息提取出来，变成更容易读取操作的信息，故
            Qsp 中包含的信息和 Block_idx, Addr_idx 等是等价的(不完全等价, 差一个对称性限制条件)，不同在存储顺序和方式
            Block_idx: is a map from idx to block info
            
            note1: 
                at 2015-8-31, I changed the def of iTensor. The previous weng's treatment
                let sum(qn of each leg) = reverse(totqn),  now change to 
                    sum(qn of each leg) = totqn
                前者其实绕了个弯，把totqn理解成了一个dummy leg 并且是conj的，完全没有必要这么做
                它会造成一定概念上的混乱。后者更加 well defined 
                
            # I removed the following code,  use np.nditer instead 
                #遍历所有的量子数组合
                if order == 'F':   #这里实际上意味着按照 fortran order 对量子数组合排序的
                    inc = True
                    i = 0
                    while inc and i<rank:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<QSp[i].nQN :
                            inc = False
                        else:
                            iQN[i] = 0
                            i = i+1
                            
                elif order == 'C':
                    inc = True
                    i = rank-1
                    while inc==True and i>= 0:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<self.QSp[i].nQN :
                            inc = False
                        else:
                            iQN[i] = 0
                            i = i-1
        
        """
        rank, QSp, totQN = self.rank, self.QSp, self.totQN
        rank_1 = rank if rank != 0 else 1 
       
        if rank == 0: 
            QSp = [totQN.qsp_class().null()]   # only use it temporarilly to generate idx
            
        
        arg = [QSp[i].nQN for i in range(rank_1)]
        #iqn = np.ndindex(*arg)  #iQN like a pointer, 用于给量子数组合编号  #iQN[i]用作leg i 上的量子数 计数
                                # if using 'C' order,  this conflicts with get_idx
        iqn = ndindex(*arg, order='F')  #I added an arg "order" for np.ndindex
        self.idx_dim = np.prod(arg, dtype=int)

        self.idx = np.ndarray(self.idx_dim, int)   #-1                
        self.idx[: ] = self.idx_dim   #self.idx[: ] = -1 
        self.Block_idx = np.ndarray((3, self.idx_dim), dtype=int, order='F')  #here use F order such that access in mem is much faster. todo: transpose Block_idx and use C order 
        self.Addr_idx = np.ndarray((rank_1, self.idx_dim), dtype=int, order='F')  # 实际使用的addr_inx的长度为 self.nidx        
        
        nidx=0
        totDim=0
        totqn = self.totQN 
        
        #todo: this loop can be made to be paralell 
        p = -1
        for iQN in iqn:
            p += 1 
            tqn = QSp[0].QNs[iQN[0]]  #这里计算了总量子数 tqni, 用于判断量子数组合是否满足指定的对称性要求, 这个不其眼的一步实际上是核心——实现了稀疏存储
            for i in range(1, rank):
                #tqn = tqn + QSp[i].QNs[iQN[i]]  
                tqn = tqn.__add__(QSp[i].QNs[iQN[i]]) 
            
            if tqn == totqn:
                d = 1
                for i in range(rank):  #计算某一block的data size 
                    d = d*QSp[i].Dims[iQN[i]]
                
                self.idx[p] = nidx  #给出了0量子数组合与所有量子数组合的序号间的关系 self.idx 和 self.Block_idx[2]互为反函数 如果总量子数为0，则idx[p] =- 1(默认值) 在self.block中都是记录不为0的量子数组合
                
                self.Block_idx[0, nidx] = totDim  #data block在self.data中的position
                self.Block_idx[1, nidx] = d   #data block变成1d数组的长度
                self.Block_idx[2, nidx] = p   #position in quantum number combinations
                
                #记录不为0的量子数组合，在所有量子数组合中的位置 Addr实为将（QN1, ..., QNn)-> ind 的映射, 将n个指标拉直了, 对每一个iTensor都定义了这个函数 Addr_idx这个二维数组的每一列实际上是所有非零block的量子数的编号(而不是量子数点值！)的组合
                self.Addr_idx[0, nidx] = 0 #for rank=0
                self.Addr_idx[0:rank, nidx] = iQN[0:rank]
                nidx += 1 
                totDim += d   #最终得到self.data 的总长度
            
        self.Addr_idx = self.Addr_idx[:, :nidx]
        self.Block_idx = self.Block_idx[:, :nidx]
        self.nidx = nidx
        self.totDim = totDim

    def set_data_entrance_parallel(self, order="F"):
        """
            实现1维稀疏存储, 对应于下式的第二个->号
            T->(D, S)->(dat, Block_ind, Addr_ind)
            this step is in effect 把Qsp 中的信息提取出来，变成更容易读取操作的信息，故
            Qsp 中包含的信息和 Block_idx, Addr_idx 等是等价的(不完全等价, 差一个对称性限制条件)，不同在存储顺序和方式
            Block_idx: is a map from idx to block info
            
            note1: 
                at 2015-8-31, I changed the def of iTensor. The previous weng's treatment
                let sum(qn of each leg) = reverse(totqn),  now change to 
                    sum(qn of each leg) = totqn
                前者其实绕了个弯，把totqn理解成了一个dummy leg 并且是conj的，完全没有必要这么做
                它会造成一定概念上的混乱。后者更加 well defined 
                
            # I removed the following code,  use np.nditer instead 
                #遍历所有的量子数组合
                if order == 'F':   #这里实际上意味着按照 fortran order 对量子数组合排序的
                    inc = True
                    i = 0
                    while inc and i<rank:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<QSp[i].nQN :
                            inc = False
                        else:
                            iQN[i] = 0
                            i = i+1
                            
                elif order == 'C':
                    inc = True
                    i = rank-1
                    while inc==True and i>= 0:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<self.QSp[i].nQN :
                            inc = False
                        else:
                            iQN[i] = 0
                            i = i-1
        
        """
        rank, QSp, totQN = self.rank, self.QSp, self.totQN
        rank_1 = rank if rank != 0 else 1 
       
        if rank == 0: 
            QSp = [totQN.qsp_class().null()]   # only use it temporarilly to generate idx
        arg = [QSp[i].nQN for i in range(rank_1)]
        iqn = np.ndindex(*arg, order='F')  #iQN like a pointer, 用于给量子数组合编号  #iQN[i]用作leg i 上的量子数 计数
        self.idx_dim = np.prod(arg, dtype=int)
        
        n = 1
        pool = Pool(processes=n, initializer=None)
        temp = pool.map(self.helper, iqn)
        temp = [t for t in enumerate(temp)if t[1]]
        size = len(temp)
        
        self.idx = np.ndarray(self.idx_dim, int)   #-1                
        self.idx[: ] = self.idx_dim   #self.idx[: ] = -1 
        self.Block_idx = np.ndarray((3, size), dtype=int, order='F')  #here use F order such that access in mem is much faster. todo: transpose Block_idx and use C order 
        self.Addr_idx = np.ndarray((rank_1, size), dtype=int, order='F')  # 实际使用的addr_inx的长度为 self.nidx        
        
        nidx=0
        totDim=0
        
        for t in temp:
            p, (iQN, d) = t
            self.idx[p] = nidx  #给出了0量子数组合与所有量子数组合的序号间的关系 self.idx 和 self.Block_idx[2]互为反函数 如果总量子数为0，则idx[p] =- 1(默认值) 在self.block中都是记录不为0的量子数组合
            
            self.Block_idx[0, nidx] = totDim  #data block在self.data中的position
            self.Block_idx[1, nidx] = d   #data block变成1d数组的长度
            self.Block_idx[2, nidx] = p   #position in quantum number combinations
            
            #记录不为0的量子数组合，在所有量子数组合中的位置 Addr实为将（QN1, ..., QNn)-> ind 的映射, 将n个指标拉直了, 对每一个iTensor都定义了这个函数 Addr_idx这个二维数组的每一列实际上是所有非零block的量子数的编号(而不是量子数点值！)的组合
            self.Addr_idx[0, nidx] = 0 #for rank=0
            self.Addr_idx[0:rank, nidx] = iQN#[0:rank]
            nidx += 1 
            totDim += d   #最终得到self.data 的总长度
        
        self.nidx = nidx
        self.totDim = totDim
    
        

    def helper(self, iQN):
        pass
        tqn = self.QSp[0].QNs[iQN[0]]  #这里计算了总量子数 tqni, 用于判断量子数组合是否满足指定的对称性要求, 这个不其眼的一步实际上是核心——实现了稀疏存储
        for i in range(1, self.rank):
            tqn = tqn + self.QSp[i].QNs[iQN[i]]  
        if tqn == self.totQN:
            d = 1
            for i in range(self.rank):  #计算某一block的data size 
                d = d*self.QSp[i].Dims[iQN[i]]
            return (iQN, d)
        else:
            return False

    def set_data_entrance_non_Abelian(self, order="F"):
        """
            ref:
                mculloch 2007 
                singh 2012
        
        """
        rank, QSp = self.rank, self.QSp
        rank_1 = rank if rank != 0 else 1 
        temp = 1   
        for i in range(rank_1):
            temp *= QSp[i].nQN
        self.idx_dim = temp   # 量子数组合 总数目 
        
        self.idx = np.ndarray(self.idx_dim, int)   #-1                
        self.idx[: ] = self.idx_dim   #self.idx[: ] = -1 
        self.Block_idx = np.ndarray((3, self.idx_dim), dtype=int, order='F')  #here use F order such that access in mem is much faster. todo: transpose Block_idx and use C order 
        self.Addr_idx = np.ndarray((rank_1, self.idx_dim), dtype=int, order='F')  # 实际使用的addr_inx的长度为 self.nidx        
        
        self.cg_coeff = np.ndarray(self.idx_dim, float)  
        
        iQN = np.zeros(rank_1, dtype=int)   #iQN 用于给量子数组合编号  #iQN[i]用作leg i 上的量子数 计数
        nidx=0
        totDim=0
        totqn = self.totQN 
        
        
        for p in range(self.idx_dim):
            match = False
            if rank == 2: 
                if QSp[0].QNs[iQN[0]]  == QSp[1].QNs[iQN[1]].conj():
                    match = True 
            elif rank == 3:
                j0, m0 = QSp[0].QNs[iQN[0]].val 
                j1, m1 = QSp[1].QNs[iQN[1]].val 
                j2, m2 = QSp[2].QNs[iQN[2]].val 
                
                #print_vars(vars(),  ['j0, j1, m0, m1', 'j2, m2'])
                try:
                    #cg = clebsch_gordan(j0, j1, j2, m0, m1, -m2)
                    cg = clebsch_gordan(j0, j1, j2, m0, m1, m2)
                except:
                    cg = 0
                match = True if cg else False
             
            elif rank == 4:
                j0, m0 = QSp[0].QNs[iQN[0]].val 
                j1, m1 = QSp[1].QNs[iQN[1]].val 
                j2, m2 = QSp[2].QNs[iQN[2]].val 
                j3, m3 = QSp[3].QNs[iQN[3]].val 
                
                j01 = j0 + j1
                cg = 0
                for m01 in range(-j01, j01+1):
                    try:
                        cg_01 =  clebsch_gordan(j0,  j1, j01, m0,  m1, -m01)
                        cg_02 =  clebsch_gordan(j01, j2, j3,  m01, m2, -m3)
                        cg += cg_01 * cg_02
                    except:
                        pass
                match = True if cg else False
                
            else:
                raise NotImplemented 
            if match:
                #print_vars(vars(),  ['cg'])
                d = 1
                for i in range(rank):  #计算某一block的data size 
                    d = d*QSp[i].Dims[iQN[i]]
                
                self.idx[p] = nidx  #给出了0量子数组合与所有量子数组合的序号间的关系 self.idx 和 self.Block_idx[2]互为反函数 如果总量子数为0，则idx[p] =- 1(默认值) 在self.block中都是记录不为0的量子数组合
                self.Block_idx[0, nidx] = totDim  #data block在self.data中的position
                self.Block_idx[1, nidx] = d   #data block变成1d数组的长度
                self.Block_idx[2, nidx] = p   #position in quantum number combinations
                
                #记录不为0的量子数组合，在所有量子数组合中的位置 Addr实为将（QN1, ..., QNn)-> ind 的映射, 将n个指标拉直了, 对每一个iTensor都定义了这个函数 Addr_idx这个二维数组的每一列实际上是所有非零block的量子数的编号(而不是量子数点值！)的组合
                self.Addr_idx[0, nidx] = 0 #for rank=0
                self.Addr_idx[0:rank, nidx] = iQN[0:rank]
                nidx += 1 
                totDim += d   #最终得到self.data 的总长度
                
                #print_vars(vars(),  ['d'])
                
                
                self.cg_coeff[nidx] = cg
            self.transvers_qn_group(iQN,  order)
       
        self.Addr_idx = self.Addr_idx[:, :nidx]
        self.Block_idx = self.Block_idx[:, :nidx]
        self.nidx = nidx
        self.totDim = totDim
        
        self.cg_coeff = self.cg_coeff[:nidx]
            
    def transvers_qn_group(self, iQN, order):
        """
            遍历所有的量子数组合
        
        """
        if order == 'F': 
            inc = True
            i = 0
            while inc and i<self.rank:  #这里实际上意味着按照 fortran order 对量子数组合排序的
                iQN[i] = iQN[i]+1
                if iQN[i]<self.QSp[i].nQN :
                    inc = 0
                else:
                    iQN[i] = 0
                    i = i+1
        elif order == 'C':
            inc = True
            i = self.rank-1
            while inc and i>= 0:
                iQN[i] = iQN[i]+1
                if iQN[i]<self.QSp[i].nQN :
                    inc = False
                else:
                    iQN[i] = 0
                    i = i-1
    
    def __getitem__(self, qn_id_tuple): 
        """
            status:  not fully tested 
            
        """
        pq = 0     
        assert len(qn_id_tuple)==self.rank 
        for i in range(self.rank-1, 0 , -1):
            pq = (pq+qn_id_tuple[i])*self.QSp[i-1].nQN
        pq = pq+qn_id_tuple[0]
        temp = self.idx[pq]
        if temp >= self.idx_dim:
            raise ValueError("error, qn_id_tuple %s not permited. valid values are \n%s "%(
                qn_id_tuple, self.Addr_idx))
        
        sh = [self.QSp[i].Dims[qn_id_tuple[i]] for i in range(self.rank)]
        self.get_block(temp)
        return self.get_block(temp), sh 
    
    def ravel_qn_id_tuple(self, qn_id_tuple):
        """
            from 量子数组合 的多维编号 to 一维编号的映射
            it may be used as follows: 
                temp = self.ravel_qn_id_tuple(qn_id_tuple)
                i = self.idx[temp]
                i is the id of the qn_id_tuple  in self.Addr_idx, self.Block_idx 
                can use it to self.get_block(i),  etc 
        """
        pq = 0     
        assert len(qn_id_tuple)==self.rank 
        for i in range(self.rank-1, 0 , -1):
            pq = (pq+qn_id_tuple[i])*self.QSp[i-1].nQN
        pq = pq+qn_id_tuple[0]
        
        return pq 
    
    def get_block(self, i, linear=True, order='F'): 
        """
            params:
                i: can either be an int or qn_id_tuple 
        """
        if not isinstance(i, int):
            i = self.get_idx(i)
        p = self.Block_idx[0, i]
        size = self.Block_idx[1, i]
        if linear:
            #print_vars(vars(),  ['type(self.data)'])
            return self.data[p: p + size] 
        else:
            sh = self.get_block_shape(i)
            return  self.data[p: p+size].reshape(sh, order=order)            
    
    def get_block_shape(self, i):
        qn_id_tuple = self.Addr_idx[:, i]
        return [self.QSp[i].Dims[qn_id_tuple[i]] for i in range(self.rank)]
    
    def set_block(self, i, data): 
        if not isinstance(i, int):
            i = self.get_idx(i)
        p = self.Block_idx[0, i]
        size = self.Block_idx[1, i]
        #self.data[p: p + size] = data
        self.data[p: p + size] = data[:]  # this checks dtype, the above dont
    
    @staticmethod 
    def example(qsp=None, rank=4, totqn=None, symmetry='Z2', dtype=np.float64, rand_seed=None): 
        """
            convenient method for testing 
        """
        if qsp is not None : 
            pass
        else:
            cls = symmetry_to_Qsp(symmetry)
            if symmetry == 'Z2' : 
                qsp = cls.easy_init( [1, -1], [2, 2]).copy_many(rank)
            elif symmetry == 'U1' : 
                #qsp = cls.easy_init( [0, 1, -1], [2, 1, 1]).copy_many(rank, ) 
                qsp = cls.easy_init( [0, 1, -1], [2, 1, 1]).copy_many(rank, reverse=list(range(rank//2, rank)))  #
            elif symmetry == 'Travial' : 
                qsp = cls.easy_init( [1], [5]).copy_many(rank)
                
        res= iTensor(QSp=qsp, dtype=dtype, totQN=totqn)
        if rand_seed is not None :
            np.random.seed(rand_seed)
        if dtype == np.float64:  
            res.data[:] = np.random.random(res.data.size)-0.5
        elif dtype == np.complex128: 
            res.data[:] = np.random.random(res.data.size)- 1j* np.random.random(res.data.size)
        else: 
            raise ValueError(dtype)
        return res
    
    def unregister_del(self):
        """
        unregister self from iTensor.DATA_BUFFER if iTensor.DATA_BUFFER was used
        """
        
        iTensor.DATA_BUFFER[self.buf_ref[0]].in_use[self.buf_ref[1]]=False
    
    def __del__(self):
        """
            task: replace __del__ with __exit__ in future

            it was said __del__ is BAD!
            但它仍然可以使用，只要__init__ 不包含 raise exception
            "If you do use __del__ make sure you are covered for any case 
            in which __init__ didn’t finish running."
        """
        if self.use_buf and self.buf_ref[0]!=-1:
            if self.use_gpu != 1: 
                iTensor.DATA_BUFFER[self.buf_ref[0]].in_use[self.buf_ref[1]]=False
            else:
                iTensor.DATA_BUFFER_GPU[self.buf_ref[0]].in_use[self.buf_ref[1]]=False

    @property
    def qsp_class(self): 
        return self.QSp[0].__class__
    
    @property 
    def qn_id(self):
        return self.qsp_class.QnClass.qn_id()
            
    @property
    def name(self):
        return self.type_name
    
    @property
    def symmetry(self): 
        return self.qsp_class.QnClass.SYMMETRY 
    
    @property
    def dtype(self):
        return self.data.dtype
    
    @property
    def nbytes(self):
        return self.data.nbytes
    
    
    def memory_use(self):
        a = self.data.nbytes/1e6
        sh = [i.totDim for  i in self.shape]
        b = np.prod(sh)*8/1e6 
        return ('sparse:', a, 'Mbyte\n', 'dense:', b, 'Mbyts')
    
    if 1: #for compatable with numpy
        @property  
        def shape(self): 
            return tuple(self.QSp)   #np.shape is a tuple so convert it 

        @property  #an even shorter name ！
        def sh(self):   
            return tuple(self.QSp)   #np.shape is a tuple so convert it 

        @property
        def size(self): 
            return self.totDim 
        
        @property
        def ndim(self):  
            return self.rank 
        
        #@property   #this is stupid. only temp use 
        def ravel(self): 
            return self.data 
        
        @property
        def T(self): 
            #assert self.rank == 2 
            return self.permutation([1, 0])


    @classmethod
    def buff_free(cls):
        for i in range(3):
            print(cls.DATA_BUFFER[i].in_use[:64])

    def copy_struct(self, other=None, has_data=True, use_buf=False):
        """
            see iTensor_CopyStruct in f90
        """
        QSp= [q.copy() for q in self.QSp]
        totQN = self.totQN.copy()
        other = iTensor(rank=self.rank, QSp=QSp, totQN=totQN, 
                dtype=self.dtype, has_data=has_data, use_buf=use_buf, 
                use_gpu=self.use_gpu)
        other.type_name = self.type_name 
        return other

    def copy_struct_new(self, other=None, index_data=False,  use_buf=False):
        """
            see iTensor_CopyStruct in f90
        """
        QSp= [q.copy() for q in self.QSp]
        totQN = self.totQN.copy()
        other = iTensor(rank=self.rank, QSp=QSp,totQN=totQN, 
                index_data=index_data, has_data=False, use_buf=use_buf)
        return other

    def shallow_copy(self):
        return self

    def copy(self, use_buf=False):
        """
            deepcopy:
                copy also data
        """
        
        other = self.copy_struct(use_buf=use_buf)
        totDim=self.totDim
        #for large array(size>10^4) slicing is faster than copy
        other.data[:totDim]=self.data[:totDim]
        return other

    def copy_new(self, buffer=None, use_buf=False):
        """
        see Tensor_ShallowCopy in f90
        """
        
        other = self.copy_struct(use_buf=use_buf)
        other.set_data_entrance(order="F")
        totDim = self.totDim
        if use_buf:   #use internal DATA_BUFFER; else use external buffer or no buffer
            buffer = other.buffer_assign(data_size=self.totDim, item_size=np.dtype(self.dtype).itemsize)
        other.data = np.ndarray(self.totDim, buffer=buffer, dtype=self.dtype, order="C")
        other.data[:totDim]=self.data[:totDim]
        return other
    
    def __str__(self, keys=None, data_format='ndarray', show_frame=1):
        """
            compare with repr,  this one is simpler 
        """
        rank = self.rank
        if rank == 0:
            rank = 0
        if keys is None:
            #keys=["rank", "type_name", 'ind_labels',  "nidx","totQN","QNs", "Dims","totDim", "Block_idx", "Addr_idx", "data"]
            keys=['rank', 'type_name', 'ind_labels', 
                    'nidx', 'totQN','QNs', 'Dims','totDim', 
                     'data']
        
        str0="----Begin iTensor----------------------------------------------\n"

        str1=""
        str2= '----End iTensor----------------------------------------------\n'
        for k in keys:
            if k == 'data':
                
                if data_format == 'ndarray' and self.size<1000:
                    t = self.to_ndarray()
                    if self.dtype != complex: 
                        temp = '(in matrix view)\n' + str(t.round(5))
                    else:
                        temp = '(in matrix view(real part))\n' + str(t.round(5).real)
                        temp += '\n(imag part)\n' + str(t.round(5).imag)
                #elif data_format == 'none' :
                #    temp = "...."
                else:
                    temp = "\n"
                    for i in range(self.nidx):
                        start = self.Block_idx[0, i]
                        d = self.Block_idx[1, i]
                        qn_id_tuple= self.Addr_idx[:, i]
                        qn_id_linear = self.Block_idx[0, i]
                        #qn_id_linear = self.ravel_qn_id_tuple(qn_id_tuple)
                        #temp += '%d '%i + str(qn_id_tuple) + ": "
                        temp += ''.join(['%d '%i, str(qn_id_tuple), 
                            '*'.join(map(str, self.get_block_shape(i))), 
                            ':']) 
                        if np.all(self.data[start:start + d]==0.0):
                            #temp += "\n all 0.0"
                            temp += ' 0\n' 
                            pass
                        else:
                            
                            if d<= 8: 
                                temp += str(self.data[start:start + d].round(5) ) + "\n"
                            else:
                                temp += str(self.data[start:start + 4].round(5) )  + "..."+ str(self.data[start+d-4:start + d].round(5) ) + " %d "%d + "  only first and last 4 entries" + "\n" 

                #temp = str(self.__dict__[k][:self.totDim].round(10)) 
            elif k == "QNs":
                #temp = [self.QSp[i].QNs[0:self.QSp[i].nQN] for i in xrange(rank)]
                temp = [str(self.QSp[i].QNs[0:self.QSp[i].nQN]) for i in range(rank)]
                #temp  = repr(self.__dict__[k][:rank]) 
                temp = str(temp)
            elif k  == "Dims" :
                temp  = repr(self.__dict__[k][:rank]) 
            else:
                temp = repr(self.__dict__[k])

            str1 += k+":\t" + temp +"\n" 
        
        if show_frame: 
            res=str0+ str1 +str2
        else: 
            #res= ''.join(['\n', str1])
            res= str1 
        return res
    
    def __repr__(self, keys=None, fewer=0, show_frame=1):
        """
        """
        #keys=["rank", "nidx","idx_dim", "Dims","totDim","data"]
        rank = self.rank
        if rank == 0:
            rank = 0
        if keys is None:
            #keys=["rank", "type_name", 'ind_labels',  "nidx","totQN","QNs", "Dims","totDim", "Block_idx", "Addr_idx", "data"]
            keys=["rank", "type_name", 'ind_labels', "nidx","totQN","QNs", "Dims","totDim", "Block_idx", "Addr_idx", "data"]

        str0="----Begin iTensor----------------------------------------------\n"

        str1=""
        str2= '----End iTensor----------------------------------------------\n'
        for k in keys:
            if k == 'data':
                if fewer:
                    temp = "...."
                else:
                    temp = "\n"
                    for i in range(self.nidx):
                        start = self.Block_idx[0, i]
                        d = self.Block_idx[1, i]
                        qn_id_tuple = self.Addr_idx[:, i]
                        qn_tuple = self.get_qn_from_qnid(qn_id_tuple)
                        qn_tuple = np.asarray([i.val for i in qn_tuple])
                        qn_id_linear = self.Block_idx[0, i]
                        #temp += ''.join(['%d '%i, str(qn_id_tuple), 
                        temp += ''.join(['%d '%i, str(qn_tuple) + ' ', 
                            '*'.join(map(str, self.get_block_shape(i))), 
                            ':']) 
                        if np.all(self.data[start:start + d]==0.0):
                            #temp += "\n all 0.0"
                            temp += ' 0\n' 
                            pass
                        else:
                            
                            if d<= 8: 
                                temp += str(self.data[start:start + d].round(5) ) + "\n"
                            else:
                                temp += str(self.data[start:start + 4].round(5) )  + "..."+ str(self.data[start+d-4:start + d].round(5) ) + " %d "%d + "  only first and last 4 entries" + "\n" 

                #temp = str(self.__dict__[k][:self.totDim].round(10)) 
            elif k == "QNs":
                #temp = [self.QSp[i].QNs[0:self.QSp[i].nQN] for i in xrange(rank)]
                temp = [str(self.QSp[i].QNs[0:self.QSp[i].nQN]) for i in range(rank)]
                #temp  = repr(self.__dict__[k][:rank]) 
                temp = str(temp)
            elif k  == "Dims" :
                pass 
                temp  = repr(self.__dict__[k][:rank]) 
            elif k == "idx":
                temp = repr(self.__dict__[k])
            elif k == "Block_idx":
                temp = '\n\t' + str(self.Block_idx[0, :self.nidx]) + "\n\t" + str(self.Block_idx[1, :self.nidx]) + "\n\t" + str(self.Block_idx[2, :self.nidx]) + "\n" 
                #temp = str(self.Block_idx[:, 0:self.nidx])
            elif k == "Addr_idx":
                temp = '\n'+str(self.__dict__[k][:rank, :self.nidx]) 
            else:
                temp = repr(self.__dict__[k])

            str1 += k+":\t" + temp +"\n" 
        
        if show_frame: 
            res=str0+ str1 +str2
        else: 
            #res= ''.join(['\n', str1])
            res= str1 
        return res
    
    def show_data(self): 
        print(self.__repr__(keys=['data'], fewer=0, show_frame=0))
    
    def show_struct(self): 
        print_vars(vars(), ['self.Addr_idx[:,:self.nidx].T', 'self.Block_idx[:,:self.nidx].T', ], 
                sep = '\n', head = 'result of show_struct:\n', 
                key_val_sep='=')
        
    @staticmethod
    def unit_tensor(rank, QSp, totQN=None, dtype=np.float64, use_gpu=0):
        """
            Q:  注意区分几种情况，
            
               --->--.-->---
                
               -->--|      
                    |
               -->--|
               
               |-->--           
               |                 
               |-->--
               
           后面两种叫单位张量吗？ 
        
        """
        totQN = totQN if totQN is not None else QSp[0].QnClass.qn_id()
        if rank%2 != 0:
            raise ValueError("rank shold be even, rank=%s"%(rank, ))
        t = iTensor(rank, QSp, totQN, dtype=dtype, use_gpu=use_gpu) 
        pTot = 1
        
        rank1 =t.rank//2 #use t.rank/2 would yeild a float 2.0
        
        for i in range(rank1):
            pTot = pTot*t.QSp[i].nQN
        #这里rank1只有总rank的一半，pTot! = self.idx_dim

        t.data[:] = 0.0    
        pos= np.ndarray(rank1, "int")
        for i in range(pTot):
            p = i+i*pTot  #this is diagonal line
            idx=t.idx[p]
            p = t.Block_idx[0,idx]
            pos[0:rank1] = t.Addr_idx[0:rank1, idx]
            d = 1
            for j in range(rank1):
                d = d*t.QSp[j].Dims[pos[j]]
            if use_gpu != 1:
                t.data[p:p+d**2] = np.identity(d,dtype=t.dtype).ravel()
            else:
                t.data[p:p+d**2] = cp.identity(d,dtype=t.dtype).ravel()
                
        return t
    
    @staticmethod
    def identity(qsp, dtype=np.float64, use_gpu=0):
        if hasattr(qsp, 'QNs'):
            qsp = qsp.copy_many(2, reverse=[1])
        return iTensor.unit_tensor(2, qsp, dtype=dtype, use_gpu=use_gpu)
    
    def zeros(qsp, dtype=np.float64, totQN=None):
        res = iTensor(QSp=qsp, dtype=dtype, totQN=totQN)
        res.data[:] = 0.0
        return res 
    
    def get_position(self, qn_ind_tuple):
        """
            map a rank-dim index to linear index 
            
            this func is used in following way: 
                p3=T3.get_position(iQN3[:T3.rank])
                idx3 = T3.idx[p3]
                p3 = T3.Block_idx[0,idx3]
                data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
            params: 
                qn_ind_tuple: 量子数编号的组合
                量子数的编号也是用的 前小后大的存储顺序
        """
        if len(qn_ind_tuple)==0: #this is for rank 0 tensor
            return 0
        
        dims=tuple(self.QSp[i].nQN for i in range(self.rank)) 
        p = np.ravel_multi_index(qn_ind_tuple, dims, order='F')
        
        return p
    
    def get_idx(self, qn_id_tuple):
        p = self.get_position(qn_id_tuple)
        if p >= self.idx_dim:  # wrong qn_id_tuple
            return -1 
        i = self.idx[p]
        return int(i)   #convert np.int64 to int 
    
    def get_qnid_from_qn(self, qn_tuple):
        """
            map a tuple of qns to a tuple of their ids
        """
        try:
            return tuple(self.QSp[i].QNs.index(qn_tuple[i]) for i in range(self.rank))
        except ValueError:
            return None
    
    def get_qn_from_qnid(self, qn_id_tuple):
        return tuple(self.QSp[i].QNs[qn_id_tuple[i]] for i in range(self.rank))
        
    
    def get_position_rev(self, index_linear):
        """ 
            map a ind_linear to a tuple 
            
            Get Position in T, get iQN from pos
        
        """
        dims=tuple(self.QSp[i].nQN for i in range(self.rank)) 
        pos = np.unravel_index(index_linear, dims, order='F')
            
        return pos

    def set_element(self, qn_id_tuple, sub_ind, element):
        """
            see iTensor_SetElement in f90
            params: 
                qn_id_tuple: block coordinate, len(qn_id_tuple)=self.rank
                    for each i, qn_id_tuple[i] in xrange(self.Qsp[i].nQN),  
                    对于Z2 symm 即 [0, 1] 
                    note 这不是量子数的val,  but numbering  
                    并且这些量子数fuse到一起后 = self.totQN
                sub_ind：data coodinate in the block,  len(sub_ind)=self.rank
                    for each i, sub_ind[i] in xrange(self.Qsp[i].Dims), 
                    对于Z2 symm 也是[0, 1] for any i
            #整个的思路是，先找到block entry，再找data position
        """
        pq = 0; pi=0
        for i in range(self.rank-1, 0 , -1):
            #print "iiiii i", i
            pq = (pq+qn_id_tuple[i])*self.QSp[i-1].nQN
            pi = (pi+sub_ind[i])*self.QSp[i-1].Dims[qn_id_tuple[i-1]]
            #print self.QSp[i-1].Dims[qn_id_tuple[i-1]]
        
        pq = pq+qn_id_tuple[0]
        pi = pi+sub_ind[0]
        temp = self.idx[pq]
        if temp >= self.idx_dim:
            raise ValueError("error, qn_id_tuple values not permited. qn_id_tuple=%s"%(qn_id_tuple, ))
        pidx = self.Block_idx[0, temp]
        block_len = self.Block_idx[1, temp]
        if pi>= block_len:
            raise ValueError("error, sub_ind values not permited, sub_ind=%s"%sub_ind)
        #print "ppppqqqq  ", pq, pi,'idx', temp, 'pidx', pidx, pidx + pi 
        #if isinstance(element, complex):
        #    raise  
        self.data[pidx+pi] = element 

    def get_element(self, qn_id_tuple, sub_ind):
        """
            params:
                qn_id_tuple: a tuple of length self.rank 
                    block coordinate
                    for each i, qn_id_tuple[i] in xrange(self.Qsp[i].nQN),  
                    对于Z2 symm 即 [0, 1] 
                    note 这不是量子数的val,  but numbering  
                    并且这些量子数fuse到一起后 = self.totQN
                sub_ind：data coodinate in the block,  len(sub_ind)=self.rank
                    for each i, sub_ind[i] in range(self.Qsp[i].Dims), 
                    对于Z2 symm 也是[0, 1] for any i
        """
        pq = 0; pi=0
        for i in range(self.rank - 1, 0, -1):
            pq = (pq+qn_id_tuple[i])*self.QSp[i-1].nQN
            pi = (pi+sub_ind[i])*self.QSp[i-1].Dims[qn_id_tuple[i-1]]
        
        pq = pq+qn_id_tuple[0]
        pi = pi+sub_ind[0]  
        idx = self.idx[pq]
        #print_vars(vars(),  ['idx'])
        pidx = self.Block_idx[0, idx]
        idx = self.idx[pq]
        if idx >= self.idx_dim:
            print("error, qn_id_tuple values not permited")
            return 
        pidx = self.Block_idx[0, idx]
        block_len = self.Block_idx[1, idx]
        if pi>= block_len:
            print("error,  sub_ind values not permited")
            return 
        
        X = self.data[pidx+pi]
        return X

    class ShapeError(Exception):
        """
        only a tiny test of exception
        """
        pass

    @staticmethod 
    def is_match(qq2, dd2, qq1, dd1, info=1): 
        """
            status: 
                not completed

        """
        res = True 
        len2 = len(qq2)
        len1 = len(qq1)
        
        match = [0 for i in range(len(qq2))]
        i = 0
        for a, q2 in enumerate(qq2): 
            if i >= len1: 
                return False 
            x = qq1[i]
            y = dd1[i]
            j = i 
            for b in range(j, len1): 
                i += 1
                if x.val == q2.val and y <= dd2[a] : 
                    if info>0: 
                        #print a, i-1, x, y, q2   
                        print_vars(vars(), ['a', 'i-1', 'x.val', 'y', 'q2.val'], sep=' ')
                    #if a != len2-1 or b == len1-1 : 
                        match[a] = 1
                        break 
                else:
                #if 1: 
                    if b  == len1-1: 
                        return False 
                    x = x + qq1[b+1]
                    y = y*dd1[b + 1]
        if info>0: print(match) 
        res= all(match)                            
        return res

    @staticmethod 
    def is_match_new(qq2, dd2, qq1, dd1, info=1): 
        """
             status: 
                not completed
       
        """ 
        res = True 
        len2 = len(qq2)
        len1 = len(qq1)
        match = [0 for i in range(len(qq2))]
        
        id1 = 0
        id2 = 0
        iq1 = 0
        iq2 = 0
        for i in range(max(len1, len2)): 
            d1 = dd1[id1]
            d2 = dd2[id2]
            q1 = qq1[iq1]
            q2 = qq2[iq2]
            if d1 == d2: 
                if q1.val != q2.val: 
                    return False 
            elif d1>d2: 
                if q1.val == q2.val: 
                    pass
                else: 
                    pass
            elif d1<d2: 
                if q1.val == q2.val: 
                    pass
                else: 
                    pass 
        if info>0: print(match) 
        res= all(match)                            
        return res

    def reshape_general_but_not_realizable(self, qsp_new): 
        """
            status: 
                not completed
            re-indexing tensor like np.reshape 
            index (or qsp) can be joined or splited 
            this function is one of the most non-travial one, many improvements needed 
            index match 应该是个（线性）方程组问题，不能局部地inspect
            howto: 
                first reshape the qn second inner block 
                最终是要得到一个映射，从self to other Addr_idx ，Block_idx, 之间的映射关系
                    如果有merge/split的guidance, Addr_idx 的映射是容易确定的, 满足量子数的约束即可
                
            
        """
        qsp = self.QSp
        q = qsp[0]
        for i in qsp[1: ]:  
            q = q.add(i)
        p = qsp_new[0]
        for i in qsp_new[1: ]:  
            p = p.add(i)
        assert q == p  
        
        qn_cls = qsp[0].QnClass
        res = iTensor(QSp=qsp_new, totQN=self.totQN)
        
        t1 = self
        t2 = res 
        rank1 = t1.rank 
        rank2 = t2.rank 
        #match_table if to find out corespondence between blocks 
        match_table = {i: [] for i in range(t2.nidx)}
        
        is_match = iTensor.is_match   
        for i in range(t2.nidx): 
            qn_id_tuple_2 = t2.Addr_idx[:, i]
            qn_tuple_2 = [t2.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_2)]
            dim_tuple_2 = [t2.QSp[i2].Dims[q] for i2, q in enumerate(qn_id_tuple_2) ]
            #print  qn_tuple_2 
            for j in range(t1.nidx): 
                qn_id_tuple_1 = t1.Addr_idx[:, j]
                qn_tuple_1 = [t1.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_1)]
                dim_tuple_1 = [t1.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_1) ]
                #if i == 1 and is_match(qn_tuple_2, qn_tuple_1):  
                #if i == 1 and j == 1 : 
                if i == 0 and j == 0 : 
                #if 1: 
                    print(i, j) 
                    print(qn_tuple_2 , dim_tuple_2)
                    print(qn_tuple_1 , dim_tuple_1) 
                    print(is_match(qn_tuple_2, dim_tuple_2,  qn_tuple_1, dim_tuple_1))
                        
                if is_match(qn_tuple_2, dim_tuple_2,  qn_tuple_1, dim_tuple_1): 
                    match_table[i].append(j)
        
        print(match_table)     
        return res    
    
    def reshape_bac(self, *qsp_new): 
        """
            since a most general purpose reshape is difficult to realize, 
            here use merge_qsp/split_qsp instead. 
            this is only for convenience and a little bit slow. 
            for production use .merge_qsp and .split_qsp directly 
            note:
                运行以下code
                    from tensor import iTensorFactory 
                    pau = pauli_mat()
                    #sx = pau['sx']
                    sx = np.asarray([[0, 1], [1, 0]])
                    sxx=np.multiply.outer(sx, sx).transpose([0, 2, 1, 3]).reshape(4, 4)
                    print sxx 
                    sx = iTensorFactory.pauli_mat('Z2')['sigma_x']
                    sxx=sx.direct_product(sx) 
                    #先reshape，后to_ndarray
                    print sxx.merge_qsp((0, 1), (2, 3)).to_ndarray(data_order='F').round(5)
                    #先to_ndarray，后reshape 
                    print sxx.to_ndarray(data_order='F').reshape(4, 4)
               输出为 
                    [[0 0 0 1]
                     [0 0 1 0]
                     [0 1 0 0]
                     [1 0 0 0]]
                    [[ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]
                     [ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]]
                    [[ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]
                     [ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]]
               之所以不同，是因为，对称张量的reshape 和 ndarray的reshape 名字相同，但
               实质有差别， iTensor.reshape 包含了将量子数指标的合并 
               
        """
        #print_vars(vars(), ['self.QSp', 'qsp_new'])
        if hasattr(qsp_new[0], '__iter__'): 
            qsp_new = qsp_new[0]
    
        if self.rank>len(qsp_new):   #merge 
            which = 'merge'
            qsp2 = qsp_new; qsp3 = self.QSp
        elif self.rank<len(qsp_new):  #split
            which = 'split' 
            qsp3 = qsp_new; qsp2 = self.QSp
        else: 
            raise ValueError('I dont know it is spilt or merge, %d, %d'%(self.rank, len(qsp_new)))
        
        #inspect leg_map
        leg_map = {}
        #merged_leg_list = []
        ii = 0
        for i, q in enumerate(qsp2): 
            leg_map[i] = ()
            temp =  self.qsp_class.null()
            count = 0
            for j, p in enumerate(qsp3[ii: ]): 
                
                temp = temp*p 
                if temp < q :
                    leg_map[i] += (ii + j, )
                    count += 1  
                elif temp  == q: 
                    leg_map[i] += (ii + j, )
                    ii += count + 1  
                    #if len(leg_map[i])>1: 
                    #    merged_leg_list.append(leg_map[i])
                    break 
                else: 
                    raise ValueError('cant be reshaped, check qsp.\n\tself.qsp=%s\n\tqsp_new=%s'%(self.QSp, qsp_new)) 
 
        #如果 p == null leg_map 有可能判断错误，此情况下, 用下面几行补救 
        last_leg_3a = leg_map[len(qsp2)-1][-1]
        last_leg_3b = len(qsp3)-1
        if last_leg_3a  == last_leg_3b:  #normal case 
            pass 
        elif last_leg_3a  + 1 == last_leg_3b: 
            if qsp3[-1] ==  self.qsp_class.null(): 
                leg_map[len(qsp2)-1] += (last_leg_3b, )
            else: 
                raise Exception("inspecting leg_map failed")
        else: 
            raise Exception("inspecting leg_map failed")
        
        if which == 'merge':
            arg = []
            for i in sorted(leg_map):
                if len(leg_map[i])>1: 
                    arg.append(leg_map[i])
            return self.merge_qsp(*tuple(arg))
        elif which == 'split' : 
            arg = []
            for k, v in leg_map.items(): 
                if len(v)>1: 
                    arg.append(k)
                    arg.append([qsp3[i] for i in v])
            return self.split_qsp(*tuple(arg))
        else: 
            raise ValueError('I dont know it is spilt or merge')
    
    def reshape(self, *qsp_new): 
        """
            since a most general purpose reshape is difficult to realize, 
            here use merge_qsp/split_qsp instead. 
            this is only for convenience and a little bit slow. 
            for production use .merge_qsp and .split_qsp directly 
            note:
                1. 
                运行以下code
                    from tensor import iTensorFactory 
                    pau = pauli_mat()
                    #sx = pau['sx']
                    sx = np.asarray([[0, 1], [1, 0]])
                    sxx=np.multiply.outer(sx, sx).transpose([0, 2, 1, 3]).reshape(4, 4)
                    print sxx 
                    sx = iTensorFactory.pauli_mat('Z2')['sigma_x']
                    sxx=sx.direct_product(sx) 
                    #先reshape，后to_ndarray
                    print sxx.merge_qsp((0, 1), (2, 3)).to_ndarray(data_order='F').round(5)
                    #先to_ndarray，后reshape 
                    print sxx.to_ndarray(data_order='F').reshape(4, 4)
               输出为 
                    [[0 0 0 1]
                     [0 0 1 0]
                     [0 1 0 0]
                     [1 0 0 0]]
                    [[ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]
                     [ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]]
                    [[ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]
                     [ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]]
               之所以不同，是因为，对称张量的reshape 和 ndarray的reshape 名字相同，但
               实质有差别， iTensor.reshape 包含了将量子数指标的合并 
               2. 
                    see the note in the doc string of merge_qsp 
               
        """
        #print_vars(vars(), ['self.QSp', 'qsp_new'])
        if hasattr(qsp_new[0], '__iter__'): 
            qsp_new = qsp_new[0]
    
        if self.rank>len(qsp_new):   #merge 
            which = 'merge'
            qsp2 = qsp_new; qsp3 = self.QSp
        elif self.rank<len(qsp_new):  #split
            which = 'split' 
            qsp3 = qsp_new; qsp2 = self.QSp
        else: 
            raise ValueError('I dont know it is spilt or merge, %d, %d'%(self.rank, len(qsp_new)))
        
        #inspect leg_map
        leg_map = {}
        #merged_leg_list = []
        ii = 0
        for i, q in enumerate(qsp2): 
            leg_map[i] = ()
            temp =  self.qsp_class.null()
            count = 0
            for j, p in enumerate(qsp3[ii: ]): 
                
                #temp = temp*p 
                temp = temp.__mul__(p)
                if temp.totDim < q.totDim :
                    leg_map[i] += (ii + j, )
                    count += 1  
                elif temp == q: 
                    leg_map[i] += (ii + j, )
                    ii += count + 1  
                    #if len(leg_map[i])>1: 
                    #    merged_leg_list.append(leg_map[i])
                    break 
                else: 
                    raise ValueError('cant be reshaped, check qsp.\n\tself.qsp=%s\n\tqsp_new=%s'%(self.QSp, qsp_new)) 
 
        #如果 p == null leg_map 有可能判断错误，此情况下, 用下面几行补救 
        last_leg_3a = leg_map[len(qsp2)-1][-1]
        last_leg_3b = len(qsp3)-1
        if last_leg_3a  == last_leg_3b:  #normal case 
            pass 
        elif last_leg_3a  + 1 == last_leg_3b: 
            if qsp3[-1] ==  self.qsp_class.null(): 
                leg_map[len(qsp2)-1] += (last_leg_3b, )
            else: 
                raise Exception("inspecting leg_map failed, reshape function is not perfect yet")
        else: 
            raise Exception("inspecting leg_map failed, reshape function is not perfect yet")
        
        if which == 'merge':
            arg = []
            for i in sorted(leg_map):
                if len(leg_map[i])>1: 
                    arg.append(leg_map[i])
            return self.merge_qsp(*tuple(arg))
        elif which == 'split' : 
            arg = []
            for k, v in leg_map.items(): 
                if len(v)>1: 
                    arg.append(k)
                    arg.append([qsp3[i] for i in v])
            return self.split_qsp(*tuple(arg))
        else: 
            raise ValueError('I dont know it is spilt or merge')
    
    def index_merge_1(self, ind):
        """
            in numpy reshape serves as merge and split of indices of a dense array, I take use of that.
        """
        n1, n2 = min(ind), max(ind)
        ind1 = range(self.rank)  #indices not changed
        for i in ind: ind1.remove(i)  #indices to be merged
        ind2 = ind
        in1 = [i for i in ind1 if i <min(ind) ] + [i-1 for i in ind1 if i> max(ind)]
        in2 = min(ind)
        print(ind1, ind2, in1, in2)
        
        qsp_  = np.array([q.copy() for q in self.QSp], object)
        qsp = np.ndarray(self.rank-1, object)
        qsp[in1] = qsp_[ind1]
        qsp[in2] = qsp_[ind2[0]].add(qsp_[ind2[1]])
        qsp = qsp.tolist()
        totqn = self.totQN.copy()
        res= iTensor(self.rank-1, qsp, totqn)
        
        for i in range(res.nidx):
            iq = res.Addr_idx[:, i]
            iq1 = res.Addr_idx[in1, i]
            iq2 = res.Addr_idx[in2, i]
            sha = [res.QSp[x].Dims[iq[x]]  for x in range(res.rank)]
            p2 = res.Block_idx[0, i]
            size2 = res.Block_idx[1, i]
            stack = []
            for j in range(self.nidx):
                iqn = self.Addr_idx[:, j]
                
                if np.all(iqn[ind1]==iq[in1]):
                    shape2 = [self.QSp[x].Dims[iqn[x]]  for x in ind2] 
                    sha_ = [x for x in sha]
                    sha_[in2] = np.prod(shape2)
                    p = self.Block_idx[0, j]
                    size = self.Block_idx[1, j]
                    temp = self.data[p:p + size].reshape(sha_) 
                    stack.append(temp)
                    print("ii", i, "jj", j, end=' ') 
                    print(iqn, iq, sha, sha_)
            xx = np.concatenate(stack)
            res.data[p2:p2 + size2] = xx.ravel()

        return res
    
    def index_merge(self, ind):
        qsp_n = [self.QSp[i].copy() for i in ind]
        for q in qsp_n: q.reverse()
        qsp_1 = qsp_n[0].copy()
        for q in qsp_n[1:]:
            qsp_1 = qsp_1.add(q)
        qsp_1.reverse()
        #r = self.rank-len(ind) + 1
        r = len(ind) + 1
        qsp = qsp_n + [qsp_1]
        merger = iTensor(r, qsp, self.totQN.copy())
        merger.data[:] = 0.0
        
        for i in range(merger.nidx):
            p = merger.Block_idx[0, i]
            l = merger.Block_idx[1, i]
            iqn = merger.Addr_idx[:, i]
            shape = [merger.QSp[x].Dims[iqn[x]] for x in range(merger.rank)]
            data = merger.data[p:p + l].reshape(shape, order="C")
            print("iii", shape)
            if 1:
                for i1 in range(shape[0]):
                    for i2 in range(shape[1]):
                        #i3 = i1 + shape[1]*i2
                        i3 = i2 + shape[1]*i1
                        data[i1, i2, i3] = 1.0 
        
        print(merger)
        res, nothing= self.contract(merger, range(self.rank), ind + [self.rank + 100])
        return res
    
    def index_merge_simple(self): 
        """
            I dont know how to implement merge, so I first make a simple one, 
            change a type (2, 2) tensor to type (1, 1) only for Z2 symm. Use it
            for combine two site in simulation 
            
            NOTE: only suit for Z2 symm. a narive implementation
        """
        
        assert self.rank == 4
        #assert self.QSp[0].QnClass.SYMMETRY == 'Z2' 
        if 0: 
            qsp = QspZ2.easy_init([1,-1], [1, 1 ]) 
            qsp=qsp.copy_many(4, reverse=[2, 3])
            #qn=QnZ2.qn_id()
            tqn = self.totQN.copy()
            res = iTensor(rank, qsp, tqn) 
        r_orig = self.rank
        lower_ind_orig = [0, 1]
        upper_ind_orig = [2, 3]
        
        if 1: 
            qsp_lower_list = [q.copy() for q in self.QSp[: 2]]
            qsp_lower = qsp_lower_list[0]
            for q in qsp_lower_list[1: ]: 
                qsp_lower = qsp_lower.tensor_prod(q)
            qsp_upper_list = [q.copy() for q in self.QSp[2: r_orig]]    
            qsp_upper = qsp_upper_list[0]
            for q in qsp_upper_list[1: ]: 
                qsp_upper = qsp_upper.tensor_prod(q)
            qsp = [qsp_lower, qsp_upper]
            qn = self.totQN.copy()
        res = iTensor(2, qsp, qn) 
        
        for j in range(self.nidx): 
            
            lower_qn_list = [self.QSp[x].QNs[self.Addr_idx[x, j]] for x in lower_ind_orig]
            upper_qn_list = [self.QSp[x].QNs[self.Addr_idx[x, j]] for x in upper_ind_orig]
            lower_qn = qsp[0].QnClass.sum(lower_qn_list)
            upper_qn = qsp[0].QnClass.sum(upper_qn_list)
            
            lower_qn_id = res.QSp[0].QNs.index(lower_qn)
            upper_qn_id = res.QSp[1].QNs.index(upper_qn)
                 
            qDims = [lower_qn_id, upper_qn_id]
            def temp(): 
                mapper = {(0, 0): 0, (1, 1): 1, (0, 1):0, (1, 0):1}
                aaa = self.Addr_idx[:, j]
                a = mapper[(aaa[0], aaa[1])]
                b = mapper[(aaa[2], aaa[3])]
                return [a, b]
            p = self.Block_idx[0, j]
            x = self.data[p]
            res.set_element(qDims=qDims, iDims=temp(), X=x)           
        #print ' -- mmerge -- '       
        return  res
   
    def index_merge_simple_2(self): 
        """
            I dont know how to implement merge, so I first make a simple one, 
            change a type (2, 2) tensor to type (1, 1) only for Z2 symm. Use it
            for combine two site in simulation 
        """
        raise NotImplemented('this is not correct')
        assert self.rank == 4
        #assert self.QSp[0].QnClass.SYMMETRY == 'Z2' 
        if 0: 
            qsp = QspZ2.easy_init([1,-1], [1, 1 ]) 
            qsp=qsp.copy_many(4, reverse=[2, 3])
            #qn=QnZ2.qn_id()
            tqn = self.totQN.copy()
            res = iTensor(rank, qsp, tqn) 
        r_orig = self.rank
        lower_ind_orig = [0, 1]
        upper_ind_orig = [2, 3]
        
        if 1: 
            qsp_lower_list = [q.copy() for q in self.QSp[: 2]]
            qsp_lower = qsp_lower_list[0]
            for q in qsp_lower_list[1: ]: 
                qsp_lower = qsp_lower.tensor_prod(q)
            qsp_upper_list = [q.copy() for q in self.QSp[2: r_orig]]    
            qsp_upper = qsp_upper_list[0]
            for q in qsp_upper_list[1: ]: 
                qsp_upper = qsp_upper.tensor_prod(q)
            qsp = [qsp_lower, qsp_upper]
            qn = self.totQN.copy()
        res = iTensor(2, qsp, qn) 
        
        if 0: 
            #先对齐block        
            for i in range(self.nidx): 
                p = self.Block_idx[0, i]
                l = self.Block_idx[1, i]
                data = self.data[p: p + l]
                print(data)
                #iqn = self.Addr_idx[:, i]
                lower_qn_list = [self.QSp[j].QNs[self.Addr_idx[j, i]] for j in lower_ind_orig]
                upper_qn_list = [self.QSp[j].QNs[self.Addr_idx[j, i]] for j in upper_ind_orig]
                lower_qn = qsp[0].QnClass.sum(lower_qn_list)
                upper_qn = qsp[0].QnClass.sum(upper_qn_list)
                
                lower_qn_id = res.QSp[0].QNs.index(lower_qn)
                upper_qn_id = res.QSp[1].QNs.index(upper_qn)
                qn_id = [lower_qn_id, upper_qn_id]
                temp = res.get_position(qn_id)
                pos = res.idx[temp]
                print(pos)
                #res.data[pos]
        if 1: 
            
            for i in range(res.nidx): 
                qn_id = res.Addr_idx[:, i]
                block_start = res.Block_idx[0, i]
                pointer = block_start
                
                for j in range(self.nidx): 
                    
                    lower_qn_list = [self.QSp[x].QNs[self.Addr_idx[x, j]] for x in lower_ind_orig]
                    upper_qn_list = [self.QSp[x].QNs[self.Addr_idx[x, j]] for x in upper_ind_orig]
                    lower_qn = qsp[0].QnClass.sum(lower_qn_list)
                    upper_qn = qsp[0].QnClass.sum(upper_qn_list)
                    
                    lower_qn_id = res.QSp[0].QNs.index(lower_qn)
                    upper_qn_id = res.QSp[1].QNs.index(upper_qn)
                    #print 'iiii', i, j
                    if qn_id[0]  == lower_qn_id and qn_id[1] == upper_qn_id: 
                        block_start_orig = self.Block_idx[0, j]
                        l = self.Block_idx[1, j]
                        data = self.data[block_start_orig: block_start_orig + l]
                        res.data[pointer: pointer + l] = data
                        pointer += l 
                        
        #print '--mmerge-2 ---'
        return  res

    def split_qsp(self, *args):
        """
            example: 
                t.split_qsp(1, [qb, qc], 3, [qe, qf, qg])
            note: 
                see the note in doc string of merge_qsp
        """
        legs_to_split = args[0::2]
        qq = args[1::2]  # a list of qsp list 
        qsp = []
        qsp_class= self.qsp_class
        qn_class= qsp_class.QnClass 
        ii = 0   # ii points to leg of t3
        l_prev = -1 
        # leg_map = {x:(x, ) for x in xrange(self.rank)}
        leg_map = {}
        for i, l in enumerate(legs_to_split): 
            if l-l_prev>1:  #legs between (l_prev, l) are not splited
                qsp.extend(self.QSp[l_prev+1:l])
                for x in range(1, l-l_prev): 
                    leg_map[l_prev+x] = (x+ii-1, )
                ii +=  l-l_prev-1
                    
            temp = qq[i]
            qsp.extend(temp)
            assert self.QSp[l] == qsp_class.prod_many(temp), (self.QSp[l],  qsp_class.prod_many(temp)) 
            length = len(temp)
            leg_map[l] = range(ii, ii + length)
            ii +=  length
            l_prev = l 
        qsp.extend(self.QSp[legs_to_split[-1]+1: ])
        for x in range(legs_to_split[-1]+1, self.rank):
            leg_map[x] = (ii, )
            ii += 1  
        
        #print_vars(vars(), ['leg_map', 'len(qsp)']); raise 
        
        res = iTensor(QSp=qsp, totQN=self.totQN.copy(), 
                dtype=self.dtype, use_gpu=self.use_gpu)
        #print_vars(vars(),  ['self.data.dtype'], '', '  ')
        t2 = self
        t3 = res 
        t3ind = list(range(t3.nidx))
        for i in range(t2.nidx): 
            qn_id_tuple_2 = t2.Addr_idx[:t2.rank, i]
            qn_tuple_2 = [t2.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_2)]
            dim_tuple_2 = [t2.QSp[i2].Dims[q] for i2, q in enumerate(qn_id_tuple_2) ]
            data_block_2 = t2.get_block(i)
            data_block_2 = data_block_2.reshape(dim_tuple_2, order='F')
            
            #block_shape = [1 if x not in legs_to_split else len(leg_map[x]) for x in range(t2.rank)]
            #print_vars(vars(), ['qn_id_tuple_2', 'dim_tuple_2'])
            start_dicts= {x: {} for x in range(t2.rank) }
            
            #for j in range(t3.nidx):
            matched_j_list = []
            for j in t3ind: 
                qn_id_tuple_3 = t3.Addr_idx[:t3.rank, j]
                qn_tuple_3 = [t3.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_3)]
                match = True 
                #for l in legs_to_split: 
                for l in leg_map: 
                    if qn_tuple_2[l] != qn_class.sum([qn_tuple_3[_i] for _i in leg_map[l]]): 
                        match = False 
                        break 
                if match:         
                    dim_tuple_3 = [t3.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_3) ]
                    dim_delta  = np.zeros(t2.rank, dtype=int)
                    
                    qn_id_tuple_grouped = [tuple([qn_id_tuple_3[y] for y in leg_map[x] ]) for x in range(t2.rank)]
                    dim_tuple_grouped = [np.prod([dim_tuple_3[y] for y in leg_map[x] ]) for x in range(t2.rank)]
                    
                    sh_start = np.zeros(t2.rank, dtype=int)
                    for k, v in enumerate(qn_id_tuple_grouped): 
                        temp = start_dicts[k]
                        if v in temp: 
                           sh_start[k] = temp[v]
                        else:
                            if 'totdim' in temp: 
                                sh_start[k] = temp['totdim']
                                temp[v] = temp['totdim']
                                temp['totdim'] += dim_tuple_grouped[k] 
                            else:
                                sh_start[k] = 0
                                temp[v] = 0
                                temp['totdim'] =  dim_tuple_grouped[k] 

                    
                    dim_delta = np.asarray(dim_tuple_grouped, dtype=int)
                    sh_end = sh_start + dim_delta 
                    #print_vars(vars(), ['i', 'j', 'qn_id_tuple_3', 'qn_id_tuple_grouped', 
                    #    'dim_tuple_grouped',  'sh_start', 'sh_end'], sep=', ')
                    sl = [slice(sh_start[iii], sh_end[iii]) for iii in range(t2.rank)]
                    data = data_block_2[tuple(sl)].ravel(order='F')
                    t3.set_block(j, data)
            
            for j in matched_j_list: 
                t3ind.remove(j)
                    
        return  t3         
    
    def merge_qsp(self, *args):
        """
            args should be list tuples of coningueous legs
            note: 
                the merge/split of qsp may be tricky at some points.
                several ATTENTION follows 
                1. sometimes, one need conj_new a leg,  before reshape it 
                2. reshape may change the memory of data
                    ***A REMARKABLE fact is that, the final order of data depend
                    on the course of reshape***: starting from a tensor t,  if
                    reshape it se several times to a tensor t1; and reshape t
                    several times to a t2 such that t1.qsp = t2.qsp. if the
                    course of the reshape are different, then the memory of
                    data of t1 and t2 may differ !!
                    
                    it seems that, in contrast, dense tensor dont depend on the course 
                    merge qsp (final qsp 更粗) 容易出此问题，split qsp （更细) 可能无此问题。
                    
                    #code example     
                        np.random.seed(1234)
                        q = QspU1.easy_init([1, -1], [1, 1])
                        qsp = q.copy_many(4) 
                        
                        t4 = iTensor.example(qsp=qsp)    # starting from a rank-4 tensor 
                        t_13 = t4.merge_qsp((0, ), (1, 2, 3))
                        t4_1 = t4.merge_qsp_all()
                        t_13_1 = t_13.merge_qsp_all()
                        t_13_1_4 = t_13_1.split_qsp(0, q.copy_many(4))
                        
                        print_vars(vars(),  ['t4.data', 't_13.data', 't4_1.data', 't_13_1.data', 
                            't_13_1_4.data' , 't4_1.shape', 't_13_1.shape' ])
                        
                    #output 
                         #t4.data=[-0.30848055  0.12210877 -0.06227226  0.28535858  0.27997581 -0.22740739]
                         #t_13.data=[-0.30848055  0.12210877  0.28535858 -0.06227226  0.27997581 -0.22740739]
                         #t4_1.data=[-0.30848055  0.12210877 -0.06227226  0.28535858  0.27997581 -0.22740739]
                         #t_13_1.data=[-0.30848055  0.12210877  0.28535858 -0.06227226  0.27997581 -0.22740739]
                         #t_13_1_4.data=[-0.30848055  0.12210877  0.28535858 -0.06227226  0.27997581 -0.22740739]
                         #t4_1.shape=([4]1+[2]4+[0]6+[-2]4+[-4]1,)
                         #t_13_1.shape=([4]1+[2]4+[0]6+[-2]4+[-4]1,)
                    
                    #t4_1 and t_13_1 have same qsp,  but their data differ 
                    
                    a drastic example of this can be found in vmps.mps.split_ksites, the order to split matters
                
        """
        
        qsp_class = self.qsp_class
        qn_class = qsp_class.QnClass 
       
        qsp = [] 
        leg_map = {}
        ii = 0  #ii points to legs of res 
        l_prev =  -1 
        leg_groups_to_merge = args
        for ll in args:
            l = ll[0]
            if l-l_prev>1: 
                qsp.extend(self.QSp[l_prev+1: l])
                for x in range(1, l-l_prev): 
                    #leg_map[l_prev+x] = (x+ii-1, )
                    leg_map[x+ii-1] =(l_prev + x, )
                ii +=  l-l_prev-1 
                
            l_prev = ll[-1]
            leg_map[ii] = range(l, l_prev + 1)
            temp = [self.QSp[l].copy() for l in leg_map[ii]]
            qsp.append(qsp_class.prod_many(temp))
            ii += 1 
        qsp.extend([self.QSp[i].copy() for i in range(args[-1][-1]+1, self.rank)] )
        for x in range(args[-1][-1]+1, self.rank):
            leg_map[ii] = (x, )
            ii += 1  
        
        #print_vars(vars(),  ['qsp'])
        
        
        t2 = iTensor(QSp=qsp, totQN=self.totQN.copy(), 
                dtype=self.dtype, use_gpu=self.use_gpu)
        t3 = self
        t3ind = list(range(t3.nidx))
        for i in range(t2.nidx): 
            qn_id_tuple_2 = t2.Addr_idx[:, i]
            qn_tuple_2 = [t2.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_2)]
            dim_tuple_2 = [t2.QSp[i2].Dims[q] for i2, q in enumerate(qn_id_tuple_2) ]
            data_block_2 = t2.get_block(i)
            data_block_2 = data_block_2.reshape(dim_tuple_2, order='F')
            #print_vars(vars(), ['qn_id_tuple_2', 'dim_tuple_2'])
            start_dicts= {x: {} for x in range(t2.rank) }
            assert data_block_2.base is t2.data   #in the following, data_block_2 should effect as a "pointer"
                
            #for j in range(t3.nidx):
            matched_j_list = []
            for j in t3ind: 
                qn_id_tuple_3 = t3.Addr_idx[:, j]
                qn_tuple_3 = [t3.QSp[_i].QNs[q] for _i, q in enumerate(qn_id_tuple_3)]
                dim_tuple_3 = [t3.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_3) ]
                
                match = True 
                for l, v in leg_map.items(): 
                    if qn_tuple_2[l] != qn_class.sum([qn_tuple_3[_i] for _i in leg_map[l]]): 
                        match = False
                        break 
                    
                if match:     
                    matched_j_list.append(j)
                    dim_tuple_3 = [t3.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_3) ]
                    dim_delta  = np.zeros(t2.rank, dtype=int)
                    
                    qn_id_tuple_grouped = [tuple([qn_id_tuple_3[y] for y in leg_map[x] ]) for x in range(t2.rank)]
                    dim_tuple_grouped = [np.prod([dim_tuple_3[y] for y in leg_map[x] ]) for x in range(t2.rank)]
                    
                    sh_start = np.zeros(t2.rank, dtype=int)
                    sh_end = np.zeros(t2.rank, dtype=int)
                    for k, v in enumerate(qn_id_tuple_grouped): 
                        temp = start_dicts[k]
                        if v in temp: 
                           sh_start[k] = temp[v]
                        else:
                            if 'totdim' in temp: 
                                sh_start[k] = temp['totdim']
                                temp[v] = temp['totdim']
                                temp['totdim'] += dim_tuple_grouped[k] 
                            else:
                                sh_start[k] = 0
                                temp[v] = 0
                                temp['totdim'] =  dim_tuple_grouped[k] 

                    
                    dim_delta = np.asarray(dim_tuple_grouped, dtype=int)
                    sh_end = sh_start + dim_delta 
                    #print_vars(vars(), ['i', 'j', 'qn_id_tuple_3', 'qn_id_tuple_grouped', 
                    #    'dim_tuple_grouped',  'sh_start', 'sh_end'], sep=', ')
                    db3 = t3.get_block(j)
                    db3 = db3.reshape(dim_delta,  order='F')
                    sl = [slice(sh_start[iii], sh_end[iii]) for iii in range(t2.rank)]
                    data_block_2[tuple(sl)] = db3 
            
            for j in matched_j_list: 
                t3ind.remove(j)
        return  t2         
    
    def merge_qsp_all(self): 
        return self.merge_qsp(list(range(self.rank)))
    
    def split_2to3(self, which, qsp_list:list): 
        """
            I dont know how to split in general, 
            so I first write a simple one. this can be used in mps 
            params:
                which: 
                    which leg to split,  take value in [0, 1]
        """
        #warnings.warn("todo: avoid duplicate compare")
        assert self.rank == 2
        if which == 0: 
            qsp = qsp_list+[self.QSp[1]]
            common2, common3 = 1, 2
        elif which == 1: 
            qsp = [self.QSp[0]] + qsp_list
            common2, common3 = 0, 0 
        #assert self.qsp_class.prod_many(qsp_list)==self.QSp[which]   #use when debug 
        t2 = self 
        t3 = iTensor(QSp=qsp, totQN=self.totQN.copy(), 
                dtype=self.dtype, use_gpu=self.use_gpu)
        #print_vars(vars(),  ['t3.device'])
        jj = list(range(t3.nidx))
        for i in range(t2.nidx): 
            qn_id_tuple_2 = t2.Addr_idx[:, i]
            matrix = t2.get_block(i)
            if which == 0: 
                n, m = (t2.QSp[i2].Dims[q] for i2, q in enumerate(qn_id_tuple_2))
                matrix = matrix.reshape((n, m), order='F')
            p = 0
            temp = []
            for j in jj: 
                qn_id_tuple_3 = t3.Addr_idx[:, j]
                if qn_id_tuple_2[common2] == qn_id_tuple_3[common3]:  #they match 
                    if which == 0:  #d2 = m
                        d0, d1, d2 = (t3.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_3) )                                           
                        data = matrix[p:p+d0*d1, :].ravel(order='F')
                        p += d0*d1  
                        t3.set_block(j, data)
                    else:   #d0 = n
                        size = t3.Block_idx[1, j]
                        data = matrix[p:p+size]
                        p += size
                        t3.set_block(j, data)
                    temp.append(j)
            for t in temp:
                jj.remove(t)
                        
        return  t3         


    def merge_3to2(self, which): 
        """
            params:
                which: can be (0, 1) or (1, 2)
        
        """
        q0, q1, q2 = self.QSp 
        if which == (0, 1) :
            qm = q0.tensor_prod(q1)
            qsp = [qm, q2]
            common2, common3 = 1, 2
        else:
            qm = q1.tensor_prod(q2)
            qsp = [q0, qm]
            common2, common3 = 0, 0
                
        t2 = iTensor(QSp=qsp, totQN=self.totQN.copy(), 
                dtype=self.dtype, use_gpu=self.use_gpu)
        t3 = self
        jj = list(range(t3.nidx))
        for i in range(t2.nidx): 
            qn_id_tuple_2 = t2.Addr_idx[:, i]
            n, m = (t2.QSp[i2].Dims[q] for i2, q in enumerate(qn_id_tuple_2) )
            matrix = t2.get_block(i)
            matrix = matrix.reshape((n, m), order='F')
            temp = []
            xx = 0; yy = 0
            for j in jj:
                qn_id_tuple_3 = t3.Addr_idx[:, j]
                if qn_id_tuple_2[common2] == qn_id_tuple_3[common3]:  #they match 
                    d0, d1, d2 = (t3.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple_3) )
                    data = t3.get_block(j)
                    if which == (0, 1): 
                        nx = d0*d1
                        matrix[xx:xx+nx, :] = data.reshape((nx, m), order='F')
                        xx += nx 
                    else: 
                        ny = d1*d2
                        matrix[:, yy:yy+ny] = data.reshape((n, ny), order='F')
                        yy += ny 
                    temp.append(j)
            for t in temp:
                jj.remove(t)
        return  t2         

    def permutation(self, P, buffer=None, use_buf=False):
        """
            permutation 是一个re-index的操作, 这一操作保持指标集整体不变，而局部置换
            todo: 
                I may write a inplace permutation, so needless QSp.copy()
                ind_labels  also need permuted
            for permutation of iTensor, first we permute the order of QNs, next we permute each block which is a nTensor
            注意这里是Fortran order
            permutatiion won't change total quantum number
                
        """
        rank = self.rank
        QSp=[self.QSp[P[i]] for i in range(rank)]   #no copy, faster
        totQN = self.totQN
        res=iTensor(rank, QSp, totQN, buffer=buffer, 
                dtype=self.dtype, use_buf=use_buf, use_gpu=self.use_gpu)
        pos=np.empty(self.rank, int)
        Dims=np.empty(self.rank, int)
        
        #permute each block
        #找到原来的block的位置与新的位置间的转换pidx<-->qidx
        #warnings.warn("using array_permutation_np")
        for n  in range(self.nidx):
            
            data = self.get_block(n)
            
            pos[0] = 0  # for rank=0
            pos[0:rank] = self.Addr_idx[0:rank,n]  # qn number index 
            pos_new = [pos[P[i]] for i in range(rank)]
            data_p = res.get_block(pos_new)
            shape = self.get_block_shape(n)
            if 0:            
                Dims= np.asarray(shape, dtype=int)
                #attention_may_be_not_efficient  可以改成inplace 
                #res.data_p[qidx:qidx+totDim]=array_permutation.array_permutation_np(data,rank,Dims,P)
                # 如果错误信息为 #error: failed in converting 5th argument `b' of array_permutation_64_ifort.array_permutation_fort_parallel to C/Fortran array
                #则检查 Dims，其中可能包含了0维 
                array_permutation.array_permutation_inplace(data, self.rank, Dims, P, data_p)
                
            else:
                #data = as_strided(data, shape)   # this seems not needed, because reshape should not allocate memory 
                data1 = data.reshape(shape, order='F')
                data2 = data1.transpose(P)
                #note: reshape and transpose won't allocate memory, which can be checked by shares_memory. So this code is enough efficient.
                #assert np.shares_memory(data, data1)
                #assert np.shares_memory(data1, data2)
                data_p[:] = data2.ravel(order='F')
                

        return res
    
    #def transpose 
    transpose = permutation   # compatible with numpy 
    
    def contract_core(self, T2, div, preserve_qsp=False, data=None, use_buf=False):
        """
            把T1，和T2的非零block 如果量子数组合相等则收缩
            todo:
                force check qsp pair are reverse of each other before contract
            locals:
                div: num. of legs to be contracted for each tensor
                buffer: use buffer to save data of T3
        """

        rank1, rank2 = self.rank, T2.rank 
        rank3 = rank1+rank2-div-div
        #tQN = self.totQN+T2.totQN
        tQN = self.totQN.__add__(T2.totQN)
        shift = rank1-div
        QSp = self.QSp[:shift]  #copy qsp is more robust, while no copy is faster
        QSp.extend(T2.QSp[div:rank2])

        if rank3==0:
            QSp = []  #QSp = [self.QSp[0].null()]
        
        
        dtype = np.result_type(self.dtype, T2.dtype)  # 自动根据 self.dtype 和 T2.dtype 推导最精准的输出类型（完美支持 32位/64位 和 实数/复数）
        
        T3 = iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, 
                dtype=dtype, use_buf=use_buf, use_gpu=self.use_gpu)
        T3.data[:]=0.0   #T3.data=np.zeros(T3.totDim)  #this is really bad, need to re-allocate space for data
        
        nidx3 = 0
        iQN1=np.empty(self.rank + 1, int)
        iQN2=np.empty(T2.rank + 1, int)        
        iQN3=np.empty(T3.rank + 1, int) # +1 to avoid T3.rank=0

        for idx2 in range(T2.nidx):
            iQN2[0] = 0  #!for rank=0
            iQN2[0:rank2] = T2.Addr_idx[0:rank2, idx2]
            if div == rank2: 
                Dim2 = 1 
            else:
                Dim2 = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div, rank2)], dtype=int)
            
            #issue: the following may be replaced by np.tensordot
            Dimc = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div)], dtype=int)  # np.prod([])=1.0, so use dtype=np.int 
            
            data2 = T2.get_block(idx2)
            data2 = data2.reshape((Dimc,Dim2), order='F')    
                
            for idx1 in range(self.nidx):
                iQN1[0] = 1 #!for rank=0
                iQN1[0:rank1]=self.Addr_idx[0:rank1,idx1]
                iseq = np.all(iQN1[shift:shift+div] == iQN2[0:div])  #注意这里写得不适当，准确地，如果是U1 symm 的话应该是T1, T2相应的量子数的值正好差个符号, 而这里是用量子数的位置处理了, 并假定....写不清楚啊
                if not iseq:  #如果量子数组合相等则收缩
                    continue
                if shift == 0: 
                    Dim1 = 1   #when shift = 1.0,  np.prod yields 1.0, should be converted to int 
                else:
                    Dim1 = np.prod([self.QSp[i].Dims[iQN1[i]] for i in range(shift)])
                
                data1 = self.get_block(idx1)
                data1 = data1.reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                
                iQN3[0:shift] = iQN1[0:shift]  #iQN3[0] = 1 #!for rank=1
                iQN3[shift:rank3] = iQN2[div:rank2]
                
                data3 = T3.get_block(iQN3[:T3.rank]).reshape((Dim1,Dim2), order='F')    
                
                #use_gpu = self.use_gpu
                #transfer_data = False
                #if self.USE_GPU_FOR_BLOCK and data3.size >= self.USE_GPU_MUL_LIM:
                #    use_gpu = True
                #    transfer_data = True
                #common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                #        dtype = dtype, use_gpu=use_gpu, transfer_data=transfer_data)
                
                use_gpu = self.use_gpu
                if self.use_gpu == 2 and (Dim1*Dim2*Dimc) < self.USE_GPU_MUL_LIM:
                    use_gpu = 0
                common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                        dtype = dtype, use_gpu=use_gpu)
        
        return T3

    
    def prepare_leg(self,T2, V1, V2, info=0):
        """
            todo: rewrite this function using C language 
            issue: this seems not efficient, improvement is needed.
                However,  for use of low rank tensors in MPS algs, this is not a big issue.
                Furthur, after putting it into tensor_player, this is not at all a problem. 
                So dont need to fix this soon.
                
        """
        if 0:
            V1 = np.array(V1) if not isinstance(V1, np.ndarray) else V1
            V2 = np.array(V2) if not isinstance(V2, np.ndarray) else V2
            #V1 = np.asarray(V1, dtype=object) if not isinstance(V1, np.ndarray) else V1
            #V2 = np.asarray(V2, dtype=object) if not isinstance(V2, np.ndarray) else V2
            print_vars(vars(),  ['V1', 'V2'])
            if self.rank<V1.size: V1 = V1[:self.rank]
            if T2.rank<V2.size: V2 = V2[:T2.rank]
            V_1n2 = np.intersect1d(V1, V2, True)
            V3 = np.setxor1d(V1, V2, True,)
            
        else:
            V1 = tuple(V1)
            V2 = tuple(V2)
            if self.rank<len(V1): V1 = V1[:self.rank]
            if T2.rank<len(V2): V2 = V2[:T2.rank]
            V_1n2 = set(V1).intersection(V2)
            #V3 = np.setxor1d(V1, V2, True,)
            #V3 = [i for i in V1 + V2 if i not in V_1n2]
            #V3 = np.ndarray(self.rank + T2.rank-len(V_1n2), dtype=object)
            V3 = []

        
        Vp1 = np.empty(self.rank, np.int32)   # order of permutated index for tensor1 
        Vp2 = np.empty(T2.rank, np.int32)
        
        k = 0  #indexing Vp1
        j = 0  #indexing Vp2
        
        for i in range(self.rank): #below calculate Vp1, Vp2
            if V1[i] not in V_1n2:
                Vp1[k] = i   # 这里之所以用k，而非直接用Vp1[k]是因为下面k的值要继续，对于l 也一样 p意指position，记录的T1中没有收缩的指标(leg)的编号
                V3.append(V1[i])   #T3 来自T1 V1的外腿
                k += 1 

        #for i in range(len(V_1n2)):  #这一段程序做了两件事：1.验证内线上维数相等，2.计算了totDim3
        for v12 in V_1n2:
            pos1 = V1.index(v12)
            pos2 = V2.index(v12)
                
            Vp1[k] = pos1
            k  += 1 #注意这里k接着上面的值了
            Vp2[j] = pos2
            j += 1   
            
            if self.Dims[pos1] != T2.Dims[pos2]:  #note this checking is still not complete. One should check qsp1 == qsp2.conj() instead  
            #if self.QSp[pos1] != T2.QSp[pos2].conj():
                #print_vars(vars(),  ['pos1', 'pos2'])
                msg ="""error, dim of index to be contracted not equal: 
                    {0.type_name}, {1.type_name}
                    ind_label_1={V1}, ind_label_2={V2}, ind_label_1n2={V_1n2}
                    dims_1={0.Dims}, dims_2={1.Dims}
                    
                    """.format(self, T2, V1=V1, V2=V2, V_1n2=V_1n2)
                #self.save('/tmp/aaa'); T2.save('/tmp/bbb')
                raise Exception(msg)
        
        for i in range(T2.rank):
            if V2[i] not in V_1n2:
                Vp2[j] = i
                V3.append(V2[i])  #T3 V3 来自T2 V2的外腿
                j += 1 
        return V_1n2, Vp1, Vp2, V3 

    def contract(self, T2, V1=None, V2=None, final_ind_labels=None, 
            out_Vc=False, data=None, use_buf=False, 
            return_v3 = False, 
            preserve_qsp=False, track_name=0,  info=0):
        """
            issue: todo: in future not return_v3 by default 
            V1,2,3 are arrays (maps) 张量指标 —> 自然数。用自然数来标记所有张量的指标
            在V1,V2中可能有相同的元素，存在s_1n2中，
            T1, T2 contract yielding T3
            params:
                Vi: of type np.ndarray(,"int")
                Vp1,先记录了T1的外腿，后记录内腿指标； Vp2先记录了内腿，后记录了外腿指标
                use_buf: use data buffer for T3
        """
        
        V1 = self.ind_labels if V1 is None else V1   #issue: when self contract with self but with different ind_labels, this cause problem
        V2 = T2.ind_labels if V2 is None else V2
        
        V_1n2, Vp1, Vp2, V3 = self.prepare_leg(T2, V1, V2)        
        T1=self
        try:
            nT1 = T1.permutation(Vp1, use_buf=use_buf)    #把T1 按照 Vp1 重排
            nT2 = T2.permutation(Vp2, use_buf=use_buf)
            T3 = nT1.contract_core(nT2, len(V_1n2), data=data, use_buf=use_buf)
        except Exception as err:
            V1 = tuple(V1)
            V2 = tuple(V2)
            msg = '\n\t'.join([
                'additional err info:  ', 
                'If the exception is like: IndexError:index 4 is out of bounds for axis 1 with size 4,it is', 
                'most likely the qsp of legs to be contracted not match,  espcially lack a reverse of qsp', 
                'names: {0.type_name}\t{1.type_name}'.format(self, T2), 
                'ind labels: {}, \t{}'.format(V1, V2) , 
                'ind to contract:%s'%(tuple(V_1n2), )] 
                #+ [ '\t'+str((i, self.shape[V1.index(i)], T2.shape[V2.index(i)])) for i in V_1n2] 
                )
            if 1:
                for i in V_1n2:
                    sh1, sh2 = self.shape[V1.index(i)],  T2.shape[V2.index(i)]
                    if sh1.conj() != sh2:
                        msg  += '\t\n\t -***-> : {}, {}, {}'.format(i, sh1, sh2)
            
           
            if not err.args: 
                       err.args=('',)
            err.args = (str(err.args[0]) + "\n"*2 + msg,)+err.args[1:]
            raise 
        
        if info>0: 
            T3.type_name = str(self.type_name) + "-" + str(T2.type_name)
        T3.ind_labels= V3
        if final_ind_labels is not None:
            ll = list(V3)  # head.ind_labels may be np.ndarray, convert it to list
            order = [ll.index(i) for i in final_ind_labels]
            T3 = T3.permutation(order)
            T3.ind_labels = final_ind_labels 
        
        if out_Vc:
            Vc= V_1n2
            return  T3, V3, Vc 
        else:
            if return_v3:
                return T3, V3
            else:
                return T3
    
    @staticmethod 
    def contract_tensor_list(tlist, final_ind_labels=None, check=False):   #def ctl
        """
            requires each t.ind_labels is not None 
            params:
                tlist can be a nested list, still it is efficient! like [t1, t2, [t3, t4]]
        """
        if check:
            assert len(temp)==len(set(temp))  #simple check,  no duplicate tensors, so that I can use ind_labels to contract        
        head = tlist[0]
        if isinstance(head, list):
            head = iTensor.contract_tensor_list(head)
        try:
            for t in tlist[1: ]:
                if isinstance(t, list):
                    t = iTensor.contract_tensor_list(t)
                head = head.contract(t)
        except Exception:
            #print_vars(vars(),  ['head.ind_labels', 't.ind_labels'])
            raise  
        if final_ind_labels is not None:
            ll = list(head.ind_labels)  # head.ind_labels may be np.ndarray, convert it to list
            order = [ll.index(i) for i in final_ind_labels]
            head = head.permutation(order)
            head.ind_labels = final_ind_labels 
        
        #print_vars(vars(),  ['head.ind_labels'])
        return head 
    ctl = CTL= contract_tensor_list 
    
    def set_ind_labels(*args):
        n = len(args)//2 
        for i in range(n):
            o = args[i*2]
            ind = args[i*2+1]
            o.ind_labels= ind 
    
    def dot(self, other, data=None,  use_buf=False): 
        """
            mainly for compatible with numpy 
            multipyly tow rank-2 tensors 
        """
        assert self.rank == 2 and other.rank == 2    
        return self.contract_core(other, div=1, data=data, use_buf=use_buf)
    
    def norm(self): 
        #a = list(range(self.rank))
        #b = list(range(self.rank))
        #temp, _= self.contract(self, a, b)
        #return math.sqrt(temp.data[0])
        return np.linalg.norm(self.data)
    
    def trace(self):
        """
            status_1_uncheck
            see Tensor_Trace in f90
            从此函数看出， 要确定iTensor的对角线首先要确定它的对角块，即出入脚量子数相等;
            然后把对角块展开成2D matrix 求其trace
            
            see also Tensor_svd.trace_rank2, which only suits for rank2 tensor 
        """
        X = 0.0
        rank = self.rank
        rank2 = rank//2
        iQN=np.empty(self.rank,"int")
        for idx  in range(self.nidx):
            iQN[0:rank]=self.Addr_idx[0:rank,idx]
            #IsDiag = iArrayEq(rank2, iQN[0], iQN[rank2+1])
            IsDiag=np.all(iQN[:rank2]==iQN[rank2:])
            #judge wether on diagnal line
            if  not  IsDiag:  continue
            d = 1
            for i  in range(rank2):
                d = d*self.QSp[i].Dims[iQN[i]]
                a=self.Block_idx[0,idx]
            X+= self.data[a:a+d*d].reshape(d,d).trace()
        return X
    
    def trace_ind_pairs(self, *ind_pairs):
        """
            this is exactly the general trace function
            I named it trace_ind_pairs just because trace is occupied by meth above this
            rename it just to trace in future. 
            params:
                ind_pairs: 
                    pairs like (0, 1), (3, 5), ('a', 'c')
                    it is actually ind label pairs
        """
        if self.ind_labels is None:
            self.ind_labels= list(range(self.rank))
        self.ind_labels= tuple(self.ind_labels)  #convert to tuple, if it is a ndarray
        id_list = []
        for i, j in ind_pairs:
            #id=iTensor.identity([self.shape[j], self.shape[i]])
            i1, j1 = self.ind_labels.index(i), self.ind_labels.index(j)
            id=iTensor.identity([self.shape[j1], self.shape[i1]])
            id.ind_labels= (i, j)
            id_list.append(id)
        id_list.insert(0, self)
        res=iTensor.contract_tensor_list(id_list)
        
        return res
            
    def partial_trace(self, site_list):
        """
            self is assumed to be a density matrix 
            
            this is a temporary implementation during calc central charge
            a full and fast implementation should not use tensor contraction, but direct sum uing index of sparse tensor
            param: 
                site_list  ranges in [0, self.rank-1]
        """
        rank = self.rank
        assert rank%2 == 0 
        rank_half= rank//2
        
        nsite = len(site_list)

        qsp = self.QSp[0].copy_many(2*nsite, reverse=list(range(nsite, 2*nsite)))
        #I = iTensor(2, qsp.copy_many(2, reverse=[1]), qsp.QnClass.qn_id())
        I = iTensor.unit_tensor(2*nsite, qsp)

        v1 = range(rank)
        #v2 = [i, rank_half + i]
        #v2 = xrange(nsite) + [rank_half  + i for i in nsite]
        v2 = site_list + [rank_half + s for s in site_list]
        res, leg= self.contract(I, v1, v2)
        warnings.warn('need reorder legs? see iTensor')
        
        #print 'leg', leg
        
        return res
    
    def conjugate_new(self, d, buffer=None, use_buf=False):
        """

            first, permute leg order from {0, 1, ..d-1, d, ..rank-1} to {d, ..rank-1, 0, 1, ..d-1 }即把朝上(下)的leg 扳到下(上)面
            second, reverse qn
        """
        ord = np.ndarray(self.rank, "int")
        for i in range(d,self.rank):
            ord[i-d] = i
        for i  in range(d):
            ord[i+self.rank-d] = i
        temp=self.permutation(ord, buffer=buffer, use_buf=use_buf)
        
        Tc = temp.copy_struct_new()
        Tc.reverse_qsp()
        Tc.set_data_entrance()
        Tc.data = temp.data
        #temp.reverse_qsp(inplace=True)
        return Tc

    def conjugate(self, d, buffer=None, use_buf=False):
        """
            see Tensor_Conjugate in f90
            first, permute leg order from {0, 1, ..d-1, d, ..rank-1} to {d, ..rank-1, 0, 1, ..d-1 }即把朝上(下)的leg 扳到下(上)面
            second, reverse qn
        """
        ord = np.ndarray(self.rank, "int")
        for i in range(d,self.rank):
            ord[i-d] = i
        for i  in range(d):
            ord[i+self.rank-d] = i
        temp=self.permutation(ord, buffer=buffer, use_buf=use_buf)
        
        temp.reverse_qsp(inplace=False)
        if self.dtype == complex:  # this line maybe not needed, numpy may check this before do the conj 
            np.conj(self.data, self.data)
        return temp
    
    def conj(self): 
        """
            reverse the qsp and conjugate the data (if it is complex)
        """
        #assert self.rank == 3  
        A = self.copy()
        A.reverse_qsp()
        np.conj(A.data, out=A.data)  #A_conj.data = A.data.conj()
        return A 
    
    def dag(self):
        """
            only for self.rank = 2
            hermite conjugate of self
        
        """
        return self.T.conj()
    
    def inv(self):   #def inverse 
        """
            the tensor is required to be square matrix
        """
        assert self.rank == 2 
        res = self.copy()
        for i in range(self.nidx):   # 遍历非零blocks
            data = self.get_block(i, linear=False, order='F')
            data_inv = scipy.linalg.inv(data)
            data_inv = data_inv.ravel(order='F')
            res.set_block(i, data_inv)
        return res 
    
    
    def diagonal_tensor_rank2(qsp): 
        """
            A tensor t satisfis:  tt* = t*t  =  I 
            Note a diagonal tensor is NOT identity tensor.
            Not any qsp can make diagonal_tensor_rank2
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

    def conj_new(self, i, use_buf=False): 
        """
            the correct way to reverse the direction of legs of iTensor 
            #improve: it is better to only provide a u tensor, that is enough 
            
            note: 
                it seems that not every itensor can do conj_new, that is
                because not able to construct for rank2-dual-tensor for any
                qsp. see note of diagonal_tensor_rank2. 
                
                this also depends on whether change totqn after conj_new
                
                so this func is not very complete
            issue:
                这里的几个函数需要整理，名字太乱，功能重复
        """
        #assert self.rank == 3  

        if self.symmetry == 'U1' : 
            q = self.QSp[i].copy()
            q.reverse()
            u = iTensor.diagonal_tensor_rank2(q)
            
            label = list(range(self.rank))
            label[i] = -1
            label_u = [-1, 1000]
            #A, _ = self.contract(u, [0, 1, 2], [2, 3],  use_buf=use_buf)
            A = self.contract(u, label, label_u, use_buf=use_buf)
            if i != self.rank-1:  #need re-order the legs 
               temp = list(range(i)) + [A.rank-1] + list(range(i, self.rank-1))
               A  = A.transpose(temp) 
        elif self.symmetry in ['Z2', 'Travial'] : 
            A = self.copy(use_buf=use_buf)
        else: 
            raise  ValueError('symmetry is %s'% self.symmetry)
        #np.conj(A.data, out=A.data)   # note should no complex conjugate here !
        return A 

    def reverse_qsp(self, inplace=True):
        """
        see iTensor_ReverseQN in f90
        what's the relation of conjugate and reverse_qsp?
        """
        if not inplace:
            QSp = [q.copy() for q in self.QSp]
            totQN = self.totQN.copy()
            self.QSp = QSp
            self.totQN = totQN
        for n in range(self.rank):
            self.QSp[n].reverse()
        self.totQN.reverse()

    def dual(self, use_buf=False):
        rank = self.rank
        qsp = [q.copy() for q in self.QSp]
        for q in qsp:
            q.reverse()
        totqn = self.totQN.copy()
        totqn.reverse()
        res= iTensor(rank, qsp, totqn, use_buf=use_buf)
        res.data[:] = self.data[:]
        return res

    def to_nTensor(self):
        """
        status_1_verified
        """
        #init_nTensor(self.rank, T.Dims, Tp)
        Tp=nTensor(self.rank,self.Dims)
        
        rank = self.rank
        Dims= np.empty(rank, 'int')
        
        #for i  in xrange(rank): Dims[i] = self.QSp[i].nQN
        #quan_num_space_size=[self.QSp[i]] for i in xrange(rank)]
        
        
        #block_size=np.empty(rank,"int")
        for idx  in range(self.nidx):
            #这一loop先确定block坐标, block的大小
            block_pos_in_data= self.Block_idx[0,idx]
            qn_comb_pos = self.Block_idx[2,idx]            
            #qn_comb_pos is the position of QN, 线性坐标与多维坐标点转换
            
            #pos is coordinate of a set of quantum number 
            qn_coord=self.get_position_rev(qn_comb_pos)
            block_size=[self.QSp[i].Dims[qn_coord[i]] for i in range(rank)]
                #block 每个指标的维数,与上面的Dims意思不同
            
            #attention_omitted_something
            Ti=nTensor(rank, block_size, shallow=False)
            print("block_pos",idx,"block_size",block_size,"idx",idx,"qn_coord",qn_coord)

            for dat_pos_in_blk  in range(self.Block_idx[1,idx]):
                print("\t dat_pos_in_blk ", dat_pos_in_blk)
                #这一loop确定block内的坐标
                data_coord_in_block=Ti.get_position_rev(dat_pos_in_blk)
                #data_coord_in_block is coordinate in the block <===  p is linear coordinate
                
                data_coord=[]
                for i  in range(rank):
                    d = 0
                    for j  in range(qn_coord[i]):
                        d = d+self.QSp[i].Dims[j]
                    
                    data_coord.append(int(d + data_coord_in_block[i]))
                
                data_pos=Tp.get_position(data_coord)
                print("\t\t data_coord",data_coord,"data_pos",data_pos, end=' ')
                
                Tp.data[data_pos] = self.data[block_pos_in_data+dat_pos_in_blk]
                print(round(self.data[block_pos_in_data+dat_pos_in_blk],10))
        return Tp
    
    def to_ndarray(self, qn_val_tuple=None,  data_order='F'): 
        """
            
            converge a symm iTensor instance to a 
            dense np.ndarray instance 
            this func should be invaluable in dev and debug 
            not for production use 
            param:
                qn_val_tuple: only select a block with given qns
            note:
                运行以下code
                    from tensor import iTensorFactory 
                    pau = pauli_mat()
                    #sx = pau['sx']
                    sx = np.asarray([[0, 1], [1, 0]])
                    sxx=np.multiply.outer(sx, sx).transpose([0, 2, 1, 3]).reshape(4, 4)
                    print sxx 
                    sx = iTensorFactory.pauli_mat('Z2')['sigma_x']
                    sxx=sx.direct_product(sx) 
                    #先reshape，后to_ndarray
                    print sxx.merge_qsp((0, 1), (2, 3)).to_ndarray(data_order='F').round(5)
                    #先to_ndarray，后reshape 
                    print sxx.to_ndarray(data_order='F').reshape(4, 4)
               输出为 
                    [[0 0 0 1]
                     [0 0 1 0]
                     [0 1 0 0]
                     [1 0 0 0]]
                    [[ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]
                     [ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]]
                    [[ 0.  0.  0.  1.]
                     [ 0.  0.  1.  0.]
                     [ 0.  1.  0.  0.]
                     [ 1.  0.  0.  0.]]
               之所以不同，是因为，对称张量的reshape 和 ndarray的reshape 名字相同，但
               实质有差别， iTensor.reshape 包含了将量子数指标的合并 
            
        """
        if qn_val_tuple is not None:
            #todo: implement this. I use this mainly for doing ED,  selecting specific sector
            raise NotImplemented 
        #res=nTensor(self.rank,self.Dims)
        #res = np.ndarray(self.Dims, dtype=self.dtype, order=order)
        res = np.ndarray(self.Dims, dtype=self.dtype) 
        res[: ] = 0.0 
        blocks= {}
       
        qn_range = [q.nQN for q in self.QSp]
        #qn_range.reverse()
        qn_id_tuple_all = itertools.product(*tuple([list(range(i)) for  i in qn_range]) )
        qn_id_tuple_all = list(qn_id_tuple_all)
        #print_vars(vars(), ['self.Dims', 'list(qn_id_tuple_all)', 'self.Addr_idx[:self.rank, :]']) 
        data_blocks= {}
        for i in range(self.nidx): 
            #qn_id_tuple_symm.append(self.Addr_idx[:self.rank, i])  
            temp=self.Addr_idx[:self.rank, i].tolist()  
            data_blocks[tuple(temp)] = self.get_block(i)
            
        start_dicts= {x: {} for x in range(self.rank) }
        for qn_id_tuple in qn_id_tuple_all: 
            block_origin = np.zeros(self.rank, dtype=int)
            dim_tuple = [self.QSp[i1].Dims[q] for i1, q in enumerate(qn_id_tuple) ]
            for k, v in enumerate(qn_id_tuple): 
                temp = start_dicts[k]
                if v in temp: 
                   block_origin[k] = temp[v]
                else:
                    if 'totdim' in temp: 
                        block_origin[k] = temp['totdim']
                        temp[v] = temp['totdim']
                        temp['totdim'] += dim_tuple[k] 
                    else:
                        block_origin[k] = 0
                        temp[v] = 0
                        temp['totdim'] =  dim_tuple[k] 
            block_end = block_origin + np.asarray(dim_tuple)
            if qn_id_tuple in data_blocks: 
                #sl = [slice(block_origin[iii], block_end[iii]) for iii in range(self.rank)]
                sl = [slice(block_origin[iii], block_end[iii]) for iii in range(self.rank)]
                res[tuple(sl)] = data_blocks[qn_id_tuple].reshape(dim_tuple, order=data_order)
            
            #print_vars(vars(), ['block_origin', 'block_end'], sep=' ')
                
            
        return res 
   
    def to_matrix(self):
        """
            only for rank-2 tensor, and taking use of to_ndarray
        """
        assert self.rank == 2 
        raise NotImplemented 
   
    def from_ndarray(self, array, qsp): 
        raise NotImplemented('todo: to be implemented')
        pass 
    
    def get_nposition(self):
        """
        see iTensor_GetNPosition
        status_1
        """
        pass
    
    def direct_product(self, T2, order='F', use_buf=False):
        """
            this function was originally defined under Tensor class, now is moved here
            parameters:
                T1^{I1}_{J1}, T2^{I2}_{J2}
            returns:
                T3^{I1I2}_{J1J2}
                i.e. index(QSp) of T3 is ordered J1, J2, I1, I2 in sequel


        """
        T1=self
        rank1 = T1.rank//2
        rank2 = T2.rank//2      
        rank3 = rank1+rank2

        QSp=[0  for i in range(T1.rank + T2.rank)] 
        for i in range(rank1):
            QSp[i] = T1.QSp[i].copy()
            QSp[i+rank3] = T1.QSp[i+rank1].copy()
        for i in range(rank2):
            QSp[i+rank1] = T2.QSp[i].copy()
            QSp[i+rank1+rank3] = T2.QSp[i+rank2].copy()
        
        #这里暗含了一个张量积标识fusion的过程
        #totQN = T1.totQN+T2.totQN
        totQN = T1.totQN.__add__(T2.totQN)
        T3= iTensor(T1.rank+T2.rank, QSp, totQN, use_buf=use_buf)
        #print 'ttt', totQN.val, T1.totQN.val, T2.totQN.val
        
        T3.data[:] = 0.0
        V1 = np.empty(self.MaxRank, "int")
        V2 = np.empty(self.MaxRank, "int")
        V3 = np.empty(self.MaxRank, "int")
        for idx1 in range(T1.nidx):
            pidx1 = T1.Block_idx[0, idx1]
            V1[0:T1.rank] = T1.Addr_idx[0:T1.rank,idx1]
            nA=1; mA=1
            for i in range(rank1):

                nA = nA*T1.QSp[i].Dims[V1[i]]
                mA = mA*T1.QSp[i+rank1].Dims[V1[i+rank1]]
                V3[i] = V1[i]
                V3[i+rank3] = V1[i+rank1]

            for idx2 in range(T2.nidx):
                pidx2 = T2.Block_idx[0, idx2]
                V2[0:T2.rank] = T2.Addr_idx[0:T2.rank, idx2]
                nB=1; mB=1
                for i in range(rank2):
                    nB = nB*T2.QSp[i].Dims[V2[i]]
                    mB = mB*T2.QSp[i+rank2].Dims[V2[i+rank2]]
                    V3[i+rank1] = V2[i]
                    V3[i+rank1+rank3] = V2[i+rank2]
                #idx3 linear position of  qn combination V3
                idx3 = T3.get_position(V3[:rank3*2])   
                
                pidx3 = T3.Block_idx[0, T3.idx[idx3]]
                    #print 'selffff.idx', T1.idx, T2.idx
                #Matrix_DirectProduct[T1.data[pidx1], nA, mA, T2.data[pidx2], nB, mB, T3.data[pidx3]]
                data1 = T1.data[pidx1:pidx1 + nA*mA].reshape((nA, mA), order=order)
                data2 = T2.data[pidx2:pidx2 + nB*mB].reshape((nB, mB), order=order)                
                if 0:
                    T3.data[pidx3:pidx3 + nA*nB*mA*mB] = common_util.matrix_direct_product(data1, data2).ravel(order="F")
                else:
                    T3.data[pidx3:pidx3 + nA*nB*mA*mB] = np.kron(data2, data1).ravel(order="F")

        return T3
    
    tensor_prod = direct_product  # def tensor_prod


    def direct_sum(self, other):
        """
            intro:
                MPS.direct_sum of two rank-3 tensors
                MPO.direct_sum of two rank-4 tensors
                but still lacks a general one
               
        """
        raise NotImplemented
        pass


    def matrix_view(self, n=None, order='C', data_order='F', round=None):
        """
            map rank (n, rank-n) tensor to a 2-d mat. 
            this method is used for debug or getting a intuitive view of the tensor

            Returns: 2-d array res
        """
        if n is None:
            n = self.rank//2
        pos= 0
        rank = self.rank
        dim1 = np.prod(self.Dims[:n])
        dim2 = np.prod(self.Dims[n:self.rank])
        res= np.zeros(np.prod(self.Dims[:self.rank]), dtype=self.dtype).reshape((dim1, dim2))
        
        ddd = [self.QSp[i].nQN for i in range(self.rank)]
        pos1 = 0 
        pos2 = 0        
        A = np.prod( [self.QSp[i].nQN for i in range(n)])
        B = np.prod( [self.QSp[i].nQN for i in range(n, self.rank)])
        iQN = np.zeros(self.rank, dtype="int")
        daaa = np.zeros(B, "int")
        for x in range(A):
            pos2 = 0
            for y in range(B):
                pos1 = daaa[y]
                p = x*A  + y
                #print 'p', p, x, y
                idx = self.idx[p]
                if idx < self.nidx:
                    d = np.prod([self.QSp[i].Dims[iQN[i]]  for i in range(rank) ])
                    d1 = np.prod([self.QSp[i].Dims[iQN[i]] for i in range(n)])
                    d2 = np.prod([self.QSp[i].Dims[iQN[i]] for i in range(n,self.rank)])            
                    start = self.Block_idx[0, idx]
                    res[pos1:pos1+d1, pos2:pos2+d2] = self.data[start:start + d].reshape(d1, d2, order=data_order)
                    #print res
                db = np.prod([self.QSp[i].Dims[iQN[i]] for i in range(n, self.rank)])
                            
                pos2 += db 
                da = np.prod([self.QSp[i].Dims[iQN[i]] for i in range(n)])
                daaa[y] += da 

                if order == "F": 
                    inc = 1
                    i = 0
                    while inc==1 and i<rank:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<self.QSp[i].nQN :
                            inc = 0
                        else:
                            iQN[i] = 0
                            i = i+1
                elif order == "C":
                    inc = True 
                    i = rank-1
                    while inc==True and i>= 0:
                        iQN[i] = iQN[i]+1
                        if iQN[i]<self.QSp[i].nQN :
                            inc = False
                        else:
                            iQN[i] = 0
                            i = i-1                
                
                
                #da = np.prod([self.QSp[i].Dims[kkk[i]] for i in xrange(n)])
        if round is None:
            return res
        else:
            return res.round(round)  
    matrix_form = matrix_view # def matrix_form 
    
    def matrix_view_rank2_z2(self):
        """
            a simple one,  only for test
        """
        shape = self.Dims 
        res = np.zeros(shape)
        
        bsh00 = [self.QSp[i].Dims[0] for i in range(2)]
        bsh11 = [self.QSp[i].Dims[1] for i in range(2)]
        data = self.get_block(0)
        res[: bsh00[0], : bsh00[1]] = data.reshape(bsh00, order='F') 
        data = self.get_block(1)
        res[bsh11[0]: , bsh11[1]:] = data.reshape(bsh11, order='F') 
        #print_vars(vars(), ['res.shape', 'bsh00', 'bsh11', 'res'])
        return res 
    
    def is_adjoint_to(self, other, precision=None, out_more=True):
        data1 = self.data
        r = other.rank//2
        other_conj = other.conjugate(r)
        data2 = other_conj.data
        dist= np.sum(np.abs(data1-data2))
        if precision is None:
            precision = 1e-14
        
        if dist < precision:
            res= True
        else:
            res= False
        
        if out_more:
            return res, "precision", precision,  "dist", dist 
        else:
            return res
    
    def is_hermite(self, precision=None, out_more=True):
        res= self.is_adjoint_to(self, precision, out_more)
        return res
    
    def is_unitary(self, tol=1e-14, show_more=True):
        assert self.rank == 2 
        #t = self.dot(self.conj())
        t = self.contract(self.conj(), (0, 1), (2, 1))
        t = t.to_ndarray()   # this is not effitient 
        d = t.shape[0]
        I = np.identity(d, dtype=t.dtype)
        #print_vars(vars(),  ['t.data-I'])
        diff = np.max(np.abs(t.data-I)) 
        #res  = np.allclose(t, I, atol=tol)
        res = True if diff<tol else False
        if show_more:
            res = (res, 'max diff:',  diff)
        return res 
    
    def is_close_to(self, other, tol=1e-15): 
        if other == 0:  #0, 1 for convenience 
            res = np.allclose(self.data, 0, atol=tol)
        elif other == 1: 
            other = iTensor.unit_tensor(self.rank, self.QSp)
            res= np.allclose(self.data, other.data, atol=tol)
        else: 
            res= np.allclose(self.data, other.data, atol=tol)
        return res 
    
    def commutator(self, other):
        assert self.rank == 2 
        return self.dot(other)-other.dot(self)
    
    def commutator_plus(self, other):
        """
            anti-commutator for fermions
        """
        assert self.rank == 2 
        return self.dot(other)+other.dot(self)
    
    if 0:
    #not completely correct
        def is_adjoint_to_2(self, other, precision=None, out_more=True):
            t1 = self.matrix_view()
            t2 = other.matrix_view()
            dist= np.sum(np.abs(t1-t2.T))
            if precision is None:
                precision = 1e-12
            
            if dist < precision:
                res= True
            else:
                res= False
            
            if out_more:
                return res, "precision", precision,  "dist", dist 
            else:
                return res

        def is_hermite_2(self, precision=None, out_more=True):
            res= self.is_adjoint_to_2(self, precision, out_more)
            return res

    def expand(self, QSp):
        """
            这个函数目前只支持 self.QSp 量子数数目  和 QSp量子数数目相等的情况
            否则会报错
        """
        for i in range(self.rank): 
            assert self.QSp[i] <=  QSp[i], (self.QSp[i], QSp[i])
        rank = self.rank
        T2 = iTensor(rank, QSp, self.totQN.copy())
        T2.data[:] = 0.0
        iQN = np.ndarray(rank, int)
        #Dims1 = np.ndarray(rank, int)
        #Dims2 = np.ndarray(rank, int)
        #transvers all dense blocks
        for idx  in range(self.nidx):
            p1 = self.Block_idx[0,idx] #-1
            iQN[0:rank] = self.Addr_idx[0:rank, idx]
            pidx2 = T2.get_position(iQN)
            idx2 = T2.idx[pidx2]
            p2 = T2.Block_idx[0,idx2] #-1
            # above find the start address p1 and p2 for corresonding blocks
            Dims1 = tuple((self.QSp[i].Dims[iQN[i]] for i in range(self.rank)))
            Dims2 = tuple((T2.QSp[i].Dims[iQN[i]] for i in range(T2.rank)))
            
            #within each dense block, transvers elements
            for ip1 in range(self.Block_idx[1,idx]):
                #ip1 is linear index, ip1->pos
                
                pos = np.unravel_index(ip1, Dims1, order='F')
                ip2 = np.ravel_multi_index(pos, Dims2, order='F')
                    
                T2.data[p2+ip2] = self.data[p1+ip1]
        return T2

    def act_on_qsp(self, qsp): 
        pass
    
    def invert_diagonal(self, use_pinv=False, pinv_tol=None): 
        """
            only for special use is idmrg 
            status: not tested 
        """
        assert self.rank == 2, "only supprted rank-2 square tensor"  
        for i in range(self.nidx): 
            qn_id_tuple = self.Addr_idx[:, i]
            p = self.Block_idx[0, i]
            size = self.Block_idx[1, i]
            data = self.data[p: p + size]
            d = self.QSp[0].Dims[qn_id_tuple[0]]
            data = data.reshape(d, d)
            
            if not use_pinv: 
                np.fill_diagonal(data, 1./data.diagonal())
            else: 
                if pinv_tol is not None : 
                    __pinv_tol = pinv_tol 
                else: 
                    __pinv_tol = float_info.epsilon * np.max(data.diagonal()) * d  
                data = np.linalg.pinv(data, __pinv_tol)
            self.data[p: p + size] = data.ravel(order='F')
    
    @staticmethod 
    def from_diagonal(qn_list, data_list): 
        pass 
    
    def set_data_to_zero(self): 
        self.data[: ] = 0.0
    
    def set_data_random(self):
        """
        """
        if self.use_gpu != 1:  
            self.data[:] = np.random.random(self.data.size)-0.5
        else:
            self.data[:] = cp.random.random(self.data.size)-0.5
        
    def reduce_1d_qsp(self, i):   # def reduce_dummy_index
        """
            suppose t has a qsp[i] whose totDim is 1. reduce it, such that 
            totqn of t is qsp[i].qn
            performance issue: 
                use contract is slow, direct operate on Addr_idx is better.
                only use this func for none production use 
        """
        if i == -1:
            i = self.rank-1
        qsp = self.QSp[i].conj()
        assert qsp.totDim == 1 
        qn = qsp.QNs[0].copy()
        #qn.reverse()
        leg = iTensor(QSp=[qsp], totQN=qn)
        leg.data[0] = 1.0 
        res = self.contract(leg, range(self.rank), [i])
        return res 
        
        remove_travial_ind = reduce_1d_qsp 
    
    def insert_1d_qsp(self, i, qn=None):  # def insert_dummy_index
        """
            performance issue: 
                see that under reduce_1d_qsp 
            params:
                i: insert the qsp before index i. 
                    if i>= self.rank, then the 1d qsp is appended as the last leg.  
        """
        if i == -1:
            i = self.rank 
        qn = qn if qn is not None else self.qsp_class.QnClass.qn_id()
        if isinstance(qn, int): 
            qn = self.qsp_class.QnClass(qn)
            
        qsp = self.qsp_class(1, [qn], [1])
        qnr = qn.copy()
        #qnr.reverse()
        leg = iTensor(QSp=[qsp], totQN=qnr)
        leg.data[0] = 1.0
        res = self.contract(leg, list(range(self.rank)), [self.rank])
        order = list(range(self.rank))
        order.insert(i, self.rank)
        res= res.permutation(order)
        return res 
        
        insert_travial_ind = insert_1d_qsp
     
    def shift_qn(self, qn_delta, qsp_id): 
        """
            shift totqn and  qn of a qsp at the same time, not changing data 
            an inplace operation:
                self.totQN -> self.totQN + qn_delta 
                self.QSp[qsp_id].QNs + qn_delta 
            this is efficient enough for production use 
            params:
                qn_delta: of type QnU1, etc.
        """
        #qsp_delta = self.qsp_class(1, [qn_delta], [1])
        #self.QSp[qsp_id] = self.QSp[qsp_id]*qsp_delta 
        self.QSp[qsp_id].shift_qn(qn_delta)
        #self.totQN = self.totQN + qn_delta
        self.totQN = self.totQN.__add__(qn_delta) 
    
    @staticmethod 
    def get_player_status():
        """
            tensor_player status
        """
        res = {}
        if 'single' in tensor_player.__module__:
            name_list = ['permutation', 'contract_core']
            for i in name_list:
                meth = getattr(iTensor, i)
                #print_vars(vars(),  ['meth.__closure__[4].cell_contents'])
                inner = meth.__closure__[1].cell_contents
                calls_tot = getattr(inner, 'calls_tot')
                res[i] = calls_tot
        
        return res  
    
    @staticmethod
    def reset_player():
        """
            when tensor_player.STATE changing 
                play    -> record   √
                record  -> record   may cause problem
                stop    -> record   may cause problem 
            in some rare and hard cases, it is not able to automatically clear
            the tape before 'record', then I may force reset tensor player
            through calling this function
            
        """
        if tensor_player.version == 'single':
            name_list = [  'set_data_entrance', 'contract_core', 'permutation' ]
            for i in name_list:
                meth = getattr(iTensor, i)
                reset = meth.__closure__[4].cell_contents
                reset.__call__()
        elif tensor_player.version == 'multiple':
            raise NotImplemented
   
    def save(self, path):
        save(self, path)
        print('itensor saved')
    
    def load(path):
        return load(path)
    
    def diagonal(self):
        """
            only suitable for rank-2 tensors. This is mainly for debug. 
        
        """
        from merapy.tensor_svd import iTensor_rank2_operation
        assert self.rank == 2 
        return iTensor_rank2_operation.diag_rank2(self)
    
    def change_device(self, which):
        assert which in ('cpu', 'gpu')
        if self.device == which:
            return 
        if which == 'gpu' :
            self.data = cp.asarray(self.data)
            self.use_gpu = 1
        elif which == 'cpu' :
            self.data = cp.asnumpy(self.data)
            self.use_gpu = 0
        else:
            raise ValueError 
    
    
    @property
    def device(self):
        return 'cpu' if isinstance(self.data, np.ndarray) else 'gpu'

class performance_iTensor(object):
    def __init__(self, symmetry):
        self.symmetry = symmetry
        pass
    
    def _permute(self, symmetry="U1", rank=8, dim=4, nqn=None,  NUM_OF_THREADS=8, iter_times=1000):
        import os
        os.environ["OMP_NUM_THREADS"] = str(NUM_OF_THREADS)
        os.environ["OMP_SCHEDULE"] = "static"#"dynamic" 
        import time
        t = iTensorFactory.simple(rank=rank, dim=[dim]*rank, symmetry=symmetry, nqn=nqn)
        buf = np.ndarray(t.data.size)
        t.data[:] = np.arange(t.data.size)
        r = rank//2
        tensor_player.STATE = "record"
        t1=t.permutation(list(range(r,rank)) + range(r), buffer=buf) 
        tensor_player.STATE = "play"
        def func(which):
            array_permutation.permute_player_fort = \
                    array_permutation.__getattribute__("permute_player_fort" + which)
            #print "do schedule dynamic"
            t0 = time.clock(); ta = time.time()
            for i in range(iter_times):
                t1=t.permutation(list(range(r,rank)) + range(r), buffer=buf)
            tb = time.time(); t1 = time.clock()
            if which == "_parallel_runtime":
                #print os.environ["OMP_SCHEDULE"] 
                which  = which +  "  " + os.environ["OMP_SCHEDULE"] 
            print(which, "\t", t1-t0, tb-ta)


        
        import common_64_ifort as c64
        c64.set_num_of_threads(NUM_OF_THREADS)

        #func("")
        #func("_parallel")        
        func("_parallel_dynamic")
        #func("_parallel_guided")
        #func("_parallel_runtime")

    def permute(self):
        """
            report:
                for all trunc_dim parallel_dynamic is always faster than non-parallel
                for trunc_dim <= 8,  parallel faster than non-parallel
        """
        for n in range(1, 8):
            self._permute(symmetry="U1", rank=8, dim=13, nqn=5, NUM_OF_THREADS=n, iter_times=10)


    def _contract(self, symmetry="U1", rank1=4, rank2=4, dim=4, 
            nqn=None,  NUM_OF_THREADS=4, iter_times=1000):
        import os
        import common_64_ifort as c64
        import time
        os.environ["OMP_NUM_THREADS"] = str(NUM_OF_THREADS)
        #c64.set_num_of_threads(NUM_OF_THREADS)
        os.environ["OMP_SCHEDULE"] = "static"#"dynamic"#"dynamic" 
        
        T1 = iTensorFactory.simple(rank=rank1, dim=[dim]*rank1, symmetry=symmetry, nqn=nqn)
        T2 = iTensorFactory.simple(rank=rank2, dim=[dim]*rank2, reverse=[0, 1], symmetry=symmetry, nqn=nqn)
        #T2.QSp[0].reverse()
        #T2.QSp[1].reverse()
        #buf = np.ndarray(t.data.size)
        #t.data[:] = np.arange(t.data.size)
        tensor_player.STATE = "record"
        T3=T1.contract_core(T2, 2)
        tensor_player.STATE = "play"
        #print "eeee"; exit()
        buff = T3.data
        
        def func(which):
            if which != "":
                common_util.contract_core_player_fort = \
                        c64.__getattribute__("contract_core_player_fort" + which)
            #print "do schedule dynamic"
            t0 = time.clock(); ta = time.time()
            for i in range(iter_times):
                T1.contract_core(T2, 2, data=buff)
            tb = time.time(); t1 = time.clock()
            if which == "_parallel_runtime":
                #print os.environ["OMP_SCHEDULE"] 
                which  = which +  "  " + os.environ["OMP_SCHEDULE"] 
            print(which, "\t", t1-t0, tb-ta)

        func("")
        #func("_paralell_ordered")
        #func("_paralell_critical")
        func("_paralell_test")
        #func("_paralell_reduction")
        #func("_paralell_reduction_1")
    
    def contract(self):
        self._contract(symmetry="U1", rank1=8, rank2=4, dim=16, nqn=3,  
                NUM_OF_THREADS=6, iter_times=10)

class Test_iTensor(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_tensor_player_single(self): 
        if tensor_player.version == 'multiple':
            return 
        for i in range(10):
            print_vars(vars(),  ['i'], '', ' ')
            rank = 4
            set_player_state_auto(iter=i, record_at=0, info=1)    
            #set_STATE_end_1(iter=i, record_at=0, stop_at=10000000, power_on=True) 
            u = iTensor.example()
            #u.contract_core(u, 2)
            u.permutation([1, 3, 2, 0])
        #tensor_player.STATE = 'stop'
        set_player_state_manual('stop')
        
        for i in range(1, 2):
            set_player_state_auto(iter=i, record_at=1, info=1)    
            t1 = iTensor.example(rank=4)
            t2 = iTensor.example(rank=4)
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3.permutation([0, 2, 3, 1, 4, 5])
            t3.permutation([0, 2, 3, 1, 4, 5])
        #tensor_player.STATE = 'stop'
        set_player_state_manual('stop')
        
        status = iTensor.get_player_status()
        print_vars(vars(),  ['status'])
        
        
        iTensor.reset_player()
        for i in range(1, 5):
            set_player_state_auto(iter=i, record_at=1, info=1)    
            t1 = iTensor.example(rank=4)
            t2 = iTensor.example(rank=4)
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3.permutation([0, 2, 3, 1, 4, 5])
            t3.permutation([0, 2, 3, 1, 4, 5])
        #tensor_player.STATE = 'stop'
        set_player_state_manual('stop')

        status = iTensor.get_player_status()
        print_vars(vars(),  ['status'])
    
    def test_tensor_player_multiple(self): 
        if tensor_player.version == 'single':
            return 
        if 1:
            print(iTensor.contract_core)
            print(iTensor.__init__)
        
        for i in range(10):
            print_vars(vars(),  ['i'], '', ' ')
            rank = 4
            set_player_state_auto(iter=i, record_at=0, info=1)    
            #set_STATE_end_1(iter=i, record_at=0, stop_at=10000000, power_on=True) 
            u = iTensor.example()
            #u.contract_core(u, 2)
            u.permutation([1, 3, 2, 0])
        #tensor_player.STATE = 'stop'
        print((type(tensor_player.the_tape)))
        print_vars(globals(),  ['tensor_player.the_tape.values()'])
        print(tensor_player.the_tape.keys())
        print(tensor_player.the_tape.get(1))
        set_player_state_manual('stop')
       
        for i in range(1, 4):
            set_player_state_auto(iter=i, record_at=1, info=1)    
            t1 = iTensor.example(rank=4)
            t2 = iTensor.example(rank=4)
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3.permutation([0, 2, 3, 1, 4, 5])
            t3.permutation([0, 2, 3, 1, 4, 5])
        #tensor_player.STATE = 'stop'
        set_player_state_manual('stop')
        
        tensor_player.the_tape.reset()   #dont remove this line, otherwise, other tests may fail
        
        #status = iTensor.get_player_status()
        #print_vars(vars(),  ['status'])
        if 0:
            pass
            if tensor_player.version == 'single':
                return 
            if 1:
                print(iTensor.contract_core)
                print(iTensor.__init__)
            
            for i in range(10):
                print_vars(vars(),  ['i'], '', ' ')
                rank = 4
                set_player_state_auto(iter=i, record_at=0, verbose=1, info=1)    
                #set_STATE_end_1(iter=i, record_at=0, stop_at=10000000, power_on=True) 
                u = iTensor.example()
                #u.contract_core(u, 2)
                u.permutation([1, 3, 2, 0])
                
            #tensor_player.STATE = 'stop'
            print((type(tensor_player.the_tape)))
            print((tensor_player.the_tape))
            raise  
            print_vars(globals(),  ['tensor_player.the_tape.values()'])
            print(tensor_player.the_tape.keys())
            print(tensor_player.the_tape.get(1))
            set_player_state_manual('stop')
            
        
        
        def test_get_element_and_set_element():
            #print w.__repr__()  #['Addr_idx']
            #这么做不对
            #w244 = w.copy()
            #w244.QSp = [QSbase(), QSbase().add(QSbase()), QSbase().add(QSbase())]
            #print "w244", w244.QSp
            qsp = [QSbase(), QSbase().add(QSbase()), QSbase().add(QSbase())]
            w244 = iTensor(3, qsp, totQN)
            
            qDims= [1, 0, 1]
            w244.set_element(qDims, [0, 1, 1], 5.)
            x = w244.get_element(qDims, [0, 1, 1])
            print(x)
            print(w244.__repr__(["data"]))
           
    def test_permutation(self): 
        t = iTensor.example()
        t.permutation([2, 1, 0, 3])
        
        #complex dtype 
        t = iTensor.example(rank=2, dtype=np.complex128)
        tp=t.permutation([1, 0])
        np.set_printoptions(4)
        print_vars(vars(),  ['t.to_ndarray()', 'tp.to_ndarray()'], '', ' ')
    
    def test_contract(self):
        from merapy import qsp_any
        if 0:
            q0 = qsp_any('U1', [1, -1], [3, 4])
            q1 = qsp_any('U1', [1, -1], [2, 1])
            qsp = [q0, q1]
            t = iTensor(QSp=qsp)
            t.data[:] = np.arange(t.data.size)
            t1 = t.conj()
            res = t.contract(t1, [0, 1, ], [0, 5])
            print_vars(vars(),  ['res.data'])
            raise  
            #res_old = np.asarray([245100., 262860., 280620., 298380., 316140., 262860., 296677., 330494., 364311., 398128.])
            #self.assertTrue(np.all(res.data[:10]==res_old))
       
        if 1:
            q = qsp_any('U1', [0, 1, -1], [5, 3, 2])
            qsp = q.copy_many(3)
            
            t = iTensor(QSp=qsp)
            t.data[:] = np.arange(t.data.size)
            t1 = t.conj()
            res = t.contract(t1, [0, 1, 3], [1, 0, 4])
            res_old = np.asarray([245100., 262860., 280620., 298380., 316140., 262860., 296677., 330494., 364311., 398128.])
            #print_vars(vars(),  ['res.data[:10]'])
            self.assertTrue(np.all(res.data[:10]==res_old))
        
    
    def test_to_ndarray(self): 
        q = QspZ2.easy_init([1, -1], [2, 2])
        qsp = q.copy_many(2)
        t = iTensor(QSp=qsp); t.data[: ]=np.arange(t.size)
        tn=t.to_ndarray()
        temp=[[ 0,  2,  0,  0,], 
             [ 1,  3,  0,  0,], 
             [ 0,  0,  4,  6,], 
             [ 0,  0,  5,  7,]]
        self.assertTrue(np.all(tn==np.asarray(temp)))
        
        if 0: 
            pau = iTensorFactory.pauli_mat( 'Z2')
            temp = ['sx', 'sy', 'sz']
            for i in temp: 
                t = pau[i]
                tn = t.to_ndarray()
                print(i) 
                print(tn) 
                print(t.matrix_view())
        
    def test_rank_zero(self): 
        print('aaaaaaaaaaaaaaaaaaaa')
        QSp = []
        t0 = iTensor(rank=0,  QSp=QSp, totQN=QspU1.QnClass.qn_id())
        print(t0) #.data 
    
    def test_rank_zero_1(self): 
        V1=[0,1,2,3]
        V2=[2,3,0,1]
        qsp = QspZ2.easy_init([1, -1], [2, 2])
        qq = qsp.copy_many(4)
        oo = iTensor(4, qq)
        rho_2 = iTensor(4, qq)
        #res, legs = oo.contract(rho_2, V1, V2, use_buf=True)  #rho and H_2 fully contracted to a scalar
        res=oo.contract_core(rho_2, 4, use_buf=True)  #rho and H_2 fully contracted to a scalar
        res= res.data[0]

        pass 
    
    def test_matrix_view(self): 
        pass 
    
    def test_split_2to3(self): 
        if 1:
            q0 = qsp_any('U1', [0], [14])
            q1 = qsp_any('U1', [0], [4])
            t2 = iTensor(QSp=[q0, q1])
            t2.data[: ] = np.arange(t2.size)
            q1a = qsp_any('U1', [1, 2, 3],    [1, 2, 3])
            q1b = qsp_any('U1', [-1, -2, -3], [1, 2, 3])
            #print_vars(vars(),  ['t2.matrix_view()'])
            t3 = t2.split_2to3(0, [q1a, q1b])            
            res_old=np.asarray([ 0., 14., 28., 42.,  1.,  2.,  3.,  4., 15., 16., 17., 18., 29., 30., 31., 32., 43., 44., 45., 46.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 19., 20., 21., 22., 23., 24., 25., 26., 27., 33., 34., 35., 36., 37., 38., 39., 40., 41., 47., 48., 49., 50., 51., 52., 53., 54., 55.])               
            self.assertTrue(np.all(t3.data==res_old))
        
        if 1:  
            q0 = QspZ2.easy_init([1, -1], [8, 8])
            q1 = QspZ2.easy_init([1, -1], [4, 2])
            t2 = iTensor(QSp=[q0, q1])
            t2.data[: ] = np.arange(t2.size)
            t2.show_data()
            c2 = t2.contract(t2, [0, 1], [0, 2])
            print_vars(vars(), ['c2.data', 'c2.matrix_view()',  't2_mat_simple.T.dot(t2_mat_simple)'])
            
            qa = QspZ2.easy_init([1, -1], [2, 2])
            qb = QspZ2.easy_init([1, -1], [2, 2])
            
            t3 = t2.split_2to3(0, [qa, qb])
            c3 = t3.contract(t3, [0, 1, 2], [0, 1, 3])
            c3.show_data()
            print(c2.data)
            print(c3.data) 
            self.assertTrue(np.all(c2.data==c3.data))
        if 1:  
            q0 = QspZ2.easy_init([1, -1], [2, 5])
            q1 = QspZ2.easy_init([1, -1], [8, 8])
            t2 = iTensor(QSp=[q0, q1])
            t2.data[: ] = np.arange(t2.size)
            t2.show_data()
            c2 = t2.contract(t2, [0, 1], [2, 1])
            print_vars(vars(), ['c2.data', 'c2.matrix_view()',  't2_mat_simple.T.dot(t2_mat_simple)'])
            
            qa = QspZ2.easy_init([1, -1], [2, 2])
            qb = QspZ2.easy_init([1, -1], [2, 2])
            
            t3 = t2.split_2to3(1, [qa, qb])
            c3 = t3.contract(t3, [5, 1, 2], [3, 1, 2])
            c2.show_data()
            c3.show_data()
            self.assertTrue(np.all(c2.data==c3.data))
    
    def test_merge_3to2(self): 
        if 1:  
            q0 = QspZ2.easy_init([1, -1], [2, 4])
            q1 = QspZ2.easy_init([1, -1], [2, 3])
            q2 = QspZ2.easy_init([1, -1], [3, 2])
            
            t3 = iTensor(QSp=[q0, q1, q2])
            t3.data[: ] = np.arange(t3.size)
            t2 = t3.merge_3to2((0, 1))
            
            c3 = t3.contract(t3, [0, 1, 2], [0, 1, 3])
            c2 = t2.contract(t2, [0, 1], [0, 2])
            
            print(c2.data)
            print(c3.data) 
            self.assertTrue(np.all(c2.data==c3.data))
        
        if 1:  
            q0 = QspU1.easy_init([1, -1], [2, 4])
            q1 = QspU1.easy_init([1, -1], [2, 3])
            #q2 = QspU1.easy_init([1, -1], [3, 2])
            q2 = q0*q1; q2.reverse()
            
            t3 = iTensor(QSp=[q0, q1, q2])
            t3.data[: ] = np.arange(t3.size)
            t2 = t3.merge_3to2((0, 1))
            u3 = t3.copy(); u3.reverse_qsp()
            c3 = t3.contract(u3, [0, 1, 2], [0, 1, 3])
            u2 = t2.copy(); u2.reverse_qsp()
            c2 = t2.contract(u2, [0, 1], [0, 2])
           
            self.assertTrue(np.all(c2.data==c3.data))
          
        if 1:  
            q0 = QspZ2.easy_init([1, -1], [2, 5])
            q1 = QspZ2.easy_init([1, -1], [2, 2])
            q2 = QspZ2.easy_init([1, -1], [2, 3])
            t3 = iTensor(QSp=[q0, q1, q2])
            t3.data[: ] = np.arange(t3.size)
            c3 = t3.contract(t3, [5, 1, 2], [3, 1, 2])
            
            t2 = t3.merge_3to2((1, 2))
            c2 = t2.contract(t2, [0, 1], [2, 1])
            c2.show_data()
            c3.show_data()
            self.assertTrue(np.all(c2.data==c3.data))
        
        if 1:  
            q0 = QspU1.easy_init([0, 1, -1], [1, 2, 5])
            q1 = QspU1.easy_init([0, 1, -1], [2, 1, 2])
            q2 = QspU1.easy_init([0, 1, -1], [3, 2, 3])
            t3 = iTensor(QSp=[q0, q1, q2])
            t3.data[: ] = np.arange(t3.size)
            u3 = t3.copy(); u3.reverse_qsp()
            c3 = t3.contract(u3, [5, 1, 2], [3, 1, 2])
            
            t2 = t3.merge_3to2((1, 2))
            u2 = t2.copy(); u2.reverse_qsp()
            c2 = t2.contract(u2, [0, 1], [2, 1])
            c2.show_data()
            c3.show_data()
            self.assertTrue(np.all(c2.data==c3.data))
        
    def test_split_qsp(self): 
        if 1: 
            qa = QspZ2.easy_init([1, -1],  [1, 1])
            qb = QspZ2.easy_init([1, -1],  [3, 3])
            qc = QspZ2.easy_init([1, -1],  [2, 2])
        else:  
            qa = QspU1.easy_init([0, 1, -1], [2, 1, 1])
            qb = QspU1.easy_init([ 1, -1], [3, 3])
            qc = QspU1.easy_init([0, 1, -1], [2, 1, 1])
        if 1: 
            t = iTensor(QSp=[qa, qb*qc]) 
            t.data[:] = np.arange(t.size)
            tt=t.split_qsp(1, [qb, qc])
            ttt=t.split_2to3(1, [qb, qc])
            print(t.data) 
            print(tt.data) 
            print(ttt.data) 
            self.assertTrue(np.all(tt.data==ttt.data))
        if 1: 
            t = iTensor(QSp=[qa*qb, qc]) 
            t.data[:] = np.arange(t.size)
            tt=t.split_qsp(0, [qa, qb])
            ttt=t.split_2to3(0, [qa, qb])
            print(t.data) 
            print(tt.data) 
            print(ttt.data) 
            self.assertTrue(np.all(tt.data==ttt.data))
        if 1: 
            t = iTensor(QSp=[qa*qb]) 
            t.split_qsp(0, t.QSp+[t.qsp_class.null()])
            t.reshape(t.QSp + [t.qsp_class.null()] )
        
    def test_split_qsp_2(self): 
        #qsp_class= QspZ2 
        qsp_class= QspU1 
        q = qsp_class.easy_init([1, -1],  [1, 1])
        qsp = q.copy_many(6)
        q0 = qsp_class.prod_many(qsp[: 3])    
        q1 = qsp_class.prod_many(qsp[3: ])    
        t = iTensor(QSp=[q0, q1])
        t.data[: ] = np.arange(t.size)
        tt = t.split_qsp(0, qsp[: 3], 1, qsp[3: ])
        print(tt.data.tolist())
        old = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 11.0, 12.0, 7.0, 8.0, 9.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0]
        self.assertTrue(np.all(tt.data==old))
    
    def test_split_qsp_by_contract(self): 
        if 1: 
            #qsp_class= QspU1
            qsp_class= QspZ2
            qa = qsp_class.easy_init([1, -1],  [2, 2])
            qb = qsp_class.easy_init([1, -1],  [3, 2])
            qc = qsp_class.easy_init([1, -1],  [2, 3])
            qd = qsp_class.easy_init([1, -1],  [2, 3])
            qe = qsp_class.easy_init([1, -1],  [2, 3])
        
        if 1:  # pass  
            t = iTensor(QSp=[qa*qb, qc*qd]); t.data[:] = np.arange(t.size)
            t2=t.split_qsp(0, [qa, qb], 1, [qc, qd])
            u = t.copy(); u.reverse_qsp();  u2 = t2.copy(); u2.reverse_qsp() 
            c2 = t.contract(u, [0, 1], [0, 1])
            c3 = t2.contract(u2, [0, 1, 2, 3], [0, 1, 2, 3])
            self.assertTrue(np.all(c2.data==c3.data))

        if 1:  #pass 
            t = iTensor(QSp=[qa*qb, qc, qd*qe]); t.data[:] = np.arange(t.size)
            t2=t.split_qsp(0, [qa, qb], 2, [qd, qe])
            c2 = t.contract(t, [0, 100, 1], [0, 1000, 1])
            c3 = t2.contract(t2, [0, 1, 100, 2, 3], [0, 1, 1000, 2, 3])
            c2.show_data()
            c3.show_data()
            self.assertTrue(np.all(c2.data==c3.data))
        
        if 1: 
            t = iTensor(QSp=[qa*qb, qc*qd, qe]); t.data[:] = np.arange(t.size)
            t2=t.split_qsp(0, [qa, qb], 1, [qc, qd])
            c2 = t.contract(t, [0, 1, 100], [0, 1, 1000])
            c3 = t2.contract(t2, [0, 1, 2, 3, 100], [0, 1, 2, 3, 1000])
            self.assertTrue(np.all(c2.data==c3.data))
       
        if 1:   #pass 
            t = iTensor(QSp=[qa*qb*qc*qd, qe]); t.data[:] = np.arange(t.size)
            t2=t.split_qsp(0, [qa, qb, qc, qd])
            c2 = t.contract(t, [0, 2], [0, 1])
            c3 = t2.contract(t2, [0, 1, 2, 3, 100], [0, 1, 2, 3, 1000])
            self.assertTrue(np.all(c2.data==c3.data))
    
    def test_merge_qsp(self): 
        if  0: 
            qsp_class= QspZ2 
            #qsp_class= QspU1 
            q = qsp_class.easy_init([1, -1],  [2, 2])
            qsp = q.copy_many(3)
            t = iTensor(QSp=qsp)
            t.data[: ] = np.arange(t.size)
            tt=t.merge_qsp((0, 1))
            #tt=t.merge_qsp((1, 2))
            print(t.data)
            print(tt.data) 
        if  1: 
            qsp_class= QspU1 
            q = qsp_class.easy_init([1, -1],  [2, 2])
            qsp = q.copy_many(4)
            t = iTensor(QSp=qsp)
            t.data[: ] = np.arange(t.size)
            tt=t.merge_qsp((0, 1))
            #tt=t.merge_qsp((1, 2))
            #t.show_data() 
            #tt.show_data()
            print_vars(vars(), ['t', 'tt'])
    
    def test_merge_qsp_2(self): 
        for qsp_class in [QspZ2, QspU1]:
            q = qsp_class.easy_init([1, -1],  [1, 1])
            qsp = q.copy_many(6)
            t = iTensor(QSp=qsp)
            t.data[: ] = np.arange(t.size)
            
            tt = t.merge_qsp((0, 1), (2, 3, 4), )
            t3 = tt.split_qsp(0, [q, q], 1, [q, q, q])
            t4 = t3.merge_qsp((0, 1), (2, 3, 4), )
            #print t.data
            #print t3.data 
            self.assertTrue(np.all(t.data==t3.data))
    
    def test_merge_qsp_by_contract(self): 
        if 1: 
            #qsp_class = QspZ2
            qsp_class = QspU1
            qa = qsp_class.easy_init([1, -1],  [2, 2])
            qb = qsp_class.easy_init([1, -1],  [3, 2])
            qc = qsp_class.easy_init([1, -1],  [2, 3])
            qd = qsp_class.easy_init([1, -1],  [2, 3])
            qe = qsp_class.easy_init([1, -1],  [2, 3])
        t3 = iTensor(QSp=[qa, qb, qc, qd, qe]); t3.data[:] = np.arange(t3.size)
        if 1:  # pass  
            t2 = t3.merge_qsp((0, 1), (3, 4))
            c2 = t2.contract(t2, [0,100,  1], [0, 1000,  1])
            c3 = t3.contract(t3, [0, 1, 100,  2, 3], [0, 1, 1000,  2, 3])
            self.assertTrue(np.all(c2.data==c3.data))

        if 1: 
            t2=t3.merge_qsp((0, 1), (2, 3))
            c2 = t2.contract(t2, [0, 1, 100], [0, 1, 1000])
            c3 = t3.contract(t3, [0, 1, 2, 3, 100], [0, 1, 2, 3, 1000])
            self.assertTrue(np.all(c2.data==c3.data))
       
        if 1:   #pass 
            t = iTensor(QSp=[qa*qb*qc*qd, qe]); t.data[:] = np.arange(t.size)
            t2=t.split_qsp(0, [qa, qb, qc, qd])
            c2 = t.contract(t, [0, 2], [0, 1])
            c3 = t2.contract(t2, [0, 1, 2, 3, 100], [0, 1, 2, 3, 1000])
            self.assertTrue(np.all(c2.data==c3.data))

    def test_reshape(self): 
        if 1: 
            qa = QspZ2.easy_init([1, -1],  [2, 2])
            qb = QspZ2.easy_init([1, -1],  [3, 2])
            qc = QspZ2.easy_init([1, -1],  [2, 3])
            qd = QspZ2.easy_init([1, -1],  [2, 3])
            qe = QspZ2.easy_init([1, -1],  [2, 3])
        if 1:    
            t3 = iTensor(QSp=[qa, qb, qc, qd, qe]); t3.data[:] = np.arange(t3.size)
            a=t3.reshape([qa*qb, qc, qd*qe] )
            b=t3.merge_qsp((0, 1), (3,4))
            self.assertTrue(np.all(a.data==b.data))
        
        if 1:    
            t2 = iTensor(QSp=[qa*qb, qc, qd*qe]); t2.data[:] = np.arange(t2.size)
            a=t2.reshape([qa, qb, qc, qd, qe] )
            b=t2.split_qsp(0, (qa, qb), 2, (qd, qe))
            self.assertTrue(np.all(a.data==b.data))
        if 1: 
            t = iTensor(QSp=[qa*qb]) 
            t.split_qsp(0, t.QSp+[t.qsp_class.null()])
            t.reshape(t.QSp + [t.qsp_class.null()] )
            
    def test_reshape_u1(self): 
        if 1:
            null = QspU1.easy_init([0], [1])
            d = QspU1.easy_init([1, -1],  [1, 1])
            D = null*d*d; D.reverse()
            qsp = [null, d.copy(), d.copy(), D]
            t3 = iTensor(QSp=qsp)
            t2 = t3.reshape(null, d**2, D.copy())
            u3 = t3.copy(); u3.reverse_qsp()
            u2 = t2.copy(); u2.reverse_qsp()
            c3 = t3.contract(u3, [1, 2, 3, 4], [8, 2, 3, 10])
            c2 = t2.contract(u2, [1, 2, 4], [8, 2, 10])
            self.assertTrue(np.all(c3.data==c2.data))   

        if 1:
            Dl = QspU1.easy_init([0], [1])
            d = QspU1.easy_init([1, -1],  [1, 1])
            Dr = Dl*d*d; Dr.reverse()
            qsp = [Dl, Dr, d**2]
            A = iTensor(QSp=qsp)
            A1 = A.reshape(Dl, Dr, d, d)
            B = A1.transpose((0, 2, 1, 3))
            #print_vars(vars(), ['A.QSp', 'A1.QSp', 'B.QSp', 'C.QSp'])
            print_vars(vars(), ['d*d', 'A.QSp', 
                'A.reshape(Dl, Dr, d, d).QSp', 
                'A.reshape(Dl, Dr, d, d).transpose((0, 2, 1, 3)).QSp', 
                'Dl*d', 'Dr*d', 
                #'A.reshape(Dl, Dr, d, d).transpose((0, 2, 1, 3)).reshape(Dl*d, Dr*d).QSp', 
                ])
            C = B.reshape(Dl*d, Dr*d)
            C = B.merge_qsp((0, 1), (2, 3))
    
    def test_conj_new(self): 
        #t = iTensor.example(symmetry='U1')
        np.random.seed(1234)
        a = QspU1.easy_init([1, -1], [1, 1])
        b = QspU1.easy_init([0, 1, -1], [1, 2, 2])
        c = QspU1.easy_init([0, 1, -1], [1, 3, 3])
        t = iTensor(QSp=[a, b, c])
        t.data = np.random.random(t.size)
        
        t1=t.conj_new(2)
        t1_old = np.array([ 0.19151945,  0.62210877,  0.43772774,  0.78535858,  0.80187218, 0.95813935,  0.87593263,  0.77997581,  0.27259261,  0.27646426])
        print_vars(vars(),  ['t.shape', 't1.shape', 'repr(t1.data)', ]) 
        self.assertTrue(np.allclose(t1.data, t1_old, atol=1e-8))
        #res_old = array([ 0.8980255 ,  0.3927732 ,  0.06395552,  0.37972657,  0.23755964, 0.84733446])
        
        t1=t.conj_new(0)
        print_vars(vars(),  ['t.shape', 't1.shape', 'repr(t1.data)', ]) 
        t1_old = np.array([ 0.19151945,  0.62210877,  0.43772774,  0.78535858,  0.77997581, 0.27259261,  0.27646426,  0.80187218,  0.95813935,  0.87593263])
        self.assertTrue(np.allclose(t1.data, t1_old, atol=1e-8))
    
    def test_reduce_and_insert_1d_qsp(self): 
        pass
        if 1: 
            q1 = QspU1.easy_init([1, -1], [1, 1]) 
            q3 = QspU1.easy_init([2], [1]) 
            qsp = q1.copy_many(2)  + [q3]  + q1.copy_many(2)
            t = iTensor(QSp= qsp)
            t.data[: ] = np.random.random(t.size)
            t1=t.reduce_1d_qsp(2)
            self.assertTrue(t1.totQN._val==-2)
        if 1: 
            t2 = t1.insert_1d_qsp(2, q3.QNs[0]._val)
            print_vars(vars(),  ['t.data', 't2.data'])
            self.assertTrue(np.allclose(t.data, t2.data))
            print_vars(vars(),  ['t.shape', 't2.shape'])
            self.assertTrue(t.shape==t2.shape)

    def test_itensor_gpu(self):
        from merapy.tensor_svd import Tensor_svd
        try:
            type(cp)
        except:
            print('cupy not installed, return')
            return 
        Dl = qsp_any('U1', qns=[0, 1, -1, ], dims=[10, 5, 5, ])
        Dr = Dl.conj()
        d = qsp_any('U1', qns=[1, -1], dims=[1, 1])
        qsp = [Dl, Dr, d]
        use_gpu = 1
        #construction
        t = iTensor(QSp=qsp, use_gpu=use_gpu)
        t.data = cp.random.random(t.data.size)-0.5
        assert t.device == 'gpu' 
        
        
        #contract
        tt = t.contract(t.conj(), (0, 1, 2), (3, 1, 2))
        
        #transpose
        tp = t.transpose((0, 2, 1))
        self.assertTrue(tp.device=='gpu')
        
        #merge and split index
        t2 = t.merge_qsp((1, 2))
        t3 = t2.split_qsp(1, (Dr, d))
        self.assertTrue(t==t3)
        
        #svd
        u, s, v =Tensor_svd.svd_rank2(tt)
        print_vars(vars(),  ['tt.sh'])
        temp = u.dot(s).dot(v)
        #print_vars(vars(),  ['temp==tt'])
        self.assertTrue(temp.is_close_to(tt))
    
    def dev_test_su2_symm(self):
        pass
            
        if 0:
            q0 = qsp_any('SU2', qns=[(0, 0), (1, -1), (1, 0), (1, 1)], dims=[2, 2, 2, 2])
            q1 = qsp_any('SU2', qns=[(0, 0), (1, -1), (1, 0), (1, 1)], dims=[2, 2, 2, 2])

            t = iTensor(QSp=[q0, q1])
            print_vars(vars(),  ['t'])
            
        if 0:  #rank-3
            q0 = qsp_any('SU2', qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[2, 2, 2, 1])
            q1 = qsp_any('SU2', qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[2, 2, 2, 1])
            q2 = qsp_any('SU2', qns=[(0, 0), (2, 0),  (2, 2)], dims=[2, 2, 1])
            if 0:
                j0, j1, m0, m1 = 2, 0, -2, 0
                j2, m2 = 2, -2 
                cg = clebsch_gordan(j0, j1, j2, m0, m1, m2)
                print_vars(vars(),  ['cg'])
                raise  

            t = iTensor(QSp=[q0, q1, q2])
            print_vars(vars(),  ['t'])
            print_vars(vars(),  ['t.cg_coeff'])
            
        if 0:  #rank-4
            q0 = qsp_any('SU2', qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[2, 2, 2, 1])
            q1 = qsp_any('SU2', qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[2, 2, 2, 1])
            q2 = qsp_any('SU2', qns=[(0, 0), (2, 0),  (2, 2)], dims=[2, 2, 1])
            q3 = qsp_any('SU2', qns=[(0, 0), (2, 0),  (2, 2)], dims=[2, 2, 1])

            t = iTensor(QSp=[q0, q1, q2, q3])
            #print_vars(vars(),  ['t'])
            print_vars(vars(),  ['t.cg_coeff'])
    
    def test_temp(self): 
        warnings.filterwarnings("ignore")
        
        from merapy.tensor_svd import Tensor_svd

        dtype= np.float32
        print_vars(vars(),  ['dtype'])

        for n in [10, 50, 100, 500, 1000, 2000]:
        #for n in [8, 9, 32, 33, 64, 65, 128, 129, 256, 257, 512, 513, 1024, 1025]:
            dims  =  [2.7, 1.2, 1]
            _dims = [int(n*i) for i in dims]
            Dl = qsp_any('U1', qns=[0, 1, -1, ], dims=_dims)

            dims  =  [3, 1.2, 1]
            _dims = [int(n*i) for i in dims]
            Dr = qsp_any('U1', qns=[0, 1, -1, ], dims=_dims)
            d = qsp_any('U1', qns=[1, -1], dims=[1, 1])
            qsp = [Dl, Dr, d]
            
            print(f'\nqsp is Dl={Dl}')
            
            
            #for use_gpu in [0,  1,  2]:
            for use_gpu in [0,  1,  ]:
                #iTensor.USE_GPU_FOR_BLOCK = use_gpu
                t = iTensor(QSp=qsp, use_gpu=use_gpu, dtype=dtype)
                t.set_data_random()
                tc = t.conj()
                t0 = time.time()
                for i in range(1):
                    t.contract(tc, (0, 1, 2), (3, 1, 2), use_buf=1)
                t1 = time.time()

                print(f'time for iTensor contraction on {use_gpu}: {t1-t0}')
            
        
            
        
        
        
    
        
            


if __name__ == "__main__":
    #warnings.filterwarnings("ignore")
    if 1: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list_iTensor = [
           #'test_permutation', 
           #'test_contract', 
           #'test_to_ndarray', 
           #'test_rank_zero', 
           #'test_rank_zero_1', 
           #
           #'test_matrix_view', 
           #
           #'test_split_2to3', 
           #'test_split_qsp', 
           #'test_split_qsp_2', 
           #'test_split_qsp_by_contract', 
           
           #'test_merge_3to2', 
           #'test_merge_qsp', 
           #'test_merge_qsp_2', 
           #'test_merge_qsp_by_contract', 
           #'test_reshape', 
           #'test_reshape_u1', 
           #
           #'test_conj_new', 
           #'test_reduce_and_insert_1d_qsp', 
           #'test_tensor_player_single', 
           #'test_tensor_player_multiple', 
           'test_itensor_gpu'
           #'dev_test_su2_symm', 
           #'test_temp', 
        ]
        
        
        for a in add_list_iTensor: 
            suite.addTest(Test_iTensor(a))
        unittest.TextTestRunner().run(suite)
       
  
