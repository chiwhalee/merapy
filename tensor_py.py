#!/usr/bin/env python
#coding=UTF8

"""  
    tensor.f90
    Improvements:
        1. 区分 leg direction
        2. 统一用ndarray 存储
        3. add label for legs

    Questions: q:
        q868: fortran program seems incorrect
        q763: deallocate, nullify difference?
        1-a. rank/2  what is rank is odd?
        2 - a. is nTensor just that normal?
        3 - a. tBuffer 暂时不管？
        4.allocated, deallocate 有必要吗
        5: -a
            在SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
                           B, LDB, BETA, C, LDC )
            中， A， B， C都是 2D矩阵，而在下面contract中， matrix_multiply(nT1.data, nT2.data, T3.data, 1.0, 0.0)
            矩阵维数不同， fortran 允许这么做？
            在unit_nTensor中同样如此
        6.-a trace means what trace?
        7. find_leg 具体是什么？
        8.-a  q144 random， how to deal with?
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
import unittest
import nose
import warnings
import pickle
import pprint
import numpy as np

#from quantum_number import *  #QuantSpace, QN_idendity, QSp_null, QSp_base
from merapy.ntensor import TensorBase, nTensor 
from merapy.quantum_number import *  #QuantumNum, QuantSpace, QN_idendity, QSp_null, QSp_base
from merapy import common_util
from merapy import array_permutation
from merapy.set1 import *
from merapy import crandom
from merapy.utilities import get_local
from merapy import decorators
from merapy.decorators import *

#import scipy.weave as weave
#from scipy import linalg
#from scipy.weave import converters        
 
#import numexpr


__all__ = ["TensorBase", "nTensor", "iTensor","iTensorFactory", "test_iTensor"]


class tBuffer(object):
    """
    the point of tBUffer is :当tensor用过自动garbage collect后，只要把tensor使用
    的buff标记in_use改为-1，而buff并未必删除，可以继续使用
    """
    def __init__(self, size, dtype):
        """
        see init_tBuffer in f90
        size: number of tensors
        """
        #self.T=[Tensor() for i in xrange(size)]
        #make self.T simpliy a place holder
        self.T=[np.ndarray(0, dtype=dtype) for i in xrange(size)]
        #for i in xrange(size): self.T[i].nullify()
        self.in_use=np.ndarray(1024,"bool")
        self.in_use[:] = False
        self.size=size   
        
    def delete(size):
        for i  in xrange(tB.size):
            tB.T[i].delete()
        self.in_use[:] = False
        sef.size = 0
        self.T = []

#meth_names= ["__init__", "set_data_entrance", "contract_core", "permutation"]
#meth_names.pop(0)

@decorate_methods(decorator=tensor_player, meth_names=None)
class iTensor(TensorBase):
    num_of_instance = 0
    #T_BUFFER=[tBuffer(size=64, dtype=TensorBase.dtype) for  i in xrange(3)]
    T_BUFFER=[tBuffer(size=100, dtype=TensorBase.dtype) for  i in xrange(4)]
    #BUFFER_ON = False
    
    def __init__(self, rank,  QSp, totQN, order="F",  
            buffer=None, use_buf=False, index_data=True, has_data=True, shallow=False):
        """
            totQN 决定了张量的变换性质 totQN = id for invariant, other for covariant
            QUESTIONG, Q:
                QSp, totQN应该在init外面copy，还是在里面copy？ copy or not this is a question
                -a.  这只是个惯例保持一致即可。但从概念上看，似乎应该在外面更好，因为初始化一个张量是一个具体的张量，那么传
                进来的QSp也应该是某"个"具体的QSp，而非某一"类"QSp。
                并且，这样做的好处是，方便在init外控制 copy，or reference

            attention_here here the default ordering of quantum number combination is changed from fortran to C 
            
            shallow: 
                only structure or with data
            self.Dims[i]:  
                在leg i 对应的空间维数
            buf_ref: 
                buf_ref[0]---which iTensor.T_BUFFER, buf_ref[1]---which tensor in the iTensor.T_BUFFER
                default -1 means not using iTensor.T_BUFFER
            issue:
                in future, QSp only searve as a way to construct iTensor, forbit it being modified, because modifying qsp alone
                is meanless and dangeous if not copyed.
                this raise a question: should qsp immutable, permanante or mutalbe?
                In this way, one only need to construct very fewer distinct QSps,  only reference, no copying
            
        """
        #comment this only for a little faster
        #TensorBase.__init__(self, rank, None)
        self.rank = rank
        self.type_name = ""
        self.index_order = order

        self.ndiv = None   #ndiv 实际是把协变反变腿分开 一个(ndiv, rank-ndiv) tensor
        self.use_buf=use_buf
        self.buf_ref = np.array([-1, -1], np.int)
        self.rank=rank
        
        #NO COPYING CONVENTION
        self.totQN = totQN   #.copy()
        self.QSp = QSp[:rank]
        #self.QSp = QSp

        if rank==0:
            #self.rank = 1
            self.Dims= np.empty(1, np.int)
            self.totDim=1
            self.QSp = QSp[:1]
            #self.QSp = QSp
            self.totQN = totQN
        else:        
            self.Dims=np.empty(rank, dtype=np.int)
        self.Dims[0] =1
        self.Dims= [QSp[i].totDim for  i in xrange(rank)]

        if index_data:
            self.set_data_entrance(order=order )

        #use_buf = False;buffer=None;warnings.warn("not use bufferrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"*100)
        if has_data:
            if buffer is not None:
                pass
            elif use_buf:   #use internal T_BUFFER; else use external buffer or no buffer
                buffer = self.buffer_assign(data_size=self.totDim)
            
            self.data = np.ndarray(self.totDim, buffer=buffer, dtype=self.dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered

    def buffer_assign(self, data_size, n=None):
        """
            see use_tBuffer in f90
            let self point to a iTensor.T_BUFFER
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
            elif self.rank <= 6:
                n = 1
            else:
                n = 2
        for i in xrange(iTensor.T_BUFFER[n].size):
            if not iTensor.T_BUFFER[n].in_use[i]:  
                self.buf_ref[0]=n
                self.buf_ref[1]=i                
                iTensor.T_BUFFER[n].in_use[i] = True 
                if iTensor.T_BUFFER[n].T[i].size<data_size:
                    iTensor.T_BUFFER[n].T[i] = np.empty(data_size, dtype=iTensor.dtype)
                
                return iTensor.T_BUFFER[n].T[i].data
        raise Exception('Error, All buffer elements are in use, stop %s\n '%(str(iTensor.T_BUFFER[0].in_use[:100], )))
    
    def set_data_entrance(self, order="F"):
        """
            实现1维稀疏存储, 对应于下式的第二个->号
            T->(D, S)->(dat, Block_ind, Addr_ind)
            this step is in effect 把Qsp 中的信息提取出来，变成更容易读取操作的信息，故
            Qsp 中包含的信息和 Block_idx, Addr_idx 等是等价的(不完全等价, 差一个对称性限制条件)，不同在存储顺序和方式
            Block_idx: is a map from idx to block info
        
        """
        rank, QSp, totQN = self.rank, self.QSp, self.totQN
        pTot =1
        for i in xrange(rank):
            pTot *= QSp[i].nQN
        self.idx_dim=pTot

        #attention_this_may_be_wrong 在python中  -1对应着最后一个元素，而Fortran中什么也不对应, 所以改成下面的
        #self.idx=np.array([int(self.idx_dim)]*self.idx_dim,"int")   #-1                
        self.idx=np.array([int(self.idx_dim)]*self.idx_dim,"int")   #-1                
        
        nidx=0
        totDim=0
        
        #iQN=[0]*self.rank
        #iQN[i]用作leg i 上的量子数 计数
        if rank == 0: 
            iQN = np.zeros(1, dtype="int")
            self.Addr_idx=np.ndarray((1,self.idx_dim),dtype="int")        
        else:
            iQN = np.zeros(self.rank, dtype="int")
            #iQN 用于给量子数组合编号
            # 实际使用的addr_inx的长度为 self.nidx
            self.Addr_idx=np.ndarray((self.rank, self.idx_dim),dtype="int")        
        
        tQN_r= self.totQN.copy()  # here must copy
        tQN_r.reverse()
        self.Block_idx=np.ndarray((3, self.idx_dim),dtype="int")

        for p in xrange(pTot):
            tQN = self.QSp[0].QNs[iQN[0]]
            #这里计算了总量子数 tQN
            
            for i in xrange(1,rank):
                tQN = tQN+self.QSp[i].QNs[iQN[i]]
            #判断量子数组合是否满足指定的对称性要求, 这个不其眼的一步实际上是核心——实现了稀疏存储
            if tQN==tQN_r:
                #q895  why not tQN == totQN?  因为operater- state duality, 从而在操作下变换正好相反？
                #其中包含协变/反变的意味
                
                d = 1
                for i in xrange(rank):
                    d = d*QSp[i].Dims[iQN[i]]
                    #print "ddddd d", i, p, d,iQN[i], len(QSp[i].Dims)
                #计算某一block的总数据量

                self.idx[p] = nidx
                #给出了0量子数组合与所有量子数组合的序号间的关系
                #self.idx 和 self.Block_idx[2]互为反函数
                #如果总量子数为0，则idx[p] =- 1(默认值) 
                #在self.block中都是记录不为0的量子数组合
                #0 is position,  1 is size, 2 is position in quantum number combinations
                self.Block_idx[0, nidx] = totDim
                #在该block在self.data中的position
                self.Block_idx[1, nidx] = d
                #block变成1d数组的长度
                self.Block_idx[2, nidx] = p
                
                #记录不为0的量子数组合，在所有量子数组合中的位置
                #Addr实为将（QN1, ..., QNn)-> ind 的映射, 将n个指标拉直了, 对没一个iTensor都定义了这个函数
                #Addr_idx这个二维数组的每一列实际上是所有非零block的量子数的编号(而不是量子数点值！)的组合
                self.Addr_idx[0,nidx] = 0 #for rank=0
                self.Addr_idx[0:rank,nidx] = iQN[0:rank]
                nidx = nidx+1                
                totDim = totDim+d
                #最终得到self.data 的总长度
            
            #遍历所有的量子数组合
            if order == "F": 
                inc = 1
                i = 0
                #attention_please  这里实际上意味着按照 fortran order 对量子数组合排序的
                while inc==1 and i<rank:
                    iQN[i] = iQN[i]+1
                    if iQN[i]<QSp[i].nQN :
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

        self.nidx = nidx
        self.totDim = totDim
    
    def unregister(self):
        """
        unregister self from iTensor.T_BUFFER if iTensor.T_BUFFER was used
        """
        #if self.use_buf and (self.buf_ref[0]!=-1):
        iTensor.T_BUFFER[self.buf_ref[0]].in_use[self.buf_ref[1]]=False
    
    def __del__(self):
        """
        task: replace __del__ with __exit__ in future

        it was said __del__ is BAD!
        但它仍然可以使用，只要__init__ 不包含 raise exception
        "If you do use __del__ make sure you are covered for any case 
        in which __init__ didn’t finish running."


        """
        if self.use_buf and (self.buf_ref[0]!=-1):
            #if iTensor != None:
            iTensor.T_BUFFER[self.buf_ref[0]].in_use[self.buf_ref[1]]=False

    @classmethod
    def buff_free(cls):
        for i in xrange(3):
            print cls.T_BUFFER[i].in_use[:64]

    def copy_struct(self, other=None, use_buf=False):
        """
            see iTensor_CopyStruct in f90
        """
        QSp= [q.copy() for q in self.QSp]
        totQN = self.totQN.copy()
        other = iTensor(rank=self.rank, QSp=QSp, totQN=totQN, use_buf=use_buf)
        return other

    def copy_struct_new(self, other=None, use_buf=False):
        """
            see iTensor_CopyStruct in f90
        """
        QSp= [q.copy() for q in self.QSp]
        totQN = self.totQN.copy()
        other = iTensor(rank=self.rank, QSp=QSp,totQN=totQN, 
                index_data=False, has_data=False, use_buf=use_buf)
        return other

    def shallow_copy(self):
        return self

    def copy(self, use_buf=False):
        """
            deepcopy
        """
        
        other = self.copy_struct(use_buf=use_buf)
        totDim=self.totDim
        #for large array(size>10^4) slicing is faster than copy
        
        other.data[:totDim]=self.data[:totDim]
        #print "cccc", 
        #try:
        #    other.data[:totDim]=self.data[:totDim]
        #except:
        #    msg = "%d\n%d\n%s"%(totDim, other.totDim, other)
        #    raise Exception(msg)
        return other

    def copy_new(self, buffer=None, use_buf=False):
        """
        see Tensor_ShallowCopy in f90
        """
        
        other = self.copy_struct(use_buf=use_buf)
        other.set_data_entrance(order="F")
        totDim = self.totDim
        if use_buf:   #use internal T_BUFFER; else use external buffer or no buffer
            buffer = other.buffer_assign(data_size=self.totDim)
        other.data = np.ndarray(self.totDim, buffer=buffer, dtype=self.dtype, order="C")
        other.data[:totDim]=self.data[:totDim]
        return other

    def __repr__(self, keys=None, fewer=0):
        """
        """
        #keys=["rank", "nidx","idx_dim", "Dims","totDim","data"]
        rank = self.rank
        if rank == 0:
            rank = 1
        if keys is None:
            keys=["rank", "type_name", "ndiv", "buf_ref","nidx","totQN","QNs", "Dims","totDim", "Block_idx", "Addr_idx", "data"]

        str0="----Begin iTensor----------------------------------------------\n"

        str1=""
        str2= '----End iTensor----------------------------------------------\n'
        for k in keys:
            if k == 'data':
                if fewer:
                    temp = "...."
                else:
                    temp = "\n"
                    for i in xrange(self.nidx):
                        start = self.Block_idx[0, i]
                        d = self.Block_idx[1, i]
                        qns= self.Addr_idx[:, i]
                        temp += str(qns) + ": "
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
                temp = [str(self.QSp[i].QNs[0:self.QSp[i].nQN]) for i in xrange(rank)]
                #temp  = repr(self.__dict__[k][:rank]) 
                temp = str(temp)
            elif k  == "Dims" :
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

        res=str0+ str1 +str2
        return res
    
    def show_data(self): 
        print self.__repr__(keys=['data'], fewer=0)
        
    @staticmethod
    def unit_tensor(rank, QSp, totQN=None):
        """
        
        """
        totQN = QSp[0].QnClass.qn_id()
        
        if rank%2 != 0:
            raise ValueError("rank shold be even, rank=%s"%(rank, ))
        t = iTensor(rank, QSp, totQN) 
        pTot = 1
        #use t.rank/2 would yeild a float 2.0
        rank1 =t.rank//2
        
        for i in xrange(rank1):
            pTot = pTot*t.QSp[i].nQN
        #这里rank1只有总rank的一半，pTot! = self.idx_dim

        t.data[:] = 0.0    
        pos= np.ndarray(rank1, "int")
        for i in xrange(pTot):
            p = i+i*pTot  #this is diagonal line
            idx=t.idx[p]
            p = t.Block_idx[0,idx]
            pos[0:rank1] = t.Addr_idx[0:rank1, idx]
            d = 1
            for j in xrange(rank1):
                d = d*t.QSp[j].Dims[pos[j]]

            t.data[p:p+d**2] = np.identity(d,dtype=t.dtype).ravel()
        return t
    
    @staticmethod
    def identity(qsp):
        return iTensor.unit_tensor(2, qsp)

    def get_position(self,iQN):
        """
            iQN: 量子数编号的组合
            量子数的编号也是用的 前小后大的存储顺序
            see iTensor_GetPosition in f90
            Get Position in T, get pos from iQN
        """
        if len(iQN)==0: #this is for rank 0 tensor
            return 0
        Dims=np.array([self.QSp[i].nQN for i in xrange(self.rank)],"int")
        #attention: this Dims is not self.Dims
        p=common_util.matrix_get_position(iQN,Dims)
        return p

    def get_position_rev(self, index_linear):
        """
            see iTensor_GetPosition_rev
            Get Position in T, get iQN from pos
        
        """
        Dims=np.array([self.QSp[i].nQN for i in xrange(self.rank)],"int")        
        #pos=np.empty(self.rank,"int")
        pos=common_util.matrix_get_position_rev(index_linear,Dims)
        return pos

    def set_element(self, qDims, iDims, X):
        """
            
            see iTensor_SetElement in f90
            params: 
                qDims: block coordinate, len(qDims)=self.rank
                    for each i, qDims[i] in xrange(self.Qsp[i].nQN),  对于Z2 symm即[0, 1]
                    并且这些量子数fuse到一起后 = self.totQN
                iDims：data coodinate in the block,  len(iDims)=self.rank
                    for each i, iDims[i] in xrange(self.Qsp[i].Dims), 对于Z2 symm 也是[0, 1] for any i
            #整个的思路是，先找到block entry，再找data position
        """
        pq = 0; pi=0
        for i in xrange(self.rank-1, 0 , -1):
            #print "iiiii i", i
            pq = (pq+qDims[i])*self.QSp[i-1].nQN
            pi = (pi+iDims[i])*self.QSp[i-1].Dims[qDims[i-1]]
            #print self.QSp[i-1].Dims[qDims[i-1]]
        
        pq = pq+qDims[0]
        pi = pi+iDims[0]
        temp = self.idx[pq]
        if temp >= self.idx_dim:
            raise ValueError("error, qdims values not permited")
        pidx = self.Block_idx[0, temp]
        block_len = self.Block_idx[1, temp]
        if pi>= block_len:
            raise ValueError("error, idims values not permited, iDims=%s"%iDims)
        #print "ppppqqqq  ", pq, pi,'idx', temp, 'pidx', pidx, pidx + pi 
        self.data[pidx+pi] = X

    def get_element(self, qDims, iDims):
        """
        status_1_verified
        see iTensor_GetElement in f90
        """
        pq = 0; pi=0
        for i in xrange(self.rank - 1, 0, -1):
            pq = (pq+qDims[i])*self.QSp[i-1].nQN
            pi = (pi+iDims[i])*self.QSp[i-1].Dims[qDims[i-1]]
        
        pq = pq+qDims[0]
        pi = pi+iDims[0]  
        pidx = self.Block_idx[0,self.idx[pq]]
        temp = self.idx[pq]
        if temp >= self.idx_dim:
            print "error, qDims values not permited"
            return 
        pidx = self.Block_idx[0, temp]
        block_len = self.Block_idx[1, temp]
        if pi>= block_len:
            print "error,  iDims values not permited"
            return 
        
        X = self.data[pidx+pi]
        return X

    def is_same_shape(self, other):
        """
        status_1_verified
        see iTensor_SameShape in f90
        """

        S = self.rank == other.rank
        if  not S:
            reason = -9
            return reason
        
        rank = self.rank
        for i in xrange( rank):
            S = self.QSp[i] == other.QSp[i]
            S= self.QSp[i]
            if  not S:
                print  'error, QSp[i]',i
                print self.QSp[i]
                print other.QSp[i]
                reason = -1
                return reason
        
        S = self.nidx  ==  other.nidx
        if  not S:
            reason = -2
            return reason
        
        S = self.totDim  ==  other.totDim
        if  not S:
            reason = -3
            return reason
        
        S = self.idx_dim  ==  other.idx_dim
        if  not S:
            reason = -4
            return reason
        
        for i in xrange( self.idx_dim):
            S = self.idx[i] == other.idx[i]
            if  not S:
                reason = -5
                return reason
        
        for i in xrange( self.nidx):
            S = self.Block_idx[0,i] == other.Block_idx[0,i]
            if  not S:
                reason = -6
                return reason
            
            S = self.Block_idx[1,i] == other.Block_idx[1,i]
            if  not S:
                reason = -7

                return reason
        reason = 1
        return reason
    
    class ShapeError(Exception):
        """
        only a tiny test of exception
        """
        pass

    def is_same_shape_new(self, other):
        """
        there is bug in it
        see iTensor_SameShape in f90
        """
        temp = ['rank', 'QSp', "nidx", "totDim", "idx_dim","idx", "Block_idx"]
        for k in temp:
            a = self.__getattribute__(k)
            b = other.__getattribute__(k)
            if isinstance(a, np.ndarray):
                if k == "idx":
                    end = self.idx_dim
                    res= np.all(a[:end]==b[:end])
                elif k == "Block_idx":
                    end == self.nidx 
                    #end == 1 
                    res= np.all(a[:2, :end]==b[:2, :end])
            else: 
                res = a == b 

            if not res:
                return False, k

        return True, None

    def index_merge_1(self, ind):
        """
            in numpy reshape serves as merge and split of indices of a dense array, I take use of that.
        """
        n1, n2 = min(ind), max(ind)
        ind1 = xrange(self.rank)  #indices not changed
        for i in ind: ind1.remove(i)  #indices to be merged
        ind2 = ind
        in1 = [i for i in ind1 if i <min(ind) ] + [i-1 for i in ind1 if i> max(ind)]
        in2 = min(ind)
        print ind1, ind2, in1, in2
        
        qsp_  = np.array([q.copy() for q in self.QSp], np.object)
        qsp = np.ndarray(self.rank-1, np.object)
        qsp[in1] = qsp_[ind1]
        qsp[in2] = qsp_[ind2[0]].add(qsp_[ind2[1]])
        qsp = qsp.tolist()
        totqn = self.totQN.copy()
        res= iTensor(self.rank-1, qsp, totqn)
        
        for i in xrange(res.nidx):
            iq = res.Addr_idx[:, i]
            iq1 = res.Addr_idx[in1, i]
            iq2 = res.Addr_idx[in2, i]
            sha = [res.QSp[x].Dims[iq[x]]  for x in xrange(res.rank)]
            p2 = res.Block_idx[0, i]
            size2 = res.Block_idx[1, i]
            stack = []
            for j in xrange(self.nidx):
                iqn = self.Addr_idx[:, j]
                
                if np.all(iqn[ind1]==iq[in1]):
                    shape2 = [self.QSp[x].Dims[iqn[x]]  for x in ind2] 
                    sha_ = [x for x in sha]
                    sha_[in2] = np.prod(shape2)
                    p = self.Block_idx[0, j]
                    size = self.Block_idx[1, j]
                    temp = self.data[p:p + size].reshape(sha_) 
                    stack.append(temp)
                    print "ii", i, "jj", j, 
                    print iqn, iq, sha, sha_
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
        
        for i in xrange(merger.nidx):
            p = merger.Block_idx[0, i]
            l = merger.Block_idx[1, i]
            iqn = merger.Addr_idx[:, i]
            shape = [merger.QSp[x].Dims[iqn[x]] for x in xrange(merger.rank)]
            data = merger.data[p:p + l].reshape(shape, order="C")
            print "iii", shape
            if 1:
                for i1 in xrange(shape[0]):
                    for i2 in xrange(shape[1]):
                        #i3 = i1 + shape[1]*i2
                        i3 = i2 + shape[1]*i1
                        data[i1, i2, i3] = 1.0 
        
        print merger
        res, nothing= self.contract(merger, xrange(self.rank), ind + [self.rank + 100])
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
        
        for j in xrange(self.nidx): 
            
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
            for i in xrange(self.nidx): 
                p = self.Block_idx[0, i]
                l = self.Block_idx[1, i]
                data = self.data[p: p + l]
                print data
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
                print pos
                #res.data[pos]
        if 1: 
            
            for i in xrange(res.nidx): 
                qn_id = res.Addr_idx[:, i]
                block_start = res.Block_idx[0, i]
                pointer = block_start
                
                for j in xrange(self.nidx): 
                    
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
        
    def index_split(self, ind, new_qsp_list):
        pass
    
    def permutation(self, P, buffer=None, use_buf=False):
        """
            permutation 是一个re-index的操作, 这一操作保持指标集整体不变，而局部置换
            todo: 
                I may write a inplace permutation, so needless QSp.copy()
                ind_labels  also need permuted
            see Tensor_Permutation in f90
            for permutation of iTensor, first we permute the order of QNs, next we permute each block which is a nTensor
            注意这里是Fortran order
            permutatiion won't change total quantum number
                
        """
        #print "in original"
        rank = self.rank
        #permute QSp, QNs
        #first is more efficient; second is more robust
        if 0:       #copy is more robust, while no copy is faster
            QSp=[self.QSp[P[i]].copy() for i in xrange(rank)]
            totQN = self.totQN.copy()
        else:
            QSp=[self.QSp[P[i]] for i in xrange(rank)]
            totQN = self.totQN

        Tp=iTensor(rank, QSp, totQN, buffer=buffer, use_buf=use_buf)
        
        #Tp.data[:] = 0.0

        pos=np.empty(self.rank, "int")
        Dims=np.empty(self.rank, "int")
        
        #permute each block
        #找到原来的block的位置与新的位置间的转换pidx<-->qidx
        #warnings.warn("using array_permutation_np")
        for n  in xrange(self.nidx):
            pidx = self.Block_idx[0,n]
            totDim = self.Block_idx[1,n]
            pos[0] = 0  # for rank=0
            pos[0:rank] = self.Addr_idx[0:rank,n]
            np1 = 0
            for i  in xrange(rank-1, 0, -1):
                np1 = (pos[P[i]]+np1)*Tp.QSp[i-1].nQN
                Dims[i] = self.QSp[i].Dims[pos[i]]
            i = 0
            np1 = np1+pos[P[i]]
            Dims[i] = self.QSp[i].Dims[pos[i]]
            qidx = Tp.Block_idx[0,Tp.idx[np1]]            
            
            temp= self.data[pidx:pidx+totDim]   #.copy()
            
            #attention_may_be_not_efficient  可以改成inplace 
            #Tp.data[qidx:qidx+totDim]=array_permutation.array_permutation_np(temp,rank,Dims,P)
            
            data = Tp.data[qidx:qidx+totDim]
            array_permutation.array_permutation_inplace(temp, self.rank, Dims, P, data)
            #Tp.data[qidx:qidx+totDim]=array_permutation.array_permutation(temp,rank,Dims,P)

        return Tp

    def contract_core(self, T2, div, preserve_qsp=False, data=None, use_buf=False):
        """
        see iTensor_Contraction2 in f90
        把T1，和T2的非零block 如果量子数组合相等则收缩
        locals:
            div: num. of legs to be contracted for each tensor
            buffer: use buffer to save data of T3
        """
        #print "ccc"

        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        tQN = self.totQN+T2.totQN
        shift = rank1-div

        #copy is needless conceptially
        if 0:       #copy is more robust, while no copy is faster
            QSp = [self.QSp[i].copy() for i in xrange(shift)]
            QSp.extend([T2.QSp[i].copy() for i in xrange(div, rank2)])
        else:
            QSp = self.QSp[:shift]
            QSp.extend(T2.QSp[div:rank2])


        if rank3==0:
            QSp = [self.QSp[0].null()]
        T3 = iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        T3.data[:]=0.0
        
        #T3.data=np.zeros(T3.totDim)  #this is really bad, need to re-allocate space for data
        
        nidx3 = 0
        alpha = 1.0; beta=1.0
        iQN1=np.empty(self.rank + 1,"int")
        iQN2=np.empty(T2.rank + 1,"int")        
        iQN3=np.empty(T3.rank + 1,"int") # +1 to avoid T3.rank=0

        for idx2 in xrange(T2.nidx):
            iQN2[0] = 0  #!for rank=0
            iQN2[0:rank2]=T2.Addr_idx[0:rank2,idx2]
            p2 = T2.Block_idx[0,idx2]
            
            Dim2 = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in xrange(div, rank2)])
            Dimc = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in xrange(div)])
            
            for idx1 in xrange(self.nidx):
                iQN1[0] = 1 #!for rank=0
                iQN1[0:rank1]=self.Addr_idx[0:rank1,idx1]
                p1 = self.Block_idx[0, idx1]
                #注意这里写得不适当，准确地，如果是
                #U1 symm. 的话应该是T1, T2相应的量子数的值正好差个符号, 而这里是用量子数的位置处理了, 并假定....写不清楚啊
                iseq = np.all(iQN1[shift:shift+div] == iQN2[0:div])
                if not iseq:
                    #如果量子数组合相等则收缩
                    continue
                Dim1 = np.prod([self.QSp[i].Dims[iQN1[i]] for i in xrange(shift)])
                
                iQN3[0] = 1 #!for rank=1
                iQN3[0:shift] = iQN1[0:shift]
                iQN3[shift:rank3] = iQN2[div:rank2]
                p3=T3.get_position(iQN3[:T3.rank])
                idx3 = T3.idx[p3]
                
                # one frequent error is the IndexError: index (9) out of xrange (0<=index<9) in dimension 1
                #dic = locals()
                #temp = ["idx3", "iQN1", "iQN2", "iQN3"]
                #print get_local(dic, temp)

                p3 = T3.Block_idx[0,idx3]
                
                data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                
                #T3.data[p3] = self.data[p1].dot(T2.data[p2])
                
                #data3=common_util.matrix_multiply(data1, data2, alpha, beta)   #alpha beta here have no effect
                
                #attention_may_be_not_efficient  这里学要为data3分配内存，能否直接在T3.data上操作？
                data3=iTensor.mul_temp(data1, data2, alpha, beta)

                #T3.data[p3:p3+Dim1*Dim2]=data3.ravel()[:]
                T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
                
                #t3data = T3.data[p3:p3+Dim1*Dim2] 
                #data3 = data3.ravel("F")
                #t3data  += data3   #attention_here  fortran order
                #t3data = numexpr.evaluate("t3data + data3")
                #print "iii",p3, data3, idx2, idx1,iQN1[:4],iQN2[:4],iQN3[:4] #,T3.data[:5].round(4)
                #print "iii",p3, data3, idx2, idx1,iQN1[:4],iQN2[:4],iQN3[:4] #,T3.data[:5].round(4)

        return T3

    def contract_core_new(self, T2, div, data=None, use_buf=False):
        """
        status_1_uncheck
        see iTensor_Contraction2 in f90
        把T1，和T2的非零block 如果量子数组合相等则收缩
        locals:
            div: num. of legs to be contracted for each tensor
            buffer: use buffer to save data of T3
        """
        
        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        tQN = self.totQN+T2.totQN
        shift = rank1-div

        #QSp = [self.QSp[i].copy() for i in xrange(shift)]
        #QSp.extend([T2.QSp[i].copy() for i in xrange(div, rank2)])
        #copy is needless conceptially
        QSp = self.QSp[:shift]
        QSp.extend(T2.QSp[div:rank2])
        
        T3= iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        T3.data[:]=0.0
        
        #print "RRRR T3 ", T3.rank, "use_buff", use_buf, T3.use_buf
        #T3.data=np.zeros(T3.totDim)  #this is really bad, need to re-allocate space for data
        #T3.data[:T3.totDim] = 0.0
        
        nidx3 = 0
        alpha = 1.0; beta=1.0
        iQN1=np.empty(self.rank,"int")
        iQN2=np.empty(self.rank,"int")        
        iQN3=np.empty(self.rank,"int")                
        for idx2 in xrange(T2.nidx):
            iQN2[0] = 0  #!for rank=0
            iQN2[0:rank2]=T2.Addr_idx[0:rank2,idx2]
            p2 = T2.Block_idx[0,idx2]
            
            Dim2 = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in xrange(div, rank2)])
            Dimc = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in xrange(div)])
            
            for idx1 in xrange(self.nidx):
                iQN1[0] = 1 #!for rank=0
                iQN1[0:rank1]=self.Addr_idx[0:rank1,idx1]
                p1 = self.Block_idx[0,idx1]
                iseq = np.all(iQN1[shift:shift+div] == iQN2[0:div])
                if not iseq:
                    #如果量子数组合相等则收缩
                    continue
                Dim1 = np.prod([self.QSp[i].Dims[iQN1[i]] for i in xrange(shift)])
                
                iQN3[0] = 1 #!for rank=1
                iQN3[0:shift] = iQN1[0:shift]
                iQN3[shift:rank3] = iQN2[div:rank2]
                p3=T3.get_position(iQN3[:T3.rank])
                idx3 = T3.idx[p3]
                p3 = T3.Block_idx[0,idx3]

                data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                
                
                #data3=T3.data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order="F")
                #below doesn't work as data3.base is NOT longer T3.data 
                #data3 = np.ndarray(Dim1*Dim2, order="F")
                data3 = np.ndarray((Dim1, Dim2), order="F")
                common_util.matrix_multiply_inplace(data1, data2, data3, alpha=1.0, beta=1.0) 
                
                
                #common_util.matrix_multiply_inplace(data1, data2, T3.data[p3:p3+Dim1*Dim2], alpha=1.0, beta=1.0) 
                



                
                #T3.data[p3] = self.data[p1].dot(T2.data[p2])
                
                #data3=common_util.matrix_multiply(data1, data2, alpha, beta)   #alpha beta here have no effect
                
                #attention_may_be_not_efficient  这里学要为data3分配内存，能否直接在T3.data上操作？
                #data3=iTensor.mul_temp(data1, data2, alpha, beta)

                #T3.data[p3:p3+Dim1*Dim2]=data3.ravel()[:]
                #T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
                
                #t3data = T3.data[p3:p3+Dim1*Dim2] 
                #data3=common_util.matrix_multiply(data1, data2, alpha, beta)   #alpha beta here have no effect
                #call dgemm('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)

                #data3 = data3.ravel("F")
                #t3data  += data3   #attention_here  fortran order
                #t3data = numexpr.evaluate("t3data + data3")

                #print "iii",p3, data3, idx2, idx1,iQN1[:4],iQN2[:4],iQN3[:4],T3.data.round(4)

        return T3

    @staticmethod
    def mul_temp(data1, data2, alpha, beta):
        """ only for profiling"""
        return common_util.matrix_multiply(data1, data2, alpha, beta) 
    
    def prepare_leg(self,T2, V1, V2, info=0):
        if not isinstance(V1, np.ndarray):
            V1=np.array(V1)
        if not isinstance(V2, np.ndarray):
            V2=np.array(V2)
        
        if self.rank<V1.size:
            V1=V1[:self.rank]
        if T2.rank<V2.size:
            V2=V2[:T2.rank]

        V_1n2 = np.intersect1d(V1, V2,True)
        V3=np.setxor1d(V1,V2,True)

        if info>0:
            print "contracted legs:\n\t", V_1n2

        k=0
        l=0
        
        Vp1=np.empty(self.rank,'i')
        Vp2=np.empty(T2.rank,'i')

        #below calculate Vp1, Vp2
        for i  in xrange(self.rank):
            if V1[i] not in V_1n2:
                Vp1[k] = i   
                # 这里之所以用k，而非直接用Vp1[k]是因为下面k的值要继续，对于l 也一样
                # p意指position，记录的T1中没有收缩的指标(leg)的编号
                V3[l] = V1[i]   #T3 来自T1 V1的外腿
                l = l+1
                k = k+1

        #这一段程序做了两件事：1.验证内线上维数相等，2.计算了totDim3
        #V_1n2=list(s_1n2)
        for i in xrange(len(V_1n2)):
            pos = np.where(V1==V_1n2[i])[0]
            pos2 = np.where(V2==V_1n2[i])[0]
            Vp1[k] = pos
             
            k = k+1  #注意这里k接着上面的值了
            if self.Dims[pos] != T2.Dims[pos2]:
                #raise Exception("Error, size of contract tensor does not match, V1=%s, V2=%s\nself=%s\nT2=%s"%(V1, V2, self, T2))
                #msg ="""error, dim of index to be contracted not equal: %s, %s, 
                #\n V1=%s, V2=%s, self.Dims=%s, T2.Dims=%s"""%(
                #        self.type_name, T2.type_name, V1, V2, self.Dims, T2.Dims)
                #msg ="""error, dim of index to be contracted not equal: {pos}""".format(locals())
                msg ="""error, dim of index to be contracted not equal: 
                    {0.type_name}, {1.type_name}
                    {V1}, {V2}, {V_1n2}
                    {0.Dims}, {1.Dims}
                    
                    """.format(self, T2, V1=V1, V2=V2, V_1n2=V_1n2)
                raise Exception(msg)
        
        k=0
        for i  in xrange(len(V_1n2)):
            #pos=V2.index(V_1n2[i])
            pos=np.where(V2==V_1n2[i])[0]
            Vp2[k] = pos
            k = k+1
       
        for i  in xrange(T2.rank):
            if V2[i] not in V_1n2:
                Vp2[k] = i
                #try:
                #    Vp2[k] = i
                #except:
                #    print Vp2, k, "len", len(Vp2)
                #    raise
                V3[l] = V2[i]  #T3 V3 来自T2 V2的外腿
                k = k+1
                l = l+1
        return V_1n2, Vp1, Vp2, V3 

    def contract(self,T2, V1, V2, order_final=None, out_Vc=False, data=None, use_buf=False, preserve_qsp=False, info=0):
        """
            see Tensor_Contraction2 in f90
            V1,2,3 are arrays (maps) 张量指标 —> 自然数。用自然数来标记所有张量的指标
            在V1,V2中可能有相同的元素，存在s_1n2中，
            
            T1, T2 contract yielding T3
            todo:
                order_final is not implememnted, when contracting the orignal order is changed, finally should be restored
            Args:
                Vi: of type np.ndarray(,"int")

                Vp1,先记录了T1的外腿，后记录内腿指标； Vp2先记录了内腿，后记录了外腿指标
        """
        if info>0:
            #\t{V1[0:rank1]}\t#{V2}[0:rank2]
            msg = """\ncontracting: {0.type_name}\t{1.type_name}
                rank: \t{0.rank}\t{1.rank}
                label: \t{V1}\t{V2} 
                dims: {0.Dims}\t{1.Dims}""".format(self, T2, V1=V1[:self.rank], V2=V2[:T2.rank])
            print msg
            if 0:
                print "\n\nstart contracting:", self.type_name, T2.type_name
                print "rank1, rank2, V1, V2:\n\t",self.rank, T2.rank,  V1[:self.rank], V2[:T2.rank]
                #print "T1.data: ", self.data[:4].round(5), "...", self.data[-4:].round(5)
                #print "T2.data: ", T2.data[:4].round(5), "...", T2.data[-4:].round(5)
                print self.data.round(3)
                print T2.data.round(3)
        if 0:
            #to impletement in the future
            self.ind_labels = dict(zip(range(self.rank), V1[:self.rank]))
            T2.ind_labels = dict(zip(range(T2.rank), V1[:T2.rank]))
        
        V_1n2, Vp1,Vp2,V3=self.prepare_leg(T2,V1,V2)        
        T1=self
        
        nT1=T1.permutation(Vp1, use_buf=use_buf)    #把T1 按照 Vp1 重排
        nT2=T2.permutation(Vp2, use_buf=use_buf)
        
        try:
            T3 = nT1.contract_core(nT2, V_1n2.size, data=data, use_buf=use_buf)
        except IndexError as err:
            msg = """ contract_core failed
            V_1n2:%s
            V1:%s\t  V2:%s
            Vp1:%s\t Vp2:%s
            """%(V_1n2, V1[:self.rank], V2[:T2.rank], Vp1[:self.rank], Vp2[:T2.rank])
            msg += """dims: {0.Dims}\t{1.Dims}""".format(self, T2, V1=V1[:self.rank], V2=V2[:T2.rank])
            msg +=  "\n%s\t %s"%(self.__repr__(keys=['QNs', 'Dims'], fewer=True), 
                    T2.__repr__(keys=['QNs', 'Dims'],fewer=True))
            #raise Exception(msg)
            print msg
            raise err
        

        T3.type_name = str(self.type_name) + "-" + str(T2.type_name)
        T3.ind_labels= V3
        
        if info>0: 
            #print "T3.data:", T3.data[range(4) + xrange(-4, 0)].round(5)        
            #print "T3.data:", T3.data[range(4)].round(5), "...", T3.data[range(-4, 0)].round(5)
            msg = """\n\tpermuted order: \t{Vp1}\t{Vp2}
                    \tyielding:  rank3={rank3}\tV3={V3}""".format( rank3=T3.rank, V3=V3, Vp1=Vp1[:self.rank], Vp2=Vp2[:T2.rank])
            print msg
        
        if out_Vc:
            Vc= V_1n2
            return  T3, V3, Vc 
        else:
            return T3, V3

    def trace(self):
        """
            status_1_uncheck
            see Tensor_Trace in f90
            从此函数看出， 要确定iTensor的对角线首先要确定它的对角块，即出入脚量子数相等;
            然后把对角块展开成2D matrix 求其trace
        """
        X = 0.0
        rank = self.rank
        rank2 = rank//2
        iQN=np.empty(self.rank,"int")
        for idx  in xrange(self.nidx):
            iQN[0:rank]=self.Addr_idx[0:rank,idx]
            #IsDiag = iArrayEq(rank2, iQN[0], iQN[rank2+1])
            IsDiag=np.all(iQN[:rank2]==iQN[rank2:])
            #judge wether on diagnal line
            if  not  IsDiag:  continue
            d = 1
            for i  in xrange(rank2):
                d = d*self.QSp[i].Dims[iQN[i]]
                a=self.Block_idx[0,idx]
            #X = X+common_util.matrix_trace(d, self.data[self.Block_idx[0,idx]])
            X+= self.data[a:a+d*d].reshape(d,d).trace()
        return X

    def partial_trace_bac(self, i):
        """
            this is a temporary implementation during calc central charge
            a full and fast implementation should not use tensor contraction, but direct sum uing index of sparse tensor
        """
        rank = self.rank
        assert rank%2 == 0 
        r_half = rank//2

        qsp = self.QSp[i].copy_many(2, reverse=[1])
        #I = iTensor(2, qsp.copy_many(2, reverse=[1]), qsp.QnClass.qn_id())
        I = iTensor.unit_tensor(2, qsp)

        v1 = xrange(rank)
        v2 = [i, r_half + i]
        res, leg= self.contract(I, v1, v2)
        #print 'leg', leg
        
        return res

    def partial_trace(self, site_list):
        """
            this is a temporary implementation during calc central charge
            a full and fast implementation should not use tensor contraction, but direct sum uing index of sparse tensor
        """
        rank = self.rank
        assert rank%2 == 0 
        rank_half= rank//2
        
        nsite = len(site_list)

        qsp = self.QSp[0].copy_many(2*nsite, reverse=range(nsite, 2*nsite))
        #I = iTensor(2, qsp.copy_many(2, reverse=[1]), qsp.QnClass.qn_id())
        I = iTensor.unit_tensor(2*nsite, qsp)

        v1 = xrange(rank)
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
        for i in xrange(d,self.rank):
            ord[i-d] = i
        for i  in xrange(d):
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
        for i in xrange(d,self.rank):
            ord[i-d] = i
        for i  in xrange(d):
            ord[i+self.rank-d] = i
        temp=self.permutation(ord, buffer=buffer, use_buf=use_buf)
        
        temp.reverse_qsp(inplace=False)
        return temp

    
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
        for n in xrange(self.rank):
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
        for idx  in xrange(self.nidx):
            #这一loop先确定block坐标, block的大小
            block_pos_in_data= self.Block_idx[0,idx]
            qn_comb_pos = self.Block_idx[2,idx]            
            #qn_comb_pos is the position of QN, 线性坐标与多维坐标点转换
            #pos=self.get_position_rev(qn_comb_pos)
            #pos is coordinate of a set of quantum number 
            qn_coord=self.get_position_rev(qn_comb_pos)
            block_size=[self.QSp[i].Dims[qn_coord[i]] for i in xrange(rank)]
                #block 每个指标的维数,与上面的Dims意思不同
            
            #attention_omitted_something
            Ti=nTensor(rank, block_size, shallow=False)
            print "block_pos",idx,"block_size",block_size,"idx",idx,"qn_coord",qn_coord

            for dat_pos_in_blk  in xrange(self.Block_idx[1,idx]):
                print "\t dat_pos_in_blk ", dat_pos_in_blk
                #这一loop确定block内的坐标
                data_coord_in_block=Ti.get_position_rev(dat_pos_in_blk)
                #data_coord_in_block is coordinate in the block <===  p is linear coordinate
                
                data_coord=[]
                for i  in xrange(rank):
                    d = 0
                    for j  in xrange(qn_coord[i]):
                        d = d+self.QSp[i].Dims[j]
                    
                    data_coord.append(int(d + data_coord_in_block[i]))
                
                    
                
                data_pos=Tp.get_position(data_coord)
                print "\t\t data_coord",data_coord,"data_pos",data_pos,
                
                Tp.data[data_pos] = self.data[block_pos_in_data+dat_pos_in_blk]
                print round(self.data[block_pos_in_data+dat_pos_in_blk],10)
        return Tp

    def get_nposition(self):
        """
        see iTensor_GetNPosition
        status_1
        """
        pass
    
    def direct_product_back_back(self, T2):
        """
        status_1_verified
        this function was originally defined under Tensor class, now is moved here
        see Direct_Product in f90

        """
        QSp=[QuantSpace() for i in xrange(self.rank)] 

        T1=self
        #attention: here the division may be incorrect
        #q: why /2?
        rank1 = T1.rank//2
        rank2 = T2.rank//2      
        rank3 = rank1+rank2
        for i in xrange(rank1):
            #q: attention: 注意这里还有下面从0开始计数,可能不对
            #attention 这里应该用deepcopy？
            QSp[i] = T1.QSp[i]
            QSp[i+rank3] = T1.QSp[i+rank1]
        
        for i in xrange(rank2):
            QSp[i+rank1] = T2.QSp[i]
            QSp[i+rank1+rank3] = T2.QSp[i+rank2]
        
        totQN = T1.totQN+T2.totQN
        #print 'ttt', totQN.val, T1.totQN.val, T2.totQN.val
        #init_Tensor(T1.rank+T2.rank, QSp, totQN, T3)
        #q: attention:  here may be incoorect
        
        T3= iTensor(T1.rank+T2.rank, QSp, totQN)
        
        T3.data[:] = 0.0
        V1 = np.empty(self.rank, "int")
        V2 = np.empty(self.rank, "int")
        V3 = np.empty(self.rank, "int")
        for idx1 in xrange(T1.nidx):
            pidx1 = T1.Block_idx[0, idx1]
            V1[0:T1.rank] = T1.Addr_idx[0:T1.rank,idx1]
            nA=1; mA=1
            for i in xrange(rank1):

                nA = nA*T1.QSp[i].Dims[V1[i]]
                mA = mA*T1.QSp[i+rank1].Dims[V1[i+rank1]]
                V3[i] = V1[i]
                V3[i+rank3] = V1[i+rank1]

            for idx2 in xrange(T2.nidx):
                pidx2 = T2.Block_idx[0, idx2]
                V2[0:T2.rank] = T2.Addr_idx[0:T2.rank, idx2]
                nB=1; mB=1
                for i in xrange(rank2):
                    nB = nB*T2.QSp[i].Dims[V2[i]]
                    mB = mB*T2.QSp[i+rank2].Dims[V2[i+rank2]]
                    V3[i+rank1] = V2[i]
                    V3[i+rank1+rank3] = V2[i+rank2]
                #idx3 linear position of  qn combination V3
                idx3 = T3.get_position(V3[:rank3*2])   
                
                pidx3 = T3.Block_idx[0, T3.idx[idx3]]
                    #print 'selffff.idx', T1.idx, T2.idx
                #Matrix_DirectProduct[T1.data[pidx1], nA, mA, T2.data[pidx2], nB, mB, T3.data[pidx3]]
                data1 = T1.data[pidx1:pidx1 + nA*mA].reshape((nA, mA))
                data2 = T2.data[pidx2:pidx2 + nB*mB].reshape((nB, mB))                
                T3.data[pidx3:pidx3 + nA*nB*mA*mB] = common_util.matrix_direct_product(data1, data2).ravel()
        return T3
    
    def direct_product_bac(self, T2, order="F", use_buf=False):
        """
        this function was originally defined under Tensor class, now is moved here
        see Direct_Product in f90
        parameters:
            T1^{I1}_{J1}, T2^{I2}_{J2}
        returns:
            T3^{I1I2}_{J1J2}
            i.e. index(QSp) of T3 is ordered J1, J2, I1, I2 in sequel


        """
        QSp=[None for i in xrange(self.MaxRank)] 

        T1=self
        #attention: here the division may be incorrect
        #q: why /2?
        rank1 = T1.rank//2
        rank2 = T2.rank//2      
        rank3 = rank1+rank2
        for i in xrange(rank1):
            QSp[i] = T1.QSp[i].copy()
            QSp[i+rank3] = T1.QSp[i+rank1].copy()
        for i in xrange(rank2):
            QSp[i+rank1] = T2.QSp[i].copy()
            QSp[i+rank1+rank3] = T2.QSp[i+rank2].copy()
        
        #这里暗含了一个张量积标识fusion的过程
        totQN = T1.totQN+T2.totQN
        T3= iTensor(T1.rank+T2.rank, QSp, totQN, use_buf=use_buf)
        #print 'ttt', totQN.val, T1.totQN.val, T2.totQN.val
        
        T3.data[:] = 0.0
        V1 = np.empty(self.MaxRank, "int")
        V2 = np.empty(self.MaxRank, "int")
        V3 = np.empty(self.MaxRank, "int")
        for idx1 in xrange(T1.nidx):
            pidx1 = T1.Block_idx[0, idx1]
            V1[0:T1.rank] = T1.Addr_idx[0:T1.rank,idx1]
            nA=1; mA=1
            for i in xrange(rank1):

                nA = nA*T1.QSp[i].Dims[V1[i]]
                mA = mA*T1.QSp[i+rank1].Dims[V1[i+rank1]]
                V3[i] = V1[i]
                V3[i+rank3] = V1[i+rank1]

            for idx2 in xrange(T2.nidx):
                pidx2 = T2.Block_idx[0, idx2]
                V2[0:T2.rank] = T2.Addr_idx[0:T2.rank, idx2]
                nB=1; mB=1
                for i in xrange(rank2):
                    nB = nB*T2.QSp[i].Dims[V2[i]]
                    mB = mB*T2.QSp[i+rank2].Dims[V2[i+rank2]]
                    V3[i+rank1] = V2[i]
                    V3[i+rank1+rank3] = V2[i+rank2]
                #idx3 linear position of  qn combination V3
                idx3 = T3.get_position(V3[:rank3*2])   
                
                pidx3 = T3.Block_idx[0, T3.idx[idx3]]
                    #print 'selffff.idx', T1.idx, T2.idx
                #Matrix_DirectProduct[T1.data[pidx1], nA, mA, T2.data[pidx2], nB, mB, T3.data[pidx3]]
                #data1 = T1.data[pidx1:pidx1 + nA*mA].reshape((nA, mA), order="F")
                #data2 = T2.data[pidx2:pidx2 + nB*mB].reshape((nB, mB), order="F")                
                data1 = T1.data[pidx1:pidx1 + nA*mA].reshape((nA, mA))
                data2 = T2.data[pidx2:pidx2 + nB*mB].reshape((nB, mB))                
                T3.data[pidx3:pidx3 + nA*nB*mA*mB] = common_util.matrix_direct_product(data1, data2).ravel(order=order)

        return T3

    def direct_product(self, T2, order=None, use_buf=False):
        """
        this function was originally defined under Tensor class, now is moved here
        see Direct_Product in f90
        parameters:
            T1^{I1}_{J1}, T2^{I2}_{J2}
        returns:
            T3^{I1I2}_{J1J2}
            i.e. index(QSp) of T3 is ordered J1, J2, I1, I2 in sequel


        """
        

        T1=self
        #attention: here the division may be incorrect
        #q: why /2?
        rank1 = T1.rank//2
        rank2 = T2.rank//2      
        rank3 = rank1+rank2

        QSp=[0  for i in xrange(T1.rank + T2.rank)] 
        for i in xrange(rank1):
            QSp[i] = T1.QSp[i].copy()
            QSp[i+rank3] = T1.QSp[i+rank1].copy()
        for i in xrange(rank2):
            QSp[i+rank1] = T2.QSp[i].copy()
            QSp[i+rank1+rank3] = T2.QSp[i+rank2].copy()
        
        #这里暗含了一个张量积标识fusion的过程
        totQN = T1.totQN+T2.totQN
        T3= iTensor(T1.rank+T2.rank, QSp, totQN, use_buf=use_buf)
        #print 'ttt', totQN.val, T1.totQN.val, T2.totQN.val
        
        T3.data[:] = 0.0
        V1 = np.empty(self.MaxRank, "int")
        V2 = np.empty(self.MaxRank, "int")
        V3 = np.empty(self.MaxRank, "int")
        for idx1 in xrange(T1.nidx):
            pidx1 = T1.Block_idx[0, idx1]
            V1[0:T1.rank] = T1.Addr_idx[0:T1.rank,idx1]
            nA=1; mA=1
            for i in xrange(rank1):

                nA = nA*T1.QSp[i].Dims[V1[i]]
                mA = mA*T1.QSp[i+rank1].Dims[V1[i+rank1]]
                V3[i] = V1[i]
                V3[i+rank3] = V1[i+rank1]

            for idx2 in xrange(T2.nidx):
                pidx2 = T2.Block_idx[0, idx2]
                V2[0:T2.rank] = T2.Addr_idx[0:T2.rank, idx2]
                nB=1; mB=1
                for i in xrange(rank2):
                    nB = nB*T2.QSp[i].Dims[V2[i]]
                    mB = mB*T2.QSp[i+rank2].Dims[V2[i+rank2]]
                    V3[i+rank1] = V2[i]
                    V3[i+rank1+rank3] = V2[i+rank2]
                #idx3 linear position of  qn combination V3
                idx3 = T3.get_position(V3[:rank3*2])   
                
                pidx3 = T3.Block_idx[0, T3.idx[idx3]]
                    #print 'selffff.idx', T1.idx, T2.idx
                #Matrix_DirectProduct[T1.data[pidx1], nA, mA, T2.data[pidx2], nB, mB, T3.data[pidx3]]
                data1 = T1.data[pidx1:pidx1 + nA*mA].reshape((nA, mA), order="F")
                data2 = T2.data[pidx2:pidx2 + nB*mB].reshape((nB, mB), order="F")                
                T3.data[pidx3:pidx3 + nA*nB*mA*mB] = common_util.matrix_direct_product(data1, data2).ravel(order="F")

        return T3


    def tensor_prod(self, T2, order="F"):
        return self.direct_product(T2, order)

    def matrix_view(self, n=None, order='C', round=None):
        """
        map rank (n, rank-n) tensor to a 2-d mat. 
        this method is used for debug or getting a intuitive view of the tensor

        Returns: 2-d array res
        """
        if n == None:
            n = self.rank//2
        pos= 0
        rank = self.rank
        dim1 = np.prod(self.Dims[:n])
        dim2 = np.prod(self.Dims[n:self.rank])
        res= np.zeros(np.prod(self.Dims[:self.rank])).reshape((dim1, dim2))
        
        ddd = [self.QSp[i].nQN for i in xrange(self.rank)]
        pos1 = 0 
        pos2 = 0        
        A = np.prod( [self.QSp[i].nQN for i in xrange(n)])
        B = np.prod( [self.QSp[i].nQN for i in xrange(n, self.rank)])
        iQN = np.zeros(self.rank, dtype="int")
        daaa = np.zeros(B, "int")
        for x in xrange(A):
            pos2 = 0
            for y in xrange(B):
                pos1 = daaa[y]
                p = x*A  + y
                
                #print 'p', p, x, y
                idx = self.idx[p]
                if idx < self.nidx:
                    d = np.prod([self.QSp[i].Dims[iQN[i]]  for i in xrange(rank) ])
                    d1 = np.prod([self.QSp[i].Dims[iQN[i]] for i in xrange(n)])
                    d2 = np.prod([self.QSp[i].Dims[iQN[i]] for i in xrange(n,self.rank)])            
                    start = self.Block_idx[0, idx]
                    res[pos1:pos1+d1, pos2:pos2+d2] = self.data[start:start + d].reshape(d1, d2)
                    #print res
                db = np.prod([self.QSp[i].Dims[iQN[i]] for i in xrange(n, self.rank)])
                            
                pos2 += db 
                da = np.prod([self.QSp[i].Dims[iQN[i]] for i in xrange(n)])
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
                #kkk= common_util.matrix_get_position_rev(p, ddd)
                #da = np.prod([self.QSp[i].Dims[kkk[i]] for i in xrange(n)])
        if round is None:
            return res
        else:
            return res.round(round)  

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
        not tested
        """
        rank = self.rank
        T2 = iTensor(rank, QSp, self.totQN.copy())
        T2.data[:] = 0.0
        iQN = np.ndarray(rank, np.int)
        Dims1 = np.ndarray(rank, np.int)
        Dims2 = np.ndarray(rank, np.int)
        #transvers all dense blocks
        for idx  in xrange(self.nidx):
            p1 = self.Block_idx[0,idx] #-1
            iQN[0:rank] = self.Addr_idx[0:rank, idx]
            pidx2 = T2.get_position(iQN)
            idx2 = T2.idx[pidx2]
            p2 = T2.Block_idx[0,idx2] #-1
            # above find the start address p1 and p2 for corresonding blocks
            for i  in xrange(rank):
                Dims1[i] = self.QSp[i].Dims[iQN[i]]
                Dims2[i] = T2.QSp[i].Dims[iQN[i]]
            
            #within each dense block, transvers elements
            for ip1 in xrange(self.Block_idx[1,idx]):
                #ip1 is linear index, ip1->pos
                pos = common_util.matrix_get_position_rev(ip1, Dims1)
                ip2 = common_util.matrix_get_position(pos, Dims2)
                T2.data[p2+ip2] = self.data[p1+ip1]
        return T2

    def act_on_qsp(self, qsp): 
        pass
    
class iTensor_new(TensorBase):
    def __init__(self,rank,  QSp, totQN, shallow=None, use_buf=None):
        """
        status_1_verified
        self.Dims[i]:  在leg i 对应的空间维数
        
        """
        TensorBase.__init__(self,rank,[0]*rank)
        use_buf0=False
        shallow0=False
        self.ndiv = None

        self.use_buf=use_buf
        self.shallow=shallow
        if shallow!= None:
            shallow0=shallow
        if use_buf!=None:
            use_buf0=use_buf
        
        if use_buf0:
        #attention_omitted_something
            self.use_buf=True
            if self.buf_ref[0] ==-1:
                self.use_Buffer(rank,QSp,totQN)
                #attention_omitted_something

        if self.use_buf:
            pass
            #attention_omitted_something

        self.rank=rank
        self.QSp=list(QSp[0:rank])
        self.totQN= totQN

        #self.quantum_number = [self.QSp[0].QNs[i] for i in xrange(2)]  #这里仅考虑Z2 symm，nQn=2
        #self.qunt_num_comb = np.ndarray((self.rank, 2))

        pTot =1
        self.Dims=np.empty(rank,dtype=self.dtype)
        self.Dims[0] =1
        for i in xrange(rank):
            pTot *= QSp[i].nQN
            self.Dims[i]= QSp[i].totDim
        self.idx_dim=pTot

        if self.shallow:
            print r"#attention_omitted_something"

        #if (not self.allocated) or (#self.max_ind_size >pTot):
        #attention_omitted_something

        
        #self.idx=np.ndarray(self.idx_dim,"int")   #-1
        #self.idx=np.array([-1L]*self.idx_dim,"int")   #-1        
        #attention_this_may_be_wrong 在python中  -1对应着最后一个元素，而Fortran中什么也不对应, 所以改成下面的
        self.idx=np.array([int(self.idx_dim)]*self.idx_dim,"int")   #-1                
        
        nidx=0
        totDim=0
        
        if rank==0:
            self.Dims[0]=1
            self.totDim=1
            self.QSp[1]= QSp_null
            self.totQN = QN_idendity  #.copy()
        
        iQN=[0]*self.rank
        #iQN[i]用作leg i 上的量子数 计数
        tQN_r= self.totQN.copy()
        tQN_r.reverse()

        temp = [self.QSp[i].nQN for i in xrange(rank)]

        self.data= np.ndarray(temp, dtype="object")

        for p in xrange(pTot):
            tQN = self.QSp[0].QNs[iQN[0]]
            for i in xrange(1,rank):
                tQN = tQN+self.QSp[i].QNs[iQN[i]]
            
            if tQN==tQN_r:
                block_shape= [QSp[i].Dims[iQN[i]] for i in xrange(rank)]
                self.data[tuple(iQN)] = np.empty(block_shape, dtype=self.dtype)
            else:
                self.data[tuple(iQN)] = 0

            inc = 1
            i = 0
            #print "iii ", iQN,  tQN==tQN_r
            while inc==1 and i<rank:
                iQN[i] = iQN[i]+1
                if iQN[i]<QSp[i].nQN :
                    inc = 0
                else:
                    iQN[i] = 0
                    i = i+1
            
        #self.nidx = nidx
        #self.totDim = totDim


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
    def isometry():
        pass
    
    @staticmethod
    def V111_1():
        pass
    
    #unit_tensor = iTensor.unit_tensor
    def unit_tensor_del_it_or_not(self, rank):
        QSp = self.qsp_base.copy_many(rank, reverse=range(rank//2, rank))
        totQN = self.qn_identity.copy()
        return iTensor.unit_tensor(rank, QSp, totQN)

    @staticmethod
    def pauli_mat_2site(symmetry, which=None):
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
                    print 'in tensor_py.py changeed '*30
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
                totqn.set_val(1)  #note this differs with sigma_z
                sigma_p1 = iTensor(rank, qsp, totqn)
                #sigma_p1 := sp×I
                sigma_p1.data[:] = 0.0
                sigma_p2 = sigma_p1.copy()
                
                # [-1,-1]->[1,-1]
                qDims[0:2] = [2,0]   #|ud><dd|
                iDims[0:2] = [0,1]
                sigma_p1.set_element(qDims, iDims, 1.0)
                
                # [-1,1]->[1,1]
                qDims[0:2] = [0,1]  #|uu><du|
                iDims[0:2] = [0,0]
                sigma_p1.set_element(qDims, iDims, 1.0)
                
                # [-1,-1]->[-1,1]
                qDims[0:2] = [2,0]  
                iDims[0:2] = [0,0]
                sigma_p2.set_element(qDims, iDims, 1.0)
                # [1,-1]->[1,1]        
                qDims[0:2] = [0,1]  
                iDims[0:2] = [1,0]
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
            
            s0, sx, sy, sz = sigma_0, sigma_x, sigma_y, sigma_z    

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
            totqn.set_val(2)  
            sigma_p = iTensor(rank, qsp, totqn)
            sigma_p.data[0] = 1.0
            sigma_m = sigma_p.conjugate(1)
            if 0:
                print sigma_z #.matrix_view()
                print sigma_p.matrix_view()
                print sigma_m.matrix_view()
            
            szz = sigma_z.tensor_prod(sigma_z)
            spm = sigma_p.tensor_prod(sigma_m)
            smp = sigma_m.tensor_prod(sigma_p)
        
        
        return vars()
    
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
            print "check using two methods to define h2 and h3"
            h2a = xx  + (-1.0)*yy
            h2b = 2.0*(pm + mp)
            print np.all(h2a.data==h2b.data)

            h3a = (-1.0)*y.direct_product(i).direct_product(y) + x.direct_product(i).direct_product(x)
            h3b = 2.0*(p.direct_product(i).direct_product(m) + m.direct_product(i).direct_product(p))
            print np.all(h3a.data==h3b.data)
        if self.symmetry == "Z2":
            pauli = self.pauli_mat()
            I_2, szz, spm, smp= pauli["sigma_0"], pauli["szz"], pauli["spm"], pauli["smp"]
            pm, mp = spm, smp
            sigma_z, sigma_p, sigma_m = pauli["sigma_z"], pauli["sigma_p"], pauli["sigma_m"]
            x, y=  pauli["sigma_x"], pauli["sigma_y"]
            zz, xx, yy=  pauli["szz"], pauli["sxx"], pauli["syy"]

            z, p, m = sigma_z, sigma_p, sigma_m
            i = I_2
            
            print "check using two methods to define h2 and h3"
            h2a = xx  + (-1.0)*yy
            h2b = 2.0*(pm + mp)
            print h2a.data
            print h2b.data

            h3a = (-1.0)*y.direct_product(i).direct_product(y) + x.direct_product(i).direct_product(x)
            h3b = 2.0*(p.direct_product(i).direct_product(m) + m.direct_product(i).direct_product(p))
            print np.all(h3a.data==h3b.data)
    
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
    
if 0:
    class Tensor(iTensor):
        """
        thie class is translated fortran, will be removed or replaced
        """
        def __init__(self,use_buf=False):
            """
            status:completed:uncheck
            see init_Tensor in f90
            """
            self.use_buf=use_buf
            self.buf_ref = np.array([-1, -1], "int")


        def nullify(self, use_buf=False):
            """
            see initialize_Tensor in f90
            """
            #self.allocated= False 
            self.shallow= False         
            self.use_buf = use_buf
            self.buf_ref[:]=-1  

            #self.max_data_size = 0
            #self.max_ind_size = 0
            self.rank = 0
            self.nidx = 0
            self.totDim = 0        

            self.idx = None
            self.Block_idx = None
            self.Addr_idx = None
            self.data = None

        def delete(self):
            """
            status_0
            deprecate it?

            """


        def shallow_copy_deprecate(self):
            """
            see Tensor_ShallowCopy in f90
            """
            
            other = self.copy_struct()
            other.idx[:] = self.idx[:]
            other.Block_idx[:, :] = self.Block_idx[:, :]
            other.Addr_idx[:, :] = self.Addr_idx[:, :]
            other.data[:] = self.data[:]
            return other

        @staticmethod 
        def find_leg(iLeg,  S):
            """
            status_1_uncheck
            see FindLeg in f90
            FindLeg is to find all the legs that are connected between two tensors
            应该是找到 腿 iLeg 在哪两个tensor里，存储在 Is[0:2]中
            returns:
                Is
            """

            #Is=2*[None]
            Is=[]
            for i in xrange(len(S)):
                #if InSet[iLeg, S[i]]>0:  
                if iLeg in S[i]:
                    #Is[k] = i
                    Is.append(i)
                    #k = k+1
            if len(Is)>2:
                print "error, see find_leg"
            return Is

        def random(self ):
            """
            status_0
            以后再处理随机数生成器
            """

        def str__(self):
            """
            print_Tensor in f90
            status_1
            """

        def unit_tensor_old(self):
            """
            status_1_deprecated,  
            this is precisely unit_itensor,  so just inhereite from iTensor
            see UnitTensor in f90
            """
            rank = self.rank
            rank1 = rank//2
            
            self.data[:] = 0.0
            for idx  in xrange(self.nidx):
                p = self.Block_idx[0,idx]
                iQN[0:rank] = self.Addr_idx[0:rank, idx]
                
                #AllEq = iArrayEq(rank1, iQN[1], iQN[rank1+1])
                AllEq = np.all(iQN[0:rank1]==iQN[rank1:rank])
                if  not AllEq:  
                    continue
                
                d=1
                for i  in xrange(rank1):
                    d = d*self.QSp[i].Dims[iQN[i]]
                #Unit_Matrix(self.data[p], d)
                self.data[p:p + d**2] = np.identity(d).ravel()

#   class Tensor_Link is moved to tensor_network.py


class test_iTensor(object):
    def __init__(self, symmetry, dim=None):
        self.qn_identity, self.qsp_base, self.qsp_null = init_System_QSp(symmetry)
        QSbase = self.qsp_base.copy
        totQN=self.qn_identity.copy
        if dim is None:
            QSbase = self.qsp_base.copy
        else:
            QSbase = self.qsp_base.__class__.max(dim).copy
        #qsp = self.qsp_base.max(dim)

        self.u=iTensor(4,[QSbase() for i in xrange(4)],totQN())
        self.w = iTensor(3, [QSbase() for i in xrange(3)], totQN())

        self.u22 = iTensor(2,[QSbase() for i in xrange(2)],totQN())
        
        self.u222=  iTensor(3,[QSbase() for i in xrange(3)],totQN())
        self.u2222=  iTensor(4,[QSbase() for i in xrange(4)],totQN())

    def dump(self):
        print "this func not work "
        import pickle
        out = open("/tmp/test", "wb")
        pickle.dump(self.u.QSp[0].QNs, out)
        out.close()
        exit()


    @classmethod
    def instance(cls, which="all"):
        u = cls.u.copy()
        w = cls.w.copy()
        u22 = cls.u22.copy()
        #t0 = cls.t0.copy()

        keys = locals().keys()

        if which == "all":
            return u, w, u22
        if which in keys:
            return locals()[which]

    @classmethod
    def init(cls):
        pass
        print cls.w
    
    @staticmethod
    def use_buf():
        u=iTensor(4,[QSp_base.copy() for i in xrange(4)],QN_idendity.copy(), use_buf=True)
        print u.data
        print iTensor.T_BUFFER[0].T[0]
        print u.data.base is iTensor.T_BUFFER[0].T[0]

        return u

    def test_copy():
        """  --- pass  """
        w1 = w.copy()
        #w1.QSp_
        w1.QSp[0] = QSbase().add(QSbase())
        print "www1 ", w1.QSp, '\n'*2
        print w.QSp
   
    @classmethod
    def is_same_shape(cls):
        """  ---pass """
        u = cls.u.copy()
        w = cls.w.copy()
        a = u.is_same_shape(u)
        print a
    
    def test_get_position_and_rev() :
        """ --- pass """
        p=w.get_position([0,0,0])
        print p
        for i in xrange(7):
            pos=w.get_position_rev(i)
            print pos
    
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
        print x
        print w244.__repr__(["data"])
    #test_get_element_and_set_element()
    def test_to_ntensor():
        """ ---pass   """
        qsp = [QSbase(), QSbase().add(QSbase()), QSbase().add(QSbase())]
        w= iTensor(3, qsp, totQN)
        
        w.set_element([0,0,0], [0, 0, 0], 3.)
        #w.set_element([1,0,1], [0, 1, 1], 4.)            
        w.set_element([0,1,1], [0, 1, 1], 5.)                        
        print w.__repr__(["data"])
        nw=w.to_nTensor()
        print nw
    #test_to_ntensor() 
    def test_trace():
        a=u.trace()
        print a
        pass
    #test_trace()
    
    @classmethod
    def contract_core(cls):
        """ --- not sure"""
        u = cls.u.copy()
        w = cls.w.copy()
        a=u.contract_core(u,2)
        print a

    @classmethod
    def contract_core_buff(cls):
        """using buffer --- pass"""
        u = cls.u.copy()
        w = cls.w.copy()
        a=u.contract_core(u,2,use_buf=True)
        print a.data
        print iTensor.T_BUFFER[0].T[0]
        print a.data.base is iTensor.T_BUFFER[0].T[0]


    @classmethod
    def contract_core1(cls):
        """ contract with a rank 0 tensor """
        t0 = cls.t0.copy()
        u = cls.u.copy()
        a=u.contract_core(t0,0)
        print a

    @classmethod
    def contract_core2(cls):
        """ contract to a rank 0 tensor """
        u = cls.u.copy()
        a=u.contract_core(u.copy(),4)
        print a

    @classmethod
    def contract_U_H2(cls):
        """   --- pass"""
        U = cls.u.copy()
        H2 = cls.u.copy()
        U.data[:] = [ 1.,  0.,  1.,  0.,  0.,  1.,  0.,  1.]
        H2.data[:] = [ 0.,  -1.,  -1.,  -0.,2 -0.,2 -1.,  -1.,  -2.]


        res= U.contract_core(H2, 2)
        #res= H2.contract_core(U,2)
        print U.data
        print H2.data

        print res.data.round(5)
    
    def contract(self):
        """    --- not sure """
        t1 = self.u2222.copy()
        t1.data[:] = np.arange(t1.totDim)
        t2 = self.u222.copy()
        t2.data[:] = np.arange(t2.totDim)
        #v1 = [1, 2, 5, 9]
        #v2 = [5, 1, 3, 6]
        v1 = ['a', 'b', 'c', 'd']
        v2 = ['e', 'c', 'a']
        t3 = t1.contract(t2, v1, v2,use_buf=True)
        print t3[0].rank
        print t3[0]
    
    @decorators.timer
    def contract_large_tensor(self, threads=10):
        import os
        os.environ["OMP_NUM_THREADS"] = str(threads)
        qn_identity, qsp_base, qsp_null = init_System_QSp("Z2")
        qsp_max , qsp_max2= qsp_base.__class__.max(16)
        QSbase = qsp_max.copy
        totQN = qn_identity.copy
        ranku = 6
        rankv = 6
        Vu = xrange(ranku)
        Vw = xrange(rankv)
        #Vw.reverse()

        for i in xrange(1):
            u=iTensor(ranku,[QSbase() for i in xrange(ranku)],totQN())
            w = iTensor(rankv, [QSbase() for i in xrange(rankv)], totQN())
            out=u.contract_core(w, 3)
            #out, leg=u.contract(w,Vu, Vw )
        print out.rank



    def decoration(self):
        if 1:
            pass
            common_util.contract_core_player_fort =\
                    common_util.contract_core_player_fort_parallel

        for i in xrange(5):
            
            rank = 4
            set_STATE_end_1(iter=i, record_at=0, stop_at=10000000, power_on=True) 
            qsp = self.qsp_base.copy_many(rank)
            totQN = self.qn_identity.copy()
            
            u = iTensor(rank=rank, QSp=qsp, totQN=totQN)
            u.data[:] = xrange(u.data.size)
            #v=u.contract(u, [0, 1, 2, 3], [8, 0, 3, 1])
            v=u.contract_core(u, 2)
            print v.data
            
            #set_STATE_end(i, 8)

    def test_permutation():
        """  ---not sure"""
        qsp = [QSbase(), QSbase().add(QSbase()), QSbase().add(QSbase())]
        w244 = iTensor(3, qsp, totQN)
        qDims= [1, 0, 1]
        w244.set_element(qDims, [0, 1, 1], 5.)
        w244.set_element(qDims, [0, 1, 0], 4.)
        w244.set_element(qDims, [0, 0, 1], 3.)            
        w244.set_element(qDims, [0, 0, 0], 2.)                        
        print w244.matrix_view(1).round(4)  
        
        W=w244.permutation([0,2,1])
        print W.matrix_view(1).round(4)   
        print W.get_element([1, 1, 0], [0, 1, 0])  #result is supposed to be 3.0
        print W.get_element([1, 1, 0], [0, 0, 1])  #result is supposed to be 4.0
    
    def test_permutation1():
        u=simple_itensor()[0]
        print u.matrix_view(2)
        u1 = u.permutation([1, 0, 2, 3])
        u2 = u.permutation([0, 1, 3, 2])
        print "\n", u1.matrix_view(2)
        print "\n", u2.matrix_view(2)

    @classmethod
    def permutation_buffon(cls):
        #u=simple_itensor()[0]
        u=cls.u.copy()
        u.data[:]=np.arange(u.totDim)
        print u.matrix_view(2)
        u1 = u.permutation([1, 0, 2, 3],use_buf=True)
        print iTensor.T_BUFFER[0].T[0]
        u2 = u.permutation([0, 1, 3, 2])
        print "\n", u1.matrix_view(2)
        print "\n", u2.matrix_view(2)



    @classmethod 
    def unit_tensor(cls):
        """   ----pass  """
        u = cls.u.copy()
        t = u.unit_tensor()
        print t.matrix_view(2).round(5)

        u22 = cls.u22.copy()
        t = u22.unit_tensor()
        #print u22
        print t.matrix_view()

    def conjugate():
        """   ---pass  perfect! just what i need"""
        u = simple_itensor()[0]
        print u.matrix_view(2)
        u1=u.conjugate(2)
        print u1.matrix_view(2)

    def test_direct_product():
        """  -----pass"""
        rank=2
        QSp=[QSp_base.copy() for i in xrange( rank )]
        QSp[0].add_to_quant_space(1, 1)
        #QSp[1].add_to_quant_space(1, 2)
        #QSp[0].add_to_quant_space(-1, 5)
        #QSp[1].add_to_quant_space(-1, 3)
        #print 'qqq', QSp[0]
        QSbase = QSp_base.copy
        totQN=QN_idendity.copy()
        u=iTensor(rank,QSp,totQN)

        u.data[:u.totDim] = np.arange(1, u.totDim + 1)
        #print u
        u1=u.matrix_view()
        print u1

        uu=u.direct_product(u)
        u2=uu.matrix_view()
        print u2
        print uu.rank

        print uu.__repr__(['data'])

    def direct_product2():
        rank = 2
        totQN = QN_idendity.copy()

        QSp = [QSp_base.copy() for i in xrange(rank)]
        
        sigma_x=iTensor(rank, QSp, totQN)
        sigma_x.data[0:2] = [1., -1.]
        
        #totQN.val = -1
        totQN.set_val(-1)
        I_2=iTensor(rank, QSp, totQN )
        I_2.data[0:2] = [1.0,1.0]

        #print sigma_x 
        #print I_2
        #print sigma_x.totDim
        #print QN_idendity
        xx= sigma_x.direct_product(I_2)
    
    def direct_prod_vs_contract(self):
        u = self.u
        self.u.data[:] = np.random.random(u.data.size)
        uu, nothing=self.u.contract(self.u, [1, 2, 3, 4], [5, 6, 7, 8])
        print uu.data[:5]
        uu = u.direct_product(u)
        print uu.data[:5]


    @staticmethod
    def direct_product_u1():
        """
        --- pass
        test direct_product for U1 symm.
        """
        rank=2
        QSp=[QSp_base.copy(),QSp_base.copy()]
        QSp[1].reverse()
        totQN=QN_idendity.copy()
        t=iTensor(rank,QSp,totQN)
        t.data[:]=np.arange(6)
        #print t
        t1=t.copy()
        print t.direct_product(t1)

        """
        --- pass
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
                 0.00E+00 0.00E+00 0.00E+00 0.10E+01 0.00E+00 0.00E+00 0.20E+01 0.30E+01 0.00E+00 0.20E+01 0.00E+00 0.30E+01 0.40E+01 0.60E+01 0.60E+01 0.90E+01
                T%Block_QN=     2     1     0     0
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     1     2     0     0
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     1     0     1     0
                 0.00E+00 0.40E+01 0.80E+01 0.12E+02
                T%Block_QN=     0     1     1     0
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     2     0     2     0
                 0.00E+00 0.50E+01 0.10E+02 0.15E+02
                T%Block_QN=     0     2     2     0
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     1     0     0     1
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     0     1     0     1
                 0.00E+00 0.40E+01 0.80E+01 0.12E+02
                T%Block_QN=     1     1     1     1
                 0.16E+02
                T%Block_QN=     0     0     2     1
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     2     1     2     1
                 0.20E+02
                T%Block_QN=     1     2     2     1
                 0.00E+00
                T%Block_QN=     2     0     0     2
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     0     2     0     2
                 0.00E+00 0.50E+01 0.10E+02 0.15E+02
                T%Block_QN=     0     0     1     2
                 0.00E+00 0.00E+00 0.00E+00 0.00E+00
                T%Block_QN=     2     1     1     2
                 0.00E+00
                T%Block_QN=     1     2     1     2
                 0.20E+02
                T%Block_QN=     2     2     2     2
                 0.25E+02
            ----End Tensor-------------------------------------------



            ----Begin iTensor----------------------------------------------
                rank:	4
                idx_dim:	81
                ndiv:	None
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
                data:	[0 0 0 0]: [ 0.  0.  0.  1.  0.  0.  2.  3.  0.  2.  0.  3.  4.  6.  6.  9.]
                [2 1 0 0]: [ 0.  0.  0.  0.]
                [1 2 0 0]: [ 0.  0.  0.  0.]
                [1 0 1 0]: [  0.   4.   8.  12.]
                [0 1 1 0]: [ 0.  0.  0.  0.]
                [2 0 2 0]: [  0.   5.  10.  15.]
                [0 2 2 0]: [ 0.  0.  0.  0.]
                [1 0 0 1]: [ 0.  0.  0.  0.]
                [0 1 0 1]: [  0.   4.   8.  12.]
                [1 1 1 1]: [ 16.]
                [0 0 2 1]: [ 0.  0.  0.  0.]
                [2 1 2 1]: [ 20.]
                [1 2 2 1]: [ 0.]
                [2 0 0 2]: [ 0.  0.  0.  0.]
                [0 2 0 2]: [  0.   5.  10.  15.]
                [0 0 1 2]: [ 0.  0.  0.  0.]
                [2 1 1 2]: [ 0.]
                [1 2 1 2]: [ 20.]
                [2 2 2 2]: [ 25.]

            ----End iTensor----------------------------------------------

        """

    @staticmethod
    def direct_product_u1_2():
        """
        another test
        test direct_product for U1 symm.
        """
        rank=2
        QSp=[QSp_base.copy(),QSp_base.copy()]
        QSp[1].reverse()
        
        totQN=QN_idendity.copy()
        totQN.set_val(1)
        tp=iTensor(rank,QSp,totQN)
        tp.data[:]=np.arange(1, tp.totDim + 1)

        rank=2
        QSp=[QSp_base.copy(),QSp_base.copy()]
        QSp[1].reverse()

        
        totQN=QN_idendity.copy()
        totQN.set_val(-1)
        tm=iTensor(rank,QSp,totQN)
        tm.data[:]=np.arange(1, tm.totDim + 1)

        print tp.direct_product(tm)

    def expand_u(self):
        pass
        #qn, qsp, null = init_System_QSp(self.symmetry)
        t = self.u.copy()
        t.data[:] = 1.0
        #print t.matrix_view()
        #qsp1 = QspU1.easy_init(qn=)
        new_qsp = [q.copy() for q in t.QSp]
        for q in new_qsp:
            q.Dims= [i*2 for i in q.Dims]
        #print new_qsp
        t2 = t.expand(new_qsp)
        #print t2.matrix_view()
        print t, t2
        print t2.data
    
    def index_merge(self):
        if 0:
            t = self.u.copy()
            t.data[:] = np.arange(t.data.size) + 1
            t1 = t.index_merge([2, 3])
            #t1.data[:] = 1.0
            #print t, t1
            print t.data, t1.data

        if 1:
            t = self.u22.copy()
            t.data[:] = np.arange(t.data.size) + 1
            t1 = t.index_merge([0, 1])
            #t1.data[:] = 1.0
            #print t, t1
            print t.data, t1.data

class performance_iTensor(object):
    pass
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
        t1=t.permutation(range(r,rank) + xrange(r), buffer=buf) 
        tensor_player.STATE = "play"
        def func(which):
            array_permutation.permute_player_fort = \
                    array_permutation.__getattribute__("permute_player_fort" + which)
            #print "do schedule dynamic"
            t0 = time.clock(); ta = time.time()
            for i in xrange(iter_times):
                t1=t.permutation(range(r,rank) + xrange(r), buffer=buf)
            tb = time.time(); t1 = time.clock()
            if which == "_parallel_runtime":
                #print os.environ["OMP_SCHEDULE"] 
                which  = which +  "  " + os.environ["OMP_SCHEDULE"] 
            print which, "\t", t1-t0, tb-ta


        
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
        for n in xrange(1, 8):
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
            for i in xrange(iter_times):
                T1.contract_core(T2, 2, data=buff)
            tb = time.time(); t1 = time.clock()
            if which == "_parallel_runtime":
                #print os.environ["OMP_SCHEDULE"] 
                which  = which +  "  " + os.environ["OMP_SCHEDULE"] 
            print which, "\t", t1-t0, tb-ta

        func("")
        #func("_paralell_ordered")
        #func("_paralell_critical")
        func("_paralell_test")
        #func("_paralell_reduction")
        #func("_paralell_reduction_1")
    
    def contract(self):
        self._contract(symmetry="U1", rank1=8, rank2=4, dim=16, nqn=3,  
                NUM_OF_THREADS=6, iter_times=10)


if __name__== "__main__":
    warnings.filterwarnings("ignore")
    if 0: 
        if 0: 
            ti = test_iTensor("Z2", 4)
            ti_e = test_iTensor("Travial")
            itf = iTensorFactory("Travial")
            pf = performance_iTensor("U1")
            #pauli_mat=itf.pauli_mat()
            #pf.permute()
            #pf.contract()
            #ti.index_merge()


            #t=iTensorFactory.simple(6, [8]*6, "U1")

            #itf.pauli_mat_1site()


            #itf.check_pauli_mat()
            #print ti_e.qsp_base
            #ti_e.contract()

            
            #print ti.instance(which="t0")
            #ti.dump()
            #ti.contract_core()
            #ti.contract_core_buff()
            #ti.contract_core1()
            #ti.contract_core2()
            #ti.permutation_buffon()
            #ti.contract()
            #ti.contract_large_tensor()
            #ti.decoration()
            #ti.direct_product_u1()   
            
            #ti.direct_product_u1_2()
            #ti.contract_U_H2()
            #ti.is_same_shape()
            #u=ti.use_buf()
            #print iTensor.buff_free()

            #tiu.init()
            #print tiu.unit_tensor()
            #ti.expand_u()
        
        if 0: 
            rank=4
            qsp = QspZ2.easy_init([1,-1], [1, 1 ]) 
            qsp = qsp.copy_many(4)
            qn = QnZ2.qn_id()
            t = iTensor(rank, qsp, qn) 
            t.data[: ] = xrange(t.data.size)
            print t.__repr__(fewer=0,keys=['data'])
            tm=t.index_merge_simple()
            print tm.__repr__(fewer=0,keys=['data'])
            print t.matrix_view()
            print tm.matrix_view()
            
        if 0: 
            symm = 'Z2'
            pau1 = iTensorFactory.pauli_mat_1site(symm)
            pau2 = iTensorFactory.pauli_mat_2site(symm)
            if 1: 
                print '----- 0 -----'
                s0 = pau1['sigma_0']
                s0 = pau1['sigma_0']
                s0i = s0.direct_product(s0)
                #s0i.show_data()
                print s0i.matrix_view()
                s0i_mer = s0i.index_merge_simple()
                #s0i_mer.show_data()
                print s0i_mer.matrix_view()
                print pau2['s00'].matrix_view()
            if 1: 
                print '----- z -----'
                sz = pau1['sigma_z']
                s0 = pau1['sigma_0']
                szi = sz.direct_product(s0)
                print szi.matrix_view()
                szi_mer = szi.index_merge_simple()
                print szi_mer.matrix_view()
                szi_mer = szi.index_merge_simple_2()
                print szi_mer.matrix_view()
               
                print pau2['szi'].matrix_view()
            if 1: 
                print '----- x -----'
                sx = pau1['sigma_x']
                s0 = pau1['sigma_0']
                sxi = sx.direct_product(s0)
                print sxi.matrix_view()
                sxi_mer = sxi.index_merge_simple()
                print sxi_mer.matrix_view()
                sxi_mer = sxi.index_merge_simple_2()
                print sxi_mer.matrix_view()
               
                print pau2['sxi'].matrix_view()
                
            if 1: 
                print '----- y -----'
                sy = pau1['sigma_y']
                s0 = pau1['sigma_0']
                syi = sy.direct_product(s0)
                print syi.matrix_view()
                syi_mer = syi.index_merge_simple()
                print syi_mer.matrix_view()
                print pau2['syi'].matrix_view()
  
    if 0: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #Test_nTensor('test_set_element'), 
           #Test_nTensor('test_permutation'), 
           #Test_nTensor('test_contract'), 
           
        ]
        for a in add_list: 
            suite.addTest(a)
        unittest.TextTestRunner().run(suite)
       
  
