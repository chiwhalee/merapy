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
from past.utils import old_div
from collections import OrderedDict
import os, sys
import unittest 
import pickle as pickle
import cProfile 
import time 
import numpy as np
import pstats 
import os.path
import warnings
#import pprocess


from scipy.linalg.blas import dgemm, sgemm, zgemm
try:
    import cupy as cp
    #from cupy.cuda import cublas
    import cupy.cublas as cublas 
except:
    pass


try:
    import triton
    import triton.language as tl
except:
    triton = None
    pass
   


import merapy.common_util as common_util

#from merapy import array_permutation
from merapy.context_util import redirect
from merapy.utilities import save, print_vars


num_of_instance = 0


__all__ = ["decorate_methods", "tensor_player", "set_STATE_end_1", "set_STATE_end_simple", "timer", "profileit", 'set_player_state_auto']
#__all__ = ["decorate_methods", "tensor_player", "set_STATE_end", "set_STATE_end_1", "timer", "profileit"]

"""
    design goal: 
        1. use new player   
        2. State Information Retention  in tape 
        3. not using many @ statement       ---done
"""

class Tape(dict):
    def __init__(self):
        self.calls= 0
        #self.tape = {}     
        self.calls = 0
        self.calls_tot = 0
        self.reach_tape_end = False
        self.tape_cleared = True
        self.STATE = 'stop'
        self.NEXT_STATE = None 
        self.PREV_STATE = None   #only for backward compatible 
        self.need_reset = False
    
    def __setitem__(self, k, v):
        self.calls += 1 
        self.calls_tot += 1 
        dict.__setitem__(self, self.calls, v)
    
    def __getitem__(self, k):
        self.calls += 1 
        return dict.__getitem__(self, self.calls)
        
    def reset(self):
        pass
        self.tape = {}   # blank tape 
        self.calls = 0   # record/play at #song 
        self.calls_tot = 0   #total num of 'songs' recorded 
        self.reach_tape_end = False 
        self.tape_cleared = False
    
    def show(self):
        res= '\ncontent of the tape:\n'
        res += '\tcalls_tot={0.calls_tot}\n'.format(self)
        res +=  '\t' + super(dict, self).__str__()[:1500]
        print(res)

        return res 
        
        
#TapeList = [Tape()]   # at least one tape 
TapeList = OrderedDict()
TapeList[0] = Tape()   #construct a default tape
        

def decorate_methods(decorator, meth_names):
    meth_names_dic = {
            'iTensor':
                [
                '__init__', 
                #'set_data_entrance', 
                'contract_core', 
                'prepare_leg', 
                'permutation', 
                ], 
            'QuantSpaceBase':
                [
                'copy'
                ], 
            #'Tensor_svd':['group_legs']   directly use tensor_player decorated see tensor_svd.py
        }
    meth_original_bac = {}
    
    meth_names_dic["QuantSpaceBase"].pop(0)

    def ClassDecorate(Class):
        class_name = Class.__name__
        meth_names = meth_names_dic[class_name]
        if 0:
            print( "%s are decorated in %s"%(meth_names, class_name))
        for attr in meth_names:
            meth_original_bac[class_name] = Class.__dict__[attr] #backup original method
            deced_meth = decorator(which=attr)(Class.__dict__[attr])
            
            if isinstance(Class.__dict__[attr], classmethod): #spcial treatment for classmethod
                deced_meth = classmethod(deced_meth)
            setattr(Class, attr, deced_meth) 
            #print "%s is decorated in %s. originally %s; now %s"%(
            #        attr, class_name, meth_original_bac[class_name], deced_meth)

        return Class
    
    return ClassDecorate

if 1: #define recorder and player
    def contract_core_recorder(func, self, T2, div, data=None, use_buf=False):
        """
            把T1，和T2的非零block 如果量子数组合相等则收缩
            locals:
                div: num. of legs to be contracted for each tensor
                buffer: use buffer to save data of T3
        """
        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        #tQN = self.totQN+T2.totQN
        tQN = self.totQN.__add__(T2.totQN)
        shift = rank1-div
        if 1:
            QSp = [self.QSp[i].copy() for i in range(shift)]
            QSp.extend([T2.QSp[i].copy() for i in range(div, rank2)])
        else:
            QSp = self.QSp[:shift]
            QSp.extend(T2.QSp[div:rank2])
        
        if rank3==0:
            #QSp = self.QSp[0].null()
            #QSp = [self.QSp[0].null()]
            QSp = []
        
        dtype = np.promote_types(self.dtype, T2.dtype)  # 自动根据 self.dtype 和 T2.dtype 推导最精准的输出类型（完美支持 32位/64位 和 实数/复数）
        if self.dtype != dtype:
            self.change_dtype(dtype)
        if T2.dtype != dtype:
            T2.change_dtype(dtype)
        
        
        T3= self.__class__(rank=rank3, QSp=QSp, totQN=tQN, 
                buffer=data, use_buf=use_buf, dtype=dtype, use_gpu=self.use_gpu)
        T3.data[:]=0.0
        
        nidx3 = 0
        iQN1=np.empty(self.rank + 1, np.int64)
        iQN2=np.empty(T2.rank + 1, np.int64)        
        iQN3=np.empty(T3.rank + 1, np.int64) # +1 to avoid T3.rank=0
        
        ind_count = 0
        
        contract_record = np.ndarray((self.nidx*T2.nidx, 6), np.int64)
        

        for idx2 in range(T2.nidx):
            iQN2[0] = 0  #!for rank=0
            iQN2[0:rank2]=T2.Addr_idx[0:rank2,idx2]
            p2 = T2.Block_idx[0,idx2]
            
            Dim2 = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div, rank2)], dtype=np.int64)
            Dimc = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div)], dtype=np.int64)
            data2 = T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F') 
            
            for idx1 in range(self.nidx):
                iQN1[0] = 1 #!for rank=0
                iQN1[0:rank1]=self.Addr_idx[0:rank1,idx1]
                p1 = self.Block_idx[0, idx1]
                iseq = np.all(iQN1[shift:shift+div] == iQN2[0:div])
                if not iseq:
                    #如果量子数组合相等则收缩
                    continue
                Dim1 = np.prod([self.QSp[i].Dims[iQN1[i]] for i in range(shift)], dtype=np.int64)
                
                iQN3[0] = 1 #!for rank=1
                iQN3[0:shift] = iQN1[0:shift]
                iQN3[shift:rank3] = iQN2[div:rank2]
                p3=T3.get_position(iQN3[:T3.rank])
                idx3 = T3.idx[p3]
                p3 = T3.Block_idx[0,idx3]
                
                contract_record[ind_count][:] = (p1, p2, p3, Dim1, Dim2, Dimc)
                ind_count  += 1  
        

                data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                
                data3 = T3.data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order='F')    
                
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
                
                
        
        tensor_player.the_tape[tensor_player.the_tape.calls] = (
                'contract', contract_record[:ind_count, :], ind_count)
        #print('xxxxxxxxxx contr', tensor_player.the_tape.calls)
        return T3

    def init_recorder(func,  self, rank=None,  QSp=None, totQN=None, order="F",  dtype=np.float64, 
            buffer=None, use_buf=False, index_data=True, has_data=True, 
            use_gpu=0):
        """
            这种做法有些过于激进, 更安全的是下面的 data_entrance_recorder/player 
            note1:
                one may pop out 'data' or not. Both can run. The draw back of the later 
                is that it will use much more memory,  so I use the former. 
        
        """
        func(self, rank=rank,  QSp=QSp, totQN=totQN, order="F",  dtype=dtype, 
            buffer=buffer, use_buf=use_buf, index_data=index_data, has_data=has_data, 
            use_gpu=use_gpu)
        if 1:  #note1
            struct_dict = self.__dict__.copy()   #here copy is NEED, or else self will lost data after pop 
            if has_data:
                struct_dict.pop('data')
            #struct_dict.pop("buf_ref")
        else:
            struct_dict = self.__dict__
        
        #print_vars(vars(),  ['self.buf_ref'])
         
        tensor_player.the_tape[tensor_player.the_tape.calls] = struct_dict
    
    #@profile
    def init_player(self, rank=None,  QSp=None, totQN=None, order="F",  dtype=np.float64, 
            buffer=None, use_buf=False, index_data=True, has_data=True, use_gpu=0):
        """
            note1:
                Although one may remove .copy, which sould still work.
                But, adding a .copy() is safer and better.  If not copy
                here, 'data' will also refered to by the tape. This
                will cause problem in a few circumstances. For example:
                1) If tape is stored, then data will be sotred also. 2)
                if make contraction recursively,  like t =
                t.contract(t),  it will yeild wrong results. So,  I use copy here.
        """
        self.__dict__ = tensor_player.the_tape[tensor_player.the_tape.calls].copy()  #note1 
        self.buf_ref = np.array([-1, -1], np.int64)
        if has_data:
            if  buffer is None and use_buf:                    
                if self.use_gpu != 1:
                    #buffer = self.buffer_assign(data_size=self.totDim if dtype==float else self.totDim*2) 
                    buffer = self.buffer_assign(data_size=self.totDim, item_size=np.dtype(dtype).itemsize) 
                    
                else:
                    #buffer = self.buffer_assign_gpu(data_size=self.totDim if dtype==float else self.totDim*2)  
                    buffer = self.buffer_assign_gpu(data_size=self.totDim, item_size=np.dtype(dtype).itemsize) 
                    
            if self.use_gpu != 1:
                self.data = np.ndarray(self.totDim, buffer=buffer, dtype=dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered
            else:
                #buffer = buffer if buffer is None else buffer.data
                if isinstance(buffer, cp.ndarray):
                    buffer = buffer.data 
                self.data = cp.ndarray(self.totDim, memptr=buffer, dtype=dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered
        
        #self.__dict__ = tensor_player.the_tape[tensor_player.the_tape.calls]
        #if self.dtype != dtype:   # It could happen that dtype change in later iterations, e.g. in TDVP algrithom
        #    self.data = np.ndarray(self.totDim, buffer=buffer, dtype=dtype, order="C")   #as a mater of fact, 1D array is both C and F ordered               
            

    def data_entrance_recorder(self, order="F"):
        rank, QSp, totQN = self.rank, self.QSp, self.totQN
        rank_1 = rank if rank != 0 else 1 
        
        if rank == 0: 
            QSp = [totQN.qsp_class().null()]   # only use it temporarilly to generate idx
       
        #self.Dims = np.array([QSp[i].totDim for  i in xrange(rank_1)])
            
        temp = 1   
        for i in range(rank_1):
            temp *= QSp[i].nQN
        self.idx_dim = temp   # 量子数组合 总数目 

        #attention_this_may_be_wrong 在python中  -1对应着最后一个元素，而Fortran中什么也不对应, 所以改成下面的
        #self.idx=np.array([int(self.idx_dim)]*self.idx_dim, int)   #-1                
        self.idx=np.ndarray((self.idx_dim, ), int)   #-1                
        self.idx[: ] = self.idx_dim
        self.Block_idx=np.ndarray((3, self.idx_dim),dtype=int)
        
        #iQN[i]用作leg i 上的量子数 计数
        iQN = np.zeros(rank_1, dtype=int)   #iQN 用于给量子数组合编号
        # 实际使用的addr_inx的长度为 self.nidx
        self.Addr_idx=np.ndarray((rank_1, self.idx_dim),dtype=int)        
            
        nidx=0
        totDim=0
            
        tQN_r= self.totQN.copy()  # here must copy
        #tQN_r.reverse()

        for p in range(self.idx_dim):
            tQN = QSp[0].QNs[iQN[0]]
            #这里计算了总量子数 tQN, 然后
            #判断量子数组合是否满足指定的对称性要求
            for i in range(1,rank):
                tQN = tQN+QSp[i].QNs[iQN[i]]
            if tQN==tQN_r:
                #q895  why not tQN == totQN?  因为operater- state duality, 从而在操作下变换正好相反？
                #其中包含协变/反变的意味
                
                d = 1
                for i in range(rank):
                    d = d*QSp[i].Dims[iQN[i]]
                    #print "ddddd d", i, p, d,iQN[i], len(QSp[i].Dims)
                #计算某一block的总数据量

                self.idx[p] = nidx
                #0 is position,  1 is size, 2 is position in quantum number combinations
                self.Block_idx[0,nidx] = totDim
                #在该block在self.data中的position
                self.Block_idx[1,nidx] = d
                #block变成1d数组的长度
                self.Block_idx[2,nidx] = p
                #记录不为0的量子数组合，在所有量子数组合中的位置
                self.Addr_idx[0,nidx] = 0 #for rank=0
                self.Addr_idx[0:rank,nidx] = iQN[0:rank]
                #Addr_idx这个二维数组的每一列实际上是所有非零block的量子数的编号(而不是量子数点值！)的组合
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
        tensor_player.the_tape[tensor_player.the_tape.calls] = ( self.Addr_idx.copy(), self.Block_idx.copy(), 
                self.idx.copy(), self.nidx, self.idx_dim, self.totDim)

    def data_entrance_player(self, order="F"):
        (self.Addr_idx, self.Block_idx, self.idx, self.nidx, 
                self.idx_dim, self.totDim)=tensor_player.the_tape[tensor_player.the_tape.calls]

    def permute_recorder(func,  self, P, buffer=None, use_buf=False):
        rank = self.rank
        #permute QSp, QNs
        #first is more efficient; second is more robust
        if 0:
            QSp=[self.QSp[P[i]] for i in range(rank)]
            totQN = self.totQN
        else:
            QSp=[self.QSp[P[i]].copy() for i in range(rank)]
            totQN = self.totQN.copy()

        #Tp=iTensor(rank, QSp, totQN, buffer=buffer, use_buf=use_buf)
        Tp=self.__class__(rank, QSp, totQN, buffer=buffer, 
                dtype=self.dtype, use_buf=use_buf, 
                use_gpu = self.use_gpu)
        pos=np.empty(self.rank,"int")
        Dims=np.empty(self.rank,"int")

        tape_ind = np.ndarray((self.nidx, 3), np.int32)
        #tape_dim = np.ndarray((self.nidx, 32), np.int32)
        #tape_ord = np.ndarray((self.nidx, 32), np.int32)
        tape_dim = np.ndarray((self.nidx, rank), np.int32)
        tape_ord = np.ndarray((self.nidx, rank), np.int32)
        
        for n  in range(self.nidx):
            pidx = self.Block_idx[0,n]
            totDim = self.Block_idx[1,n]
            pos[0] = 0  # for rank=0
            pos[0:rank] = self.Addr_idx[0:rank,n]
            np1 = 0
            for i  in range(rank-1, 0, -1):
                np1 = (pos[P[i]]+np1)*Tp.QSp[i-1].nQN
                Dims[i] = self.QSp[i].Dims[pos[i]]
            i = 0
            np1 = np1+pos[P[i]]
            Dims[i] = self.QSp[i].Dims[pos[i]]
            qidx = Tp.Block_idx[0,Tp.idx[np1]]            
            
            data = self.data[pidx:pidx+totDim]   #.copy()
            
            if 0:
                #data = Tp.data[qidx:qidx+totDim]
                #array_permutation.array_permutation_inplace(temp, self.rank, Dims, P, data)
                Tp.data[qidx:qidx+totDim]=array_permutation.array_permutation(
                        data, rank,Dims,P)
            else:
                data = data.reshape(Dims, order='F')
                #data = np.transpose(data, axes=P)
                data = data.transpose(P)
                Tp.data[qidx:qidx+totDim]= data.ravel(order='F')
                
            
            tape_ind[n][0:3] = (pidx, qidx, totDim)
            tape_dim[n][:] = Dims[:rank]
            tape_ord[n][:] = P[:rank]
        if 0:
            tape_sort = tape_ind[:, 2].argsort()   ##I may implement this later, sort the tape according to totDim
            #tape_sort.sort()
            tensor_player.the_tape[tensor_player.the_tape.calls] = ('permute', self.nidx, tape_ind[tape_sort], tape_dim[tape_sort], tape_ord[tape_sort])
        else:
            tensor_player.the_tape[tensor_player.the_tape.calls] = ('permute', self.nidx, tape_ind, tape_dim, tape_ord)  #here I insert 'permute' for easier debuging
            
        return Tp

    #@profile
    def permute_player(self, P, buffer=None, use_buf=False):
        Tp=self.__class__(None, None, None, dtype=self.dtype, 
                buffer=buffer, use_buf=use_buf, use_gpu=self.use_gpu)
        
        _, nidx, tape_ind, tape_dim, tape_ord = tensor_player.the_tape[tensor_player.the_tape.calls]
        
        if 0:
            if self.data.dtype == np.float64: 
                array_permutation.permute_player_fort(
                            tape_ind, tape_dim, tape_ord, self.data, Tp.data)
            else:  #complex 
                array_permutation.complex_permute_player_fort(
                            tape_ind, tape_dim, tape_ord, self.data, Tp.data)
        if 1:
            for ind in range(nidx):
                pidx, qidx, totDim  = tape_ind[ind]
                dims = tape_dim[ind]
                P = tape_ord[ind]
                data = self.data[pidx:pidx+totDim]
                data = data.reshape(dims, order='F')
                #data = np.transpose(data, axes=P)
                data = data.transpose(P)
                Tp.data[qidx:qidx+totDim] = data.ravel(order='F')

        return Tp

    #@profile
    def contract_core_player(self, T2, div, data=None, use_buf=False):
        """
            I have tried to make the following parallel, using either python or fortran code.
                python:
                    contract_core_player_parallel_1, 2, 3
                    they are all slow
                fortran:
                    contract_core_player_fort_paralell_critical,  result correct
                    contract_core_player_fort_paralell_reduction,  result correct
                    contract_core_player_fort_paralell_ordered,  result correct
                    
                    The results are all correct, but I don't remember why I
                    did not use these parallel versions. I may check these out in future. 
            UPDATE: actually,  in practice the parallel version may be slower,  so not needed. 
                The reason is that the sizes of the data blocks are
                very inhomogenous, so the lagest data block would be
                bottle neck. 
            
        """
            
        
        dtype = np.promote_types(self.dtype, T2.dtype)  # 自动根据 self.dtype 和 T2.dtype 推导最精准的输出类型（完美支持 32位/64位 和 实数/复数）
        if self.dtype != dtype:
            self.change_dtype(dtype)
        if T2.dtype != dtype:
            T2.change_dtype(dtype)
           
        T3 = self.__class__(rank=None, QSp=None, totQN=None, 
                buffer=data, dtype=dtype, use_buf=use_buf)
        
        T3.data[:]=0.0
        self_data = self.data
        T2_data = T2.data
        T3_data = T3.data
        #print_vars(vars(),  ['self.use_gpu'])
        #if self.use_gpu:
        #    self_data = cp.asarray(self_data)
        #    T2_data = cp.asarray(T2_data)
        #    T3_data = cp.asarray(T3_data)
        
        _, rec, num_rec = tensor_player.the_tape[tensor_player.the_tape.calls]
        
        if 0:
            if dtype == np.float64:  
                common_util.contract_core_player_fort(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                #common_util.contract_core_player_fort_paralell_critical(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                #common_util.contract_core_player_fort_paralell_ordered(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                #common_util.contract_core_player_fort_paralell_reduction(self.data, T2.data, T3.data, rec, num_rec=num_rec)
            else: 
                common_util.contract_core_player_fort_complex(self.data, T2.data, T3.data, rec, num_rec=num_rec)
        
        if 0:
             for ind in range(num_rec):
                p1, p2, p3, Dim1, Dim2, Dimc = rec[ind]
                data1 = self_data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2 = T2_data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                data3 = T3_data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order='F')    
                
                use_gpu = self.use_gpu
                if self.use_gpu == 2 and (Dim1*Dim2*Dimc) < self.USE_GPU_MUL_LIM:
                    use_gpu = 0
                    
                common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                        dtype = dtype, use_gpu=use_gpu)
                
        if self.use_gpu ==0:
             for ind in range(num_rec):
                p1, p2, p3, Dim1, Dim2, Dimc = rec[ind]
                data1 = self_data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2 = T2_data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                data3 = T3_data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order='F')    
                    
                common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                        dtype = dtype, use_gpu=0)
                
        elif self.use_gpu == 1:
            if not self.USE_BATCHED_GEMM_FOR_GPU:
                for ind in range(num_rec):
                    p1, p2, p3, Dim1, Dim2, Dimc = rec[ind]
                    data1 = self_data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                    data2 = T2_data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                    data3 = T3_data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order='F')    
                        
                    common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                            dtype = dtype, use_gpu=1)
            else:
                #contract_core_player_cublas_batched(self,  T2,  T3,  rec,  num_rec)
                contract_core_player_triton(self,  T2,  T3,  rec,  num_rec)       
                
        elif self.use_gpu == 2:
             for ind in range(num_rec):
                p1, p2, p3, Dim1, Dim2, Dimc = rec[ind]
                data1 = self_data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2 = T2_data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                data3 = T3_data[p3:p3+Dim1*Dim2].reshape((Dim1,Dim2), order='F')    
                
                if self.use_gpu == 2 and (Dim1*Dim2*Dimc) < self.USE_GPU_MUL_LIM:
                    use_gpu = 0
                else:
                    use_gpu = 2
                common_util.gemm_all(data1, data2, data3, alpha=1.0, beta=1.0,
                        dtype = dtype, use_gpu=use_gpu)
                
        else:
            raise  ValueError(self.use_gpu)
                
        
        #if self.use_gpu:
        #    #T3.data = cp.asnumpy(T3_data) 
        #    T3.data = T3_data.get()
            
        return T3
    
    def contract_core_player_cublas_batched(self, T2, T3, rec, num_rec):
        """
            基于 cuBLAS 指针数组(Pointer Array)的非均匀块 Batched GEMM 实现
            消除 for ind in range(num_rec) 的 Python 循环延迟
        """
        raise NotImplemented  #会报错，这个办法根本行不通
        
        if num_rec == 0:
            return T3

        # 提取三者显存的一维裸指针物理基地址 (十六进制首地址)
        #base_A = self.data.data.ptr
        #base_B = T2.data.data.ptr
        #base_C = T3.data.data.ptr
        
        base_A = int(self.data.data.ptr)
        base_B = int(T2.data.data.ptr)
        base_C = int(T3.data.data.ptr)                
        
        # 获取数据类型的字节大小
        itemsize = self.data.itemsize

        # 3. 利用 CuPy 的高级切片，直接在 GPU 侧把磁带里的相对偏移 p1, p2, p3 转换为绝对显存物理地址数组
        # rec 每一行是: [p1, p2, p3, Dim1, Dim2, Dimc]
        
        #p1_arr = rec[:num_rec, 0]
        #p2_arr = rec[:num_rec, 1]
        #p3_arr = rec[:num_rec, 2]
        
        p1_arr = cp.array(rec[:num_rec, 0], dtype=cp.int64) # 内存偏移量应该用int64
        p2_arr = cp.array(rec[:num_rec, 1], dtype=cp.int64)
        p3_arr = cp.array(rec[:num_rec, 2], dtype=cp.int64)                
        
        A_offsets = base_A + p1_arr * itemsize
        B_offsets = base_B + p2_arr * itemsize
        C_offsets = base_C + p3_arr * itemsize

        #  准备 cuBLAS 所需要的元数据数组 (M, N, K)
        # cuBLAS 接收 C 类型的 int* 数组
        #m_arr = rec[:num_rec, 3].astype(cp.int32)
        #n_arr = rec[:num_rec, 4].astype(cp.int32)
        #k_arr = rec[:num_rec, 5].astype(cp.int32)
        
        m_arr = cp.array(rec[:num_rec, 3], dtype=cp.int32)  # 注意 gemm接口要求 int32
        n_arr = cp.array(rec[:num_rec, 4], dtype=cp.int32)
        k_arr = cp.array(rec[:num_rec, 5], dtype=cp.int32)                
      
        # 针对不规则矩阵，cuBLAS 要求每一块的领先维度（Leading Dimension）也必须作为数组传入
        # 因为是一维扁平紧凑排布的，LD 刚好等于各自的行数
        lda_arr = m_arr.copy()
        ldb_arr = k_arr.copy()
        ldc_arr = m_arr.copy()

        #调用 CuPy 低级的 cuBLAS 句柄（Handle）进行批处理发射
        #handle = cublas.get_handle()
        handle = cp.cuda.device.get_cublas_handle()
        
        # 根据数据类型动态路由到具体的 C++ 算子
        # 注意：cuBLAS 默认是 Fortran 连续（列优先），如果语义是 C 优先（行优先）
        # 可以利用数学公式 A @ B = (B^T @ A^T)^T，或者直接在这里调整 'N' 和 'T' 的转置参数
        # 这里以默认的内存排布为准：
        transA = cublas.CUBLAS_OP_N
        transB = cublas.CUBLAS_OP_N
        
        # 设置阿尔法和贝塔系数的指针
        alpha = np.array([1.0], dtype=self.data.dtype)
        beta  = np.array([0.0], dtype=self.data.dtype)
        
        # 根据 dtype 分流到底层的汇编多态算子
        if self.data.dtype == cp.float64:
            cublas.dgemmBatched(
                handle, transA, transB, 
                m_arr.data.ptr, n_arr.data.ptr, k_arr.data.ptr,
                alpha.ctypes.data, A_offsets.data.ptr, lda_arr.data.ptr,
                B_offsets.data.ptr, ldb_arr.data.ptr,
                beta.ctypes.data, C_offsets.data.ptr, ldc_arr.data.ptr, 
                num_rec
            )
        elif self.data.dtype == cp.complex128:
            cublas.zgemmBatched(
                handle, transA, transB, 
                m_arr.data.ptr, n_arr.data.ptr, k_arr.data.ptr,
                alpha.ctypes.data, A_offsets.data.ptr, lda_arr.data.ptr,
                B_offsets.data.ptr, ldb_arr.data.ptr,
                beta.ctypes.data, C_offsets.data.ptr, ldc_arr.data.ptr, 
                num_rec
            )
        elif self.data.dtype == cp.float32:
            cublas.sgemmBatched(
                handle, transA, transB, 
                m_arr.data.ptr, n_arr.data.ptr, k_arr.data.ptr,
                alpha.ctypes.data, A_offsets.data.ptr, lda_arr.data.ptr,
                B_offsets.data.ptr, ldb_arr.data.ptr,
                beta.ctypes.data, C_offsets.data.ptr, ldc_arr.data.ptr, 
                num_rec
                    ) 
        elif self.data.dtype == cp.complex64:
            raise NotImplemented
            #cublas.cgemmBatched(...) # 同理分流

        return T3


    # 1. 编写 Triton 硬件核心 Kernel
    @triton.jit
    def _triton_symmetric_gemm_kernel_bac(
        A_ptr, B_ptr, C_ptr, rec_ptr,
        stride_rec_row, # 磁带二维表的步长
        alpha,
        BLOCK_SIZE_M: tl.constexpr, 
        BLOCK_SIZE_N: tl.constexpr, 
        BLOCK_SIZE_K: tl.constexpr):
        # 🚀 绝妙之处：每个 Thread Block 处理磁带录制好的一个独立 Block 任务！
        task_id = tl.program_id(0)
        
        # 定位到当前任务在磁带表里的整行首地址
        current_rec_ptr = rec_ptr + task_id * stride_rec_row
        
        # 像读取汇编寄存器一样，瞬间捞出该 Block 的所有物理元数据
        p1 = tl.load(current_rec_ptr + 0)
        p2 = tl.load(current_rec_ptr + 1)
        p3 = tl.load(current_rec_ptr + 2)
        M  = tl.load(current_rec_ptr + 3)
        N  = tl.load(current_rec_ptr + 4)
        K  = tl.load(current_rec_ptr + 5)
        
        # 计算当前 Block 在三串大一维显存中的绝对首地址指针
        local_A_ptr = A_ptr + p1
        local_B_ptr = B_ptr + p2
        local_C_ptr = C_ptr + p3
        
        # 生成当前 Thread Block 内部的局域并行掩码网格 (Tiling)
        offs_m = tl.arange(0, BLOCK_SIZE_M)
        offs_n = tl.arange(0, BLOCK_SIZE_N)
        offs_k = tl.arange(0, BLOCK_SIZE_K)
        
        # 针对不规则非均匀块，生成越界守护掩码（Mask），杜绝 ILLEGAL_ADDRESS
        mask_m = offs_m < M
        mask_n = offs_n < N
        
        # 在 GPU 寄存器（Register）里开辟一块干净的累加清空区
        accumulator = tl.zeros((BLOCK_SIZE_M, BLOCK_SIZE_N), dtype=tl.float64)
        
        # 沿着收缩键合维度 K 进行分块迭代
        for k in tl.range(0, tl.cdiv(K, BLOCK_SIZE_K)):
            k_remaining = K - k * BLOCK_SIZE_K
            mask_k = offs_k < k_remaining
            
            # 计算局域指针矩阵
            # 假设内存是标准 C 连续排布
            a_ptrs = local_A_ptr + (offs_m[:, None] * K + (k * BLOCK_SIZE_K + offs_k[None, :]))
            b_ptrs = local_B_ptr + ((k * BLOCK_SIZE_K + offs_k[:, None]) * N + offs_n[None, :])
            
            # 协同搬运：同步从全局显存加载数据到 Shared Memory / Register
            a = tl.load(a_ptrs, mask=(mask_m[:, None] & mask_k[None, :]), other=0.0)
            b = tl.load(b_ptrs, mask=(mask_k[:, None] & mask_n[None, :]), other=0.0)
            
            # 核心：调用 Tensor Core 硬件乘法器发射局域矩阵乘法
            accumulator += tl.dot(a, b)
            
        # 🌟 Triton 的大招：算子融合 Epilogue（在寄存器里直接做完物理缩放，再写回显存）
        accumulator = accumulator * alpha
        
        # 终点着陆：安全写回 C 矩阵的一维扁平显存
        c_ptrs = local_C_ptr + (offs_m[:, None] * N + offs_n[None, :])
        tl.store(c_ptrs, accumulator, mask=(mask_m[:, None] & mask_n[None, :]))
        

    # 2. 在 Player 中包装并调用 Triton 编译器入口
    def contract_core_player_triton_bac(self, T2, T3, rec, num_rec):
        """
        Triton 版本的图回放执行器
        """
        import cupy as cp
        if num_rec == 0: return T3
        
        # Triton 必须接收标准的 PyTorch/Triton 认识的指针或者通过 CuPy 桥接
        # 确保你的磁带 rec 本身已经作为 Tensor 驻留在 GPU 上
        if not isinstance(rec, cp.ndarray):
            rec_gpu = cp.array(rec[:num_rec], dtype=cp.int64)
        else:
            rec_gpu = rec[:num_rec].astype(cp.int64)

        # 动态分析当前最大 Block 的上限，来决定 Triton 线程块的编译 Tiling 大小
        # 选择最邻近的 2 的幂次方以适配硬件对齐
        max_m = int(cp.max(rec_gpu[:, 3]).item())
        max_n = int(cp.max(rec_gpu[:, 4]).item())
        max_k = int(cp.max(rec_gpu[:, 5]).item())
        
        def next_power_of_2(x):
            return 1 if x == 0 else 2**(x - 1).bit_length()

        # 编译期常量，决定每个 SM 开辟多少共享内存
        BLOCK_M = max(16, min(128, next_power_of_2(max_m)))
        BLOCK_N = max(16, min(128, next_power_of_2(max_n)))
        BLOCK_K = max(16, min(128, next_power_of_2(max_k)))

        # 🚀 发射：Grid 大小就是总的任务数，实现完全的 Block 间硬件并行
        grid = (num_rec, )
        
        _triton_symmetric_gemm_kernel[grid](
            self.data.data.ptr, 
            T2.data.data.ptr, 
            T3.data.data.ptr, 
            rec_gpu.data.ptr,
            rec_gpu.strides[0], # 传入磁带表的行步长
            1.0, # alpha
            BLOCK_SIZE_M=BLOCK_M, 
            BLOCK_SIZE_N=BLOCK_N, 
            BLOCK_SIZE_K=BLOCK_K
        )
        
        return T3


    @triton.jit
    def _triton_symmetric_gemm_float32_kernel(
        A_ptr, B_ptr, C_ptr, rec_ptr,
        stride_rec_row, 
        alpha,
        BLOCK_SIZE_M: tl.constexpr, 
        BLOCK_SIZE_N: tl.constexpr, 
        BLOCK_SIZE_K: tl.constexpr
        ):
        # 🌟 每一个 Thread Block 独立认领磁带（rec_ptr）中的一行任务
        task_id = tl.program_id(0)
        
        # 1. 构造标准的 2 的幂（8元素）一维向量网格与掩码，规避 3.1 前端 AST 审查
        offs_rec = tl.arange(0, 8)
        mask_rec = offs_rec < 6
        
        # 2. 算出一维绝对字节地址（因为 rec_ptr 是 int64 磁带，每个元素占 8 字节）
        rec_offsets_int = task_id * stride_rec_row * 8 + offs_rec * 8
        rec_abs_addr = rec_ptr + rec_offsets_int
        
        # 3. 强转为指针向量，一枪把当前任务的 6 个元数据捞进寄存器
        rec_ptrs = tl.cast(rec_abs_addr, tl.pointer_type(tl.int64))
        rec_data = tl.load(rec_ptrs, mask=mask_rec, other=0)
        
        # 无损提取标量（在 3.1 极度挑食的 PassManager 下最稳健的硬件条件选通写法）
        p1 = tl.sum(tl.where(offs_rec == 0, rec_data, 0))
        p2 = tl.sum(tl.where(offs_rec == 1, rec_data, 0))
        p3 = tl.sum(tl.where(offs_rec == 2, rec_data, 0))
        M  = tl.sum(tl.where(offs_rec == 3, rec_data, 0))
        N  = tl.sum(tl.where(offs_rec == 4, rec_data, 0))
        K  = tl.sum(tl.where(offs_rec == 5, rec_data, 0))
        
        # 4. 🌟 像素级物理对账：因为换成了 float32，大张量里的每个元素占 4 字节！
        # 所以基地址变址要乘以 4
        local_A_ptr = tl.cast(A_ptr, tl.pointer_type(tl.float32)) + p1
        local_B_ptr = tl.cast(B_ptr, tl.pointer_type(tl.float32)) + p2
        local_C_ptr = tl.cast(C_ptr, tl.pointer_type(tl.float32)) + p3
        
        # 5. 经典 2D Tiling 切片网格
        offs_m = tl.arange(0, BLOCK_SIZE_M)
        offs_n = tl.arange(0, BLOCK_SIZE_N)
        offs_k = tl.arange(0, BLOCK_SIZE_K)
        
        mask_m = offs_m < M
        mask_n = offs_n < N
        
        # 初始化单精度累加器
        accumulator = tl.zeros((BLOCK_SIZE_M, BLOCK_SIZE_N), dtype=tl.float32)
        
        # 沿着 K 维度进行经典分块迭代
        for k in tl.range(0, tl.cdiv(K, BLOCK_SIZE_K)):
            k_remaining = K - k * BLOCK_SIZE_K
            mask_k = offs_k < k_remaining
            
            # 二维滑窗寻址（严格遵循 Row-Major 原理）
            a_ptrs = local_A_ptr + (offs_m[:, None] * K + (k * BLOCK_SIZE_K + offs_k[None, :]))
            b_ptrs = local_B_ptr + ((k * BLOCK_SIZE_K + offs_k[:, None]) * N + offs_n[None, :])
            
            a = tl.load(a_ptrs, mask=(mask_m[:, None] & mask_k[None, :]), other=0.0)
            b = tl.load(b_ptrs, mask=(mask_k[:, None] & mask_n[None, :]), other=0.0)
            
            # 🌟 炸裂输出：彻底激活 RTX 3090 Tensor Core 的 Ampere mma.sync 核心管线！
            
            # 强制让 PyTorch / CUDA 后端关掉 TF32 欺骗，回归真正的单精度物理硬件管线
            accumulator += tl.dot(a, b,  allow_tf32=False)
            
        accumulator = accumulator * alpha
        
        # 写回最终单精度矩阵
        c_ptrs = local_C_ptr + (offs_m[:, None] * N + offs_n[None, :])
        tl.store(c_ptrs, accumulator, mask=(mask_m[:, None] & mask_n[None, :]))


    def contract_core_player_triton(self, T2, T3, rec, num_rec):
        """
        【满血物理对账版】Triton 发射器：强制显存连续化，斩断一切变址串位
        """
        import cupy as cp
        if num_rec == 0: return T3
        
        # 🌟 1. 核心拯救行动：强制对大张量进行物理显存规整（C_CONTIGUOUS）
        # 如果它们是转置后的 View，这几行会强制在显存里重新开辟连续空间并拷贝，确保一维变址绝对正确
        if hasattr(self.data, 'flags') and not self.data.flags['C_CONTIGUOUS']:
            self.data = cp.ascontiguousarray(self.data)
        if hasattr(T2, 'flags') and not T2.flags['C_CONTIGUOUS']:
            T2 = cp.ascontiguousarray(T2)
        # T3 作为输出缓冲区，也必须保证在物理上是严格连续的一维空间
        if hasattr(T3, 'flags') and not T3.flags['C_CONTIGUOUS']:
            T3 = cp.ascontiguousarray(T3)

        # 2. 确保磁带本身规整
        if not isinstance(rec, cp.ndarray):
            rec_gpu = cp.array(rec[:num_rec], dtype=cp.int64)
        else:
            rec_gpu = rec[:num_rec].astype(cp.int64)

        # 3. 动态自适应推算最优 Tiling 尺寸
        max_m = int(cp.max(rec_gpu[:, 3]).item())
        max_n = int(cp.max(rec_gpu[:, 4]).item())
        max_k = int(cp.max(rec_gpu[:, 5]).item())
        
        def next_power_of_2(x):
            return 1 if x == 0 else 2**(x - 1).bit_length()

        BLOCK_M = max(16, min(128, next_power_of_2(max_m)))
        BLOCK_N = max(16, min(128, next_power_of_2(max_n)))
        BLOCK_K = max(16, min(128, next_power_of_2(max_k)))

        grid = (num_rec, )
        stride_rec_row = int(rec_gpu.strides[0] // 8)

        # 🌟 4. 判别当前正在运行的数据类型，精准动态匹配显存字节步长 (float32=4字节, float64=8字节)
        # 之前如果写死了 8 而上层传入了 float32，或者写死了 4 而上层是 float64，寻址直接会错乱 2 倍或 0.5 倍！
        element_size = self.data.dtype.itemsize  # 自动获取 4 或 8

        # 5. 根据真实数据类型，调用对应精度的内核
        # 如果上层是 float64，建议使用之前为你准备的带有 `tl.dot(..., allow_tf32=False)` 或广播归约的核函数
        # 这里以自适应寻址为例，在内核中配合正确的 element_size：
        
        # 重新核对：内核中计算指针时，由于 Triton 内部 `+ p1` 是基于指针类型的（Pointer Arithmetic）
        # 在 Triton 内部：指针 + 1 自动代表跳过一个元素（即自动乘以了 sizeof(dtype)）！
        # 🌟 致命对账：如果你之前在内核里写了 `local_A_ptr = A_ptr + p1 * 8`，那就是重复乘了 8 字节，指针直接飞到了九霄云外！
        
        _triton_symmetric_gemm_float32_kernel[grid](
            int(self.data.data.ptr), 
            int(T2.data.data.ptr), 
            int(T3.data.data.ptr), 
            int(rec_gpu.data.ptr),
            stride_rec_row,
            1.0, 
            BLOCK_SIZE_M=BLOCK_M, 
            BLOCK_SIZE_N=BLOCK_N, 
            BLOCK_SIZE_K=BLOCK_K
        )
        
        return T3


    def prepare_leg_recorder(func, self,T2, V1, V2, info=0):
        """
        """
        res = func(self, T2, V1, V2, info)
        tensor_player.the_tape[tensor_player.the_tape.calls] = ('prepare', res)
        return res

    def prepare_leg_player(self,T2, V1, V2, info=0):
        """
        """
        return tensor_player.the_tape[tensor_player.the_tape.calls][1] 

    def contract_core_player_parallel_1(self, T2, div, data=None, use_buf=False):
        """
        not work
        a parallel version
        """
        if 0:
            from . import tensor_py
            iTensor = tensor_py.iTensor
        
        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        tQN = self.totQN+T2.totQN
        shift = rank1-div
        
        QSp = self.QSp[:shift]
        QSp.extend(T2.QSp[div:rank2])
        
        #T3= iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        if rank3==0:
            #QSp = self.QSp[0].null()
            QSp = [self.QSp[0].null()]
        T3= self.__class__(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        T3.data[:]=0.0
        
        nidx3 = 0
        alpha = 1.0; beta=1.0

        self.num_of_contract = 0

        def calc(ind):
                self.num_of_contract += 1 
                #global T3
                p1, p2, p3, Dim1, Dim2, Dimc = tensor_player.the_tape[tensor_player.the_tape.calls][ind]
                #print Dim1, Dimc, self.data[p1:p1+Dim1*Dimc]
                #print tensor_player.the_tape[tensor_player.the_tape.calls][idx2, idx1]
                data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                #data3=iTensor.mul_temp(data1, data2, alpha, beta)
                data3=self.__class__.mul_temp(data1, data2, alpha, beta)
                T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
                print(ind, self.num_of_contract,  T3.data[p3:p3+Dim1*Dim2], data3)

        if 1:
            for ind in tensor_player.the_tape[tensor_player.the_tape.calls]:
                calc(ind)
        limit = 20
        #pprocess.pmap(calc, tensor_player.the_tape[tensor_player.the_tape.calls].keys(), limit=limit)
        print(T3.data)
        return T3

    def contract_core_player_parallel_2(self, T2, div, data=None, use_buf=False):
        """
        this one worked but much slower
        see iTensor_Contraction2 in f90
        把T1，和T2的非零block 如果量子数组合相等则收缩
        locals:
            div: num. of legs to be contracted for each tensor
            buffer: use buffer to save data of T3
        """
        if 0:
            from . import tensor_py
            iTensor = tensor_py.iTensor
        
        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        tQN = self.totQN+T2.totQN
        shift = rank1-div
        
        QSp = self.QSp[:shift]
        QSp.extend(T2.QSp[div:rank2])
        
        #T3= iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        if rank3==0:
            #QSp = self.QSp[0].null()
            QSp = [self.QSp[0].null()]
        T3= self.__class__(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        T3.data[:]=0.0
        
        nidx3 = 0
        alpha = 1.0; beta=1.0
        
        def calculate(ind):
            p1, p2, p3, Dim1, Dim2, Dimc = tensor_player.the_tape[tensor_player.the_tape.calls][ind]
            data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
            data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
            #data3=iTensor.mul_temp(data1, data2, alpha, beta)
            data3=self.__class__.mul_temp(data1, data2, alpha, beta) 
            return data3

        limit = 10
        results = pprocess.Map(limit=limit, reuse=1)

        # Wrap the calculate function and manage it.
        calc = results.manage(pprocess.MakeReusable(calculate))
        #calc = results.manage(pprocess.MakeParallel(calculate))
        ind_list = list(tensor_player.the_tape[tensor_player.the_tape.calls].keys())
        for ind in ind_list:
            calc(ind)
        
        for ind in ind_list:
            n = ind_list.index(ind)
            data3 = results[n]
            p1, p2, p3, Dim1, Dim2, Dimc = tensor_player.the_tape[tensor_player.the_tape.calls][ind]
            T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
        return T3

    def contract_core_player_parallel_3(self, T2, div, data=None, use_buf=False):
        """
        this one worked but much slower
        see iTensor_Contraction2 in f90
        把T1，和T2的非零block 如果量子数组合相等则收缩
        locals:
            div: num. of legs to be contracted for each tensor
            buffer: use buffer to save data of T3
        """
        if 0:
            from . import tensor_py
            iTensor = tensor_py.iTensor
        
        rank1 = self.rank
        rank2=T2.rank
        rank3=rank1+rank2-div-div
        tQN = self.totQN+T2.totQN
        shift = rank1-div
        
        QSp = self.QSp[:shift]
        QSp.extend(T2.QSp[div:rank2])
        
        #T3= iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        if rank3==0:
            #QSp = self.QSp[0].null()
            QSp = [self.QSp[0].null()]
        T3= self.__class__(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
        T3.data[:]=0.0
        
        nidx3 = 0
        alpha = 1.0; beta=1.0
        
        def calculate(ind):
            p1, p2, p3, Dim1, Dim2, Dimc = tensor_player.the_tape[tensor_player.the_tape.calls][ind]
            data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
            data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
            #data3=iTensor.mul_temp(data1, data2, alpha, beta)
            data3=self.__class__.mul_temp(data1, data2, alpha, beta) 
            return data3

        ind_list = list(tensor_player.the_tape[tensor_player.the_tape.calls].keys())
        limit = 5
        results = pprocess.pmap(calculate, ind_list, limit=limit)
        
        for ind in ind_list:
            n = ind_list.index(ind)
            data3 = results[n]
            p1, p2, p3, Dim1, Dim2, Dimc = tensor_player.the_tape[tensor_player.the_tape.calls][ind]
            T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
        return T3

    def group_legs_recorder(cls,itensor, ndiv=None):
        """
        此函数的目的是得到：
            量子数组合标记 idx 与 它在前后分块 标记p1, p2间的转换关系： 
                cls.QN_Group[0,idx] = gidx   # 在goup后的量子数中的编号
                cls.QN_Group[1,idx] = p1   #idx 在前一块中的编号
                cls.QN_Group[2,idx] = p2   #idx 在后一块中的编号            
                
                cls.QNG_Addr1[0,p1] = gidx #两个编号间的关系
                cls.QNG_Addr1[2,p1] = dim1  
        see iTensor_GroupLegs in f90
        """
        qsp = itensor.QSp[0].__class__
        cls.reset(QuantSpace=qsp)

        div = itensor.ndiv
        if ndiv is not None:  
            div = ndiv
        if div not in list(range(1, itensor.rank)):
            print("error, div exceeds range")
            print("div, rank", div, itensor.rank, itensor.ndiv)
            exit()
        rank = itensor.rank        

        dim1 = np.prod([itensor.QSp[i].nQN for i in range(div)])
        dim2 = np.prod([itensor.QSp[i].nQN for i in range(div, rank)])

        #if (.noself.allocated(QN_Group)) then
        #attention_omitted_something
        cls.QN_Group_Size = dim1*dim2  # 实际上它就是 self.idx_dim

        cls.QNG_Addr1 = np.ndarray((3, cls.QN_Group_Size), "int")
        cls.QNG_Addr2 = np.ndarray((3, cls.QN_Group_Size), "int")
        cls.QN_Group = np.ndarray((3, cls.QN_Group_Size), "int")
        
        #if (QN_Group_Size.lself.dim1*dim2) then
        #attention_omitted_something

        #if  not allocated[V_buf]:  
        #attention_omitted_something
        cls.V_size = 1024
        #cls.V_buff = np.ndarray((cls.V_size, cls.V_size), dtype=itensor.dtype)

        
        cls.QNG_Addr1[:, :]=-1
        cls.QNG_Addr2[:, :]=-1
        cls.QSp_Group1.nQN = 0
        cls.QSp_Group2.nQN = 0
        
        iQN = np.empty(itensor.MaxRank, "int")
        cls.group1_size = 0
        cls.group2_size = 0
        
        for idx in range(itensor.nidx):
            iQN[0:rank] = itensor.Addr_idx[0:rank, idx]
            #iQN[i]是leg i上的量子数计数器

            #这是前一半
            i = 0
            p1 = iQN[i]
            QN1 = itensor.QSp[i].QNs[iQN[i]]
            dim1 = itensor.QSp[i].Dims[iQN[i]]
            #do i = 2, div
            for i in range(1, div):
                p1 = p1*itensor.QSp[i].nQN+iQN[i]
                QN1 = QN1+itensor.QSp[i].QNs[iQN[i]]
                dim1 = dim1*itensor.QSp[i].Dims[iQN[i]]
            #p1: idx在前div个量子数组合中的线性坐标
           
            #这是后一半
            i = div 
            p2 = iQN[i]
            QN2 = itensor.QSp[i].QNs[iQN[i]]
            dim2 = itensor.QSp[i].Dims[iQN[i]]
            #do i = div+2, rank
            for i in range(div + 1, rank):
                p2 = p2*itensor.QSp[i].nQN+iQN[i]
                QN2 = QN2+itensor.QSp[i].QNs[iQN[i]]
                dim2 = dim2*itensor.QSp[i].Dims[iQN[i]]
            #p2 = p2

            gidx = cls.QSp_Group1.has_quant_num(QN1)
            #print "iii", 'idx', idx, 'p1', p1, 'p2', p2,'QN1', QN1, 'gidx', gidx, 'addr1', cls.QNG_Addr1[1, p1], 'addr2', cls.QNG_Addr1[1, p2]
            #newly_added

            if gidx < 0:  
                cls.QNG_Addr1[1,p1] = 0
                cls.QNG_Addr2[1,p2] = 0
                cls.QSp_Group1.add_to_quant_space(QN1, dim1)
                cls.QSp_Group2.add_to_quant_space(QN2, dim2)
                cls.group1_size += 1                 
                cls.group2_size += 1                                 
            else:
                if cls.QNG_Addr1[1,p1] < 0:  
                    cls.QNG_Addr1[1,p1] = cls.QSp_Group1.Dims[gidx]
                    cls.QSp_Group1.add_to_quant_space(QN1, dim1)
                    cls.group1_size += 1 
                if cls.QNG_Addr2[1,p2] < 0:  
                    cls.QNG_Addr2[1,p2] = cls.QSp_Group2.Dims[gidx]
                    cls.QSp_Group2.add_to_quant_space(QN2, dim2)
                    cls.group2_size += 1                     
            #attention_this_may_be_wrong
            #the following line is added by lzh
            if gidx<0: gidx = cls.QSp_Group1.has_quant_num(QN1) 

            cls.QN_Group[0,idx] = gidx   # 在goup后的量子数中的编号
            cls.QN_Group[1,idx] = p1   #idx 在前一块中的编号
            cls.QN_Group[2,idx] = p2   #idx 在后一块中的编号

            cls.QNG_Addr1[0,p1] = gidx #两个编号间的关系
            cls.QNG_Addr2[0,p2] = gidx

            cls.QNG_Addr1[2,p1] = dim1  #相当于之前的self.Block_idx[1,nidx] = d
            cls.QNG_Addr2[2,p2] = dim2
        
        tensor_player.the_tape[tensor_player.the_tape.calls] = ()
        
        temp = ("V_size", "V_buff", "QN_Group_Size", "QN_Group", 
                "QSp_Group1", "QSp_Group2", "QNG_Addr1", "QNG_Addr2")
        for t in temp:
            tensor_player.the_tape[tensor_player.the_tape.calls] += (getattr(cls, t),) 
            #tensor_player.the_tape[tensor_player.the_tape.calls] += (getattr(cls, t).copy(),) 
            warnings.warn("copy needed here?")

    def group_legs_player(cls,itensor, ndiv=None):
        temp = ("V_size", "V_buff", "QN_Group_Size", "QN_Group", 
                "QSp_Group1", "QSp_Group2", "QNG_Addr1", "QNG_Addr2")
        for i in range(len(temp)):
            #cls.__setattr__(temp[i], tensor_player.the_tape[tensor_player.the_tape.calls][i])
            setattr(cls, temp[i], tensor_player.the_tape[tensor_player.the_tape.calls][i])
    
    Qsp_copy_recorder = None
    
    def Qsp_copy_player(self, other=None):
        return self
        if other is None:
            return None
        else:
            self.nQN = other.nQN
            self.totDim=other.totDim
            self.RefQN=other.RefQN.copy()
            nqn=other.nQN
            self.QNs[:nqn]=[other.QNs[i].copy() for i in range(nqn)]
            #self.QNs=other.QNs.copy()
            self.Dims=other.Dims.copy()
            self.Addr=other.Addr.copy()



def tensor_player(which):
    """
        draw back  of this approach: debug is difficult
        this is the closure or factory function coding patten using nested funcitons, 
        as Mark Lutz said (on P.421) " 
            Although classes (described in Part VI of this book) are usually best at remembering state
            because they make it explicit with attribute assignments, such functions provide an
            alternative. ..... Moreover, function nesting is commonly used for decorators
            (explored in Chapter 38)—in some cases, it’s the most reasonable coding pattern"
    
    """
    #tensor_player.the_tape = TapeList[0]   #load default tape 
    
    def inner(func):
        recorder_dic = {'__init__':init_recorder, 
                        #'set_data_entrance':data_entrance_recorder, 
                        'contract_core':contract_core_recorder, 
                        'prepare_leg':prepare_leg_recorder,
                        'permutation':permute_recorder, 
                        #'group_legs':group_legs_recorder, 
                        #'Qsp_copy':Qsp_copy_recorder
                        }
        player_dic = {  
                        '__init__':init_player, 
                        #'set_data_entrance':data_entrance_player, 
                        'contract_core': contract_core_player, 
                        'permutation':  permute_player, 
                        'prepare_leg':prepare_leg_player,
                        #'group_legs':group_legs_player, 
                        #'Qsp_copy':Qsp_copy_player
                        }
        
        try:
            recorder = recorder_dic[which]
            player = player_dic[which]
        except:
            raise KeyError("wrong recorder/player name:%s"%which)
        
        def wrapper(*args, **kargs):
            """
                intro:
                    This defines the logics of the tensor palyer, it behaves exactly
                    as a normal walkman player. 
                
                tensor_player can have foure states:
                    ['play', 'record', 'stop', 'pause']
            
            """
            the_tape = tensor_player.the_tape
            if the_tape.STATE == "play":
                if the_tape.calls == the_tape.calls_tot: 
                    the_tape.calls = 0
                    the_tape.reach_tape_end = True
                    the_tape.tape_cleared = False

                #the_tape.calls += 1 
                return player(*args, **kargs)

            elif the_tape.STATE == "record":
                if the_tape.reach_tape_end: 
                    #print('tape of %s is full, its length is %d, clear tape before record.'%( func.__name__, len(the_tape.tape))) 
                    the_tape.tape = {}     
                    the_tape.calls = 0
                    the_tape.calls_tot = 0
                    the_tape.reach_tape_end = False
                
                if the_tape.PREV_STATE  == 'record':
                    if not the_tape.tape_cleared: 
                        the_tape.tape = {}     
                        the_tape.calls = 0
                        the_tape.calls_tot = 0
                        the_tape.reach_tape_end = False
                        the_tape.tape_cleared = True

                #the_tape.calls += 1 
                #the_tape.calls_tot += 1   #calls_tot only increases in record
                return recorder(func, *args, **kargs)
            
            elif the_tape.STATE == "stop": # return original method
                if the_tape.need_reset:  #issue: this actually has a problem 
                    the_tape.reset()
                    the_tape.need_reset = False #issue: this can only reset one method 
                    
                return func(*args, **kargs)
            elif the_tape.STATE == 'pause':
                return func(*args, **kargs)
            else:
                print('ttttttttttttt', type(the_tape.STATE))
                print(the_tape.STATE)
                raise ValueError
        return wrapper
    return inner


#tensor_player.STATE = 'stop'   # this line is needed, because sometiems tensor_player.STATE is used but tensor_player has not been called as decorator 
tensor_player.version = 'multiple'
tensor_player.the_tape = TapeList[0]

def set_player_state_auto(iter, record_at, tape_id=0, stop_at=10000000, verbose=False, info=0,  power_on=True):
    """
        arg verbose will be deprecated 
        issue: todo:  其实，这应该弄成context manager 
        
    """
    #global tensor_player
    the_tape = TapeList[tape_id]
    tensor_player.the_tape = the_tape
    the_tape.PREV_STATE = the_tape.STATE
    if power_on:
        if the_tape.NEXT_STATE is not None:  #NEXT_STATE 就是 手动指定，而非按照下面自动判断STATE, 
            the_tape.STATE = the_tape.NEXT_STATE
            the_tape.NEXT_STATE = None
        elif iter == record_at: 
            the_tape.STATE = "record"
            
        elif iter < record_at : #or iter>= stop_at:
            the_tape.STATE = "stop"
        else:
            the_tape.STATE = "play"
            #the_tape.calls= 0
        if info>-1:  #default display this msg 
            #if info>1:
            #    print('iiii', iter, the_tape.PREV_STATE)
            if the_tape.PREV_STATE != the_tape.STATE: 
                print("STATE of the_tape is changed from '%s' to '%s' at iter=%d"%(
                        the_tape.PREV_STATE, the_tape.STATE, iter))
            if the_tape.PREV_STATE == 'record' and the_tape.STATE == 'record':   # in very rare curcumstances, this may happen; add this line for robustness
                print('attention, PREV_STATE and current STATE of the_tape are both "record"')
    else:
        #the_tape.STATE = "stop"
        for i in TapeList.values():
            i.STATE = 'stop'
            i.reset()
        tensor_player.the_tape = TapeList[0]
        if info>0:
            print("the_tape is stopped")

def set_player_state_manual(state, tape_id=None, info=0):
    """
        params:
            state: in ['record', 'play', 'stop', 'pause']
    """
    tape_id = tape_id if tape_id is not None else 0
    TapeList[tape_id].STATE  = state
    if info>0:
        print('set state of tape %d to %s'%(tape_id, state))

def use_tensor_player(func):   #a decorator 
    #to implement 
    pass

def get_player_state(tape_id=None):
    tape_id = tape_id if tape_id is not None else 0
    return TapeList[tape_id].STATE

set_STATE_end_1 = set_player_state_auto 



def set_STATE_end_simple(iter,  record_at=0, resume=False, power_on=True):
    the_tape = TapeList[0]
    if power_on:
        if iter == record_at:
            the_tape.STATE = "record"
        else:
            the_tape.STATE = "play"
        print("\nset player to %s"%the_tape.STATE)
        #print "iter = %d"%iter, "set player to STATE= %s"%tensor_player.STATE
    else:
        the_tape.STATE = "stop"
        print("tensor_player is stopped")


class TestIt(unittest.TestCase): 
    
    @classmethod 
    def setUpClass(cls):
        from merapy.tensor_py import iTensor, qsp_any
        #the following line is required 
        iTensor=decorate_methods(decorator=tensor_player, meth_names=None)(iTensor)
        cls.iTensor = iTensor 
        #tape = Tape()
        #TapeList['ttt'] = tape 
        #cls.tape = tape 
        #print(TapeList[0])
    
    def setUp(self): 
        #TapeList[0].reset()
        #self.__class__.tape.reset()
        set_player_state_manual('stop', tape_id=0)
        #set_player_state_manual('stop', tape_id='ttt')
        tensor_player.STATE = 'stop'
    
    def tearDown(self):
        #self.__class__.tape.reset()
        set_player_state_manual('stop', tape_id=0)
        #set_player_state_manual('stop', tape_id='ttt')
        tensor_player.STATE = 'stop'
 
    def test_tensor_player(self): 
        from merapy import qsp_any
        #iTensor = TestIt.iTensor
        iTensor = self.__class__.iTensor
        
        q0 = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[4, 2, 2, 1, 1])
        q1 = q0.conj()
        qsp = [q0, q1]
        t = iTensor(QSp=qsp)
        np.random.seed(3333)
        t.data = 2*(np.random.random(t.size) -0.5)
        
        N = 5
        t0 = time.time()
        for i in range(N):
            print('i=', i)
            set_player_state_auto(iter=i, record_at=0, info=1)
            t=t.dot(t)
        t1 = time.time()
        print_vars(vars(),  ['t.norm()'])
        self.assertAlmostEqual(t.norm(), 0.38622236988000797, 10)
        print(tensor_player.the_tape.keys())
        set_player_state_manual('stop', tape_id=0)
        tensor_player.STATE = 'stop'
        tape = TapeList[0]
        print_vars(vars(),  ['tape.keys()'])
        tape.show()
            
        TapeList[0].reset()

    def test_tensor_player_gpu(self): 
        from merapy import qsp_any
        #iTensor = TestIt.iTensor
        iTensor = self.__class__.iTensor
        try:
            type(cp)
        except:
            print('cupy not installed, return')
            return 
        
        Dl = qsp_any('U1', qns=[0, 1, -1, ], 
                dims=[200, 130, 130, ])
        #Dl = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[4, 2, 2, 1, 1])
        Dr = Dl.conj()
        d = qsp_any('U1', qns=[1, -1], dims=[1, 1])
        qsp = [Dl, Dr, d]
        
        use_gpu = 1
        t = iTensor(QSp=qsp, use_gpu=use_gpu)
        tc = t.conj()
        
        N = 10
        for i in range(N):
            print_vars(vars(),  ['i'])
            set_player_state_auto(iter=i, record_at=0, tape_id=0,  info=0)
            res=t.contract(tc, [0, 1, 2], [1, 3, 2])
            #t.transpose((0, 2, 1))
        assert res.device == 'gpu'  

    def test_tensor_player_2(self):
        from merapy import qsp_any
        iTensor = self.__class__.iTensor
        #q0 = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[20, 10, 10, 5, 5])
        q0 = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[4, 2, 2, 1, 1])
        q1 = q0.conj()
        qsp = [q0, q0.copy(), q0.conj(), q0.conj()]
        t = iTensor(QSp=qsp)
        
        N = 10
        
        t0 = time.time()
        for i in range(N):
            #print('i=', i)
            t.transpose((1, 0, 3, 2))
            t.contract(t, [0, 1, 2, 3], [2, 3, 4, 5])
        t1 = time.time()
        
        
        for i in range(N):
            #print('i=', i)
            #set_player_state_auto(iter=i, record_at=1, info=1)    
            set_player_state_auto(iter=i, record_at=0, tape_id=0,  info=0)
            t.transpose((1, 0, 3, 2))
            t.contract(t, [0, 1, 2, 3], [2, 3, 4, 5])
            #t.contract(t, [0, 1], [1, 2])
        t2 = time.time()
            
        print_vars(vars(),  ['t1-t0'])
        print_vars(vars(),  ['t2-t1'])
        #tensor_player.STATE = 'stop'
        #print(TapeList[0])
        #print(tensor_player.the_tape.keys())
        
        set_player_state_manual('stop', tape_id=0)
            
           
        #tensor_player.STATE = 'stop'
    
    def test_tensor_player_performance(self):
        if 1:
            n, m = 150, 100
            a=cp.random.random((n, m))
            #a = a + a.T 
            ag = a.copy()
            N = 2

            t0 = time.time()
            for i in range(N):
                a.dot(a.T)
                #cp.linalg.eigh(a)
            t1 = time.time()
            print('time by GPU:', t1-t0)
        
        from merapy import qsp_any
        iTensor = self.__class__.iTensor

        #print_vars(vars(),  ['TapeList[0]'])
        #TapeList[0].reset()
        print(TapeList[0])
        

        # the following line is required 
        #iTensor=decorate_methods(decorator=tensor_player, meth_names=None)(iTensor)
        #Dl = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[400, 100, 100, 50, 50])
        Dl = qsp_any('U1', qns=[0, 1, -1, ], 
                dims=[500, 130, 130, ])
        #Dl = qsp_any('U1', qns=[0, 1, -1, 2, -2], dims=[4, 2, 2, 1, 1])
        Dr = Dl.conj()
        d = qsp_any('U1', qns=[1, -1], dims=[1, 1])
        qsp = [Dl, Dr]
        
        use_gpu = 1
        
        t = iTensor(QSp=qsp, use_gpu=use_gpu)
        tc = t.copy()
        N = 10
        
        data = cp.ndarray(t.data.size) if t.use_gpu ==1 else np.ndarray(t.data.size) 
        
        t1 = time.time()
        for i in range(N):
            set_player_state_auto(iter=i, record_at=0, tape_id=0,  info=0)
            t.dot(tc, data=data.data)
            #t.dot(tc)
            #t.dot(tc, use_buf=1)
        t2 = time.time()
        
        
        
        print_vars(vars(),  ['t.data.size'])
        print_vars(vars(),  ['t.use_gpu'])
            
        print_vars(vars(),  ['t2-t1'])
        print(TapeList[0].calls)
      
        tensor_player.STATE = 'stop'
        print(tensor_player.the_tape.keys())
        set_player_state_manual('stop', tape_id=0)
            
        TapeList[0].reset()
           
        #tensor_player.STATE = 'stop'
    
    def test_temp(self): 
        print(dir(cublas))
        
        dgemm  =  cublas.gemm
        a = np.ndarray((3, 3))
        b = np.ndarray((3, 3))
        c = np.ndarray((3, 3))
        a = cp.asarray(a)
        b = cp.asarray(b)
        c = cp.asarray(c)
        #dgemm(a)
        print_vars(vars(),  ['dgemm.__doc__'])
        raise  
        
        from merapy.tensor_py import iTensor 
        if 1:
            print(iTensor.contract_core)
            print(iTensor.__init__)
        
        #iTensor.reset_player()
        print('iiiiiiiiiii', id(tensor_player))
        for i in range(1, 10):
            print('i=', i)
            set_player_state_auto(iter=i, record_at=1, info=1)    
            t1 = iTensor.example(rank=4)
            t2 = iTensor.example(rank=4)
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3 = t1.contract(t2, [0, 1, 2, 3], [4, 2, 5, 6])
            t3.permutation([0, 2, 3, 1, 4, 5])
            t3.permutation([0, 2, 3, 1, 4, 5])
        #tensor_player.STATE = 'stop'
        print(TapeList[0])
        print(tensor_player.the_tape.keys())

        #status = iTensor.get_player_status()
        #print_vars(vars(),  ['status'])


if __name__ == "__main__":


    warnings.filterwarnings('ignore')
    if 0:
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
        #'test_tensor_player_performance', 
        'test_tensor_player', 
        #'test_tensor_player_gpu', 
        #'test_tensor_player_2', 
        #'test_temp', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))

        unittest.TextTestRunner(verbosity=0).run(suite)

        

