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
import os, sys
import unittest 
import pickle as pickle
import cProfile 
import numpy as np
import warnings

#from merapy import common_util
import merapy.common_util as common_util

from merapy import array_permutation
from merapy.context_util import redirect



"""
    design goal: 
        1. use new player   
        2. State Information Retention  in tape 
        3. not using many @ statement       ---done
"""


def decorate_methods(decorator, meth_names):
    meth_names_dic = {
            "iTensor":[
                "__init__", 
                "set_data_entrance", 
                "contract_core", 
                "permutation"
                ], 
            "QuantSpaceBase":["copy"], 
            #"Tensor_svd":["group_legs"]   directly use tensor_player decorated see tensor_svd.py
            }
    meth_original_bac = {}
    
    #print 'ccccccc'
    #meth_names_dic["iTensor"].remove('contract_core')

    meth_names_dic["iTensor"].pop(0) 
    meth_names_dic["QuantSpaceBase"].pop(0)

    def DecoDecorate(Class):
        class_name = Class.__name__
        meth_names = meth_names_dic[class_name]
        #if __debug__:   #this is annoying
        #    print "%s are decorated in %s"%(meth_names, class_name)
        for attr in meth_names:
            meth_original_bac[class_name] = Class.__dict__[attr] #backup original method
            deced_meth = decorator(which=attr)(Class.__dict__[attr])
            
            if isinstance(Class.__dict__[attr], classmethod): #spcial treatment for classmethod
                deced_meth = classmethod(deced_meth)
            setattr(Class, attr, deced_meth) 
            #print "%s is decorated in %s. originally %s; now %s"%(
            #        attr, class_name, meth_original_bac[class_name], deced_meth)

        return Class
    
    return DecoDecorate

#permute_tape = np.ndarray((700000, 40), np.int32)

def tensor_player(which):
    """
        draw back  of this approach: debug is difficult
        this is the closure or factory function coding patten using nested funcitons, 
        as Mark Lutz said (on P.421) " 
            Although classes (described in Part VI of this book) are usually best at remembering state
            because they make it explicit with attribute assignments, such functions provide an
            alternative. ..... Moreover, function nesting is commonly used for decorators
            (explored in Chapter 38)—in some cases, it’s the most reasonable coding pattern"

        which = "init", "contract", "permute"
    """
    mapper = {"__init__":"init", 
            "set_data_entrance":"data_entrance", 
            "permutation":"permute", 
            "contract_core":"contract", 
            "copy":"Qsp_copy"
            }
    if which in mapper:
        which = mapper[which]
    
    tensor_player.STATE = "stop"   #default not using player
    tensor_player.PREV_STATE = 'stop'
    tensor_player.NEXT_STATE = None
    tensor_player.reset = False 
    
    tensor_player.save_tape = False
    tensor_player.load_tape = False
    tensor_player.tape_prefix = ""
    
    def inner(func):
        """
            this trick works only in that every inner returned by tensor_player
            are distinct objs       
            note:
                1.  I have used  'inner.tape' to make 'tape' be able to used in wrapper, 
                    another way to realize this is to use 'nonlocal' statement in wrapper
                2.  inner.attr can be inspected in the following way:
                    #t is an iTensor instance 
                    inner = t.permutation.__closure__[1]
                    print( inner.cell_contents.calls_tot)
        
        """
        
        def reset():
            inner.tape = {}   # blank tape 
            inner.calls = 0   # record/play at #song 
            inner.calls_tot = 0   #total num of 'songs' recorded 
            inner.reach_tape_end = False 
            inner.tape_cleared = False
        
        reset()
        
        inner.tape_saved = False
        inner.tape_loaded = False

        if 1:   #define recorder and player
            def init_recorder(self, rank,  QSp, totQN, order="F",  
                    buffer=None, use_buf=False, index_data=True, has_data=True, shallow=False):
                """
                    这种做法有些过于激进, 更安全的是下面的 data_entrance_recorder/player 
                
                """
                func(self, rank,  QSp, totQN, order="F",  
                    buffer=buffer, use_buf=use_buf, index_data=index_data, has_data=has_data, shallow=shallow)
                
                struct_dict = self.__dict__.copy()
                if has_data:
                    struct_dict.pop("data")
                struct_dict.pop("buf_ref")
                inner.tape[inner.calls] = struct_dict
            
            #@profile
            def init_player(self, rank,  QSp, totQN, order="F",  
                    buffer=None, use_buf=False, index_data=True, has_data=True, shallow=False):
    
                #following is wrong,  because change self wold affect inner.tape,  e.g. append data attr
                #self.__dict__ = inner.tape[inner.calls]     
                
                #following is not deepcopy, self and tape share same value, but NOT share key
                self.__dict__ = inner.tape[inner.calls].copy()
                
                self.buf_ref = np.array([-1, -1], np.int)
                
                #print 'ssss', self.dtype 
                if has_data:
                    if use_buf:   #use internal T_BUFFER; else use external buffer or no buffer
                        buffer = self.buffer_assign(data_size=self.totDim)
                    self.data = np.ndarray(self.totDim, buffer=buffer, dtype=self.dtype, order="C")

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
                inner.tape[inner.calls] = ( self.Addr_idx.copy(), self.Block_idx.copy(), 
                        self.idx.copy(), self.nidx, self.idx_dim, self.totDim)

            def data_entrance_player(self, order="F"):
                ( self.Addr_idx, self.Block_idx, self.idx, self.nidx, 
                        self.idx_dim, self.totDim)=inner.tape[inner.calls]

            def permute_recorder(self, P, buffer=None, use_buf=False):
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
                Tp=self.__class__(rank, QSp, totQN, buffer=buffer, dtype=self.dtype, use_buf=use_buf)
                pos=np.empty(self.rank,"int")
                Dims=np.empty(self.rank,"int")

                tape_ind = np.ndarray((self.nidx, 3), np.int32)
                tape_dim = np.ndarray((self.nidx, 32), np.int32)
                tape_ord = np.ndarray((self.nidx, 32), np.int32)
                
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
                    
                    temp= self.data[pidx:pidx+totDim]   #.copy()
                    
                    #data = Tp.data[qidx:qidx+totDim]
                    #array_permutation.array_permutation_inplace(temp, self.rank, Dims, P, data)
                    Tp.data[qidx:qidx+totDim]=array_permutation.array_permutation(temp,rank,Dims,P)
                    
                    #record[n] = (qidx, pidx, totDim, Dims.copy())
                    tape_ind[n][0:3] = [pidx, qidx, totDim]
                    tape_dim[n][:rank] = Dims[:rank]
                    tape_ord[n][:rank] = P[:rank]
                if 0:
                    #I may implement this later, sort the tape according to totDim
                    print("ssssort")
                    print(self.nidx)
                    tape_sort = tape_ind[:, 2].argsort()
                    #tape_sort.sort()
                    print(tape_ind[tape_sort[-1::-1]])
                    inner.tape[inner.calls] = (self.nidx, tape_ind[tape_sort], tape_dim[tape_sort], tape_ord[tape_sort])
                    #inner.tape[inner.calls] = (self.nidx, tape_ind[tape_sort[-1::-1]],tape_dim[tape_sort[-1::-1]], tape_ord[tape_sort[-1::-1]])
                else:
                    inner.tape[inner.calls] = (self.nidx, tape_ind, tape_dim, tape_ord)
                return Tp

            def permute_player_bac(self, P, buffer=None, use_buf=False):
                
                rank = self.rank
                #permute QSp, QNs
                if 0:
                    QSp=[self.QSp[P[i]].copy() for i in range(rank)]
                    totQN = self.totQN.copy()
                else:
                    QSp=[self.QSp[P[i]] for i in range(rank)]
                    totQN = self.totQN
                
                #Tp=self.__class__(rank, QSp, totQN,  buffer=buffer, use_buf=use_buf)
                Tp=self.__class__(rank, QSp, totQN, dtype=self.dtype, buffer=buffer, use_buf=use_buf)
                
                nidx, tape_ind, tape_dim, tape_ord = inner.tape[inner.calls]
                if self.data.dtype == float: 
                    #print('ppppp', self.data.size, Tp.data.size)
                    array_permutation.permute_player_fort(self.rank, 
                                tape_ind, tape_dim, tape_ord, self.data, Tp.data, nidx, Tp.data.size)
                else:  #complex 
                    array_permutation.complex_permute_player_fort(self.rank, 
                                tape_ind, tape_dim, tape_ord, self.data, Tp.data, nidx, Tp.data.size)

                return Tp

            #@profile
            def permute_player(self, P, buffer=None, use_buf=False):
                #since __init__ has been decorated, just do in following way                
                Tp=self.__class__(rank=None, QSp=None, totQN=None, buffer=buffer, use_buf=use_buf)
                
                nidx, tape_ind, tape_dim, tape_ord = inner.tape[inner.calls]
                array_permutation.permute_player_fort(self.rank, 
                            tape_ind, tape_dim, tape_ord, self.data, Tp.data, nidx, Tp.data.size)

                return Tp
            
            def contract_core_recorder(self, T2, div, data=None, use_buf=False):
                """
                    see iTensor_Contraction2 in f90
                    把T1，和T2的非零block 如果量子数组合相等则收缩
                    locals:
                        div: num. of legs to be contracted for each tensor
                        buffer: use buffer to save data of T3
                """
                
                #contract_record_1 = {}
                contract_record_1 = np.ndarray((self.nidx*T2.nidx, 6), np.int)
                ind_count = 0
                
                rank1 = self.rank
                rank2=T2.rank
                rank3=rank1+rank2-div-div
                tQN = self.totQN+T2.totQN
                shift = rank1-div
                if 1:
                    QSp = [self.QSp[i].copy() for i in range(shift)]
                    QSp.extend([T2.QSp[i].copy() for i in range(div, rank2)])
                else:
                    QSp = self.QSp[:shift]
                    QSp.extend(T2.QSp[div:rank2])
                
                #T3= iTensor(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf)
                if rank3==0:
                    #QSp = self.QSp[0].null()
                    #QSp = [self.QSp[0].null()]
                    QSp = []

                dtype = complex if self.dtype == complex or T2.dtype == complex else float 
                
                T3= self.__class__(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, use_buf=use_buf, dtype=dtype)
                T3.data[:]=0.0
                
                nidx3 = 0
                alpha = 1.0; beta=1.0
                iQN1=np.empty(self.rank + 1, np.int)
                iQN2=np.empty(T2.rank + 1, np.int)        
                iQN3=np.empty(T3.rank + 1, np.int) # +1 to avoid T3.rank=0

                for idx2 in range(T2.nidx):
                    iQN2[0] = 0  #!for rank=0
                    iQN2[0:rank2]=T2.Addr_idx[0:rank2,idx2]
                    p2 = T2.Block_idx[0,idx2]
                    
                    Dim2 = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div, rank2)], dtype=np.int)
                    Dimc = np.prod([T2.QSp[i].Dims[iQN2[i]] for i in range(div)], dtype=np.int)
                    
                    for idx1 in range(self.nidx):
                        iQN1[0] = 1 #!for rank=0
                        iQN1[0:rank1]=self.Addr_idx[0:rank1,idx1]
                        p1 = self.Block_idx[0, idx1]
                        iseq = np.all(iQN1[shift:shift+div] == iQN2[0:div])
                        if not iseq:
                            #如果量子数组合相等则收缩
                            continue
                        Dim1 = np.prod([self.QSp[i].Dims[iQN1[i]] for i in range(shift)], dtype=np.int)
                        
                        iQN3[0] = 1 #!for rank=1
                        iQN3[0:shift] = iQN1[0:shift]
                        iQN3[shift:rank3] = iQN2[div:rank2]
                        p3=T3.get_position(iQN3[:T3.rank])
                        idx3 = T3.idx[p3]
                        p3 = T3.Block_idx[0,idx3]
                        
                        contract_record_1[ind_count][:] = (p1, p2, p3, Dim1, Dim2, Dimc)
                        ind_count  += 1  
                        #contract_record_1[idx2, idx1] = np.array((p1, p2, p3, Dim1, Dim2, Dimc), np.int)
                        #contract_record_1[idx2, idx1] = (p1, p2, p3, Dim1, Dim2, Dimc)
                

                        data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                        data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F') 
                        data3=self.__class__.mul_temp(data1, data2, alpha, beta, dtype=dtype)
                        T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
                
                inner.tape[inner.calls] = contract_record_1[:ind_count, :], ind_count
                return T3

            def contract_core_player_bac(self, T2, div, data=None, use_buf=False):
                """
                    see iTensor_Contraction2 in f90
                    把T1，和T2的非零block 如果量子数组合相等则收缩
                    locals:
                        div: num. of legs to be contracted for each tensor
                        buffer: use buffer to save data of T3
                """
                rank1 = self.rank
                rank2 = T2.rank
                rank3 = rank1+rank2-div-div
                tQN = self.totQN+T2.totQN
                shift = rank1-div
                
                QSp = self.QSp[:shift]
                QSp.extend(T2.QSp[div:rank2])
                
                if rank3==0:
                    QSp = [self.QSp[0].null()]
                    
                dtype = complex if self.dtype == complex or T2.dtype == complex else float 
                T3 = self.__class__(rank=rank3, QSp=QSp, totQN=tQN, buffer=data, dtype=dtype, use_buf=use_buf)
                rec, num_rec = inner.tape[inner.calls]
                if dtype == float:  
                    common_util.contract_core_player_fort(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                else: 
                    common_util.contract_core_player_fort_complex(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                    #common_util.contract_core_player_fort(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                    
                return T3
            
            #@profile
            def contract_core_player(self, T2, div, data=None, use_buf=False):
                """
                    see iTensor_Contraction2 in f90
                    把T1，和T2的非零block 如果量子数组合相等则收缩
                    locals:
                        div: num. of legs to be contracted for each tensor
                        buffer: use buffer to save data of T3
                """
                T3 = self.__class__(rank=None, QSp=None, totQN=None, buffer=data, use_buf=use_buf)
                
                rec, num_rec = inner.tape[inner.calls]
                common_util.contract_core_player_fort(self.data, T2.data, T3.data, rec, num_rec=num_rec)
                return T3
            
            def contract_core_player_parallel_1(self, T2, div, data=None, use_buf=False):
                """
                not work
                a parallel version
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

                self.num_of_contract = 0

                def calc(ind):
                        self.num_of_contract += 1 
                        #global T3
                        p1, p2, p3, Dim1, Dim2, Dimc = inner.tape[inner.calls][ind]
                        #print Dim1, Dimc, self.data[p1:p1+Dim1*Dimc]
                        #print inner.tape[inner.calls][idx2, idx1]
                        data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                        data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                        #data3=iTensor.mul_temp(data1, data2, alpha, beta)
                        data3=self.__class__.mul_temp(data1, data2, alpha, beta)
                        T3.data[p3:p3+Dim1*Dim2]  += data3.ravel('F')[:]   #attention_here  fortran order
                        print(ind, self.num_of_contract,  T3.data[p3:p3+Dim1*Dim2], data3)

                if 1:
                    for ind in inner.tape[inner.calls]:
                        calc(ind)
                limit = 20
                #pprocess.pmap(calc, inner.tape[inner.calls].keys(), limit=limit)
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
                    p1, p2, p3, Dim1, Dim2, Dimc = inner.tape[inner.calls][ind]
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
                ind_list = list(inner.tape[inner.calls].keys())
                for ind in ind_list:
                    calc(ind)
                
                for ind in ind_list:
                    n = ind_list.index(ind)
                    data3 = results[n]
                    p1, p2, p3, Dim1, Dim2, Dimc = inner.tape[inner.calls][ind]
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
                    p1, p2, p3, Dim1, Dim2, Dimc = inner.tape[inner.calls][ind]
                    data1=self.data[p1:p1+Dim1*Dimc].reshape((Dim1,Dimc), order='F')    #attention_here fortran order
                    data2=T2.data[p2:p2+Dim2*Dimc].reshape((Dimc,Dim2), order='F')    
                    #data3=iTensor.mul_temp(data1, data2, alpha, beta)
                    data3=self.__class__.mul_temp(data1, data2, alpha, beta) 
                    return data3

                ind_list = list(inner.tape[inner.calls].keys())
                limit = 5
                results = pprocess.pmap(calculate, ind_list, limit=limit)
                
                for ind in ind_list:
                    n = ind_list.index(ind)
                    data3 = results[n]
                    p1, p2, p3, Dim1, Dim2, Dimc = inner.tape[inner.calls][ind]
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
                
                inner.tape[inner.calls] = ()
                
                temp = ("V_size", "V_buff", "QN_Group_Size", "QN_Group", 
                        "QSp_Group1", "QSp_Group2", "QNG_Addr1", "QNG_Addr2")
                for t in temp:
                    inner.tape[inner.calls] += (getattr(cls, t),) 
                    #inner.tape[inner.calls] += (getattr(cls, t).copy(),) 
                    warnings.warn("copy needed here?")

            def group_legs_player(cls,itensor, ndiv=None):
                temp = ("V_size", "V_buff", "QN_Group_Size", "QN_Group", 
                        "QSp_Group1", "QSp_Group2", "QNG_Addr1", "QNG_Addr2")
                for i in range(len(temp)):
                    #cls.__setattr__(temp[i], inner.tape[inner.calls][i])
                    setattr(cls, temp[i], inner.tape[inner.calls][i])
            
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
        
            recorder_dic = {"init":init_recorder, 
                            "data_entrance":data_entrance_recorder, 
                            "contract":contract_core_recorder, 
                            "permute":permute_recorder, 
                            "group_legs":group_legs_recorder, 
                            "Qsp_copy":Qsp_copy_recorder
                            }
            player_dic = {  "data_entrance":data_entrance_player, 
                            "contract": contract_core_player_bac, #contract_core_player_bac, 
                            "permute":  permute_player_bac, #permute_player_bac, 
                            "group_legs":group_legs_player, 
                            "init":init_player, 
                            "Qsp_copy":Qsp_copy_player
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
                    ['play', 'record', 'stop']
                    I may add a 'pause'
            
            """
            #print "STATEEE", tensor_player.STATE, func.func_name , 
            #print "STATEEE", tensor_player.STATE #, tensor_player.__code__#, func.func_name 
            #print('wraperrrrrrrrrrrrr', id(tensor_player))
            #print('ffffffffff', func.__name__, inner.calls, inner.calls_tot)
            
            if tensor_player.STATE == "play":
                if inner.calls == 0:
                    if tensor_player.load_tape:
                        if not inner.tape_loaded:
                            fn = tensor_player.tape_prefix + "__" +  func.__name__ + ".tape"
                            inn = open(fn, "rb")
                            tape = pickle.load(inn)
                            inn.close()
                            inner.tape = tape
                            inner.calls_tot = max(tape.keys())
                            inner.reach_tape_end = True
                            print("%s loaded successfully. directly play afterwards."%fn)
                            
                            inner.tape_loaded = True  #only need to load once

                if inner.calls == inner.calls_tot: 
                    if tensor_player.save_tape:  #this is when just after record
                        if not inner.tape_saved:
                            fn = tensor_player.tape_prefix + "__" +  func.__name__ + ".tape"
                            out = open(fn, "wb")
                            pickle.dump(inner.tape, out)
                            out.close()
                            print("%s saved successfully."%fn)
                            inner.tape_saved = True
                    inner.calls = 0
                    inner.reach_tape_end = True
                    inner.tape_cleared = False

                inner.calls += 1 
                return player(*args, **kargs)

            elif tensor_player.STATE == "record":
                if inner.reach_tape_end: 
                    #print('tape of %s is full, its length is %d, clear tape before record.'%( func.__name__, len(inner.tape))) 
                    inner.tape = {}     
                    inner.calls = 0
                    inner.calls_tot = 0
                    inner.reach_tape_end = False
                    
                #print('ppppppppp', tensor_player.PREV_STATE)
                if tensor_player.PREV_STATE  == 'record':
                    if not inner.tape_cleared: 
                        inner.tape = {}     
                        inner.calls = 0
                        inner.calls_tot = 0
                        inner.reach_tape_end = False
                        inner.tape_cleared = True

                inner.calls += 1 
                inner.calls_tot += 1   #calls_tot only increases at here 
                
                return recorder(*args, **kargs)
            
            elif tensor_player.STATE == "stop": # return original method
                #inner.tape_cleared = False
                if tensor_player.reset:  #issue: this actually has a problem 
                    reset()
                    
                    inner.tape_saved = False
                    inner.tape_loaded = False
                    tensor_player.load_tape = False
                    tensor_player.save_tape = False
                    
                    tensor_player.reset = False #issue: this can only reset one method 
                    
                return func(*args, **kargs)

        return wrapper
    return inner


tensor_player.version = 'single'
tensor_player.STATE = 'stop'   # this line is needed, because sometiems tensor_player.STATE is used but tensor_player has not been called as decorator 
#tensor_player.NEXT_STATE = None 


def set_STATE_end_simple(iter,  record_at=0, resume=False, power_on=True):
    if power_on:
        if iter == record_at:
            tensor_player.STATE = "record"
        else:
            tensor_player.STATE = "play"
        print("\nset player to %s"%tensor_player.STATE)
        #print "iter = %d"%iter, "set player to STATE= %s"%tensor_player.STATE
    else:
        tensor_player.STATE = "stop"
        print("tensor_player is stopped")

#issue: todo:  其实，这应该弄成context manager 
def set_player_state_auto(iter, record_at, stop_at=10000000, tape_id=None, verbose=False, info=0,  power_on=True):
    """
        arg verbose will be deprecated 
        
    """
    tensor_player.PREV_STATE = tensor_player.STATE
    #print('setttttttttttt', id(tensor_player))
    if power_on:
        if tensor_player.NEXT_STATE is not None:  #NEXT_STATE 就是 手动指定，而非按照下面自动判断STATE, 
            tensor_player.STATE = tensor_player.NEXT_STATE
            tensor_player.NEXT_STATE = None
        elif iter == record_at: 
            tensor_player.STATE = "record"
            if tensor_player.load_tape:
                path = tensor_player.tape_prefix + '__set_data_entrance.tape'
                msg = 'test existence of tape %s...'%path
                if os.path.exists(path): 
                    tensor_player.STATE = "play"
                    tensor_player.save_tape = False
                    msg  += 'tape found. set tensor_player.save_tape=Flase'
                else: 
                    msg  += ' but not found. set tensor_player.load_tape=False'
                    tensor_player.load_tape = False
                print(msg) 
            
        elif iter < record_at : #or iter>= stop_at:
            tensor_player.STATE = "stop"
        else:
            tensor_player.STATE = "play"
        if info>-1:  #default display this msg 
            #if info>1:
            #    print('iiii', iter, tensor_player.PREV_STATE)
            if tensor_player.PREV_STATE != tensor_player.STATE: 
                print("STATE of tensor_player is changed from '%s' to '%s' at iter=%d"%(
                        tensor_player.PREV_STATE, tensor_player.STATE, iter))
            if tensor_player.PREV_STATE == 'record' and tensor_player.STATE == 'record':   # in very rare curcumstances, this may happen; add this line for robustness
                print('attention, PREV_STATE and current STATE of tensor_player are both "record"')
        
    else:
        tensor_player.STATE = "stop"
        if info>0:
            print("tensor_player is stopped")

def set_player_state_manual(state, tape_id=None, info=0):
    """
        tape_id is not used
    """
    tensor_player.STATE = state
    if info>0:
        print('set tensor_player.STATE to %s'%(state, ))

def get_player_state(tape_id=None):
    """
        tape_id is not used
    """
    return tensor_player.STATE 


set_STATE_end_1 = set_player_state_auto 

def change_tape(t, state):
    raise NotImplemented("not completed")
    state_bac = tensor_player.STATE
    meth_names= ["set_data_entrance", "contract_core", "permutation"]
    iTensor = t.__class__   # t is just used for passing one hook in
    meth_bac = {}
    for attr in meth_names:
        meth_bac[attr] = iTensor.__dict__[attr] #backup original method
        meth_deced = tensor_player(which=attr)(iTensor.__dict__[attr])
        setattr(iTensor, attr, meth_deced) 
    return meth_bac, satte_bac

def restore_tape(t, meth_bac, state_bac):
    raise NotImplemented("not completed")
    if 1:
        iTensor = t.__class__
        for i in meth_bac:
            setattr(iTensor, i, meth_bac[i])
        tensor_player.STATE = state_bac


