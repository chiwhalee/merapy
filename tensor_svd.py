#coding=utf8
#include "header.f90"

import unittest 
import numpy as np
import numpy.linalg as linalg
import scipy 
from scipy.sparse.linalg import eigs


from merapy.tensor import iTensor
#from tensor_py import tensor_player
from merapy.decorators import *
from merapy.quantum_number import * #import QuantSpace
from merapy.utilities import print_vars 
import merapy.common_util as common_util 
from merapy.utilities import print_vars 

""" 
issue225: 
    error occor when use ifort compiled matrix_svd_unitilize and matrix size is too large.
    I use scipy.linalg.svd instead of matrix_svd_unitilize. 
    but I am not sure whether scipy.linalg.svd is as fast as fortran code
    and there would be NOT converge error
"""


print "call dsyev what is lwork? in common4py.f90 this value may be improper"

class Tensor_svd(object):
    """
    this implementation without inheriete from iTensor
    """
    @classmethod
    def reset(cls, QuantSpace):
        """
        QuantSpace is a class !
        
        """
        cls.V_size = 0 
        cls.V_buff = None
        cls.QN_Group_Size = 0 
        cls.QN_Group = None
        #cls.QSp_Group1 = QuantSpace(n=0, qns=[], dims=[])
        #cls.QSp_Group2 = QuantSpace(n=0, qns=[], dims=[])
        cls.QSp_Group1 = QuantSpace.empty()#()
        cls.QSp_Group2 = QuantSpace.empty() #()

        cls.QNG_Addr1 = None
        cls.QNG_Addr2 = None

    @classmethod
    @tensor_player(which="group_legs")
    def group_legs(cls,itensor, ndiv=None):
        """
            reshape of iTensor, i.e. group its QSps
            其中的要点是，要知道在原来的张量中的某一个data block 在reshape后的张量的地址在哪、长度是多少
            把量子数组合分成两段 (q_1, q_2, ..., q_k;  q_{k+1}, ..., q_{rank}). 
            第1段对应于下指标，第二段是上指标
            此函数的目的是得到：
                量子数组合标记 idx 与 它在前后两段 标记p1, p2间的转换关系： 
                    gidx: 
                        group后的gidx, for Travial in [0]; for Z2 in [0,1]; 
                        for U1 in [0,1,2]
                    cls.QN_Group[0,idx] = gidx   # 在goup后的量子数中的编号
                    cls.QN_Group[1,idx] = p1   #idx 在前一块中的编号
                    cls.QN_Group[2,idx] = p2   #idx 在后一块中的编号            
                    
                    cls.QNG_Addr1[0,p1] = gidx #两个编号间的关系
                    cls.QNG_Addr1[2,p1] = dim1  
                    QSp_Group1, QSp_Group2: 将指标QSp merge成1，2

            examples: 
                group_legs(u, div=3), where u is a rank 4 ,Z2 tensor
                    self.Addr_idx:
                    [[0 1 1 0 1 0 0 1]
                     [0 1 0 1 0 1 0 1]
                     [0 0 1 1 0 0 1 1]
                                            #reshape 即在这里把量子数组合切成两段
                     [0 0 0 0 1 1 1 1]]
                    self.Block_idx:
                    [[ 0  1  2  3  4  5  6  7]
                     [ 1  1  1  1  1  1  1  1]
                     [ 0  3  5  6  9 10 12 15]]

                    QN_Group:               #see self.Addr_idx
                    [[0 0 0 0 1 1 1 1] 
                     [0 6 5 3 4 2 1 7]
                     [0 0 0 0 1 1 1 1]]
                    QNG_Addr1:
                    [[0 1 1 0 1 0 0 1]
                     [0 2 1 3 0 2 1 3]
                     [1 1 1 1 1 1 1 1]]     #len of data block
                    QNG_Addr2:
                    [[0 1]
                     [0 0]
                     [1 1]]                 #len of data block
                    
                    QSp_Group1:
                        class:	QspZ2
                    totDim:	8
                    Dims:	[4 4]

                    QSp_Group2:
                        class:	QspZ2
                    totDim:	2
                    Dims:	[1 1]
                    
                group_legs(u, div=3), where u is a rank 3 ,U1 tensor            
                    self.Addr_idx:
                    [[0 2 1 2 0 1 0]
                     [0 1 2 0 2 0 1]
                     [0 0 0 1 1 2 2]]
                    self.Block_idx:
                    [[ 0  8 10 12 14 16 18]
                     [ 8  2  2  2  2  2  2]
                     [ 0  5  7 11 15 19 21]]

                    QN_Group:
                    [[0 0 0 1 1 2 2]
                     [0 7 5 6 2 3 1]
                     [0 0 0 1 1 2 2]]
                    QNG_Addr1:
                    [[ 0  2  1  2 -1  0  1]
                     [ 0  2  2  0 -1  5  0]
                     [ 4  2  2  2 -1  1  2]]
                    QNG_Addr2:
                    [[0 1 2]
                     [0 0 0]
                     [2 1 1]]
                    
                    QSp_Group1:
                        class:	QspU1
                    totDim:	14
                    Dims:	[6 4 4]
                    
                    QSp_Group2:
                        class:	QspU1
                    totDim:	4
                    Dims:	[2 1 1]

        """
        qsp = itensor.QSp[0].__class__ #() #.empty()
        cls.reset(QuantSpace=qsp)

        div = itensor.ndiv
        if ndiv is not None:  
            div = ndiv
        assert 0<div<itensor.rank, "error, div exceeds range, div=%s, itensor.rank=%s, itensor.ndiv=%s"%(
                div, itensor.rank, itensor.ndiv)
        #if div not in range(1, itensor.rank):
        #    raise ValueError("error, div exceeds range, div=%s, itensor.rank=%s, itensor.ndiv=%s"
        #            %(div, itensor.rank, itensor.ndiv))

        rank = itensor.rank        

        dim1 = np.prod([itensor.QSp[i].nQN for i in range(div)])
        dim2 = np.prod([itensor.QSp[i].nQN for i in range(div, rank)])
        cls.QN_Group_Size = dim1*dim2  # 实际上它就是 self.idx_dim

        cls.QNG_Addr1 = np.ndarray((3, cls.QN_Group_Size), "int")
        cls.QNG_Addr2 = np.ndarray((3, cls.QN_Group_Size), "int")
        cls.QN_Group = np.ndarray((3, cls.QN_Group_Size), "int")
        

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
    
    @classmethod
    def show_group(cls, itensor):
        """
        newly_added
        """
        print 'self.Addr_idx:\n', itensor.Addr_idx[:,:itensor.nidx]
        print 'self.Block_idx:\n', itensor.Block_idx[:, :itensor.nidx]
        print "\nQN_Group:\n", cls.QN_Group[:, 0:itensor.nidx]
        print "QNG_Addr1:\n", cls.QNG_Addr1[:, 0:cls.group1_size]
        print "QNG_Addr2:\n", cls.QNG_Addr2[:, 0:cls.group2_size]
        print 'QSp_Group1:\n\t', cls.QSp_Group1.__repr__(exclude=['RefQN', 'QNs', 'nQN'])
        print 'QSp_Group2:\n\t', cls.QSp_Group2.__repr__(exclude=['RefQN', 'QNs', 'nQN'])        
        
    @classmethod
    def svd(cls, itensor, ndiv=None, nSd=None, pr=None):
        """
            svd of u*s*v.T =  t, while finally replace t with t = u*v.T
        """
        div = itensor.ndiv
        if ndiv is not None: 
            div = ndiv
        #iTensor_group_legs(T,div)
        cls.group_legs(itensor, div)
        
        nvSd = 0
        Sd = np.empty(2048, dtype=itensor.dtype)
        vSd = np.empty(2048, dtype=itensor.dtype)
        #print "sss", self.matrix_view(div)
        for gidx  in range(cls.QSp_Group1.nQN):
            V_buf = cls.get_block1(itensor, gidx=gidx)
            #print 'ggg ', gidx, V_buf
            nV, mV = V_buf.shape
            
            if 1:
                #issue225
                V_buf, vSd =common_util.matrix_svd_unitilize(V_buf)
            else:
                # below 2 lines are not needed any more at present, but keep it only for future use
                u, s, v=scipy.linalg.svd(V_buf, full_matrices=False)
                V_buf=common_util.matrix_multiply(u,v,1.0,0.0)

            
            nvSd = min(nV,mV)
            if pr:
                print cls.QSp_Group1.QNs[gidx].val, vSd[0:nvSd]
            if nSd is not None:  
                Sd[nSd+0:nSd+nvSd] = vSd[0:nvSd]
                nSd = nSd+nvSd
            cls.set_block1(itensor, gidx, V_buf)
        
        return Sd
    
    @classmethod
    def svd_new(cls, itensor, ndiv=None, nSd=None, pr=None):
        """
            the svd above is merely for mera minimization, this one 
            is the general one 
            status: not completed 
        """
        
        div = itensor.ndiv
        if ndiv is not None: 
            div = ndiv

        cls.group_legs(itensor, div)
        
        #QSp[div].QNs[0] = cls.QSp_Group1.QNs[gidx]
        #QSp[div].reverse()
        #
        #U = iTensor(div+1, QSp, tQN)
        
        qsp_l = [itensor.QSp[i].copy() for i in range(div)]
        qsp_r = [itensor.QSp[i].copy() for i in range(div,itensor.rank)]
        qsp_c = qsp_l[0].__class__.prod_many(qsp_l )
        qsp_c_rev = qsp_c_rev.copy(); qsp_c_rev.reverse( )
        U = iTensor(QSp=qsp_l + [qsp_c])
        S= iTensor(QSp=[qsp_c, qsp_c_rev])
        V = iTensor( QSp=[qsp_c_rev + qsp_r])
 
        Sd = np.empty(2048, dtype=itensor.dtype)
        vSd = np.empty(2048, dtype=itensor.dtype)
        for gidx  in range(cls.QSp_Group1.nQN):
            V_buf = cls.get_block1(itensor, gidx=gidx)
            #print 'ggg ', gidx, V_buf
            nV, mV = V_buf.shape
            
            if 1:
                #issue225
                u, s, v =common_util.matrix_svd(V_buf)
            else:
                # below 2 lines are not needed any more at present, but keep it only for future use
                u, s, v=scipy.linalg.svd(V_buf, full_matrices=False)
                V_buf=common_util.matrix_multiply(u,v,1.0,0.0)

            nvSd = min(nV,mV)
            if pr:
                print cls.QSp_Group1.QNs[gidx].val, vSd[0:nvSd]
            if nSd is not None:  
                Sd[nSd+0:nSd+nvSd] = vSd[0:nvSd]
                nSd = nSd+nvSd
            
            #cls.set_block1(itensor, gidx, V_buf)
            #the set_block1 below should be defined 
            cls.set_block1(U, u)
            cls.set_block1(S, s)
            cls.set_block1(V, v)
        
        return Sd
    
    @staticmethod 
    def svd_rank2(tensor, trunc_dim=10000, full_matrices=False, return_trunc_err=False, 
        normalize_singular_val=False):     
        """
            status: 
                tested  but not quit thorough 
            itensor should be prepare into a rank 2 tensor 
        """
        tt = tensor  # a shorter name 
        uu, ss, vv = {}, {}, {}
        dim = []
        qns = []
        for i in range(tt.nidx): 
            qn_id_tuple = tt.Addr_idx[:, i]
            qn0, qn1 = qn_id_tuple
            dl = tt.QSp[0].Dims[qn0]; dr = tt.QSp[1].Dims[qn1]
            p  = tt.Block_idx[0, i]
            size  = tt.Block_idx[1, i]
            assert dl*dr == size, (dl, dr, size) 
            mat = tt.data[p: p+dl*dr].reshape(dl, dr, order='F')
            if 0:
                #issue225
                #issue: if mat is of shape (1, 1), 那么经过svd，mat的值会被修改成 1.0 ！！！ 还不知道为什么，可能是f2py的bug
                #暂时先不用它，而改用numpy
                u, s, v = common_util.matrix_svd(min(dl, dr), mat)
                
            else:
                u, s, v=scipy.linalg.svd(mat, full_matrices=False)
            uu[i] = u
            ss[i] = s
            vv[i] = v
            dim.append(s.size)
            qn = tt.QSp[1].QNs[qn1]
            qns.append(qn.copy())
        
        #issue:  may happen d = 0, some qn is fully truncated
        if trunc_dim < np.sum(dim):  
            temp = np.zeros(np.sum(dim), dtype=float)
            d0 = 0
            for i, d in enumerate(dim): 
                temp[d0: d0 + d] = ss[i]
                d0 += d
            
            temp.sort()
            temp = temp[: : -1]
            val_min = temp[trunc_dim]
            norm = np.linalg.norm(temp[: trunc_dim])
            #trunc_err = np.sum(temp[trunc_dim: ])
            trunc_err = 1-norm**2  
            
            for i in range(tt.nidx): 
                s= ss[i]
                s= s[np.where(s>val_min)]
                ss[i] = s/norm 
                d = s.size 
                dim[i] = d
                uu[i] = uu[i][:, :d]
                vv[i] = vv[i][:d, :]
        else:
            trunc_err = 0.0
        qsp_sr = tt.QSp[0].__class__(len(dim), qns, dim)
        qsp_sl = qsp_sr.copy(); qsp_sl.reverse()
        U = iTensor(QSp=[tt.QSp[0], qsp_sr])
        S = iTensor(QSp=[qsp_sl, qsp_sr])
        V = iTensor(QSp=[qsp_sl, tt.QSp[1]])
        
        #print_vars(vars(), ['U', 'S','V' ])
        for i in xrange(U.nidx): 
            p  = U.Block_idx[0, i]
            size  = U.Block_idx[1, i]
            U.data[p: p + size] = uu[i].ravel(order='F')
            
            p  = V.Block_idx[0, i]
            size  = V.Block_idx[1, i]
            V.data[p: p + size] = vv[i].ravel(order='F')
            
            p  = S.Block_idx[0, i]
            size = S.Block_idx[1, i]
            S.data[p: p + size] = np.diag(ss[i]).ravel(order='F')
        if not return_trunc_err:
            return U, S, V 
        else:
            return U, S, V , trunc_err 
 
    @classmethod
    def random_unit_tensor(cls,itensor, d):
        itensor.randomize_data()  #here method randomize_data inheriets from TensorBase
        cls.svd(itensor, ndiv=d)
    
    @classmethod
    def eig(cls, itensor, totQN=None):
        """
            see iTensor_Eig in f90
            thif func is only used in calc central_charge
            params:
                totQN:  can be used to select excitation states
            return:
                Vg: eigen vec of itensor with greatest eigval
                Eg: eigen energy 
        """
        rank = itensor.rank
        div = rank//2
        cls.group_legs(itensor, div)
        E_size=1024*64
        #lwork = 5*E_size
       
        QSp = [itensor.QSp[i].copy() for i in range(div)]
        if totQN is None:
            totQN = QSp[0].QnClass.qn_id()
        QSp.append(QSp[0].null())
        tQN = QSp[0].QNs[0].qn_id()

        #E = np.empty(E_size)
        E = np.ndarray(E_size, order="F")
        pos = np.empty(div+1, np.int)
 
        for gidx  in range(cls.QSp_Group1.nQN):
            #print "ggg", gidx, 
            if  cls.QSp_Group1.QNs[gidx] !=  totQN:  
                continue
            V_buf = cls.get_block1(itensor, gidx)
            
            common_util.matrix_eigen_vector(V_buf, E)
            #E, V_buf = np.linalg.eig(V_buf)
            #need sort !!!
            #print E[:3]/2.0
            
            V_buf = V_buf.ravel(order="F")

            Eg = E[0]
            QSp[div].QNs[0] = cls.QSp_Group1.QNs[gidx]
            QSp[div].reverse()
            Vg = iTensor(div+1, QSp, tQN)
            Vg.data[:] = 0.0
            
            for idx in range(itensor.nidx):
                if cls.QN_Group[0,idx] != gidx:  
                    continue
                p1 = cls.QN_Group[1,idx]; p2=cls.QN_Group[2,idx]
                nT = cls.QNG_Addr1[2,p1]; mT = cls.QNG_Addr2[2,p2]
                x = cls.QNG_Addr1[1,p1]; y = cls.QNG_Addr2[1,p2]
                pos[0:div] = itensor.Addr_idx[0:div,idx]
                pos[div] = 0

                pV = Vg.get_position(pos)
                p = Vg.Block_idx[0, Vg.idx[pV]]
                #Vg.data[p:p+nT-1] = V_buf[x+1:x+nT]
                Vg.data[p:p+nT] = V_buf[x:x+nT]
        return Vg, E

    @classmethod
    def eig_sparse(cls, itensor, qn_list, k=10, tol=None):
        """
            eigenvalue decomposition of sparse tensor
        """
        #tol = 1e-18 if tol is None else tol
        tol = tol if tol is not None else 1e-12
        
        rank = itensor.rank
        div = rank//2
        cls.group_legs(itensor, div)
        
        res= {}
        #print "QSp_Group input is:\n%s"%cls.QSp_Group1
        #print "QSp_Group output is:\n%s"%cls.QSp_Group2
        #print "selected input qn are %s"%qn_list
        for gidx in range(cls.QSp_Group1.nQN):
            #print "ggg", gidx, cls.QSp_Group1.QNs[gidx]
            qn = cls.QSp_Group1.QNs[gidx]  
            if qn in qn_list:  
                V_buf = cls.get_block1(itensor, gidx)
                d = V_buf.shape[0]
                #print 'dddd', d
                if d <=  20:
                    vals, vecs= linalg.eig(a=V_buf)
                elif k>d:
                    k = d-2
                    vals, vecs= eigs(A=V_buf, k=k, tol=tol )
                else:
                    vals, vecs= eigs(A=V_buf, k=k, tol=tol )
                #res[gidx] = qn, vals
                res[qn._val] = vals

            else:
                print ("qn %s not found in qn_list %s, QSp_Group1.QNs are %s"%(qn,  qn_list, cls.QSp_Group1.QNs))

        return res
    
    @classmethod
    def eigs(cls, itensor): 
        pass 
    
    @classmethod
    def get_block(cls, div, target_QN, gidx= -1, need_group=True):
        """ 
        """
        #gidx = -1

        if need_group:
            cls.group_legs(div)
        
        if gidx == -1: 
            for i in range(cls.QSp_Group1.nQN):
                if cls.QSp_Group1.QNs[i] == target_QN:  
                    gidx = i
                    break  
        if gidx == -1 :  
            print 'Can not find target QN'
            exit()
        else:
            return cls.get_block1(gidx)
    
    @classmethod
    def get_block1(cls, itensor, gidx):
        """ 
        """
        
        nV = cls.QSp_Group1.Dims[gidx]
        mV = cls.QSp_Group2.Dims[gidx]        
        VV = np.empty((nV, mV), dtype=itensor.dtype, order='F')
        
        for idx  in range(itensor.nidx):
            if cls.QN_Group[0,idx] != gidx:  
                continue
            p1 = cls.QN_Group[1,idx]; p2=cls.QN_Group[2,idx]
            nT = cls.QNG_Addr1[2,p1]; mT = cls.QNG_Addr2[2,p2]
            x = cls.QNG_Addr1[1,p1]; y = cls.QNG_Addr2[1,p2]
            p = itensor.Block_idx[0,idx]
            
            if nT == 1: 
                temp = np.ndarray((nT, mT), order='F')
                temp[:,:] = itensor.data[p:p+nT*mT].reshape((nT, mT)) #[:, :]
            else:
                temp = itensor.data[p:p+nT*mT].reshape((nT, mT), order='F')
            #common_util.set_matrix(VV,temp , x,y,  True)
            VV[x:x + nT, y:y + mT] = temp
        return VV
    
    @classmethod
    def set_block(cls, div, target_QN, VV):
        """
            differ with get_block by set_matrix
        """
        cls.group_legs(div)

        for gidx  in range(cls.QSp_Group1.nQN):
            if cls.QSp_Group1.QNs[gidx] == target_QN:  
                break
        
        if gidx > cls.QSp_Group1.nQN:  
            print 'Can not find target QN'
            exit()
        cls.set_block1(gidx, VV)
       
    @classmethod
    def set_block1_bac(cls, itensor, gidx, VV):
        nV = cls.QSp_Group1.Dims[gidx]
        mV = cls.QSp_Group2.Dims[gidx]        
        
        for idx  in range(itensor.nidx):
            if cls.QN_Group[0,idx] != gidx:  
                continue
            p1 = cls.QN_Group[1,idx]; p2=cls.QN_Group[2,idx]
            nT = cls.QNG_Addr1[2,p1]; mT = cls.QNG_Addr2[2,p2]
            x = cls.QNG_Addr1[1,p1]; y = cls.QNG_Addr2[1,p2]
            p = itensor.Block_idx[0,idx]
            temp = np.ndarray((nT, mT), order='F')
            if 0:
                common_util.set_matrix(VV, temp, x,y,  False )
            else:
                #common_util.set_matrix(VV, temp, x,y,  False )
                temp[:, :] = VV[x:x + nT, y:y + mT]

            #attention_here  fortran order must be used
            itensor.data[p:p+nT*mT] = temp.ravel('F')

    @classmethod
    def set_block1(cls, itensor, gidx, VV):
        nV = cls.QSp_Group1.Dims[gidx]
        mV = cls.QSp_Group2.Dims[gidx]        
        
        for idx  in range(itensor.nidx):
            if cls.QN_Group[0,idx] != gidx:  
                continue
            p1 = cls.QN_Group[1,idx]; p2=cls.QN_Group[2,idx]
            nT = cls.QNG_Addr1[2,p1]; mT = cls.QNG_Addr2[2,p2]
            x = cls.QNG_Addr1[1,p1]; y = cls.QNG_Addr2[1,p2]
            p = itensor.Block_idx[0,idx]
            temp = np.ndarray((nT, mT), order='F')
            if 0:
                common_util.set_matrix(VV, temp, x,y,  False )
            else:
                #common_util.set_matrix(VV, temp, x,y,  False )
                temp[:, :] = VV[x:x + nT, y:y + mT]

            #attention_here  fortran order must be used
            itensor.data[p:p+nT*mT] = temp.ravel('F')


class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
        from tensor import  test_iTensor
        #u,  w, u2234= simple_itensor()
        #u = test_iTensor.instance("u")
        #w = test_iTensor.instance("w")
        u = test_iTensor(symmetry='U1').u
        self.u = u 
        self.qn_identity, self.qsp_base, self.qsp_null = init_System_QSp(symmetry='Z2')
    
    def test_temp(self): 
        pass 
    
    def test_group_legs(self, rank=3):
        """
            --- pass 
        """
        u = iTensor(rank, self.qsp_base.copy_many(rank), self.qn_identity)

        Tensor_svd.group_legs(u, 2)
        Tensor_svd.show_group(u)
        
    def test_eig(self):
        rank = 4
        u = iTensor(rank, self.qsp_base.copy_many(rank), self.qn_identity)
        u.data[:] = 1.0
        print "tensor to be eiged", u
        totqn = self.qn_identity.copy()
        print "tttt", type(totqn)
        a, b=Tensor_svd.eig(u, totqn)
        print "a", a, "b", b
        pass
   
    def test_svd(self):
        """  
            --- not sure
        """
        #u = cls.u.copy()
        u = self.u 
        u.data[:] = 0
        Tensor_svd.svd(u, ndiv=2)
        #u.QSp[0].reverse()
        #u.QSp[1].reverse()
        #u.QSp[2].reverse()
        print u.data.round(5)
    
    def test_svd_rank2(self): 
        np.random.seed(1234)
        t = iTensor.example(rank=2, symmetry='U1')
        t.data[: ] = np.random.random(t.totDim)
        U, S, V=Tensor_svd.svd_rank2(t)
        print S.data 
        
        U, S, V = Tensor_svd.svd_rank2(t, trunc_dim=5)
        print  S.data  
        data = [ 1.51112271, 0.,0., 0.34817902, 0.35781727, 0.50099513]
        
        self.assertTrue(np.allclose(S.data, data, atol=1e-8))
        print_vars(vars(), ['U.shape', 'V.shape'])
        
    def test_svd_rank2_2(self): 
        if 1: 
            np.random.seed(1234)
            qsp = QspU1.easy_init([0, 1, -1], [8, 5, 5])
            qsp = qsp.copy_many(2, reverse=[1])
            t = iTensor.example(qsp=qsp, rank=2, symmetry='U1')
            data = np.random.random(t.totDim)
            t.data[: ] = data
            t1 = t
        
        if 1: 
            np.random.seed(1234)
            qsp = QspZ2.easy_init([1, -1], [8, 8])
            qsp = qsp.copy_many(2, reverse=[1])
            t = iTensor.example(qsp=qsp, rank=2 ) 
            data = np.random.random(t.totDim)
            t.data[: ] = data
            t2 = t 
        if 0:     
            for t in [t1, t2]: 
                U, S, V=Tensor_svd.svd_rank2(t, trunc_dim=100)
                vv=V.dot(V.conjugate(1))
                uu = U.conjugate(1).dot(U)
                self.assertTrue(uu.is_close_to(1))
                self.assertTrue(vv.is_close_to(1))
                a = U.dot(S).dot(V)
                #print_vars(vars(), ['a.to_ndarray().shape', 't.to_ndarray().shape'])
                #print_vars(vars(), ['a.to_ndarray()[:3][:3]', 't.to_ndarray()[:3][:3]'])
                self.assertTrue(np.allclose(a.to_ndarray(), t.to_ndarray(), 1e-14)) 
        
        
        if 1:     
            for ii, t in enumerate([t1, t2]): 
            #for t in [t1]: 
                t.data = t.data/t.norm()
                U, S, V=Tensor_svd.svd_rank2(t, trunc_dim=12)
                vv=V.dot(V.conjugate(1))
                uu = U.conjugate(1).dot(U)
                self.assertTrue(uu.is_close_to(1))
                self.assertTrue(vv.is_close_to(1))
                a = U.dot(S).dot(V)
                #print_vars(vars(), ['a.to_ndarray().shape', 't.to_ndarray().shape'])
                #print_vars(vars(), ['a.to_ndarray()[:3][:3]', 't.to_ndarray()[:3][:3]'])
                #self.assertTrue(np.allclose(a.to_ndarray(), t.to_ndarray(), 1e-14)) 
                overlap = t.contract(a, [0, 1], [0, 1])[0].data[0]
                print_vars(vars(), ['t.norm()', 't.to_ndarray().shape', 'overlap'])
                old = {0:0.996667569187076, 1: 0.99875849173495423}[ii]
                self.assertAlmostEquals(overlap, old, 12)
                
    def random_unit_tensor(cls):
        u = self.u
        Tensor_svd.random_unit_tensor(u, 2)
        print u.data
    
    @classmethod
    def random_unit_tensor_large(cls):
        #import os
        #os.environ["OMP_NUM_THREADS"] = str(threads)
        qn_identity, qsp_base, qsp_null = init_System_QSp("Travial")
        qsp_max , qsp_max2= qsp_base.__class__.max(4)
        QSbase = qsp_max.copy
        totQN = qn_identity.copy
        ranku = 8
        u=iTensor(ranku,[QSbase() for i in xrange(ranku)], totQN())
        Tensor_svd.random_unit_tensor(u, 2)
   
    
if __name__ == "__main__":
    if 0: 
        tt = test_Tensor_svd(symmetry="Travial")
        #tt.group_legs(rank=3)
        tt.svd()
        #tt.eig()
        #tt.random_unit_tensor()
        #tt.random_unit_tensor_large()
    if 0: #examine
        if 0: 
            from concurrencytest import ConcurrentTestSuite, fork_for_tests
            loader = unittest.TestLoader()
            suite = []
            for cls in [TestIt]: 
                temp = loader.loadTestsFromTestCase(cls)
                suite.append(temp)
            suite = unittest.TestSuite(suite)
            
            suite = ConcurrentTestSuite(suite, fork_for_tests(2))
            unittest.TextTestRunner(verbosity=0).run(suite)
        else: 
                
            #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
            #unittest.TextTestRunner(verbosity=0).run(suite)    
            TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
            unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_temp', 
           #'test_eig', 
           #'test_svd', 
           #'test_svd_rank2', 
           #'test_svd_rank2_2', 
           #'test_group_legs', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)
    
    




