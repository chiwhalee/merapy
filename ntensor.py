#!/usr/bin/env python
#coding=UTF8


"""
    the name of nTensor  can be interpreted as normal ternsor or tensor that inhereitated
    from np.ndarray 
"""
import unittest 
import numpy as np 
import cPickle as pickle 

from merapy import crandom
from merapy.utilities import print_vars

class TensorBase(object):
    """
    an abstract tensor base class 
    this class is newly added 
    """
    MaxRank= 32  # maximal rank of tensors supported
    dtype="float"

    def __init__(self,rank=None, ind_labels=None,use_buf=False):
        """
        
        ind_labels: dict
        type_name: may be "U", "V", "V2" etc
        """
        self.rank=rank
        self.ind_labels= {i:None for i in xrange(rank)}
        if ind_labels is not None:
            #print 'iii', ind_labels
            
            for i in xrange(rank):
                print 'iii', i, ind_labels
                self.ind_labels[i] = ind_labels[i]
        self.data = None
        self.type_name = ""

    def set_ind_labels(self, ind_labels):
        """
        ind_labels:  a dict whose keys are from [0,rank-1], values arbitrary
        """
        self.ind_labels.update(ind_labels)
    
    def get_label_list(self):
        res= [self.ind_labels[i] for i in xrange(self.rank)]
        return res
    def _dims(self):
        return self.Dims[:self.rank]

    def copy_deprecated(self):
        from copy import deepcopy
        res= deepcopy(self)
        return res
    
    def shallow_copy_deprecate(self):
        """
        attention:  still deep copy, I wrote this only for debuging 
        """
        from copy import deepcopy
        res= deepcopy(self)
        return res

    def matrix_view(self, p, order='F'):
        """
        newly_added
        map rank (p, q) tensor to a 2-d mat. 
        """
        dim1 = np.prod(self.Dims[:p])
        dim2 = np.prod(self.Dims[p:self.rank])
        res=self.data.reshape((dim1, dim2), order=order)
        return res

    def randomize_data(self):
        """
        status_1_with_problem 
        see Random_Tensor in f90
        this method is not reasonable,  since totDim and data are not defined for TensorBase
        """
        #try:
            #self.data[:] = np.random.rand(self.totDim, dtype=self.dtype)
            #attention_this_may_be_wrong  rand does not support dtype
            #self.data[:self.totDim] = np.random.random(self.totDim)

        self.data[:self.totDim] = [crandom.rand()-0.5 for i in xrange(self.totDim)]

        #except: print "error, self.data is not materialized"

    def scale(self,alpha,inplace=False):
        """
        status_1_verified
        see Tensor_Scale in f90
        """
        
        if inplace:  
            #for i  in xrange(self.totDim): self.data[i] = alpha*self.data[i]            
            #self.data[:]= self.data[:]*alpha
            self.data = self.data*alpha
        else:
            T2 = iTensor(self.rank, self.QSp, self.totQN)
            #for i in xrange(self.totDim):  T2.data[i] = alpha*self.data[i]
            T2.data = self.data*alpha
            return T2
    
    def add_with_scale(self,T2, alpha,beta):
        """
        see Tensor_Add_WithScale
        """
        reason= self.is_same_shape(T2)
        if  reason<0:  
            print 'Error in Tensor_Add_WithScale, self,other different shape:', reason
            #print 'self.' + reason + ":", 
            #print self.__getattribute__(reason)
            #print 'other.' + reason + ":", 
            #print T2.__getattribute__(reason)
            #exit()
            raise iTensor.ShapeError

        #T3 = Tensor()
        #T3.init_Tensor(self.rank, self.QSp, self.totQN)
        T3 = self.copy()
        #T3.data[:] = alpha*self.data[:] + beta*T2.data[:]
        T3.data = alpha*self.data + beta*T2.data
        return T3

    def __iadd__(self, other):
        self.data += other.data 
        return self

    def __add__(self, other):
        res= self.copy()
        res.data = self.data + other.data 
        return res
    
    def __radd__(self, other):
        res= self.copy()
        if isinstance(other, float): 
            res.data = self.data + other
        else: 
            res.data = self.data + other.data 
        return res


    def __sub__(self, other):
        res= self.copy()
        res.data = self.data-other.data
        return res

    def __rmul__(self, scalar):
        #self.data *= scalar  # this is wrong! 
        #return self
        res= self.copy()
        res.data *= scalar
        return res 
    
    
    def __imul__(self, scalar):
        self.data *= scalar
        return self

    def add_with_scale_new(self,T2, alpha,beta):
        """
        bug in it
        see Tensor_Add_WithScale
        """
        same, reason= self.is_same_shape(T2)
        if  not same:  
            print 'Error in Tensor_Add_WithScale, self,other different shape:', reason
            print 'self.' + reason + ":", 
            print self.__getattribute__(reason)
            print 'other.' + reason + ":", 
            print T2.__getattribute__(reason)
            exit()
        T3 = Tensor()
        T3.init_Tensor(self.rank, self.QSp, self.totQN)
        #for i in xrange(T3.totDim):
        T3.data[:] = alpha*self.data[:] + beta*T2.data[:]
        return T3
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

    def save(self, fn):
        out = open(fn, "wb")
        pickle.dump(self, out)
        out.close()
        
    @staticmethod
    def load(fn):
        inn = open(fn, "rb")
        res  = pickle.load(inn)
        return res

# DONT DELETE THE OLD nTensor,  MAY STILL NEEDED 
#below is deprecated,  as it is not quite compatible with np.ndarray; but it may be of use later
if 0:
    pass 
    #class nTensor(TensorBase):
    #    def __init__(self, rank, Dims,ind_labels=None, shallow=None, use_buf=None):
    #        """
    #            normal, non-symmetric Tensors
    #        """
    #        
    #        TensorBase.__init__(self, rank=rank, ind_labels=ind_labels)

    #        self.allocated=False
    #        self.use_buf=use_buf
    #        self.shallow=shallow
    #        use_buf0 = False

    #        #self.max_data_size=0

    #        if use_buf is not None:  
    #            use_buf0 = use_buf
    #        #$ Use buffer        
    #        if use_buf0 or self.use_buf:
    #            if self.buf_ref[0] == -1:
    #            #$use_tBuffer(rank, Dims, T)
    #                return
    #        #$ Already use buffer, but rank change, use different one
    #        if self.use_buf:
    #            if  ((rank<=6) and (self.rank == 6) )or ((rank == 6) and (self.rank<=6)) :
    #                  #$delete_Tensor(T)
    #                #$use_tBuffer(rank, Dims, T)
    #                return
    #            
    #        if rank > self.MaxRank:
    #            print "tensor rank exceeds bounds"
    #            print "['Tensor rank=',I3, '> MaxRank', I3]", rank, self.MaxRank
    #            exit()
    #            #attention_this_may_be_wrong
    #            #stop
    #        
    #        if self.shallow:
    #            if  shallow==None:
    #                print  "Error, can not initialized again a shallow tensor"
    #                #stop
    #                exit()
    #        
    #        #self.rank = rank
    #        #attention_this_may_be_wrong  I modified self.MaxRank into self.rank for self.Dims
    #        
    #        if rank == 0:
    #            self.Dims=np.array([1],"int")   #[0] = 1
    #            self.totDim = 1
    #        elif rank>0:        
    #            self.Dims=np.ndarray(self.MaxRank,"int")
    #            self.Dims[:rank] = Dims[0:rank]
    #        
    #        self.totDim = np.prod(Dims[:rank])

    #        self.data=np.empty(self.totDim,dtype=self.dtype)
    #        
    #        if shallow is not None : 
    #            self.shallow = shallow
    #            self.allocated = False
    #            #nullify[self.data]
    #            if shallow:
    #                self.allocated = True
    #                return
    #            
    #        self.shallow = False
    #        #if (self.totDim > #self.max_data_size) or ( not self.allocated):
    #        #attention_omitted_something
    #        if 0:
    #            #self.max_data_size = self.totDim
    #            if self.allocated:
    #                #deallocate[self.data]
    #                pass
    #                #attention_this_may_be_wrong
    #            
    #            self.allocated  == True
    #            #allocate[self.data[self.totDim]]
    #            self.data = 0.0
    #            if self.use_buf:
    #                pass
    #        #$iTensor.T_BUFFER[self.buf_ref[1]].T[self.buf_ref[2]].data => self.data
    #            
    #    def permutation(self,P, perm_order='F'):
    #        """
    #        see nTensor_Permutation in f90
    #        P: 是张量 leg 编号的重排
    #        这里主要是调用array_permutation来实现

    #        out: 
    #            Tp: 和self rank相同，但每个leg的维数不同（置换了）

    #        nshape[i]:
    #        """
    #        nshape = np.ndarray(self.rank, "int")
    #        norder = np.ndarray(self.rank, "int")
    #        for i  in xrange(self.rank):
    #            nshape[i] = self.Dims[P[i]]
    #            norder[P[i]] = i            
    #        
    #        #如果rank=1,travial 情况...
    #        Tp= nTensor(self.rank, nshape)
    #        if self.rank<=1:
    #            Tp.data[0:Tp.totDim] = self.data[0:self.totDim]
    #            return Tp
    #        
    #        #判断次序是否变化，如果leg次序和原来的相同，则直接复制self给Tp，return
    #        hasChange = False
    #        for i  in xrange(self.rank):
    #            if norder[i] !=i:
    #                hasChange = True
    #        if  not hasChange:
    #            Tp.data[0:Tp.totDim] = self.data[0:self.totDim]
    #            return  Tp


    #        #如果次序变化，进行如下操作
    #        if perm_order=="F":
    #            Tp.data=array_permutation.array_permutation(self.data,self.rank, self.Dims, P  )
    #        elif perm_order=="C":
    #            Tp.data=array_permutation.array_permutation_C(self.data,self.rank, self.Dims, P  )

    #        return Tp
    #    
    #    def get_position_rev(self, p):
    #        pos=common_util.matrix_get_position_rev(p,self.Dims)
    #        return pos

    #    def get_position(self, pos):
    #        """
    #        see nTensor_GetPosition in f90       
    #        """
    #        p=common_util.matrix_get_position(pos,self.Dims)
    #        return p
    #        
    #    def __repr__(self, print_data=None):
    #        """
    #        """
    #        
    #        keys=["rank", "Dims","totDim","data"]

    #        str0="----Begin Tensor----------------------------------------------\n"

    #        str1=""
    #        str2= '----End    Tensor----------------------------------------------\n'
    #        for k in keys:
    #            temp = repr(self.__dict__[k])
    #            if k == 'data':
    #                temp = repr(self.__dict__[k][:self.totDim].round(10)) 
    #            elif k == 'Dims':
    #                temp  = repr(self.__dict__[k][:self.rank]) 
    #            str1 += k+":\n\t" + temp +"\n" 

    #        res=str0+ str1 +str2
    #        
    #        return res

    #    def svd(self, d,nsd=None):
    #        """
    #            return:
    #                singlular values
    #        """
    #        dim1 = np.prod(self.Dims[:d])
    #        dim2 = np.prod(self.Dims[d:self.rank])
    #        
    #        min_dim = min(dim1, dim2)
    #        #allocate[U[dim1,min_dim], VT[min_dim, dim2], S[min_dim]]
    #        #U = np.ndarray((dim1, min_dim), self.dtype)
    #        #VT = np.ndarray((min_dim, dim2), self.dtype)
    #        a=self.data.reshape(dim1,dim2)
    #        U, S, VT = common_util.matrix_svd(min_dim, a)                
    #        
    #        #T2 = nTensor(self.rank, self.Dims)
    #        #common_util.matrix_multiply(dim1, min_dim, dim2, U, VT, T2.data, 1.0, 0.0)
    #        #T2.data=common_util.matrix_multiply(U.reshape(dim1,min_dim),VT.reshape(min_dim,dim2),1.0,0.0)
    #        #T2.data=T2.data.ravel()
    #        if nsd:
    #            Sd = S[:nsd]
    #        else:
    #            Sd=S
    #        return  Sd
    #        #if nsd!=None: nsd = min_dim
    #        #if Sd!=None: Sd[0:min_dim] = S[0:min_dim]
    #        #deallocate[U,VT,S]
    #        #attention_omitted_something

    #    def get_element(self, ind_vals):
    #        """
    #        parameters:
    #            ind_val: a list of integers specify values of each indices
    #        指标前低后高
    #        """
    #        p = 0
    #        for i in xrange(self.rank-1,-1, -1):
    #            p = p*self.Dims[i]+ind_vals[i]
    #        #p = p+1
    #        X = self.data[p]
    #        return X

    #    def set_element(self, ind_vals, X):
    #        """
    #        see nTensor_SetElement in f90
    #        这里ind_vals= [i0, i1, ..., in-1]
    #        前面的指标在低位，后高；不同于通常的惯例前高后低
    #        Dims: an iterable
    #        """
    #        #p: position
    #        p = 0  
    #        for i in xrange(self.rank-1, -1, -1):
    #            #p = p*self.Dims[i]+ind_vals[i]-1
    #            p = p*(self.Dims[i]) + ind_vals[i]  
    #        #p = p+1  
    #        #注意这里Fortran 和py的区别
    #        self.data[p] = X

    #    def unit_tensor(self):
    #        """
    #        status_0p9
    #        see unit_nTensor in_f90
    #        attention_this_may_be_wrong
    #        maybe this method shouldn't be in a class
    #        """
    #        d = 1
    #        for i  in xrange(self.rank/2):
    #            d = d*self.Dims[i]
    #        self.data=np.identity(d)

    #    def contract(self,T2, V1, V2,  out_Vc=False, info=0,  use_buf=None):
    #        """
    #        V1,2,3 are arrays (maps) 张量指标 —> 自然数。用自然数来标记所有张量的指标
    #        在V1,V2中可能有相同的元素，存在s_1n2中，
    #        
    #        T1, T2 contract yielding T3
    #        Args:
    #            Vi: of type np.ndarray(,"int")

    #            Vp1,先记录了T1的外腿，后记录内腿指标； Vp2先记录了内腿，后记录了外腿指标
    #        """
    #        T1=self
    #        #s1_list=list(S1)
    #        #s2_list=list(S2)
    #        V1=V1[:self.rank]
    #        V2=V2[:T2.rank]
    #        V1=list(V1)
    #        V2=list(V2)
    #        
    #        V_1n2= [i for i in V2 if i in V1]   #V1 \cap V2
    #        V_1m2= [i for i in V1 if i not in V2]  #V1-V2
    #        V_2m1=[i for i in V2 if i not in V1]   #V2-V1

    #        if info>0:
    #            print "rank1, rank2, V1, V2:\n\t",self.rank, T2.rank,  V1, V2
    #            print "contracted legs:\n\t", V_1n2

    #        #rank3 = len(s_1u2 - s_1n2)
    #        rank3=len(OrderSet.symm_diff(V1,V2))
    #        V3=[None]*rank3

    #        k=0
    #        l=0
    #        totDim1=1   #totDim1,T1外腿总维数,在下面MatrixMultiply中要用到
    #        Vp1=np.empty(self.MaxRank,'i')
    #        Vp2=np.empty(self.MaxRank,'i')
    #        Dims=np.empty(self.MaxRank,'i')
    #        for i  in xrange(T1.rank):
    #            if V1[i] not in V_1n2:
    #                Vp1[k] = i   
    #                # 这里之所以用k，而非直接用Vp1[k]是因为下面k的值要继续，对于l 也一样
    #                # p意指position，记录的T1中没有收缩的指标(leg)的编号
    #                Dims[l] = T1.Dims[i]
    #                totDim1 = totDim1*T1.Dims[i]
    #                V3[l] = V1[i]   #T3 来自T1 V1的外腿
    #                l = l+1
    #                k = k+1

    #        totDim3 = 1
    #        #这一段程序做了两件事：1.验证内线上维数相等，2.计算了totDim3
    #        #V_1n2=list(s_1n2)
    #        for i in xrange(len(V_1n2)):
    #            pos =  V1.index(V_1n2[i])
    #            pos2 = V2.index(V_1n2[i])
    #            Vp1[k] = pos
    #            totDim3 = totDim3*T1.Dims[pos]
    #            k = k+1  #注意这里k接着上面的值了
    #            
    #            if T1.Dims[pos] != T2.Dims[pos2]:
    #                print "['Error, size of contract tensor does not match']"
    #                print "T1.rank, T2.rank = ", T1.rank, T2.rank
    #                print V1, V2
    #                raise Exception
    #        #print "totDim3:\n\t",totDim3
    #        k=0
    #        for i  in xrange(len(V_1n2)):
    #            pos=V2.index(V_1n2[i])
    #            Vp2[k] = pos   
    #            k = k+1
    #       
    #        totDim2 = 1   #dim of external leg  of T2
    #        for i  in xrange(T2.rank):
    #            if V2[i] not in V_1n2:
    #                Vp2[k] = i
    #                Dims[l] = T2.Dims[i]
    #                totDim2 = totDim2*T2.Dims[i]
    #                V3[l] = V2[i]  #T3 V3 来自T2 V2的外腿
    #                k = k+1
    #                l = l+1
    #         
    #        T3=nTensor(rank=rank3,Dims=Dims,use_buf=use_buf)
    #        
    #        if info>0: 
    #            print "Vp1Vp2:\n\t", Vp1[:len(V1)],Vp2[:len(V2)]

    #        nT1=T1.permutation(Vp1)    #把T1 按照 Vp1 重排
    #        nT2=T2.permutation(Vp2)

    #        
    #        data1=nT1.data.reshape((totDim1,-1))
    #        data2=nT2.data.reshape((-1,totDim2))

    #        #print "data", type(data1), type(data2), data1.shape, data2.shape
    #        #print data1.dtype, data2.dtype

    #        data3=common_util.matrix_multiply(data1, data2, 1.0, 0.0)
    #        T3.data=data3.ravel()
    #        #T3.data=nT1.data*nT2.data

    #        
    #        #nT1.delete()
    #        #nT2.delete()
    #        #attention_this_may_be_wrong nT1,nT2 应该会自动 garbage collection
    #        if out_Vc:
    #            Vc= V_1n2
    #            return  T3, V3, Vc 
    #        else:
    #            return T3, V3

    #    def contract_new(self, T2, perm_order="F", mul_algorithm="lapack", info=0, nVc=None, Vc=None,use_buf=None):
    #        """
    #        V1,2,3 are arrays (maps) 张量指标 —> 自然数。用自然数来标记所有张量的指标
    #        在V1,V2中可能有相同的元素，存在s_1n2中，
    #        
    #        T1, T2 contract yielding T3
    #        Args:
    #            Vi: of type np.ndarray(,"int")

    #            Vp1,先记录了T1的外腿，后记录内腿指标； Vp2先记录了内腿，后记录了外腿指标
    #        """
    #        V1=self.get_label_list()
    #        V2=T2.get_label_list()
    #        V3=[None]*32
    #        

    #        T1=self
    #        
    #        V_1n2= [i for i in V2 if i in V1]   #V1 \cap V2
    #        V_1m2= [i for i in V1 if i not in V2]  #V1-V2
    #        V_2m1=[i for i in V2 if i not in V1]   #V2-V1

    #        if nVc !=None:
    #            nVc = len(S)
    #            Vc=np.array(list(S),"int")
    #        rank3 =len(V_1m2) + len(V_2m1)    #len(s_1u2 - s_1n2)

    #        k=0
    #        l=0
    #        totDim1=1   #totDim1,T1外腿总维数,在下面MatrixMultiply中要用到
    #        Vp1=np.empty(self.MaxRank,'i')
    #        Vp2=np.empty(self.MaxRank,'i')
    #        Dims=np.empty(self.MaxRank,'i')
    #        for i  in xrange(T1.rank):
    #            if V1[i] not in V_1n2:
    #                Vp1[k] = i   
    #                # 这里之所以用k，而非直接用Vp1[k]是因为下面k的值要继续，对于l 也一样
    #                # p意指position，记录的T1中没有收缩的指标(leg)的编号
    #                Dims[l] = T1.Dims[i]
    #                totDim1 = totDim1*T1.Dims[i]
    #                V3[l] = V1[i]   #T3 来自T1 V1的外腿
    #                l = l+1
    #                k = k+1

    #        totDim3 = 1
    #        #这一段程序做了两件事：1.验证内线上维数相等，2.计算了totDim3
    #        for i in xrange(len(V_1n2)):
    #            pos =  V1.index(V_1n2[i])
    #            pos2 = V2.index(V_1n2[i])
    #            Vp1[k] = pos
    #            totDim3 = totDim3*T1.Dims[pos]
    #            k = k+1  #注意这里k接着上面的值了
    #            
    #            if T1.Dims[pos] != T2.Dims[pos2]:
    #                print "['Error, size of contract tensor does not match']"
    #                print "T1.rank, T2.rank = ", T1.rank, T2.rank
    #                raise Exception
    #        k=0
    #        for i  in xrange(len(V_1n2)):
    #            pos=V2.index(V_1n2[i])
    #            Vp2[k] = pos   
    #            k = k+1
    #       
    #        totDim2 = 1   #dim of external leg  of T2
    #        for i  in xrange(T2.rank):
    #            if V2[i] not in V_1n2:
    #                Vp2[k] = i
    #                Dims[l] = T2.Dims[i]
    #                totDim2 = totDim2*T2.Dims[i]
    #                V3[l] = V2[i]  #T3 V3 来自T2 V2的外腿
    #                k = k+1
    #                l = l+1
    #            
    #        T3=nTensor(rank3,Dims,use_buf)
    #        
    #        #print "new order of legs Vp1Vp2:\n", Vp1[:len(V1)],Vp2[:len(V2)]
    #        
    #        nT1=T1.permutation(Vp1,perm_order=perm_order)    #把T1 按照 Vp1 重排
    #        nT2=T2.permutation(Vp2,perm_order=perm_order)

    #        

    #        if mul_algorithm=="lapack":
    #            data1=nT1.data.reshape((totDim1,-1))
    #            data2=nT2.data.reshape((-1,totDim2))
    #            data3=common_util.matrix_multiply(data1, data2, 1.0, 0.0)
    #        
    #        elif mul_algorithm== "numpy":
    #            data1=nT1.data.reshape(nT1.Dims[:nT1.rank])
    #            data2=nT2.data.reshape(nT2.Dims[:nT2.rank])
    #            data3=data1.dot(data2)

    #        T3.data=data3.ravel()
    #        if info>0:
    #            print "V1, V2", V1, V2
    #            print "Vp1, Vp2", Vp1[:T1.rank], Vp2[:T2.rank]
    #            print "contracted legs:", V_1n2
    #            print "leg order of new tensor:", V3[:T3.rank]
    #        
    #        #T3.data=nT1.data*nT2.data
    #        
    #        #nT1.delete()
    #        #nT2.delete()
    #        #attention_this_may_be_wrong nT1,nT2 应该会自动 garbage collection
    #        return T3

    #    def contract_np(self, T2):
    #        """
    #        another implementation of contract using nataive numpy.tensordot
    #        """
    #        l1= self.get_label_list()
    #        l2= T2.get_label_list()
    #        l_1n2=[i  for  i in l1 if i in l2]
    #        n=len(l_1n2)
    #        axis1=[l1.index(l_1n2[i]) for i in  xrange(n)]
    #        axis2=[l2.index(l_1n2[i]) for i in  xrange(n)]
    #        
    #        #print "labels:\t",l1,l2
    #        #print "axis to contract:\t",axis1, axis2
    #        data1= self.data.reshape(self.Dims[:self.rank],order='C')
    #        data2=T2.data.reshape(T2.Dims[:T2.rank],order='C')

    #        data3=np.tensordot(data1, data2 ,(axis1,axis2) )

    #        data3=data3.ravel("C")

    #        return data3


    #    def direct_product(self,T2):
    #        """
    #        status_1
    #        see nDirect_Product in f90
    #        """
    #        T1=self
    #        rank1 = int(T1.rank/2)
    #        rank2 = int(T2.rank/2)
    #        #attention_this_may_be_wrong
    #        rank3 = rank1+rank2
    #        nA=1; 
    #        mA=1
    #        for i  in xrange(rank1):
    #            Dims[i] = T1.Dims[i]    #for covariant legs, likewise nA below is xrange of covariant index of T1
    #            Dims[i+rank3] = T1.Dims[i+rank1]    #for contravariant legs,  likewise mA is ...
    #            nA=nA*T1.Dims[i]    
    #            mA=mA*T1.Dims[i+rank1]
    #        
    #        nB=1; mB=1
    #        for i  in xrange(rank2):
    #            Dims[i+rank1] = T2.Dims[i]   
    #            Dims[i+rank1+rank3] = T2.Dims[i+rank2]  
    #            nB=nB*T2.Dims[i]
    #            mB=mB*T2.Dims[i+rank2]
    #        T3=nTensor(rank3+rank3, Dims)

    #        Matrix_DirectProduct(T1.data, nA, mA, T2.data, nB, mB, T3.data)
    #        #attention_omitted_something
    #        return T3
    #class Test_nTensor(unittest.TestCase):         
    #    def setUp(self): 
    #        t_r1_a=nTensor(1,[5],['e'])
    #        t_r1_a.data=np.arange(t_r1_a.totDim,dtype='d')
    #        t_r1_b=nTensor(1,[5],['e']);        
    #        t_r1_b.data=np.arange(t_r1_b.totDim,dtype='d')
    #
    #        t_r2 = nTensor(2,[2,4], ind_labels=['a','d'])
    #        t_r2.data[:t_r2.totDim] = np.arange(np.prod(t_r2.Dims[:2]))
    #
    #        t_r3=nTensor(3, [2,3,6], ind_labels=['a','b','c'])
    #        t_r3.data[:t_r3.totDim] = np.arange(t_r3.totDim)    
    #        self.t_r3 = t_r3
    #
    #        t_r4=nTensor(4, [6,5,4,2], ind_labels=['c','e','d','a'])
    #        t_r4.data[:t_r4.totDim] = np.arange(t_r4.totDim)    
    #        self.t_r4 = t_r4
    #        
    #        t3=  nTensor(1,[3], ind_labels=['b']);  t3.data=np.arange(3,dtype='d')
    #        t4=  nTensor(1,[4], ind_labels=['d']);  t4.data=np.arange(4,dtype='d')
    #        t23=  nTensor(2,[2,3], ind_labels=['a','b']);  t23.data=np.arange(6,dtype='d')
    #        t24=  nTensor(2,[2,4], ind_labels=['a','d']);  t24.data=np.arange(8,dtype='d')
    #        t45=  nTensor(2,[4,5], ind_labels=['d','e']);  t45.data=np.arange(20,dtype='d')
    #        t234=  nTensor(3,[2,3,4], ind_labels=['a','b','d']);  t234.data=np.arange(24,dtype='d')
    #        t435=  nTensor(3,[4,3,5], ind_labels=['d','b','e']);  t435.data=np.arange(60,dtype='d')
    #        t435=  nTensor(3,[4,3,5], ind_labels=['d','b','e']);  t435.data=np.arange(60,dtype='d')
    #
    #    def test_nTensor_init():
    #        T=Tensor()
    #        Dims= [1,2,4,4]
    #        rank=4
    #        T=nTensor(rank=rank,Dims=Dims)
    #        print T
    #    
    #    def test_set_element(self):
    #        print "=============================== ="*3+"\n"
    #        print "test nTensor.set_element &.get_element  --- pass\n"
    #        rank = 2
    #        Dims= [3, 4]
    #        """A_{i1, i2, ... in} = A()   """
    #        T = nTensor(rank=rank, Dims=Dims)
    #        
    #        print T.data, '\n', T.data.reshape(Dims)
    #        ind1 = [2, 3]; ind2 = [1, 2] ; ind3 = [0, 0]
    #        T.set_element(ind1, 6)
    #        T.set_element(ind2, 2)
    #        T.set_element(ind3, 1)
    #        print T.data,  '\n', T.data.reshape(Dims)
    #        print T.get_element(ind1), T.get_element(ind2), T.get_element(ind3)
    #
    #        print " -------------------------------------------------------- "
    #        print "这里set, get都是前低后高, 而numpy reshape时相反, 是前高后低.下面illustrate this"
    #        
    #        rank = 3
    #        Dims= [3, 4, 5]
    #        t = nTensor(rank, Dims)
    #        t.data = np.arange(np.prod(Dims))
    #        print t.data[:t.totDim]
    #        data1 = np.arange(np.prod(Dims)).reshape(Dims)
    #        #print data1
    #        print "Dims:", Dims
    #        ind_vals = (1, 3, 2); ind_vals_1=ind_vals[::-1]
    #        print "the index is:\t", ind_vals ,  "\treverse ind:\t", ind_vals_1
    #        print "t.get_element([ind_vals]):", t.get_element(ind_vals)
    #        #print "t.data.reshape(Dims)"
    #        data2 = data1.T
    #        data3 = t.data.reshape(Dims[::-1])
    #        #print "reverse ind:\t", ind_vals_1
    #        print "shape of data2:", data2.shape  
    #        print "data2[ind_vals]:\t", data2[ind_vals], " data1[ind_vals_1]:\t", data1[ind_vals_1]
    #        print "data2[ind_vals_1]:\t", data2[ind_vals_1], "data1[ind_val]\t", data1[ind_vals]
    #        print "shape of data3:\t", data3.shape
    #        print "data3[ind_vals_1]:\t", data3[ind_vals_1]
    #        print "so,  the ansewer is : t.get_element(ind_val) = self.data.reshape(Dims[::-1])[ind_vals[::-1]]"
    #        print "while transpose of t.data.reshape is not relevent!!"
    #        #T.data = np
    #
    #    def test_permutation(self):
    #        print "=============================== ="*3+"\n"
    #        print "test nTensor.permutation ---pass\n"
    #        rank = 2
    #        Dims= [2, 4]
    #        t = nTensor(rank,Dims)
    #        t.data[:t.totDim] = np.arange(np.prod(Dims))
    #        data1 = np.arange(np.prod(Dims)).reshape(Dims)
    #        data2=np.swapaxes(data1, 0, 1)
    #        print 'data1:', data1.shape
    #        print "for a rank 2 tensor, is swapaxes just transpose? the answer is: ", 
    #        print (data1 == data2.T ).all()
    #        
    #        self.assertTrue((data1 == data2.T ).all())
    #        #t.Dims= [6, 4]
    #        print data1
    #        #print t.Dims[0:t.rank], '\t', t.data[0:t.totDim]
    #        print t 
    #        #p = np.random.permutation(t.rank)
    #        p = [1, 0]
    #        #p=[2,1]
    #        tp = t.permutation(p)
    #        print tp
    #    
    #    def test_contract(self):
    #        """a gross test ---pass """
    #        #V1=[50,60,30]
    #        #V2=[30,50, 80,100]
    #        #V3=np.empty(32,'i')
    #        V1=['a','b','c']
    #        V2=['c','e','d','a']
    #        V3=[None]*32  #['b','d','e']
    #        t_r3 = self.t_r3
    #        t_r4 = self.t_r4
    #        tc = t_r3.contract(t_r4, V1, V2, V3,  use_buf=None)
    #        print tc[0].data
    #    
    #    def test_contract_new_and_np():
    #        """test contract of two rank1 tensors---pass"""
    #        tc=t_r1_a.contract_new(t_r1_b)
    #        tc_np=t_r1_a.contract_np(t_r1_b)
    #        print "tc:\n",tc.data
    #        print "tc_np:\n",tc_np    #, type(tc_np.data), tc_np.data.shape
    #
    #        print "---pass"
    #        tc= t24.contract_new(t45)
    #        tc1= t24.contract_np(t45)
    #        print tc.data
    #        print tc1
    #
    #        print "---pass"
    #        tc= t24.contract_new(t435)
    #        tc1= t24.contract_np(t435)
    #        print tc.data
    #        print tc1
    #
    #        def subtest_mul_algorithm():
    #            print "test mul_algorithm lapack or numpy.dot yield same results--- pass"
    #            tc= t24.contract_new(t45)
    #            tca=t24.contract_new(t45,mul_algorithm="numpy")
    #            print tc.data
    #            print tca.data
    #        #subtest_mul_algorithm()
    #    
    #    def test_contract1():
    #        tc=t435.contract_new(t234,perm_order="C")
    #        tc1= t435.contract_np(t234)
    #        print tc.data
    #        print tc1
    #   
    #    def test_svd():
    #        """  --- pass"""
    #        m=5; n=3
    #        t2=nTensor(rank=2,Dims=[m,n])
    #        t2.data=np.random.rand(m*n).ravel()
    #        A=t2.data.reshape(m,n)
    #        print np.linalg.svd(A)[1]
    #        print t2.data
    #        print t2.svd(1)
    #
    #def test_nTensor():
    #    def test_Tensor():
    #        rank=4
    #        QSp=[QSp_base.copy()]*rank
    #        QSbase = QSp_base.copy
    #        totQN=QN_idendity.copy()
    #        u = Tensor()
    #        u.init_Tensor(rank,QSp,totQN)
    #        u.data[:] = np.arange(u.totDim)
    #        #print u
    #
    #        
    #        def test_scale():
    #            """-----------pass """
    #            print u.matrix_view(2)
    #            u1 = u.scale(3)
    #            print u1.matrix_view(2)
    #            u.scale(3, True)
    #            print u.matrix_view(2)
    #        
    #        def add_with_scale():
    #            """  ---pass """
    #            print u.matrix_view()
    #            w = test_iTensor.w.copy()
    #            uu= u.add_with_scale(u, 2, 1)
    #            print uu.matrix_view()
    #       
    #    def performance():
    #        import time
    #
    #        def performance_ntensor_contract(dim,n):
    #            """ 
    #            --- incomplete
    #            compare the perfomance of  f2py contract and np.tensordot
    #
    #            report:
    #
    #                当矩阵维数 >= 2**14时，matrix_multiply or np.dot all break down
    #            """
    #            #Dims=[8,3,4]
    #            Dims=[dim]*n
    #            rank=len(Dims)
    #            label1=range(rank)
    #            T=  nTensor(rank, Dims,  ind_labels=label1);  
    #            #label2 = xrange(rank-1, -1, -1)
    #            label2 = xrange(rank)
    #            T1= nTensor(rank, Dims,  ind_labels=label2);  
    #
    #            T.data=np.random.rand(np.prod(Dims))            
    #            T1.data=np.random.rand(np.prod(Dims))
    #
    #            t1=time.time()
    #            a=T.contract_new(T1,mul_algorithm="lapack", info=0)
    #            #a=3
    #            t2=time.time()
    #            print t2-t1
    #
    #            t3=time.time()
    #            b=T.contract_np(T1)
    #            #b=2
    #            t4=time.time()
    #            
    #            print "results:", a.data,b
    #            t21=t2-t1; t43=t4-t3
    #            print "fort: %2.5f"%t21, "numpy: %2.5f"%t43,'\t', int(float(t21)/(float(t43)+1E-50))
    #        #performance_ntensor_contract(2, 10)
    #
    #        def performance_ntensor_contract_1():
    #            """  
    #            every index has dim 2
    #            """
    #            for n in xrange(2,18,2):
    #                print "rank:", n, "\n","\t",
    #                performance_ntensor_contract(dim=2,n=n)
    #        #performance_ntensor_contract_1()
    #
    #        def performance_ntensor_contract_2():
    #            """  
    #            every index has dim 4
    #            """
    #            for n in xrange(2,14,2):
    #                print "rank:", n, "\n","\t",
    #                performance_ntensor_contract(dim=4,n=n)
    #        #performance_ntensor_contract_2()
    #
    #    def test_TensorLink():
    #        def test_init():
    #            tl = TensorLink(2)
    #            print tl.__dict__
    #    
    #    def test_iTensor_new():
    #        rank=3
    #        temp = QSp_base.copy
    #        QSp=[temp()]*rank
    #        QSbase = QSp_base.copy
    #        totQN=QN_idendity.copy()
    #        
    #        u=iTensor_new(rank,QSp,totQN)
    #        print u
    #        print u.data, 'shape', u.data.shape
    #        print np.where(u.data==0)
"""
about __new__: 
    __new__ is static class method, while __init__ is instance method.
    __new__ has to create the instance first, so __init__ can initialize it. 
        which can be a new one (typically that task is delegated to type.__new__), 
        an existing one (to implement singletons, "recycle" instances from a pool, and so on), 
        or even one that's not an instance of the class. If __new__ returns an instance of the class (new or existing), __init__ then gets called on it; if __new__ returns an object that's not an instance of the class, then __init__ is not called.    
    Note that __init__ takes self as parameter. Until you create instance there is no self
    
    T.__new__(S, ...) -> a new object with type S, a subtype of T
"""

def generalized_matrix_multiply(A, B, mul, A_nzb=None, B_nzb=None): 
    """
    """
    M = A.shape[0] 
    K = A.shape[1]
    N = B.shape[1]
    res = np.ndarray((M, N), dtype=np.object)
    if A_nzb is None: 
        for i in xrange(M): 
            for j in xrange(N): 
                t = 0.0
                for k in xrange(K): 
                    a = A[i, k]; b = B[k, j]
                    if isinstance(a, float) or isinstance(b, float): 
                        pass
                    else: 
                        #t  += updateCleft_new(self.Cl[0, j], A1.conj(), mpo4[1][j, i], A1)
                        t += mul(a, b) 
                res[i, j] = t 
    else: 
        raise NotImplemented('below is still not satisfactory')
        res[:, :] = 0.0
        assert B_nzb is not None 
        for i, k1 in A_nzb: 
            for k2, j in B_nzb: 
                if k1 == k2: 
                    res[i, j] += mul(A[i, k1],  B[k2, j]) 
    return res 

def contract_tensors(X, numindX, indX, Y, numindY, indY, out=None):
    """
        this func can be completely replaced by np.tensordot
        the out param is not useful yet
    """
    
    Xsize = X.shape
    Ysize = Y.shape
    
    indXl = range(numindX)
    indYr = range(numindY)
    #print_vars(vars(),  ['indX'])
    for i in indX:
        indXl.remove(i)
    for i in indY:
        indYr.remove(i)
        
    
    sizeXl = [Xsize[i] for i in indXl]
    sizeX = [Xsize[i] for i in indX]
    sizeYr = [Ysize[i] for i in indYr]
    sizeY = [Ysize[i] for i in indY]

    #names=["sizeXl","sizeX","sizeYr","sizeY","indXl","indYr"]
    #for n in names: print n, locals()[n] 
    if 0: 
        varables= vars()
        varables.update(Xshape=X.shape, Yshape=Y.shape)
        msg = """
        X.shape=%(Xshape)s, Y.shape=%(Yshape)s 
        sizeX=%(sizeX)s, sizeY=%(sizeY)s
        indX=%(indX)s, indY=%(indY)s
        """%varables
        print msg

    #if np.prod(sizeX) != np.prod(sizeY): 
    if not np.all(sizeX==sizeY): 
        varables= vars()
        varables.update(Xshape=X.shape, Yshape=Y.shape)
        msg = """
        error: indX and indY are not of same dimension: 
        X.shape=%(Xshape)s, Y.shape=%(Yshape)s 
        sizeX=%(sizeX)s, sizeY=%(sizeY)s
        indX=%(indX)s, indY=%(indY)s
        """%varables
        raise Exception(msg)
        

    if len(indYr)==0:  # if Y is completely contracted
        
        if len(indXl)==0: # if X is also completeley contracted 
            
            X=X.transpose(indX)
            X=X.reshape(-1, np.prod(sizeX))
            Y=Y.transpose(indY)
            Y=Y.reshape(np.prod(sizeY), -1)
            if out is None: 
                #Z = X.dot(Y) #contract to a scalar 
                Z = np.dot(X, Y)
                Zsize = 1
                return Z#, Zsize
            else: 
                np.dot(X, Y, out)
                return out#, Zsize
        else:
            
            X=X.transpose(indXl+indX)
            X=X.reshape((np.prod(sizeXl), np.prod(sizeX)))
            Y=Y.transpose(indY )
            Y=Y.reshape(np.prod(sizeY), np.prod(sizeYr))
            Zsize =  sizeXl
            if out is None: 
                Z = X.dot(Y)
                Z=Z.reshape(Zsize)
                #Z.ind_labels= indXl + indYr
                return Z
            else: 
                np.dot(X, Y, out)
                res = out.reshape((Zsize, 1)) 
                #res.ind_labels= indXl + indYr
                return res
    
    #print X.shape, indXl, indX,  indXl + indX
    X=X.transpose(indXl+indX)
    X=X.reshape((np.prod(sizeXl), np.prod(sizeX)))
    Y=Y.transpose(indY+indYr)
    Y=Y.reshape((np.prod(sizeY),np.prod(sizeYr)))
    #print X.shape,Y.shape

    Zsize = sizeXl + sizeYr
    if out is None: 
        Z = X.dot(Y)
        Z = Z.reshape(Zsize)
        #numindX = len(Z.shape)
        return Z  
    else: 
        np.dot(X, Y, out)
        res = out.reshape(Zsize)
        #res.ind_labels = indXl + indYr
        return res

def contract_tensors_new(X, Y, indX, indY): 
    """
        a wrapper and better interface for contract_tensors;
        this will be the standard one in future
    """
    if not hasattr(X, 'QSp'): 
        return contract_tensors(X, X.ndim, indX, Y, Y.ndim, indY)
    else:
        # transform of function interface 
        n = X.rank 
        m = Y.rank 
        V1 = np.arange(n, dtype=int)
        V2 = np.arange(n, n + m, dtype=int)
        a = - np.arange(1, 1 + len(indX))
        V1[indX] = a
        V2[indY] = a
        #return X.contract(Y, V1, V2, use_buf=True)[0]
        return X.contract(Y, V1, V2)[0]

class nTensor(np.ndarray, TensorBase):
    def __new__(cls, shape, dtype=float, buffer=None, offset=0,
          strides=None, order=None, info=None):
        # Create the ndarray instance of our type, given the usual
        # ndarray input arguments.  This will call the standard
        # ndarray constructor, but return an object of our type.
        # It also triggers a call to InfoArray.__array_finalize__
        obj = np.ndarray.__new__(cls, shape, dtype, buffer, offset, strides,
                         order)
        # set the new 'info' attribute to the value passed
        obj.info = info
        # Finally, we must return the newly created object:
        return obj
    
    def __init__(self, shape, dtype=float, buffer=None, offset=0,
          strides=None, order=None, info=None):
        pass 
        #TensorBase.__init__(self, self.ndim)
    
    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from
        # ndarray.__new__(InfoArray, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it -
        # i.e. those of a standard ndarray.
        #
        # We could have got to the ndarray.__new__ call in 3 ways:
        # From an explicit constructor - e.g. InfoArray():
        #    obj is None
        #    (we're in the middle of the InfoArray.__new__
        #    constructor, and self.info will be set when we return to
        #    InfoArray.__new__)
        if obj is None: return
        # From view casting - e.g arr.view(InfoArray):
        #    obj is arr
        #    (type(obj) can be InfoArray)
        # From new-from-template - e.g infoarr[:3]
        #    type(obj) is InfoArray
        #
        # Note that it is here, rather than in the __new__ method,
        # that we set the default value for 'info', because this
        # method sees all creation of default objects - with the
        # InfoArray.__new__ constructor, but also with
        # arr.view(InfoArray).
        self.info = getattr(obj, 'info', None)
        # We do not need to return anything

    def contract(self, other, ind_list_1, ind_list_2): 
        """
            althouth I defined ndarray in contract_tensors func, it returns an instantce of nTensor 
        """
        return contract_tensors(self, self.ndim, ind_list_1, other, other.ndim, ind_list_2)

class Test_nTensor(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_temp(self): 
        a = nTensor((3, 3, 3), int)
        print(a)
        print  super(a.__class__) 
    
    def test_contract(self): 
        a = nTensor((2, 2, 2))
        b = nTensor((2, 2, 2))
        c=a.contract(b, [0, 1, 2], [0, 1, 2])
        print  type(c)
        
        

def test_TensorBase():
    t = TensorBase()
    def test_init():
        t=TensorBase()
    def test_randomize_data():
        #mrandom.csrand(1)
        t.totDim = 10
        t.data = np.arange(t.totDim, dtype=t.dtype)
        t.randomize_data()
        print t.data
    #test_randomize_data()

        
if __name__ == '__main__' : 


    if 0: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #Test_nTensor('test_temp'), 
           Test_nTensor('test_contract'), 
           
        ]
        for a in add_list: 
            suite.addTest(a)
        unittest.TextTestRunner().run(suite)
       

