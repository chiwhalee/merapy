#!/usr/bin/env python
#coding=UTF8


"""
    the name of nTensor  can be interpreted as normal ternsor or tensor that inhereitated
    from np.ndarray 
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
from builtins import object
import unittest 
import numpy as np 
import numbers 
from numbers import Number
import pickle as pickle 

#from merapy import crandom
from merapy.utilities import print_vars

class TensorBase(object):
    """
    an abstract tensor base class 
    this class is newly added 
    """
    MaxRank= 32  # maximal rank of tensors supported
    #dtype = float 

    def __init__(self,rank=None, ind_labels=None,use_buf=False):
        """
        
        ind_labels: dict
        type_name: may be "U", "V", "V2" etc
        """
        self.rank=rank
        self.ind_labels= {i:None for i in range(rank)}
        if ind_labels is not None:
            #print 'iii', ind_labels
            
            for i in range(rank):
                print('iii', i, ind_labels)
                self.ind_labels[i] = ind_labels[i]
        self.data = None
        self.type_name = ""

    def set_ind_labels(self, ind_labels):
        """
        ind_labels:  a dict whose keys are from [0,rank-1], values arbitrary
        """
        self.ind_labels.update(ind_labels)
    
    def get_label_list(self):
        res= [self.ind_labels[i] for i in range(self.rank)]
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
        #reason= self.is_same_shape(T2)
        reason = 1
        if  reason<0:  
            print('Error in Tensor_Add_WithScale, self,other different shape:', reason)
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

    def __neg__(self):
        res= self.copy() 
        res.data = -self.data 
        return res 
    
    def __iadd__(self, other):
        self.data += other.data 
        return self

    def __add__(self, other):
        res= self.copy()
        res.data = self.data + other.data 
        return res
    
    def __radd__(self, other):
        res = self.copy()
        #if isinstance(other, float): 
        if isinstance(other, Number): 
            res.data = self.data + other
        else: 
            res.data = self.data + other.data 
        return res


    def __sub__(self, other):
        res = self.copy()
        #if isinstance(other, Number): #this is erro pron
        #    raise  
        #    res.data = self.data - other
        #else:
        #    res.data = self.data - other.data
        res.data = self.data - other.data
        return res

    def __rmul__(self, scalar):
        """
            calculate scalar*self 
        """
        #self.data *= scalar  # this is wrong!,  because this is not intend to be inplace operation
        #return self
        res = self.copy_struct(has_data=False)
        res.data = scalar*self.data
        # I dont use the following,  because if scalar is of type complex, numpy raise if self.data is not of type complex 
        #res = self.copy()
        #res.data *= scalar  
        return res 
    
    def __truediv__(self, scalar):
        res = self.copy()
        res.data /= scalar 
        return res 
    
    def __itruediv__(self, scalar):
        self.data /= scalar 
        return self
    
    def __div__(self, scalar):
        res = self.copy()
        res.data /= scalar 
        return res 
    
    def __imul__(self, scalar):
        self.data *= scalar
        return self

    def add_with_scale_new(self,T2, alpha,beta):
        """
        bug in it
        see Tensor_Add_WithScale
        """
        #same, reason= self.is_same_shape(T2)
        same = 1
        
        if  not same:  
            print('Error in Tensor_Add_WithScale, self,other different shape:', reason)
            print('self.' + reason + ":", end=' ') 
            print(self.__getattribute__(reason))
            print('other.' + reason + ":", end=' ') 
            print(T2.__getattribute__(reason))
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
        for i in range(len(S)):
            #if InSet[iLeg, S[i]]>0:  
            if iLeg in S[i]:
                #Is[k] = i
                Is.append(i)
                #k = k+1
        if len(Is)>2:
            print("error, see find_leg")
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

def generalized_matrix_multiply(A, B, mul, null=0.0,  A_nzb=None, B_nzb=None): 
    """
        used to multiply two mpo's which are block_wise
        params:
            mul: a generalized 'multiply' operation, which can be defined to be any 
                two element operations!
        returns: 
            C = A*B 
    
    """
    M = A.shape[0] 
    K = A.shape[1]
    N = B.shape[1]
    res = np.ndarray((M, N), dtype=np.object)
    if A_nzb is None: 
        for i in range(M): 
            for j in range(N): 
                t = null 
                for k in range(K): 
                    a = A[i, k]; b = B[k, j]
                    if isinstance(a, float) or isinstance(b, float): 
                        pass
                    else: 
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
    
    indXl = list(range(numindX))
    indYr = list(range(numindY))
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
        print(msg)

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
            Y=Y.reshape(np.prod(sizeY, dtype=np.int), np.prod(sizeYr, dtype=np.int))
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
    X=X.reshape((np.prod(sizeXl, dtype=np.int), np.prod(sizeX, dtype=np.int)))
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

def contract_tensors_new(X, Y, indX, indY, use_buff=False): 
    """
        a wrapper and better interface for contract_tensors;
        this will be the standard one in future
    """
    if not hasattr(X, 'QSp'): 
        #raise  # remove below 
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
        return X.contract(Y, V1, V2, return_v3=True, use_buf=use_buff, )[0]

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
        print(super(a.__class__)) 
    
    def test_contract(self): 
        a = nTensor((2, 2, 2))
        b = nTensor((2, 2, 2))
        c=a.contract(b, [0, 1, 2], [0, 1, 2])
        print(type(c))
        
        

def test_TensorBase():
    t = TensorBase()
    def test_init():
        t=TensorBase()
    def test_randomize_data():
        #mrandom.csrand(1)
        t.totDim = 10
        t.data = np.arange(t.totDim, dtype=t.dtype)
        #t.randomize_data()
        print(t.data)
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
       

