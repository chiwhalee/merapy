#! /usr/bin/env python
#coding=UTF8

"""
    issue:
        performance issue: may need to change qn val to np.int16 or np.int8, this may be faster? 
        
    some math
        In representation theory, QnZ2, QnU1 etc defined below are in fact the
        so called "WEIGHT SPACE" \Lambda={\lambda_1, \lambda_2, ...}.
            这些类中的operator实际上给出了 "FUSION
            RULE"，即给出多个矢量空间的张量积的量子数的运算法则(i.e. define an
            algebra)
        
        While QspZ2,  QSpU1 etc define a vector space V, which is direct sum
        decomposition according to the weights: V = \osum_{i}{blablabla...}. 
        note this osum is ORDERED
    
    todo:
        1. 用descriptor
        2. MaxQNNum 是fortran的做法，在python中可以改成可变长度的
        3. Qn 和QSp 都应该 弄成 immutable type, 像 integer，tuple那样。
            mutable type使得不得不 到处都得copy，很麻烦。这是长时间使用的经验感受.
            imutable 就可一大胆放心、随意第使用.
            这要求，所有对其的操作都是 not-inplaced 
            
    Qsp
        a Qsp is a (reducible) L-module for lie algebra L
    QuantumNumber 定义了对于给定李代数表示的root(or weight?)集合，以及fusion rule

    Z2: spin parity invariance

    as for Z2 (parity) symmetry, there is no difference between a vector and its dual, 
    that is to say need not distinguish upper and lower legs
"""

#from past.utils import old_div
from builtins import object
import unittest 
#import numpy as np
import warnings
from abc import (ABCMeta, ABC, abstractmethod)
from collections.abc import Iterable, Sized
from typing import (Any, Tuple, List, Type, Callable,  Generic, TypeVar, Callable, get_type_hints, 
        Optional, ClassVar, no_type_check, overload, Union)

#from py3nj import (clebsch_gordan, wigner, wigner3j, wigner6j, wigner9j)


#from merapy.utilities import print_vars
#from merapy.decorators import decorate_methods, tensor_player 

__all__=[  "QuantSpaceBase",  
    "QnU1", "QnZ2", "QnZ3", "QnTravial", 
    "QspU1", "QspZ2", "QspZ3", "QspTravial", 'make_qsp', 'qsp_any',
     "symmetry_to_Qn", "symmetry_to_Qsp", 'symmetry_to_QspClass', "qn_factory"]

GROUP_NAMES = ["Travial", "Z2", "Z3", "U1"]


Ta = TypeVar('Ta', bound='A') 
class A:
    def __init__(self, i:int)->None:
        self.i = i
    
    def __add__(self:Ta, other:Ta)->int:
        return self.i + other.i


Tqn = TypeVar('Tqn', bound='QnBase') 

#Q = TypeVar('Q', bound=QnBase)

#class QnBase(object):
class QnBase(ABC):   #with object new python cannot load old pickled files,  I dont know why 
    #__metaclass__ = ABCMeta
    """
        abstract base class, not meant for direct use
    """
    
    #SYMMETRY:ClassVar[str] = "NaN"   #add this line only for passing compile of mypyc
    IS_Abelian = True
    QnId:ClassVar[Union[int, None]] = None
    SYMMETRY:ClassVar[Union[str, None]] = None
    #def __init__(self):
    #    pass
    #def __init__(self, value):
    #    pass
    #def __getstate__1(self) -> int:
    #    """ not workable"""
    #    return self._val
    #def __setstate__1(self, state) -> int:
    #    self._val = state["_val"]
        
    @property
    def val(self) -> int:
        return self._val

    #@val.setter #it seems that setter effects only when QnBase inherite from object
    #def val(self, new_val):
    #    #print 
    #    raise NotImplemented("this don't work, because of olb-style class")
    #    self._val = new_val
    
    def __repr__(self)->str:
        """
        newly added
        """
        #res= repr(self.val)
        #res="(" + "%5d"*self.num_of_symm +")"  % tuple(self.val)
        #if isinstance(self.val, np.ndarray):
        #if isinstance(self._val, Iterable):
        #    res=("(" + "%5d"*len(self._val)+")")  % tuple(self._val)
        if 0:
            pass
        elif isinstance(self._val, int):
            res= ("(" + "%5d"")")  % self._val
        elif isinstance(self, QnSU2):
            res = '({}, {})'.format(self._val[0], self._val[1])
        else:
            raise ValueError(type(self._val))
        return res
    
    def __eq__(self, other):
        return self._val == other._val 
    
    def __ne__(self, other):
        return self._val != other._val 
    
    def __lt__(self, other): 
        return self._val < other._val 
    
    def __hash__(self):
        return self._val 
    
    def set_val(self,val):
        """
        newly added, assign val to self.val
        val: an instance of QuantumNum  or   an integer or an iterable
        """
        #if isinstance(val, iterable):
        if hasattr(val, "__iter__"):
            self._val = val[0]
        elif isinstance(val, int):
            self._val = val
        else:
            raise ValueError("type of val not correct: %s"%type(val))
    
    @classmethod
    def qn_id(cls):
        return cls(cls.QnId)
    
    @classmethod
    def sum(cls, qn_list): 
        res = qn_list[0]
        for i in qn_list[1: ]: 
            #res= res + i
            res= res.__add__(i)
        return res
    
    def qsp_class(self): 
        return symmetry_to_Qsp(self.SYMMETRY)
    
    def conj(self): 
        res = self.copy()
        res.reverse()
        return res 
    
    def copy(self):   
        return self.__class__(self._val)
    
class QnTravial(QnBase):
    SYMMETRY:ClassVar[str] = "Travial"
    NUM_OF_SYMM = 1  #this may be an issue NUM_OF_SYMM should be 0?
    QNS = (1, )
    QnId:ClassVar[int] = 1
    def __init__(self, value=1) -> None:
        """
            the param value is just for consistency with other QnClass 
        """
        self._val = 1

    def set_val(self, val):
        """
            _val always is 1
        """
        self._val = 1

    def copy(self):
        return QnTravial()
    
    def reverse(self):
        pass # do nothing
    
    def __add__(self, other:'QnTravial')->'QnTravial':
        return self 
    
    def __eq__(self, other):
        """
            always True
        """
        if self._val != other._val:
            raise Exception
        return True
    
    def __ne__(self, other):
        """
            always False since only one qn 1
        """
        if self._val != other._val:
            raise Exception
        return False

    @classmethod
    def qn_id(cls):
        return QnTravial()

class QnZ2(QnBase):
    """ 
        自旋朝下粒子数的即偶性
        spectrum  of down-spin-parity operator (-1)^(-sz) is {1, -1}
        the spectrum itself forms an algebra isometric to Z2
        this calss define the operations for this algebra
        这一定义，即采用
        CONVENTION: 
            对于单个自旋,1, -1 分别对应于朝上和朝下
            
    """
    SYMMETRY:ClassVar[str] = "Z2"   #precisely, this is used for spin-parity symm
    NUM_OF_SYMM = 1
    QNS = (1, -1)   #all the elements in the algebra, 
    QnId:ClassVar[int] = 1 # identity in the algebra
    
    def __init__(self, value:int) -> None:
        """
            把val统一写成np.ndarray类型, 这样有利有弊: 利在可用array的各种运算，比如*法；弊在赋值麻烦些
            num_of_symm: newly added atribute by lzh, 
        """
        
        self._val = value
    
    def copy(self):
        return QnZ2(self._val)  #if use the following line, program will grow slower and slower, I don't know why yet
        
    def reverse(self):
        """
            reverse refers to conjugate repr. of a group
        """
        pass


    def __add__(self, other:'QnZ2')->'QnZ2':
        """
            see QN_Add in f90
            q: 所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnZ2(self._val*other._val)

    @classmethod
    def qn_id(cls):
        return QnZ2(1)

class QnZ3(QnBase):
    SYMMETRY:ClassVar[str] = "Z3"
    NUM_OF_SYMM = 1
    #QNS = (1, -1)
    QNS = (0, 1, 2)
    QnId:ClassVar[int] = 0 # identity in the algebra
    
    def __init__(self, value:int) -> None:
        """
            把val统一写成np.ndarray类型, 这样有利有弊: 利在可用array的各种运算，比如*法；弊在赋值麻烦些
            num_of_symm: newly added atribute by lzh, 
        """
        self._val = value
    
    def copy(self):
        return QnZ3(self._val)
        
    def reverse(self):
        """
            reverse refers to conjugate repr. of a group
            see QN_Reverse in f90
            0-->0
            1-->2
            2-->1
        """
        pass
        self._val = (3-self._val)%3

    def __add__(self, other:'QnZ3')->'QnZ3':
        """
            see QN_Add in f90
            q: 所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnZ3((self._val + other._val)%3)
    
    def __reduce__del(self, other):
        return QnZ3(self.val*other.val)

class QnU1(QnBase):
    SYMMETRY:ClassVar[str]= "U1"
    NUM_OF_SYMM = 1
    QNS = tuple(range(-10, 11))
    QnId:ClassVar[int] = 0 # identity in the algebra
    
    def __init__(self, value:int) -> None:
        """
            把val统一写成np.ndarray类型, 这样有利有弊: 利在可用array的各种运算，比如*法；弊在赋值麻烦些
            num_of_symm: newly added atribute by lzh, 
        """
        self._val = value
    
    #def __new__x(cls, value):
    #    self = super(QnU1, cls).__new__(cls, value)
    #    self._val = value
    #    return value
    
    def copy(self):
        return  QnU1(self._val)
        
    def reverse(self):
        """
            reverse refers to conjugate repr. of a group
        """
        self._val = -self._val
    
    
    def __add__(self, other:'QnU1')->'QnU1':   #def __add__(self, other:'QnU1')->Type['QnU1']:
        """
            所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnU1(self._val + other._val)

    def add(self, other:'QnU1')->'QnU1':
        """
            所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnU1(self._val + other._val)

    

    def conj(self): 
        return QnU1(-self._val)
    
    @classmethod
    def qn_id(cls):
        return QnU1(0)

class QnSU2(QnBase):
    SYMMETRY:ClassVar[str] = 'SU2'
    NUM_OF_SYMM = 1
    QnId:ClassVar = (0, 0)    # (j, m)
    
    IS_Abelian = False
    #def __init__(self, value):
    def __init__(self, value:Tuple[int, int]) -> None:
        """
        """
        self._val = value
    def conj(self):
        return QnSU2((self._val[0], -self._val[1]))
    
    def __add__(self, other):
        """
        """
        return QnSU2((self._val[0] + other._val[0], self._val[1] + other._val[1]))

    def reverse(self):
        self._val = (self._val[0], -self._val[1])
    
    #commented only for debuging mypy
    #def cg_coeff(qn0, qn1, qn2):
    #    return clebsch_gordan(qn0._val[0],  qn1._val[0], qn2._val[0], 
    #            qn0._val[1],  qn1._val[1], qn2._val[1])

type_QnU1 = TypeVar('type_QnU1', bound=QnU1)

T = TypeVar('T', bound='QuantSpaceBase')


#meth_names= ["copy"]
#@decorate_methods(tensor_player, meth_names)
class QuantSpaceBase(object):
    """
        abstract base class
        so MaxQNNum only defined in subclassed
        issue:
            method may be slow:
                __eq__, Dims , update, tensor_prod
        discuss:
            to define self._dims as a list of np.ndarray ?
            It appears that the former is better, reason:
                1. The number of qn is not limited now,  
                    Then no need let all qn cetered around qn_id which
                    reqire me must shift qn when defining a mps with nontrival
                    totqn. I can now expand range of qn as wide as I want. This is quite critical for a  simpler code. 
                    so it is much more flexible and convenient. 
                2. may be even faster, as the tests show! Because no need to allocate memory when construct a qsp!
                3. save memory
            
    """
    IS_Abelian:ClassVar[bool] = True
    #QnClass:Type = QnBase
    #QnClass:ClassVar[Type['QnBase']] = QnBase
    QnClass:ClassVar
    #def __init__(self, n, qns, dims) -> None:
    def __init__(self, n:int, qns:List, dims:List[int])->None:
        """
            QSp真正有用的只有 QNs, _dims 两个属性
            _dims: 每个量子数对应子空间的维数
            
            note:
                1. is it better to change self._dims to be varying length 
                2. is it better to store self._dim in numpy array or list?
                some refs:
                   http://stackoverflow.com/questions/46860970/why-use-numpy-over-list-based-on-speed0       
 
        """
        self.nQN = n
        #self._dims = np.empty(self.MaxQNNum, int) #MaxQNNum only defined in subclassed
        self._dims = dims
        self._totDim = -1 
        self.QNs = qns
    
    @property
    def totDim(self) -> int:
        """
            lazy evaluation
        """
        #if hasattr(self, "_totDim"):
        #    return self._totDim
        if self._totDim != -1 :
            return self._totDim
        else:
            #self._totDim = np.sum(self._dims[:self.nQN])  
            self._totDim = sum(self._dims[:self.nQN])     #np.sum is no faster than builtin sum when acting on normal list of ints
            return self._totDim

    @property
    def tot_dim(self) -> int:  # replace tot_Dim in future 
        """
            lazy evaluation
        """
        #if hasattr(self, "_totDim"):
        #    return self._totDim
        if self._totDim is not None:
            return self._totDim
        else:
            #self._totDim = np.sum(self._dims[:self.nQN])  
            self._totDim = sum(self._dims[:self.nQN])     #np.sum is no faster than builtin sum when acting on normal list of ints
            return self._totDim
    
    @property
    @abstractmethod
    def symmetry(self) -> Optional[str]: 
        return self.QnClass.SYMMETRY
    
    @property
    def Dims(self) -> List[int]:
        return self._dims
    
    @Dims.setter
    def Dims(self, new_dims):
        self._dims[:len(new_dims)] = new_dims
        self._totDim = sum(self._dims)

    def __ge__(self, other) -> bool: 
        if self.QnClass != other.QnClass: 
            return False
        if self.nQN < other.nQN: 
            return False
        return not self < other
        
    def __le__(self, other) -> bool:
        info=0  #mypyc not allow pass 3 args to __le__
        if self.QnClass != other.QnClass: 
            if info>0: print('__le__ reason: ', 1)
            return False
        if not self.nQN <= other.nQN:
            if info>0: print('__le__ reason: ', 2)
            return False
        #if not np.all(self._dims[:self.nQN]<=other._dims[:self.nQN]):
        if not all(self._dims[i]<=other._dims[i] for i in range(self.nQN)):
            if info>0: print('__le__ reason: ', 3, self._dims[:self.nQN], other._dims[:self.nQN])
            
            return False
        return True
    
    def __lt__(self, other): 
        if self.QnClass != other.QnClass: 
            return False
        if not (self <= other and self != other):  
            return False
        return True
      
    def __gt__(self, other): 
        """
            由 le 派生
        """
        if self.QnClass != other.QnClass: 
            return False
        if self <= other:
            return False
        if self.nQN != other.nQN:   #neither completely <, >,  = ; not comparable
            return False
        return True

    def __eq__(self, other) -> bool:
        """
        """
        if self.nQN != other.nQN:
            return False
        if self.totDim != other.totDim:
            return False
            
        nQN = self.nQN
        #if not np.all(self.QNs[:nQN]==other.QNs[:nQN]):
        if not all([self.QNs[i] == other.QNs[i] for i in range(nQN)]):
            return False 
        #if not np.all(self._dims[:nQN]==other._dims[:nQN]):
        if not all([self._dims[i] == other._dims[i] for i in range(nQN)]):
            return False 
        return True

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def __repr__(self) -> str:
        """ 
       
        """
        keys= ["class", 'nQN', 'QNs', '_dims', ]
        res= ""
        n = self.nQN
        temp = []
        for i in range(n): 
            val = self.QNs[i]._val
            dim = self._dims[i]
            s= '[%s]%d'%(val,dim)
            temp.append(s)
        res = '+'.join(temp)

        return res
    
    def __mul__(self, other):    #not supported by mypyc, so I comment it
        #return self.add(other)
        return self.tensor_prod(other)
    
    #def __div__(self, other): 
    #    """
    #        q1=QspZ2.easy_init([1, -1], [2, 2])
    #        q2=QspZ2.easy_init([1, -1], [1, 3])
    #        print  q1*q2   == q1*q1 
    #        #they are equal 
    #    """
    #    msg = 'divide of qsp cant be defined,  as it is not unique. see doc string'
    #    raise NotImplemented(msg)
    
    def __pow__(self, n:int):  #not supported by mypyc, so I comment it
        res= self.__class__.null()
        for i in range(n): 
            res= res.__mul__(self)
        return res 
    
    @classmethod
    #def easy_init(cls:Type[T], qns=None, dims=None)->T:
    def easy_init(cls:Type[T], qns=Optional[List[Tqn]], dims=Optional[List[int]])->T:
        """ a slow but easy init """
        symm =  cls.QnClass.SYMMETRY 
        if symm== 'Travial':
            qns1 = [QnTravial()]
        else:
            qns1 = [cls.QnClass(i) for i in qns]
        n = len(qns1)
        dims = list(dims)
        return cls(n=n, qns=qns1, dims=dims) 

    def reverse(self)->None:
        """
        """
        for i in range(self.nQN):
            self.QNs[i].reverse()
    
    def conj(self):
        """
            not inplace reverse 
        """
        res= self.copy()
        res.reverse()
        return res 
    
    #@tensor_player(which="Qsp_copy")
    def copy(self:T, reverse=False)->T:  # see here for the neccisity of type hint T https://www.python.org/dev/peps/pep-0484/#id34
        qns= [q.copy() for q in self.QNs[:self.nQN]]
        #other = self.__class__(n=self.nQN, qns=qns, dims=self._dims)
        other = self.__class__(n=self.nQN, qns=qns, dims=self._dims.copy())
        other._totDim = self.totDim
        if reverse: 
            other.reverse()
        return other
    
    def copy_from(self, other):
        self.nQN = other.nQN
        self._totDim=other.totDim
        nqn=other.nQN
        if nqn != 0:
            self.QNs[:nqn]=[other.QNs[i].copy() for i in range(nqn)]
        self._dims=other._dims.copy()

    def update(self, qsp_max=None):
        """
            this method only has limited use, and in non critical block 
        """
        if qsp_max is not None:
            nqn = 0
            QNs= []
            _dims= [0 for i in range(self.MaxQNNum)]
            
            for n in range(self.nQN):
                i= qsp_max.has_quant_num(self.QNs[n])
                if i>=0:
                    QNs.append(self.QNs[n])
                    _dims[nqn]=min(self._dims[n], qsp_max._dims[i] )
                    nqn += 1 
            temp = self.__class__(nqn, QNs, _dims)
            self.copy_from(temp)

    def copy_many(self, n, reverse=None):
        res=[self.copy() for i in range(n)]
        if reverse is not None:
            for i in reverse:
                res[i].reverse()
        return res

    def has_quant_num(self, qn):
        """
            qn: instance of Quantum Num  or a tuple
        """
        try:
            return self.QNs.index(qn)
        except:
            return -1 
    has_qn = has_quant_num   #def has_qn

    def add_to_quant_space(self, qn, d):
        """
            this is in fact direct sum of vector spaces subject to symmetry
            note that direct sum is not a commutable operation
        """
        i = self.has_quant_num(qn)
        assert i<self.MaxQNNum, self.MaxQNNum
        if i<0:  #when qn not in self.QNs
            self.QNs.append(qn)
            self._dims.append(d)
            self.nQN += 1
        else:  # when qn in self.QNs 
            self._dims[i] += d  
        self._totDim = self.totDim  + d

    def tensor_prod(self:T, other:T)->T:
        """ 
             although named add, acturally tensorprod of self and other
            todo: change the name to prod in future 
            note:
                after tensor prod, the value of the qn in the qsp is not ordered. 
        """
        if self.nQN==0:
            return other.copy()    #need copy here?
        if other.nQN==0:
            return self.copy()     #need copy here?
    
        res = self.__class__(n=0, qns=[], dims=[])
        res._totDim = 0

        for i in range(self.nQN):
            for j in range(other.nQN):
                #qn = self.QNs[i] + other.QNs[j]
                qn = self.QNs[i].__add__(other.QNs[j])
                d = self._dims[i] * other._dims[j]
                res.add_to_quant_space(qn, d)
        return res
    
    #add = tensor_prod  #def add
    
    @staticmethod
    def prod_many( qsp_list): 
        q = qsp_list[0]
        for a in qsp_list[1: ]: 
            q = q.tensor_prod(a)
        return q 
    
    @classmethod
    #def null(cls) ->Type[QuantSpaceBase]:
    #def null(cls)->Union['QspZ2', 'QspU1', 'QspTravial', 'QuantSpaceBase']:  #invaled
    #def null(cls) -> 'QuantSpaceBase':
    def null(cls:Type[T]) -> T:
        """
            a trivial (1D) vector space
        """
        qn = cls.QnClass.qn_id()
        return  cls(n=1, qns=[qn], dims=[1])
    
    @classmethod
    def empty(cls):
        return cls(n=0, qns=[], dims=[])
    
    def expand_totdim(self, totdim): 
        qns = self.QNs
        qns = [q.val for q in qns]
        dims = self.Dims[:self.nQN].copy() 
        
        #ratial = dims/float(self.totDim)
        ratial = [d/float(self.totDim) for d in dims ]
        #dims = totdim * ratial
        dims= [totdim*d for d in ratial]
        dims = [int(i) for i in dims]
        diff = sum(dims)-totdim
        if diff >0 or diff >self.nQN: 
            raise
        else: 
            for i in range(-diff): 
                dims[i] += 1 
        
        res= self.__class__.easy_init(qns, dims)
        return res 
    
    def expand(self, other):
        for i in range(other.nQN):
            qn = other.QNs[i]
            dim = other.Dims[i]
            self.add_to_quant_space(qn, dim)
    
    @classmethod
    def qn_id(cls:Type[T])->Tqn:
        return cls.QnClass.qn_id()
    

class QspTravial(QuantSpaceBase):
    MaxQNNum = 1
    #QnClass:Type[QnBase]= QnTravial
    QnClass:ClassVar = QnTravial
    def __init__(self, n, qns, dims) -> None:
        #self.nQN = n
        #self._dims[:n] = dims[:n]
        #self._totDim = None
        #self.QNs = qns
        QuantSpaceBase.__init__(self, n, qns, dims)

    @classmethod
    def set_base(cls, dim=None):
        if dim is None:
            dim = 2
        qsp_base = QspTravial(n=1, qns=[cls.QnClass.qn_id()], dims=[dim])
        
        #qn_identity = QnTravial()
        qn_identity = cls.QnClass.qn_id()
        
        qsp_null = QspTravial.null()
        return qn_identity, qsp_base, qsp_null 
    
    @classmethod
    def max(cls, trunc_dim, nqn=None):
        qsp_max = cls(n=1, qns=None, dims=[trunc_dim])
        if 1:
            qsp_max2 = cls(n=1, qns=None, dims=[trunc_dim])
        return qsp_max
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        assert len(dims)==1
        return cls(n=1, qns=[QnTravial()], dims=dims) 
    
    def shift_qn(self, qn=None):
        """
            just do nothing 
        """
        pass
    
class QspZ2(QuantSpaceBase):
#class QspZ2(Generic[T]):
    MaxQNNum = 2
    #QnClass:Type = QnZ2
    QnClass:ClassVar[Type['QnZ2']] = QnZ2
    def __init__(self, n, qns, dims) -> None:
        """

        """
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
    
    #def __div__(self, other): 
    #    """
    #        solve this linear eqn of (d0, d1): 
    #            d0*b0 + d1*b1 = c0
    #            d0*b1 + d1*b0 = c1
    #        -->
    #            d0*b0*b1 + d1*b1^2 = c0*b1
    #            d0*b1*b0 + d1*b0^2 = c1*b0
    #        -->
    #            d1 = (c0*b1 - c1*b0)/(b1^2 - b0^2)
    #            d0 = (c1*b0 - c0*b1)/(b1^2 - b0^2)
    #            
    #            d0, d1 may be not unique 
    #            
    #    """
    #    raise NotImplemented
    #    
    #    #b0 = other._dims[0]
    #    #b1 = other._dims[1]
    #    #
    #    #c0 = self._dims[0]
    #    #c1 = self._dims[1]
    #    #res= self.__class__.easy_init([1, -1], [d0, d1])
    #    
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        n = len(dims)
        assert dims is not None 
        if qns is None: 
            qns = [1, -1]
            assert len(qns)>= len(dims)
            qns = qns[: len(dims)]
        qns1 = [cls.QnClass(i) for i in qns]
        return cls(n=n, qns=qns1, dims=dims)

    @classmethod
    def set_base(cls, dim=None):
        if dim is None:
            dim = 2
        #qn_identity = QnZ2(1)
        qn_identity = cls.QnClass.qn_id()

        dim_half = dim//2
        qns1 = [1, -1]
        qns = [QnZ2(i) for i in qns1]
        dims = [dim_half, dim_half]
        qsp_base = QspZ2(n=2, qns=qns, dims=dims)

        qsp_null = QspZ2.null()
        
        if 0:
            num_good_qn = 2
            goodqns= [QnZ2(1) for i in range(num_good_qn)]
            for i in range(num_good_qn):
                goodqns[i].val[0] = [1,-1][i]

        return qn_identity, qsp_base, qsp_null
    
    @classmethod
    def max(cls, trunc_dim, nqn=None):
        if trunc_dim%2 != 0:
            raise ValueError("trunc_dim should be multipy of 2.  trunc_dim=%d"%trunc_dim)
        dim_half = trunc_dim//2
        dims = [dim_half, dim_half]  
        qns1 = [1, -1]
        qns=[QnZ2(i) for i in qns1]
        qsp_max = QspZ2(n=2, qns=qns, dims=dims)

        if 1:
            dims= [2, 2]
            qns=[QnZ2(i) for i in qns1]
            qsp_max2 = QspZ2(n=2, qns=qns, dims=dims)
            #qsp_max2.RefQN[0,0:2] = [2,1]
            #qsp_max2.RefQN[1,0:2] = [0,1]
            qsp_max2.update()
        #return qsp_max, qsp_max2
        return qsp_max
   
    #@classmethod
    #def qn_id(cls)->QnZ2:
    #    return QnZ2(1)
    
    #@classmethod
    #def null(cls) -> 'QspZ2': #def null(cls) ->Type[QuantSpaceBase]:  #def null(cls)->Union['QspZ2', 'QspU1', 'QspTravial', 'QuantSpaceBase']:  #invaled
    #    """
    #        a trivial (1D) vector space
    #    """
    #    qn = QnZ2.qn_id()
    #    return  cls(n=1, qns=[qn], dims=[1])
   
class QspZ3(QuantSpaceBase):
    MaxQNNum = 3
    QnClass:ClassVar = QnZ3
    def __init__(self, n, qns, dims) -> None:
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)

    @classmethod
    def set_base(cls, dim=None):
        if dim is None:
            dim = 3
        #qn_identity = QnZ3(1)
        qn_identity = cls.QnClass.qn_id()

        dim_over_3 = dim//3
        qns1 = [0, 1, 2]
        qns = [QnZ3(i) for i in qns1]
        dims = [dim_over_3]*3
        qsp_base = QspZ3(n=3, qns=qns, dims=dims)

        qsp_null = QspZ3.null()
        
        return qn_identity, qsp_base, qsp_null
    
    @classmethod
    def max(cls, trunc_dim, nqn=None):
        if trunc_dim%3 != 0:
            raise ValueError("trunc_dim should be multipy of 3.  trunc_dim=%d"%trunc_dim)
        dim_over_3 = trunc_dim//3
        dims = [dim_over_3]*3  
        qns1 = [0, 1, 2]
        qns=[QnZ3(i) for i in qns1]
        qsp_max = QspZ3(n=3, qns=qns, dims=dims)

        
        return qsp_max

class QspU1(QuantSpaceBase):
    #MaxQNNum shouldn't to small, because two qn can fuse into a new qn so that nQN can increase MaxQNNum is related to max rank of tensors in the tensor net
    #issue: 数目太多慢，太少不够，需要改进
    #MaxQNNum = 30 
    #for fermion hubbard model, if nu=0.1, it requeres this number larger
    #MaxQNNum:ClassVar = 20  
    MaxQNNum:ClassVar[int] = 100
    
    #QnClass:Type[QnBase] = QnU1
    #QnClass:Type = QnU1
    #QnClass:type = QnU1
    QnClass:ClassVar = QnU1
    def __init__(self, n:int, qns:List, dims:List[int])->None:
        """
        """
        #QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
        
        self.nQN = n
        self._dims = dims
        self._totDim:int = -1 
        self.QNs = qns
        
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        n = len(dims)
        assert dims is not None 
        if qns is None: 
            qns = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, -8, 8, -9, 9, -10, 10]
            assert len(qns)>= len(dims)
            qns = qns[: len(dims)]
        qns1 = [cls.QnClass(i) for i in qns]
        dims = list(dims)
        return cls(n=n, qns=qns1, dims=dims)

    @classmethod
    def set_base(cls, dim=None):
        """
            for U1 symm. the space can be identified with two spin-1/2 as follows 
            1  <-->    |uu>  
            -1 <-->    |dd>
            0a <-->    |du>    #this is bad, 0a--|ud>, 0b--|du>  should be more conventianal
            0b <-->    |ud>
             OR ??
            1  <-->    <uu|  
            -1 <-->    <dd|
            0a <-->    <ud| 
            0b <-->    <du|
        """

        qn_identity = QnU1(0)
        #qn_identity = cls.QnClass.qn_id()
        dims= [2, 1, 1]
        qns1 = (0, 1, -1)
        qns = [QnU1(i) for i in qns1] 
        qsp_base = QspU1(n=3, qns=qns, dims=dims)
        if 1:
            qsp_null = QspU1.null()
        return qn_identity, qsp_base, qsp_null
    
    @classmethod
    def set_base_1site(cls, dim=None):
        qn_identity = cls.QnClass.qn_id()
        dims= (1, 1)
        qns1 = (1, -1)
        qns = [QnU1(i) for i in qns1] 
        qsp_base = QspU1(n=2, qns=qns, dims=dims)
        qsp_null = QspU1.null()
        
        return qn_identity, qsp_base, qsp_null
    
    @classmethod
    def max(cls, trunc_dim, nqn=None):
        nqn = nqn if nqn is not None else 3

        qsp_max = {}
        for i in range(1, 10):
            qsp_max[(3, i*4)] = QspU1.easy_init((0, 1, -1), [2*i, 1*i, 1*i])
            #qsp_max[(5, i*4  + 2)] = QspU1.easy_init((0, 1, -1, 2, -2), i*np.array((2, 1, 1, 1, 1)))
        
        #qsp_max[(3, 14)] = QspU1.easy_init((0, 1, -1), np.array((8, 3, 3)))
        
        qn5 = (0, 1, -1, 2, -2)
        qsp_max[(5, 6)] = QspU1.easy_init( qn5, [2, 1, 1, 1, 1])
        qsp_max[(5, 10)] = QspU1.easy_init(qn5, [4, 2, 2, 1, 1])
        qsp_max[(5, 13)] = QspU1.easy_init(qn5, [5, 3, 3, 1, 1])
        qsp_max[(5, 14)] = QspU1.easy_init(qn5, [6, 3, 3, 1, 1])
        qsp_max[(5, 16)] = QspU1.easy_init(qn5, [6, 4, 4, 1, 1])
        qsp_max[(5, 17)] = QspU1.easy_init(qn5, [7, 4, 4, 1, 1])
        qsp_max[(5, 18)] = QspU1.easy_init(qn5, [8, 4, 4, 1, 1])
        qsp_max[(5, 21)] = QspU1.easy_init(qn5, [9, 5, 5, 1, 1])
        qsp_max[(5, 24)] = QspU1.easy_init(qn5, [10, 5, 5, 2, 2])
        qsp_max[(5, 27)] = QspU1.easy_init(qn5, [11, 6, 6, 2, 2])
        qsp_max[(5, 30)] = QspU1.easy_init(qn5, [12, 7, 7, 2, 2])
        qsp_max[(5, 40)] = QspU1.easy_init(qn5, [16, 9, 9, 3, 3])
        qsp_max[(5, 50)] = QspU1.easy_init(qn5, [20, 11, 11, 4, 4])
        try:
            return qsp_max[(nqn, trunc_dim)]
        except KeyError as err:
            keys= list(qsp_max.keys())
            key3 = [x for x in keys if x[0]==3]; key3.sort()
            key5 = [x for x in keys if x[0]==5]; key5.sort()
            raise Exception("key (%d, %d) is not found. \n allowed keys are \n\t%s\n\t%s"%(trunc_dim, nqn, key3, key5))
    
    def shift_qn(self, qn_delta):
        """
           change all the qns by qn_delta, not changing dim; an inplace operation
        """
        n= self.nQN
        if not isinstance(qn_delta, int):
            qn_delta = qn_delta._val
        for i in range(n):
            self.QNs[i]._val += qn_delta 
    
    def __mul__(self, other:'QspU1')->'QspU1': 
        #return self.add(other)
        return self.tensor_prod(other)
   
        
   
class QspU1_half(QuantSpaceBase):
    """
        half integer or odd U1 qn
    """
    pass

class QspSU2(QuantSpaceBase):
    MaxQNNum = 20 
    QnClass = QnSU2
    IS_Abelian = False
    def __init__(self, n, qns, dims, ) -> None:
        """
        """
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        n = len(dims)
        assert dims is not None 
        if qns is None: 
            qns = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, -8, 8, -9, 9, -10, 10]
            assert len(qns)>= len(dims)
            qns = qns[: len(dims)]
        qns1 = [cls.QnClass(i) for i in qns]
        dims = list(dims)
        return cls(n=n, qns=qns1, dims=dims)

    def tensor_prod_test(self, other):
        """ 
            refs:
                singh and vidal 2012 
            
            todo: change the name to prod in future 
            note:
                after tensor prod, the value of the qn in the qsp is not ordered. 
        """
        if self.nQN==0:
            return other.copy()    #need copy here?
        if other.nQN==0:
            return self.copy()     #need copy here?
    
        res = self.__class__(n=0, qns=[], dims=[])
        res._totDim = 0

        for i in range(self.nQN):
            for j in range(other.nQN):
                d0, d1 = self._dims[i], other._dims[j]
                if 0:
                    j0, m0 = self.QNs[i].val 
                    j1, m1 = other.QNs[j].val 
                    
                    m = m0 + m1
                    j = j0 + j1 
                    
                    #these 2 lines are comment temprarily when debug mypy 
                    #jj, val = wigner.drc3jj(j0, j1, m0, m1)
                    #arg = np.where(val)
                    
                    #print_vars(vars(),  ['j0, j1', 'm0, m1', 'd0, d1'])
                    #print_vars(vars(),  ['val', 'arg'])
                
                #res.add_to_quant_space(qn, d)
        return res




def symmetry_to_Qn(symmetry):
    temp = {"Travial":QnTravial, "Z2":QnZ2, "Z3":QnZ3, "U1":QnU1}
    return temp[symmetry]

def symmetry_to_Qsp(symmetry):
    dic = {"Travial":QspTravial, "Z2":QspZ2, "Z3":QspZ3, "U1":QspU1, 'SU2':QspSU2}
    return  dic[symmetry]

symmetry_to_QspClass= symmetry_to_Qsp

def make_qsp(symmetry, qns=None, dims=None):
    """
        inefficient, not for production use 
    """
    if qns is not None:
        if hasattr(qns[0], '_val'):
            qns= [q._val for q in qns]
    cls= symmetry_to_Qsp(symmetry)
    return cls.easy_init(qns, dims)

qsp_any  = any_qsp = make_qsp  #def qsp_any   def any_qsp

def qn_factory(symmetry, val):
    if symmetry == "Travial":
        return QnTravial()
    else:
        return symmetry_to_Qn(symmetry)(val)

qn_any = qn_factory

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    
    
    def test_shift_qn(self): 
        q = QspU1.easy_init(None, [2, 1, 1])
        print(q)
        q.shift_qn(2)
        q1 = QspU1.easy_init([2, 3, 1], [2, 1, 1])
        print(q1)
        self.assertTrue(q==q1)
      
    
    def test_old_all(self): 
        if 1: 
            q1 = QspU1.max(16, 3)
            q2 = QspU1.max(17, 5)
            
            assert (not q1<q2)  or (not q1==q2) or (not q1<q2)
            print(q2)
            print(q2 >= q1) 
            
        if 1: 
            q1=QspU1.easy_init( qns=(0, 1, -1), dims=(2, 2, 2) )
            q2=QspU1.easy_init( qns=(1, 2, 0), dims=(2, 2, 2) )
            print(q1,  q2)
            print(q1.tensor_prod(q2))
    
    def test_compare(self): 
        q1 = QspU1.max(8)
        q2 = QspU1.max(14, 5)
        
        #assert (not q1<q2)  or (not q1==q2) or (not q1<q2)
        print(q2)
        print(q1) 
        #print q2  ==   q1 
        self.assertTrue(q1 < q2)
        self.assertTrue(q2 >=  q1)
  
    def test_power(self):
        pass 
        q = QspZ2.easy_init( [1, -1], [2, 2])
        print(q)
        print(q**2)
        print(q*q) 
        self.assertTrue(q**3 == q*q*q) 

        q = QspU1.easy_init( [0, 1, -1], [2, 1, 1])
        print(q)
        print(q**2)
        print(q*q) 
        self.assertTrue(q**3 == q*q*q) 
        print(q**0) 

    def test_tensor_prod(self): 
        if 0:
            a = QspU1.easy_init([1, -1], [2, 4])
            b = QspU1.easy_init([1, -1], [4, 2])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            #print_vars(vars(), ['a', 'b', 'c'])
            #print_vars(vars(), ['a*b', 'b*a', 'a*c', 'c*a'])
            #print_vars(vars(), ['(a*b).QNs', '(b*a).QNs', '(a*c).QNs', '(c*a).QNs'])
        if 0:
            a = QspU1.easy_init([1, -1], [2, 4])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            #print_vars(vars(), ['a', 'b', 'c'])
            #print_vars(vars(), ['(a*b).QNs', '(b*a).QNs', '(a*c).QNs', '(c*a).QNs'])

        if 0:
            a = QspU1.easy_init([0, 1, 2], [1, 1, 1])
            #print_vars(vars(), ['a', 'a**2', 'a**3'])
            print(a.QNs , type(a.QNs), a.QNs.index(QnU1(2)))
        if 0:
            a = QspU1.easy_init([1, -1], [2, 4])
            b = QspU1.easy_init([1, -1], [4, 2])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            #print_vars(vars(), ['a*b*c', 'a*(b*c)'])
            
        if 1:
            a = qsp_any('Travial', qns=[1], dims=[2])
            b = qsp_any('Travial', qns=[1], dims=[2])
            c = a.tensor_prod(b)
            c1 = qsp_any('Travial', qns=[1], dims=[4])
            self.assertTrue(c==c1)
            
            
    def xtest_expand_totdim(self): 
        if 1: 
            q=QspZ2.easy_init([1, -1], [2, 3])
            qe = q.expand_totdim(11)
            #print_vars(vars(),  ['qe'])
            self.assertTrue(qe==QspZ2.easy_init([1, -1], [5, 6]))

        q=QspU1.easy_init([0, 1, -1], [4, 2, 3])
        qe = q.expand_totdim(26)
        #print_vars(vars(),  ['qe'])
        self.assertTrue(qe==QspU1.easy_init([0, 1, -1], [12, 6, 8]))

    def test_temp(self): 
        
        #a = A(1)
        #b = A(2)
        #c = a + b
        #print(c.i)
        
        a = QnU1(1)
        print(a)
        print(QspU1.MaxQNNum)
        pass
        #if 0:
        #    j0, j1, m0, m1 = 10, 10, -2, 2
        #    jj, val = wigner.drc3jj(j0, j1, m0, m1)
        #    arg = np.where(val)
        #    print_vars(vars(),  ['j0, j1', 'm0, m1', 'd0, d1'])
        #    print_vars(vars(),  ['val', 'arg'])
        #    
        #if 1:
        #    q0 = QnSU2((1, -1))
        #    q1 = QnSU2((1, 1))
        #    q2 = QnSU2((2, 0))
        #    print_vars(vars(),  ['q0.cg_coeff(q1, q2)'])
        #    print(clebsch_gordan(2, 2, 0 , -2, 2, 0))
        #    
        #if 0:
        #    qn = QnSU2((1, 1))
        #    qn1 = qn.conj()
        #    print_vars(vars(),  ['qn', 'qn1'])
        #    
        #    q = QspSU2.easy_init(qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[4, 4, 4, 4])
        #    print_vars(vars(),  ['q'])
        #    q = qsp_any('SU2', qns=[(0, 0), (2, -2), (2, 0), (2, 2)], dims=[4, 4, 4, 4])
        #    
        #    qq = q.tensor_prod(q)
        #    print_vars(vars(),  ['qq.totDim', 'q.totDim'])
        #    print_vars(vars(),  ['qq'])

        #if 1:
        #    q = QspSU2.easy_init(qns=[(1, 1), (1, -1)], dims=[1, 1])
        #    q = QspSU2.easy_init(qns=[(1, 1), (1, -1)], dims=[1, 1])
        #    qq = q.tensor_prod(q)
        #    print_vars(vars(),  ['qq.totDim', 'q.totDim'])
        #    print_vars(vars(),  ['qq'])


if __name__ == "__main__":
        
    if 0:
        #TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        del TestIt.test_temp
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_shift_qn', 
           #'test_old_all', 
           #'test_compare', 
           #'test_power', 
           #'test_tensor_prod', 
           #'test_expand_totdim', 
           'test_temp', 
            ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)
       

