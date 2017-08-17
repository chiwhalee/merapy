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
import unittest 
import numpy as np
import warnings

from merapy.utilities import print_vars
from merapy.decorators import decorate_methods, tensor_player 
from merapy.utilities import print_vars 

#__all__=["QuantumNum", "QuantSpace", "QN_idendity", "QSp_null", "QSp_base", "init_System_QSp", "QspU1", "QspZ2", "reset_System_QSp"]
__all__=["QuantumNum", "QuantSpace", "QuantSpaceBase",  "init_System_QSp", 
    "QnU1", "QnZ2", "QnZ3", "QnTravial", 
    "QspU1", "QspZ2", "QspZ3", "QspTravial", 
    "reset_System_QSp", "symmetry_to_Qn", "symmetry_to_Qsp", 'symmetry_to_QspClass', "qn_factory"]



#self.MaxQNNum= 16


#group_names=["Travial", "Z2", "Z2Z2", "U1", "U1Z2", "U1U1"]
GROUP_NAMES = ["Travial", "Z2", "Z3", "U1"]


class QnBase():
    """
    abstract base class, not meant for direct use
    """
    
    def __getstate__1(self):
        """ not workable"""
        #print "pickled"
        return self._val
    def __setstate__1(self, state):
        self._val = state["_val"]
        
    @property
    def val(self):
        return self._val

    @val.setter #it seems that setter effects only when QnBase inherite from object
    def val(self, new_val):
        #print 
        raise NotImplemented("this don't work, because of olb-style class")
        self._val = new_val
    
    def __repr__(self):
        """
        newly added
        """
        #res= repr(self.val)
        #res="(" + "%5d"*self.num_of_symm +")"  % tuple(self.val)
        if isinstance(self.val, np.ndarray):
            #res=("(" + "%5d"*QuantumNum.NUM_OF_SYMM+")")  % tuple(self.val)
            res=("(" + "%5d"*self.val.shape[0]+")")  % tuple(self.val)
            #print "in __repr______"
            #res= str(self.val)
        elif isinstance(self.val, int):
            res= ("(" + "%5d"")")  % self.val
            #res=("(" + "%5d"*self.val.shape[0]+")")  % tuple(self.val)
            #res=("(" + "%5d"*self.NUM_OF_SYMM+")")  % tuple(self.val)
        else:
            raise ValueError(self.val)
        return res
    
    def set_val_back(self,val):
        """
        newly added, assign val to self.val
        val: an instance of QuantumNum  or   an integer or an iterable
        """
        if val.__class__ == self.__class__:
            self.val[:]=val.val[:]
            return 
        if type(val)==int:
            self.val=np.array([val])  #(val,)
            return
        print "vvv", type(val), self.__class__
        print "self", type(self.val), "other", type(val)
        self.val[:]=val[:]
    
    def __eq__(self, other):
        #return np.all(self.val == other.val)
        return self._val == other._val 
    
    def __ne__(self, other):
        #return np.all(self.val != other.val)
        return self._val != other._val 

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
            res= res + i
        return res
    
    def qsp_class(self): 
        return symmetry_to_Qsp(self.SYMMETRY)
    
    def conj(self): 
        res = self.copy()
        res.reverse()
        return res 
    
class QnTravial(QnBase):
    SYMMETRY = "Travial"
    #this may be an issue NUM_OF_SYMM should be 0?
    NUM_OF_SYMM = 1
    QNS = (1, )
    QnId = 1
    def __init__(self, value=None):
        """
            the param value is just for consistency with other QnClass 
        """
        #self._val = np.array([1], np.int)
        self._val = 1

    def set_val(self, val):
        """
        _val always is 1
        """
        #print "ssss", type(val)
        #self._val = val
        self._val = 1

    def copy(self):
        return QnTravial()
    def reverse(self):
        pass # do nothing
    def __add__(self, other):
        return QnTravial()
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
        #print "qqqq", type(self), type(self._val), type(other._val)
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
    SYMMETRY = "Z2"   #precisely, this is used for spin-parity symm
    NUM_OF_SYMM = 1
    QNS = (1, -1)   #all the elements in the algebra, 
    QnId = 1   # identity in the algebra
    
    def __init__(self, value):
        """
            把val统一写成np.ndarray类型, 这样有利有弊: 利在可用array的各种运算，比如*法；弊在赋值麻烦些
            num_of_symm: newly added atribute by lzh, 
        """
        
        #self._val = np.array(value, dtype=np.int, ndmin=1)
        #self._val = np.array([value], "int")
        self._val = value
    
    def copy(self):
        #if use the following line, program will grow slower and slower, I don't know why yet
        return QnZ2(self._val)
        #return QnZ2(self._val[0])
        
    def reverse(self):
        """
        reverse refers to conjugate repr. of a group
        """
        pass

    def __add__(self, other):
        """
            see QN_Add in f90
            q: 所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnZ2(self._val*other._val)
    
    def __reduce__(self,other):
        return QnZ2(self.val*other.val)

class QnZ3(QnBase):
    SYMMETRY = "Z3"
    NUM_OF_SYMM = 1
    #QNS = (1, -1)
    QNS = (0, 1, 2)
    QnId = 0
    def __init__(self, value):
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

    def __add__(self, other):
        """
            see QN_Add in f90
            q: 所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
            故变成加法
        """
        return QnZ3((self._val + other._val)%3)
    
    def __reduce__(self, other):
        return QnZ3(self.val*other.val)

class QnU1(QnBase):
    SYMMETRY = "U1"
    NUM_OF_SYMM = 1
    QNS = tuple(range(-10, 11))
    QnId = 0

    def __init__(self, value):
        """
        把val统一写成np.ndarray类型, 这样有利有弊: 利在可用array的各种运算，比如*法；弊在赋值麻烦些
        num_of_symm: newly added atribute by lzh, 
        """
        
        #self._val = np.array(value, dtype=np.int, ndmin=1)
        self._val = value
    
    def __new__1(cls, value):
        self = super(QnU1, cls).__new__(cls, value)
        self._val = value
        return value
    
    def copy(self):
        #print self.val, self._val[0]
        #return  QnU1(self.val)
        return  QnU1(self._val)
        
    def reverse(self):
        """
        reverse refers to conjugate repr. of a group
        see QN_Reverse in f90
        """
        #self.val = -self.val
        self._val = -self._val


    def __add__(self, other):
        """
        q: 所谓 add 应该是该对称算子乘积的特征值, 因为是Abelian群，
        故变成加法
        """

        #return QnU1(self.val+other.val)
        return QnU1(self._val + other._val)
    
    def __reduce__(self,other):
        """
        see QN_Minus in f90
        status_1
        """
        qn3= QuantumNum()

        if self.SYMMETRY=="Z2" or "Z2Z2":
            qn3.val = self.val*other.val
            qn3.val = self.val*other.val
        if self.SYMMETRY=="U1":
            qn3.val = self.val-other.val
        if self.SYMMETRY=="U1Z2":
            qn3.val[0] = self.val[0]-other.val[0]
            qn3.val[1] = self.val[1]*other.val[1]
        return qn3
    

class QuantumNum(QnZ2, QnU1):
    SYMMETRY = "Z2"
    pass


meth_names= ["copy"]
#meth_names.pop(0)
@decorate_methods(tensor_player, meth_names)
class QuantSpaceBase(object):
    """
        abstract base class
        so MaxQNNum only defined in subclassed
    """
        
    def __init__(self, n, qns, dims):
        """
            QSp真正有用的只有 QNs, _dims 两个属性
            _dims: 每个量子数对应子空间的维数
        """
        self.nQN = n
        self._dims=np.empty(self.MaxQNNum, int) #MaxQNNum only defined in subclassed
        #self._dims=np.empty(n, int) 
        self._dims[:n] = dims[:n]
        #self.totDim = np.sum(self._dims[:n])  #for performance, this is not pre

        #self.QNs= np.ndarray(self.MaxQNNum, np.object)
        #self.QNs[:n]= qns[:n]
        self.QNs= qns
    
    @property
    def totDim(self):
        """
            lazy evaluation
        """
        if hasattr(self, "_totDim"):
            return self._totDim
        else:
            self._totDim = np.sum(self._dims[:self.nQN])  
            return self._totDim
    
    @property
    def symmetry(self): 
        return self.QnClass.SYMMETRY
    
    @property
    def Dims(self):
        return self._dims
    
    @Dims.setter
    def Dims(self, new_dims):
        self._dims[:len(new_dims)] = new_dims
        self._totDim = sum(self._dims)

    def __ge__(self, other): 
        if self.QnClass != other.QnClass: 
            return False
        if self.nQN < other.nQN: 
            return False
        return not self < other
        
    def __le__(self, other, info=0):
        if self.QnClass != other.QnClass: 
            if info>0: print '__le__ reason: ', 1
            return False
        if not self.nQN <= other.nQN:
            if info>0: print '__le__ reason: ', 2
            return False
        if not np.all(self._dims[:self.nQN]<=other._dims[:self.nQN]):
            if info>0: print '__le__ reason: ', 3, self._dims[:self.nQN], other._dims[:self.nQN]
            
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

    def __eq__(self, other):
        """
            improving: in futrue use key_property to implement this
        """
        key1 = ["nQN", "totDim"]
        for i in key1:
            res = self.__getattribute__(i)==other.__getattribute__(i)
            if not res:
                #print self.__getattribute__(i)
                #print other.__getattribute__(i)
                return res
        key2 = ["QNs", "_dims"]
        nQN = self.nQN
        for i in key2:
            res = self.__getattribute__(i)[:nQN]==other.__getattribute__(i)[:nQN]
            res= np.all(res)
            if not res:
                #print self.__getattribute__(i)
                #print other.__getattribute__(i)
                return res
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__old(self,exclude=["RefQN"]):
        """ 
       
        """
        keys= ["class", 'nQN', 'QNs', '_dims', 'RefQN']
        if exclude:
            for k in exclude:
                keys.remove(k)
        res= ""
        n = self.nQN
        for k in keys:
            if k == 'QNs': 
                temp = str(self.QNs[:n])
            elif k == '_dims': 
                temp = str(self._dims[:n])
            elif k == 'RefQN': 
                temp = str(self.RefQN[:, :n])
            elif k == "class":
                temp = self.__class__.__name__
            else:
                temp = str(self.__dict__[k])

            #res += k+":\n\t" + temp +"\n"
            res += k+":\t" + temp +"\n"

        return res
    
    def __repr__(self,exclude=["RefQN"]):
        """ 
       
        """
        keys= ["class", 'nQN', 'QNs', '_dims', 'RefQN']
        if exclude:
            for k in exclude:
                keys.remove(k)
        res= ""
        n = self.nQN
        temp = []
        for i in range(n): 
            val = self.QNs[i]._val
            dim = self._dims[i]
            s= '[%d]%d'%(val,dim)
            temp.append(s)
        res = '+'.join(temp)

        return res
    
    def __mul__(self, other): 
        return self.add(other)
    
    def __div__(self, other): 
        """
            q1=QspZ2.easy_init([1, -1], [2, 2])
            q2=QspZ2.easy_init([1, -1], [1, 3])
            print  q1*q2   == q1*q1 
            #they are equal 
        """
        msg = 'divide of qsp cant be defined,  as it is not unique. see doc string'
        raise NotImplemented(msg)
        
    
    def __pow__(self, n): 
        res= self.__class__.null()
        for i in range(n): 
            res= res*self
        return res 
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        n = len(qns)
        symm =  cls.QnClass.SYMMETRY 
        if symm== 'Travial':
            qns1 = None
        else:
            qns1 = [cls.QnClass(i) for i in qns]
        return cls(n=n, qns=qns1, dims=dims) 

    def reverse(self):
        """
        see QSp_Reverse in f90"""
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
    def copy(self, reverse=False):
        qns= [q.copy() for q in self.QNs[:self.nQN]]
        other=self.__class__(n=self.nQN, qns=qns, dims=self._dims)
        other._totDim = self.totDim
        if hasattr(self, "RefQN"):
            other.RefQN=self.RefQN.copy()
        if reverse: 
            other.reverse()
        return other
    
    def copy_from(self, other):
        self.nQN = other.nQN
        self._totDim=other.totDim
        if hasattr(self, "RefQN"):
            self.RefQN=other.RefQN.copy()
        nqn=other.nQN
        if nqn != 0:
            #print other.QNs
            self.QNs[:nqn]=[other.QNs[i].copy() for i in range(nqn)]
        self._dims=other._dims.copy()

    def update(self, qsp_max=None):
        """
        calculate totDim, RefQN
        see QSp_update in f90
        """
        D_Max0 = np.ndarray(self.MaxQNNum, int)
        D_Max0[:]=10000000   

        #temp=self.__class__.__call__()
        #temp = self.__class__.__call__(n=0, qns=[], dims=[])
        temp = self.__class__.empty()
        if qsp_max is not None:
            #temp.RefQN[:,:]=0
            #temp.nQN=0
            nqn = 0
            for n in range(self.nQN):
                i= qsp_max.has_quant_num(self.QNs[n])
                if i>=0:
                    temp.QNs.append(self.QNs[n])
                    #temp.QNs[nqn]=self.QNs[n]
                    #temp.nQN=temp.nQN+1
                    temp._dims[nqn]=min(self._dims[n], qsp_max._dims[i] )
                    nqn += 1 
                    #temp.RefQN[0,i]=min(self.RefQN[0,n], qsp_max.RefQN[0,i])
                    #temp.RefQN[1,i]=min(self.RefQN[1,n], qsp_max.RefQN[1,i])
            temp.nQN = nqn
            self.copy_from(temp)
            
            #self.__dict__.update(temp.__dict__)
        if 0:
            n=0
            #self.Addr[n]=0
            totDim=self._dims[0]
            
            for n in range(1, self.nQN):
                #self.Addr[n]=self.Addr[n-1]+self._dims[n-1]
                if self._dims[n]>D_Max0[n]:
                    self._dims[n]=D_Max0[n]
                totDim += self._dims[n]
            n=self.nQN+1
            
            #attention_this_may_be_wrong
            #self.Addr[n]=self.Addr[n-1] + self._dims[n-1]
            self._totDim =totDim

    def copy_many(self, n, reverse=None):
        res=[self.copy() for i in xrange(n)]
        if reverse is not None:
            for i in reverse:
                res[i].reverse()
        return res


    def link(self):
        return self

    def is_equal(self,other):
        """
        deprecated, use __eq__ instead
        status_1
        see QSp_IsEq in f90
        """
        #keys=self.__dict__.keys()
        keys= ["nQN", "totDim", "QNs","_dims"]
        
        for k in keys:
            a=self.__dict__[k]; 
            b=other.__dict__[k]
            res= a == b
            #if type(a) !=np.ndarray:
            if not isinstance(a, np.ndarray):
                if res != True:
                    print k+" not equal " , a, b
                    return res
            else:
                res=res.all()
                if res != True:
                    print k+" not equal, their difference is "
                    #ind = np.where(res1==False)
                    #print ind
                    print a-b
                    return res
        return True


    def key_property(self):
        qns = tuple([q._val for q in self.QNs]) 
        dims= tuple(self._dims[:self.nQN].tolist())
        res = {'qns':qns, 'dims':dims, 'nqn':len(qns)}
        return res

    def square(self):
        """
        see CombineQSp2 in f90
        reture another instace of QuantSpace that is the squre of self
        """
        RefQN=np.ndarray(2,self.MaxQNNum)
        res= QuantSpace()
        for i in range(self.nQN):
            for j in range(self.nQN):
                qn=self.QNs[i]+self.QNs[j]
                n,m= self._dims[i], self._dims[j]
                d=n*m
                if i==j:
                    res.add_to_quant_space(qn,d)
                else:
                    res.add_to_quant_space(qn,d+d)
                if pidx>= self.MaxQNNum:
                    continue 
                if i==j:
                    RefQN[0,pidx] = RefQN[0,pidx]+n*(n+1)/2
                    RefQN[1,pidx] = RefQN[1,pidx]+n*(n-1)/2
                else:
                    RefQN[0,pidx] = RefQN[0,pidx]+d
                    RefQN[1,pidx] = RefQN[1,pidx]+d
        res.Addr=0
        res.RefQN=RefQN
        return res
    def square3(self):
        """
        see CombineQSp3 in f90
        reture another instace of QuantSpace that is the squre of self
        """
        RefQN=np.ndarray(2,self.MaxQNNum)
        res= QuantSpace()
        for i in range(self.nQN):
            for j in range(self.nQN):
                qn=self.QNs[i]+self.QNs[j]
                n,m= self._dims[i], self._dims[j]
                d=n*m
                pidx=res.add_to_quant_space(qn,d)
                if pidx>= self.MaxQNNum:
                    continue 
                if i==j:
                    RefQN[0,pidx] = RefQN[0,pidx]+n*(n+1)/2
                    RefQN[1,pidx] = RefQN[1,pidx]+n*(n-1)/2
                else:
                    RefQN[0,pidx] = RefQN[0,pidx]+d
                    RefQN[1,pidx] = RefQN[1,pidx]+d
        res.Addr=0
        res.RefQN=RefQN
        return res

    def has_quant_num(self, qn):
        """
            see InQuantSpace in f90
            qn: instance of QuantumNum  or a tuple
        
        """
        if self.nQN ==0:  #empty quantum space
            res = -1    # -1 indicates self.nQN=0
            return res
        
        for i in xrange(self.nQN):
            if qn == self.QNs[i]:
                return i 
        res = -1
        return res
    
    def get_qn_id(self, qn):
        #try:
        #    return self.QNs.index(qn)
        #except ValueError:
        #    return -1 
        return self.QNs.index(qn) 

    def add_to_quant_space(self, qn, d):
        """
            this is in fact direct sum of vector spaces subject to symmetry
            note that direct sum is not a commutable operation
        """
        i = self.has_quant_num(qn)
        ind= abs(i)
        if ind > self.MaxQNNum: return None

        #when qn not in self.QNs
        if i<0:
            i=self.nQN
            self.QNs.append(qn)

            self.nQN+=1
            self._dims[i]=0

        # when qn in self.QNs 的情况
        self._dims[i] += d
        if not hasattr(self, "_totDim"):
            self._totDim = self.totDim
        self._totDim += d
        return ind
   
    def tensor_prod(self, other):
        """ 
            acturally tensorprod of self and other,  though named add
            q571, RefQN is not added!, RefQN[:,:] is arbitrary valued
            RefQN 写起来比较麻烦
            todo: change the name to prod in future 
        """
        if self.nQN==0:
            return other.copy()    #need copy here?
        if other.nQN==0:
            return self.copy()     #need copy here?
    
        res = self.__class__(n=0, qns=[], dims=[])
        res._totDim = 0

        for i in xrange(self.nQN):
            for j in xrange(other.nQN):
                qn =self.QNs[i] + other.QNs[j]
                d= self._dims[i]*other._dims[j]
                res.add_to_quant_space(qn, d)
        return res
    
    add = tensor_prod 
    
    @staticmethod
    def prod_many( qsp_list): 
        q = qsp_list[0]
        for a in qsp_list[1: ]: 
            q = q.tensor_prod(a)
        return q 
    
    @classmethod
    def null(cls):
        qn_identity = cls.QnClass.qn_id()
        qsp_null = cls(n=1, qns=[qn_identity], dims=[1])
        #qsp_null.update()
        return qsp_null
    
    @classmethod
    def empty(cls):
        return cls(n=0, qns=[], dims=[])
    
    def qn_to_ind(self, qn): 
        """
            reverse of self.QNs[ind]
        """
        #return self.QNs.index(qn)
        for i in range(len(self.QNs)): 
            if self.QNs[i]._val == qn._val: 
                return i
    
    def expand_totdim(self, totdim): 
        qns = self.QNs
        qns = [q.val for q in qns]
        dims = self.Dims[:self.nQN].copy() 
        if 0: 
            n = totdim/self.totDim
            residual = totdim%self.totDim 
            print_vars(vars(),  ['qns', 'dims', 'self.totDim', 'n'])
            dims *= n 
            #residual = totdim-()
            if residual>0: 
                #for i, j in enumerate(xrange(totdim-dim.size)): 
                for i in xrange(residual):
                    print_vars(vars(),  ['residual', 'i'])
                    dims[i] += 1
        
        ratial = dims/float(self.totDim)
        #print_vars(vars(),  ['ratial'])
        dims = totdim * ratial
        dims= dims.astype(int)
        diff = np.sum(dims)-totdim
        if diff >0 or diff >self.nQN: 
            raise
        else: 
            for i in xrange(-diff): 
                dims[i] += 1 
        
        res= self.__class__.easy_init(qns, dims)
        return res 
                
                

class QspTravial(QuantSpaceBase):
    MaxQNNum = 1
    QnClass= QnTravial
    def __init__(self, n, qns, dims):
        self.nQN = n
        self._dims = np.empty(self.MaxQNNum, np.int) #MaxQNNum only defined in subclassed
        self._dims[:n] = dims[:n]
        #self._totDim = self._dims[0]
        self.QNs = [QnTravial()]
        #self.QNs = np.array([QnTravial()])
        self.RefQN = np.empty((2, 1))
        #self.Addr=np.empty(self.MaxQNNum+1,"int")


    @classmethod
    def set_base(cls, dim=None):
        if dim is None:
            dim = 2
        qsp_base = QspTravial(n=1, qns=None, dims=[dim])
        
        #qn_identity = QnTravial()
        qn_identity = cls.QnClass.qn_id()
        
        qsp_null = QspTravial.null()
        return qn_identity, qsp_base, qsp_null 
    
    @classmethod
    def max(cls, trunc_dim, nqn=None):
        #if trunc_dim%2 != 0:
        #    raise ValueError("trunc_dim should be multipy of 2.  trunc_dim=%d"%trunc_dim)
        qsp_max = cls(n=1, qns=None, dims=[trunc_dim])
        if 1:
            qsp_max2 = cls(n=1, qns=None, dims=[trunc_dim])
            #qsp_max2.RefQN[0,0:2] = [2,1]
            #qsp_max2.RefQN[1,0:2] = [0,1]

        #return qsp_max, qsp_max2 
        return qsp_max

class QspZ2(QuantSpaceBase):
    MaxQNNum = 2
    QnClass = QnZ2
    def __init__(self, n, qns, dims,  RefQN=None):
        """

        """
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
        self.RefQN=np.ndarray((2, self.MaxQNNum), np.int)
    
    def __div__(self, other): 
        """
            solve this linear eqn of (d0, d1): 
                d0*b0 + d1*b1 = c0
                d0*b1 + d1*b0 = c1
            -->
                d0*b0*b1 + d1*b1^2 = c0*b1
                d0*b1*b0 + d1*b0^2 = c1*b0
            -->
                d1 = (c0*b1 - c1*b0)/(b1^2 - b0^2)
                d0 = (c1*b0 - c0*b1)/(b1^2 - b0^2)
                
                d0, d1 may be not unique 
                
        """
        raise NotImplemented
        
        b0 = other._dims[0]
        b1 = other._dims[1]
        
        c0 = self._dims[0]
        c1 = self._dims[1]
        res= self.__class__.easy_init([1, -1], [d0, d1])
        return res 

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
            qsp_max2.RefQN[0,0:2] = [2,1]
            qsp_max2.RefQN[1,0:2] = [0,1]
            qsp_max2.update()
        #return qsp_max, qsp_max2
        return qsp_max
    
class QspZ3(QuantSpaceBase):
    MaxQNNum = 3
    QnClass = QnZ3
    def __init__(self, n, qns, dims,  RefQN=None):
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
        self.RefQN=np.ndarray((2, self.MaxQNNum), np.int)

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

        if 0:
            dims= [3, 3]
            qns=[QnZ3(i) for i in qns1]
            qsp_max2 = QspZ3(n=3, qns=qns, dims=dims)
            #qsp_max2.RefQN[0,0:2] = [2,1]
            #qsp_max2.RefQN[1,0:2] = [0,1]
            qsp_max2.update()
        #return qsp_max, qsp_max2
        return qsp_max

class QspU1(QuantSpaceBase):
    #MaxQNNum shouldn't to small, because two qn can fuse into a new qn so that nQN can increase MaxQNNum is related to max rank of tensors in the tensor net
    #issue: 数目太多慢，太少不够，需要改进
    #MaxQNNum = 30 
    #MaxQNNum = 15
    MaxQNNum = 20 #for fermion hubbard model, if nu=0.1, it requeres this number larger
    QnClass = QnU1
    def __init__(self, n, qns, dims, RefQN=None):
        """
        """
        QuantSpaceBase.__init__(self, n=n, qns=qns, dims=dims)
    
    @classmethod
    def easy_init(cls, qns=None, dims=None):
        """ a slow but easy init """
        n = len(dims)
        assert dims is not None 
        if qns is None: 
            qns = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 
                    6, -6, 7, -7, -8, 8, -9, 9, -10, 10]
            assert len(qns)>= len(dims)
            qns = qns[: len(dims)]
        qns1 = [cls.QnClass(i) for i in qns]
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

        #qn_identity = QnU1(0)
        qn_identity = cls.QnClass.qn_id()
        dims= (2, 1, 1)
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
            qsp_max[(3, i*4)] = QspU1.easy_init((0, 1, -1), i*np.array((2, 1, 1)))
            #qsp_max[(5, i*4  + 2)] = QspU1.easy_init((0, 1, -1, 2, -2), i*np.array((2, 1, 1, 1, 1)))
        
        #qsp_max[(3, 14)] = QspU1.easy_init((0, 1, -1), np.array((8, 3, 3)))
        
        qn5 = (0, 1, -1, 2, -2)
        qsp_max[(5, 6)] = QspU1.easy_init( qn5, np.array((2, 1, 1, 1, 1)))
        qsp_max[(5, 10)] = QspU1.easy_init(qn5, np.array((4, 2, 2, 1, 1)))
        qsp_max[(5, 13)] = QspU1.easy_init(qn5, np.array((5, 3, 3, 1, 1)))
        qsp_max[(5, 14)] = QspU1.easy_init(qn5, np.array((6, 3, 3, 1, 1)))
        #qsp_max[(5, 16)] = QspU1.easy_init(qn5, np.array((8, 3, 3, 1, 1)))
        qsp_max[(5, 16)] = QspU1.easy_init(qn5, np.array((6, 4, 4, 1, 1)))
        qsp_max[(5, 17)] = QspU1.easy_init(qn5, np.array((7, 4, 4, 1, 1)))
        qsp_max[(5, 18)] = QspU1.easy_init(qn5, np.array((8, 4, 4, 1, 1)))
        #qsp_max[(5, 20)] = QspU1.easy_init(qn5, np.array((8, 4, 4, 2, 2)))
        qsp_max[(5, 21)] = QspU1.easy_init(qn5, np.array((9, 5, 5, 1, 1)))
        qsp_max[(5, 24)] = QspU1.easy_init(qn5, np.array((10, 5, 5, 2, 2)))
        qsp_max[(5, 27)] = QspU1.easy_init(qn5, np.array((11, 6, 6, 2, 2)))
        qsp_max[(5, 30)] = QspU1.easy_init(qn5, np.array((12, 7, 7, 2, 2)))
        qsp_max[(5, 40)] = QspU1.easy_init(qn5, np.array((16, 9, 9, 3, 3)))
        qsp_max[(5, 50)] = QspU1.easy_init(qn5, np.array((20, 11, 11, 4, 4)))
        try:
            return qsp_max[(nqn, trunc_dim)]
        except KeyError as err:
            keys= qsp_max.keys()
            key3 = filter(lambda x:x[0]==3, keys); key3.sort()
            key5 = filter(lambda x:x[0]==5, keys); key5.sort()
            raise Exception("key (%d, %d) is not found. \n allowed keys are \n\t%s\n\t%s"%(trunc_dim, nqn, key3, key5))
    
class QspU1_half(QuantSpaceBase):
    """
        half integer or odd U1 qn
    """
    pass

#class QuantSpace(QspZ2, QspU1):  pass

class QuantSpace(QuantSpaceBase):  pass


def symmetry_to_Qn(symmetry):
    temp = {"Travial":QnTravial, "Z2":QnZ2, "Z3":QnZ3, "U1":QnU1}
    return temp[symmetry]

def symmetry_to_Qsp(symmetry):
    
    return {"Travial":QspTravial, "Z2":QspZ2, "Z3":QspZ3, "U1":QspU1}[symmetry]

symmetry_to_QspClass= symmetry_to_Qsp

def make_qsp(symmetry, qns=None, dims=None):
    """
        inefficient, not for production use 
    """
    cls= symmetry_to_Qsp(symmetry)
    return cls.easy_init(qns, dims)

def qn_factory(symmetry, val):
    if symmetry == "Travial":
        return QnTravial()
    else:
        return symmetry_to_Qn(symmetry)(val)

def qsp_factory(symmetry, qns, dims):
    pass

#mapper = {"Z2":QspZ2, "Z3":QspZ3, "U1":QspU1, "Travial":QspTravial, "Trav":QspTravial}

def init_System_QSp(symmetry="Z2", dim=None, combine_2site=True):
    try:
        #cls = mapper[symmetry]
        cls= symmetry_to_Qsp(symmetry)
    except:
        raise KeyError(symmetry)
    
    QuantumNum.SYMMETRY = symmetry
    qn_identity, qsp_base, qsp_null = cls.set_base(dim=dim)
    
    #if symmetry == "Travial" and combine_2site:
    if symmetry in ["Travial", 'Z2'] and combine_2site:
        #qn_identity, qsp_base, qsp_null = QspTravial.set_base(dim=4)
        qn_identity, qsp_base, qsp_null = cls.set_base(dim=4)

    if symmetry == "U1" and not combine_2site:
        qn_identity, qsp_base, qsp_null = QspU1.set_base_1site(dim=dim)

    
    return qn_identity, qsp_base, qsp_null


def reset_System_QSp(symmetry):
    #global QN_idendity, QSp_base, QSp_null
    warnings.warn("QuantSpace has set base class %s"%symmetry)
    try:
        cls= mapper[symmetry]
    except:
        raise KeyError(symmetry)
    
    QuantSpace.__bases__ = (cls, )

    #QN_idendity, QSp_base , QSp_null= init_System_QSp(symmetry=symmetry)
    #this dont work, because after assignment,  they become local,  although they are assumed global
    return init_System_QSp(symmetry=symmetry)


if 0:
    def swap(qn1,qn2,isZ2):
        """
        q: 不理解意思
        status_1
        """
        if isZ2:
            res= qn1==-1 and qn2 == -1
        else:
            res= (abs(qn1)//2)==1 and (abs(qn2//2)==1)
        return res

    def swap1(QN1, QN2):
        if QuantumNum.SYMMETRY =="Z2":
            res= swap(QN1.val, QN2.val, True)
        elif QuantumNum.SYMMETRY =="Z2Z2":
            res= swap(QN1.val[1],QN2.val[0], True)
        elif QuantumNum.SYMMETRY=="U1":
            res= swap(QN1.val, QN2.val, False)
        elif QuantumNum.SYMMETRY=="U1U1":
            res= swap(QN1.val[1],qn2.val[0], False)
        return res

    def swap2(QN1, QN2):
        if QuantumNum.SYMMETRY =="Z2":
            res= swap(QN1.val, QN2.val, True)
        elif QuantumNum.SYMMETRY =="Z2Z2":
            res= swap(QN1.val[0],QN2.val[1], True)
        elif QuantumNum.SYMMETRY=="U1":
            res= swap(QN1.val, QN2.val, False)
        elif QuantumNum.SYMMETRY=="U1U1":
            res= swap(QN1.val[0], QN2.val[1], False)
        return res
    def qspace_get_dim(n,QSp_list,nd,_dims):
        """
        status_1_uncheck
        see QSpace_GetDim in f90
        q: still don't understand this function
        """
        pTot=1
        for i in range(n):
            pTot*=QSp_list[i].nQN
        nd=pTot

        _dims=[]

        iQN =[1L for i in range(n)]

        for p in range(pTot):
            d=1
            for i in range(n):
                d *= QSp_list[i]._dims([iQN[i]])
                _dims.append(d)
                inc =1
                
                i=0
                while inc==1 and i <=n:
                    iQN[i]+= 1
                    if iQN[i] <= QSp_list[i].nQN:
                        inc =0
                    else:
                        iQN[i]=1
                        i+=1

        #attention_this_may_be_wrong
        return nd, _dims


# ====================================================================================== 

class test_qn(object):
    def __init__(self, symmetry):
        self.symmetry = symmetry
        self.qn_identity, self.qsp_base,  self.qsp_null = init_System_QSp(symmetry)

        pass
    def add_to_quant_space(self):
        qb= self.qsp_base.copy()
        qn = self.qn_identity.copy()
        qb.add_to_quant_space(qn, 3)
        print qb
    
    def reverse(self):
        qn = self.qn_identity.copy()
        QSp_base = self.qsp_base
        qn.val = -1
        qn.reverse()
        print qn

        qs= QSp_base.copy()
        qs.reverse()
        print qs
    
    def has_quant_num(self):
        QSp_base = self.qsp_base.copy()
        print "\ntest of has_quant_num: ", QSp_base.has_quant_num(QSp_base.QNs[0]) #QSp_null.QNs[0])
        
    def add(self):
        print "\ntest of self.add(): "
        
        QN_idendity, QSp_base , QSp_null= init_System_QSp(symmetry=self.symmetry) 
        #print QSp_base
        temp= QSp_base.add(QSp_base) #.add(QSp_base)
        #print QSp_base
        print temp
    
    def update(self):
        print "uuuuu"*10
        #QN_idendity, QSp_base , QSp_null= init_System_QSp(symmetry=self.symmetry, combine_2site=False) 
        
        if self.symmetry == "U1" :
            QSp_base = QspU1.easy_init((1, -1), (1, 1))
            print QSp_base
            temp= QSp_base.add(QSp_base).add(QSp_base)
            print temp
            ##print QSp_base.__class__.max(4)
            qsp_max = QspU1.easy_init((1, -1, 3, -3), (2, 2, 1, 1))
            print "max", qsp_max
            temp.update(qsp_max)
            print temp
    
    def update_2(self):
        print "uuuuu"*10
        #QN_idendity, QSp_base , QSp_null= init_System_QSp(symmetry=self.symmetry, combine_2site=False) 
        
        if self.symmetry == "U1" :
            QSp_base = QspU1.easy_init((1, -1, 3, -3), (2, 2, 1, 1)) 
            print QSp_base
            temp= QSp_base.add(QSp_base).add(QSp_base)
            print temp
            ##print QSp_base.__class__.max(4)
            qsp_max = QspU1.easy_init((1, -1, 3, -3), (2, 2, 1, 1))
            print "max", qsp_max
            temp.update(qsp_max)
            print temp

    @staticmethod
    def add_u1():
        QN_idendity, QSp_base=init_System_QSp("U1")
        qsp_max=QspU1()
        qsp_max.nQN = 3
        temp=[0,1,-1,2,-2]
        for i in range(5):
            qsp_max.QNs[i].val[0] = temp[i] 
        qsp_max._dims[0:5] = [2,1,1,1,1]

        qsp_max.update()
        qsp_max2=qsp_max.copy()  # not used
        
        a=QSp_base.add(QSp_base).add(QSp_base)
        print QSp_base
        print a
        a.update(qsp_max=qsp_max)
        print a
    
    def property(self):
        print "property not work at present"
        #qsp_base = self.qsp_base.copy()
        qn = QnU1(1)
        qn._val = -2
        assert qn.val  == -2
        if 0:
            qn = QnU1(0)
            qn.val = -2
            print qn._val
            assert qn._val == -2, "qn._val = %d"%qn._val
            #assert qn._val == -2 

    def pickle(self):
        #import pickle as cPickle
        import cPickle
        QN_idendity = QnU1(0)
        out = open("/tmp/test", "wb")
        cPickle.dump(QN_idendity, out)
        out.close()
        #exit()

        qq = QspU1()

        out = open("/tmp/test", "wb")
        cPickle.dump(qq, out)
        out.close()
        #exit()

    def le(self):
        a = QspU1.max_1((3, 4))
        b = QspU1.max_1((5, 14))
        print a <= b
        exit()

    def eq(self):
        pass
        a = QspU1.max(12, 3)
        b = QspU1.max(12, 3)
        c = a.copy()
        print a == b
        
class _performence(object):
    @classmethod
    def copy_qn(cls):
        import time
        t1=time.clock()
        for i in range(20000):
            a=QN_idendity.copy()
        t2=time.time()
        for i in range(20000):
            a=QN_idendity.copy1()
        t3=time.clock()
        print t2-t1
        print t3-t2
    @classmethod      
    def copy_qsp(cls):
        import time
        t1=time.clock()
        for i in range(100):
            QSp_base.copy()
        t2=time.clock()
        for i in range(100):
            QSp_base.copy_old()
        t3=time.clock()
        print "copy",t2-t1
        print "copy1",t3-t2

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_temp(self): 
        q = QspU1.easy_init(None, [2, 1, 1])
        print(q)
      
    
    def test_old_all(self): 
        for symm in GROUP_NAMES:
        #for symm in ["U1"]:
        #for symm in []:
            print symm*10
            tq=test_qn(symmetry=symm)
            print tq.qsp_base
            print tq.qsp_base.key_property()
            pf=_performence
            #reset_System_QSp(symmetry="Trav")
            tq.has_quant_num()
            tq.add()
            print "1"*100
            tq.update()
            tq.reverse()
            #tq.add_to_quant_space()
            #tq.property()
            #tq.pickle()
            tq.eq()
            #pf.copy_qn()
            #pf.copy_qsp()
        if 1: 
            q1 = QspU1.max(16, 3)
            q2 = QspU1.max(17, 5)
            
            assert (not q1<q2)  or (not q1==q2) or (not q1<q2)
            print q2
            print q2 >= q1 
            
        if 1: 
            q1=QspU1.easy_init( qns=(0, 1, -1), dims=(2, 2, 2) )
            q2=QspU1.easy_init( qns=(1, 2, 0), dims=(2, 2, 2) )
            print q1,  q2
            print q1.tensor_prod(q2)
    
    def test_compare(self): 
        q1 = QspU1.max(8)
        q2 = QspU1.max(14, 5)
        
        #assert (not q1<q2)  or (not q1==q2) or (not q1<q2)
        print q2
        print q1 
        #print q2  ==   q1 
        self.assertTrue(q1 < q2)
        self.assertTrue(q2 >=  q1)
  
    def test_power(self):
        pass 
        q = QspZ2.easy_init( [1, -1], [2, 2])
        print  q
        print q**2
        print q*q 
        self.assertTrue(q**3 == q*q*q) 

        q = QspU1.easy_init( [0, 1, -1], [2, 1, 1])
        print  q
        print q**2
        print q*q 
        self.assertTrue(q**3 == q*q*q) 
        print  q**0 

    def test_tensor_prod(self): 
        if 0:
            a = QspU1.easy_init([1, -1], [2, 4])
            b = QspU1.easy_init([1, -1], [4, 2])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            print_vars(vars(), ['a', 'b', 'c'])
            print_vars(vars(), ['a*b', 'b*a', 'a*c', 'c*a'])
            print_vars(vars(), ['(a*b).QNs', '(b*a).QNs', '(a*c).QNs', '(c*a).QNs'])
        if 0:
            a = QspU1.easy_init([1, -1], [2, 4])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            print_vars(vars(), ['a', 'b', 'c'])
            print_vars(vars(), ['(a*b).QNs', '(b*a).QNs', '(a*c).QNs', '(c*a).QNs'])

        if 0:
            a = QspU1.easy_init([0, 1, 2], [1, 1, 1])
            print_vars(vars(), ['a', 'a**2', 'a**3'])
            print a.QNs , type(a.QNs), a.QNs.index(QnU1(2))
        if 1:
            a = QspU1.easy_init([1, -1], [2, 4])
            b = QspU1.easy_init([1, -1], [4, 2])
            c = QspU1.easy_init([1, 0, -1], [4, 3, 2])
            print_vars(vars(), ['a*b*c', 'a*(b*c)'])
            
    def test_expand_totdim(self): 
        if 1: 
            q=QspZ2.easy_init([1, -1], [2, 3])
            qe = q.expand_totdim(11)
            print_vars(vars(),  ['qe'])
            self.assertTrue(qe==QspZ2.easy_init([1, -1], [5, 6]))

        q=QspU1.easy_init([0, 1, -1], [4, 2, 3])
        qe = q.expand_totdim(26)
        print_vars(vars(),  ['qe'])
        self.assertTrue(qe==QspU1.easy_init([0, 1, -1], [12, 6, 8]))
        
    
if __name__ == "__main__":
        
    if 0:
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
           'test_temp', 
           #'test_old_all', 
           #'test_compare', 
           #'test_power', 
           #'test_tensor_prod', 
           #'test_expand_totdim', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)
       

