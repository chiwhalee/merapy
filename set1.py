#coding=utf8
"""
    see Set.f90

    名为set，实际上是ordered set！！

    直接用 builtin set，功能足够and even more

    let A be of type set:
    对应如下：
        attributes:
            A%N--> len(A)
            A%data--> no corresponding, note that set in python is a type, not a class?
        methods:
            Set2Array -->  np.array(list(set))
            Array2Set -->  set(A)
            CommInSet -->  s1.intersection(b)
            remove  -->  s1.remove(a)
            MinusSet -->  s1-s2
            JointSet --> s1.union(s2)

"""
import numpy as np
from copy import deepcopy

class Set_back(object):
    MaxLength=32
    def __init__(self,N=0,data=np.ndarray(MaxLength,"int")):
        self.N=N
        self.data=data
    def array_to_set(self,n,A):
        """
        see Array2Set in f90
        """
        self.N = 0
        for i in range(n):
            self.insert_new(A[i],self)

    def insert_new(self):
        
        Insert = InSet(m, S)
        if Insert < 0:  
            self.N = self.N+1
            self.data[self.N]=m
            Insert=self.N
    def to_array(self):
        """
        status_1
        see Set2Array in f90
        Returns: a np.ndarray
        """
        arr= np.ndarray(self.data[:self.N],"int")
        return arr

class Set(object):

    MaxLength=32

    def __init__():
        self.n=0
        self.data=[0]*self.MaxLength

    def copy(self):
        """ 
        see AssignSet in f90
        """
        return deepcopy(self)

    def has_elemnet(self,m):
        """
        see InSet in f90
        """
        if m in self.data:
            return self.data.index(m)
        else:
            return -1
    def inset(self,m):
        """
        see Insert2Set in f90
        """
        if m not in self.data:
            self.data.append(m)
            self.n=self.n+1

    def remove(self,m):
        """
        see Remove in f90
        """
        if m in self.data:
            self.data.remove(m)
            self.n=self.n-1

    def intersect(self,other):
        """
        see CommInSet in f90
        """
        res= Set()



class OrderSet():

    @staticmethod
    def intersect(a, b):
        return [i for i in a if i in b]

    @staticmethod
    def symm_diff(a, b):
        """
        symmetric difference of lists a and b, that is,  a\cup b  - (a\cap b)
        """
        return [i for i in a if i not in b] + [i for i in b if i not in a]
    
    @staticmethod
    def diff(a, b):
        """
        a\b
        """
        return [i for i in a if i not in b]
    




if __name__ =="__main__":
    ss=Set()
    print ss.MaxLength
    print type(ss.to_array()), ss.to_array()



