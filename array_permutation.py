#coding=utf8
"""
see ArrayPermutation.f90


issue:
performance:
        
        order_[:rank]=order[:rank]
        Dims_[:rank]=Dims[:rank],  could also be removed
        out_buff if not good choice

"""
import numpy as np
import platform
import warnings
import unittest
import socket
from merapy.utilities import print_vars

arch=  platform.architecture()[0]   #32 or 64 
os1 = platform.system()
hostname = socket.gethostname()

try:            
    if hostname == 'VirtualBox-Lab':
        from merapy.lib.array_permutation_64_ifort_virtual import (
               array_permutation_fort_parallel , 
               #I have several versions of permute_player, the parallel dynamic is supposed to be fattest
               #but it is not thoroughly tested against too mainy cases
               permute_player_fort_parallel_dynamic as permute_player_fort  
                )
    else:
        if os1 == 'Linux':
            from merapy.lib.array_permutation_64_ifort import (
                   array_permutation_fort_parallel , 
                   array_permutation_fort_parallel_complex, 
                   #I have several versions of permute_player, the parallel dynamic is supposed to be fattest
                   #but it is not thoroughly tested against too mainy cases
                   permute_player_fort_parallel_dynamic as permute_player_fort, 
                   complex_permute_player_fort_parallel_dynamic as complex_permute_player_fort
                    )
        elif os1  == 'Windows':
            from merapy.lib.win.array_permutation_gfort import (
                   array_permutation_fort_parallel , 
                   array_permutation_fort_parallel_complex, 
                   #I have several versions of permute_player, the parallel dynamic is supposed to be fattest
                   #but it is not thoroughly tested against too mainy cases
                   permute_player_fort_parallel_dynamic as permute_player_fort, 
                   complex_permute_player_fort_parallel_dynamic as complex_permute_player_fort
                    )


except ImportError:
    raise
    #from array_permutation_64 import array_permutation as array_permutation_fort



__all__ = ["array_permutation", "array_permutation_np"]

#for better performance
Dims_ = np.ndarray(32, np.int)
order_ = np.ndarray(32, np.int)
MaxSize = 22**4     #support at most 6-rank tensor 22dim for each index
#out_buff = np.ndarray(MaxSize, np.float64)


#@profile  #used for line_profiler
def array_permutation(a, rank, Dims, order, out=None):
    """  
        this wraps array_permutation through f2py, a little adjustment has been made
        a: an 1D array
        rank:  rank of a
        order: new order of a's indices
    """
    #
        #在下面的array_permutation中，按照array从1开始计数，要根据情况转换
        #if index_start==0:
        #    order=map(lambda x: x+1, order)    
        #把这一句话改成下面numpy的做法后，真个mera程序速度提高1倍!!
        
        #在f90中，dims，order都是固定长度为32的array，so I make the following
        #Dims_=np.ndarray(32,np.int)
        #Dims_=np.ndarray(rank,np.int)
    Dims_[:rank]=Dims[:rank]
    
    #order_=np.ndarray(rank, np.int)
    #order_=np.ndarray(32, np.int)
    order_[:rank]=order[:rank]
    
    if out is None: 
        
        out = np.ndarray(a.size, dtype=a.dtype)
        #out = np.ndarray(a.size,dtype=a.dtype, buffer=out_buff)
        if a.dtype == float: 
            array_permutation_fort_parallel(rank, Dims_, order_, a, out)
        else: 
            array_permutation_fort_parallel_complex(rank, Dims_, order_, a, out)
            
        return out 
    else: 
        if a.dtype == float: 
            array_permutation_fort_parallel(rank, Dims_, order_, a, out)
        else: 
            array_permutation_fort_parallel_complex(rank, Dims_, order_, a, out)


#@profile  
def array_permutation_inplace(a, rank, Dims, order, b):
    """  
        this is important for performace critical code !!
        this wraps array_permutation through f2py, a little adjustment has been made
        a: an 1D array
        rank:  rank of a
        order: new order of a's indices
    """
    #在下面的array_permutation中，按照array从1开始计数，要根据情况转换
    #if index_start==0:
    #    order=map(lambda x: x+1, order)    
    #把这一句话改成下面numpy的做法后，真个mera程序速度提高1倍!!
    
    #在f90中，dims，order都是固定长度为32的array，so I make the following
    #Dims_=np.ndarray(32,np.int)
    #Dims_=np.ndarray(rank,np.int)
    Dims_[:rank]=Dims[:rank]
    
    #order_=np.ndarray(rank, np.int)
    #order_=np.ndarray(32, np.int)
    order_[:rank]=order[:rank]

    #b=np.ndarray(a.size, dtype=a.dtype)
    #b=np.ndarray(a.size,dtype=a.dtype, buffer=out_buff)
    #print_vars(vars(),  ['a.dtype', 'b.dtype'], sep=' ')
    if a.dtype == float:  
        array_permutation_fort_parallel(rank, Dims_, order_, a, b)
    else: 
        array_permutation_fort_parallel_complex(rank, Dims_, order_, a, b)
        
def array_permutation_C(a, rank, Dims, order, index_start=0):
    """
    make array_permutation work for C order 
    """
    a1= a.reshape(Dims[:rank],order="C")
    a2= a1.ravel("F")
    a3=array_permutation(a2,rank,Dims,order)

    new_shape=[Dims[i] for i in order[:rank]]

    a4=a3.reshape(new_shape,order="F")
    a5=a4.ravel("C")
    return a5

def array_permutation_np(a, rank, Dims, order, index_start=0):
    Dims= Dims[:rank]
    order = order[:rank]
    #print "order", order,"Dims", Dims, "shape", a.shape,"size", a.size
    a=a.reshape(Dims[:rank], order="F")
    res=a.transpose(order)
    res= res.ravel(order="F")
    return res


def swapaxis_F(a,shape,ax1,ax2):
    """a is a 1-d array 
    I make numpy.swapaxes to adapt to Fortran order, inorder to compare the result with ArrayPermutation above
    
    """
    a_nd_f = a.reshape(tuple(shape),order="F")
    #a_1d_f = a_nd_f.flatten("F")
    #a_nd_f=a_1d_f.reshape(a.shape, order="F")
    
    a_swap_nd_f= a_nd_f.swapaxes(ax1,ax2)
    #a_swap_1d_f= a_swap_nd_f.flatten("F")
    a_swap_1d_f= np.ravel(a_swap_nd_f,order="F")
    return a_swap_1d_f    


#import sage
#import sage.all
"""
from sage.combinat.permutation import Permutation
print Permutation.__module__
p=Permutation([2, 3, 1, 4])
print p
v1=[3,1,4,2]
print p.reduced_word()
v1,p*v1
"""

def array_permutation_1(rank, Dims,order, totDim,A, B):
    """
    written in pure, python,  deprecated, 得到的结果不对
    see ArrayPermutation in f90

    rank是腿的个数，Dims是每条腿的维度，
    order是permutation的次序，
    A,B是存放数据的array, 
    totDim是总维度(各条腿维度之积，fortran子程序需要指定array的长度）
    
    Args:
        rorder：应该是order的逆函数
    
    """

    #ifdef __USE_OMP__
    #attention_omitted_something
    #endif    
    
    #a few checks to see whether dims match
    ntotDim = 1
    nDims = np.ndarray(rank, "int")
    rorder = np.ndarray(rank, "int")
    for i  in range(rank):
        ntotDim = ntotDim*Dims[i]
        nDims[i] = Dims[order[i]]
        rorder[order[i]] = i
    if ntotDim != totDim:  print  'error in ArrayPermutation, size not match'
    
    print "nDims",type(nDims),nDims.shape,  nDims[:rank], "\n", "Dims", Dims[:rank]

    step = np.ones(32, "int")
    rstep = np.ones(32, "int")

    for i  in range(rank):
        for j  in range(rorder[i]-1):
            step[i] = step[i]*nDims[j]
        rstep[i] = (Dims[i]-1)*step[i]
    
    print "step", step[:rank], "rstep", rstep[:rank]

    ##$    write[fmt,*] "[2[I5,1x,'| '],", rank, "[1x,I5], ' | ',", &
    ##$          rank, "[1x,I5], ' | ', 200[1x,I5]]"
    
    if 0:
        #ifdef __USE_OMP__    
        #$OMP PARALLEL private[p,q,i,npp,inc,idx,p1,p2,pp,res]    
        nths = OMP_GET_NUM_THREADS()
        nth = OMP_GET_THREAD_NUM()
        
        npp = totDim/nths
        p1 = nth*npp+1
        p2 = [nth+1]*npp
        if nth == nths-1:  p2 = totDim
        
        pp = p1-1
        for i  in range(rank-1):
            res = pp/Dims[i]
            idx[i] = pp-res*Dims[i]+1        
            pp = res
        idx[rank] = pp+1    
        
        q = idx[order[rank]]-1
        for i  in range(rank-1, 1,-1):
            q = nDims[i]*q+idx[order[i]]-1
        q = q+1

        #else:
    #p1=0   #p1 = 1 originally 
    #p2 = totDim    
    q=0     #q=1 originally
    #idx[0:rank] = 1
    idx = np.zeros(32, "int")
    #endif    
    for p  in range(totDim):
        print 'p, q, i', p, q, i
        B[q] = A[p]
        
        inc = True
        i=0 #i = 1  originally
        #下面要确定q的值
        while inc  and (i < rank): 
            idx[i] = idx[i]+1
            #if abs(q) > 100: print 'Error i q ',i, q

            if idx[i] <= Dims[i]: 
                q = q+step[i]
                inc = False
            else:
                q = q-rstep[i]
                idx[i] = 1
                i = i+1                


class Test_array_permute(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_temp(self): 
        pass 
    
    def test_array_permutation_inplace(self): 
        print array_permutation_fort_parallel.__doc__ 
        print array_permutation_fort_parallel_complex.__doc__ 
        a = np.random.random((3, 4, 2, 3))
        b = 1j*np.random.random((3, 4, 2, 3))
        out = np.ndarray((3, 4, 2, 3), dtype=float, order='F').ravel()
        array_permutation_inplace(a, 4, [3, 4, 2, 3], [0, 2, 1, 3], out)
        
        #test complext dtype 
        ab = a + b 
        out_c = np.ndarray((3, 4, 2, 3), dtype=complex, order='F').ravel()
        array_permutation_inplace(ab, 4, [3, 4, 2, 3], [0, 2, 1, 3], out_c)
        
    
    def test_0(self):
        """
            test for rank 2 and see
            wether array_permutation for rank2 is just a transposetion
        
        """
        
        rank = 2 
        Dims = [2, 4 ]
        order=[1,0]
        totDim  = np.prod(Dims[:rank])
        a = np.ndarray(totDim, np.float64);    
        a[:]=range(totDim)        
        print a.size, a.shape
    
    def test_1(self):
        """
            test for rank 2 and see
            wether array_permutation for rank2 is just a transposetion
        
        """
        #print test1.__doc__
        rank = 2 
        Dims_ = [2, 4 ]
        shape=tuple(Dims_)
        Dims=np.ndarray(32,'i')
        Dims[:rank]=Dims_[:rank]
        #order_=[2,1]
        order_=[1,0]
        
        order=np.ndarray(32,'i')
        order[:rank]=order_[:rank]

        totDim  = np.prod(Dims[:rank])
        #a = np.ndarray((1,totDim),'d');    
        a = np.ndarray(totDim,'d');    
        a[:]=range(totDim)        
        af= np.ndarray(totDim,'d',order="F")
        af[:]=range(totDim)
        #b=np.ndarray((1,totDim),'d',order="FORTRAN")
        #b=np.ndarray(totDim,'d',order="FORTRAN")
        b=np.ndarray(totDim,'d',order="F")
        b1=np.ndarray((1,totDim),'d',order="FORTRAN")
        
        b=array_permutation(a,rank, Dims, order)
        
        print "let me start an array a:\n",a
        print "af:\n", af
        print "there seems no apearing difference with the two defs"
        print "shape of a:\n",a.shape
        print "I want to reshape it into a matrix with dims", Dims_
        print "and compare between C and F order"
        a_24c= a.reshape(Dims_)
        #a_24f=a.reshape(Dims_,'F')
        a_24f=a.reshape((2,4),order='F')
        a_42f=a.reshape(Dims_[::-1])
        af_24c=af.reshape(Dims_, order='F')
        af_24f=af.reshape(Dims_, order='F')
        print "reshape of a in C order a_24c:\n" ,a_24c
        print "reshape of af in F order af_24c:\n" ,af_24f        
        print "reshape of af in F order af_24f:\n" ,af_24f                
        print "reshape of a in F order a_24f:\n" ,a_24f
        print "reshape of a in F order a_42f:\n" ,a_42f
        print u"So, the truth is  reshape 与定义a的order无关，而与reshape(,order)有关"

        print "----------------------------------"
        print "array_permutation a--> b in this order: " , order[:rank]        
        print "b\n",b
        print "now I reshape b"        
        b_24c=b.reshape(Dims_,order='C')
        b_24f=b.reshape(Dims_,order='F')
        print "b_24c:\n",b_24c
        print "b_24f:\n",b_24f
        #print "b1:\n",b1.reshape((2,4))
        #print "b1:\n",b1.reshape((4,2))

        print "-----------------------------------"
        print "now use numpy.swapaxes do the samething"
        a_24c_swap= a_24c.swapaxes(0,1)
        print "a_24c_swap:\n", a_24c_swap
        print "which is just transpose of a_24c"
        print "a_24c.T:\n", a_24c.T
        print "now convert a_24c_swap to 1-D, yield a_24c_swap_linear:"
        a_24c_swap_linear= a_24c_swap.reshape(totDim)
        print a_24c_swap_linear
        print "so it not eqal b, next see what if I a_24c_swap is reshaped in F order?"
        a_24c_swap_linear_f= a_24c_swap.reshape(totDim,order="F")
        print "a_24c_swap_linear_f:\n",  a_24c_swap_linear_f
        print "it's totally not what I want, So I faild, but next try:"


        print "--------------------------------------------"
        a_24f_swap= a_24f.swapaxes(0,1)
        print "a_24f_swap:\n", a_24f_swap
        a_24f_swap_linear= a_24f_swap.reshape(totDim)
        print a_24f_swap_linear
        print "so it not eqal b, next see what if I a_24f_swap is reshaped in F order?"
        a_24f_swap_linear_f= a_24f_swap.reshape(totDim,order="F")
        print "a_24f_swap_linear_f:\n",  a_24f_swap_linear_f
        print "is just b.  So I GET IT!  "
        print "-==================================="
        print r""" 用np.swapaxes得到和array_permutation 相同结果点方法是：
        a_1d= np.ndarray(totDim,order=)
        b_1d= array_permutation(a_1d,norder=[ i<-->j ])
        
        则，可以用如下方法得到相同结果：
        a_ndf= a_1d.reshape(Dims, order="F")
        b_ndc
        a_swap_nd_c= a_ndf.swapaxes(axi-1,axj-1)
        b_1df
        a_swap_1d_f=b_ndc.reshape(totDim,order='F')
        则
        b_1d= a_swap_1d_f
        
        总的来说，上面操作相当于定义了一个oder=“F”的swap
        def swapaxis_F(a,shape,ax1,ax2):
           a is a 1-d array 
            a_nd_f = a.reshape(tuple(shape),order="F")
            #a_1d_f = a_nd_f.flatten("F")
            #a_nd_f=a_1d_f.reshape(a.shape, order="F")
            
            a_swap_nd_f= a_nd_f.swapaxes(ax1,ax2)
            #a_swap_1d_f= a_swap_nd_f.flatten("F")
            a_swap_1d_f= np.ravel(a_swap_nd_f,order="F")
            return a_swap_1d_f    
        in the following I test this function
        """
        print swapaxis_F(a,Dims_,0,1)

    def test_2(self):
        """
            a test for rank 4 case ---pass
        """

        #rank = 4 
        Dims_ = [2, 3, 4,2 ]
        rank=len(Dims_)
        Dims=np.ndarray(32,'i')
        Dims[:rank]=Dims_[:rank]
        #order_ = np.random.permutation(rank)
        #order_=[1,3,2,4]
        #order_=[2,1,4,3]
        #order_=[2,3,1,4]
        order_=[1,2,0,3]
        order=np.ndarray(32,'i')
        order[:rank]=order_[:rank]
        totDim  = np.prod(Dims[:rank])
        a = np.ndarray(totDim,'d');    
        a[:]=range(totDim)
        b=np.ndarray(totDim,'d',order="F")
        #b=array_permutation(a,rank, Dims, order,   index_start=1)
        b=array_permutation(a,rank, Dims, order)
        print "a:\n",a.shape 
        print "a_2342\n",a.reshape(Dims[:rank])
        print "order", order[:rank]        
        print "b:\n",b
        dim1=[3,2,4,2]
        b1= swapaxis_F(a,Dims_,0,1)
        b1=swapaxis_F(b1,dim1,1,2)
        print "b1:\n", b1
        print np.all(b==b1)
    
    @unittest.skip('this func is very old, there is some problem in it,  I have no time to figure out')
    def test_3(self):
        """ 
            I want to get a array_permutation for C order, 
            make use of original array_permutation for F order
        """
        shape = (2,5)
        a = np.arange(np.prod(shape))
        print "a", a
        norder = [1,0]
        a25 = a.reshape(shape)
        b = a25.swapaxes(0,1)
        print "after swapaxes of a:\n", b.ravel()
        print "next to see which of the following array_permutation yields b"
        
        print '=='*20
        print "start from a"
        a_ap=array_permutation(a, 2, shape,norder)
        print "after array_permutation, a_ap:\n", a_ap
        nshape=(5,2)
        print "a_ap in new shape", a_ap.reshape(nshape)
        print "a_ap in new shape and F order:", a_ap.reshape(nshape,order='F')
        print a_ap.reshape(nshape,order='F').ravel("F")
        
        print "--"*20
        print "start from a25fc"
        a25f=a.reshape(shape,order="F")
        xx=a25f.ravel("C")
        yy=array_permutation(xx,2,shape,norder)
        yyF=yy.reshape(nshape,order="F")
        yyC=yy.reshape(nshape,order="C")
        print xx,"\n", yy,"\n", yyF,"\n", yyC
        print yyF.ravel("C")
        print yyC.ravel("F")

        print "--"*20
        print "start from a25cf"
        a25c=a.reshape(shape)
        xx=a25c.ravel("F")
        yy=array_permutation(xx,2,shape,norder)
        yyF=yy.reshape(nshape,order="F")
        yyC=yy.reshape(nshape,order="C")
        print xx,"\n", yy,"\n", yyF,"\n", yyC
        print yyF.ravel("C")

        print "So, I get it!, it is yyF, yyF=b"


        """
            let me summarize:
            a--> a.reshape(shape)-->a.ravel("F") -->array_permutation
            reshape(nshape,order="F")-->ravel("C")
        
        """
    
    def test_4(self):
        print "test of function array_permutation_C  ---pass"
        def sub_test1():
            shape=(4,5)
            a=np.arange(np.prod(shape),dtype="d")
            norder=[1,0]
            a25=a.reshape(shape)
            b=a25.swapaxes(0,1)
            print "after swapaxes of a:\n", b.ravel()

            c=array_permutation_C(a,2,shape,norder)
            print c
        def sub_test2():
            shape=(4,5,2)
            a=np.arange(np.prod(shape),dtype="d")
            norder=[1,0,2]
            a25=a.reshape(shape)
            b=a25.swapaxes(0,1)
            print "after swapaxes of a:\n", b.ravel()

            c=array_permutation_C(a,len(shape),shape,norder)
            print c
   
    def test_np(self):
        """
        np.transpose is equivalent to array_permutation !!----pass
        """
        from merapy.lib.array_permutation_64_ifort import array_permutation_fort
        print array_permutation_fort.__doc__
        #print array_permutation_player_fort.__doc__
        dims= [5, 3, 4, 1]
        rank = len(dims)
        order = [3, 0, 2, 1]
        a = np.ndarray(np.prod(dims), order="F")
        a[:] = np.arange(np.prod(dims))
        a1 = a.copy()
        print "original array"
        print a
        b=array_permutation(a,rank, dims, order)
        c = array_permutation_np(a1.reshape(dims, order="F"),rank, dims, order,  index_start=0)
        print b
        print c
    

if __name__ == "__main__":

    
    if 0: 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        
        suite.addTest(Test_array_permute('test_temp'))
        #suite.addTest(Test_array_permute('test_0'))
        #suite.addTest(Test_array_permute('test_1'))
        #suite.addTest(Test_array_permute('test_2'))
        #suite.addTest(Test_array_permute('test_3'))

        unittest.TextTestRunner().run(suite)








