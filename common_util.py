#coding=utf8

"""
q:
    q131, svd results could be different each time running?
"""
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
from past.utils import old_div
from builtins import object
import os, unittest
import sys 
import pickle as pickle
import numpy as np
from scipy import linalg
import platform
import socket
import mkl 

from scipy.linalg.blas import dgemm, sgemm, zgemm
try:
    import cupy as cp
    #from cupy.cuda import cublas
    import cupy.cublas as cublas 
except:
    pass


from merapy.utilities import print_vars

#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties

arch=  platform.architecture()[0]   #32 or 64 
hostname = socket.gethostname()
os1 = platform.system()
IS_PY3 = sys.version_info.major>2

if 1:
    if os1 == 'Linux':
        if IS_PY3:
            #from merapy.lib.linux_py3.common_64_ifort import *
            try:
                from merapy.lib.linux_py3.common import *
            except ImportError:
                from merapy.lib.linux_py3.common_gfort import *
        else:
            from merapy.lib.linux_py2.common_64_ifort import *
        #from merapy.lib.common_64_ifort import *
    elif os1 == 'Windows':
        from merapy.lib.win.common_gfort import *
    
    _set_num_of_threads = set_num_of_threads
else:
    _set_num_of_threads is None

def set_num_of_threads(n, info=1):
    """
        a simple wrapper 
        refs:
            https://stackoverflow.com/questions/29559338/set-max-number-of-threads-at-runtime-on-numpy-openblas
            https://github.com/numpy/numpy/issues/11826
                    It depends on the BLAS that is being used, but in
                    most cases, yes it is possible. Using ctypes it
                    shouldn't be hard, especially if you know which
                    BLAS you are using. Putting an as reliable as
                    possible utility function into numpy has its tricky
                    parts (plus for ATLAS and maybe Accelerate it just
                    is not possible to adjust threads).  But, most
                    users end up using either OpenBLAS or MKL (via
                    Anaconda) and if you do not aim for perfect
                    reliability I think creating a function that works
                    for your system is very straight forward.
                    
            https://github.com/joblib/threadpoolctl
                    threadpoolctl does not attempt to limit the size of
                    Python multiprocessing pools (threads or processes)
                    or set operating system-level CPU affinity
                    constraints: threadpoolctl only interacts with
                    native libraries via their public runtime APIs.            
                
        other ways:
            import mkl
            mkl.set_num_threads(n)
    """
    if info:
        print('set_num_of_threads to %d'%n)
    if _set_num_of_threads is not None: 
        _set_num_of_threads(n)
    else:
        mkl.set_num_threads(n)


def gemm_all(a, b, c, alpha=1.0, beta=0.0, dtype=float,  use_gpu=False, ):
    
    if not use_gpu:
        if dtype == float:
            dgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
        else:
            zgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
    else:
        cublas.gemm('N', 'N', a, b, out=c, alpha=1.0, beta=beta) 



class test_common(object):
    def __init__(self):
        pass
    def test_add_inplace(self):
        print(test_add_inplace.__doc__)
        data1 = np.ndarray(10)
        data = np.ndarray((2, 5))

        data1[:] = 1.0
        data[:] = 0.0
        #data3[:] = 5.0
        #data2 = data[10:]
        data2 = data.reshape((1, 10), order="F")
        print(data2.base is data)
        
        test_add_inplace(data1, data2)
        print(data2)
        print(data)
    #@timer
    def contract_core_player_fort(self):
        #print contract_core_player_fort.__doc__
        data1 = np.ndarray(8)
        data2 = np.ndarray(8)
        data3 = np.ndarray(8)
        data1[:4] = 2; data1[4:] = 3; 
        data2[:4] = np.identity(2).ravel(); data2[4:] = np.identity(2).ravel()

        data3[:] = 0.0
        #rec = np.arange(6)
        #rec = np.array([1, [1, 1, 1, 2, 2, 5]])
        n = 12
        rec = np.ndarray((n, 6), np.int)
        #contract_record_1[ind_count][:] = (p1, p2, p3, Dim1, Dim2, Dimc)
        rec[:n//2][:] = [0, 0, 0, 2, 2, 2]
        rec[n//2:][:] = [4, 4, 4, 2, 2, 2]
        #print rec
        import common_64_ifort as c64
        #from common_64_ifort import contract_core_player_fort_paralell_ordered, contract_core_player_fort_paralell_bac
        if 1:
            for i in range(1):
                contract_core_player_fort(data1, data2, data3, rec)
                print(data3)
        else:
            for i in range(10000):
                #c64.contract_core_player_fort_paralell_test(data1, data2, data3, rec)
                c64.contract_core_player_fort_paralell_reduction_1(data1, data2, data3, rec)
                print(data3)
    
    def set_num_of_threads(self, n):
        import common_64_ifort as c64
        c64.set_num_of_threads(n)

    @staticmethod
    def matrix_svd_unitilize():
        """   
             --- not sure 
        q131
        """
        import numpy as np
        from scipy.linalg import svd
        print(matrix_svd_unitilize.__doc__)
        a = np.arange(30).reshape(6, 5)

        a1, s= matrix_svd_unitilize(a)
        print('a1', a1)
        print('s', s)

        print("compare results with numpy svd")
        print("--"*20)
        u, s, v=svd(a, full_matrices=False)
        print("uv", u.dot(v))
        print('s', s)

    @staticmethod
    def matrix_svd_unitilize_large_mat():
        """   
             --- not sure 
        q131
        """
        import numpy as np
        from scipy.linalg import svd
        print(matrix_svd_unitilize.__doc__)
        n = 16; m = 4
        a = np.arange(n*m).reshape(n, m)

        if 1:
            a1, s= matrix_svd_unitilize(a)
            print('a1', a1)
            print('s', s)

        print("compare results with numpy svd")
        print("--"*20)
        #u, s, v=svd(a, full_matrices=False)
        #print "uv", u.dot(v)
        #print 's', s

    @staticmethod
    def matrix_svd_unitilize_1():
        """   
             --- not sure_ 
        q131
        """
        from scipy.linalg import svd
        import numpy as np
        print(matrix_svd_unitilize.__doc__)
        a = np.arange(4).reshape(2, 2)
        #a = np.arange(4).reshape(4, 1)
        a[:] = 0.0
        for i in range(5):
            a1, s= matrix_svd_unitilize(a)
            print('a1', a1)
            print('s', s)

        print("compare results with numpy svd")
        print("--"*20)
        u, s, v=svd(a, full_matrices=False)
        print("uv", u.dot(v))
        print('s', s)

    @staticmethod
    def matrix_multiply():
        import numpy as np
        print(matrix_multiply.__doc__)
        print("test matrix_multiply  --- pass")
        #a=np.arange(3*4,dtype='d').reshape((3,4))
        #b=np.arange(4*5,dtype='d').reshape((4,5))
        #c=np.empty((3,5),dtype='d')
        
        #attention must reshape it to 2D
        for i in range(1):
            m=1000;  k= 500;  n = 1000
            a=np.arange(m*k, dtype='d').reshape((m,k))
            
            b=np.arange(k*n,dtype='d').reshape((k,n))
            c=np.empty(1,dtype='d')
            #c_out=matrix_multiply(a,b,c,1.0,0.0)
            c_out=matrix_multiply(a,b,1.0,0.0)
            c_np=a.dot(b)
        print(c_out.round(10), "\n")
        print(c_np, "\n")
        print(c.round(10))

        print(np.all(c_out==c_np))
    

    @staticmethod
    def matrix_multiply_inplace_1():
        import numpy as np
        print(matrix_multiply_inplace.__doc__)
        print("test matrix_multiply  --- pass")
        #a=np.arange(3*4,dtype='d').reshape((3,4))
        #b=np.arange(4*5,dtype='d').reshape((4,5))
        #c=np.empty((3,5),dtype='d')
        
        #attention must reshape it to 2D
        m=100
        a=np.ndarray(m,dtype='d',order="F")#.reshape((100,1),order="F")
        n= 50
        b=np.ndarray(n,dtype='d',order="F").reshape((50,1),order="F")
        
        Cin=np.ndarray(200,order="F")
        
        matrix_multiply_inplace(a,b,Cin,1.0,0.0)
        c_np=a.dot(b)
        #print c_np, "\n"
        #print c.round(10)

        print(np.all(Cin==c_np))
    
    
    def test_set_matrix(self):
        """----pass """
        print(set_matrix.__doc__)
        a = np.zeros((5, 5), order='F')
        b = np.ones((3, 3), order='F')
        x, y = 1, 1
        #set_matrix(a, b, 1, 1, forward=False)
        set_matrix(a, b, x, y, forward=False)
        print(a)
        print(b)
        print("compare with np  =================")
        a1 = np.zeros((5, 5), order='F')
        b1 = np.ones((3, 3), order='F')
        #x1, y1 = b.shape
        #a1[x:x+x1, y:y+y1] = b


class TestCommon(unittest.TestCase):
    def setUp(self): 
        pass
    def test_temp(self) : 
        #print array_permutation_fort_parallel.__doc__ 
        print(matrix_multiply_complex.__doc__) 

    
    def test_matrix_svd_1by1(self) : 
        #print(matrix_svd.__doc__) 
        a = np.random.random((1, 1)) 
        a_orig = a 
        a_copy = a.copy()
        print_vars(vars(),  ['a', 'a_copy', 'a_orig'])
        print_vars(vars(),  ['id(a_orig)', 'id(a_copy)'])
        
        u, s, v=matrix_svd(1, a) 
        print_vars(vars(),  ['u', 's', 'v'])
        print(a, id(a))
        print_vars(vars(),  ['a', 'a_copy', 'a_orig'])
        print_vars(vars(),  ['id(a_orig)', 'id(a_copy)'])
        
        print("下面结果会是false，为了提醒 maxtrix_svd 对于1by1 矩阵有问题")
         
        if IS_PY3:
            #update,  it seems the bug is fixed for py3 
            #self.assertTrue(a_orig[0, 0]==a_copy[0, 0])
            pass
        else:
            self.assertFalse(a_orig[0, 0]==a_copy[0, 0])
            
       
    def test_get_num_of_threads(self): 
        set_num_of_threads(5)
        print(get_num_of_threads.__doc__)
        print(get_max_threads.__doc__)
        print('num_of_threads', get_num_of_threads())
        print('max_threads', get_max_threads())
        print('num_procs', get_num_of_procs())
        import multiprocessing 
        print(multiprocessing.cpu_count())
        import threading
        print(threading.active_count())

    def test_get_position_and_rev(self):
        """test both position and position_rev"""
        
        mgp = matrix_get_position  
        mpgr = matrix_get_position_rev  
         
        rank=3
        #dims=np.array([2,3,4])
        dims=[2,2,2]
        pos=[1,1,1]
        pos=[0,0,1]

        posr=np.ndarray(rank)
        print(mpgr.__doc__)
        for i in range(8):
            posrrr=mpgr(i, dims)
            print('iii', i, posrrr, mgp(posrrr,dims))
    
    def test_matrix_svd(self): 
        print(matrix_svd.__doc__) 
        #n, m = 200, 80
        n, m = 100, 80
        a = np.random.random((n, m))
        u, s, v=matrix_svd(m, a)
        print(s) 
        print(u.shape, s.shape, v.shape)
        print((u*s).dot(v) - a)
        
        U, S, V=np.linalg.svd(a, full_matrices=0)
        print(S) 
        val, vec = np.linalg.eigh(a.T.dot(a)) 
        #print S**2-val[::-1]
        self.assertTrue(np.allclose(s, S, atol=1e-15))

    def test_matrix_eigen_vector(self):
        import numpy as np
        from scipy.linalg import eig
        print(matrix_eigen_vector.__doc__)
        n =  5
        vec = np.ndarray(5, order="F")
        #vec = np.asfortranarray(vec)
        a = np.random.random((n, n))
        b = np.random.random((n, n))
        a = a + a.T
        b[:, :] = a
        a = np.asfortranarray(a)
        
        matrix_eigen_vector(a, vec)

        print("a", a, vec[:n])

        print("compare with np ==========")
        val, vec=eig(b)
        print(val)
        print(vec)
        
    def test_matrix_multiply(self):
        if 1: 
            print(matrix_multiply.__doc__)
            #attention must reshape it to 2D
            m=15;  n= 20 
            a=np.arange(m,dtype='d').reshape((3,5))
            b=np.arange(n,dtype='d').reshape((5,4))
            c=np.empty(1,dtype='d')
            #c_out=matrix_multiply(a,b,c,1.0,0.0)
            c_out=matrix_multiply(a,b,1.0,0.0)
            c_np=a.dot(b)
            print(c_out.round(10), "\n")
            print(c_np, "\n")
            self.assertTrue(np.all(c_out==c_np))
        
        if 1:  #dtype = complex 
            print(matrix_multiply_complex.__doc__)
            #attention must reshape it to 2D
            m=15; n= 20 
            a=np.arange(m,dtype=complex).reshape((3,5))   +  1j
            b=np.arange(n,dtype=complex).reshape((5,4))  +  1j 
            #a=np.arange(m,dtype=float).reshape((3,5))
            #b=np.arange(n,dtype=float).reshape((5,4)) 
            c=np.empty(1,dtype=complex)
            #c_out=matrix_multiply(a,b,c,1.0,0.0)
            c_out=matrix_multiply_complex(a,b,1.0,0.0)
            c_np=a.dot(b)
            print(c_out.round(10), "\n")
            print(c_np, "\n")
            print(c.round(10))
            self.assertTrue(np.all(c_out==c_np))
        
    def test_matrix_multiply_inplace(self):
        print(matrix_multiply_inplace.__doc__)
        if 1:  
            #attention must reshape it to 2D
            m=3
            k=5
            n=8
            a=np.arange(m*k,dtype='d').reshape((m,k),order="F")
            b=np.arange(n*k,dtype='d').reshape((k,n),order="F")
            
            Cin=np.empty((m,n),order="F")
            
            
            matrix_multiply_inplace(a, b,Cin,1.0,0.0  )  #,3,3,1)
            c_np=a.dot(b)
            print(Cin[:])
            self.assertTrue(np.all(Cin==c_np))
        
        if 1:  #inplace version   
            #attention must reshape it to 2D
            print(matrix_multiply_inplace_complex.__doc__)
            m=7
            k=5
            n=3
            a=np.arange(m*k,dtype=complex).reshape((m,k),order="F")
            a += 1j*a  
            b=np.arange(n*k,dtype=complex).reshape((k,n),order="F")
            b  += -0.5j*b  
            
            Cin=np.empty((m,n),order="F", dtype=complex)
            
            matrix_multiply_inplace_complex(a, b,Cin,1.0 + 0j, 0.0  )  #,3,3,1)
            c_np=a.dot(b)
            print(Cin[:])

            self.assertTrue(np.all(Cin==c_np))


if __name__=="__main__":
    
    def test_common_func():
        import numpy as np

        def test_unit_matrix():
            from common_32 import unit_matrix
            #print c32.__doc__
            print("test on unit_matrix: ---pass")
            print(unit_matrix.__doc__)
            a=np.zeros((4,4),"d",order="FORTRAN")
            #a=np.zeros((4,4),"d")
            unit_matrix(a)
            print(a)
        #test_unit_matrix()

        def test_transpose():
            from common_32 import transpose4py
            print(transpose4py.__doc__)
            
            a=np.arange(16).reshape((4,4))
            b=np.zeros((4,4),'d')
            b=transpose4py(a,b)
            print(a,"\n",b)
        #test_transpose()

        def test_matrix_direct_product():
            from common_32 import matrix_direct_product 
            print(matrix_direct_product.__doc__)
            print("test matrix_direct_product: ---pass")
            a=np.arange(3).reshape((3,1))
            b=np.arange(5).reshape((1,5))
            c=matrix_direct_product(a,b)
            print(c)

            a = np.arange(15)
            b = np.arange(8)
            #following won't work, we must reshape it 
            #c = matrix_direct_product(a, b, 3, 5, 2, 4)
            #print c.shape
        #test_matrix_direct_product()


        def test_matrix_trace():
            """ ---pass """
            print(matrix_trace.__doc__)
            a=np.arange(100).reshape(10,10)
            b=matrix_trace(a)
            print(b, a.trace())
        #test_matrix_trace()
        def test_svd():
            """  ---pass """
            print(matrix_svd.__doc__)
            m=8; n=70  ; mn=min(m,n)
            a=np.random.rand(m*n).reshape(m,n)
            b=matrix_svd(mn,a)
            c= np.linalg.svd(a, full_matrices=True)
            print(b[1].round(10))
            print(c[1].round(10))
        #test_svd()
        def test_matrix_svd_unitilize():
            """   
                 --- not sure 
            q131
            """
            from scipy.linalg import svd
            print(matrix_svd_unitilize.__doc__)
            a = np.arange(30).reshape(6, 5)

            a1, s= matrix_svd_unitilize(a)
            print('a1', a1)
            print('s', s)

            print("compare results with numpy svd")
            print("--"*20)
            u, s, v=svd(a, full_matrices=False)
            print("uv", u.dot(v))
            print('s', s)
        #test_matrix_svd_unitilize()

        def test_iShift():
            """  ---not sure"""
            a=np.arange(10,dtype="int")
            #print ishift.__doc__
            ishift(a,4,8)
            print(a)
        #test_iShift()


    if 0: #examine
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCommon)
        unittest.TextTestRunner(verbosity=0).run(suite)    
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
            #'test_temp', 
            #'test_matrix_multiply', 
            #'test_matrix_multiply_inplace', 
            #'test_get_num_of_threads', 
            'test_matrix_svd_1by1', 
            #'test_matrix_svd', 
            #'test_matrix_eigen_vector', 
        ]
        for a in add_list: 
            suite.addTest(TestCommon(a))
        
        unittest.TextTestRunner().run(suite)
       






