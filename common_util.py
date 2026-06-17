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

from scipy.linalg.blas import dgemm, sgemm, zgemm, cgemm 
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

if 0:
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
    _set_num_of_threads  = None

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


            
def gemm_all(a, b, c, alpha=1.0, beta=0.0, dtype=np.float64,  use_gpu=0, transfer_data=False):

    if use_gpu == 0:
        if dtype == np.float64:
            dgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
        elif dtype == np.complex128:
            zgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
        elif dtype == np.float32:
            sgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
        elif dtype == np.complex64:
            cgemm(alpha, a, b, beta=beta, c=c, overwrite_c=True)
        else:
            raise TypeError(f"不支持当前的物理数据类型: {dtype}")
            
    elif use_gpu == 2:
        a = cp.asarray(a, dtype=dtype)
        b = cp.asarray(b, dtype=dtype)
        out = cublas.gemm('N', 'N', a, b, alpha=alpha, beta=beta) 
        #out.get(order='F', out=c)
        c += out.get(order='F') 
    elif use_gpu == 1:
        #assert isinstance(a, cp.ndarray), "use_gpu=1 时，输入必须已经是 CuPy 显存数组"
        #print_vars(vars(),  ['a.dtype', 'b.dtype', 'c.dtype'])
        cublas.gemm('N', 'N', a, b, out=c, alpha=alpha, beta=beta) 
    else:
        raise ValueError(use_gpu) 

def examine_cupy_array(t, ):
    """ 
        Exammine a cupy ndarray memory
    """
    import cupy as cp
    if not isinstance(t.data, cp.ndarray):
        return
        
    # 1. 提取低级 CUDA 内存指针对象
    memptr = t.data.data  # cupy.cuda.MemoryPointer
    mem_block = memptr.mem  # cupy.cuda.Memory
    
    ptr_start = memptr.ptr        # 数据的起始绝对物理地址
    element_size = t.data.itemsize # 每个元素占的字节数 (complex64=8, complex128=16)
    view_bytes = t.data.size * element_size # 视图期望的字节大小
    
    print(f"\n==== 🔍 显存审计报告: {t.name} ====")
    print(f"数据类型 (dtype)  : {t.data.dtype}")
    print(f"形状与步长(Shape/Strides): {t.data.shape} / {t.data.strides}")
    print(f"元素总数 (size)   : {t.data.size}")
    print(f"视图占据字节(Bytes) : {view_bytes} 字节")
    print(f"物理首地址 (Pointer): {hex(ptr_start)}")
    print(f"底层内存块总大小    : {mem_block.size} 字节")
    print(f"内存物理连续性标志  : C_Contig={t.data.flags['C_CONTIGUOUS']}, F_Contig={t.data.flags['F_CONTIGUOUS']}")
    
    # 核心物理对账：你承诺的视图大小，绝对不能超过底层物理块的大小
    if view_bytes > mem_block.size:
        print("🚨 [硬件级越界] 严重警告：该视图需要的字节数超过了底层 CUDA 物理分配块！")


class TestCommon(unittest.TestCase):
    def setUp(self): 
        pass
    
    def test_temp(self) : 
        #print array_permutation_fort_parallel.__doc__ 
        pass
       
    def xtest_get_num_of_threads(self): 
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


    if 1: #examine
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCommon)
        unittest.TextTestRunner(verbosity=0).run(suite)    
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
            #'test_temp', 
        ]
        for a in add_list: 
            suite.addTest(TestCommon(a))
        
        unittest.TextTestRunner().run(suite)
       






