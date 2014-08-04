#coding=utf8

"""
q:
    q131, svd results could be different each time running?
"""
import os, unittest
import cPickle as pickle
import numpy as np
import platform
import socket
#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties

arch=  platform.architecture()[0]

hostname = socket.gethostname()


#if hostname == 'VirtualBox-Lab' and 0:
#    from merapy.lib.common_64_ifort_virtual import *
#else:
#    from merapy.lib.common_64_ifort import *

from merapy.lib.common_64_ifort import *


del arch
del platform


"""
sometimes fort func don't require passing in "F" ordered arrays, such as matrix_direct_product; while someitmes it do  why?

"""

class PickleDb:
    """
        not implemented
    """
    def __init__(self, fn=None, path=None):

        if fn is not None:
            path = os.path.abspath(fn); path = fn[-1::-1];  n = path.find("/");   
            path = path[n:]; path = path[-1::-1]
            self.fn_only = path[:n]
        elif path is not None:
            x = self.find_max_dim_file(path)
            print "file with max layer and trunx_dim is %s"%x
            #fn=path + x
            self.fn_only = x
            fn = "%s%s"%(path, x) if path[-1]=="/" else "%s/%s"%(path, x)
        else:
            raise Exception
        self.path = path
        #print fn, path
        inn = open(fn, "rb")
        self.obj = pickle.load(inn)
    
    def find_max_dim_file(self, path):
        nlist = os.listdir(path)
        #print nlist
        pickle_files = [i for i in nlist if "pickle" in i]
        print "all pickle files under specified path are %s"%pickle_files
        
        def _parse_name(fn):
            if "transition" in fn: #'12-transition.pickle'
                dim = "0"
            elif "_" in fn:  #like "5_13.pickle or 5_13-4lay.pickle"
                dim = fn.split("-")[0].split("_")[1].split(".")[0]
            else:
                if "lay" in fn: #"like "8-4lay.pickle""
                    dim = fn.split("-")[0]
                else:   
                     #like "8.pickle"
                    dim = fn.split(".")[0]
            
            i = fn.find("lay")
            layer = fn[i-1] if i>0 else "3"
            #print "ddd", [dim, layer]
            return eval(dim), eval(layer)

        temp = []
        for i in pickle_files:
            #print i, _parse_name(i)
            try:
                dim, layer = _parse_name(i)
            except:
                raise Exception(i)
            temp.append((layer, dim, i))
        temp.sort(key=lambda x:x[:2], reverse=True)
        res= temp[0][2]
        return res
        
    def _plot(self, *args, **kwargs):
        plt.plot(*args, **kwargs)
        plt.grid(True)
        plt.show()
    
    def _scaling_dim_exact(self):
        exact = {}

        exact["xx"] = {
                0:(0.0, ) + (1.0, )*4   +  (2.0, )*8, 
                1:(0.25,) + (1.25, )*4  +  (2.25, )*8 , 
                2:          (1.0, )*4  + (2.0, )*8}
        
        exact["Ising"] = {
                1:(0, ) + (1, )   +  (2, )*4, 
                -1:(1/8.,) + (1 + 1/8., )*2  +  (2+1/8., )*3 
                }
        exact["heisbg_long"] = {
                0:(0, ) + (0.5, )*2  + (1.0, )*2  + (1.5, )*2, 
                1:(0.5, ) + (1, )*2  + (1.5, )*2, 
                2:(2.0, )*3, 
                }
        exact["heisbg"] = exact['heisbg_long']

        ee = self.obj.energy_exact
        mapper = {-1.2732:"xx", -1.6077:"heisbg_NNN", -1.6449:"heisbg_long", -1.7725:'heisbg'}
        for i in mapper:
            if abs(i-ee)<1e-2:
                model = mapper[i]
        exact_val = exact[model]
        if self.obj.model == "Ising": 
            exact_val = exact["Ising"]

        return exact_val

    def plot_scaling_dim(self, qns=None, num_of_fields=None, calc_err=False):
        """
            one rec is like this:
            {1025:
                [(0, [0.0, 1.0, 1.16326228, 1.19447544, 1.19447544, 1.23058869, 1.6277161, 1.6277161, 1.82188553, 1.82188553]),
                 (1, [0.87008947, 1.08398164, 1.12192732, 1.50196967, 1.63179961, 1.74482795, 1.88535985, 2.0160677, 2.13870443, 2.13870443]),
                 (2, [1.0722275, 1.78260076, 2.37668506, 2.60884736, 2.60884736, 2.800528, 2.800528, 2.85390176, 2.93886244, 2.93886244])
                 ]
             } or
             {1025:
                [(1, [-0.0, 1.00009277, 1.99996219, 1.99997006, 2.00000392, 2.00289768, 2.72124425, 2.72124425, 2.72266868, 2.72266868]), 
                [(-1, [0.1250008, 1.12503901, 1.12508933, 2.12244868, 2.12717037, 2.12830034, 2.35310428, 2.35310428, 2.35354695, 2.35354695])}
            
        """
        num_of_fields = num_of_fields if num_of_fields is not None else 2
        try:
            rec = self.obj.scaling_dim_record
        except AttributeError:
            rec = self.obj.conformal_dim_record
        
        if calc_err:
            exact_val = self._scaling_dim_exact()
            print exact_val

        x = rec.keys()
        x.sort()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax1 = fig.add_subplot(111)
        
        xy = self.get_energy_err()
        ax1.plot(*xy, linestyle="dashed")

        xy = self.get_energy_diff()
        ax1.plot(*xy, linestyle="dashed")
        #ax1.xlim(xmin=x[0])
        
        if qns is None:
            if self.obj.symmetry == "U1":
                qns = range(3)
            elif self.obj.symmetry == "Z2":
                qns= [1, -1]
            else:
                qns= [1]
        fontP = FontProperties()
        fontP.set_size('small')
   
        legend = ["eng_err"]
        color_map = {0:"red", 1:"blue", 2:"green"}
        line_sty_map = {0:"solid", 1:"dashed", 2:"dotted"}

        for qn in qns:
            for i in range(num_of_fields):
                qn_ = qns.index(qn)
                if qn_ == 0 and num_of_fields== 1: i = 1
                linewidth = 3.5 if (qn_, i) in [(0, 1), (1, 0), (2, 0)] else 1.0
                y1 = np.array([(rec[ii][qn][i]) for ii in x])
                
                if calc_err:
                    try:
                        y1 = y1-exact_val[qn][i]
                        y1 = np.log10(np.abs(y1))
                    except IndexError:
                        print "index qn=%d, i=%d not found in exact_val, omit it "%(qn, i)
                        continue
                ax.plot(x, y1, "o", linewidth=linewidth,  
                        color=color_map[qn_], linestyle="solid", marker=None)
                
                #if (qn_, i) in [(0, 1), (1, 0), (2, 0)]:
                legend.append("%d, %d"%(qn, i))
        box = ax.get_position()
        box1 = ax1.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax1.set_position([box.x0, box1.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(legend, loc='center left', bbox_to_anchor=(1, 0.5), prop=fontP)
        plt.xlim(xmin=x[0])
        #plt.ylim(ymin=-0.1)
        #plt.tick_params()
        #plt.minor_tickson()
        #plt.legend(legend, loc=3)
        plt.grid(True)
        
        print "ssss\n", rec[x[-1]]
        if hasattr(self.obj, "log"):
            print "llll\n", self.obj.log
        x = "_err" if calc_err else ""
        y = self.fn_only.split(".")[0]
        fn=self.path + "scaling_dim-%s--%s.png"%(x, y)
        plt.savefig(fn ,bbox_inches='tight') # save as png
        os.popen("display " + fn)

    def plot_energy(self):
        """
            one record is like this:
            {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        """
        recs = self.obj.energy_record
        ee = self.obj.energy_exact
        x = recs.keys()
        x.sort()
        y = [recs[i][0] for i in x]
        y = np.array(y)
        y = np.log10(y-ee)
        #print recs[1]
        self._plot(x, y,"o", color='red',linestyle='dashed')

    def get_energy_err(self):
        """
            one record is like this:
            {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        """
        recs = self.obj.energy_record
        ee = self.obj.energy_exact
        x = recs.keys()
        x.sort()
        y = [recs[i][0] for i in x]
        y = np.array(y)
        y = np.log10(y-ee)
        #print recs[1]
        #self._plot(x, y,"o", color='red',linestyle='dashed')
        return x, y
    
    def plot_energy_diff(self):
        """
        one record is like this:
        {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        
        """
        recs= self.obj.energy_record
        x = recs.keys()
        x.sort()
        dy = np.array([recs[i+1][0] - recs[i][0] for i in x[:-2]])
        dy = np.log10(-dy)
        #print recs[1]
        self._plot(x[:-2], dy,"o", color='red',linestyle='dashed')

    def get_energy_diff(self):
        """
        one record is like this:
        {1:(-0.093284875794309885, 1366290999.381305, 529.42)}
        
        """
        recs= self.obj.energy_record
        x = recs.keys()
        x.sort()
        #dy = np.array([recs[i+1][0] - recs[i][0] for i in x[:-2]])
        dy = np.array([recs[x[i+1]][0] - recs[x[i]][0] for i in range(len(x)-2)])
        dy = np.log10(-dy)
        #print recs[1]
        #self._plot(x[:-2], dy,"o", color='red',linestyle='dashed')
        
        return x[:-2], dy




def set_matrix_np(a, b, x, y, forward):
    x1, y1 = b.shape
    if forward:
        a[x:x+x1, y:y+y1] = b
    else:
        b[:, :] = a[x:x + x1, y:y + y1]

def timer(func):
    import time
    from timeit import timeit
    def wraper(*args, **kargs):
        t1=time.clock()
        ta = time.time()

        res= func(*args, **kargs)
        t2=time.clock()
        tb = time.time()
        q_iter = 1
        print "cpu time", (t2-t1)/q_iter,  "\t wall time", (tb-ta)/q_iter, "\t", func.func_name
        return res 
    return wraper


class test_common():
    def __init__(self):
        pass
    def test_add_inplace(self):
        print test_add_inplace.__doc__
        data1 = np.ndarray(10)
        data = np.ndarray((2, 5))

        data1[:] = 1.0
        data[:] = 0.0
        #data3[:] = 5.0
        #data2 = data[10:]
        data2 = data.reshape((1, 10), order="F")
        print data2.base is data
        
        test_add_inplace(data1, data2)
        print data2
        print data
    @timer
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
                print data3
        else:
            for i in range(10000):
                #c64.contract_core_player_fort_paralell_test(data1, data2, data3, rec)
                c64.contract_core_player_fort_paralell_reduction_1(data1, data2, data3, rec)
                print data3
    
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
        print matrix_svd_unitilize.__doc__
        a = np.arange(30).reshape(6, 5)

        a1, s= matrix_svd_unitilize(a)
        print 'a1', a1
        print 's', s

        print "compare results with numpy svd"
        print "--"*20
        u, s, v=svd(a, full_matrices=False)
        print "uv", u.dot(v)
        print 's', s

    @staticmethod
    def matrix_svd_unitilize_large_mat():
        """   
             --- not sure 
        q131
        """
        import numpy as np
        from scipy.linalg import svd
        print matrix_svd_unitilize.__doc__
        n = 16; m = 4
        a = np.arange(n*m).reshape(n, m)

        if 1:
            a1, s= matrix_svd_unitilize(a)
            print 'a1', a1
            print 's', s

        print "compare results with numpy svd"
        print "--"*20
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
        print matrix_svd_unitilize.__doc__
        a = np.arange(4).reshape(2, 2)
        #a = np.arange(4).reshape(4, 1)
        a[:] = 0.0
        for i in range(5):
            a1, s= matrix_svd_unitilize(a)
            print 'a1', a1
            print 's', s

        print "compare results with numpy svd"
        print "--"*20
        u, s, v=svd(a, full_matrices=False)
        print "uv", u.dot(v)
        print 's', s

    @staticmethod
    def matrix_multiply():
        import numpy as np
        print matrix_multiply.__doc__
        print "test matrix_multiply  --- pass"
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
        print c_out.round(10), "\n"
        print c_np, "\n"
        print c.round(10)

        print np.all(c_out==c_np)
    
    @staticmethod
    def matrix_multiply_inplace():
        import numpy as np
        from scipy import linalg
        print matrix_multiply_inplace.__doc__
        print "test matrix_multiply  --- pass"
        #a=np.arange(3*4,dtype='d').reshape((3,4))
        #b=np.arange(4*5,dtype='d').reshape((4,5))
        #c=np.empty((3,5),dtype='d')
        
        #attention must reshape it to 2D
        m=1
        k=5
        n=1
        a=np.arange(m*k,dtype='d').reshape((m,k),order="F")
        b=np.arange(n*k,dtype='d').reshape((k,n),order="F")
        
        Cin=np.empty((m,n),order="F")
        
        print Cin[:]
        matrix_multiply_inplace(a, b,Cin,1.0,0.0  )  #,3,3,1)
        c_np=a.dot(b)
        print c_np[:], "\n"
        print Cin[:]

        print np.all(Cin==c_np)

    @staticmethod
    def matrix_multiply_inplace_1():
        import numpy as np
        print matrix_multiply_inplace.__doc__
        print "test matrix_multiply  --- pass"
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

        print np.all(Cin==c_np)
    
    def matrix_eigen_vector(self):
        import numpy as np
        from scipy.linalg import eig
        print matrix_eigen_vector.__doc__
        n =  5
        vec = np.ndarray(65536, order="F")
        #vec = np.asfortranarray(vec)
        a = np.random.random((n, n))
        b = np.random.random((n, n))
        a = a + a.T
        b[:, :] = a
        a = np.asfortranarray(a)
        
        matrix_eigen_vector(a, vec)

        print "a", a, vec[:n]

        print "compare with np =========="
        val, vec=eig(b)
        print val
        print vec




        pass

    def test_set_matrix(self):
        """----pass """
        print set_matrix.__doc__
        a = np.zeros((5, 5), order='F')
        b = np.ones((3, 3), order='F')
        x, y = 1, 1
        #set_matrix(a, b, 1, 1, forward=False)
        set_matrix(a, b, x, y, forward=False)
        print a
        print b
        print "compare with np  ================="
        a1 = np.zeros((5, 5), order='F')
        b1 = np.ones((3, 3), order='F')
        #x1, y1 = b.shape
        #a1[x:x+x1, y:y+y1] = b
        set_matrix_np(a1, b1, x, y, forward=False)
        print a1
        print b1

class TestCommon(unittest.TestCase):
    def setUp(self): 
        pass
    def test_temp(self) : 
        pass 
    
    def test_matrix_svd_1by1(self) : 
        print matrix_svd.__doc__ 
        a = np.random.random((1, 1)) 
        a_orig = a 
        a_copy = a.copy()
        print a, id(a)
        u, s, v=matrix_svd(min(a.shape), a) 
        print a, id(a)
        print  a_orig, a_copy
        print "下面结果会是false，为了提醒 maxtrix_svd 对于1by1 矩阵有问题"
        self.assertTrue(a_orig[0, 0]==a_copy[0, 0])
       
    def test_get_num_of_threads(self): 
        set_num_of_threads(5)
        print get_num_of_threads.__doc__
        print get_max_threads.__doc__
        print 'num_of_threads', get_num_of_threads()
        print 'max_threads', get_max_threads()
        print 'num_procs', get_num_of_procs()
        import multiprocessing 
        print multiprocessing.cpu_count()
        import threading
        print  threading.active_count()

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
        print mpgr.__doc__
        for i in range(8):
            posrrr=mpgr(i, dims)
            print  'iii', i, posrrr, mgp(posrrr,dims)
    
    def test_matrix_svd(self): 
        print matrix_svd.__doc__ 
        a = np.random.random((10, 8))
        u, s, v=matrix_svd(8, a)
        print  s 
        U, S, V=np.linalg.svd(a, full_matrices=0)
        print  S 
        self.assertTrue(np.allclose(s, S, atol=1e-15))
    
if __name__=="__main__":

    def test_common_func():
        import numpy as np

        def test_unit_matrix():
            from common_32 import unit_matrix
            #print c32.__doc__
            print "test on unit_matrix: ---pass"
            print unit_matrix.__doc__
            a=np.zeros((4,4),"d",order="FORTRAN")
            #a=np.zeros((4,4),"d")
            unit_matrix(a)
            print a
        #test_unit_matrix()

        def test_transpose():
            from common_32 import transpose4py
            print transpose4py.__doc__
            
            a=np.arange(16).reshape((4,4))
            b=np.zeros((4,4),'d')
            b=transpose4py(a,b)
            print a,"\n",b
        #test_transpose()

        def test_matrix_direct_product():
            from common_32 import matrix_direct_product 
            print matrix_direct_product.__doc__
            print "test matrix_direct_product: ---pass"
            a=np.arange(3).reshape((3,1))
            b=np.arange(5).reshape((1,5))
            c=matrix_direct_product(a,b)
            print c

            a = np.arange(15)
            b = np.arange(8)
            #following won't work, we must reshape it 
            #c = matrix_direct_product(a, b, 3, 5, 2, 4)
            #print c.shape
        #test_matrix_direct_product()

        def test_matrix_multiply():
            print matrix_multiply.__doc__
            print "test matrix_multiply  --- pass"
            #a=np.arange(3*4,dtype='d').reshape((3,4))
            #b=np.arange(4*5,dtype='d').reshape((4,5))
            #c=np.empty((3,5),dtype='d')
            
            #attention must reshape it to 2D
            m=100
            a=np.arange(m,dtype='d').reshape((20,5))
            n= 50
            b=np.arange(n,dtype='d').reshape((5,10))
            c=np.empty(1,dtype='d')
            #c_out=matrix_multiply(a,b,c,1.0,0.0)
            c_out=matrix_multiply(a,b,1.0,0.0)
            c_np=a.dot(b)
            print c_out.round(10), "\n"
            print c_np, "\n"
            print c.round(10)

            print np.all(c_out==c_np)
        #test_matrix_multiply()

        def test_matrix_trace():
            """ ---pass """
            print matrix_trace.__doc__
            a=np.arange(100).reshape(10,10)
            b=matrix_trace(a)
            print b, a.trace()
        #test_matrix_trace()
        def test_svd():
            """  ---pass """
            print matrix_svd.__doc__
            m=8; n=70  ; mn=min(m,n)
            a=np.random.rand(m*n).reshape(m,n)
            b=matrix_svd(mn,a)
            c= np.linalg.svd(a, full_matrices=True)
            print b[1].round(10)
            print c[1].round(10)
        #test_svd()
        def test_matrix_svd_unitilize():
            """   
                 --- not sure 
            q131
            """
            from scipy.linalg import svd
            print matrix_svd_unitilize.__doc__
            a = np.arange(30).reshape(6, 5)

            a1, s= matrix_svd_unitilize(a)
            print 'a1', a1
            print 's', s

            print "compare results with numpy svd"
            print "--"*20
            u, s, v=svd(a, full_matrices=False)
            print "uv", u.dot(v)
            print 's', s
        #test_matrix_svd_unitilize()

        def test_iShift():
            """  ---not sure"""
            a=np.arange(10,dtype="int")
            #print ishift.__doc__
            ishift(a,4,8)
            print a
        #test_iShift()


    if 0: #examine
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCommon)
        unittest.TextTestRunner(verbosity=0).run(suite)    
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
            #TestCommon('test_get_num_of_threads'), 
            #TestCommon('test_temp'), 
            TestCommon('test_matrix_svd_1by1'), 
            #TestCommon('test_matrix_svd'), 
       
        ]
        for a in add_list: 
            suite.addTest(a)
        
        unittest.TextTestRunner().run(suite)
       



    #tc = test_common()
    #tc.contract_core_player_fort()
    #tc.set_num_of_threads(2)
    #tc.test_add_inplace()
    #tc.matrix_svd_unitilize_1()
    #tc.matrix_multiply_inplace()
    #tc.matrix_multiply()
    #tc.matrix_svd_unitilize_large_mat()
    #tc.matrix_eigen_vector()
    #tc.test_set_matrix()




