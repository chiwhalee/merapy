#coding=utf8
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
from builtins import object
from past.utils import old_div
import os, sys
import unittest 
import pickle as pickle
import cProfile 
import numpy as np
import pstats 
import os.path
import warnings
#import pprocess

#from merapy import common_util
import merapy.common_util as common_util

from merapy import array_permutation
from merapy.context_util import redirect

#import tensor_py
if 0:
    from merapy.tensor_player_single import (decorate_methods, tensor_player, 
            set_STATE_end_1, set_STATE_end_simple, set_player_state_auto, 
            set_player_state_manual, get_player_state)
else:
    from merapy.tensor_player_multiple import (decorate_methods, tensor_player, 
            set_player_state_auto, set_player_state_manual, get_player_state, 
            set_STATE_end_1, set_STATE_end_simple)

num_of_instance = 0


__all__ = ["decorate_methods", "tensor_player", "set_STATE_end_1", "set_STATE_end_simple", "timer", "profileit", 'set_player_state_auto']
#__all__ = ["decorate_methods", "tensor_player", "set_STATE_end", "set_STATE_end_1", "timer", "profileit"]



if 1:
    def timer(func, func_name=None):
        import time
        from timeit import timeit
        if hasattr(func, '__name__'): 
            func_name = func.__name__
        def wraper(*args, **kargs):
            t1=time.process_time()
            ta = time.time()

            res= func(*args, **kargs)
            t2=time.process_time()
            tb = time.time()
            q_iter = 1
            print("cpu time", old_div((t2-t1),q_iter),  "\t wall time", old_div((tb-ta),q_iter), "\t", func_name)
            return res 
        return wraper

    def count(func):
        def wrapper(*args, **kargs):
            global num_of_instance
            num_of_instance  += 1
            print("num_of_instance", num_of_instance)
            #print "aaaaaaa"*100, args
            return func(*args, **kargs)
        return wrapper


    def count_instance(func):
        def wrapper(*args, **kargs):
            global num_of_instance
            num_of_instance += 1
            print("num_of_instance", num_of_instance)
            return func(*args, **kargs)
        return wrapper


    class counter(object):
        num_of_instance = 0
        def __init__(self, aClass):
            self.num_of_instance += 1 
            self.aClass= aClass
            print("num is: ",  self.num_of_instance)
        def __call__(self,  rank,  QSp, totQN, order="F",  buffer=None, use_buf=False, shallow=False):
            return self.aClass(rank=rank,  QSp=QSp, totQN=totQN, order=order,  buffer=buffer, use_buf=use_buf, shallow=shallow)


if 1:  #profilers
    def profiled(path, multi=True):
        """
        Decorator to allow individual functions to be profiled, without profiling
        the whole program.  This allows for much more targeted profiling, which is
        necessary for threaded programs or if performance becomes an issue.

        multi: if True, adds a sequential number to each profile run to prevent
                        name collisions.
               if False, the last invocation wins.
        """
        # This extra layer of indirection is so the decorator can accept arguments
        # When the user doesn't provide arguments, the first arg is the function,
        # so detect that and show an error.
        if not isinstance(path, str):
            raise Exception("This decorator takes a path argument")
        d = os.path.dirname(path)
        assert os.path.exists(d)
        p, q = os.path.splitext(path)
        i = [0]

        def decorator(func):
            def newfunc(*args, **kwargs):
                pr = cProfile.Profile()
                ret = pr.runcall(func, *args, **kwargs)
                if multi:
                    fn = "%s.%d%s" % (p, i[0], q)
                else:
                    fn = path
                i[0] += 1
                pr.dump_stats(fn)
                return ret
            # Be well-behaved
            newfunc.__name__ = func.__name__
            newfunc.__doc__ = func.__doc__
            newfunc.__dict__.update(func.__dict__)
            return newfunc
        return decorator

    def profile_simple(func):
        def wrapper(*args, **kwargs):
            datafn = func.__name__ + ".profile" # Name the data file sensibly
            prof = cProfile.Profile()
            retval = prof.runcall(func, *args, **kwargs)
            prof.dump_stats(datafn)
            return retval
        return wrapper

    def profileit(path, fn=None, stats=True):
        """
            example usage:  @profileit("profile_for_func1_001") 
        """
        fn = fn if fn is not None else '/tmp/profile'
        def inner(func):
            def wrapper(*args, **kwargs):
                prof = cProfile.Profile()
                retval = prof.runcall(func, *args, **kwargs)
                # Note use of name from outer scope
                if stats:
                    inn=open(path, 'a')
                    bac = sys.stdout
                    sys.stdout = inn
                    prof.dump_stats(fn)
                    
                    p = pstats.Stats(fn) 
                    p.sort_stats("cumulative").print_stats(60) 
                    inn.close()
                    sys.stdout = bac
                
                print('\n\nprofiling result is saved in %s'%(path, ))
                return retval
            return wrapper
        return inner



def func(x, y=2):
    a = 3
    res= a + x + y
    return res

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass 
    
    
    def test_timer_count(self): 
        pass 
    
        #@profileit(fn="profile.tmp")
        @count
        def func(x, y):
            return x**2+ y
        
        class A(object):
            @timer
            @count
            def __init__(self, x,z,  y=3):
                self.x = x + z
                print("x", x*y)

        for i in range(4):
            func(i, i)
            print(num_of_instance)

        A(5, 22);A(6, 2)

        #print tensor_player.STATE

    def test_temp(self): 
        pass
            
    

if __name__ == "__main__":
    if 0: 
        #@profileit(fn="profile.tmp")
        @count
        def func(x, y):
            return x**2+ y
        
        class A(object):
            @timer
            @count
            def __init__(self, x,z,  y=3):
                self.x = x + z
                print("x", x*y)

        for i in range(4):
            func(i, i)
            print(num_of_instance)

        A(5, 22);A(6, 2)

        #print tensor_player.STATE


    #warnings.filterwarnings('ignore')
    if 0:
        TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
        #'test_temp', 
        #'test_timer_count', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))

        unittest.TextTestRunner(verbosity=0).run(suite)

        

