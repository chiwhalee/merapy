#!/usr/bin/env python
#coding=utf8

import unittest

from merapy.graphics import TestIt as Test_graphics 
from merapy.run import TestIt as Test_run 
from merapy.utilities import TestIt as Test_util 
from merapy.context_util import TestIt as Test_context 
from merapy.tensor_svd import TestIt as Test_svd 
from merapy.common_util import TestCommon 
from merapy.tensor_network import TestIt as Test_tensor_network 
from merapy.top_level import TestIt as Test_top_level 

#suite = unittest.TestLoader().loadTestsFromTestCase(TestTensor)
#python -m unittest discover --pattern=*.py

if 1:  #skip all test_temp  
    all_cls= []
    dic=locals()    
    for k, v in dic.items(): 
        #if hasattr(v, '__base__') and hasattr(v, 'test_temp'): 
        if hasattr(v, '__base__'):  
            if v.__base__.__name__ ==  'TestCase': 
                all_cls.append(v)
                if hasattr(v, 'test_temp'): 
                    temp =  unittest.skip("skip test_temp")(getattr(v, 'test_temp')) 
                    setattr(v, 'test_temp', temp)

if 0: 
    unittest.main()
else: 
    from concurrencytest import ConcurrentTestSuite, fork_for_tests
    loader = unittest.TestLoader()
    suite = []
    for cls in all_cls: 
        temp = loader.loadTestsFromTestCase(cls)
        suite.append(temp)
    suite = unittest.TestSuite(suite)
    suite = ConcurrentTestSuite(suite, fork_for_tests(8))
    unittest.TextTestRunner(verbosity=0).run(suite)


#
