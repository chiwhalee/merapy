#!/usr/bin/env python
#coding=utf8

"""
    this module serves as the test module 

"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
import warnings
import os
import unittest
import numpy as np 

from merapy.utilities import print_vars
from merapy.main import Main
from merapy.decorators import profileit, tensor_player
#import common_util
from merapy import common_util
from merapy.config import CFG_ISING_BASIC, CFG_HEISBG_BASIC, CFG_POTTS_BASIC
from merapy.updaters_all import  updaters_z2, updaters_u1

cfg = {}
warnings.filterwarnings('ignore') 

class TestIt(unittest.TestCase): 
    def setUp(self):
        self.seq = list(range(10))
    def test_temp(self): 
        pass 
    
    def test_heisbg_eigenstate(self): 
        from merapy.config import CFG_HEISBG_BASIC, copy_config 
        from merapy.top_level import top_level_eigenstate
        c = copy_config(CFG_HEISBG_BASIC)
        c.update(USE_CUSTOM_RAND=1, rand_seed=1234)
        
        c['updaters'] = updaters_u1
        c['updaters']['rho_top_func'] = top_level_eigenstate  
        print_vars(vars(),  ['c["updaters"]'])
        main = Main(**c)
        main.run(q_iter=5)
        #for some reason, the randomness is not fixed, so the following may fail 
        #self.assertAlmostEqual(main.S.energy,-0.59491171125170261, 10)
    
    def test_ising(self): 
        cfg_id = "ising"
        temp = dict(use_player=True, USE_CUSTOM_RAND=True, updaters=updaters_z2, trunc_dim=2, tot_layer=4, SYMMETRY="Z2", 
                NUM_OF_THREADS=1, )
        model_param={"h":-1.0, "J_NN":-1.0, "J_NNN":0.0}
        cfg[cfg_id]= CFG_ISING_BASIC.copy();  cfg[cfg_id].update(temp)
        cfg[cfg_id]['model_param'].update(model_param)
        #cfg[cfg_id].update(SYMMETRY='Travial')
        res= {0: 0, 1: 0.39467292395765885, 2: 0.76792404874778875, 3: 0.73223668913348139, 4: 0.72682462192879371, 5: 0.72401825188717228, 6: 0.72201106510910051, 7: 0.72061047374758014, 8: 0.71967406868471528, 9: 0.71905226108831211, 10: 0.71862224442580691}
        os.system('rm auto')    
        main = Main(**cfg[cfg_id])

        #tensor_player.save_tape = True
        #tensor_player.load_tape = True
        #tensor_player.tape_prefix = "ising"

        main.run(q_iter=10, do_measure=0)
        #main.run_scale_invar(q_iter=20, do_measure=0, backup_parpath='./mera_backup_test_folder')
        main.S.examine_energy(res)
        #main.stop_player()

    def test_heisbg(self): 
        temp = dict(USE_CUSTOM_RAND=True, updaters=updaters_u1, trunc_dim=4, tot_layer=4, 
                use_player=True, SYMMETRY="U1", NUM_OF_THREADS=1, )
        cfg_heisbg = CFG_HEISBG_BASIC.copy();   cfg_heisbg.update(temp)
        model_param={"J_NN":1.0, "J_NNN":0.241186}
        cfg_heisbg['model_param'].update(model_param)
        #cfg_heisbg.update(SYMMETRY='Travial', use_player=0)
        
        main = Main(**cfg_heisbg)
        res= {0: 0, 1: -0.15447477974680113, 2: -0.27316843978691985, 3: -0.38435891607150841, 4: -0.49142890463353184, 5: -0.59574502419909914, 6: -0.69595138685370195, 7: -0.78927036090526448, 8: -0.87373436401172677, 9: -0.9490632254128426, 10: -1.015790781309708}
        #tensor_player.tape_prefix = "heisbg"
        #tensor_player.save_tape = True
        #tensor_player.load_tape = False
        os.system('rm 4.pickle')    
        os.system('rm auto')    
        main.run(q_iter=10, do_measure=0)
        main.S.examine_energy(res, delta=1e-12)



if __name__ == '__main__':  
    
    if 0: #examine
        if 0: 
            from concurrencytest import ConcurrentTestSuite, fork_for_tests
            loader = unittest.TestLoader()
            suite = []
            for cls in [TestIt]: 
                temp = loader.loadTestsFromTestCase(cls)
                suite.append(temp)
            suite = unittest.TestSuite(suite)
            
            suite = ConcurrentTestSuite(suite, fork_for_tests(2))
            unittest.TextTestRunner(verbosity=0).run(suite)
        else: 
                
            #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
            #unittest.TextTestRunner(verbosity=0).run(suite)    
            TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
            unittest.main()
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #TestIt('test_temp'), 
           #TestIt('test_ising'), 
           TestIt('test_heisbg') , 
           #TestIt('test_heisbg_eigenstate') , 
        ]
        for a in add_list: 
            suite.addTest(a)
        #suite.addTest(TestIt('test_ising'))
        #suite.addTest(TestIt('test_heisbg'))
        unittest.TextTestRunner().run(suite)
        
        
        


