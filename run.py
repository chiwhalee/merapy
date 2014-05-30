#!/usr/bin/env python
#coding=utf8

"""
    this module serves as the test module 

"""

import warnings
import os
import unittest

from merapy.main import Main
from merapy.decorators import profileit, tensor_player
#import common_util
from merapy import common_util
from merapy.config import CFG_ISING_BASIC, CFG_HEISBG_BASIC, CFG_POTTS_BASIC ,  updaters_z2, updaters_u1

cfg = {}
warnings.filterwarnings('ignore') 

class TestIt(unittest.TestCase): 
    def setUp(self):
        self.seq = range(10)
    def test_temp(self): 
        pass
    def test_ising(self): 
        cfg_id = "ising"
        temp = dict(use_player=True, USE_CUSTOM_RAND=True, updaters=updaters_z2, trunc_dim=2, tot_layer=4, SYMMETRY="Z2", 
                NUM_OF_THREADS=1, )
        model_param={"h":-1.0, "J_NN":-1.0, "J_NNN":0.0}
        cfg[cfg_id]= CFG_ISING_BASIC.copy();  cfg[cfg_id].update(temp)
        cfg[cfg_id]['model_param'].update(model_param)
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
        
        main = Main(**cfg_heisbg)
        res= {0: 0, 1: -0.15447477974680113, 2: -0.27316843978691985, 3: -0.38435891607150841, 4: -0.49142890463353184, 5: -0.59574502419909914, 6: -0.69595138685370195, 7: -0.78927036090526448, 8: -0.87373436401172677, 9: -0.9490632254128426, 10: -1.015790781309708}
        #tensor_player.tape_prefix = "heisbg"
        #tensor_player.save_tape = True
        #tensor_player.load_tape = False
        os.system('rm 4.pickle')    
        os.system('rm auto')    
        main.run(q_iter=10, do_measure=0)
        main.S.examine_energy(res)



if __name__ == '__main__':  
    
    if 1: #examine
        if 1: 
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
           TestIt('test_ising'), 
           TestIt('test_heisbg') 
        ]
        for a in add_list: 
            suite.addTest(a)
        #suite.addTest(TestIt('test_ising'))
        #suite.addTest(TestIt('test_heisbg'))
        unittest.TextTestRunner().run(suite)
       

#else:
    if 0:  #heisbg test
        cfg_id = 0
        temp = dict(updaters=updaters_u1, trunc_dim=4, tot_layer=4, use_player=True, filename=None, backup_fn=None, 
                SYMMETRY="Travial", NUM_OF_THREADS=0, USE_CUSTOM_RAND=True, 
                model_param={"J_NN":1.0, "J_NNN":0.241186}, info=0)
        cfg[cfg_id] = CFG_HEISBG_BASIC.copy(); cfg[cfg_id].update(temp)
        res = {0: 0, 1: -0.15447477974680091, 2: -0.27316843978691718, 3: -0.38435891607150996, 4: -0.49142890463353295, 5: -0.59574502419910091}
        
        backup_fn = None #"tmp.pickle"
        backup_fn1 = None
        main= Main(**cfg[cfg_id])
        main.run(q_iter=5, backup_fn=backup_fn)
        main.run_scale_invar(q_iter=100, backup_fn=backup_fn)
        main.expand_layer(5)
        main.run()
        exit()
        if 0:
            import os
            os.environ["OMP_SCHEDULE"] = "dynamic"#"static"#"dynamic" 
            main= Main(**cfg[cfg_id])
            tensor_player.save_tape= True
            tensor_player.load_tape = True
            tensor_player.tape_prefix = ""
            main.run(q_iter=10)
        
        
        main.stop_player()
        if 0:
            pass
            #from common_64_ifort import contract_core_player_fort_paralell_test
            import common_64_ifort as c64
            common_util.contract_core_player_fort = c64.contract_core_player_fort_paralell_critical
            main.run(q_iter=10)
            print "reduction "*5
            common_util.contract_core_player_fort = c64.contract_core_player_fort_paralell_reduction_1
            main.run(q_iter=10)
     
    if 0:  #ising test
        
        temp = dict(updaters=updaters_u1, trunc_dim=4, tot_layer=4, use_player=True, filename=None, backup_fn=None, 
                SYMMETRY="U1", NUM_OF_THREADS=1, USE_CUSTOM_RAND=True, info=0)
        cfg = CFG_ISING_BASIC.copy(); cfg.update(temp)
        model_param={"h":-1.0, "J_NN":-1.0, "J_NNN":0.0} 
        cfg['model_param'].update(model_param)

        cfg.update({'use_player': 1, 'trunc_dim':2, 'SYMMETRY':'Travial', "combine_2site":False, "only_NN":True})
        #cfg[0].update({"combine_2site":False, "only_NN":True})
        
        from merapy.top_level import top_level_eigenstate
        cfg['updaters'].update(rho_top_func=top_level_eigenstate)
       
        main = Main(**cfg)
        
        #main.run_scale_invar(q_iter=0, backup_fn='auto', backup_parpath='./mera_backup_test_folder/', do_measure=0)
        
        import config, schedule
        backup_parpath = config.MERA_BACKUP_TENSOR_DIR  + '/mera_backup_test_folder' 
        os.system('rm /home/zhli/Documents/mera_backup_tensor/mera_backup_test_folder/*' )    
        
        if 1: 
            main.run_schedule(schedule=schedule.schedule_scale_invar, do_measure=0, backup_fn='auto', 
                    backup_parpath=backup_parpath, 
                    q_iter_max = 10, 
                    mera_shape_min = (4, 4), mera_shape_max=(4, 4), info=0, 
                    #dim_diff_remap = {4:1e-1, 8:1e-4}
                    )
            #main.run(q_iter=10, backup_parpath='./mera_backup_test_folder', do_measure=0)
            #main.run_scale_invar(q_iter=20, backup_parpath='./mera_backup_test_folder', do_measure=0)
        

