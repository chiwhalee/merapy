#coding=utf8
"""
    schedules for Main.run_schedule

"""
from collections import OrderedDict

#schedule = OrderedDict()
if 1: 
    schedule = [
        {'ansatz':'prod_state',  'shape': (4, 4),  'q_iter':10,  'energy_diff_min': 1e-10},  
        {'ansatz':'scale_invar', 'shape': (4, 4),  'q_iter':10,  'energy_diff_min': 1e-10}, 
        {'ansatz':'scale_invar', 'shape': (8, 4),  'q_iter':10,  'energy_diff_min': 1e-10}, 
        {'ansatz':'scale_invar', 'shape': (12,4),  'q_iter':10,  'energy_diff_min': 1e-10}, 
        ]

    
    schedule_prod_state = [
        {'ansatz':'prod_state',  'shape': (4, 2), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        {'ansatz':'prod_state',  'shape': (8, 2), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        {'ansatz':'prod_state',  'shape': (12,2), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        
        {'ansatz':'prod_state',  'shape': (4,3), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        {'ansatz':'prod_state',  'shape': (8,3), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        {'ansatz':'prod_state',  'shape': (12,3), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        
        {'ansatz':'prod_state',  'shape': (4,4), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        #{'ansatz':'prod_state',  'shape': (4,5), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        {'ansatz':'prod_state',  'shape': (4,6), 'q_iter':50000,  'energy_diff_min': 1e-10},  
        
        {'ansatz':'prod_state', 'shape': (8, 4),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        {'ansatz':'prod_state', 'shape': (8, 5),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        {'ansatz':'prod_state', 'shape': (8, 6),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        {'ansatz':'prod_state', 'shape': (8, 7),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        
        {'ansatz':'prod_state', 'shape': (12,4),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (12,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (12,6),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (12,7),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (12,8),  'q_iter':10000,  'energy_diff_min': 1e-9}, 

        
        {'ansatz':'prod_state', 'shape': (14,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (14,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        {'ansatz':'prod_state', 'shape': (16,4),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (16,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (16,6),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        {'ansatz':'prod_state', 'shape': (17,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (17,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 

        {'ansatz':'prod_state', 'shape': (20,4),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (20,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (20,6),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
 
        {'ansatz':'prod_state', 'shape': (21,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (21,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        {'ansatz':'prod_state', 'shape': (27,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'prod_state', 'shape': (27,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
      
        ]

    schedule_scale_invar = [
        {'ansatz':'prod_state',   'shape': (4, 4),  'q_iter':1000,  'energy_diff_min': 1e-10},  
        {'ansatz':'scale_invar', 'shape': (4, 4),  'q_iter':50000,  'energy_diff_min': 1e-10},  
        
        {'ansatz':'scale_invar', 'shape': (8, 4),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        {'ansatz':'scale_invar', 'shape': (8, 5),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
                                                                 
        {'ansatz':'scale_invar', 'shape': (12,4),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,5),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,6),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,7),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,8),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        
        # note (14, 5)=(6, 3, 3, 1, 1)
        {'ansatz':'scale_invar', 'shape': (14,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (14,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 

        # note (16, 5)=(6, 4, 4, 1, 1)
        {'ansatz':'scale_invar', 'shape': (16,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        {'ansatz':'scale_invar', 'shape': (16,4,3),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (16,5,3),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (16,6,3),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
       
        # note (17, 5)=(7, 4, 4, 1, 1)
        {'ansatz':'scale_invar', 'shape': (17,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (17,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        # note (18, 5)=(8, 4, 4, 1, 1)
        {'ansatz':'scale_invar', 'shape': (18,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (18,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 

        # note (21, 5)=(9, 5, 5, 1, 1)
        {'ansatz':'scale_invar', 'shape': (21,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (21,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        
        {'ansatz':'scale_invar', 'shape': (20,4),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (20,5),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (20,6),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
       
        {'ansatz':'scale_invar', 'shape': (24,4),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (24,5),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
       
        {'ansatz':'scale_invar', 'shape': (27,4,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (27,5,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
    
        #{'ansatz':'scale_invar', 'shape': (30,4),    'q_iter':10000,  'energy_diff_min': 1e-9}, 
        #{'ansatz':'scale_invar', 'shape': (30,5),    'q_iter':10000,  'energy_diff_min': 1e-9}, 

        ]
    
    
    schedule_scale_invar_3qn = [
        {'ansatz':'prod_state',   'shape': (4, 4),  'q_iter':1000,  'energy_diff_min': 1e-10},  
        {'ansatz':'scale_invar', 'shape': (4, 4),  'q_iter':50000,  'energy_diff_min': 1e-10},  
        
        {'ansatz':'scale_invar', 'shape': (8, 4),  'q_iter':50000,  'energy_diff_min': 1e-10}, 
        {'ansatz':'scale_invar', 'shape': (8, 5),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
                                                                 
        {'ansatz':'scale_invar', 'shape': (12,4),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,5),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,6),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,7),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (12,8),  'q_iter':50000,  'energy_diff_min': 1e-9}, 
                                                                     
        {'ansatz':'scale_invar', 'shape': (16,4),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
        {'ansatz':'scale_invar', 'shape': (16,5),  'q_iter':10000,  'energy_diff_min': 1e-9}, 
   
        ]
    
   
    
def copy_schedule(schedule, skip_list=[]): 
    #res= [i for i in schedule]
    res= []
    if len(skip_list)>0: 
        msg = 'skipped shapes in schedule are %s'%(skip_list, )
        print(msg)
    for s in schedule:
        if not s['shape'] in skip_list: 
            res.append(dict(s))
     
    return res


#---------------------------------------------------
def test_copy_schedule(): 
    s=copy_schedule( schedule_scale_invar, skip_list=[(4, 4)])
    print s
    

if __name__ == '__main__' : 
    pass
        
    
