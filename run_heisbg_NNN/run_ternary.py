#!/usr/bin/env python
#coding=UTF8
#from __future__ import absolute_import
import numpy as np 

from merapy.main import Main 
from merapy.run_heisbg_NNN.config import make_config



J2J2 = [2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1]
J2J2 = np.arange(0, 1.1, 0.1).round(3).tolist()

J2J2 += [0.241186] 

root1 = ['temp']
surfix = [''] 
grid = Main.make_grid(surfix, J2J2, root1)


config_group = []
for i in grid: 
    sur, J2, r = i
    c = make_config(J2, root1=r, surfix=sur) 
    c['schedule']['dim_diff_remap'] = {4:5e-9, 8:5e-8, 12:8e-8, 14:9e-8, 17:1e-7}   
    config_group.append(c)
    print c['backup_parpath']
    print c['schedule']['dim_diff_remap']
    
num_threads = 1
for c in config_group: 
    c['schedule'].update(mera_shape_min=(4, 4), mera_shape_max=(12, 4), skip_list=[(8, 5), (8, 10)])
    c['NUM_OF_THREADS'] = num_threads
   
    if 0:  # for testing  
        c['only_NN'] = 1
        c['backup_parpath'] = mkdtemp()
        c['schedule'].update(q_iter_max=5, mera_shape_max=(4, 4), do_measure=0)
        print  c['backup_parpath']

Main.run_many(config_group, nproc=12, parallel=1)


    
if __name__=="__main__":          

    pass



