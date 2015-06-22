#!/usr/bin/env python
#coding=utf8
import numpy as np
from tempfile import mkdtemp
from merapy.main import Main
from merapy.context_util import make_temp_dir

from config import make_config


hh = np.arange(0.8, 1.2, 0.02).round(3).tolist()

#root1 = ['test/eigen_state']
root1 = ['main']
surfix = [''] 
grid = Main.make_grid(surfix, hh, root1)


config_group = []
for i in grid: 
    sur, h, r = i
    c = make_config(h, root1=r, surfix=sur, which_top_state='scale_invar_state') 
    c['schedule']['dim_diff_remap'] = {4:1e-9, 8:3e-8, 12:8e-8, 14:9e-8, 17:1e-7}   
    c.update(register_job=1, job_group_name='mera1')
    config_group.append(c)
    print c['backup_parpath']

for c in config_group: 
    c['NUM_OF_THREADS'] = 1
    c['schedule'].update(mera_shape_min=(4, 4), 
            mera_shape_max=(8, 4), skip_list=[(8, 5), (8, 10)])

#Main.run_many(config_group, nproc=None, parallel=1)


