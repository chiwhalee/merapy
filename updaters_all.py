#!/usr/bin/env python
#coding=utf8


from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
from merapy.ascending import ascending_ham
from merapy.descending import descending_ham
from merapy.iteration import iterative_optimize_all
from merapy.all_in_once import update_all_in_once
from merapy.top_level import top_level_product_state, top_level_product_state_u1, top_level_eigenstate
from merapy.finite_site import finite_site, finite_site_u1
#import merapy.finite_site as finite_site_module
import merapy.updaters_binary as upbin
import merapy.schedule as schedule_module


updaters_u1 = {
        "ascending_func":ascending_ham, 
        "descending_func":descending_ham, 
        "update_mera_func":iterative_optimize_all, 
        "finite_range_func":finite_site_u1, 
        "rho_top_func":top_level_product_state_u1}

updaters_binary = updaters_u1.copy()
updaters_binary.update({
    "ascending_func":upbin.ascending_ham, 
    "descending_func":upbin.descending_ham, 
    "update_mera_func":upbin.iterative_optimize_all})


#updaters_z2= [
#        ascending_ham, 
#        None, 
#        update_all_in_once,
#        finite_site, 
#        top_level_product_state]

updaters_z2 = {
        "ascending_func":ascending_ham, 
        "descending_func":None, 
        "update_mera_func":update_all_in_once, 
        "finite_range_func":finite_site, 
        "rho_top_func":top_level_product_state}



