#coding=utf8

from matplotlib import pyplot as plt 
from merapy.measure_and_analysis.analysis import (RESULTDB_ROOT, 
        Analysis_mera , Analysis_idmrg, Analysis_vmps, Analysis_proj_qmc)

#heisbg_root = root = '/home/zhli/backup_tensor_dir/run-heisbg'  

heisbg_root = root = '/'.join([RESULTDB_ROOT , 'run-heisbg']) 


if 0: 
    an_mera = Analysis_mera()
    an_mera.alpha_parpath_dict[1.0]= heisbg_root + '/ternary'
    an_mera.alpha_parpath_dict['prod-state']= heisbg_root + '/ternary/prod-state'
    an_mera.alpha_parpath_dict['backup']= heisbg_root + '/ternary/backup'
    an_mera.alpha_list = [1.0, 'prod-state', 'backup']
    an = an_mera
if 1:     
    an_mera = Analysis_mera(local_root='/'.join([root, 'mera']), param_list=['Jzz']) 
    an_mera.set_parpath_dict()
    temp = [
            {'sub_dir': 'main'}, 
            {'sub_dir': 'eigen_state'}, 
            {'sub_dir': 'prod_state'}, 
            ]
    an_mera.add_sub_analysis(temp)
    an = an_mera

if 1: 
    an_vmps = Analysis_vmps(local_root='/'.join([root, 'vmps']), param_list=['Jzz']) 
    an_vmps.set_parpath_dict()
    #temp = [
    #        {'sub_dir': 'main'}, 
    #        {'sub_dir': 'main_symm'}, 
    #        {'sub_dir': 'test'}, 
    #        {'sub_dir': 'temp'}, 
    #        ]
    #an_vmps.add_sub_analysis(temp)
    an_vmps.add_sub_analysis('auto')
    an_vmps.an_test.add_sub_analysis('auto')
    

if 1: 
    an_idmrg = Analysis_idmrg()
    an_idmrg.alpha_parpath_dict[1.0] = '/'.join([root, 'idmrg'])
    an_idmrg.alpha_list = [1.0]
    an.ani = an_idmrg

 
    an_idmrg_psi= Analysis_idmrg(local_root='/'.join([root, 'idmrg-psi_guess']), param_list=['Jzz'])
    an_idmrg_psi.set_parpath_dict()
    #temp = [
    #        {'sub_dir': 'main'}, 
    #        {'sub_dir': 'main_symm'}, 
    #        {'sub_dir': 'main_fix_err_1em8'}, 
    #        {'sub_dir': 'performance'}, 
    #        {'sub_dir': 'test'}, 
    #        ]
    #an_idmrg_psi.add_sub_analysis(temp)
    an_idmrg_psi.add_sub_analysis('auto')
    an_idmrg_psi.an_main_nu_fix_err.param_list = ['Jzz', 'nu']
    an_idmrg_psi.an_main_filling_nu.param_list = ['Jzz', 'nu']
    
    #an_idmrg_psi.an_main_nu_fix_err.param_list = ['Jzz', 'nu']
    
    #an_idmrg_psi.add_sub_analysis([{'sub_dir':'test'}])
    an_idmrg_psi.an_test.add_sub_analysis('auto')
    an_idmrg_psi.an_test.an_xxz_v1v2.param_list = ['v1', 'v2', 'nu']
    
    an_idmrg_lam = Analysis_idmrg(local_root='/'.join([root, 'idmrg-lam_guess']), param_list=['Jzz'])
    an_idmrg_lam.set_parpath_dict()
    temp = [
            {'sub_dir': 'main'}, 
            {'sub_dir': 'main_symm'}, 
            {'sub_dir': 'performance'}, 
            {'sub_dir': 'test'}, 
            ]
    an_idmrg_lam.add_sub_analysis(temp)
    an_idmrg_lam.an_test.add_sub_analysis('auto')

if 1:
    an_qmc = Analysis_proj_qmc(local_root='/'.join([root, 'proj_qmc']), param_list=['Jzz']) 
    an_qmc.set_parpath_dict()




if __name__ == '__main__' : 
   
    xx = an_idmrg_lam 
    xx.outline_file()
    if 0: 
        #xx.measure_all(aa=[1.0], which=['entanglement'], force=0)
        xx.plot_entanglement_vs_chi( aa=[1.0, (1.0, 'a'), (1.0, 'ins_sites_new')], 
                marker='o', sh_min=(0, 6), sh_max=(0, 90))
        plt.show()
   
    
    

