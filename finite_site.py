#coding=utf8
"""

"""
from merapy.minimize import finite_site, finite_site_u1


__all__ = ["finite_site", "finite_site_u1"]



if __name__ == "__main__":

    import profile
    from merapy.mera import test_Mera
    from merapy.all_in_once import update_all_in_once
    from current.iteration import iterative_optimize_all
    from current.ascending import ascending_ham
    from current.descending import descending_ham


    def test(trunc_dim=2, tot_layer=4, symmetry="Z2", model="Ising", only_NN=True):
        q_iter=1
        q_one=1
        q_lay=1
        LayStart=0
        M = test_Mera.instance(trunc_dim, tot_layer, symmetry=symmetry)
        sys= test_System.instance(M, symmetry= symmetry, model=model, only_NN=only_NN)


        SystemParam.h=-1.0
        SystemParam.J_NN = -1.0
        sys.init_Hamiltonian(M)
        for iLayer in range(M.num_of_layer-1):
            M.init_layer(iLayer, True)
        #print M.__repr__(which=["U_dag", "V_dag"])        
        
        
        #finite_site(M,sys, q_one=q_one, q_lay=q_lay, q_iter=q_iter, lstart = 0, lstep=0, cont=False, info=0)
        finite_site(M,sys, ascending=ascending_ham, descending=descending_ham, update_mera=iterative_optimize_all, q_one=q_one, q_lay=q_lay, q_iter=q_iter, lstart = LayStart, lstep=0, cont=SystemParam.ReAssumed, info=0)
        M.show_param()
    
    test(trunc_dim=4, tot_layer=4, symmetry="U1", model="Heisenberg", only_NN=True)
    #test(trunc_dim=2, tot_layer=4, symmetry="Z2", model="Ising")




