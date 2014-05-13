#encoding=utf8

import warnings

from merapy.system_parameter import SystemParam
#from merapy.tensor import *
from merapy.tensor_network import *
from merapy.mera import *
from merapy.hamiltonian import *

#from init_mera_graph import *  #init_mera_graph, init_g2
#import init_mera_graph


def descending_ham(M,S,ilayer,tau=None, info=0):
    """
        ilayer: ranges from M.num_of_layer-2, 1
        calculate rho_2 for ilayer-1
        this is the average desending operator  \bar D =1/3(D_R + R_C + R_L) for ternary 
        act on S.rho_2[ilayer].O[j], obtain S.H_2[ilayer-1].O[j].data
        注意这里与f90中不同，f90中是把 ilayer+1 desending 到 ilayer 
    """
    G_2_2 = S.G_2_2

    if info>0:
        print "START desending_ham"
    
    if tau is None:  
        tau = ilayer-1
    ilayer_m1 = ilayer-1
    
    if 1:
        name = "oo"  #note that rho_2 is OO, oo is H_2
        S.rho_2[ilayer_m1].O[0].data[:] = 0.0
        temp = S.add_env_many(M, ilayer_m1, name, G_2_2, key_list=None, info=info-1)
        S.rho_2[ilayer_m1].O[0].data += temp.data


    if not S.only_NN:
        G_3_2 = S.G_3_2
        G_3_3 = S.G_3_3

        name = "ooo"
        S.rho_3[ilayer_m1].O[0].data[:] = 0.0
        for G in [G_3_2, G_3_3]:
            #temp = S.add_env_many(M, ilayer_m1, name, G, key_list=None, info=info-1)
            #S.rho_3[ilayer_m1].O[0].data += temp.data
            S.add_env_many(M, ilayer_m1, name, G, add_tensor=S.rho_3[ilayer_m1].O[0], key_list=None, info=info-1)



if __name__ == "__main__": 
    from merapy.mera import test_Mera
    if 1:
        M = test_Mera.instance(trunc_dim=2, tot_layer=4, symmetry="Z2")
        sys = test_System.instance(M, model="Ising", symmetry="Z2", only_NN=False)
    if 0:
        pass
        M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1")
        sys= test_System.instance(M, model="Heisenberg", symmetry="U1", only_NN=False)



    def reset_mera(M):
        import numpy as np
        for i in range(M.num_of_layer-1):
            t = M.V[i].tensor[0]
            rank = t.rank
            qDims= np.empty(rank, "int")
            iDims= np.empty(rank, "int")
            qDims[:]= [0, 0, 0, 0]
            iDims[:] = 0
            t.data[:] = 0.0
            t.set_element(qDims, iDims, 1.0)


    from merapy.top_level import *
    top_level_product_state_u1(M, sys)
    print sys
    for i in range(M.num_of_layer-1, 0, -1):
        print "iii", i
        descending_ham(M, sys, ilayer=i, info=5)
        exit()
        #print sys.H_2[i+1].O[0].data.round(5)
    #print M.V[M.num_of_layer-1].tensor[0]
    #print sys
    #print sys.H_2[1].O[0].data
    





