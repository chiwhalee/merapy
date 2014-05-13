#!/usr/bin/env python
#coding=utf8

"""
make compatable with binary graph
"""

from merapy.tensor_svd import Tensor_svd


def ascending_ham(M, S, ilayer, tau=None, info=0):
    """
        ilayer: ranges from 0 to M.num_of_layer-2
        calculate H_2 for ilayer+1, see Vidal 2007 fig.14
        this is the average ascending operator  \bar A =1/3(A_R + R_C + R_L) for ternary 
        act on S.H_2[ilayer].O[j], obtain S.H_2[ilayer+1].O[j].data
        注意这里与f90中不同，f90中是把 ilayer-1 ascending 到 ilayer
        
        tau: tau is used only when calculate coup coeff for long range models
    """
    #G_2_2 = init_mera_graph.G_2_2
    G_2_2 = S.G_2_2
    G_2_3 = S.G_2_3
    G_3_3 = S.G_3_3
    
    if ilayer<0 or ilayer>M.num_of_layer-2:
        raise ValueError("ilayer should be in [0, %d]. ilayer=%s"%(M.num_of_layer-2, ilayer))

    if tau is None:
        tau = ilayer
    
    if info>0:
        print "START ascending_ham"
    
    i=0
    j = 0
    ilayer_p1 = ilayer + 1
    
    name="OO"   
    S.H_2[ilayer_p1].O[j].data[:] =0.0
    for G in [G_2_2]:
        S.add_env_many(M, ilayer, name, G, add_tensor=S.H_2[ilayer_p1][0], branch=0, info=info-1)
    
    if not S.only_NN:
        name = "OOO"
        S.H_3[ilayer_p1].O[j].data[:] = 0.0
        for G in [G_2_3, G_3_3]:
            S.add_env_many(M, ilayer, name, G, add_tensor=S.H_3[ilayer_p1][0], branch=0, info=info-1)


    if info>0:
        print "H_2 at iLayer+1,",ilayer_p1
        print S.H_2[ilayer_p1].O[j].data[:4].round(3),"....", S.H_2[ilayer_p1].O[j].data[-4:].round(3) 
        print "END ascending_ham"

def descending_ham(M,S,ilayer,tau=None, info=0):
    """
        ilayer: ranges from M.num_of_layer-2, 1
        calculate rho_2 for ilayer-1
        this is the average desending operator  \bar D =1/3(D_R + R_C + R_L) for ternary 
        act on S.rho_2[ilayer].O[j], obtain S.H_2[ilayer-1].O[j].data
        注意这里与f90中不同，f90中是把 ilayer+1 desending 到 ilayer 
    """
    G_2_2 = S.G_2_2
    G_2_3 = S.G_2_3
    G_3_3 = S.G_3_3

    if info>0:
        print "START desending_ham"
    
    if tau is None:  
        tau = ilayer-1
    ilayer_m1 = ilayer-1
    
    name = "oo"  #note that rho_2 is OO, oo is H_2
    S.rho_2[ilayer_m1].O[0].data[:] = 0.0
    for G in [G_2_2, G_2_3]:
        S.add_env_many(M, ilayer_m1, name, G, add_tensor=S.rho_2[ilayer_m1][0], key_list=None, info=info-1)

    name = "ooo"
    S.rho_3[ilayer_m1].O[0].data[:] = 0.0
    for G in [G_3_3]:
        S.add_env_many(M, ilayer_m1, name, G, add_tensor=S.rho_3[ilayer_m1].O[0], key_list=None, info=info-1)

def iterative_optimize_all(M, S, layer, j=0, tau=None, info=0):
    """
        update U, V, H_2, rho_2 in layer layer
        首先计算这几个算子的evironment， 然后svd等处理
        layer: ranges from 1 to M.num_of_layer-2,  why not -1??
    """
    G_2_2 = S.G_2_2
    G_2_3 = S.G_2_3
    G_3_3 = S.G_3_3
    
    
    if info>0:
        print "START all_in_once"
    # Translational invariant         
    
    Names= ["V", "U"]
    #update V and U first, then V_dag and U_dag are just their conjugate resp.
    
    ilayer = layer
    #envs= TensorLink(size=8, use_buf=True) #contain envirenment tensors:   o, rho, h, w, u, u^+, w, w^+  
    envs= [None]*8
    num_of_tensor = 2
    
    envs[0] = M.V[ilayer][j].copy(use_buf=True)
    envs[1] = M.U[ilayer][j].copy(use_buf=True)

    for n  in range(2):
        envs[n].reverse_qsp()
    for n  in range(num_of_tensor):
        envs[n].data[:] = 0.0
    
    for G in [G_2_2, G_2_3, G_3_3]:
        weight = {g:-G[g].weight for g in G}
        for n  in range(2):
            name = Names[n]
            S.add_env_many(M, ilayer, name, add_tensor=envs[n], G_dic=G, weight_dic=weight, info=info-1)

                    

    #---------------------------------------------------------------------
    #---------------------------------------------------------------------    
    # Update Disentangler and isometry
    for n in range(num_of_tensor):
        if 1:
            envs[n].reverse_qsp()
        else:
            envs[n] = envs[n].dual(use_buf=True)
        name = Names[n]
        t = M.__getattribute__(name)
        ndiv = t[layer][0].rank//2 if name == "U" else t[layer][0].rank-1 
        Tensor_svd.svd(envs[n], ndiv=ndiv, nSd=0)   #this corresponds to vidal paper p.16 update w and u
        
        t[layer][j] = envs[n]
        
        t_dag = M.__getattribute__(name + "_dag")
        data = t_dag[layer][j].data.data
        t_dag[ilayer][j] = t[layer][j].conjugate(ndiv, buffer=data)

        if info>0:
            print "U at layer", ilayer, M.U[layer][j].data[range(5) + range(-5, 0)]
            print "END all_in_once"

