#coding=utf8

from merapy.tensor_svd import *

__all__ = ["iterative_optimize_all"]

       
def iterative_optimize_all(M, S, layer, j=0, tau=None, info=0):
    """
        update U, V, H_2, rho_2 in layer layer
        首先计算这几个算子的evironment， 然后svd等处理
        tau: not used at all in f90; #tau = ilayer-1
        layer: ranges from 1 to M.num_of_layer-2,  why not -1??
        top V is not updated here
    """
    G_2_2 = S.G_2_2
    
    if info>0:
        print "START iterative_optimize"
    # Translational invariant         

    ilayer = layer
    if tau is None:
        tau = ilayer

    Names = ["V", "U"]   #update V and U first, then V_dag and U_dag are just their conjugate resp.
    
    num_of_tensor = len(Names)
    envs = [None]*num_of_tensor     #contain envirenment tensors:   o, rho, h, w, u, u^+, w, w^+  

    for n in range(num_of_tensor):
        name = Names[n]
        if 1:
            #try:
            #    envs[n] = M.__getattribute__(name)[ilayer][j].copy(use_buf=True)
            #except:
            #    msg = """ %s\n%s\n%s
            #    """%(name, ilayer, M.__getattribute__(name)[ilayer][j])
            #    raise Exception(msg)
            envs[n] = M.__getattribute__(name)[ilayer][j].copy(use_buf=True)
            envs[n].reverse_qsp()
        else:
            envs[n] = M.__getattribute__(name)[ilayer][j].dual(use_buf=True)
        envs[n].data[:] = 0.0
    
    for g in sorted(G_2_2):
        graph  = G_2_2[g]
        #for G in G_2_2.itervalues():

        order= graph.contract_order
        weight = - graph.weight  # in fact, weight does not mater here, since we will do svd on envs[n] later
        restart =  True 
        for n  in range(num_of_tensor):
            name = Names[n]
            S.add_env(M, ilayer, graph, order, name, weight, add_tensor=envs[n], restart=restart, info=info-1)
            restart =  False 
    
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

if __name__  == "__main__":
    pass
    from mera import test_Mera
    from hamiltonian import test_System
    M = test_Mera.instance(trunc_dim=4, symmetry="U1")
    sys= test_System.instance(M, "U1", model="Heisenberg")


    iterative_optimize_all(M, sys,layer=0, info=2)
    
    #print "\nafter all_in_once"
    print M.__repr__(layers=[0], which=["U","V","U_dag","V_dag"])
    #print sys.__repr__(layers=[0,1], which=["H_2","rho_2"])



