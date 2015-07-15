#coding=utf8
"""
q:  
    q69 why no   +1?
    q41  why revers qn here?,  because it is env which just means is conjugate? 
    q59 why *0.5
    

"""

from tensor import *
from tensor_svd import *
from tensor_network import *
from mera import *
from hamiltonian import *
from tensor_reflection import TensorReflect

#from init_mera_graph import G_2_2
from merapy.diagrams.V31.graph_ternary import G_2_2
#G3, order3, weight3, GTop, order_top  = init_mera_graph()
#G2, weight2, order2= init_g2()

__all__ = ["update_all_in_once"]


def update_all_in_once(M,S,layer,j=0, tau=None, info=0):
    """
    update U, V, H_2, rho_2 in layer layer
    首先计算这几个算子的evironment， 然后svd等处理
    tau: not used at all in f90; #tau = ilayer-1
    layer: ranges from 1 to M.num_of_layer-2,  why not -1??
    顶层的V tensor的更新方法不同，因为对应于不同的图
    """
    
    if info>0:
        print "START all_in_once"
    # Translational invariant         
    jEnv = 0
    iEnv = 0
    Names= ["V", "U", "OO", "oo"]
    #update V and U first, then V_dag and U_dag are just their conjugate resp.
    
    ilayer = layer
    
    #Ops= TensorLink(size=2, use_buf=True)  #not used in f90
    #envs= TensorLink(size=8, use_buf=True) #contain envirenment tensors:   o, rho, h, w, u, u^+, w, w^+  
    envs= [None]*8
    
    num_of_tensor = 4
    ndiv = [3, 2, 2, 2]
    #ndiv = [2, 1, 2, 2]
    
    
    envs[0] = M.V[ilayer][j].copy(use_buf=True)
    envs[1] = M.U[ilayer][j].copy(use_buf=True)
    envs[2] = S.H_2[ilayer+1][j].copy(use_buf=True)   # note above Names[2]="OO"
    envs[3] = S.rho_2[ilayer][j].copy(use_buf=True)

    #q41
    for n  in range(2):
        envs[n].reverse_qsp()
    
    for n  in range(num_of_tensor):
        envs[n].data[:] = 0.0

    #for g in [-1, 0]:
    for g in [(1, 2), (2, 3)]: 
        graph = G_2_2[g]
        order = graph.contract_order
        weight = -graph.weight
        #if g == -1:  weight = weight*0.5
        if g == (1, 2):  weight = weight*0.5
        restart =  True 
        for n  in range(num_of_tensor):
            name = Names[n]
            temp=S.add_env(M,ilayer, graph, order, name, weight, restart=restart, info=info-1)
            envs[n].data += temp.data
            #here lonly change data,  but not QSp and etc.
            restart =  False 
    
    if S.USE_REFLECTION:   #note this is defferent with f90 
        for n in range(2, num_of_tensor):
            #print "iii", ilayer, n
            #print envs[n].data 
            TensorReflect.reflect_symmetrize(envs[n], ndiv=ndiv[n])
            #print envs[n].data,"\n" 
    
    
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------    
    #update H_2 and rho_2
    #S.H_2[ilayer+1].O[j]=envs[2].scale(-1.0)
    #S.H_2[ilayer+1].O[j].data = -1.0*envs[2].data
    S.H_2[ilayer+1][j].data[:] = -1.0*envs[2].data[:]
    
    X = envs[3].trace()
#    print 'X=', X
    if abs(X) >= 1.E-8:  
        X = 1.0/X
    else:
        X = -1.0
    
    #print '\nTrace[rho_2]=', X
    
    #S.rho_2[ilayer-1].O[j]=envs[3].scale(X)
    #S.rho_2[ilayer].O[j]=envs[3].scale(X)
    S.rho_2[ilayer][j].data[:] = X*envs[3].data[:]
    
    
    
    # Update Disentangler and isometry
    
    for n  in range(2):
        envs[n].reverse_qsp()
            #q: why not reverse qsp?
    
    if S.USE_REFLECTION:    
        for n  in range(2):
            TensorReflect.reflect_symmetrize(envs[n], ndiv=ndiv[n])
    
    for n  in range(2):
        #this corresponds to vidal paper p.16 update w and u
        nsd = 0
        #name = Names[n]        
        #envs[n] = Tensor_svd(envs[n])
        #Sd1=envs[n].svd(ndiv=ndiv[n], nSd=nsd)
        Tensor_svd.svd(envs[n],ndiv=ndiv[n], nSd=nsd)

    
    n = 0  #V tensor
    M.V[layer][j] = envs[n]
    data = M.V_dag[ilayer][j].data.data
    M.V_dag[ilayer][j]=M.V[ilayer][j].conjugate(ndiv[n], buffer=data)
    
    n = 1  #U  tensor
    M.U[layer][j] = envs[n]
    data = M.U_dag[ilayer][j].data.data
    M.U_dag[ilayer][j]=M.U[ilayer][j].conjugate(ndiv[n], buffer=data)

    if info>0:
        print "U at layer", ilayer, M.U[layer][j].data[range(5) + range(-5, 0)]
        print "END all_in_once"


if __name__ == "__main__":
    """  ---pass  """
    #M = test_Mera.instance()    
    #sys=  test_System.instance(M)
     
    M = test_Mera.instance(trunc_dim=2, symmetry="U1")
    sys= test_System.instance(M, model="Heisenberg")


    update_all_in_once(M, sys,layer=0, info=2)
    
    #print "\nafter all_in_once"
    print M.__repr__(layers=[0], which=["U","V","U_dag","V_dag"])
    #print sys.__repr__(layers=[0,1], which=["H_2","rho_2"])


    
      


