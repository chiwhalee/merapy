#coding=utf8
__all__ = ["ascending_ham"]

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
    
    if ilayer<0 or ilayer>M.num_of_layer-2:
        raise ValueError("ilayer should be in [0, %d]. ilayer=%s"%(M.num_of_layer-2, ilayer))

    if tau is None:
        tau = ilayer
    
    i=0
    j = 0
    name="OO"   #hamiltonian operator

    ilayer_p1 = ilayer + 1
    S.H_2[ilayer_p1].O[j].data[:] =0.0
    if info>0:
        print "START ascending_ham"
    for g in sorted(G_2_2):
        G = G_2_2[g]
        order = G.contract_order
        weight = G.weight
        S.add_env(M, ilayer, G, order, name, weight, add_tensor=S.H_2[ilayer_p1][0], info=info-1)
        
    
    if info>0:
        print "H_2 at iLayer+1,",ilayer_p1
        print S.H_2[ilayer_p1].O[j].data[:4].round(3),"....", S.H_2[ilayer_p1].O[j].data[-4:].round(3) 
        print "END ascending_ham"

def ascending_ham_1site(M, S, ilayer, tau=None, info=0):
    """
        ilayer: ranges from 0 to M.num_of_layer-2
        calculate H_2 for ilayer+1, see Vidal 2007 fig.14
        this is the average ascending operator  \bar A =1/3(A_R + R_C + R_L) for ternary 
        act on S.H_2[ilayer].O[j], obtain S.H_2[ilayer+1].O[j].data
        注意这里与f90中不同，f90中是把 ilayer-1 ascending 到 ilayer
        
        tau: tau is used only when calculate coup coeff for long range models
    """
    
    G_1_1 = init_mera_graph.G_1_1
    
    if ilayer<0 or ilayer>M.num_of_layer-2:
        raise ValueError("ilayer should be in [0, %d]. ilayer=%s"%(M.num_of_layer-2, ilayer))

    if tau is None:
        tau = ilayer
    
    i=0
    j = 0
    name="O"   #hamiltonian operator

    ilayer_p1 = ilayer + 1
    S.H_1[ilayer_p1].O[j].data[:] =0.0
    if info>0:
        print "START ascending_ham"
    for g in sorted(G_1_1):
        G = G_1_1[g]
    #for G in G_1_1.itervalues():
        order = G.contract_order
        weight= G.weight
        S.add_env(M, ilayer, G, order, name, weight, add_tensor=S.H_1[ilayer_p1][0], info=info-1)
        
    
    if info>0:
        print "H_1 at iLayer+1,",ilayer_p1
        print S.H_1[ilayer_p1].O[j].data[:4].round(3),"....", S.H_1[ilayer_p1].O[j].data[-4:].round(3) 
        print "END ascending_ham"

def ascending_generic(S): 
    """
        #not completed
    """
    pass


if __name__ == "__main__":

    from mera import test_Mera
    from hamiltonian import test_System
    if 0:
        symm = "Travial"
        M = test_Mera.instance(trunc_dim=2, tot_layer=5, symmetry=symm)
        sys = test_System.instance(M, model="Ising", symmetry=symm, only_NN=False)
    if 1:
        symm = "U1"
        M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry=symm)
        sys= test_System.instance(M, model="Heisenberg", symmetry=symm, only_NN=False, only_NNN=True)
    
    #set_ham_to_identity(sys)
    #from decorators import set_STATE_end_simple
    import decorators 
    for i in range(M.num_of_layer-1):
        #ilayer bellow 0 and >=M.num_of_layer-1 are not allowed
        decorators.set_STATE_end_simple(i, M.num_of_layer-1, iter0=0)
        ascending_ham(M, sys, ilayer=i, info=0)
    
    
    





