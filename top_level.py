#coding=utf8

import numpy as np
import unittest 

from mera import *
from hamiltonian import *
from tensor_svd import Tensor_svd
#from init_mera_graph import init_gtop
from merapy.diagrams.V31.graph_ternary import GTop, order_top
import warnings

from merapy.utilities import print_vars

"""
q:
    q82,  why everage

"""

__all__ = ["top_level_eigenstate", "top_level_product_state", "top_level_product_state_u1", 
        "top_level_product_state_u1_1site"]

#GTop, order_top = init_gtop()

if 0:
    def top_level_ham2_Ternary(M,S, topQN):
        iLayer = M.num_of_layer-2
        j=0


        ord1=[1,0,3,2]
        H1 = S.H_2[iLayer].O[j].permutation(ord1, use_buf=True)
        Tensor_Add_WithScale(H1, 1.0, S.H_2[iLayer].O[j],1.0, H2)
        H2 = H1.data[:]
        H1 = H2
        
        tQN = topQN    
        iTensor_Eig(H1, V1, E1, tQN)
        
        tQN.val = -tQN.val
        iTensor_Eig(H2, V2, E2, tQN)
        
    #    print 'E1,E2=', E1,E2
        if E2 <= E0:  
            print 'E2=', E2
            M.V[iLayer+1][j] = V2
    #        print_iTensor(M.V[iLayer+1][j], True )
            topQN = tQN
        else:
            M.V[iLayer+1][j] = V1
    #    iTensor_Eig(H2, M.V[iLayer+1][j], Eg, topQN)
        
        Tensor_Conjugate(M.V[iLayer+1][j].rank-1, M.V[iLayer+1][j], M.V_dag[iLayer+1][j])
        
    def top_level_eigenstate_Ternary(M,S):
        """
        see Ternary/TopLevel.f90
        """
        topQN = M.topQN.copy()
        top_level_ham2(M,S, topQN)
        
        G = GTop
        order = order_top
        iLayer = M.num_of_layer-1
        j=0
        
        Names = "OO"
        #init_TLink(ops, 2, use_buf= True )
        ops= [None]*2
        
        if topQN.val == 0:  
            #Tensor_ShallowCopy(rho_top, ops[1])
            ops[0] = rho_top
        else:
            Tensor_ShallowCopy(rho_top_odd, ops[1])
        
        System_Contraction(M, S, iLayer, G, 1, Names, ops, S.rho_2[iLayer-1].O[j], order, exception=1, restart= True )
        
    #    if topQN.val == -0:  
    #    print_iTensor(S.rho_2[iLayer-1].O[j], True )
        delete_TLink(ops)


def top_level_ham2(M, S, topQN):
    """
        将H_2[M.num_of_layer-2] 对角化，得到M.V[M.num_of_layer-1]
    
    """
    iLayer = M.num_of_layer-1
    j=0
    ord1=[1,0,3,2]
    H1 = S.H_2[iLayer][j].permutation(ord1, use_buf=True)
    
    #average 
    #q82
    H2 = H1 + S.H_2[iLayer][j]
    #H2 should have divide 2, but it doesn't mater, since we are doing eig on it
    tQN = topQN.copy()

    ilayer_p1 = iLayer #+ 1
    
    M.V[ilayer_p1][j], E = Tensor_svd.eig(H2, tQN)
    M.V_dag[ilayer_p1][j] = M.V[ilayer_p1][j].conjugate(M.V[ilayer_p1][j].rank-1, use_buf=True)


def top_level_eigenstate(M, S):
    """
        it is in fact a desending operation: map from rho_top of 
        rank (1, 1) to rho[M.num_of_layer-1] of rank (2, 2)
        todo: replace this function by descending 
    """
    #print_vars(vars(),  ["repr('use top_level_eigenstate')"])
    topQN = M.topQN.copy()
   
    top_level_ham2(M, S, topQN)
    
    G = GTop
    order = order_top
    iLayer = M.num_of_layer-1
    j=0
    
    name = "oo"
    weight = 0.5
    S.rho_2[iLayer][j].data[:] = 0.0
    
    
    S.rho_2[iLayer][j]=S.add_env(M, iLayer, G, order, name, weight, info=0)
    
    #average
    order1 = [1, 0, 3, 2]
    rho = S.rho_2[iLayer][j].permutation(order1, use_buf=True)
    S.rho_2[iLayer][j] += rho



"""
top_level_product_state 有3个用处，
1. 在最上层与ascending上去的H收缩求能量, 2.更新上面的V tensor，3.是对rho descending的起点
"""


def top_level_product_state(M,S):
    """
        定以倒数第二层两个site的纯态密度矩阵
        在finite_site.py 中调用过一次
        既然是product state viz. 非纠缠态, 即需要满足 tr(rho^2)=1, 再加上对称性限制, rho的形式就限制很大
        这里讲乘机态要注意是哪几部分之间的乘积态
        return:
           rho = 1/2(|uu><uu| + |dd><dd|)  这是个纠缠态, 并且是混态
    """
    iLayer = M.num_of_layer-1
    j=0
    S.rho_2[iLayer].O[j].data[:] = 0.0
    
    #rank = (2, 2)
    rank = S.rho_2[iLayer][j].rank
    qDims= np.empty(rank, "int")
    iDims= np.empty(rank, "int")
    
    # 0.5|00><00|
    qDims[:] = 0
    iDims[:] = 0
    S.rho_2[iLayer][j].set_element(qDims, iDims, 0.5)
    
    # 0.5|11><11|    
    qDims[:] = 1
    iDims[:] = 0
    S.rho_2[iLayer][j].set_element(qDims, iDims, 0.5)
    
    if not S.only_NN:
        warnings.warn("here buffer may be needed; rho_3 thus define correct?")
        
        temp = M.V[iLayer-1][j]
        QSp = temp.QSp[-1].copy_many(2)
        QSp[1].reverse()
        
        #topQN = 0
        #I = 1.0*|0><0|
        topQN = M.topQN
        I = iTensor(2, QSp, topQN)

        QSp = temp.QSp[-1].copy_many(2)
        QSp[1].reverse()
        iqn = QSp[0].has_quant_num(topQN)
        if iqn<0:
            raise Exception("Do not find suitable topQN")
        qDims= np.empty(rank, "int")
        iDims= np.empty(rank, "int")
        qDims[:] = iqn
        iDims[:] = 0
        I.data[:] = 0.0
        I.set_element(qDims, iDims, 1.0)
        S.rho_3[iLayer].O[j]=I.direct_product(I, use_buf=True).direct_product(I, use_buf=True)
    
def top_level_product_state_u1(M,S):
    """
        实际上是把次顶层V tensor的朝上的leg约化为1d，设置了次顶层V tensor的top quant number
        structure of rho_2 top is defined in mera.py,  this func assign its data
        在finite_site.py 中调用过一次
        既然是product state viz. 非纠缠态, 即需要满足 tr(rho^2)=1, 再加上对称性限制, rho的形式就限制很大
        这里讲乘机态要注意是哪几部分之间的乘积态
        
        对于finite range算法，这样选取的rho实质上起到一种投影算子的作用(在QM中，求算子期望值时，rho本身就是投影算子)，
        把次最上层的上面的腿都投影到|0>态，这样顶层就是Vidal文中所讲的 |0>^{\otimes N'} tensor product state
        return:
            1.0*|00><00|
    """
    iLayer = M.num_of_layer-1
    j=0
    #print "rrr", iLayer, S.rho_2[iLayer]
    #exit()
    S.rho_2[iLayer][j].data[:] = 0.0
    temp = M.V[iLayer-1][j]
    
    #QSp = temp.QSp[-1].copy_many(2)
    #QSp[1].reverse()
    QSp = temp.QSp[-1].copy_many(2, reverse=[1])
    topQN = M.topQN     #topQN = qn_identity 
    rank = 2
    I = iTensor(rank, QSp, topQN)   #I = 1.0*|0><0|
    
    #iqn = QSp[0].has_quant_num(topQN)
    #wordaround for decorators of copy
    iqn = temp.QSp[-1].has_quant_num(topQN)

    if iqn<0:
        tmpl= """
                not able to find suitable topQN
                topQN is: %(topQN)s
                QSp is %(QSp)s:
            """
        msg = tmpl%{"topQN":topQN, "QSp":temp.QSp[-1]}
        raise Exception(msg)
    qDims= np.empty(rank, "int")
    iDims= np.empty(rank, "int")
    qDims[:] = iqn
    iDims[:] = 0
    I.data[:] = 0.0
    I.set_element(qDims, iDims, 1.0)
    
    #rho_2 = 1.0*|00><00| 
    S.rho_2[iLayer].O[j] = I.direct_product(I, use_buf=True)
    
    if not S.only_NN:
        S.rho_3[iLayer].O[j]=I.direct_product(S.rho_2[iLayer].O[j], use_buf=True)

def top_level_product_state_u1_1site(M,S):
    """
    实际上是把次顶层V tensor的朝上的leg约化为1d，设置了次顶层V tensor的top quant number
    structure of rho_2 top is defined in mera.py,  this func assign its data
    在finite_site.py 中调用过一次
    既然是product state viz. 非纠缠态, 即需要满足 tr(rho^2)=1, 再加上对称性限制, rho的形式就限制很大
    这里讲乘机态要注意是哪几部分之间的乘积态
    
    对于finite range算法，这样选取的rho实质上起到一种投影算子的作用(在QM中，求算子期望值时，rho本身就是投影算子)，
    把次最上层的上面的腿都投影到|0>态，这样顶层就是Vidal文中所讲的 |0>^{\otimes N'} tensor product state
    return:
        1.0*|00><00|
    """
    warnings.warn("this func is not completely right, there is actually no correct topqn for U1-1site")
    iLayer = M.num_of_layer-1
    j=0
    rank = 2
    S.rho_2[iLayer].O[j].data[:] = 0.0
    temp = M.V[iLayer-1][j]
    if 0:
        t = iTensor(0, [QspU1.easy_init([1], [2])], QnU1(1))
        print t
        exit()
    
    
    S.rho_2[iLayer].O[j].data[:] = 0.0
    
    #rank = (2, 2)
    rank = S.rho_2[iLayer].O[j].rank
    qDims = np.empty(rank, "int")
    iDims = np.empty(rank, "int")
    
    # 0.5|00><00|
    qDims[:] = 1 #[0, 1, 1, 0]#[0, 1, 0, 1]
    iDims[:] = 0
    S.rho_2[iLayer][0].set_element(qDims, iDims, 0.5)
    
    # 0.5|11><11|    
    qDims[:] = 2 #[1, 0, 0, 1]#[1, 0, 1, 0]
    iDims[:] = 0
    S.rho_2[iLayer][0].set_element(qDims, iDims, 0.5)
    print S.rho_2[iLayer][0]
    #exit()
    
    
    if not S.only_NN:
        pass
        #S.rho_3[iLayer].O[j]=I.direct_product(S.rho_2[iLayer].O[j], use_buf=True)


class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    
    def test_eigen_state(self): 
        pass 
        symm = 'Travial'
        symm = 'U1'
        M = Mera.example(trunc_dim=4, tot_layer=4, symmetry=symm); 
        S= System.example(M, symm, 'Heisenberg')
        #top_level_eigenstate(M, S)
        print_vars(vars(),  ['M', 'S'])
        import decorators
        from merapy.ascending import ascending_ham 
        #np.set_printoptions(precision=3)
        decorators.tensor_player.STATE = 'stop'
        if 1: 
            for i in range(M.num_of_layer-1):
                #ilayer bellow 0 and >=M.num_of_layer-1 are not allowed
                decorators.set_STATE_end_simple(i, M.num_of_layer-1, iter0=0)
                ascending_ham(M, S, ilayer=i, info=0)
        decorators.tensor_player.STATE = 'stop'
        print_vars(vars(),  ['S'])
        top_level_eigenstate(M, S)
        #print_vars(vars(),  ['S'])
        #self.assertAlmostEqual()
        old=np.array([ 0.00184772, -0.01285457, -0.01285457,  0.03879706, -0.01285457, 0.08942899,  0.08942899, -0.2699104 , -0.01285457,  0.08942899])
        self.assertTrue(np.allclose(S.rho_2[3][0].data[: 10], old, atol=1e-8))
       
    def test_temp(self): 
        pass 
        symm = 'Travial'
        #symm = 'U1'
        M = Mera.example(trunc_dim=4, tot_layer=4, symmetry=symm); 
        S= System.example(M, symm, 'Heisenberg', info=2)
        #top_level_eigenstate(M, S)
        #print_vars(vars(),  ['M', 'S'])
        import decorators
        from merapy.ascending import ascending_ham 
        #np.set_printoptions(precision=3)
        decorators.tensor_player.STATE = 'stop'
        S.expand_layer(5)
        if 1: 
            for i in range(M.num_of_layer-1):
                #ilayer bellow 0 and >=M.num_of_layer-1 are not allowed
                decorators.set_STATE_end_simple(i, M.num_of_layer-1, iter0=0, )
                ascending_ham(M, S, ilayer=i, info=-1)
        decorators.tensor_player.STATE = 'stop'
        print_vars(vars(),  ['S'])
        top_level_eigenstate(M, S)
        
        #print_vars(vars(),  ['S'])
        #self.assertAlmostEqual()
        old=np.array([ 0.00184772, -0.01285457, -0.01285457,  0.03879706, -0.01285457, 0.08942899,  0.08942899, -0.2699104 , -0.01285457,  0.08942899])
        #self.assertTrue(np.allclose(S.rho_2[3][0].data[: 10], old, atol=1e-8))
       
    

if __name__ == "__main__":
    if 0: 
        from mera import test_Mera
        if 0:
            M = test_Mera.instance(trunc_dim=2, tot_layer=4, symmetry="Z2")
            sys= test_System.instance(M, model="Ising", symmetry="Z2", only_NN=False)
            
            #top_level_ham2(M, sys, M.topQN)
            print M.V[3][0]
            top_level_eigenstate(M, sys)
            #top_level_product_state(M, sys)
            print sys.rho_2[3].O[0]
            print sys.rho_3[3].O[0]
        if 1:
            M = test_Mera.instance(trunc_dim=3, tot_layer=4, symmetry="Z3")
            sys = test_System.instance(M, model="Potts", symmetry="Z3", only_NN=False)
            
            #top_level_ham2(M, sys, M.topQN)
            #print M.V[3][0]
            #top_level_eigenstate(M, sys)
            top_level_product_state_u1(M, sys)
            print sys.rho_2[3][0]
            print sys.rho_3[3][0]

        if 0:
            M = test_Mera.instance(trunc_dim=4, tot_layer=4, symmetry="U1")
            sys= test_System.instance(M, model="Heisenberg", symmetry="U1", only_NN=True)
            top_level_product_state_u1_1site(M, sys)
    
    if 0: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_eigen_state', 
           'test_temp', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)
       

        
        



