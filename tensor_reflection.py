#coding=UTF8

"""
MingQi: QSp_max = QSp_0
QSp_max%Dims(1:2) = 2
QSp_max%RefQN(1,1:2) = (/2,1/)
QSp_max%RefQN(2,1:2) = (/0,1/)

QSp_max%RefQN 几个数加起来要等于总维数

RefQN is for the reflection symmetry
10:46 AM 
real space reflection symmetry

MingQi: RefQN是反演对称性的维度
  
反演对称性有两个量子数，+-1
3:29 PM 
需要给出对应于 Z2对称性下反演对称性的维数。例如QSp_max%RefQN(1,1:2) = (/2,1/) 表示反演对称性为+1时，Z2为+1和-1的维度分别为2和1
  
QSp_max%RefQN(2,1:2) = (/0,1/)表示反演为-1时，Z2为+1和-1的维度分别是0，1

MingQi: Reflection_Sign是空间反演的符号
7:54 PM 
Spin_Reflection是处理自旋反演
  
不记得在Ising model中是否使用了Spin Reflection
7:55 PM 
在没有外磁场时，S->-S（or Sz->-Sz）体系不

state space  marked by two qns, spin parity p and reflection r, i.e.  |p,r>
"""



from quantum_number import *
from tensor import *
import numpy as np
from graphics  import Graphics
#import hamiltonian 
from tensor_network import TensorNetwork

class TensorReflect(object):
    @staticmethod
    def reflect(T, div, buffer=None):
        """
        this is the real space reflect of the spin chain, 
        specifically,  permute the legs in following way
        [0,1,...,div-1,   div,rank-1] ==> [div-1,div-2,...,0,   rank-1,rank-2,...,div]
        """
        rank = T.rank
        P = np.ndarray(rank,dtype="int")
        for i in range(div):
            P[i] = div-1-i
        for i in range(div, rank):
            P[i] = rank-1+div-i
        
        #print "pppp", P
        #tensor_Permutation(T, T2, P)
        T2 = T.permutation(P, buffer=buffer)
        return T2

    @staticmethod    
    def reflect_sign1(QSp0):
        """
        根据QSp0构造一个type(1,1)张量 A^{r'p'}_{rp}, 这个张量的作用 
        为下面的reflect_symmetrize做准备.
        A= A{r'p'}_{rp}|r',p'><p,r|, s.t. 
           A{r'p'}_{rp} = \delta_{p',p}\delta_{r',r}, if r=+1 
           A{r'p'}_{rp} = -\delta_{p',p}\delta_{r',r}, if r=-1
        实际上，也可以不定义这样的张量，直接在另外的张量上进行操作，进行refelct           
        """
        totQN = QN_idendity.copy()
        QSp = [QSp0.copy(),QSp0.copy()]
        QSp[0].reverse()
        #note here reverse qsp[0], while not qsp[1], 
        #this will be made clear once see the graph defiend in reflect_symmetrize below

        #this is in fact an operation on  matrix
        T=iTensor(2, QSp, totQN, use_buf=True)
        T.data[:] = 0.0  #this line is redundent since data is set to 0 in iTensor.__init__
        
        #从p的块，到r的块
        for i  in range(QSp[0].nQN):
            #get entrance for none zero data blocks
            #only blocks numberd i*QSp[0].nQN+i is none zero, i.e. q=q' 
            pidx=T.idx[i*QSp[0].nQN+i]
            pidx=T.Block_idx[0, pidx]
            for j  in range(QSp[0].Dims[i]):
                #get entrance for a data entry in the block
                #similarly only r=r' is none zero
                p = j*QSp[0].Dims[i]+j
                #print "ppp", p
                if j < QSp[0].RefQN[0,i]:
                    #reflect qn=+1
                    T.data[pidx+p] = 1.0
                else:
                    #reflect qn=-1
                    T.data[pidx+p] = -1.0
        return T
    @staticmethod
    def reflect_sign_old(T, T_sign):
        rank = T.rank
        for i  in range(T.rank):
            TensorReflect.reflect_sign1(T.QSp[i], T_sign.tensor[i])

    @staticmethod
    def reflect_sign(T):
        rank = T.rank
        T_sign = []
        for i  in range(T.rank):
            T_sign.append(TensorReflect.reflect_sign1(T.QSp[i]))
        return T_sign

    @staticmethod
    def reflect_symmetrize(T, T2=None, ndiv=None):
        div = T.ndiv
        if ndiv!=None:  
            div = ndiv
        rank = T.rank
        P = np.ndarray(rank,dtype="int")
        for i in range(div):
            P[i] = div-1-i
        for i in range(div, rank):
            P[i] = rank-1+div-i

        #initialize_Tensor(T1,use_buf= True )
        #initialize_Tensor(T3,use_buf= True )
        
        #delete_TensorNetwork(TNet_sys)
        #init_TensorNetwork(TNet_sys)
        
        #define a graph
        #each nodes have two edges
        G = Graphics()
        order = np.ndarray(rank+1, "int")
        forder = np.ndarray(rank, "int")
        G.size = rank+1
        
        for i  in range(rank):
            G.nodes[i] = 2
            G.edges[0:2,i] = [i,-i-1]  #not [i, -i], simply to avoid the case of [0, 0]             
            order[i+1] = i
            forder[i] = -i-1
        i = rank
        G.nodes[i] = rank
        for j  in range(rank):
            G.edges[j,i] = j
        order[0] = rank  #a ring
        
        #hamiltonian.System.TNet_sys = TensorNetwork() 
        #TNet_sys = hamiltonian.System.TNet_sys
        #attention_low_efficiency  
        #construct TensorNetwork each time may be  bad
        TNet_sys= TensorNetwork()
        TNet_sys.tlink.tensor[:rank] = TensorReflect.reflect_sign(T) 
        TNet_sys.tlink.tensor[rank]=T.shallow_copy()
        TNet_sys.TG = G
        
        #对T的每条腿进行空间reflection操作，this amongts to a contraction operation
        Tp = TNet_sys.contract_except(n=rank+1, order=order, 
                final_order=forder, restart=True ,info = 0)
        
        if T2!=None:  
            T2=Tp.permutation(P, buffer=T2.data.data)
        else:
            T1=Tp.permutation(P, use_buf=True)
            T3 = T
            #T.data[i] = 0.5d0*[T1.data[i]+T3.data[i]]
            T.data[:] = T1.data[:]+T3.data[:]                
        
    
    def reflect_Symmetrize_Mapping(T, RefMap, ndiv):
        initialize_Tensor(T2, use_buf= True )
        #forall[i=0:T.totDim] T.data[i]=dble[i]
        
        if ndiv!=None:  
            reflect_Symmetrize(T,T2, ndiv)
        else:
            reflect_Symmetrize(T,T2)
        
        #forall[i=0:T.totDim] RefMap[i]=nint[T2.data[i]]
        
        delete_Tensor(T2)
    
    def Check_reflect(T, T2, ndiv):
        pass
    
    def Check_Hem(T):
        rank = T.rank
        div = rank/2
        #forall[i=0:div] P[i] = div+i
        #forall[i=div+0:rank] P[i] = i-div
        
        tensor_Permutation(T, Tp, P)
        for i  in range(T.totDim):
            X = T.data[i]-Tp.data[i]
            if abs[X] >= 1.0-10:  
                print "['i=',I0,'  X=', 200[1x,E15.8]]", i, T.data[i], Tp.data[i], X
        delete_Tensor(Tp)
    
    def spin_reflect(T, Tp):
        """
        the reflection in the internal spin space, i.e., reflect
        between spin up and down
        """
        rank = T.rank
        QSp[0:rank] = T.QSp[0:rank]
        totQN = T.totQN.copy()        
        P = np.empty((ran, rank), "int")
        #P[i, r], i refers to a qn, r a leg
        #把每条腿的自旋反向, 记录在P中
        for r  in range(rank):
            MaxNQN = T.QSp[r].nQN
            if T.QSp[r].QNs[0].val == 0:  
                P[0,r] = 1
                #do i = 2, MaxNQN, 2
                for i in range(1,MaxNQN,2):
                    P[i,r] = i+1
                    P[i+1,r] = i
            else:
                for i  in range(0, MaxNQN, 2):
                    P[i,r] = i+1
                    P[i+1,r] = i
        
        #init_Tensor(rank, QSp, totQN, Tp)
        Tp = iTensor(rank, QSp, totQN)
        #Tp.data = 0.0
        iQN1 = np.ndarray(rank, "int")
        iQN2 = np.ndarray(rank, "int")
        for idx1 in range(T.nidx):
            iQN1[0:rank] = T.Addr_idx[0:rank, idx1]
            p1 = T.Block_idx[0,idx1]
            d1 = T.Block_idx[1,idx1]
            
            for i  in range(rank):
                #新的量子数组合
                iQN2[i] = P[iQN1[i],i]
            
            #iTensor_GetPosition(idx2, iQN2, Tp)
            idx2 = Tp.get_position(iQN2)

            idx2 = Tp.idx[idx2]
            p2 = Tp.Block_idx[0, idx2]
            d2 = Tp.Block_idx[1, idx2]
            #新的p2, d2
            if d1 != d2:  
                print 'Error, size not match', idx1, idx2, d1, d2
            #新的数据, over
            Tp.data[p2:p2+d2-1] = T.data[p1+d1-1:p1:-1]
    
    def Check_Spin_reflect(T):
        initialize_Tensor(Tp, use_buf= True )
        Spin_reflect(T,Tp)
        for i  in range(T.totDim):
            X = T.data[i]-Tp.data[i]
            if abs[X] >= 1.0-10:  
                print 'i,X=', i, T.data[i], X
        delete_Tensor(Tp)
    
    def spin_reflect_symmetrize(T):
        Tp = spin_reflect(T)
        T.data[:]=0.50*(T.data[:]+Tp.data[:])
    
    def Test_reflect():
        QSp[1].nQN = 3
        QSp[1].QNs[0:3].val = [0,1,-1]
        QSp[1].Dims[0:3] = [2,1,1]
        QSp[1].RefQN[1,0:3] = [1,1,1]
        QSp[1].RefQN[2,0:3] = [1,0,0]
        QSp_Update(QSp[1])
                
        QSp[2:3] = QSp[1]
        QSp_Reverse(QSp[3])
        QSp[4] = QSp[3]
        totQN = QN_idendity
        rank = 4
        init_Tensor(rank, QSp, totQN, T)
        #forall[i=0:T.totDim] T.data[i]=dble[i]
        
        #allocate[RefMap[T.totDim]]
        T.ndiv = 2
        reflect_Symmetrize_Mapping(T, RefMap)
        Random_Tensor(T)
        
        reflect_Symmetrize(T)
#        T2 = T
#        for i  in range(T.totDim):
#            T2.data[i] = 0.5d0*[T.data[i]+T1.data[i]]
        
#        T2.ndiv = 2
        iTensor_SVD(T)
        reflect_Symmetrize(T,T1)
        print_iTensor2(T,T1, True )


class test_reflect():
    @staticmethod
    def test1():
        from tensor import test_iTensor 
        u= test_iTensor.instance("u")
        #print u
        tr=TensorReflect
        #ur=tr.reflect(u,2)
        #print ur
        #us=tr.reflect_sign1(QSp_base)
        u.data[:] = np.arange(u.totDim)
        print u.data
        tr.reflect_symmetrize(u, ndiv=1)
        print u.data
    @staticmethod
    def test2():
        TR=TensorReflect
        from mera import test_Mera
        from hamiltonian import test_System
        from tensor_svd import Tensor_svd_new

        M=test_Mera.instance(trunc_dim=2)
        sys= test_System.instance(M)

        #print M
        print sys
        
        u0=M.U[0].tensor[0]
        h0 = sys.H_2[0].O[0]
        v0 = M.V[0].tensor[0] 
        L = locals()
        for i in ["u0", "h0", "v0"]:
            print i
            print L[i].data
            ndiv = 2
            if i == "v0": ndiv = 3 
            TR.reflect_symmetrize(L[i], ndiv=ndiv)
            Tensor_svd_new.svd(L[i], ndiv=ndiv)
            print L[i].data

    @staticmethod
    def test3():
        TR=TensorReflect
        from mera import test_Mera
        from hamiltonian import test_System
        from tensor_svd import Tensor_svd_new

        M=test_Mera.instance(trunc_dim=2)
        sys= test_System.instance(M)

        #print M
        #print sys
        
        v0 = M.V[0].tensor[0] 
        v1 = M.V[1].tensor[0] 
        v2 = M.V[2].tensor[0]
        v0.data[:]=v1.data[:]
        L = locals()

        q=v0.QSp[0]
        
        print "#"*100
        print q
        qs=q.add(q) #.add(q)
        #问题处在refqn 不应故等于0
        print qs
        qs.update(QSp_max=M.QSp_max)
        print qs
        if 0:    
            for i in range(4):
                print v0.QSp[i].__repr__(exclude=["QNs","Dims","nQN"])
            print "\n\n"
            for i in range(4):
                print v2.QSp[i].__repr__(exclude=["QNs","Dims","nQN"])

        if 0:
            for i in ["v0", "v1"]:
                print i, 
                print "before\n", L[i].data
                ndiv = 3
                TR.reflect_symmetrize(L[i], ndiv=ndiv)
                #Tensor_svd_new.svd(L[i], ndiv=ndiv)
                print i, "after\n", L[i].data
                print "\n"



if __name__ == "__main__":

    tr=test_reflect
    #tr.test2()
    tr.test3()
    if 0:
        a = QSp_base.add(QSp_base).add(QSp_base)
        print a
        a.update(QSp_max=QSp_base)
        print a


    """
        --- pass
        [ 0.  1.  2.  3.  4.  5.  6.  7.]
        [  0.   5.   4.   9.   5.  10.   9.  14.]
    """





    #results
    """
        using pure python implementation of iTensor
        using pure python implementation of iTensor
        v0 before
        ----Begin iTensor----------------------------------------------
        rank:	4
        idx_dim:	16
        ndiv:	3
        nidx:	8
        totQN:	(    1)
        QNs:	['[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]']
        Dims:	array([ 2.,  2.,  2.,  2.])
        totDim:	8
        idx:	array([ 0, 16, 16,  1, 16,  2,  3, 16, 16,  4,  5, 16,  6, 16, 16,  7])
        Block_idx:	[0 1 2 3 4 5 6 7]
            [1 1 1 1 1 1 1 1]
            [ 0  3  5  6  9 10 12 15]

        Addr_idx:	[[0 1 1 0 1 0 0 1]
         [0 1 0 1 0 1 0 1]
         [0 0 1 1 0 0 1 1]
         [0 0 0 0 1 1 1 1]]
        data:	[0 0 0 0]: [-0.68966]
        [1 1 0 0]: [-0.2565]
        [1 0 1 0]: [-0.11678]
        [0 1 1 0]: [ 0.66704]
        [1 0 0 1]: [-0.11111]
        [0 1 0 1]: [-0.51846]
        [0 0 1 1]: [-0.76133]
        [1 1 1 1]: [ 0.37315]

        ----End iTensor----------------------------------------------

        v0 after
        ----Begin iTensor----------------------------------------------
        rank:	4
        idx_dim:	16
        ndiv:	3
        nidx:	8
        totQN:	(    1)
        QNs:	['[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]']
        Dims:	array([ 2.,  2.,  2.,  2.])
        totDim:	8
        idx:	array([ 0, 16, 16,  1, 16,  2,  3, 16, 16,  4,  5, 16,  6, 16, 16,  7])
        Block_idx:	[0 1 2 3 4 5 6 7]
            [1 1 1 1 1 1 1 1]
            [ 0  3  5  6  9 10 12 15]

        Addr_idx:	[[0 1 1 0 1 0 0 1]
         [0 1 0 1 0 1 0 1]
         [0 0 1 1 0 0 1 1]
         [0 0 0 0 1 1 1 1]]
        data:	[0 0 0 0]: [ 0.]
        [1 1 0 0]: [-0.92354]
        [1 0 1 0]: [ 0.]
        [0 1 1 0]: [ 0.92354]
        [1 0 0 1]: [ 0.65021]
        [0 1 0 1]: [ 0.]
        [0 0 1 1]: [-0.65021]
        [1 1 1 1]: [ 0.]

        ----End iTensor----------------------------------------------

        v1 before
        ----Begin iTensor----------------------------------------------
            rank:	4
            idx_dim:	16
            ndiv:	3
            nidx:	8
            totQN:	(    1)
            QNs:	['[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]']
            Dims:	array([ 2.,  2.,  2.,  2.])
            totDim:	8
            idx:	array([ 0, 16, 16,  1, 16,  2,  3, 16, 16,  4,  5, 16,  6, 16, 16,  7])
            Block_idx:	[0 1 2 3 4 5 6 7]
                [1 1 1 1 1 1 1 1]
                [ 0  3  5  6  9 10 12 15]

            Addr_idx:	[[0 1 1 0 1 0 0 1]
             [0 1 0 1 0 1 0 1]
             [0 0 1 1 0 0 1 1]
             [0 0 0 0 1 1 1 1]]
            data:	[0 0 0 0]: [-0.42095]
            [1 1 0 0]: [-0.15087]
            [1 0 1 0]: [-0.43825]
            [0 1 1 0]: [ 0.77973]
            [1 0 0 1]: [-0.16802]
            [0 1 0 1]: [ 0.30662]
            [0 0 1 1]: [ 0.8928]
            [1 1 1 1]: [ 0.28401]

        ----End iTensor----------------------------------------------

        v1 after
        ----Begin iTensor----------------------------------------------
            rank:	4
            idx_dim:	16
            ndiv:	3
            nidx:	8
            totQN:	(    1)
            QNs:	['[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]', '[(    1) (   -1)]']
            Dims:	array([ 2.,  2.,  2.,  2.])
            totDim:	8
            idx:	array([ 0, 16, 16,  1, 16,  2,  3, 16, 16,  4,  5, 16,  6, 16, 16,  7])
            Block_idx:	[0 1 2 3 4 5 6 7]
                [1 1 1 1 1 1 1 1]
                [ 0  3  5  6  9 10 12 15]

            Addr_idx:	[[0 1 1 0 1 0 0 1]
             [0 1 0 1 0 1 0 1]
             [0 0 1 1 0 0 1 1]
             [0 0 0 0 1 1 1 1]]
            data:	[0 0 0 0]: [-0.84191]
            [1 1 0 0]: [ 0.62886]
            [1 0 1 0]: [-0.8765]
            [0 1 1 0]: [ 0.62886]
            [1 0 0 1]: [ 0.72478]
            [0 1 0 1]: [ 0.61324]
            [0 0 1 1]: [ 0.72478]
            [1 1 1 1]: [ 0.56802]

        lll           1
            T%data= ,-.4210E+00 ,-.1509E+00 ,-.4383E+00 ,0.7797E+00 ,-.1680E+00 ,0.7797E+00 ,-.1680E+00 ,0.3066E+00 ,0.8928E+00 ,0.2840E+00
            T%data= ,0.0000E+00 ,-.9306E+00 ,0.0000E+00 ,0.9306E+00 ,-.1061E+01 ,0.9306E+00 ,-.1061E+01 ,0.0000E+00 ,0.1061E+01 ,0.0000E+00

        ----End iTensor----------------------------------------------


    """

