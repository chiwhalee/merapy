#! /usr/bin/env python
#coding=UTF8
"""
    see tensornetwork.f90
    q:
        1.MaxLeg 是图中所有的leg？ too few? 
        2. q-140  why permutation in the end?
    todo: 
        support for contract by devide graph to several groups
"""
from tensor import *
from graphics import *
from set1 import OrderSet

__all__ = ["TensorLink", "TensorNetwork"]

class TensorLink():
    """
    tlink is basically pointers or references to a set of tensors
    this class was orginally defined in tensors.f90, I moved it here
    """
    def __init__(self, size, use_buf=False):
        """
        status:complete

        just construct n Tensors, no link yet

        see Init_TLink in f90
        use_buf: bool 

        """
        self.size= size 
        self.use_buf=use_buf
        #if use_buf:self.tensor=[TensorBase(use_buf=use_buf) for i in range(size)  ]
        #else:self.tensor=[TensorBase() for i in range(size)  ]
        self.tensor=[None for i in range(size)]

    def delete(self):
        """
        status:half
        see Delete_Tlink in f90
        """

    def __getitem__(self, i):
        return self.tensor[i]

    def __setitem__(self, i, x):
        self.tensor[i] = x

    def Init_TLink(tlink, n, use_buf):
        self.size = n
        #allocate[self.tensor[n]]
        if use_buf!=None:
            for i  in range(n):
                pass
                #initialize_Tensor[self.tensor(i], use_buf=use_buf)
                #attention_omitted_something
            
        else:
            for i  in range(n):
                pass
    
    def delete(self):
        for i  in range(self.size):
            self.tensor[i].delete()


class TensorNetwork(object):
    MaxLeg = 64   #see quantum_number.py
    DATA_BUFFER = None
    BUFFER_ON = False
    
    @classmethod
    def init_buffer(cls, size, num):
        cls.DATA_BUFFER = {}
        ranks= [4, 6]
        for r in ranks:
            cls.DATA_BUFFER[r] = [np.empty(size) for n in range(num)]
        cls.BUFFER_ON = True

    def __init__(self):
        """
        status_1_verified
        see init_TensorNetwork in f90
        """
        self.nLeft = 0
        self.nRight = 0
        self.tlink = TensorLink(size=self.MaxLeg)
        self.tlink_left = TensorLink(size=self.MaxLeg)
        self.tlink_right= TensorLink(size=self.MaxLeg)
        self.edges_left=np.ndarray((Graphics.MaxEdge,Graphics.MaxNode), "int")
        self.edges_right=np.ndarray((Graphics.MaxEdge,Graphics.MaxNode), "int")
        self.TG = None

    def delete(self):
        """
        deprecated
        see delete_TensorNetwork in f90
        """
        pass

    @staticmethod
    def contract_many(n,T,V, nord, order, nford, ford, info=0):
        """
        see Tensor_Contraction in f90
        this function is no more than contracting many a tensors according to a certain order
        Contract a set of tensor according to V and ord    
        n: number of tensors
        T: tensor list
        Tout output tensor
        V leg list
        ord: array for contraction order
        nord: is the dimension of ord array        
        ford: is the final order of legs which are not contracted
        """
        S=[]
        for i in range(n):
            S.append( V[:T.tensor[i].rank, i])
        #S保存所有腿的标记

        Sford = ford
        
        NumOfTensor = n
        while len(order)!=0:    #do contract with order shrinking until empty
            #print "aaa ", len(order), "order", order, "S",S 
            iLeg = order[0]   #选取stack中的第一条腿

            Is=iTensor.find_leg(iLeg, S)   #找到连接这条腿的两个tensor
            i1 = Is[0]    #i1, i2是两个被iLeg连接的tensor的编号
            i2= Is[1]
            temp, V3, Vc = T.tensor[i1].contract(T.tensor[i2], S[i1],  S[i2],
                     out_Vc=True, info=info-1)
            V[:temp.rank,NumOfTensor] = V3[:]
            T.tensor.append(temp)
            S1 = Vc  # Vc[:nVc]

            T.tensor.pop(i1)
            T.tensor.pop(i2-1)

            S.pop(i1)
            S.pop(i2-1)
            S.append(V3)
            S2 = OrderSet.diff(order, S1)    #把所有收缩掉的腿都去掉

            order = S2   #ord shrinks,  a new loop starts
        
        #attention_omitted_something
        #for i in range(Tout.rank):
            #Vc[InSet[ford[i], S[NumOfTensor]]] = i
        #Tout=T.tensor[NumOfTensor].permutation(Vc)
        #q140
        Tout=T.tensor[0]

        return Tout
        
    def contract_except(self, n,  order, exception=None, final_order=None, use_buf=True, restart=None, info=0):
        """
            always use buff
            see Contract_TensorNetwork in f90
            Contract n-node tensor network except the node `exception'.其效果是把exception节点剪掉，
            把剩余点所有节点收缩成一个节点，即左后输出的张量
            according to order `order'
            exception: 存储的是某个节点的编号，如果<0, 则不考虑
            n: number of tensors, or equivalently number of nodes in the graph
            returns:
                Tout: 
        """
        #for i in range(G.size):  print self.tlink.tensor[n].type_name

        TNet=self
        G = TNet.TG.copy()   #note the copy is necessary!

        order=order[:n]

        restart0 =  True 
        #attention_omitted_something
        order = list(order)
        if exception is not None:  
            if exception>=0:  
                nfinal = G.nodes[exception]   #最终节点 的腿的数目
                final_leg_ord =  G.edges[0:nfinal, exception]   #最终腿的顺序
                i = order.index(exception)
            #这里把exception node的两边分成两部分ord1,  ord2
            if i>= 0:  #判断 exception 在不在图中
                ord1=order[:i]
                ord2=order[i+1:n]
                ord2.reverse()
                #note here ord2 is reversed
            else:
                ord1=order
        else:
            ord1=order
            ord2 = []
        if info>0:
            print "ord", order
            print "ord1, ord2", ord1, ord2
            print [G.names[i] for i in ord1], [G.names[i] for i in ord2]


        nLeft = 0
        if ord1 != []:  
            TNet.edges_left[:,0] = G.edges[:, ord1[0]]
            #把leg编号 从 graph 赋值给 TensorNetwork
        while len(ord1)>1:   # contract until only one tensor left
            nLeft = nLeft+1  # 左边tensor的数目
            G.contract_nodes(ord1[0], ord1[1])   # 这里G.size  加了1
            if nLeft  >  TNet.nLeft:  
                TNet.nLeft = nLeft
                #合并出的tensor的序号为 nLeft
                v1 = TNet.edges_left[:, nLeft-1]
                v2 = G.edges[:, ord1[1]]
                temp, legs= TNet.tlink.tensor[ord1[0]].contract(TNet.tlink.tensor[ord1[1]], V1=v1, V2=v2, out_Vc=False, use_buf= use_buf, info=info-1)
                TNet.tlink_left.tensor[nLeft]=temp 
                rank = TNet.tlink_left.tensor[nLeft].rank
                TNet.edges_left[0:rank, nLeft] = legs[0:rank]
            
            TNet.tlink.tensor[G.size-1] = TNet.tlink_left.tensor[nLeft].shallow_copy()

            #ord1.insert(2,G.size-1)
            #print "ggg1", G
            ord1.pop(0)
            ord1.pop(0)
            ord1.insert(0, G.size-1)
        if info-1>0:
            print "left side is completed contracted\n"

        nRight = 0
        if ord2 != []:  
            TNet.edges_right[:,0] = G.edges[:, ord2[0]]
        while len(ord2)>1:
            nRight = nRight+1
            G.contract_nodes(ord2[0], ord2[1])
            if nRight  >=  TNet.nRight:  
                TNet.nRight = nRight
                v1 = TNet.edges_right[:, nRight-1]
                v2 = G.edges[:, ord2[1]]
                
                temp, legs= TNet.tlink.tensor[ord2[0]].contract(TNet.tlink.tensor[ord2[1]], v1, v2, 
                        out_Vc=None, use_buf= use_buf, info=info-1)
                TNet.tlink_right.tensor[nRight]=temp 
                rank = TNet.tlink_right.tensor[nRight].rank
                TNet.edges_right[0:rank, nRight] = legs[0:rank]
            TNet.tlink.tensor[G.size-1] = TNet.tlink_right.tensor[nRight].shallow_copy()
            
            #common_util.ishift(ord2, 1, -1)

            #ord1.insert(2,G.size-1)
            ord2.pop(0)
            ord2.pop(0)            
            ord2.insert(0, G.size-1)
            #print "ggg2", G
        
        if info-1>0:
            print "right side is completed contracted\n"
            
      
        #把左右两边点收缩掉
        if ord1 and ord2:  
            G.contract_nodes( ord1[0], ord2[0])
            #print "ttt ", type(TNet.tlink.tensor[ord1[0]]), ord1[0]
            v1 = TNet.edges_left[:, nLeft]
            v2 = TNet.edges_right[:, nRight]

            TNet.tlink.tensor[G.size-1], legs = TNet.tlink.tensor[ord1[0]].contract(TNet.tlink.tensor[ord2[0]], v1, v2, 
                    out_Vc=False, use_buf= use_buf, info=info-1)
            rank = TNet.tlink.tensor[G.size-1].rank
            #TNet.graph.edges[0:rank, G.size] = legs[0:rank]
        
        #rearange 即把最终的张量重新排一下
        need_rearrange =  False 

        if final_order is not None:              
            rank = G.nodes[G.size-1] #number of legs  of final tensor
            P=np.ndarray(rank, np.int)
            for i in range(rank): #TNet.graph.nodes[TNet.graph.size]):
                #j=final_order.index(legs[i])
                j = np.where(final_order==legs[i])
                P[j] = i
            need_rearrange =  True 
        elif exception is not None:  
            if exception >= 0:  
                rank = G.nodes[G.size-1]  
                P=np.ndarray(rank, np.int)
                for i in range(rank):
                    j=np.where(G.edges[:rank, exception]==legs[i])
                    P[j] = i
                need_rearrange =  True 

        if need_rearrange:  
            Tout=TNet.tlink.tensor[G.size-1].permutation(P, use_buf=use_buf)
        else:
            Tout = TNet.tlink.tensor[G.size-1]
        if info>0:
            print "\nentire tensor net is contracted, yielding:"
            #print "type_name: ", Tout.type_name, "\trank: ", Tout.rank, "\tleg label:", P[:Tout.rank] #G.edges[:Tout.rank, G.size-1]
            print "type_name:{0.type_name}\trank:{0.rank}\tlebel{label}\tdims:{0.Dims}".format(
                    Tout, label=P[:Tout.rank])
            print Tout.data[:4].round(5),"...", Tout.data[-4:].round(5)
        
        return Tout
    
    @staticmethod
    def make_tensor_link(S, G, trans_invar=0): 
        res = [None]*G.size
        for node_id in range(G.size): 
            node_name = G.names[node_id]
            name, lay, pos = node_name.split('_')
            if trans_invar: pos = 0
            if 1: 
                if name == 'rho': 
                    op = S.rho_2[lay][pos].copy()
                if name == 'V': 
                    op = S.mera.V[int(lay)][pos].copy()
                if name == 'Vp': 
                    op = S.mera.V_dag[int(lay)][pos].copy()
                if name == 'U': 
                    op = S.mera.U[int(lay)][pos].copy()
                if name == 'Up': 
                    op = S.mera.U_dag[int(lay)][pos].copy()
            
            op.type_name = node_name
            res[node_id] = op
        return res
       

class test_TensorNetwork(object):
    g3=Graphics()
    g5=Graphics()
    tn=TensorNetwork()
    #test_TensorNetwork.set_value_g3()
    @classmethod
    def set_value_g3(cls):
        g3 = cls.g3
        g3.size = 3
        g3.nodes[:]= -1
        g3.edges[:, :] = -1 

        g3.nodes[0] = 3
        g3.nodes[1] = 3
        g3.nodes[2] = 4        

        g3.edges[:g3.nodes[0], 0] = [0, 2, 4]
        g3.edges[:g3.nodes[1], 1] = [0, 1, 3]
        g3.edges[:g3.nodes[2], 2] = [1, 2,3,4]
        #print 'ggg', cls.g3
    
    @classmethod
    def set_value_g5(cls):
        g5 = cls.g5
        g5.size = 5
        g5.nodes[:]= -1
        g5.edges[:, :] = -1 

        #g5.nodes[0] = 3
        g5.nodes[:]=2

        g5.edges[:g5.nodes[0], 0] = [0, 4]
        g5.edges[:g5.nodes[1], 1] = [0, 1]
        g5.edges[:g5.nodes[2], 2] = [1, 2]
        g5.edges[:g5.nodes[3], 3] = [2, 3]
        g5.edges[:g5.nodes[4], 4] = [3, 4]
    
    @classmethod
    def contract_many(cls):
        """   --- seems correct """
        g3 = cls.g3
        tl=TensorLink(size=3)
        #t3a=  nTensor(1,[3], ind_labels=['b']);  t3a.data=np.arange(3,dtype='d')
        #t3b=  nTensor(1,[3], ind_labels=['b']);  t3b.data=np.arange(3,dtype='d')
        #t4=  nTensor(1,[4], ind_labels=['d']);  t4.data=np.arange(4,dtype='d')

        t222a= nTensor(3,[2,2,2]); t222a.data=np.arange(8,dtype='d')
        t222b= nTensor(3,[2,2,2]);  t222a.data=np.arange(8,dtype='d')            
        t2222= nTensor(4,[2,2,2,2]);  t2222.data=np.arange(16,dtype='d')            

        tl.tensor=[t222a,t222b,t2222]
        n=3
        T= tl
        V= g3.edges
        #print g3
        #ord=[0, 2, 1]
        ord=[0,1,2,3,4]
        nord=len(ord)
        ford=[5,6,7]
        nford=len(ford)

        to = TensorNetwork.contract_many(n,T,V, nord, ord, nford, ford)
        print to
    #test_contract_many()
    
    @classmethod
    def contract_except(cls):
        g5 = cls.g5
        tn = TensorNetwork()
        tn.TG=g5
        n=g5.size
        order=[0,2,4,5,1]
        exception= 4

        t22=[]
        for i in range(n):
            a= nTensor(2,[2,2]); a.data=np.arange(4,dtype='d')
            t22.append(a.copy())

        tn.tlink.tensor[:n]=t22[:n]
        to = tn.contract_except(n, order, exception)

        print to

    @classmethod
    def contract_except_itensor(cls):
        import tensor
        g5 = cls.g5
        tn = TensorNetwork()
        tn.TG=g5
        n=g5.size
        #order=[0,2,4,5,1]
        order = range(5)
        exception= 0

        t22=[]
        for i in range(n):
            #a= nTensor(2,[2,2]); a.data=np.arange(4,dtype='d')

            a = tensor.test_iTensor.instance(which="u22")
            t22.append(a.copy())
        tn.tlink.tensor[:n]=t22[:n]
        for i in range(n):
            tn.tlink.tensor[i].type_name = i

        #print "ttt", type(tn.tlink.tensor[0])
        to = tn.contract_except(n, order, exception, info=1)

        #print to


if __name__ == "__main__":

    
    ttn = test_TensorNetwork
    ttn.set_value_g3()
    ttn.set_value_g5()

    print 'ttt', ttn.g5
    #ttn.contract_except()
    ttn.contract_except_itensor()
    #iTensor.buff_free()





        

