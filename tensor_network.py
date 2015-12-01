#! /usr/bin/env python
#coding=UTF8
"""
    see tensornetwork.f90
    todo: 
        support for contract by devide graph to several groups
"""
import unittest 
import numpy as np 

from merapy.utilities import print_vars
from tensor import iTensor 
from graphics import Graphics 
from set1 import OrderSet

__all__ = ["TensorLink", "TensorNetwork"]

class TensorLink():
    """
    tlink is basically pointers or references to a set of tensors
    this class was orginally defined in tensors.f90, I moved it here
    """
    def __init__(self, size, use_buf=False):
        """
            just construct n Tensors, no link yet
            use_buf: bool 

        """
        self.size= size 
        self.use_buf=use_buf
        #if use_buf:self.tensor=[TensorBase(use_buf=use_buf) for i in range(size)  ]
        #else:self.tensor=[TensorBase() for i in range(size)  ]
        self.tensor=[None for i in range(size)]

    def __getitem__(self, i):
        return self.tensor[i]

    def __setitem__(self, i, x):
        self.tensor[i] = x

class TensorNetwork(object):
    """
        conceptially this is a network of tensors \hat{G} mapped from a graph G 
    """
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
        """
        self.nLeft = 0
        self.nRight = 0
        self.tlink = [None]*self.MaxLeg 
        self.tlink_left = [None]*self.MaxLeg 
        self.tlink_right = [None]*self.MaxLeg 
        self.edges_left = np.ndarray((Graphics.MaxEdge,Graphics.MaxNode), int)
        self.edges_right = np.ndarray((Graphics.MaxEdge,Graphics.MaxNode), int)
        self.edges_lists = []
        for i in range(3):  #support at most 3 sub group 
            self.edges_lists.append(np.ndarray((Graphics.MaxEdge,Graphics.MaxNode), int))
        
        self.TG = None

    @staticmethod
    def contract_tensor_list_old(T, V, order, nford, ford, info=0):
        """
            see Tensor_Contraction in f90
            this function is no more than contracting many a tensors according to a certain order
            Contract a set of tensor according to V and ord    
            params: 
                T: tensor list
                Tout:  output tensor
                V:  leg list
                ord: array for contraction order
                nord: is the dimension of ord array        
                ford: is the final order of legs which are not contracted
        """
        n = len(T)
        S = []
        for i in range(n):
            S.append( V[:T[i].rank, i])
        #S保存所有腿的标记

        length = n
        while len(order)!=0:    #do contract with order shrinking until empty
            iLeg = order[0]   #选取stack中的第一条腿

            Is=iTensor.find_leg(iLeg, S)   #找到连接这条腿的两个tensor
            i1 = Is[0]    #i1, i2是两个被iLeg连接的tensor的编号
            i2= Is[1]
            temp, V3, Vc = T[i1].contract(T[i2], S[i1],  S[i2],
                     out_Vc=True, info=info-1)
            V[:temp.rank,length] = V3[:]
            print_vars(vars(),  ['V3', 'Vc'])
            T.append(temp)
            S1 = Vc  # Vc[:nVc]

            T.pop(i1)
            T.pop(i2-1)

            S.pop(i1)
            S.pop(i2-1)
            S.append(V3)
            S2 = OrderSet.diff(order, S1)    #把所有收缩掉的腿都去掉

            order = S2   #ord shrinks,  a new loop starts
        
        #attention_omitted_something
        #for i in range(Tout.rank):
            #Vc[InSet[ford[i], S[length]]] = i
        #Tout=T[length].permutation(Vc)
        #q140
        Tout=T[0]

        return Tout
    
    @staticmethod
    def contract_tensor_list(tensor_list, labels_list, ford=None, use_buf=False, info=0):
        """
            see tensor_listensor_Contraction in f90
            this function is no more than contracting many a tensors according to a certain order
            Contract a set of tensor according to V and ord    
            params: 
                tensor_list: tensor list
                tensor_listout:  output tensor
                V:  leg list
                ord: array for contraction order
                nord: is the dimension of ord array        
                ford: is the final order of legs which are not contracted
        """
        t0 = tensor_list[0]
        label0 = labels_list[0]
        
        i = 0
        for t1 in tensor_list[1: ]: 
            i += 1 
            label1 = labels_list[i]
            t0, label0 = t0.contract(t1, label0,  label1, use_buf=use_buf, info=info-1)
        
        #attention_omitted_something
        #for i in range(tensor_listout.rank):
            #Vc[InSet[ford[i], S[length]]] = i
        #tensor_listout=tensor_list[length].permutation(Vc)
        #q140
        #tensor_listout=tensor_list[0]

        return t0, label0 
        
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
        #for i in range(G.size):  print self.tlink[n].type_name

        TNet=self
        G = TNet.TG.copy()   #note the copy is necessary!
        order=order[:n]

        
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
                temp, legs= TNet.tlink[ord1[0]].contract(TNet.tlink[ord1[1]], V1=v1, V2=v2, out_Vc=False, use_buf= use_buf, info=info-1)
                TNet.tlink_left[nLeft]=temp 
                rank = TNet.tlink_left[nLeft].rank
                TNet.edges_left[0:rank, nLeft] = legs[0:rank]
            
            TNet.tlink[G.size-1] = TNet.tlink_left[nLeft].shallow_copy()

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
                
                temp, legs= TNet.tlink[ord2[0]].contract(TNet.tlink[ord2[1]], v1, v2, 
                        out_Vc=None, use_buf= use_buf, info=info-1)
                TNet.tlink_right[nRight]=temp 
                rank = TNet.tlink_right[nRight].rank
                TNet.edges_right[0:rank, nRight] = legs[0:rank]
            TNet.tlink[G.size-1] = TNet.tlink_right[nRight].shallow_copy()
            
            #ord1.insert(2,G.size-1)
            ord2.pop(0)
            ord2.pop(0)            
            ord2.insert(0, G.size-1)

        
        if info-1>0:
            print "right side is completed contracted\n"
            
      
        #把左右两边点收缩掉
        if ord1 and ord2:  
            G.contract_nodes( ord1[0], ord2[0])
            #print "ttt ", type(TNet.tlink[ord1[0]]), ord1[0]
            v1 = TNet.edges_left[:, nLeft]
            v2 = TNet.edges_right[:, nRight]

            TNet.tlink[G.size-1], legs = TNet.tlink[ord1[0]].contract(TNet.tlink[ord2[0]], v1, v2, 
                    out_Vc=False, use_buf= use_buf, info=info-1)
            rank = TNet.tlink[G.size-1].rank
            #TNet.graph.edges[0:rank, G.size] = legs[0:rank]
        
        #rearange 即把最终的张量重新排一下
        need_rearrange =  False 
        P = None 
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
            Tout=TNet.tlink[G.size-1].permutation(P, use_buf=use_buf)
        else:
            Tout = TNet.tlink[G.size-1]
        if info>0:
            print "\nentire tensor net is contracted, yielding:"
            #print "type_name: ", Tout.type_name, "\trank: ", Tout.rank, "\tleg label:", P[:Tout.rank] #G.edges[:Tout.rank, G.size-1]
            print "type_name:{0.type_name}\trank:{0.rank}\tlebel{label}\tdims:{0.Dims}".format(
                    Tout, label=P)
            print Tout.data[:4].round(5),"...", Tout.data[-4:].round(5)
        
        return Tout

    def contract_except_test(self, n, order, exception=None, final_order=None, use_buf=True, restart=None, info=0):
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
        #for i in range(G.size):  print self.tlink[n].type_name
        
        if exception is None: 
            exception = []
        elif not hasattr(exception, '__iter__'): 
            exception = [exception]
            
        TNet = self
        G = TNet.TG.copy()   #note the copy is necessary!
        order=order[:n]
        #size = G.size - len(exception)
        #length = G.size - 1 
        length = G.size 
        order = list(order)
       
        #issue:  下面第一种收缩顺序不是最优的，
        #第二种, 当num_exc = 1时是最优的，但其他情况未必
        if 0: 
            for i in exception: 
                order.remove(i)
            tensor_list = self.tlink 
            o = order[0]
            head = tensor_list[o] 
            head_label = G.edges[: head.rank, o]
            for ii in order[1: ]: 
                t = tensor_list[ii]
                t_label = G.edges[: t.rank, ii]
                head, head_label = head.contract(t, head_label, t_label, info=info-1)
        else: 
            order = np.array(order, dtype=int)
            nested_order_list = []
            nested_label_list = []
            temp = [np.where(order==i)[0][0] for i in exception]
            #temp.append(length-1)
            temp1 = list(temp)
            temp.insert(0, -1)
            temp1.append(length)
            
            num_exc = len(exception)
            for i in range(num_exc + 1): 
                start = temp[i]
                end = temp1[i]
                if end-start>1: 
                    nested_order_list.append(order[start+1: end])
            
            tt = []
            ll = []
            for oo in nested_order_list: 
                if info>0: 
                    print 'contracting tensor_list {} ... '.format(oo)
                tlist = [self.tlink[i] for i in oo]
                llist = [G.edges[:self.tlink[i].rank, i] for i in oo]
                t, l = TensorNetwork.contract_tensor_list(tlist, llist, use_buf=use_buf, info=info)
                tt.append(t)
                ll.append(l)
            
            res, res_label = TensorNetwork.contract_tensor_list(tt, ll, use_buf=use_buf, info=info)
            head = res
            head_label = res_label 
        if 1: 
            #rearange 即把最终的张量重新排一下
            need_rearrange =  False 
            if final_order is not None:              
                raise NotImplemented 
                rank = G.nodes[G.size-1] #number of legs  of final tensor
                P = np.ndarray(rank, np.int)
                for i in range(rank): #TNet.graph.nodes[TNet.graph.size]):
                    #j=final_order.index(legs[i])
                    j = np.where(final_order==legs[i])
                    P[j] = i
                need_rearrange =  True 
                
            elif exception: 
                rank  = head.rank 
                final_label = np.ndarray(rank, dtype=int)
                p = np.ndarray(rank, dtype=int)
                i0 = 0
                for e in exception: 
                    r = self.tlink[e].rank 
                    final_label[i0: i0 + r] =  G.edges[:r, e]
                    i0 = r 
                #final_label = final_label[: head.r]
                for i in range(rank):
                    j=np.where(final_label==head_label[i])
                    p[j] = i
                need_rearrange = True 
            if need_rearrange:  
                res=head.permutation(p, use_buf=use_buf)
            else:
                res = head 
            if info>0:
                print "\nentire tensor net is contracted, yielding:"
                #print "type_name: ", res.type_name, "\trank: ", res.rank, "\tleg label:", P[:res.rank] #G.edges[:res.rank, G.size-1]
                print "type_name:{0.type_name}\trank:{0.rank}\tlebel{label}\tdims:{0.Dims}".format(
                        res, label=p)
                print res.data[:4].round(5),"...", res.data[-4:].round(5)
        
        return res 

    def contract_with_removals(self, n,  order, exception=None, final_order=None, use_buf=True, restart=None, info=0):
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
        
        exception = exception if exception is not None  else []

        TNet=self
        G = TNet.TG.copy()   #note the copy is necessary!
        order=order[:n]
        
        order = list(order)
        tensor_lists = []
        order_lists= []
        if 0: 
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
        
        if len(exception)>0: 
            start = 0
            for e in exception: 
                pass 
                #if exception>=0:  
                sub_list = []
                assert exception >= 0 
                if 1: 
                    nfinal = G.nodes[e]   #最终节点 的腿的数目
                    final_leg_ord =  G.edges[0:nfinal, e]   #最终腿的顺序
                    i = order.index(e)
                #这里把exception node的两边分成两部分ord1,  ord2
                    assert i >= 0
                    if start<i: 
                        #tensor_lists.append()
                        print_vars(vars(),  ['start'])
                        order_lists.append(order[start: i])
                        
                    start = i+1
                #if i>= 0:  #判断 exception 在不在图中
                #    ord1=order[:i]
                #    ord2=order[i+1:n]
                #    ord2.reverse()
                #    #note here ord2 is reversed
            
            if start<n: 
                order_lists.append(order[start: ])
            print_vars(vars(),  ['tensor_lists', 'order_lists', 'n'])
            
        if 0: 
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
                    temp, legs= TNet.tlink[ord1[0]].contract(TNet.tlink[ord1[1]], V1=v1, V2=v2, out_Vc=False, use_buf= use_buf, info=info-1)
                    TNet.tlink_left[nLeft]=temp 
                    rank = TNet.tlink_left[nLeft].rank
                    TNet.edges_left[0:rank, nLeft] = legs[0:rank]
                
                TNet.tlink[G.size-1] = TNet.tlink_left[nLeft].shallow_copy()

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
                    
                    temp, legs= TNet.tlink[ord2[0]].contract(TNet.tlink[ord2[1]], v1, v2, 
                            out_Vc=None, use_buf= use_buf, info=info-1)
                    TNet.tlink_right[nRight]=temp 
                    rank = TNet.tlink_right[nRight].rank
                    TNet.edges_right[0:rank, nRight] = legs[0:rank]
                TNet.tlink[G.size-1] = TNet.tlink_right[nRight].shallow_copy()
                
                #common_util.ishift(ord2, 1, -1)

                #ord1.insert(2,G.size-1)
                ord2.pop(0)
                ord2.pop(0)            
                ord2.insert(0, G.size-1)
                #print "ggg2", G
            
            if info-1>0:
                print "right side is completed contracted\n"
        
        #if 1:
        
        for ii, sub_list in enumerate(order_lists): 
            #self.edges_lists = np.ndarray((Graphics.MaxEdge, Graphics.MaxNode), dtype=int)
            sub_edge_list = self.edges_lists[ii]
            num_tensor = 0  # num of tensor contracted 
            TNet.num_tensor = 0
            if sub_list != []:  
                sub_edge_list[:,0] = G.edges[:, sub_list[0]]
                #把leg编号 从 graph 赋值给 TensorNetwork
            while len(sub_list)>1:   # contract until only one tensor left
                num_tensor = num_tensor+1  # 左边tensor的数目
                G.contract_nodes(sub_list[0], sub_list[1])   # 这里G.size  加了1
                if num_tensor  >  TNet.num_tensor:  
                    TNet.num_tensor = num_tensor
                    #合并出的tensor的序号为 num_tensor
                    v1 = sub_edge_list[:, num_tensor-1]
                    v2 = G.edges[:, sub_list[1]]
                    temp, legs= TNet.tlink[sub_list[0]].contract(TNet.tlink[sub_list[1]], V1=v1, V2=v2, out_Vc=False, use_buf= use_buf, info=info-1)
                    TNet.tlink_left[num_tensor]=temp 
                    rank = TNet.tlink_left[num_tensor].rank
                    sub_edge_list[0:rank, num_tensor] = legs[0:rank]
                
                TNet.tlink[G.size-1] = TNet.tlink_left[num_tensor].shallow_copy()
            
                #sub_list.insert(2,G.size-1)
                #print "ggg1", G
                sub_list.pop(0)
                sub_list.pop(0)
                sub_list.insert(0, G.size-1)
            
                 
            if info-1>0:
                print "left side is completed contracted\n"
            print_vars(vars(),  ['ii', 'aaaaaaaaaaaaaaaaa', 'sub_list']) 
        #把左右两边点收缩掉
        raise 
        if ord1 and ord2:  
            G.contract_nodes( ord1[0], ord2[0])
            #print "ttt ", type(TNet.tlink[ord1[0]]), ord1[0]
            v1 = TNet.edges_left[:, nLeft]
            v2 = TNet.edges_right[:, nRight]

            TNet.tlink[G.size-1], legs = TNet.tlink[ord1[0]].contract(TNet.tlink[ord2[0]], v1, v2, 
                    out_Vc=False, use_buf= use_buf, info=info-1)
            rank = TNet.tlink[G.size-1].rank
            #TNet.graph.edges[0:rank, G.size] = legs[0:rank]
        
        #rearange 即把最终的张量重新排一下
        need_rearrange =  False 
        P = None 
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
            Tout=TNet.tlink[G.size-1].permutation(P, use_buf=use_buf)
        else:
            Tout = TNet.tlink[G.size-1]
        if info>0:
            print "\nentire tensor net is contracted, yielding:"
            #print "type_name: ", Tout.type_name, "\trank: ", Tout.rank, "\tleg label:", P[:Tout.rank] #G.edges[:Tout.rank, G.size-1]
            print "type_name:{0.type_name}\trank:{0.rank}\tlebel{label}\tdims:{0.Dims}".format(
                    Tout, label=P)
            print Tout.data[:4].round(5),"...", Tout.data[-4:].round(5)
        
        return Tout

    def contract_labeled_tensor_list(self, tlist): 
        raise NotImplemented('todo') 
    
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
    
class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass 
        tn=TensorNetwork()
        if 1: 
            g3=Graphics()
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
        self.g3 = g3         
        if 1: 
            g5=Graphics()
            g5.size = 5
            g5.nodes[:]= -1
            g5.edges[:, :] = -1 

            #g5.nodes[0] = 3
            g5.nodes[:]=2
            if 1: 
                g5.edges[:g5.nodes[0], 0] = [4, 0]
                g5.edges[:g5.nodes[1], 1] = [0, 1]
                g5.edges[:g5.nodes[2], 2] = [1, 2]
                g5.edges[:g5.nodes[3], 3] = [2, 3]
                g5.edges[:g5.nodes[4], 4] = [3, 4]
        self.g5 = g5
        

        if 1: 
            g10=Graphics()
            g10.size = 10
            g10.nodes[:]= -1
            g10.edges[:, :] = -1 
            #g10.nodes[0] = 3
            g10.nodes[:]=2
            if 1: 
                g10.edges[:g10.nodes[0], 0] = [9, 0]
                for i in range(1, 10): 
                    g10.edges[:g10.nodes[i], i] = [i-1, i]
                
        self.g10 = g10
            
    def test_contract_tensor_list(self):
        if 1: 
            g3 = self.g3.copy()
            t222a = iTensor.example(rank=3, symmetry='Travial', rand_seed=1234)
            t222b = iTensor.example(rank=3, symmetry='Travial', rand_seed=1234)
            t2222 = iTensor.example(rank=4, symmetry='Travial', rand_seed=1234)
            
            T = [t222a,t222b,t2222]
            V= g3.edges
            V = []
            for i in range(g3.size): 
                V.append(g3.edges[: g3.nodes[i], i])
            print_vars(vars(),  ['g3', 'V'])
            #ord=[0, 2, 1]
            #nord=len(ord)
            ford=[5,6,7]
            
            to, label = TensorNetwork.contract_tensor_list(T, V,ford)
            print 
            self.assertAlmostEqual(to.data[0], 0.84861864605708592, 10)
            
    def test_contract_except(self):
        #g5 = cls.g5
        if 1: 
            g5 = self.g5.copy()
            tn = TensorNetwork()
            tn.TG=g5
            n=g5.size
            #order=[0,2,4,5,1]
            order=[0,1, 2, 3, 4]
            exception= 2
            t22=[]
            for i in range(n):
                a = iTensor.example(rank=2, symmetry='Travial', rand_seed=i)
                t22.append(a.copy())
            
            tn.tlink[:n]=t22[:n]
            to = tn.contract_except(n, order, exception, info=1)
            self.assertAlmostEqual(to.data[10],-0.0064333958499588409, 10)

        if 1: 
            g5 = self.g5.copy()
            tn = TensorNetwork()
            tn.TG=g5
            n=g5.size
            #order=[0,2,4,5,3]
            order=[0,1, 2, 3, 4]
            exception= 2
            
            t22=[]
            for i in range(n):
                a = iTensor.example(rank=2, symmetry='U1', rand_seed=i)
                a.type_name = str(i)
                t22.append(a.copy())
            
            tn.tlink[:n]=t22[:n]
            to = tn.contract_except(n, order, exception, info=5)
            self.assertAlmostEqual(to.data[5], 0.0066916824757841942, 10)

    def test_contract_except_new(self):
        #g5 = cls.g5
        if 1: 
            g10 = self.g10.copy()
            tn = TensorNetwork()
            tn.TG=g10
            n=g10.size
            #order=[0,2,4,5,1]
            #order=[0, 1, 2, 3, 4]
            order = range(10)
            exception= 2
            t22=[]
            for i in range(n):
                a = iTensor.example(rank=2, symmetry='Travial', rand_seed=i)
                a.type_name = str(i)
                #a = iTensor.example(rank=2, symmetry='U1', rand_seed=i)
                t22.append(a.copy())
            
            tn.tlink[:n]=t22[:n]
            exception = [2, 5]
            #to = tn.contract_with_removals(n, order, exception, info=0)
            to = tn.contract_except_test(n, order, exception, info=1)
            #print_vars(vars(),  ['to.data'])

    def test_temp(self):
        pass 
        if 1: 
            g5 = self.g5.copy()
            tn = TensorNetwork()
            tn.TG=g5
            n=g5.size
            #order=[0,2,4,5,1]
            order=[0, 1, 2, 3, 4]
            exception= 2
            t22=[]
            for i in range(n):
                a = iTensor.example(rank=2, symmetry='Travial', rand_seed=i)
                a.type_name = str(i)
                t22.append(a.copy())
            
            tn.tlink[:n]=t22[:n]
            to = tn.contract_except_test(n, order, exception, info=2)
            #to = tn.contract_except(n, order, exception, info=2)
            print_vars(vars(),  ['to.data'])
            self.assertAlmostEqual(to.data[10],-0.0064333958499588409, 10)


if __name__ == "__main__":


    if 0: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_temp', 
           #'test_contract_tensor_list', 
           #'test_contract_except', 
           'test_contract_except_new', 
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)





        

