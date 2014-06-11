#! /usr/bin/env python
#coding=UTF8
import unittest

from collections import OrderedDict
from math import log
#from Set import Set

""" 
    the graph(or network) used in TNS algorithms are by definition directed multigraphs.
"""
import numpy as np

GRAPH_MODULE_NAME_ALL = ['binary', 'modified_binary', 'ternary', 'quaternary', 'quinary', 'septenary']

class Graphics():
    MaxEdge = 64
    MaxNode = 64

    #def __init__(self,names=[""]*MaxNode, size=0, nodes=np.empty(MaxNode, "int"), edges=np.empty((MaxEdge,MaxNode),"int")):
    #VERY STRANGE, above init cause serious problem
    def __init__(self, size=0):
        """ see init_Graphics in f90
            Names：a list of name str
            size：节点数
            nodes: 表示这个节点有几条腿, -1, 则说明该节点不存在
            edges[:, n1]: Nodes n1的连接情况, 例如，edges[:G.nodes[3], 3] = [1, 3, 5, -2]
            表示连接在3号节点的所有边的编号为1，3，5，-2(负数is ok)
        
        """
        self.graph_name =""
        self.names= [""]*self.MaxNode
        self.size = size   
        self.nodes = np.empty(self.MaxNode, "int")
        self.nodes[:] = -1
        self.edges = np.empty((self.MaxEdge,self.MaxNode),"int")
        self.edges[:, :]=-1

        #for contraction of tensor net
        #self.contract_order = None # np.ndarray(Graphics.MaxNode, "int")        
        self.contract_order = [None]*size
        self.final_order = None  #arange final order of legs
        self.weight = None
        self.comp_complexity = None # todo: define this
    
    @staticmethod
    def init_from_dict(graph_def, pos_offset=1): 
        """
            graph_def is something e.g.
            G_calc_EE = {
                'V_0_0':  {1: ('Vp_0_0', 2), 2: (None, None), 3: (None, None),  4: ('V_1_0', 2)  } , 
                'V_0_1':  {1: (None, None),  2: (None, None), 3: ('Vp_0_1', 4), 4: ('V_1_1', 2)} ,

                'V_1_0':  {1: ('Vp_1_0', 2), 2: ('V_0_0', 4), 3: (None, None),  4: ('rho', 1)}, 
                'V_1_1':  {1: (None, None),  2: ('V_0_1', 4), 3: ('Vp_1_1', 4), 4: ('rho', 2)},
                    
                'Vp_0_0': {1: ('Vp_1_0', 3),                  3: (None, None),  4: (None, None)}, 
                'Vp_0_1': {1: ('Vp_1_1', 3), 2: (None, None), 3: (None, None)                   }, 
                
                'Vp_1_0': {1: ('rho', 3), 4: (None, None)}, 
                'Vp_1_1': {1: ('rho', 4), 2: (None, None)}
            }
        
        """
        #graph_def = sorted(graph_def)
        all_node= []
        all_edge = []
        #for k, val in graph_def.items(): 
        for k in sorted(graph_def): 
            val = graph_def[k]
            n1 = k
            if k not in all_node: 
                all_node.append(k)
            for pos1, v in val.items(): 
                n2 = v[0]
                pos2 = v[1]
                if v[0] not in all_node: 
                    all_node.append(v[0])
                if (n1, pos1, n2, pos2) not in all_edge and (n2, pos2, n1, pos1) not in all_edge: 
                    all_edge.append((n1, pos1, n2, pos2))
                    
        #print all_node
        if None in all_node: all_node.remove(None)
        
        #print all_edge
        
        size = len(all_node)
        res = Graphics(size=size)
        res.nodes[: ] = -1
        res.edges[:, :] = -1
        
        #edge_id = 0
        edge_id = -1 
        for  n1, pos1, n2, pos2 in all_edge:
            edge_id  += 1  
            for node, pos in [(n1, pos1), (n2, pos2)]: 
                if node is not None: 
                    node_id = all_node.index(node)
                    res.names[node_id] = node
                    if res.nodes[node_id] == -1: 
                        res.nodes[node_id] = 1
                    else:
                        res.nodes[node_id] += 1 
                    if 0: 
                        if node == 'V_1_1':  
                            print n1, pos1, n2, pos2, '\t','node_id', node_id, 'edge_id', edge_id, 
                            print 'num nodes', res.nodes[node_id]
                    res.edges[pos-pos_offset][node_id] = edge_id
                    #res.edges[pos][node_id] = edge_id
        return res
    
    from_dict = init_from_dict
    
    def set_contract_order(self, order):
        if isinstance(order[0], str): 
            self.contract_order = [self.find_node(n) for n in order]
        elif isinstance(order[0], int): 
            self.contract_order = order
    
    def set_final_order(self, order, ind_offset=1): 
        """
            case 1: order is like this: 
                [('V_0_0', 1), ('V_0_0', 2), ('U_0_0', 1), ('U_0_0', 2), ('V_0_1', 2), ('V_0_1', 3), 
                , ('V_0_0', 4), ('V_0_1', 4)]

        """
        if len(order[0])>1: 
            final_order = []
            for nname, pos in order: 
                node_id = self.find_node(nname)
                e=self.edges[pos-ind_offset, node_id][0]
                final_order.append(e)
        else: 
            raise
        self.final_order = final_order
    
    def copy(self):
        """
        see Graphics_Assignment in f90
        """
        other = Graphics()
        size = self.size
        other.size = size
        other.nodes= self.nodes.copy()
        other.edges= self.edges.copy()
        other.names[:size] = self.names[:size]
        other.weight = self.weight
        #other.contract_order = self.contract_order.copy()
        other.contract_order = [i for i in self.contract_order]
        return other

    def __copy__(self):
        return self.copy()

    def find_node(self, name):
        """
            1 name may corresponding to n nodes
        """
        
        res= [i for i in range(self.size) if name == self.names[i]]
        return res

    def __repr__old(self):
        """
        status_1_verified
        """
        from utilities import repr_format
        dic = self.__dict__
        keys= ["size", "nodes", "edges"]
        res= ""
        for k in keys:
            value_str = str(dic[k])
            if k == "nodes":
                value_str = str(dic[k][:self.size])
            if k == "edges":
                value_str = ""
                for i in range(self.size):
                    if self.nodes[i]>0:
                        value_str  += "%5d\t"%i + str(dic[k][:self.nodes[i], i]) + "\n\t"
            res+=  repr_format(k, value_str)
        return res

    def __repr__(self):
        """
        """
        #from utilities import repr_format
        dic = self.__dict__
        keys= ["size", "nodes", "edges", "contract_order", "final_order", "weight"]
        res= "\n"
        for k in keys:
            temp = str(dic[k])
            if k == "contract_order":
                temp = str(dic[k][:self.size])
            if k == "nodes":
                temp = str(dic[k][:self.size])
            if k == "edges":
                temp = ""
                for i in range(self.size):
                    if self.nodes[i]>0:
                        temp += "%5d\t"%i  + self.names[i] + "\t"+ str(dic[k][:self.nodes[i], i]) + "\n"
            res+=  k + ":\n\t" + temp + "\n"
        return res

    def find_common_leg(self,n1,n2):
        """
        status_1_verified
        see FindCommonLeg in f90
        """
        a=self.edges[:self.nodes[n1], n1].tolist()
        b = self.edges[:self.nodes[n2], n2].tolist()
        anb = [i for i in a if i in b]
        ab = [i for i in a if i not in anb ] + [i for i in b if i not in anb]
        ne = len(anb)
        #edges= np.array([-1L]*self.MaxEdge, "int")
        return anb

    def contract_nodes(self, n1,n2 ):
        """
        status_1_verified
        see Node_Contraction in f90
        Contract all edges between nodes n1 and n2 in G   
        optional argument ne, edges denote the number and id of edges 
        that have been contracted
        """
        a = self.edges[:self.nodes[n1], n1].tolist()
        b = self.edges[:self.nodes[n2], n2].tolist()
        anb = [i for i in a if i in b]
        ab = [i for i in a if i not in anb ] + [i for i in b if i not in anb]
        ne = len(anb)
        edges= anb

        self.nodes[n1] = -1
        self.nodes[n2] = -1

        self.size = self.size + 1
        #print self.size,  len(ab)
        self.nodes[self.size-1] = len(ab)
        #print ab
        if len(ab)>0:
            self.edges[:len(ab), self.size-1] = ab[:]

    def contract_nodes_new(self, n1,n2 ):
        """
        这个函数有问题，edges不好处理
        see Node_Contraction in f90
        Contract all edges between nodes n1 and n2 in G   
        optional argument ne, edges denote the number and id of edges 
        that have been contracted
        """
        a = self.edges[:self.nodes[n1], n1].tolist()
        b = self.edges[:self.nodes[n2], n2].tolist()
        anb = [i for i in a if i in b]
        ab = [i for i in a if i not in anb ] + [i for i in b if i not in anb]
        ne = len(anb)
        edges= anb

        #self.nodes[n1] = -1
        #self.nodes[n2] = -1
        self.nodes.pop(n1)
        if n1<n2:
            self.nodes.pop(n2-1)
        else:
            self.nodes.pop(n2)
        self.size = self.size -1

        self.nodes[self.size-1] = len(ab)
        #print ab
        if len(ab)>0:
            self.edges[:len(ab), self.size-1] = ab[:]

    def remove_one_node(self, node_name):
        """
            not completed
        """
        node = self.find_node(node_name)[0]
        self.nodes[node] = -1
        self.node_name[node] = "null"
    
    def to_dict(self): 
        res = {}
        for node_id in range(self.size):  
            node_name = self.names[node_id]
            edge_dic = {}
            nleg = self.nodes[node_id]
            for l in range(nleg): 
                
                leg_id = self.edges[l, node_id]
                temp=np.where(self.edges==leg_id)
                if 0: 
                    if node_name == 'U_0_1':  
                        print  l,leg_id, temp
                        if l == 3: exit() 
                    
                if len(temp[0])==1: 
                    
                    edge_dic[l] = (None, None)
                else: 
                    for i in [0, 1]: 
                        if temp[1][i] != node_id: 
                            n2 = temp[1][i]
                            name2 = self.names[n2]
                            pos2 = temp[0][i]
                            edge_dic[l] = (name2, pos2)
                
            res[node_name] = edge_dic
        return res
        
    def to_diagram(self): 
        from merapy.diagrams.diagram import Diagram
        dic_def = self.to_dict()
        G = {}
        for n in sorted(dic_def): 
            conn = dic_def[n]
            val = [(conn[i][0], conn[i][1]) for i in sorted(conn)]
            G[n] = val
        G = Diagram(G)
        return G
            
        
    def to_networkx(self): 
        pass
    

if 0:
    class TensorGraph():
        pass

def calc_ascending_site(n, n_min=1, lay=None): 
    """
    
        a site is at n, where is it at coarsed layer lay?
        Howto: 
        1. 首先存在一个系统的标记的问题, 和约定的问题 
            notation: 
                包括site的编号: (l, n) 在第l层的第n个位置
                和V tensor的编号: [l, m] 在第l层的第m个位置
            convention: 
            1) 在第0层site，第0个site，接在编号为0的Vtensor的中间的leg上
                由此对于(l, n) , 接在 [l, n//3] 上，if n%3==0：中间腿；n%3==1: 右边腿； n%3==2: 左边腿
            2) l层的编号为0的Vtensor朝上的腿，标记为第l+1层的0th site, 这样一次类推
        2. 按照上述convention，再结合mera的结构，即可。具体规律为: 
            1), 先把(l, n)，变成2site的(l, (n-1, n))
            2)  (l, (n-1, n)) 被映射到 (l+1,  )
            
        
        3.例子：
            e.g.
                    (0, 40) =>(0, 39-40)    
                                首先预处理一下，把1site变成 2site
                                40 is on the right leg of 13th V tensor
                --> (0, 13-14), 13 is on the right leg of 4th V tensor, 14 is on the left leg of 5th V tensor
                --> (0, 4-5),   4  is on the right leg of 1th V tensor, 5  is on the left leg of 2th V tensor
                --> (0-1-2)
                故 总的来说，要判断在哪(两个)个V上
    """
        
    tube = [(n-1, n)]  #auto complete to two-site ops
    
    while 0: 
        print 'hhhhh'
        n_mod3 = n%3
        n_div3 = n//3
       
        if   n_mod3 == 0: 
            res = (n_div3-1, n_div3)
        elif n_mod3 == 1: 
            res = (n_div3,   n_div3 + 1)
        elif n_mod3 == 2: 
            res = (n_div3,   n_div3 + 1)
        n = res[-1]   
        
        tube.append(res)    
        if n <=  n_min: 
            break
    
    #(4, 5)->(1, 2), (3, 4)->(1, 2), (2, 3)->(0, 1), (1, 2)->(0, 1)
    while 1: 
        if n <=  n_min: 
            break
       
        n_mod3 = n%3
        n_div3 = n//3
       
        if   n_mod3 == 0: 
            res = (n_div3-1, n_div3)
        elif n_mod3 == 1: 
            res = (n_div3,   n_div3 + 1)
        elif n_mod3 == 2: 
            res = (n_div3,   n_div3 + 1)
        n = res[-1]   
        if n >= n_min: 
            tube.append(res)    


    return tube


def coarse_grain_map(X, tau, n=1, T=None, return_full_list=0):
    """
        n: act on X by R for n times
        T: total layer  of mera
    """
    #if T is None: 
    T =log(max(X), 3)
    base = 3
    causal_cone = []
    for i in range(n): 
        X_new = []
        if tau + i>= T: 
            causal_cone.append([0])
            break 
        for x in X:
            if x%base ==  0 : 
                x1 = [x/3]
            else: 
                x1 = [x//3, x//3  + 1]
            for a in x1: 
                if a not in X_new: 
                    X_new.append(a)
        #print X, X_new
        X = list(X_new)
        causal_cone.append(tuple(X_new))
    #print causal_cone
    if return_full_list: 
        return causal_cone
    else: 
        return causal_cone[-1]

def fine_grain_map(X, tau, n): 
    """
        n: num acted by fine map 
    """
    temp = [(i, i+1) for i in range(1000)]
    res= [i for i in temp if coarse_grain_map(i, 0, n)==X]
    #res = calc_ascending_site(27, n_min=1, lay=1 )
    print res
    #res= []  
   

def calc_ascending_site_binary(n, n_min=1, lay=None): 
    """
    
        a site is at n, where is it at coarsed layer lay?
        Howto: 
        1. 首先存在一个系统的标记的问题, 和约定的问题 
            notation: 
                包括site的编号: (l, n) 在第l层的第n个位置
                和V tensor的编号: [l, m] 在第l层的第m个位置
            convention: 
            1) 在第0层site，第0个site，接在编号为0的Vtensor的中间的leg上
                由此对于(l, n) , 接在 [l, n//3] 上，if n%3==0：中间腿；n%3==1: 右边腿； n%3==2: 左边腿
            2) l层的编号为0的Vtensor朝上的腿，标记为第l+1层的0th site, 这样一次类推
        2. 按照上述convention，再结合mera的结构，即可。具体规律为: 
            1), 先把(l, n)，变成2site的(l, (n-1, n))
            2)  (l, (n-1, n)) 被映射到 (l+1,  )
            
        
        3.例子：
            e.g.
                    (0, 40) =>(0, 39-40)    
                                首先预处理一下，把1site变成 2site
                                40 is on the right leg of 13th V tensor
                --> (0, 13-14), 13 is on the right leg of 4th V tensor, 14 is on the left leg of 5th V tensor
                --> (0, 4-5),   4  is on the right leg of 1th V tensor, 5  is on the left leg of 2th V tensor
                --> (0-1-2)
                故 总的来说，要判断在哪(两个)个V上
    """
    raise NotImplemented('not completed')    
    tube = [(n-1, n)]  #auto complete to two-site ops
    
    #rescale_factor = 2
    b = 2  # rescale_factor
    
    #(4, 5)->(1, 2), (3, 4)->(1, 2), (2, 3)->(0, 1), (1, 2)->(0, 1)
    while 1: 
        if n <=  n_min: 
            break
       
        n_mod_b = n%b
        n_div_b = n//b
       
        if n_mod_b == 0: 
            res = (n_div_b-1, n_div_b)
        elif n_mod_b == 1: 
            res = (n_div_b,   n_div_b + 1)
        elif n_mod_b == 2: 
            res = (n_div_b,   n_div_b + 1)
        n = res[-1]   
        if n >= n_min: 
            tube.append(res)    


    return tube


def mera_full_graph(): 
    site = OrderedDict()
    #site at each layer, number
    num_of_layer = 6
    N0 = 2*3**(num_of_layer/2)
    Nlay = OrderedDict()
    node_U = OrderedDict()
    node_V = OrderedDict()
    
    Nlay[0] = N0
    for l in range(1, num_of_layer/2): 
        Nlay[l] = Nlay[l-1]/3
        for n in range(Nlay[l]): 
            #site[l, n] = {'in':     , 'out': {'arrow': ' ', 'pos':    }}
            n_mod3 = n%3 
            l_mod2 = l%2
            out_arrow_name = 'V' if n_mod3 == 0 else 'U'
            in_arrow_name  = 'V' if n_mod3 == 0 else 'U' 
            site[l, n] = {
                    'in_arrow':  {'name':  in_arrow_name, 'lay': l-1, 'pos': (n-n_mod3)/3}, 
                    'out_arrow': {'name': out_arrow_name, 'lay': l+1, 'pos': (n-n_mod3)/3}
                    }
    print site   

class SupportSet(object): 
    """  
        an ordered set 
        including effective support set, 
        denoted by [X;tau]
    
    """

class SupportPatern(object): 
    def __init__(self, N): 
        #assert 
        self.lattice =  tuple(range(N))
        
    def eff_lattice(self, level): 
        pass

class TestGraphics(unittest.TestCase): 
    def setUp(self): 
        pass
    def test_aaa(self): 
        pass
    def test_coarse_grain_map(self): 
        X = [55]; tau = 0; n = 1
        res=coarse_grain_map(X, tau, n)
        self.assertTrue(res==(18, 19)) 
        
   

#mera_full_graph(); exit()    



if __name__ == '__main__':        
    

    if 1: 
        g=Graphics()
        def set_value():
            g.size=8
            g.nodes[:]=-1
            g.edges[:, :]=-1
            g.names[1]="OO" #OO
            g.nodes[1]=4
            g.edges[0:4,1]=[1,2,3,4]
            g.names[2]="Up" #Up_1
            g.nodes[2]=4
            g.edges[0:4,2]=[5,6,7,8]
            g.names[3]="U" #U_1
            g.nodes[3]=4
            g.edges[0:4,3]=[7,9,10,11]
            g.names[4]="Vp" #Vp_0
            g.nodes[4]=4
            g.edges[0:4,4]=[1,12,13,5]
            g.names[5]="Vp" #Vp_1
            g.nodes[5]=4
            g.edges[0:4,5]=[2,6,14,15]
            g.names[6]="oo" #oo_3_4
            g.nodes[6]=4
            g.edges[0:4,6]=[8,14,9,16]
            g.names[7]="V" #V_1
            g.nodes[7]=4
            g.edges[0:4,7]=[11,16,15,4]
            g.names[8]="V" #V_0
            g.nodes[8]=4
            g.edges[0:4,8]=[12,13,10,3]
        #set_value()
        def set_value1():
            g.size = 3
            g.nodes[:]= -1
            g.edges[:, :] = -1 

            g.nodes[0] = 3
            g.nodes[1] = 3
            g.nodes[2] = 2        

            g.edges[:g.nodes[0], 0] = [0, 2, 3]
            g.edges[:g.nodes[1], 1] = [1, 0, 3]
            g.edges[:g.nodes[2], 2] = [2, 1]        
        set_value1()

        def test_find_common_leg():
            print g.find_common_leg(0, 1)
            print g.find_common_leg(0, 2)        
        #test_find_common_leg()

        def test_node_contract():
            print g
            g.contract_nodes(0, 1)
            print g
        #test_node_contract()


        #print g
        g1 = Graphics()
        g1.size = 4
        g2 = Graphics()
        g2.size = 2
        print g1.size
        print g    

    if 1: #examine
        suite = unittest.TestLoader().loadTestsFromTestCase(TestGraphics)
        unittest.TextTestRunner(verbosity=0).run(suite)    
       
   
    
    #fine_grain_map()    
    #print  all([1, 2]==(1, 2))
    X = (3, 4); tau = 0; n = 1
    fine_grain_map(X, tau, n)
    
    
    
    

