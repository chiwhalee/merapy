#! /usr/bin/env python
#coding=UTF8

"""
see mera.f90
q:
    q299, why svd
    q300,  what is unit_isometry
    

"""
import unittest
import numpy as np
import warnings
from math import log

from quantum_number import  * #QN_idendity, QSp_base, QuantSpace, QSp_null
from graphics import Graphics
from tensor import *
from tensor_svd import *
#from random_64 import mrandom
import crandom
from tensor_reflection import TensorReflect


__all__ = ["Mera"]

class UTensor(object):
    """
    this class may be deprecated, use list instead
    """
    def __init__(self):
        self.ilayer= None
        self.nSite= None
        self.tensor= None  #[Tensor()]
    
    def __getitem__(self, i):
        return self.tensor[i]
    def __setitem__(self, i, x):
        self.tensor[i] = x

    def shallow_copy(self):
        pass
    def strorage(self):
        pass
    def restore(self):
        pass
    def unit_disentangler(self):
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

def coarse_grain_map(X, tau=0, n=1, T=None, return_full_list=0):
    """ 
        params: 
            X:  tuple of list  
            n: act on X by R for n times
            T: total layer  of mera
    """
    if T is None: 
        #print_vars(vars(),  ['X'])
        #T =log(max(X), 3)
        T = 5 
    base = 3
    causal_cone = []
    for i in range(n): 
        X_new = []
        if tau + i>= T:  #already reach top  
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
    def __init__(self): 
        self.connectivity = None   # num of connected intervals can be divided into 
    
    @staticmethod 
    def split_disconnect(X): 
        """
            X:  must be a list of integers 
        """
        
        res = []
        iprev = X[0]
        temp = (iprev, )
        for i in X[1: ]:
            if i == iprev + 1: 
                temp  = temp + (i, )
            else: 
                res.append(temp)
                iprev = i
                temp = (iprev, )
            iprev = i 
        res.append(temp)
        res= tuple(res)
        if len(res)==1:  #make e.g. ((0, 1), ) into (0, 1)
            res= res[0]
        #print_vars(vars(),  ['X', 'res'])

        return res 
    
    @staticmethod 
    def join_disconnect(X): 
        res= []
        for i in X: 
            if isinstance(i, int): 
                res.append(i)
            else: 
                res.extend(list(i))
        return res 

class SupportPaternEff(object): 
    """
        effective support pattern 
    """
    def __init__(self, N): 
        #assert 
        self.lattice =  tuple(range(N))
        
    def eff_lattice(self, level): 
        pass
    
    @staticmethod 
    def show(spe): 
        s= '**** Support Patern ***\n'
        for i in spe: 
            p = spe[i]
            if len(p)<8: 
                s += 'l={} \t{}\n'.format(i, p)
            else: 
                s += 'l={} \t{} \n\t... {}\n'.format(i, p[:4], p[-4: ])
        print s 
        
SPE = SupportPaternEff 

    

class Mera(object):
    """
        Multi-scale entanglement renornalization ansatz representing the wave function
    """
    MaxLayer = 16

    def __init__(self, layer, nTop, qsp_0, qsp_max, qsp_max2, 
            topQN, qsp_null, qn_identity, tensor_defs=None, 
            rand_seed=None, USE_CUSTOM_RAND=False, USE_REFLECTION=False,   
            unitary_init=None):
        """
            see init_Mera in f90
            给所有的tensor allocate space

            layer是总共有多少层, 
            ilayer是层数的标志
            nSite是有多少个site
            nNeigh有多少个邻居，应该是已经不用了
            
            qsp_0是最底层的quantum space,即U tensor进来腿的维数
            Qsp_max/qsp_max2是最大的Quantum space
            Qsp_max是U tensor出去的腿的Quantum space
            qsp_max2是V tensor的
        """
        #self.nNeigh = nNeigh
        #self.NSize = nNeigh**layer
        self.nTop = nTop  #this param is not used at all 
        self.qsp_0 = qsp_0
        self.qsp_max = qsp_max #.copy()
        self.qsp_max2 = qsp_max2 #.copy()
        self.topQN = topQN
        self.qsp_null = qsp_null
        self.qn_identity = qn_identity 
        self.symmetry = self.qsp_0.QNs[0].SYMMETRY 
        self.USE_REFLECTION = USE_REFLECTION 

        self.tensor_defs= {
                "U":        {"type":(2, 2),    "type_1":"unitary",	        "dual":False,}, 
                "U_dag":    {"type":(2, 2),    "type_1":"unitary",	        "dual":True, }, 
                "V":        {"type":(3, 1),    "type_1":"disentangler",	    "dual":False,}, 
                "V_dag":    {"type":(1, 3),    "type_1":"disentanger",	    "dual":True, }, 
                }
        
        if tensor_defs is not None:
            for k, v in tensor_defs.iteritems():
                if self.tensor_defs.has_key(k):
                    self.tensor_defs[k].update(v)
                else:
                    self.tensor_defs[k] = v

        self.period = self.tensor_defs["V"]["type"][0]  #this can essecially distinguish ternary, quaternary etc

        #self.tensor_defs= ["U", "V", "U_dag", "V_dag"]
        for i in self.tensor_defs:
            setattr(self, i, [])

        #给每一层self.U[ilayer][:] allocate space
        self.num_of_layer = 0
        for ilayer in range(layer-1):  #init all layers except the top
            self.add_one_layer(unitary_init=unitary_init)
        
        self.add_top_layer(rand_init=True)
        
        self.num_of_input = 2 * 3**self.num_of_layer 
        #self.layer_eff = tuple(range(layer))
        self.layer_keys = tuple(range(layer))
        self.top_layer_key = layer - 1
        
        #if USE_CUSTOM_RAND:
        #    crandom.mrandom.csrand(rand_seed)
        #else:
        #    np.random.seed(rand_seed)
        #    crandom.rand = np.random.random
        

    def __eq__(self, other):
        """
            test two meras FORMALLY equal
        """
        res= True
        temp = ["num_of_layer", "qsp_0", "qsp_max", "topQN", "period"]
        for i in temp:
            a = self.__getattribute__(i)
            b = other.__getattribute__(i)
            if  a != b :
                print i, a, b
                res = False
                break
        return res
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    @classmethod 
    def example(cls, trunc_dim=2, tot_layer=4, symmetry="Z2", combine_2site=False):
        """
        
        """
        from quantum_number import init_System_QSp

        #from quantum_number import QSp_base, QN_idendity
        from tensor_network import TensorNetwork

        QN_idendity, QSp_base , QSp_null = init_System_QSp(symmetry=symmetry)


        #tensor.init_Buffer()

        qsp_0= QSp_base.copy()
        qsp_0.update()

        
        if symmetry == "Z2" :
            qsp_max = qsp_0.copy()
            print_vars(vars(),  ['trunc_dim'])
            qsp_max.Dims[0:2] = trunc_dim/2
            qsp_max.RefQN[0,0:2] = [1,1]
            qsp_max.RefQN[1,0:2] = [0,0]
            qsp_max.update()


            qsp_max2 = qsp_0.copy()
            qsp_max2.Dims[0:2] = 2
            qsp_max2.RefQN[0,0:2] = [2,1]
            qsp_max2.RefQN[1,0:2] = [0,1]
            qsp_max2.update()
        elif symmetry == "U1":  
            qsp_max =QspU1.max(trunc_dim)
            qsp_max2 = qsp_max.copy()
            
        elif symmetry == "Travial" :
            qsp_max= QSp_base.max(trunc_dim)
            qsp_max2 = qsp_max.copy()
        elif symmetry == "Z3":
            qsp_max= QSp_base.max(trunc_dim)
            qsp_max2 = qsp_max.copy()

        else:
            raise
        
        topQN= QN_idendity.copy()
        
        crandom.mrandom.csrand(1234)
        nTop = 2
        M= Mera(tot_layer, nTop, qsp_0, qsp_max, qsp_max2, topQN, QSp_null, QN_idendity)
        #sys.init_Hamiltonian(M)
        return M
   
    def renormalization_map(self, i): 
        raise NotImplemented
    coarse_grain_map = renormalization_map
    
    def calc_support_patten(self, support_patten_phys): 
        raise NotImplemented 
    
    def calc_asd_ops(self, support_patten_lay): 
        raise NotImplemented 
    
    def key_property(self):
        qsp_max_property = self.qsp_max.key_property()
        trunc_dim = self.qsp_max.totDim
        assert trunc_dim == sum(qsp_max_property['dims']) 
        res= {'num_of_layer':self.num_of_layer, 'trunc_dim':trunc_dim}
        res.update(qsp_max_property)
        return res


    def add_one_layer(self, unitary_init=None, rand_init=True):
        """
            see Mera_Addlayer in f90
            if true [default] randominze initial value, otherwise 
            use value of lower layer        
        """
        #ilayer = self.num_of_layer+1
        ilayer = self.num_of_layer
        self.init_layer(ilayer, unitary_init=unitary_init, rand_init=rand_init) #, rand_init=rand_init)
        #总层数增加
        self.num_of_layer = self.num_of_layer+1
    
    def init_layer(self, ilayer, unitary_init=None, rand_init=True):
        """
            assign EVERY attributes: QSp, data...  for tensors in the layer
            对U，V 进行实质性的赋初值
            rand_init:
                if true [default] randominze initial value, otherwise 
                use value of lower layer        
        """

        if unitary_init is None:
            if self.symmetry  in ["U1"]: 
                unitary_init = "random_unit_tensor"
            elif self.symmetry in ["Z2", "Z3", "Travial"]:
                unitary_init = "unit_tensor"
            else:
                raise
        self.unitary_init = unitary_init

        #这里的设置决定了mera中张量的腿的维数
        qsp_0 = self.qsp_0   #.copy() 
        qsp_max = self.qsp_max  
        qsp_max2 = self.qsp_max2  
        
        j=0
        QSp_U = [None]*4
        QSp_V = [None]*4
        
        if ilayer == 0:  
            qsp = qsp_0.copy()
        else:
            # U quantspace 接前一层V的
            temp = self.V[ilayer-1][0]
            qsp = temp.QSp[temp.rank-1].copy()   #here must copy!
            qsp.reverse()
        
        #assign init value for legs of U and V tensors
        #U的2，3两条腿朝上
        down, up = self.tensor_defs["U"]["type"]
        QSp_U = qsp.copy_many(down) + qsp.copy_many(up, reverse=range(up))
        down = self.tensor_defs["V"]["type"][0]
        QSp_V = qsp.copy_many(down)
        
        temp0 = qsp.copy()
        temp = qsp.copy()
        for i in range(down-1):
            temp = temp.add(temp0)
            temp.update(qsp_max = qsp_max)
        #attention_this_may_be_wrong   reflect_symmetrize is not compatible with qsp.add
        #temp.update(qsp_max=qsp_max)
        temp.reverse()
        QSp_V.append(temp)
        
        #just make some place holders
        for i in self.tensor_defs:
            self.__getattribute__(i).append([None])

        #给self.U[layer][:]赋值
        #rank = 4
        rank = sum(self.tensor_defs["U"]["type"])
        totQN = self.qn_identity.copy()
        for j in range(0, 1):
            #self.U[ilayer][j] = iTensor(rank, QSp_U, totQN)
            self.U[ilayer][j] = iTensor(rank, QSp_U, totQN)
            ndiv = self.U[ilayer][j].rank//2
            self.U[ilayer][j].ndiv = ndiv
            #initialize data of tensors
            if rand_init:                  
                #seems that rand_init always true
                temp = self.U[ilayer][j]
                if unitary_init == "unit_tensor": 
                    self.U[ilayer][j] = iTensor.unit_tensor(temp.rank, temp.QSp, temp.totQN)
                elif unitary_init == "random_unit_tensor":
                    Tensor_svd.random_unit_tensor(self.U[ilayer][j],2)
                else:
                    raise ValueError("unitary_init is: %s"%unitary_init)
                #self.U[ilayer][j] = unit_tensor(temp.rank, temp.QSp, temp.totQN)
            else:
                self.U[ilayer][j]=self.U[ilayer-1][j].copy()

            self.U_dag[ilayer][j]=self.U[ilayer][j].conjugate(2 )

        #给self.V[layer][:]赋值        
        #rank = 4
        rank = sum(self.tensor_defs["V"]["type"])
        for j in range(0, 1):
            self.V[ilayer][j]=iTensor(rank, QSp_V, totQN)
            ndiv = self.V[ilayer][j].rank-1
            self.V[ilayer][j].ndiv = ndiv
            if rand_init:  
                #attention_low_efficiency   random_unit_tensor may use buffer
                Tensor_svd.random_unit_tensor(self.V[ilayer][j], ndiv)
            else:
                self.V[ilayer][j] = self.V[ilayer-1][j].copy()
            
            if self.USE_REFLECTION:
                TensorReflect.reflect_symmetrize(self.V[ilayer][j], ndiv=ndiv)
                Tensor_svd.svd(self.V[ilayer][j], ndiv=ndiv)
            
            self.V_dag[ilayer][j] = self.V[ilayer][j].conjugate(ndiv)

    def add_top_layer(self,rand_init=True):
        """
            这个 top layer V, 在Finite range mera中并没有用到
            这只是为了形式完整，或备作它用
            add the top V tensor
            最上层的V特殊一些，是个rank(2, 0) type tensor
        """
        #print "adding top layer"
        ilayer = self.num_of_layer
        self.num_of_layer += 1 
        
        #attention_this_may_be_wrong
        #实际这里init_layer 是多余的，其不良后果是多增加了一个U layer，本不存在
        self.init_layer(ilayer, rand_init=rand_init, unitary_init=self.unitary_init)
        
        QSp = [None for i in range(3)]
        
        j = 0
        #topQN = self.topQN
        ten = self.V[ilayer-1][j]
        #print "lll", len(QSp), len(self.V), ilayer-1, j, ten.rank        
        QSp[0] =ten.QSp[ten.rank-1].copy()   #下面一层V tensor超上的那条腿
        QSp[0].reverse()
        QSp[1] = QSp[0].copy()

        #QSp[2] = QSp_null.copy()   #top tensor 最上面那条腿实际上是1维的，等同于没有
        
        QSp[2] = self.qsp_null #.copy()   #top tensor 最上面那条腿实际上是1维的，等同于没有
        
        #这里实际上多余 topQN 本来就是QN_idendity
        QSp[2].QNs[0]  = self.topQN
        QSp[2].update()
        QSp[2].reverse()
        
        
        totQN = self.qn_identity.copy()  
        
        rank = 3
        self.V[ilayer][j]=iTensor(rank, QSp[:rank], totQN )

        if rand_init:  
            Tensor_svd.random_unit_tensor(self.V[ilayer][j], self.V[ilayer][j].rank-1)
        else:
            self.V[ilayer][j] = self.V[ilayer-1][j].copy()
        
        self.V_dag[ilayer][j]=self.V[ilayer][j].conjugate(self.V[ilayer][j].rank-1)
    
    def expand_layer(self, rand_init):
        """
        """
        #iLayer = M%Layer+2
        #call Mera_InitLayer(M, iLayer, rand_init=.FALSE.)
        #if (present(rand_init)) then
        #   call Mera_AddLayer(M,rand_init)!, rand_init=.FALSE.)
        #else
        #   call Mera_AddLayer(M)
        #end if

        ilayer = self.num_of_layer+2
        self.init_layer( ilayer, rand_init= False )
        if rand_init!=None:  
            self.add_one_layer(rand_init)#, rand_init= False )
        else:
            self.add_one_layer()

    def expand_dim(self, qsp_max):
        M = self
        qsp_0 = M.qsp_0
        qsp_max2 = qsp_max.copy()
        for j  in range(qsp_max2.nQN):
            qsp_max2.Dims[j] = qsp_max.Dims[j]*2
        
        M.qsp_max = qsp_max.copy()
        
        for ilayer  in range(M.num_of_layer-1):
            rank_U = sum(self.tensor_defs["U"]["type"])
            rank_V = sum(self.tensor_defs["V"]["type"])
            QSp_U = [None]*rank_U
            j = 0
            if ilayer == 0:  
                QSp_U[0] = qsp_0.copy()
                #Du_in = D_0
            else:
                temp = M.V[ilayer-1][0]
                QSp_U[0] = temp.QSp[temp.rank-1].copy()   #here must copy!
                QSp_U[0].reverse()
            
            QSp_U[1:4] = QSp_U[0].copy_many(3, reverse=[1, 2])
            
            QSp_V = [None]*rank_V
            QSp_V[:rank_V-1] = QSp_U[0].copy_many(rank_V-1)
            
            temp0 = QSp_V[0].copy()
            temp = QSp_V[0].copy()
            for i in range(rank_V-2):
                temp = temp.add(temp0)
                temp.update(qsp_max)
            #temp.update(qsp_max)
            temp.reverse() 
            QSp_V[rank_V-1] = temp

            totQN = self.qn_identity.copy()
            
            for j in range(1):
                M.U[ilayer][j] = M.U[ilayer][j].expand(QSp_U) 
                self.U_dag[ilayer][j] = self.U[ilayer][j].conjugate(2)
                
                M.V[ilayer][j] = M.V[ilayer][j].expand(QSp_V) 
                self.V_dag[ilayer][j] = self.V[ilayer][j].conjugate(self.V[ilayer][j].rank-1)
        
        if 1:
            ilayer = M.num_of_layer - 1
            
            #实际这里init_layer 是多余的，其不良后果是多增加了一个U layer，本不存在
            #M.init_layer(ilayer, rand_init=rand_init, unitary_init = M.unitary_init)
            
            QSp = [None for i in range(3)]
            
            j = 0
            #topQN = self.topQN
            ten = M.V[ilayer-1][j]
            #print "lll", len(QSp), len(M.V), ilayer-1, j, ten.rank        
            QSp[0] =ten.QSp[ten.rank-1].copy()   #下面一层V tensor超上的那条腿
            QSp[0].reverse()
            QSp[1] = QSp[0].copy()

            #QSp[2] = QSp_null.copy()   #top tensor 最上面那条腿实际上是1维的，等同于没有
            QSp[2] = M.qsp_null #.copy()   #top tensor 最上面那条腿实际上是1维的，等同于没有
            
            #这里实际上多余 topQN 本来就是QN_idendity
            QSp[2].QNs[0]  = M.topQN
            QSp[2].update()
            QSp[2].reverse()
            M.V[ilayer][j]=M.V[ilayer][j].expand(QSp)
            M.V_dag[ilayer][j]=M.V[ilayer][j].conjugate(M.V[ilayer][j].rank-1)

    def Shallow_Copy_Utensor(U1,U2):
        U2.ilayer = U1.ilayer
        U2.nSite = U1.nSite
        if  not allocated[U2.tensor]:  
            allocate[U2[0:0]]#U2.nSite-1]]
            
    def init_UnitMera(M):
        pass
    
    def Add2Set(m,n,A):
        for i in range(n):
            if m == A[i]:  
                Add2Set = i
                return
        Add2Set = n+1
        A[n+1] = m
    
    def GetNeigh(i, neigh, nSite):
        GetNeigh=mod[i+neigh+nSite+nSite,nSite]
    
    def storage_UTensor(fn, UT):
        #write[fn]UT.ilayer, UT.nSite
        for i  in range(0):
            Storage_Tensor(fn, UT[i])
    
    def restore_UTensor(fn, UT):
        #read[fn]UT.ilayer, UT.nSite
        for i  in range(0):
            Restore_Tensor(fn, UT[i])

    @staticmethod
    def unit_disentanger(U):
        """
        status_1_uncheck
        """
        if U.rank  !=  4:  
            print 'Not a disentanger'
            exit()
        return iTensor.unit_tensor(U)

    @staticmethod
    def unit_isometry(W):
        """
        this method is never used
        """
        iQN = np.ndarray(W.rank, "int")
        Dims= np.ndarray(W.rank, np.int)
        W.data[:] = 0.0
        rank = W.rank
        for idx  in range(W.nidx):
            p = W.Block_idx[0,idx]
            size = W.Block_idx[1,idx]
            iQN[0:rank] = W.Addr_idx[0:rank,idx]

            # Fix me, only work for parity symmetry            
            # Final leg is odd parity            
            if iQN[rank-1] ==  1:  
                for i  in range(size):
                    #W.data[p-1+i] = np.random.rand()-0.5
                    W.data[p-1+i] = crandom.rand()-0.5
                    pass
                    continue 
            # at least one leg is not parity
            if iQN[0]  ==  1 or iQN[1]  ==   1:  
                #forall[i=0:rank] Dims[i]=W.QSp[i].Dims[iQN[i]]
                for i in range(rank):
                    Dims[i]=W.QSp[i].Dims[iQN[i]]
                for i  in range(size):
                    #here is to judge wether has odd parity legs
                    pos = common_util.matrix_get_position_rev(i,Dims)
                    if 1 in pos:
                        X = 0.0
                    else:
                        pass
                        #X = np.random.rand()-0.5
                        X = crandom.rand()-0.5
                    W.data[p+i] = X                    
                continue
            
            # all of QNs are even parity
            if iQN[0] == 0 and iQN[1] == 0:  
                #forall[i=0:rank] Dims[i]=W.QSp[i].Dims[iQN[i]]
                for i in range(rank):
                    Dims[i]=W.QSp[i].Dims[iQN[i]]
                for i  in range(size):
                    pos= common_util.matrix_get_position_rev(i, Dims)
                    if pos[rank-1] == 0:  
                        X = 0.0
                    else:
                        #X = np.random.rand()-0.5
                        X = crandom.rand()-0.5
                    W.data[p+i] = X                    
        W.data[0] = 0.0
        #q299
        #W = Tensor_svd(W)
        #W.svd(2)
        Tensor_svd.svd(W, 2)
        return W
    
    def __repr__old(self, keys=[]):
        from pprint import pprint
        
        if keys== [] :
            #keys=["num_of_layer","qsp_0","qsp_max","qsp_max2"]
            keys=["num_of_layer"]
        dic=self.__dict__
        res=""
        for k in keys:
            temp = repr(self.__dict__[k])
            #if k in ["num_of_layer"]:
            res+= k+":\n\t" + temp +"\n" 
        return res

    def __repr__(self,layers=[], which=None, fewer=True, round=5):
        """  which  = ["U", "V", "U_dag", "V_dag"]"""
        
        if layers==[] :
            layers=range(self.num_of_layer)
        res= ""
        name = ["U", "V"]
        if which != None :
            name = which
        for n in name:
            #print n
            res += "%s\n"%(n)
            for i in layers:
                tensor = self.__getattribute__(n)[i][0]
                if fewer:
                    #res += "\t%d %s   %s   %d  buff \n"%(i, tensor.data[:4].round(round), tensor.data[-4:].round(round), tensor.totDim)
                    res += "\t%d %s   %s   %d  buff %s\n"%(i, tensor.data[:4].round(round), tensor.data[-4:].round(round), tensor.totDim, tensor.buf_ref)
                    #print "\t", i, tensor.data[:4].round(digit),"  ",tensor.data[-4:].round(digit), "  ", tensor.totDim
                else:
                    res +=  "\t%d%s\n"%(i, tensor.data.round(round))
                    #print "\t", i, tensor.data.round(digit)
        return res

    def to_mps(self): 
        """
            this is an idea:  of mine:  can one change mps to mera and vise versa?
            mix the two algorithms
        """
    
    def reset_to_identity(self): 
        """
            reset each tensor in mera to identify map 
            this is used in e.g test J for long range 
        """
        M = self 
        for i in range(M.num_of_layer-1):
            t = M.V[i][0]
            rank = t.rank
            qDims= np.empty(rank, "int")
            iDims= np.empty(rank, "int")
            qDims[:]= 0#[0, 0, 0, 0]
            iDims[:] = 0
            t.data[:] = 0.0
            t.set_element(qDims, iDims, 1.0)
            
            t = M.V_dag[i][0]
            rank = t.rank
            qDims= np.empty(rank, "int")
            iDims= np.empty(rank, "int")
            qDims[:]= 0 #[0, 0, 0, 0]
            iDims[:] = 0
            t.data[:] = 0.0
            t.set_element(qDims, iDims, 1.0)
        
    
def reset_mera(M):
    """
        
    """
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

class test_Mera(object):
    @classmethod
    def instance(cls, trunc_dim=2, tot_layer=4, symmetry="Z2", combine_2site=False):
        """
        borrowed from main.py for debuging
        """
        from quantum_number import init_System_QSp

        #from quantum_number import QSp_base, QN_idendity
        from tensor_network import TensorNetwork

        QN_idendity, QSp_base , QSp_null = init_System_QSp(symmetry=symmetry)


        #tensor.init_Buffer()

        qsp_0= QSp_base.copy()
        qsp_0.update()

        
        if symmetry == "Z2" :
            qsp_max = qsp_0.copy()
            qsp_max.Dims[0:2] = trunc_dim/2
            qsp_max.RefQN[0,0:2] = [1,1]
            qsp_max.RefQN[1,0:2] = [0,0]
            qsp_max.update()


            qsp_max2 = qsp_0.copy()
            qsp_max2.Dims[0:2] = 2
            qsp_max2.RefQN[0,0:2] = [2,1]
            qsp_max2.RefQN[1,0:2] = [0,1]
            qsp_max2.update()
        elif symmetry == "U1":  
            qsp_max =QspU1.max(trunc_dim)
            qsp_max2 = qsp_max.copy()
            
        elif symmetry == "Travial" :
            qsp_max= QSp_base.max(trunc_dim)
            qsp_max2 = qsp_max.copy()
        elif symmetry == "Z3":
            qsp_max= QSp_base.max(trunc_dim)
            qsp_max2 = qsp_max.copy()

        else:
            raise
        
        topQN= QN_idendity.copy()
        
        crandom.mrandom.csrand(1234)
        nTop = 2
        M= Mera(tot_layer, nTop, qsp_0, qsp_max, qsp_max2, topQN, QSp_null, QN_idendity)
        #sys.init_Hamiltonian(M)
        return M

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    def test_example(self): 
        m = Mera.example() 
        print_vars(vars(),  ['m.num_of_input'])
    
    def test_coarse_grain_map(self): 
        
        X = [55]; tau = 0; n = 1
        res=coarse_grain_map(X, tau, n)
        self.assertTrue(res==(18, 19)) 
        
        if 1: 
            X = [1, 2, 3]; tau = 0; n = 1
            for i in range(6): 
                X = range(i, i + 3)
                res = coarse_grain_map(X, tau, n) 
                print X, res 
        
        res = SupportSet.split_disconnect(((1, 2, 3,  5, 6, 7, 10, 11)))
        self.assertEqual(res, ((1, 2, 3), (5, 6, 7), (10, 11)))

    def test_temp(self): 
        
        X = [1, 2, 3]; tau = 0; n = 1
        for i in range(6): 
            #X = range(i, i + 3)
            X = [0, 1, 4, 5,  18, 19, 22]
            res = coarse_grain_map(X, tau, n) 
            print X, res 
        
        
    
if __name__=='__main__':
    """   
        ---- pass
        init val of mera, 
        trunc dim =4
            4
            {'h': 1.0, 'tot_Layer': 4, 'J_NN': -1.0, 'gamma': 1.0}
            U
                0 [ 1.  0.  1.  0.  0.]    [ 0.  0.  1.  0.  1.]    8
                1 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
                2 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
                3 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
            V
                0 [-0.74031 -0.31417 -0.0997   0.74803 -0.11271]    [ 0.86359 -0.13002  0.3045   0.80994  0.34728]    16
                1 [-0.18833  0.00821 -0.24003 -0.02472 -0.20739]    [-0.03192 -0.18542  0.29677 -0.14335  0.19232]    128
                2 [-0.05783 -0.0474  -0.26687 -0.23953 -0.22357]    [ 0.02352 -0.29011  0.1772   0.26828 -0.14946]    128
                3 [-0.51121  0.0521  -0.25331 -0.22125  0.23093]    [-0.22125  0.23093  0.41798  0.39995  0.48461]    8



             U at layer           1           8
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01 ,0.0000E+00 ,0.1000E+01
             U at layer           2         128
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01
             U at layer           3         128
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01

             V at layer            1          16
            T%data= ,-.7403E+00 ,-.3142E+00 ,-.9970E-01 ,0.7480E+00 ,-.1127E+00 ,0.8636E+00 ,-.1300E+00 ,0.3045E+00 ,0.8099E+00 ,0.3473E+00
             V at layer            2         128
            T%data= ,-.1883E+00 ,0.8213E-02 ,-.2400E+00 ,-.2472E-01 ,-.2074E+00 ,-.3192E-01 ,-.1854E+00 ,0.2968E+00 ,-.1434E+00 ,0.1923E+00
             V at layer            3         128
            T%data= ,-.5783E-01 ,-.4740E-01 ,-.2669E+00 ,-.2395E+00 ,-.2236E+00 ,0.2352E-01 ,-.2901E+00 ,0.1772E+00 ,0.2683E+00 ,-.1495E+00
             V at layer            4           8
            T%data= ,-.5112E+00 ,0.5210E-01 ,-.2533E+00 ,-.2212E+00 ,0.2309E+00 ,-.2212E+00 ,0.2309E+00 ,0.4180E+00 ,0.4000E+00 ,0.4846E+00

            {'h': 1.0, 'tot_Layer': 4, 'J_NN': -1.0, 'gamma': 1.0}
            4
            U_dag
                0 [ 1.  0.  1.  0.  0.]    [ 0.  0.  1.  0.  1.]    8
                1 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
                2 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
                3 [ 1.  0.  0.  0.  0.]    [ 0.  0.  0.  0.  1.]    128
            V_dag
                0 [-0.74031 -0.31417 -0.42745 -0.20222 -0.37996]    [-0.48608 -0.65521  0.32476  0.80994  0.34728]    16
                1 [-0.18833 -0.21406  0.00821 -0.1828  -0.24003]    [ 0.29677 -0.2769  -0.14335 -0.02159  0.19232]    128
                2 [-0.05783  0.01005 -0.0474  -0.00325 -0.26687]    [ 0.1772   0.09246  0.26828  0.24036 -0.14946]    128
                3 [-0.51121  0.0521  -0.25331 -0.22125  0.23093]    [-0.22125  0.23093  0.41798  0.39995  0.48461]    8


             U_dag at layer           1           8
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01 ,0.0000E+00 ,0.1000E+01
             U_dag at layer           2         128
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01
             U_dag at layer           3         128
            T%data= ,0.1000E+01 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.1000E+01

             V_dag at layer            1           8
            T%data= ,-.7403E+00 ,-.3142E+00 ,-.4275E+00 ,-.2022E+00 ,-.3800E+00 ,-.4861E+00 ,-.6552E+00 ,0.3248E+00 ,0.8099E+00 ,0.3473E+00
             V_dag at layer            2         128
            T%data= ,-.1883E+00 ,-.2141E+00 ,0.8213E-02 ,-.1828E+00 ,-.2400E+00 ,0.2968E+00 ,-.2769E+00 ,-.1434E+00 ,-.2159E-01 ,0.1923E+00
             V_dag at layer            3         128
            T%data= ,-.5783E-01 ,0.1005E-01 ,-.4740E-01 ,-.3247E-02 ,-.2669E+00 ,0.1772E+00 ,0.9246E-01 ,0.2683E+00 ,0.2404E+00 ,-.1495E+00
             V_dag at layer            4         128
            T%data= ,-.5112E+00 ,0.5210E-01 ,-.2533E+00 ,-.2212E+00 ,0.2309E+00 ,-.2212E+00 ,0.2309E+00 ,0.4180E+00 ,0.4000E+00 ,0.4846E+00
        U1 symm.
        tot_layer = 4
             U at layer           1
             -.85E+00 -.92E-01 -.12E+00 0.26E+00 0.43E+00      -.10E+00 0.26E+00 0.89E+00 0.54E-01 -.10E+01          70
             U at layer           2
             0.11E+00 0.60E+00 -.19E+00 0.33E+00 -.53E+00      -.54E+00 -.31E+00 0.52E+00 -.25E+00 0.10E+01          70
             U at layer           3
             -.22E+00 0.58E+00 0.43E+00 -.27E+00 -.45E+00      0.66E+00 -.21E+00 -.20E+00 0.24E+00 -.10E+01          70

             V at layer            1
             0.25E+00 -.33E+00 0.30E+00 0.94E-01 0.14E+00      -.45E-01 0.17E+00 -.42E+00 0.84E-01 -.57E-01          70
             V at layer            2
             0.12E+00 0.14E+00 -.38E+00 0.20E-01 0.14E+00      -.63E-03 0.12E+00 -.17E+00 -.15E+00 -.26E+00          70
             V at layer            3
             -.28E+00 -.29E+00 0.32E+00 -.11E+00 -.12E+00      -.25E+00 -.30E+00 0.24E+00 -.33E+00 0.22E-01          70
             V at layer            4
             0.40E-01 -.39E+00 0.60E+00 0.10E+00 0.46E-01      -.39E+00 0.60E+00 0.10E+00 0.46E-01 -.69E+00           6
        
        tot_layer = 5
                 U at layer           1          70
                T%data= ,-.8480E+00 ,-.9178E-01 ,-.1177E+00 ,0.2608E+00 ,0.4251E+00 ,-.1033E+00 ,0.2565E+00 ,0.8916E+00 ,0.5407E-01 ,-.1000E+01
                 U at layer           2        1120
                T%data= ,0.1081E+00 ,0.5980E+00 ,0.0000E+00 ,0.0000E+00 ,-.1922E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                 U at layer           3        1120
                T%data= ,-.2175E+00 ,0.5756E+00 ,0.0000E+00 ,0.0000E+00 ,0.4282E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                 U at layer           4        1120
                T%data= ,-.1314E+00 ,0.1278E+00 ,0.0000E+00 ,0.0000E+00 ,-.6447E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00

                 V at layer            1         140
                T%data= ,0.2511E+00 ,-.3287E+00 ,0.2961E+00 ,0.9387E-01 ,0.1428E+00 ,0.0000E+00 ,0.8431E-01 ,0.0000E+00 ,-.5694E-01 ,0.0000E+00
                 V at layer            2        1120
                T%data= ,0.1179E+00 ,0.1395E+00 ,0.0000E+00 ,0.0000E+00 ,-.3763E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                 V at layer            3        1120
                T%data= ,-.2838E+00 ,-.2928E+00 ,0.0000E+00 ,0.0000E+00 ,0.3215E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                 V at layer            4        1120
                T%data= ,0.3145E+00 ,0.4537E-01 ,0.0000E+00 ,0.0000E+00 ,-.2488E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                 V at layer            5          24
                T%data= ,0.2437E+00 ,-.3740E+00 ,0.0000E+00 ,0.0000E+00 ,0.7859E+00 ,0.0000E+00 ,0.1305E+00 ,0.0000E+00 ,0.0000E+00 ,0.0000E+00
                5
                U
                    0 [-0.84801 -0.09178 -0.11775  0.26075  0.42509]    [-0.10332  0.25653  0.89162  0.05407 -1.     ]    70
                    1 [ 0.10813  0.59799 -0.19221  0.33068 -0.52952]    [-0.53539 -0.30699  0.5178  -0.24912  1.     ]    70
                    2 [-0.2175   0.57558  0.42816 -0.26694 -0.44595]    [ 0.66    -0.21431 -0.19962  0.23821 -1.     ]    70
                    3 [-0.13142  0.12784 -0.64466  0.50232  0.67435]    [ 0.30165 -0.18316  0.34321  0.59646 -1.     ]    70
                    4 [ 0.25343 -0.60954  0.48866  0.03674 -0.41219]    [-0.53241 -0.34343 -0.42782 -0.07452  1.     ]    70
                V
                    0 [ 0.2511  -0.32872  0.29614  0.09387  0.14276]    [-0.04481  0.17229 -0.42421  0.08431 -0.05694]    70
                    1 [ 0.11786  0.13954 -0.37626  0.0199   0.13874]    [-0.00063  0.12485 -0.17229 -0.1512  -0.26182]    70
                    2 [-0.28381 -0.2928   0.32147 -0.11308 -0.11924]    [-0.25076 -0.29513  0.24379 -0.32932  0.02217]    70
                    3 [ 0.31454  0.04537 -0.24882 -0.25035  0.25795]    [ 0.31649  0.253    0.06034  0.19135  0.20701]    70
                    4 [ 0.24367 -0.37398  0.78586 -0.27088  0.30458]    [-0.37398  0.78586 -0.27088  0.30458  0.13054]    6
    """ 

    pass
    if 0: 
        tm = test_Mera
        #M = tm.instance(trunc_dim=4, tot_layer=4, symmetry="U1")
        #M = tm.instance(trunc_dim=4, tot_layer=4, symmetry="Travial")
        M = tm.instance(trunc_dim=6, tot_layer=4, symmetry="Z3")
        #M = tm.instance(trunc_dim=11, tot_layer=4, symmetry="Travial")
        print M
        print type(M.topQN.val)
        #M.show_param()

        
        #print M.__repr__(which=["U_dag", "V_dag"],  fewer=1)  
        #print M.__repr__(which=["U", "V"],  fewer=1)  
        

        #print M.V[0][0]
        #print M.V[M.num_of_layer-1][0]
    
    if 0: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
        
           #'test_example',  
           #'test_coarse_grain_map',  
           'test_temp',  
           
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))
        unittest.TextTestRunner().run(suite)

   







