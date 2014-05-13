#!/usr/bin/env python
#coding=utf8
import pprint
import sys
import math
import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs , eigsh

from merapy.ascending import ascending_ham
from merapy.tensor_svd import Tensor_svd
from merapy.hamiltonian import System
from merapy.tensor_py import iTensor
from merapy.quantum_number import *
from merapy.decorators import timer, set_STATE_end_1, tensor_player, decorate_methods

__all__ = ["calc_scaling_dim_1site", "calc_scaling_dim_2site"]


func = lambda x, b: round(-math.log(abs(x), b), 8)

from merapy.graphics import Graphics
def init_g2_modified():
    """
        insert Vg into original G2
    """
    G2 = {}
    order2 = {}
    for i in [-1, 0, 1]:
        G2[i] = Graphics()    
        G2[i].graph_name = "G2[" + str(i) + "]"
        order2[i] = np.ndarray(Graphics.MaxNode, "int")        
    #order2 = np.ndarray((Graphics.MaxNode, 3), "int")        
    
    weight2 = {}
    for k in G2.keys():
        weight2[k] = 1/3.0    
    
    if 1:#G2[R] insert two Vg between V_0 and V_p0
        G2[-1].size=10
        #G2[-1].nodes=-1
        #G2[-1].edges=-1
        G2[-1].names[0]="OO" #OO
        G2[-1].nodes[0]=4
        G2[-1].edges[0:4,0]=[0,1,2,3]

        G2[-1].names[1]="Up" #Up_1
        G2[-1].nodes[1]=4
        G2[-1].edges[0:4,1]=[4,5,6,7]

        G2[-1].names[2]="oo" #oo_2_3
        G2[-1].nodes[2]=4
        G2[-1].edges[0:4,2]=[6,7,8,9]

        G2[-1].names[3]="U" #U_1
        G2[-1].nodes[3]=4
        G2[-1].edges[0:4,3]=[8,9,10,11]

        G2[-1].names[4]="Vp" #Vp_0
        G2[-1].nodes[4]=4
        G2[-1].edges[0:4,4]=[0,100,200,4]

        G2[-1].names[5]="Vp" #Vp_1
        G2[-1].nodes[5]=4
        G2[-1].edges[0:4,5]=[1,5,14,15]

        G2[-1].names[6]="V" #V_1
        G2[-1].nodes[6]=4
        G2[-1].edges[0:4,6]=[11,14,15,3]

        G2[-1].names[7]="V" #V_0
        G2[-1].nodes[7]=4
        G2[-1].edges[0:4,7]=[12,13,10,2]
        #$ 6/7            2.030001E+14
        
        #insert Vg
        G2[-1].names[8] = "Vg"
        G2[-1].nodes[8]=2
        G2[-1].edges[0:2,8]=[100,12]
        
        G2[-1].names[9] = "Vg"
        G2[-1].nodes[9]=2
        G2[-1].edges[0:2,9]=[200,13]

        order2[-1][0:10]=[3,2,1,5,6,0,4, 8, 9,  7]
    
    if 1: #G2[L] insert 3  Vg
        G2[0].size=11
        #G2[0].nodes=-1
        #G2[0].edges=-1
        G2[0].names[0]="OO" #OO
        G2[0].nodes[0]=4
        G2[0].edges[0:4,0]=[0,1,2,3]

        G2[0].names[1]="Up" #Up_1
        G2[0].nodes[1]=4
        G2[0].edges[0:4,1]=[4,5,300,7]

        G2[0].names[2]="U" #U_1
        G2[0].nodes[2]=4
        G2[0].edges[0:4,2]=[6,8,9,10]

        G2[0].names[4]="Vp" #Vp_1
        G2[0].nodes[4]=4
        G2[0].edges[0:4,4]=[1,5,13,14]

        G2[0].names[5]="oo" #oo_3_4
        G2[0].nodes[5]=4
        G2[0].edges[0:4,5]=[7,13,8,15]

        G2[0].names[6]="V" #V_1
        G2[0].nodes[6]=4
        G2[0].edges[0:4,6]=[10,15,14,3]

        G2[0].names[7]="V" #V_0
        G2[0].nodes[7]=4
        G2[0].edges[0:4,7]=[11,12,9,2]

        G2[0].names[3]="Vp" #Vp_0
        G2[0].nodes[3]=4
        G2[0].edges[0:4,3]=[0,100, 200,4]

        G2[0].names[8]="Vg" 
        G2[0].nodes[8]=2
        G2[0].edges[0:2,8]=[100, 11]

        G2[0].names[9]="Vg" 
        G2[0].nodes[9]=2
        G2[0].edges[0:2,9]=[200, 12]

        G2[0].names[10]="Vg"
        G2[0].nodes[10]=2
        G2[0].edges[0:2,10]=[300, 6]

        #$ 6/8         2.01020001E+16
        order2[0][0:11]=[2, 10, 1, 5, 6, 4, 0, 3, 8, 9, 7]
    
    if 1: #G2[C] insert 1 Vg
        G2[1].size=9
        #G2[1].nodes=-1
        #G2[1].edges=-1
        G2[1].names[0]="Up" #Up_2
        G2[1].nodes[0]=4
        G2[1].edges[0:4,0]=[0,1,2,3]

        G2[1].names[1]="OO" #OO
        G2[1].nodes[1]=4
        G2[1].edges[0:4,1]=[4,5,6,7]
        G2[1].names[2]="U" #U_2
        G2[1].nodes[2]=4
        G2[1].edges[0:4,2]=[8,3,9,10]

        G2[1].names[3]="Vp" #Vp_1
        G2[1].nodes[3]=4
        G2[1].edges[0:4,3]=[4,100,12,0]

        G2[1].names[4]="Vp" #Vp_2
        G2[1].nodes[4]=4
        G2[1].edges[0:4,4]=[5,1,13,14]
        G2[1].names[5]="oo" #oo_4_5
        G2[1].nodes[5]=4
        G2[1].edges[0:4,5]=[12,2,15,8]

        G2[1].names[6]="V" #V_1
        G2[1].nodes[6]=4
        G2[1].edges[0:4,6]=[11,15,9,6]

        G2[1].names[7]="V" #V_2
        G2[1].nodes[7]=4
        G2[1].edges[0:4,7]=[10,13,14,7]
        #$ 6/8         2.01020001E+16
        
        #insert Vg
        G2[1].names[8]="Vg" #V_1
        G2[1].nodes[8]=2
        G2[1].edges[0:2,8]=[100, 11]

        order2[1][0:9]=[2,0,5,6, 8, 3,1,4,7]

    return G2, weight2, order2

G2, weight2, order2 = init_g2_modified()

def ascending_ham_modified(M, S, ilayer, tau=None, info=0):
    """
        ilayer: ranges from 0 to M.num_of_layer-2
        calculate H_2 for ilayer+1, see Vidal 2007 fig.14
        this is the average ascending operator  \bar A =1/3(A_R + R_C + R_L) for ternary 
        act on S.H_2[ilayer].O[j], obtain S.H_2[ilayer+1].O[j].data
        注意这里与f90中不同，f90中是把 ilayer-1 ascending 到 ilayer
        
        tau: tau is used only when calculate coup coeff for long range models
    """
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
    for g in [-1, 0, 1]:
        G = G2[g]
        order = order2[g]
        weight=weight2[g]
        temp = S.add_env(M, ilayer, G, order, name, weight, info=info-1)
        S.H_2[ilayer_p1][j].data += temp.data  #*weight
        
    
    if info>0:
        print "H_2 at iLayer+1,",ilayer_p1
        print S.H_2[ilayer_p1].O[j].data[:4].round(3),"....", S.H_2[ilayer_p1].O[j].data[-4:].round(3) 
        print "END ascending_ham"

def asd_op_1site(V):
    """
    open legs are: 200, 500, 300, 400
    """
    r = V.rank
    V_dag = V.conjugate(r-1 )
    #lab1 = [0, 1, 2, 3]
    #lab2 = [4, 0, 5, 2]
    lab1 = range(r) 
    lab1[1] = 100; lab1[-1] = 300
    lab2 = [400] + range(r-1)
    lab2[2] = 500
    op, lab12 = V.contract(V_dag, lab1, lab2)
    return op

def Vg(qsp):
    """
        matrix repr of element g of group Z2, U1, etc
    """
    symmetry = qsp.QnClass.SYMMETRY
    qn = qsp.QnClass.qn_id()
    qsp = qsp.copy_many(2, reverse=[1])
    oo = iTensor(2, qsp, qn)
    oo.data[:] = 0.0
    if symmetry == "Z2":
        d = qsp[0].totDim//2
        oo.data[:oo.totDim//2] = np.identity(d).ravel()
        oo.data[oo.totDim//2:] = -np.identity(d).ravel()
    elif symmetry == "Travial": 
        #this wont give non-local operator, write this only avoid breaking program
        #without symmetry it is not easy to implement that. The point is that one dont know how to identify different subspaces for each quantum number
        #qsp = V.QSp[0].copy_many(2, reverse=[1])
        dim = qsp[0].totDim
        oo.data[:] = np.identity(dim).ravel()
    elif symmetry == "U1":
        for i in range(oo.nidx):
            p = oo.Block_idx[0, i]
            size = oo.Block_idx[1, i]
            dim = int(math.sqrt(size))
            if dim**2 != size: 
                raise Exception("%s"%oo)
            qsp = oo.QSp[0]
            a = oo.Addr_idx[:, i][0]
            m = qsp.QNs[a].val
            #print i, q
            sign = 1 if m%2== 0 else -1     #under rotation \theta->theta + pi, state |e^{-im\theta}> produce a minus sign if m is odd
            oo.data[p:p + size] = sign*np.identity(dim).ravel()
    return oo

def asd_op_1site_nonlocal(V):
    V_dag = V.conjugate(V.rank-1 )
    lab1 = [0, 1, 2, 3]
    lab2 = [6, 0]
    lab3 = [4, 6, 5, 2]
    symmetry = V.QSp[0].QnClass.SYMMETRY
    oo = Vg(V.QSp[0].copy())
    V_oo , lab_V_oo = V.contract(oo, lab1, lab2)
    V_oo_Vd, lab_V_oo_Vd= V_oo.contract(V_dag, lab_V_oo, lab3)
    
    return V_oo_Vd

def asd_op_1site_another(V, V_dag, H1, Vin, iter=None):
    H1.data = Vin
    lab1 = [0, 1, 2, 3]
    lab2 = [4, 0, 5, 2]
    lab3 = [5, 1]
    if iter.iter is not None: 
        set_STATE_end_1(iter=iter.iter, record_at=0, verbose=False)    
        iter.iter += 1 
    V_Vd, lab_v_vd= V.contract(V_dag, lab1, lab2)
    op, lab_op = V_Vd.contract(H1, lab_v_vd, lab3)
    
    Vout = op.data
    return Vout

def calc_scaling_dim_1site(path=None, S=None, iterate=False, local=True, qn_list=None, tol=None):
    state_bac = tensor_player.STATE
    tensor_player.STATE = "stop"
    tol =  tol if tol is not None else 0
    if S is None:
        S = System.load(path)
    #msg = "energy = %f, symmetry = %s"%(S.energy, S.symmetry)
    energy_err = (S.energy - S.energy_exact)
    print "energy = %1.15f, eng_err = %1.1e, symmetry = %s, trunc_dim = %d"%(
            S.energy, energy_err, S.symmetry, S.mera.qsp_max.totDim)
    if 1:
        symmetry = S.symmetry
        if qn_list is None:
            qn_list = symmetry_to_Qn(symmetry).QNS
            if symmetry == "U1":
                qn_list = range(3)
        qn_instance_list = [qn_factory(symmetry, i) for i in qn_list]
    
    #period = S.mera.period

    SIlayer = S.mera.num_of_layer-2
    V = S.mera.V[SIlayer][0].copy()
    rescale_factor = V.rank-1
    func = lambda x: round(-math.log(x, rescale_factor), 8)
    if not  iterate:
        if local:
            op = asd_op_1site(V)
        else:
            op = asd_op_1site_nonlocal(V)
        op = op.permutation([1, 2, 0, 3])

        #if not hasattr(op.QSp[0], "empty"):   #only temporarily for compatbile
        #    op.QSp[0].__class__.empty = empty
        temp = Tensor_svd.eig_sparse(op, qn_instance_list, tol=tol)
        
        #print "qns\tconformal dims"
        if 0:
            for i in temp:
                qn = temp[i][0]
                vals = abs(temp[i][1]) 
                vals.sort()
                vals= vals[-1::-1]
                try:
                    print qn, map(func, vals[:])
                except:
                    #raise Exception(str(vals))
                    msg = "vals are %s"%str(vals)
                    print msg
        else:
            for i, vals in temp.items():
                vals = np.abs(vals)
                vals.sort()
                vals = vals[-1::-1]
                print i, map(func, vals[:])

    else:
        print "iteratively eigs"

        V_dag = S.mera.V_dag[SIlayer][0]
        if 1:
            ilayer = SIlayer
            S.rho_2[ilayer][0].__class__= iTensor
            S.mera.V[ilayer][0].__class__= iTensor
            S.mera.V_dag[ilayer][0].__class__= iTensor
            S.mera.U[ilayer][0].__class__= iTensor
            S.mera.U_dag[ilayer][0].__class__= iTensor
        
        msg = ""
        qsp0 = V.QSp[0].copy()
        for qn in qn_instance_list:
            #print "%s"%qn, 
            msg += "\t%s "%qn 
            H1 = iTensor(2, qsp0.copy_many(2, reverse=[1]), qn)
            dim = H1.data.size
            calc_scaling_dim_1site.iter = 0
            temp =lambda vec:asd_op_1site_another(V, V_dag, H1, vec, iter=calc_scaling_dim_1site)
            A = LinearOperator((dim, dim), matvec=temp) 
            k = 15
            v0 = np.random.random(dim)-0.5
            if dim-1> k:
                vals, vecs = eigs(A=A, k=k, v0=v0, tol=tol)
            elif dim<3:
                #vals, vecs = eigs(A=A, k=1, tol=tol)
                print "dim too low, continue with next qn. dim = %d"%dim
                continue
            else:
                vals, vecs = eigs(A=A, k=dim-2, v0=v0, tol=tol)
            
            set_STATE_end_1(iter=None, record_at=None, stop_at=None, verbose=False, power_on=False)
            #print "num of iter = %d"%calc_scaling_dim_1site.iter, 
            vals= np.abs(vals)
            vals.sort()
            vals= vals[-1::-1]
            msg += "%s\n"%str(map(func, vals[:15])) 
        print msg

    tensor_player.STATE = state_bac

def asd_op_2site(ascending_func, S, Vin, local=True, info=0, iter=None):
    M = S.mera
    ilayer = M.num_of_layer-2
    ilayer_p1 = ilayer + 1
    if not local:
        if not S.network_tensor_dic.has_key("Vg"):
            S.network_tensor_dic.update({"Vg":("Vg", 0)})
        if not hasattr(S, "Vg"):
            temp = [[None]]*ilayer_p1
            temp[ilayer][0] = Vg(S.mera.V[ilayer][0].QSp[0].copy())
            S.__dict__["Vg"] = temp

    set_STATE_end_1(iter=iter.iter, record_at=0, verbose=False)    
    iter.iter += 1 
    S.H_2[ilayer][0].data = Vin
    if local:
        ascending_func(M, S, ilayer, info=info)
    else:
        ascending_ham_modified(M, S, ilayer, info=info)
    Vout = S.H_2[ilayer_p1][0].data
    #print "vvv", Vout.round(5)[:20]
    return Vout 

#@timer
def calc_scaling_dim_2site_bac(path=None, ascending_func=ascending_ham, local=True, S=None, qn_list=None, tol=None, k=None, verbose=True):
    """
        I compareed the results of NN heisbg model with J_nn = 0.241186 with that of Heisenberg/CritiDimention.f90 they coincied for different qns for all digits
    """
    
    tol = tol if tol is not None else 0
    is_loaded = False
    if S is None: 
        S = System.load(path)
    if 1:
        record = S.record
        try:
            #print 'try to fetch scaling_dim from S.record...  ', 
            msg = 'try to fetch scaling_dim from S.record...  '
            res = record[S.iter]['scaling_dim']
        except KeyError as err:
            #print 'record not found'
            #print err
            msg += 'record not found:  '  + err.message
            if verbose:  print msg
            is_loaded = True
        else:
            #print 'record found'
            msg  += 'record found'
            if verbose:  print msg
            return res
        
    #print "energy = %f, symmetry = %s"%(S.energy, S.symmetry)
    
    M = S.mera
    ilayer = M.num_of_layer-2
    ilayer_p1 = ilayer + 1
    data_bac = {}
    data_bac[ilayer] = S.H_2[ilayer][0];   data_bac[ilayer_p1] = S.H_2[ilayer_p1][0]
    temp1 = S.H_2[ilayer][0].data[:10].copy()
    use_player = True
    if use_player:
        state_bac = tensor_player.STATE
        meth_names = ["set_data_entrance", "contract_core", "permutation"]
        meth_bac = {}
        meth_deced = {}
        for i in meth_names:
            meth_bac[i] = iTensor.__dict__[i]
            meth_deced[i] = tensor_player(which=i)(iTensor.__dict__[i])
            setattr(iTensor, i, meth_deced[i]) 
        #iTensor_temp = decorate_methods(decorator=tensor_player, meth_names=meth_names)(iTensor)
        
            #S.rho_2[ilayer][0].__class__= iTensor_temp
            #S.mera.V[ilayer][0].__class__= iTensor_temp
            #S.mera.V_dag[ilayer][0].__class__= iTensor_temp
            #S.mera.U[ilayer][0].__class__= iTensor_temp
            #S.mera.U_dag[ilayer][0].__class__= iTensor_temp
    
    symmetry = S.symmetry
    if qn_list is None:
        qn_list = symmetry_to_Qn(symmetry).QNS
        if symmetry == "U1":
            qn_list = range(3)
    qn_instance_list = [qn_factory(symmetry, i) for i in qn_list]

    qsp0 = S.mera.V[ilayer][0].QSp[0].copy()
    rescale_factor = S.mera.V[ilayer][0].rank-1
    func = lambda x: round(-math.log(x, rescale_factor), 8)
    
    print "2site, b = %d"%rescale_factor
    res= {}
    for qn in qn_instance_list:
        #print "%s"%qn, 
        S.H_2[ilayer][0] = iTensor(4, qsp0.copy_many(4, reverse=[2, 3]), qn.copy())
        S.H_2[ilayer_p1][0] = iTensor(4, qsp0.copy_many(4, reverse=[2, 3]), qn.copy())

        only_NN = S.only_NN
        S.only_NN = True
        calc_scaling_dim_2site.iter = 0

        temp = lambda V: asd_op_2site(ascending_func, S, V, local=local, info=0, iter=calc_scaling_dim_2site)
        dim = S.H_2[ilayer][0].totDim
        A = LinearOperator((dim, dim), matvec=temp, dtype=np.float) 

        v0 = np.random.random(dim)-0.5

        k = k if k is not None else 15
        if dim-1> k:
            vals, vecs = eigs(A=A, k=k, v0=v0, tol=tol)
        elif dim<3:
            print "dim too low, continue with next qn. dim = %d"%dim
            continue
        else:
            vals, vecs = eigs(A=A, k=dim-2, v0=v0, tol=tol)

        tensor_player.STATE = "stop"
        
        vals= np.abs(vals)
        #print "vv", vals
        vals.sort()
        vals= vals[-1::-1]
        try:
            #res.append((qn._val, map(func, vals[:k])))
            res[qn._val] = tuple(map(func, vals[:k]))
        except ValueError as err:
            print err
            del S.G_2_2
            print 'S.G_2_2 is deleted, retry it'
            break
        
        S.only_NN = only_NN
        S.H_2[ilayer][0] = data_bac[ilayer]
        S.H_2[ilayer_p1][0] = data_bac[ilayer_p1]
        temp2 = S.H_2[ilayer][0].data[:10].copy()
        assert np.all(temp1==temp2), "H_2 not succesfully restored"
    if use_player:
        for i in meth_names:
            setattr(iTensor, i, meth_bac[i])
        tensor_player.STATE = state_bac
    
    if is_loaded:
        if not S.record.has_key(S.iter):  
            S.record[S.iter]={}
        S.record[S.iter].update({'scaling_dim':res})
        print 'scaling_dima at step %d is saved'%S.iter
        S._save(path)
        
    return res

def calc_scaling_dim_2site(S, ascending_func=ascending_ham, local=True,  qn_list=None, tol=None, k=None, verbose=True, **kwargs):
    """
        I compareed the results of NN heisbg model with J_nn = 0.241186 with that of Heisenberg/CritiDimention.f90 they coincied for different qns for all digits
    """
    
    tol = tol if tol is not None else 0
    is_loaded = False
    if 0:
        record = S.record
        try:
            #print 'try to fetch scaling_dim from S.record...  ', 
            msg = 'try to fetch scaling_dim from S.record...  '
            res = record[S.iter]['scaling_dim']
        except KeyError as err:
            #print 'record not found'
            #print err
            msg += 'record not found:  '  + err.message
            if verbose:  print msg
            is_loaded = True
        else:
            #print 'record found'
            msg  += 'record found'
            if verbose:  print msg
            return res
        
    #print "energy = %f, symmetry = %s"%(S.energy, S.symmetry)
    
    M = S.mera
    ilayer = M.num_of_layer-2
    ilayer_p1 = ilayer + 1
    data_bac = {}
    data_bac[ilayer] = S.H_2[ilayer][0];   data_bac[ilayer_p1] = S.H_2[ilayer_p1][0]
    temp1 = S.H_2[ilayer][0].data[:10].copy()
    use_player = True
    if use_player:
        state_bac = tensor_player.STATE
        meth_names = ["set_data_entrance", "contract_core", "permutation"]
        meth_bac = {}
        meth_deced = {}
        for i in meth_names:
            meth_bac[i] = iTensor.__dict__[i]
            meth_deced[i] = tensor_player(which=i)(iTensor.__dict__[i])
            setattr(iTensor, i, meth_deced[i]) 
        #iTensor_temp = decorate_methods(decorator=tensor_player, meth_names=meth_names)(iTensor)
        
            #S.rho_2[ilayer][0].__class__= iTensor_temp
            #S.mera.V[ilayer][0].__class__= iTensor_temp
            #S.mera.V_dag[ilayer][0].__class__= iTensor_temp
            #S.mera.U[ilayer][0].__class__= iTensor_temp
            #S.mera.U_dag[ilayer][0].__class__= iTensor_temp
    
    symmetry = S.symmetry
    if qn_list is None:
        qn_list = symmetry_to_Qn(symmetry).QNS
        if symmetry == "U1":
            qn_list = range(3)
    qn_instance_list = [qn_factory(symmetry, i) for i in qn_list]

    qsp0 = S.mera.V[ilayer][0].QSp[0].copy()
    rescale_factor = S.mera.V[ilayer][0].rank-1
    func = lambda x: round(-math.log(x, rescale_factor), 8)
    
    #print "2site, b = %d"%rescale_factor
    res= {}
    for qn in qn_instance_list:
        #print "%s"%qn, 
        S.H_2[ilayer][0] = iTensor(4, qsp0.copy_many(4, reverse=[2, 3]), qn.copy())
        S.H_2[ilayer_p1][0] = iTensor(4, qsp0.copy_many(4, reverse=[2, 3]), qn.copy())

        only_NN = S.only_NN
        S.only_NN = True
        calc_scaling_dim_2site.iter = 0

        temp = lambda V: asd_op_2site(ascending_func, S, V, local=local, info=0, iter=calc_scaling_dim_2site)
        dim = S.H_2[ilayer][0].totDim
        A = LinearOperator((dim, dim), matvec=temp, dtype=np.float) 

        v0 = np.random.random(dim)-0.5

        k = k if k is not None else 15
        if dim-1> k:
            vals, vecs = eigs(A=A, k=k, v0=v0, tol=tol)
        elif dim<3:
            print "dim too low, continue with next qn. dim = %d"%dim
            continue
        else:
            vals, vecs = eigs(A=A, k=dim-2, v0=v0, tol=tol)

        tensor_player.STATE = "stop"
        
        vals= np.abs(vals)
        #print "vv", vals
        vals.sort()
        vals= vals[-1::-1]
        try:
            #res.append((qn._val, map(func, vals[:k])))
            res[qn._val] = tuple(map(func, vals[:k]))
        except ValueError as err:
            print err
            del S.G_2_2
            print 'S.G_2_2 is deleted, retry it'
            break
        
        S.only_NN = only_NN
        S.H_2[ilayer][0] = data_bac[ilayer]
        S.H_2[ilayer_p1][0] = data_bac[ilayer_p1]
        temp2 = S.H_2[ilayer][0].data[:10].copy()
        assert np.all(temp1==temp2), "H_2 not succesfully restored"
    if use_player:
        for i in meth_names:
            setattr(iTensor, i, meth_bac[i])
        tensor_player.STATE = state_bac
    
        
    return res



if __name__ == "__main__":
    import os
    from merapy.common_util import PickleDb as PDB 

    if len(sys.argv)>1:
        fn = sys.argv[1]   
        is_plot = False
        try:
            if sys.argv[2] == "-p": is_plot = True
        except:
            pass
        if is_plot:
            pdb = PDB(fn=fn, path="./")
            pdb.plot_scaling_dim()
            pdb.plot_scaling_dim(calc_err=True)
            sys.exit(0)
        qn_list = None
        calc_scaling_dim_1site(fn, iterate=False, local=True, qn_list=qn_list)
        res = calc_scaling_dim_2site(fn, local=True, S=None, qn_list=None, tol=1e-8, k=10); print res
        sys.exit(0)
    #decorators.tensor_player.NEXT_STATE = "record"

    path = "./run-ising/" 
    path = "./run-xx-symm=trav/"
    path = "./run-xx/"
    path = "./run-heisbg-NNN/"
    path = "../current/run-long-better/"
    path = "../current/run-long-better-quat/"
    
    
    if 1:
        pass
        fn = path + "8.pickle"
    else:
        fn = None
    qn_list = [1] #[0, 1, 2]#[0, 1, 2, 3]

   
    if 0:
        pdb = PDB(fn=fn, path=path)
        pdb.plot_scaling_dim()
        pdb.plot_scaling_dim(calc_err=True)
        #pdb.plot_energy()
        #pdb.plot_energy_diff()
    else:
        fn = '/home/zhli/Documents/mera_backup_tensor/run-ising/4.pickle.ternary'
        os.environ["OMP_NUM_THREADS"] = str(2)
        #calc_scaling_dim_1site(fn, iterate=False, local=True, qn_list=qn_list)
        if 0: 
            pass
            #calc_scaling_dim_1site(fn, iterate=False, local=False, qn_list=qn_list)
            #calc_scaling_dim_1site(fn, iterate=True, local=True, qn_list=qn_list)
            #res = calc_scaling_dim_2site(fn, local=True, S=None, qn_list=qn_list, tol=1e-8, k=12)
        S= System.load(fn)
        res = calc_scaling_dim_2site(S=S, local=True, qn_list=qn_list, tol=1e-8, k=12); print res
        




