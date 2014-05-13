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
from merapy.graphics import Graphics

__all__ = ["calc_scaling_dim_1site", "calc_scaling_dim_3site"]


def ascending_ham_3site(M, S, ilayer, tau=None, info=0):
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
    G_2_3 = S.G_2_3
    G_3_3 = S.G_3_3
    
    if ilayer<0 or ilayer>M.num_of_layer-2:
        raise ValueError("ilayer should be in [0, %d]. ilayer=%s"%(M.num_of_layer-2, ilayer))

    if tau is None:
        tau = ilayer
    
    if info>0:
        print "START ascending_ham"
    
    i=0
    j = 0
    ilayer_p1 = ilayer + 1
    
    if 1:
        name = "OOO"
        S.H_3[ilayer_p1].O[j].data[:] = 0.0
        for G in [G_3_3]:
            S.add_env_many(M, ilayer, name, G, add_tensor=S.H_3[ilayer_p1][0], branch=0, info=info-1)


    if info>0:
        print "H_2 at iLayer+1,",ilayer_p1
        print S.H_2[ilayer_p1].O[j].data[:4].round(3),"....", S.H_2[ilayer_p1].O[j].data[-4:].round(3) 
        print "END ascending_ham"

def asd_op_3site(ascending_func, S, Vin, local=True, info=0, iter=None):
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
    #print "turned off player"
    set_STATE_end_1(iter=iter.iter, record_at=0, verbose=False)    
    iter.iter += 1 
    S.H_3[ilayer][0].data = Vin
    if local:
        ascending_ham_3site(M, S, ilayer, info=info)
    else:
        ascending_ham_modified(M, S, ilayer, info=info)
    Vout = S.H_3[ilayer_p1][0].data
    return Vout 


@timer
def calc_scaling_dim_3site(fn=None, ascending_func=ascending_ham_3site, local=True, 
        S=None, qn_list=None, tol=None, k=None):
    """
        I compareed the results of NN heisbg model with J_nn = 0.241186 with that of Heisenberg/CritiDimention.f90 they coincied for different qns for all digits
    """
    state_bac = tensor_player.STATE
    meth_names= ["set_data_entrance", "contract_core", "permutation"]
    meth_bac = {}
    for i in meth_names:
        meth_bac[i] = iTensor.__dict__[i]
    iTensor_temp = decorate_methods(decorator=tensor_player, meth_names=meth_names)(iTensor)
    
    tol = tol if tol is not None else 0
    if S is None:
        S = System.load(fn)
    #print "energy = %f, symmetry = %s"%(S.energy, S.symmetry)
    energy_err = (S.energy - S.energy_exact)
    print "energy = %1.15f, eng_err = %1.1e, symmetry = %s, trunc_dim = %d"%(
            S.energy, energy_err, S.symmetry, S.mera.qsp_max.totDim)

    print "3site"
    M = S.mera
    ilayer = M.num_of_layer-2
    ilayer_p1 = ilayer + 1
    data_bac = {}
    data_bac[ilayer] = S.H_3[ilayer][0];   data_bac[ilayer_p1] = S.H_3[ilayer_p1][0]
    temp1 = S.H_3[ilayer][0].data[:10].copy()
    
    S.rho_3[ilayer][0].__class__= iTensor_temp
    S.mera.V[ilayer][0].__class__= iTensor_temp
    S.mera.V_dag[ilayer][0].__class__= iTensor_temp
    S.mera.U[ilayer][0].__class__= iTensor_temp
    S.mera.U_dag[ilayer][0].__class__= iTensor_temp
    
    symmetry = S.symmetry
    if qn_list is None:
        qn_list = symmetry_to_Qn(symmetry).QNS
        if symmetry == "U1":
            qn_list = range(3)
    qn_instance_list = [qn_factory(symmetry, i) for i in qn_list]

    qsp0 = S.mera.V[ilayer][0].QSp[0].copy()
    rescale_factor = S.mera.V[ilayer][0].rank-1
    func = lambda x: round(-math.log(x, rescale_factor), 8)

    res= []
    for qn in qn_instance_list:
        #print "%s"%qn, 
        S.H_3[ilayer][0] = iTensor_temp(6, qsp0.copy_many(6, reverse=[3, 4, 5]), qn.copy())
        S.H_3[ilayer_p1][0] = iTensor_temp(6, qsp0.copy_many(6, reverse=[3, 4, 5]), qn.copy())

        only_NN = S.only_NN
        S.only_NN = True
        calc_scaling_dim_3site.iter = 0

        temp = lambda V: asd_op_3site(ascending_ham_3site, S, V, local=local, info=0, iter=calc_scaling_dim_3site)
        dim = S.H_3[ilayer][0].totDim
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

        #print "num of iter = %d"%calc_scaling_dim_3site.iter
        #set_STATE_end_1(iter=None, record_at=None, stop_at=None, verbose=False, power_on=False)
        tensor_player.STATE = "stop"
        
        vals= np.abs(vals)
        #print "vv", vals
        vals.sort()
        vals= vals[-1::-1]
        
        res.append((qn._val, map(func, vals[:k])))
        
        S.only_NN = only_NN
        #print S.H_3[ilayer][0].data[:10]
        S.H_3[ilayer][0] = data_bac[ilayer]
        S.H_3[ilayer_p1][0] = data_bac[ilayer_p1]
        temp2 = S.H_3[ilayer][0].data[:10].copy()
        if not np.all(temp1==temp2):
            raise Exception("H_3 not succesfully restored")

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
        #calc_scaling_dim_1site(fn, iterate=False, local=True, qn_list=qn_list)
        res = calc_scaling_dim_3site(fn, local=True, S=None, qn_list=qn_list, tol=1e-8, k=12); print res
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
    qn_list = None #[0, 1, 2]#[0, 1, 2, 3]

   
    if 0:
        pdb = PDB(fn=fn, path=path)
        pdb.plot_scaling_dim()
        pdb.plot_scaling_dim(calc_err=True)
        #pdb.plot_energy()
        #pdb.plot_energy_diff()
    else:
        os.environ["OMP_NUM_THREADS"] = str(4)
        calc_scaling_dim_1site(fn, iterate=False, local=True, qn_list=qn_list)
        #calc_scaling_dim_1site(fn, iterate=False, local=False, qn_list=qn_list)
        #calc_scaling_dim_1site(fn, iterate=True, local=True, qn_list=qn_list)
        #res = calc_scaling_dim_3site(fn, local=True, S=None, qn_list=qn_list, tol=1e-8, k=12)
        #for i in res:
        #    print i
        #print res
        res = calc_scaling_dim_3site(fn, local=True, S=None, qn_list=qn_list, tol=1e-8, k=12); print res
        #for i in res:
        #    print i
        




