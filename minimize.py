#!/usr/bin/env python
#coding=utf8
import unittest
import time
from math import log10
import warnings
import time
from scipy.sparse.linalg.eigen.arpack.arpack import ArpackNoConvergence

from merapy import print_vars
from merapy.decorators import set_STATE_end_1, tensor_player
from merapy.measure_and_analysis.scaling_dimension import  calc_scaling_dim_1site, calc_scaling_dim_2site 
from merapy.measure_and_analysis.measurement import measure_S
from merapy.top_level import top_level_product_state, top_level_product_state_u1
from merapy.hamiltonian import System


class Iteration(object): 
    def __init__(self, S, updaters, q_iter, q_one, q_lay, filename="auto", backup_fn='auto', 
            run_time_max=None, energy_diff_min=None, q_iter_relative=False, use_player=False, 
            info=0, **kwargs):
        
        #self.updaters= updaters
        if isinstance(updaters, list):  #backward compatable
            self.ascending_ham, self.descending_ham, self.iterative_optimize_all = updaters[:3]
        elif isinstance(updaters, dict):
            xx = ["ascending_func", "descending_func", "update_mera_func"]
            yy = ["ascending_ham", "descending_ham", "iterative_optimize_all"]
            for i in range(3):
                setattr(self, yy[i], updaters[xx[i]])

        self.q_iter = q_iter
        self.q_one = q_one
        self.q_lay = q_lay
        self.use_player = use_player
        self.run_time_max = run_time_max if run_time_max is not None else 86400*20  # unit is second
        self.energy_diff_min = energy_diff_min if energy_diff_min is not None else 0.0
        self.info = info
        
        self.filename = filename
        self.backup_fn = backup_fn
            
        self.resume = resume
        self.auto_resume = auto_resume if auto_resume is not None else True
        
        self.calc_scaling_dim = calc_scaling_dim
        self.q_iter_relative = q_iter_relative
        self.__dict__.update(**kwargs)
   
    def measure_and_control(self, S, iter, iter0,  pace):
        if S.iter%pace == 0 or iter-iter0<pace: 
            
            self.eval_energy(iter, iter0, pace)

        if S.iter%200 == 0 and self.calc_scaling_dim:    #measure scaling dim
            calc_scaling_dim_1site(S=S, iterate=False)
            #res = calc_scaling_dim_2site(ascending_func=self.ascending_ham, S=S, k=10, tol=1e-6)
            try:
                res = calc_scaling_dim_2site(ascending_func=self.ascending_ham, S=S, k=10, tol=1e-6)
                for i in res: print i, res[i];  
                print '\n'
                S.scaling_dim_record[S.iter] = res
            except Exception as err:
                print err
            
        if self.auto_resume and self.backup_fn is not None:
            self.resume_func(S, S.mera)
    
    def resume_func(self, S, M):
        try:
            S0= System.load(self.backup_fn)
        except IOError:
            print "try to resume but file %s not found"%self.backup_fn
        else:
            e = S.energy; e0 = S0.energy
            if e> e0:
                S.__dict__ = S0.__dict__
                M.__dict__ = S0.mera.__dict__
                self.S = S;  self.S.mera = M
                
                S.iter = S0.iter
                S.iter1 = S0.iter1
                print "\nAOTU_REASSUMED FROM %s AT A LOWER ENEGY.\nstored E=%2.15f, current E=%2.15f"%(
                    self.backup_fn, e0, e)
        finally:
            self.auto_resume = False
   
if 0: 
    def Rand_Seq():
        iSeq = -1

        V=int(np.random.rand()/0.25)+1
        iSeq[1] = V
        pp=2
        while pp<5: 
            V=int(np.random.rand()/0.25) + 1
            for qq  in range(pp-1):
                if V == iSeq[qq]:  
                    break
            if qq == pp:  
                iSeq[pp] = V
                pp = pp+1
   

    def finite_site_u1_new(M,S, ascending, descending, update_mera, rho_top_func, 
            q_one, q_lay, q_iter, use_player=False, filename=None, 
            lstart=1, lstep=0, resume=False, info=0):
        """
        similar to finite_site_u1 below but changed update order: calculate all rho first
        the result seems correct, but need more tests, ESPECIALLY AGAINST TOP_LEVEL_EIGENSTATE
        """
        raise NotImplemented

        #err_tol = 1.E-7
        
        iter0 = S.iter
        
        #note here start is no longer controled by the passed in param start; remove that param in future
        if iter0  == 0:
            if rho_top_func.__name__ in ["top_level_product_state", "top_level_product_state_u1"]:
                rho_top_func(M, S)     #rho_top is never changed when using top_level_product_state
        
        for iter  in range(iter0, q_iter+iter0):
            if info>0:  
                print str(iter) + "th iteratiion"
            
            #update rho_2[M.num_of_layer-1]        
            if rho_top_func.__name__ == "top_level_eigenstate": 
                rho_top_func(M, S)

            #这里似乎不一定在top layer收缩，可以任意layer
            S.eng_ham(M.num_of_layer-1)
            S.display(iter, filename)
            
            #这一步更新rho[num_of_layer-2] ,..., [1]
            for ilayer  in range(M.num_of_layer-1, lstart + 1, -1):
                #rho[0]可以求出来，但它在更新U，V时用不到，我们也不用它来收缩算能量，故不求它
                for iter_lay in range(q_lay):
                    #calculate reduced density mat. for each layer
                    descending(M,S,ilayer,  info=info-1)

            #第一遍：从下到上把H_2求出来
            #注意对于finite range，在fortran的算法里，最上层的V张量是个摆设，不是实际的最上层!, 它的值没有被更新过
            for ilayer in range(lstart, M.num_of_layer-1):
                #note the update order here, start from H then M then H:  H[0]->M[0]->H[1]->M[1]->....
                #H[0] is inited before finite_site
                for iter_lay in range(q_lay):
                    #iterative_optimize_all(M, S, ilayer,j=0, info=info-1)
                    update_mera(M,S,ilayer,j=0, info=info-1)    #update U[ilayer] etc
                ascending(M, S, ilayer, info=info-1)   #calculate H[ilayer+1] for ilayer
            

            #start =  False 
            #ilayer = M.num_of_layer
            
            if iter%100 == 0 and iter>iter0:
                if backup_fn is not None:
                    S.save(fn=backup_fn)

            S.iter += 1 
            if use_player:
                set_STATE_end(iter, q_iter, iter0=iter0, resume=resume, power_on=True)

def finite_site(M, S, ascending, descending, update_mera,  rho_top_func,
        q_one, q_lay, q_iter,  use_player=False, 
        filename=None, backup_fn=None, resume=False, 
        auto_resume=True,  lstart=1, lstep=0, start=True, info=0, **kwargs):
    """
        总的来说，要更新h2, rho, mera 这三部分，要合理安排三类算子，以及每类算子的每一曾的更新次序。
        更新h2, 必须要有mera；
        更新mera，必须要有h2和rho；
        更新rho，实际上rho不是独立的，完全由mera决定
        只有h2[0]和mera是实际的矩阵，而其他都是收缩的中间过程, 是虚的

        conduct a complete optimization for the system
        注意对于finite range，在fortran的算法里，最上层的V张量是个摆设，不是实际的最上层!, 它的值没有被更新过

        q_iter: wee vidal p.18 qiter, number to sweep  optimizing all u and w in whole mera
        q_lay: see vidal p.18 qlay, the number to repeat optimizing all u and w within that layer
        q_one: see vidal p.16 qone, number to iterate to optimize a single w or u
    """
    #include "interfaces.f90"
    
    #    Names[0:4] = ["U", "W1", "W2", "V"]
    Names= ["U", "V"]
    

    if info>0:
        print "START finite_site"

    err_tol = 1.E-7
    
    # iter=0 is init, original value
    #S.iter precisely means the state of S after iter times update, e.g. S.iter = 0 means original state
    iter0 = S.iter 
    if iter0  == 1:
        if rho_top_func.__name__ in ["top_level_product_state", "top_level_product_state_u1"]:
            rho_top_func(M, S)     #rho_top is never changed when using top_level_product_state
        start = True 
    elif iter0 > 1:
        start = False

    # tho rho is never changed
    rho_top_func(M, S)
    
    for iter  in range(iter0 + 1, q_iter+iter0 + 1):   #for iter  in range(q_iter):
        
        if info>0:
            print str(iter) + "th iteratiion"
        if use_player:
            if iter0 == 0:
                set_STATE_end_1(iter, record_at=iter0 + 2, power_on=True)
            else:
                if resume:
                    set_STATE_end_1(iter, record_at=iter0 + 1, power_on=True)
                else:
                    set_STATE_end_1(iter, record_at=-1, power_on=True)  #no more record, start play immediatly

        #第一遍：从下到上把H_2求出来
        for ilayer in range(lstart, M.num_of_layer-1):
            #注意对于finite range，在fortran的算法里，最上层的V张量是个摆设，不是实际的最上层!, 它的值没有被更新过
            if info>1:
                print "\nupdating layer", ilayer
                
            if start:   
                ascending(M,S,ilayer,info=info-1)
            else:
                for iter_lay in range(q_lay):
                    update_mera(M,S,ilayer,j=0, info=info-1)
        
        #rho_top_func(M, S)

        for ilayer  in range(M.num_of_layer-2, lstart-1, -1):
            if info>1:
                print "updating layer", ilayer
            for iter_lay in range(q_lay):
                update_mera(M,S,ilayer,j=0,  info=info-1)

        start =  False 
        #S.iter = iter
        S.iter += 1 
        #print "iiii", S.iter 
        S.eng_ham(M.num_of_layer-1)
        S.save_eng(S.iter)
        S.display(S.iter, filename)
        
    if info>0:
        print "END finite_site"

def finite_site_u1(M, S, ascending, descending, update_mera, rho_top_func, 
        q_one, q_lay, q_iter, q_iter_relative=False, use_player=False, 
        filename=None, backup_fn=None, 
        lstart=1, lstep=0, resume=False, auto_resume=True, before_scaleinvar=False, 
        run_time_max=None, energy_diff_min=0.0, info=0, **kwargs):
    """
        ref:
            参考 Evenbly, G., and G. Vidal. "Algorithms for entanglement renormalization." Physical Review B 79.14 (2009): 144108.
            p.17 A1-A4.  
            与Vidal文中稍有不同的是，see my note n5 at page 17
        todo:
            应该还用Vidal文中的迭代顺序，即先把所有的rho求出来，这样更舒服


        总的来说，要更新h2, rho, mera 这三部分，要合理安排三类算子，以及每类算子的每一曾的更新次序。
        更新h2, 必须要有mera；
        更新mera，必须要有h2和rho；
        更新rho，实际上rho不是独立的，完全由mera决定
        只有h2[0]和mera是实际的矩阵，而其他都是收缩的中间过程, 是虚的
        conduct a complete optimization for the system

        q_iter: wee vidal p.18 qiter, number to sweep  optimizing all u and w in whole mera
        q_lay: see vidal p.18 qlay, the number to repeat optimizing all u and w within that layer
        q_one: see vidal p.16 qone, number to iterate to optimize a single w or u

    """

    use_local_storage = kwargs.get('use_local_storage', False)
    last_save_time = time.time()
    if resume:    
        auto_resume = False

    backup_fn_local = kwargs.get('backup_fn_local', None)
    iter0 = S.iter1    #S.iter precisely means the state of S after iter times update, e.g. S.iter = 0 means original state

    
    if S.iter == 0:  #note here start is no longer controled by the passed in param start; remove that param in future
        if rho_top_func.__name__ in ["top_level_product_state", "top_level_product_state_u1", "top_level_product_state_u1_1site"]:
            rho_top_func(M, S)     #rho_top is never changed when using top_level_product_state
        start = True 
    else:  
        start = False
    
    if q_iter_relative:   
        q_iter = q_iter + iter0 + 1
    q_iter = q_iter + iter0 + 1 if q_iter_relative else  q_iter + 1
    is_resumed = False

    for iter  in xrange(iter0 + 1, q_iter):
        if info>0:
            print str(iter) + "th iteratiion"
        
        if use_player:
            verbose = 1 if iter - iter0< 10 else 0
            if iter0 == 0:
                set_STATE_end_1(iter, record_at=iter0 + 2, power_on=True, verbose=verbose)
            else:
                if resume:
                    set_STATE_end_1(iter, record_at=iter0 + 1, power_on=True, verbose=verbose)
                else:
                    set_STATE_end_1(iter, record_at=-1, power_on=True, verbose=verbose)  #no more record, start play immediatly
        
        if S.use_pinning_term:
            S.use_pinning_term_func()
            
        
        #第一遍：从下到上把H_2求出来
        #注意对于finite range，在fortran的算法里，最上层的V张量是个摆设，不是实际的最上层!, 它的值没有被更新过
        for ilayer in range(lstart, M.num_of_layer-1):
            #note the update order here, start from H then M then H:  H[0]->M[0]->H[1]->M[1]->....
            #H[0] is inited before finite_site
            if not start:   
                for iter_lay in range(q_lay):
                    #iterative_optimize_all(M,S,ilayer,j=0, info=info-1)
                    update_mera(M,S,ilayer,j=0, info=info-1)    #update U[ilayer] etc
            ascending(M, S, ilayer, info=info-1)   #calculate H[ilayer+1] for ilayer
        
        
        if rho_top_func.__name__ == "top_level_eigenstate":   #update rho_2[M.num_of_layer-1]        
            rho_top_func(M, S)

        #这里似乎不一定在top layer收缩，可以任意layer
        
        S.iter += 1 
        S.iter1 += 1 
        S.eng_ham(M.num_of_layer-1)
        S.save_eng(S.iter)
        steps, dt_real= S.display(S.iter, filename)

        if iter-iter0>10:
            if run_time_max is not None: 
                if steps*dt_real>run_time_max: #(10 days)
                    print "\nENERGY_DIFF IS TOO SMALL, BREAK AT HERE. steps= %(steps)d\n"%vars()
                    break
            if abs(S.energy_diff_std) < energy_diff_min:
                print "energy_diff lower than energy_diff_min = %1.2e. BREAK AT HERE."%energy_diff_min
                #self.is_break = True
                break
           
            if auto_resume: 
                path = backup_fn if not use_local_storage else backup_fn_local
                try:
                    S0= System.load(path, use_local_storage)
                except IOError:
                    print "try to resume but file %s not found"%path
                else:
                    e = S.energy; e0 = S0.energy
                    if e> e0:
                        #S= S0
                        #M = S0.mera
                        S.__dict__ = S0.__dict__
                        M.__dict__ = S0.mera.__dict__
                        S.iter = S0.iter
                        print "AOTU_REASSUMED FROM %s AT A LOWER ENERGY E=%2.10f LOWER THAN CURRENT E=%2.10f\n"%(
                                path, e0, e)
                finally:
                    auto_resume = False


        #这一步更新rho[num_of_layer-2] ,..., [1]
        for ilayer  in range(M.num_of_layer-1, lstart + 1, -1):
            #start from rho[M.num_of_layer-1] desendeing -> rho[M.num_of_layer-2]->...->rho[1] 
            #rho[0]可以求出来，但它在更新U，V时用不到，我们也不用它来收缩算能量，故不求它
            for iter_lay in range(q_lay):
                #calculate reduced density mat. for each layer
                descending(M,S,ilayer,  info=info-1)

        start =  False 
        
        if S.iter%50 == 0: 
            if time.time()-last_save_time>3600: 
                last_save_time = time.time()
                if backup_fn is not None:
                    if not use_local_storage: 
                        S.save(fn=backup_fn)
                    else: 
                        S.save(fn=backup_fn_local, use_local_storage=use_local_storage)
            
    iter_tot = S.iter1-iter0
    msg = '\ttotal number of iteration is %d'%(iter_tot); print msg
    if iter_tot>0 : 
        temp = backup_fn if not S.use_local_storage else backup_fn_local
        if temp is not None: 
            print "\nS has been changed, save at final iter"
            S.save(fn=temp, use_local_storage=S.use_local_storage)


def FiniteRangeIter(Iteration): 
    pass

class ScaleInvar():
    """
        tested by compareing  with Heisenberg/scaleinvariant.f90, the results coincide for all digits
    """
    def __init__(self, S, updaters, num_of_SIlayer, q_iter, q_one, q_lay, 
            resume=False, auto_resume=True, filename="auto", backup_fn=None, 
            run_time_max=None, energy_diff_min=None, q_iter_relative=True, use_player=False, 
            calc_scaling_dim=False, info=0, **kwargs):
        #print 'qqqqqqq', q_iter
        self.S= S
        #self.updaters= updaters
        if isinstance(updaters, list):  #backward compatable
            self.ascending_ham, self.descending_ham, self.iterative_optimize_all = updaters[:3]
        elif isinstance(updaters, dict):
            xx = ["ascending_func", "descending_func", "update_mera_func"]
            yy = ["ascending_ham", "descending_ham", "iterative_optimize_all"]
            for i in range(3):
                setattr(self, yy[i], updaters[xx[i]])

        self.num_of_SIlayer = num_of_SIlayer
        self.q_iter = q_iter
        self.q_one = q_one
        self.q_lay = q_lay
        self.use_player = use_player
        #self.run_time_max = run_time_max if run_time_max is not None else 86400*20  # unit is second
        self.run_time_max = run_time_max
        self.energy_diff_min = energy_diff_min if energy_diff_min is not None else 0.0
        self.info = info
        self.SIlayer = self.S.mera.num_of_layer-2
        self.filename = filename
        self.backup_fn = backup_fn
            
        self.resume = resume
        #self.auto_resume = auto_resume 
        self.is_resumed = False
        self.ham_op_names = S.fetch_op_names(type="H")
        self.calc_scaling_dim = calc_scaling_dim
        self.q_iter_relative = q_iter_relative
        temp = ['use_local_storage', 'backup_fn_local']
        for t in temp: 
            setattr(self, t, kwargs.get(t))
        
        self.__dict__.update(**kwargs)
        
    def scale_invariant(self):
        """
            for lay in [0, M.num_of_layer-2] update normally
            for lay in   [M.num_of_layer-2, M.num_of_layer-+ num_of_SIlayer] update scale invariantly 
        """
        self.is_break = False
        last_save_time = time.time()

        S= self.S
        M = S.mera
        info = self.info

        if self.filename == "auto":
            self.filename = self.output_fn_auto()
            print "output_fn is set to %s"%self.filename
        if 1:
            aa = "\nSTART SCALE_INVARIANT ITERATION\n"
            print aa 
            if self.filename is not None:
                out = open(self.filename, "a"); out.write(aa); out.close()
            #print "use following ham ops: %s\n"%self.ham_op_names.keys()
        
       
        iter0 = S.iter1
        
        q_iter = self.q_iter + iter0 + 1 if self.q_iter_relative else  self.q_iter + 1
        

        if not self.resume:
            S.log[iter0] = S.get_info()
        
        for iter in range(iter0 + 1, q_iter):  #for iter in range(iter0 + 1, q_iter + iter0 + 1):
            if self.use_player:
                record_at = iter0+2 if self.resume else iter0 + 1
                set_STATE_end_1(iter, record_at=record_at, power_on=True)
            
            if S.use_pinning_term:
                S.use_pinning_term_func()

            #for ilayer in range(M.num_of_layer-2):
            for ilayer in range(self.SIlayer):
                self.iterative_optimize_all(M, S, ilayer, j=0, info=info-1)
                self.ascending_ham(M, S, ilayer, info=info-1)

            self.optimize_SIlayer()
            S.iter += 1 
            S.iter1 += 1 
            
            self.measure_and_control_new(S, iter , iter0, pace=20)
            if self.is_break: break
                
            #for ilayer in range(M.num_of_layer-2, 1, -1):
            for ilayer in range(self.SIlayer, 1, -1):
                self.descending_ham(M, S, ilayer)
            
            if S.iter%50 == 0 and iter-iter0>10:
                if time.time()-last_save_time>3600: 
                    last_save_time = time.time()
                    temp = self.backup_fn if not self.use_local_storage else self.backup_fn_local
                    if temp is not None: 
                        S.save(fn=temp, use_local_storage=self.use_local_storage)
        
        iter_tot = S.iter1-iter0
        msg = '\ttotal number of iteration is %d'%(iter_tot); print msg
        if iter_tot>0 : 
            temp = self.backup_fn if not self.use_local_storage else self.backup_fn_local
            if temp is not None: 
                print "\nS has been changed, save at final iter"
                S.save(fn=temp, use_local_storage=self.use_local_storage)
            

    if 1: pass
        #def scale_invariant_new(self):
        #    """
        #    combine scale_invariant and optimize_SIlayer
        #        for lay in [0, M.num_of_layer-2] update normally
        #        for lay in   [M.num_of_layer-2, M.num_of_layer-+ num_of_SIlayer] update scale invariantly 
        #    """
        #    if self.resume:    
        #        auto_resume = False
        #    self.is_break = False

        #    S= self.S
        #    M = S.mera
        #    info = self.info
        #    q_iter = self.q_iter


        #    if self.filename == "auto":
        #        self.filename = self.output_fn_auto()
        #        print "output_fn is set to %s"%self.filename
        #    if 1:
        #        aa = "\nSTART SCALE_INVARIANT ITERATION\n"
        #        print aa 
        #        if self.filename is not None:
        #            out = open(self.filename, "a"); out.write(aa); out.close()
        #        print "use following ham ops: %s\n"%self.ham_op_names.keys()
        #    
        #    #iter0 = S.iter
        #    iter0 = S.iter1
        #    q_iter = self.q_iter + iter0 + 1 if self.q_iter_relative else  self.q_iter + 1
        #    
        #    SIlayer = self.SIlayer   # = 2
        #    SIlayer_m1 = SIlayer - 1   
        #    SIlayer_p1 = SIlayer + 1   

        #    if not self.resume:
        #        S.log[iter0] = S.get_info()
        #    #for iter in range(iter0 + 1, q_iter + iter0 + 1):
        #    for iter in range(iter0 + 1, q_iter):
        #        if self.use_player:
        #            record_at = iter0+2 if self.resume else iter0 + 1
        #            set_STATE_end_1(iter, record_at=record_at, power_on=True)
        #        
        #        S.rho_2[SIlayer_p1][0].data[:] = S.rho_2[SIlayer][0].data[:]
        #        if not S.only_NN:  S.rho_3[SIlayer_p1][0].data[:] = S.rho_3[SIlayer][0].data[:]
        #        for ilayer in range(self.SIlayer, 1, -1):
        #            self.descending_ham(M, S, ilayer)

        #        for ilayer in range(self.SIlayer):
        #            self.iterative_optimize_all(M, S, ilayer, j=0, info=info-1)
        #            self.ascending_ham(M, S, ilayer, info=info-1)
        #        
        #        if 1:
        #            h_bac = {}
        #            h_aver = {}
        #            for o in self.ham_op_names:
        #                h_bac[o]= S.__getattribute__(o)[SIlayer][0].copy(use_buf=True)
        #                h_aver[o]= S.__getattribute__(o)[SIlayer][0].copy(use_buf=True)

        #            for i in range(self.num_of_SIlayer):  #calculate \bar h
        #                if i>0:
        #                    for o in self.ham_op_names:
        #                        S.__getattribute__(o)[SIlayer][0].data[:] = S.__getattribute__(o)[SIlayer_p1][0].data[:]
        #                self.ascending_ham(M, S, SIlayer, tau=SIlayer + i)
        #                x = 1./float(S.mera.period)**(i + 1)

        #                for o in self.ham_op_names:
        #                    h_aver[o].data += x*S.__getattribute__(o)[SIlayer_p1][0].data
        #            for o in self.ham_op_names:
        #                S.__getattribute__(o)[SIlayer][0].data[:] = h_aver[o].data[:]
        #            
        #            self.iterative_optimize_all(M, S, SIlayer, tau=SIlayer, j=0, info=info-1)    #upate mera at SIlayer

        #        self.optimize_SIlayer()
        #        S.iter += 1 
        #        S.iter1 += 1 
        #        
        #        #self.measure_and_control(S, iter, iter0)
        #        self.measure_and_control_new(S, iter , iter0, pace=20)
        #        if self.is_break: break
        #            
        #        
        #        if S.iter%50 == 0 and iter-iter0>10:
        #            if self.backup_fn is not None:
        #                S.save(fn=self.backup_fn)
        #    
        #    if self.backup_fn is not None: #always save after exit 
        #        print "\nS has been changed, save at final iter..."
        #        S.save(fn=self.backup_fn)

        #def optimize_SIlayer_old(self):
        #    S= self.S
        #    M = S.mera
        #    info = self.info
        #    
        #    ilayer = self.SIlayer   # = 2
        #    ilayer_m1 = ilayer - 1   
        #    ilayer_p1 = ilayer + 1   
        #    
        #    #backup mera at SIlayer
        #    h_bac = {}
        #    for o in self.ham_op_names:
        #        h_bac[o]= S.__getattribute__(o)[ilayer][0].copy(use_buf=True)
        #    
        #    S.rho_2[ilayer_p1][0].data[:] = S.rho_2[ilayer][0].data[:]
        #    if not S.only_NN:
        #        #S.rho_3[ilayer_p1][0] = S.rho_3[ilayer][0].copy(use_buf=True)
        #        S.rho_3[ilayer_p1][0].data[:] = S.rho_3[ilayer][0].data[:]
        #    #note the following are removed, this differs with Heisenberg/scaleInvariant.f90
        #    
        #    if 0:
        #        #calculate \bar h
        #        tau = ilayer - 2
        #        x = 1.0/S.mera.period

        #        for i in range(self.num_of_SIlayer):
        #            tau += 1
        #            self.ascending_ham(M, S, ilayer, tau=tau)
        #            S.H_2[ilayer][0].data += x*S.H_2[ilayer_p1][0].data 

        #    #update SI layers of U, V 
        #    tau = ilayer 
        #    for iter_layer in range(self.q_lay):
        #        
        #        #update ham ops
        #        if iter_layer>0:
        #            for o in self.ham_op_names:
        #                #S.H_2[ilayer][0].data[:] = h_bac["H_2"].data[:]  #this sentens has no effect when q_lay = 1
        #                S.__getattribute__(o)[ilayer][0].data[:] = h_bac[o].data[:]
        #        self.ascending_ham(M, S, ilayer, tau=tau)
        #        #following two line see ref1
        #        x = 1.0/S.mera.period
        #        for o in self.ham_op_names:
        #            #下面一行 + 的是h_bac, 而非 H[ilayer], 是考虑到q_lay>1的情况
        #            S.__getattribute__(o)[ilayer][0].data  = x*S.__getattribute__(o)[ilayer_p1][0].data + h_bac[o].data 
        #        
        #        #upate mera
        #        self.iterative_optimize_all(M, S, ilayer, j=0, info=info-1)
        #        
        #        #update rho2, rho3 at both SIlayer and SIlayer+1
        #        self.descending_ham(M, S, ilayer_p1, info=info-1)
        #        if 1:
        #            S.rho_2[ilayer_p1][0].data[:] = S.rho_2[ilayer][0].data[:]
        #            if not S.only_NN:
        #                S.rho_3[ilayer_p1][0].data[:] = S.rho_3[ilayer][0].data[:]

        #    for o in self.ham_op_names:
        #        S.__getattribute__(o)[ilayer][0].data[:] = h_bac[o].data[:]

    def optimize_SIlayer(self):
        """ 
            improve: in future optimize_SIlayer should be combined with scale_invariant, this is easy, just put descending_ham before ascending_ham
            taken SIlayer as scale invar layer
            更新U, V etc in scale invar layer, to update them, we need
            h_diamond
            tau: tau is used only when calculate coup coeff for long range models
            when num_of_layer = 4, SIlayer = 2
            taken SIlayer and above as SI layer
            ref1:
                see P.36 eq.(A32) in Evenbly, Glen, and Guifre Vidal. 
                "Quantum Criticality with the Multi-scale Entanglement Renormalization Ansatz." 
                arXiv preprint arXiv:1109.5334 (2011). 
        """
        S= self.S;   M = S.mera;     info = self.info
        
        SIlayer = self.SIlayer   # = 2
        SIlayer_m1 = SIlayer - 1   
        SIlayer_p1 = SIlayer + 1   
        
        #below optimize SIlayer. first we need rho[SIlayer + 1], simply copy from sublayer
        #since SI layers are SI, rho_2[tau>=SIlayer] are all the same, so just copy it from SIlayer to SIlayer_p1
        S.rho_2[SIlayer_p1][0].data[:] = S.rho_2[SIlayer][0].data[:]
        if not S.only_NN:  S.rho_3[SIlayer_p1][0].data[:] = S.rho_3[SIlayer][0].data[:]
            
        
        #backup mera at SIlayer
        h_bac = {}
        h_aver = {}
        for o in self.ham_op_names:
            h_bac[o]= S.__getattribute__(o)[SIlayer][0].copy(use_buf=True)
            h_aver[o]= S.__getattribute__(o)[SIlayer][0].copy(use_buf=True)

        if 1:
            for i in range(self.num_of_SIlayer):  #calculate \bar h
                if i>0:
                    for o in self.ham_op_names:
                        S.__getattribute__(o)[SIlayer][0].data[:] = S.__getattribute__(o)[SIlayer_p1][0].data[:]
                self.ascending_ham(M, S, SIlayer, tau=SIlayer + i)
                x = 1./float(S.mera.period)**(i + 1)

                for o in self.ham_op_names:
                    h_aver[o].data += x*S.__getattribute__(o)[SIlayer_p1][0].data
            for o in self.ham_op_names:
                S.__getattribute__(o)[SIlayer][0].data[:] = h_aver[o].data[:]
        else:
            tau = SIlayer 
            self.ascending_ham(M, S, SIlayer, tau=tau)
            #following two line see ref1
            x = 1.0/S.mera.period
            for o in self.ham_op_names:
                S.__getattribute__(o)[SIlayer][0].data   += x*S.__getattribute__(o)[SIlayer_p1][0].data #+ h_bac[o].data 
        
        self.iterative_optimize_all(M, S, SIlayer, tau=SIlayer, j=0, info=info-1)    #upate mera at SIlayer
        self.descending_ham(M, S, SIlayer_p1, info=info-1)   #update rho2, rho3 at SIlayer 
        
        #why restore H? because what happened above SIlayer should not affect H[for all lay<SIlayer]
        #这一步不完全必要，因为之后ascending会自动覆盖现在的H, 下面两句可能只是为了接着eval_energy or eng_ham 提供正确H[SIlayer]
        for o in self.ham_op_names:
            S.__getattribute__(o)[SIlayer][0].data[:] = h_bac[o].data[:]

    def _eval_energy(self):
        """
            evaluate energy at some si layer
            the method is meant to eval a lower ground energy, not optimeze mera
                tau:  used only when calculate coup coeff for long range models
        """
        S = self.S
        M = S.mera
        SIlayer = self.SIlayer   # = 2
        SIlayer_p1 = SIlayer + 1

        for i in range(self.num_of_SIlayer):
        #这里是反复地从 SIlayer ascending 到 SIlayer_p1,  WITH tau changing
            if i>0:
                for o in self.ham_op_names:
                    S.__getattribute__(o)[SIlayer][0].data[:] = S.__getattribute__(o)[SIlayer_p1][0].data[:]
            self.ascending_ham(M, S, SIlayer, tau=SIlayer + i)
        #above only update H, no need update rho, as the later is scale-invar not changed in each layer
        
        if 1:
            S.rho_2[SIlayer_p1][0].data[:] = S.rho_2[SIlayer][0].data[:]
            if not S.only_NN:
                S.rho_3[SIlayer_p1][0].data[:] = S.rho_3[SIlayer][0].data[:]
        #in the end eval eng at the layer above the first SIlayer
        S.eng_ham(SIlayer_p1)

    def eval_energy(self, iter, iter0, pace):
        """
            eval_energy with its own tensor player
        
        """
        S= self.S
        state_bac = tensor_player.STATE
        meth_names = ["set_data_entrance", "contract_core", "permutation"]
        iTensor = S.H_2[0][0].__class__   # t is just used for passing one hook in
        if not hasattr(self, 'meth_deced'):
            self.meth_deced = {}
            self.meth_bac = {}
            for i in meth_names:
                self.meth_bac[i] = iTensor.__dict__[i] #backup original method
                self.meth_deced[i] = tensor_player(which=i)(iTensor.__dict__[i])
        for i in meth_names:
            setattr(iTensor, i, self.meth_deced[i]) 

        tensor_player.STATE = 'record' if iter == iter0+1 else 'play'

        #print 'iii', iter, tensor_player.STATE
        self._eval_energy()

        for i in self.meth_bac:
            setattr(iTensor, i, self.meth_bac[i])
        tensor_player.STATE = state_bac

        if 1:
            #S.eng_ham(ilayer=self.SIlayer)
            S.save_eng(S.iter)
            if iter-iter0<pace:   pace = 1
        
            steps, dt_real= S.display(S.iter, self.filename, pace=pace)
            #is_break = False
            if iter-iter0>10:
                if self.run_time_max is not None : 
                    if steps*dt_real>self.run_time_max : #(1days*?)
                        print "\nENERGY_DIFF IS TOO SMALL, BREAK AT HERE. total steps run = ", iter-iter0
                        self.is_break = True
                #if abs(S.energy_diff) < self.energy_diff_min:
                if abs(S.energy_diff_std) < self.energy_diff_min:
                    print "energy_diff lower than energy_diff_min = %1.2e. BREAK AT HERE."%self.energy_diff_min
                    self.is_break = True
            
    def resume_func(self, S, M, use_local_storage):
        path = self.backup_fn_local if use_local_storage else self.backup_fn
        
        try:
            S0= System.load(path, use_local_storage=use_local_storage)
        except IOError:
            print "try to resume but file %s not found"%path
        else:
            e = S.energy; e0 = S0.energy
            if e> e0:
                S.__dict__ = S0.__dict__
                M.__dict__ = S0.mera.__dict__
                self.S = S;  self.S.mera = M
                
                S.iter = S0.iter
                S.iter1 = S0.iter1
                print "\nAOTU_RESUMED FROM %s AT A LOWER ENEGY.\nstored E=%2.15f, current E=%2.15f"%(
                    path, e0, e)
        finally:
            self.is_resumed = True
            pass
    
    def measure_and_control_new(self, S, iter, iter0,  pace):
        if S.iter%pace == 0 or iter-iter0<pace: 
            
            self.eval_energy(iter, iter0, pace)

        if S.iter%200 == 0 and self.calc_scaling_dim:    #measure scaling dim
            calc_scaling_dim_1site(S=S, iterate=False)
            #res = calc_scaling_dim_2site(ascending_func=self.ascending_ham, S=S, k=10, tol=1e-6)
            try:
                res = calc_scaling_dim_2site(ascending_func=self.ascending_ham, S=S, k=10, tol=1e-6)
                for i in res: print i, res[i];  
                print '\n'
                S.scaling_dim_record[S.iter] = res
            except Exception as err:
                print err
            
        #if self.auto_resume and self.backup_fn is not None:
        if (not self.is_resumed) and self.backup_fn is not None:
            self.resume_func(S, S.mera, self.use_local_storage)
    
    def _fn_auto(self):
        dim = self.S.mera.qsp_max.totDim
        nqn = self.S.mera.qsp_max.nQN
        nqn = "5_" if nqn == 5 else "" 
        lay = self.S.mera.num_of_layer-1
        lay = "-%dlay"%lay if lay != 3 else ""
        return dim, nqn, lay

    def output_fn_auto(self):
        dim, nqn, lay = self._fn_auto()
        fn = "res-%(nqn)s%(dim)d%(lay)s.dat"%vars()
        return fn



class TestScaleInvar(unittest.TestCase): 
    def setUp(self): 
        pass
    def test_it(self): 
        from ascending import ascending_ham
        from descending import descending_ham
        from iteration import iterative_optimize_all
        #from all_in_once import 
        from finite_site import finite_site_u1
        from top_level import top_level_product_state_u1

        from main import  Main

        import warnings
        warnings.filterwarnings('ignore') 

        updaters_u1= [ascending_ham, descending_ham, iterative_optimize_all, finite_site_u1, top_level_product_state_u1]
        updaters_u1 = {"ascending_func":ascending_ham, "descending_func":descending_ham, "update_mera_func":iterative_optimize_all, 
            "finite_range_func":finite_site_u1, "rho_top_func":top_level_product_state_u1}

        from decorators import profileit

        #main=profileit(fn="profile.tmp")(main)
        ising_model_param = {"h":-1.0, "J_NN":-1.0, "J_NNN":0.0}
        

        #def autoexe():
        if 1: 
            filename = "res-heisbg.dat"
            filename = None
            J_NNN = 0.241186
            #updaters_u1[-1] = top_level_eigenstate
            heisbg_model_param = {"J_NN":1.0, "J_NNN":J_NNN, 'Jzz': 1.0}
            main = Main(updaters=updaters_u1, trunc_dim=4, tot_layer=4, q_iter=3, q_lay=1, info=0, 
                    MODEL="Heisenberg", SYMMETRY="Travial", only_NN=True, only_NNN=False, USE_REFLECTION=False, 
                    filename=filename, NUM_OF_THREADS=1, use_player=False, combine_2site=False, 
                    resume=False,  
                    unitary_init="random_unit_tensor", model_param = heisbg_model_param, message = None)
            from merapy.utilities import random_str
            import os
            backup_fn = random_str( )
            S=main.run(q_iter=3, backup_fn=None, do_measure=0)
            S=main.run_scale_invar(q_iter=10,  backup_fn=backup_fn, filename=None, do_measure=0)
            res= {0: 0, 1: -0.14306831775876239, 2: -0.28302724810836422, 3: -0.37443146281546302, 4: -0.4519070431176464, 5: -0.52841563236836264, 6: -0.60836365005944959, 7: -0.69261165207729269, 8: -0.77380310057594426, 9: -0.84634118594060803, 10: -0.91177543745081024}
            S.examine_energy(res, delta=2e-14)
            
            print  'see if resumed successfully'
            S=main.run_scale_invar(q_iter=10,  backup_fn=backup_fn, filename=None, do_measure=0)
            
            os.remove(backup_fn)


if __name__ == "__main__":

    if 1: 
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           #'test_save', 
           #'test_load', 
           #'test_example',  
           #'test_heisbg',  
        ]
        for a in add_list: 
            suite.addTest(TestScaleInvar(a))
        unittest.TextTestRunner().run(suite)

       
