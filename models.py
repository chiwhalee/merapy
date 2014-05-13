import numpy as np

from merapy.tensor_py import iTensor
from merapy.quantum_number import *

__all__ = ['potts']


class Model(object):
    def __init__(self):
        pass
        self.energy_exact = None

def potts(**kwargs):
    """
        the ham H = sum  {  Z - X_Xd - Xd_X  }
        so 
            h2 = 0.5*(ZI + IZ) - X_Xd - Xd_X
    """
    symmetry = kwargs['symmetry']
    h = kwargs['h'];  J_NNN = kwargs['J_NNN']
    qsp_base = kwargs['qsp_base']
    if symmetry == "Travial":      
        #_Z = np.array([[-1, 0, 0], [0, 2, 0], [0, 0, -1]], dtype=np.float64)
        _Z = np.array([[2, 0, 0], [0, -1, 0], [0, 0, -1]], dtype=np.float64)
        _X = np.array([[ 0, 1, 0], [0, 0, 1], [1, 0,  0]], dtype=np.float64)
        _I = np.eye(3)
        tqn = qsp_base.QnClass.qn_id()        
        Z = iTensor(2, qsp_base.copy_many(2), tqn)
        X = iTensor(2, qsp_base.copy_many(2), tqn)
        I = iTensor(2, qsp_base.copy_many(2), tqn)

        #qsp = QspTravial.easy_init(qns=[1], dims=[3])
        qsp = qsp_base.copy()
        tqn = QnTravial.qn_id()
        Z.data[:] = _Z.ravel()
        X.data[:] = _X.ravel()
        X_dag = X.conjugate(1)
        I.data[:] = _I.ravel()
    
    elif symmetry == 'Z3':
        tqn=qsp_base.QnClass.qn_id()
        Z = iTensor(2, qsp_base.copy_many(2, reverse=[1]), tqn)
        Z.data[:] = [2.0, -1.0, -1.0]
        #print Z
        X = iTensor(2, qsp_base.copy_many(2, reverse=[1]), QnZ3(1))
        X.data[:] = 1.0
        X_dag = X.conjugate(1)
        I = iTensor.unit_tensor(2, qsp_base.copy_many(2, reverse=[1]))

    #h2 = None
    h0 = {}
    if 1:
        X_Xd = X.direct_product(X_dag)
        Xd_X = X_dag.direct_product(X)
        ZI = Z.direct_product(I)
        IZ = I.direct_product(Z)
        
        #h2 = 0.5*(ZI + IZ) - X_Xd - Xd_X
        h0[2] = -0.5*(ZI + IZ) - X_Xd - Xd_X
        h0[2] +=  -1.0*I.direct_product(I)
        

        #h0[3] = X.direct_product(I).direct_product(X_dag) + X_dag.direct_product(I).direct_product(X)
        #h0[3].data *= J_NNN 
        h0[3] = I.direct_product(I).direct_product(I)  # it is better to leave it to zero
        h0[3].data[:] = 0.0
    
    if not h0[2].is_hermite(): raise

    return h0




if __name__ == '__main__':
    pass
    from merapy.quantum_number_py import *
    
    #potts(symmetry='Travial', qsp_base=QspTravial.easy_init(qns=[1], dims=[3]), h=1.0)
    qsp_base=QspZ3.easy_init(qns=[0, 1, 2], dims=[1, 1, 1])

    h0=potts(symmetry='Z3', qsp_base=qsp_base, h=1.0, J_NNN=0.0)
    print h0[2].data




    
    

    


