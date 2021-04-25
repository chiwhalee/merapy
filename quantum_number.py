from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

#from system_parameter import lang_tensor

from future import standard_library
standard_library.install_aliases()
from builtins import *
#which = 'mypy'
which = 'py'

if which == "py":
    from merapy.quantum_number_py import *
elif which == 'mypy':
    from merapy.lib.quantum_number_mypy import *
    
elif which == "cython":
    print("using cython implementation of quantum_number"*10)
    #import pyximport 
    #pyximport.install()  # above two lines are necessary, but why?

    from quantum_number_pyx  import *
    #from tensor_pyx import test_iTensor
else:
    print("wrong module to import, see quantum_number.py")
    exit()

__all__ = [
        'QspZ2', 'QspU1', 'QspTravial', 'QnZ2', 'QnU1', 'make_qsp', 'qsp_any', 
        'symmetry_to_Qn', 'symmetry_to_Qsp', 
        ]

