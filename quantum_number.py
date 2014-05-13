
#from system_parameter import lang_tensor
import system_parameter
which = system_parameter.SystemParam.lang_tensor

#which = "cython"

if which == "py":
    from quantum_number_py import *
elif which == "cython":
    print "using cython implementation of quantum_number"*10
    #import pyximport 
    #pyximport.install()  # above two lines are necessary, but why?

    from quantum_number_pyx  import *
    #from tensor_pyx import test_iTensor
else:
    print "wrong module to import, see quantum_number.py"
    exit()


