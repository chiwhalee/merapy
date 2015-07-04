


#import system_parameter
#which = system_parameter.SystemParam.lang_tensor
which = 'py'


if which == "py":
    from tensor_py import *
elif which == "cython":
    print "using cython implementation of iTensor"
    import pyximport 
    pyximport.install()  # above two lines are necessary, but why?

    from tensor_pyx import *
    #from tensor_pyx import test_iTensor
else:
    print "wrong module to import, see tensor.py"
    exit()


