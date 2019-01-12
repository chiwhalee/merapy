from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import



#import system_parameter
#which = system_parameter.SystemParam.lang_tensor
from future import standard_library
standard_library.install_aliases()
from builtins import *
which = 'py'


if which == "py":
    #from .tensor_py import *
    from merapy.tensor_py import *
elif which == "cython":
    print("using cython implementation of iTensor")
    import pyximport 
    pyximport.install()  # above two lines are necessary, but why?

    from tensor_pyx import *
    #from tensor_pyx import test_iTensor
else:
    print("wrong module to import, see tensor.py")
    exit()


