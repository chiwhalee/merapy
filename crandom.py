from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import *
import platform
import socket
import sys 

arch=  platform.architecture()[0]   #32 or 64 
hostname = socket.gethostname()
os1 = platform.system()

is_py3 = sys.version_info>(3, 0)

if hostname == 'VirtualBox-Lab':
    from merapy.lib.random_64_ifort_virtual import mrandom
else:
    if os1 == 'Linux':
        if is_py3:
            try:
                from merapy.lib.linux_py3.random2 import mrandom
            except:
                from merapy.lib.linux_py3.random2_gfort import mrandom
                
        else:
            from merapy.lib.linux_py2.random_64_ifort import mrandom
            
    elif os1 == 'Windows':
        from merapy.lib.win.random_gfort import mrandom
    

rand = mrandom.crand

