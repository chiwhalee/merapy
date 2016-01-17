import platform
import socket


arch=  platform.architecture()[0]   #32 or 64 
hostname = socket.gethostname()
os1 = platform.system()

if hostname == 'VirtualBox-Lab':
    from merapy.lib.random_64_ifort_virtual import mrandom
else:
    if os1 == 'Linux':
        from merapy.lib.random_64_ifort import mrandom
    elif os1 == 'Windows':
        from merapy.lib.win.random_gfort import mrandom
    

rand = mrandom.crand

