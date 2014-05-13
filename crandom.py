import platform
import socket
hostname = socket.gethostname()
arch=platform.architecture()[0]
if hostname == 'VirtualBox-Lab':
    from merapy.lib.random_64_ifort_virtual import mrandom
else:
    from merapy.lib.random_64_ifort import mrandom

rand = mrandom.crand

