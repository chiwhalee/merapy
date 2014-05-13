
import unittest
import cPickle as pickle
from contextlib import contextmanager
import nose
import sys
import os
#from tempfile import mkdtemp
import tempfile
import unittest
import shutil
import rpyc
from plumbum import SshMachine
import rpyc
from rpyc.utils.zerodeploy import DeployedServer
import socket 


@contextmanager
def working_directory(path):
    current_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(current_dir)

if 0: 
    with working_directory("data/stuff"):
        # do something within data/stuff
        pass
    # here I am back again in the original working directory

@contextmanager
def redirect(**kwds):
    """
        example: 
            with redirect(stdout=open("log.txt", "a")): # these print statements will go to /tmp/log.txt 
                print "Test entry 1"
                print "Test entry 2"
            print "Back to normal stdout again"

    """    
    stream_names = ["stdin", "stdout", "stderr"]
    old_streams = {}
    try:
        for sname in stream_names:
            stream = kwds.get(sname, None)
            if stream is not None and stream != getattr(sys, sname):
                old_streams[sname] = getattr(sys, sname)
                setattr(sys, sname, stream)
        yield
    finally:
        for sname, stream in old_streams.iteritems():
            setattr(sys, sname, stream)

@contextmanager
def make_temp_dir(dir=None):
    temp_dir = tempfile.mkdtemp(dir=dir)
    try: 
        yield temp_dir
    finally: 
        shutil.rmtree(temp_dir)
   
@contextmanager
def make_temp_filename(suffix='', prefix='tmp', dir=None, delete=True):
    """This context creates a visible temporary file with size 0 in the file system.
    This file is removed by default when the context exits (delete = True).
    Unlike the tempfile module, this context only offers a filename (and not an open file) to
    the code that uses it.
    """
    fd, name = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir)
    try:
        os.close(fd)
        yield name
    finally:
        if delete:
            try:
                os.remove(name)
            except OSError:
                pass   


# back to the normal stdout
def get_ip_address(): 
    ip = socket.gethostbyname(socket.gethostname())
    #print 'ip is %s'%ip 
    
    try: 
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("gmail.com",80))
        #print('ip address is %s'%s.getsockname()[0])
        ip = s.getsockname()[0] 
        s.close()
    except: 
        print 'get ip address using gmail failed'
        
    return ip

def ssh_connect(hostname): 
    """
        a convenient func
    """
    
    if hostname in ['QTG-WS1-ubuntu', 'local']: 
        #ip, user = '210.45.117.30', 'zhli'
        args= dict(host='210.45.117.30', user='zhli', connect_timeout=100)
    elif hostname == 'sugon': 
        ip, user = '211.86.151.102', 'zhihuali'
        args= dict(host='211.86.151.102', user='zhihuali', connect_timeout=3600)
    else: 
        raise ValueError(hostname)
    
    ip = get_ip_address()
    msg = 'try ssh connecting from %s to %s'%(ip, args['host'])
    ssh = SshMachine(**args)
    msg +=   '  ... succeed'
    print(msg)
    return ssh

@contextmanager
def rpyc_conn_zerodeploy(hostname): 
    try:
        ssh = ssh_connect(hostname)
        server = DeployedServer(ssh)
        conn = server.classic_connect()
        path = ["/home/zhli/dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python"]
        conn.modules.sys.path.extend(path)
        # see this link, may be can find a way out
        #http://stackoverflow.com/questions/856116/changing-ld-library-path-at-runtime-for-ctypes
        
        yield conn
    finally: 
        conn.close()
        server.close()
        ssh.close()

@contextmanager
def rpyc_conn_local_zerodeploy(): 
    try: 
        get_ip_address()
        ssh = SshMachine('210.45.117.30', user='zhli', keyfile=None)  
        server = DeployedServer(ssh)
        conn = server.classic_connect()
        path = ["/home/zhli/dropbox/My-documents/My-code/quantum-many-body/mera-algorithms/python"]
        conn.modules.sys.path.extend(path)
        
        # search this https://www.google.com.hk/search?newwindow=1&safe=off&q=python+change+ld_library_path+at+runtime&oq=python+change+ld_library_path+at+runtime&gs_l=serp.3...3045.4184.0.4379.7.7.0.0.0.1.327.327.3-1.1.0....0...1c.1.39.serp..7.0.0.umAiWXHnYYs
        # see this link, may be can find a way out
        # http://stackoverflow.com/questions/1178094/change-current-process-environment
        #http://stackoverflow.com/questions/856116/changing-ld-library-path-at-runtime-for-ctypes
        if 0: 
            conn.modules.os.environ['VENV_LD_LIBRARY_PATH'] = os.environ['VENV_LD_LIBRARY_PATH']
            conn.modules.os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']
            conn.modules.os.environ['LIBRARY_PATH'] = os.environ['LIBRARY_PATH']
            conn.modules.os.environ['PATH'] = os.environ['PATH']
            conn.modules.os.environ['PYTHONPATH'] = os.environ['PYTHONPATH']
            conn.modules.os.environ['NLSPATH'] = os.environ['NLSPATH']
            lib =  ["/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64/libifport.so.5", 
                    "/opt/intel/composer_xe_2013.1.117/compiler/lib/mic/libimf.so"]
           
            conn.execute('import ctypes')
            for l in lib: 
                conn.execute('ctypes.cdll.LoadLibrary("%s")'%l)
            #conn.modules.os.environ.update(os.environ) 
        yield conn
    finally: 
        conn.close()
        server.close()
        ssh.close()

@contextmanager
def rpyc_conn(hostname): 
    """
        EOFError
    """
    try: 
        #ssh = SshMachine('210.45.117.30', user='zhli', keyfile=None)  
        ssh = ssh_connect(hostname)
        conn = rpyc.classic.ssh_connect(ssh, 17013)
        yield conn
    finally:
        try: 
            conn.close()
        except : 
            #raise Exception('connection failed')
            print 'rpyc connection failed'
            pass
        try: 
            ssh.close() 
        except: 
            print 'ssh connection failed'
            pass

@contextmanager
def rpyc_conn_local(): 
    try: 
        #EOFError:
        get_ip_address()
        
        ssh = SshMachine('210.45.117.30', user='zhli', keyfile=None)  
        conn = rpyc.classic.ssh_connect(ssh, 17013)
        yield conn
    finally:
        try: 
            conn.close()
        except : 
            #raise Exception('connection failed')
            print 'rpyc connection failed'
            pass
        try: 
            ssh.close() 
        except: 
            print 'ssh connection failed'
            pass

def rpyc_load(path, use_local_storage=False): 
    if not use_local_storage: 
        inn = open(path, "rb")
        res = pickle.load(inn)
        inn.close()
    else:
        if 0: # previouse trials and errors, keep this, from it can better understand rpyc
            if 0: # cannot load,  reason is that inn and load are not in the same namespace    
                with rpyc_conn_local() as conn: 
                    pass
                    print conn
                    inn = conn.builtin.open(path, "rb")
                    res = None
                    res = cPickle.load(inn)
                    inn.close()
                    res = conn.modules['merapy.hamiltonian'].System.load(path) 
            
            if 0: # can load,  but cannot return; reason is that, when conn closed,  S is not available
                with rpyc_conn_local() as conn: 
                    pass
                    res= None
                    path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/4.pickle' 
                    inn = conn.builtin.open(path, "rb")
                    res = conn.modules.cPickle.load(inn) 
                    #res = pickle.load(inn) 
                    res = conn.modules['merapy.hamiltonian'].System.update_S(res) 
                    print res
                    inn.close()
            if 0:  #this worked , but below is better 
                with rpyc_conn_local() as conn: 
                    #path = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/4.pickle' 
                    inn = conn.builtin.open(path, "rb")
                    res = conn.modules.cPickle.load(inn) 
                    res = conn.modules['merapy.hamiltonian'].System.update_S(res) 
                    res = rpyc.classic.obtain(res)
                    inn.close()
            
        with rpyc_conn_local() as conn: 
        #with rpyc_conn_local_zerodeploy() as conn:
            load_local = conn.modules['merapy.context_util'].rpyc_load
            res= load_local(path)
            res = rpyc.classic.obtain(res)
    return res

def rpyc_save(path, obj, use_local_storage=False): 
    if not use_local_storage: 
        out = open(path, "wb")
        pickle.dump(obj, out)
        out.close()
    else: 
        with rpyc_conn_local_zerodeploy() as conn:
            out = conn.builtin.open(path, 'wb')
            conn.modules.cPickle.dump(obj, out)
            out.close()

class DistCompute(object): 
    pass
    

class TestIt(unittest.TestCase): 
    def setUp(self): 
        pass
    def test_temp(self): 
           
        with make_temp_dir() as dd: 
            fn = dd + '/aaa'
            obj = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaa'
            rpyc_save(fn, obj, use_local_storage=1)
            a=rpyc_load(fn, use_local_storage=1)
            self.assertTrue(a==obj)
            
            fn ='/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0/4.pickle' 
            a=rpyc_load(fn, use_local_storage=1)
    
    def xtest_make_temp_dir(self): 
        with make_temp_dir() as dd: 
            print dd

    def test_redirect(self):
        with redirect(stdout=open("log.txt", "a")): # these print statements will go to /tmp/log.txt 
            print "Test entry 1"
            print "Test entry 2"
        print "Back to normal stdout again"
    
    def test_rpyc_conn_local(self): 
        with rpyc_conn_local() as conn: 
            print conn.modules['os']
            print conn.modules.merapy
    
    def test_rpyc_conn_local_zero(self): 
        with rpyc_conn_local_zerodeploy() as conn: 
            print conn.modules['os']

        with rpyc_conn_zerodeploy(hostname='local') as conn: 
            print conn.modules['os']
        
        with rpyc_conn_zerodeploy(hostname='sugon') as conn: 
            print conn.modules['os']
        

if __name__ == '__main__' : 
    pass
    #nose.run()

    if 0: #examine
        #suite = unittest.TestLoader().loadTestsFromTestCase(TestIt)
        #unittest.TextTestRunner(verbosity=0).run(suite)    
        unittest.main()
        
    else: 
        suite = unittest.TestSuite()
        add_list = [
           TestIt('test_temp'), 
           TestIt('test_rpyc_conn_local'), 
           TestIt('test_rpyc_conn_local_zero'), 
        ]
        for a in add_list: 
            suite.addTest(a)
        unittest.TextTestRunner().run(suite)
       

