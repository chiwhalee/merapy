#!/usr/bin/env python        
# coding:utf-8

"""
    RPyC is a seamless library for integrating Python processes on many machines/processes
    to fully utilize computational resource, support communication between each hosts

"""
import unittest 
import cPickle as pickle
import cloud
pickle_any = cloud.serialization.cloudpickle
import os

from collections import OrderedDict
from rpyc import Service
from rpyc.utils.server import ThreadedServer
import rpyc
from plumbum import SshMachine

class TaskManager(Service):
    """
        POOL def: 
            primary key: job_id 
            field_names: 
                backup_parpath_local
                project
                param
    """
    POOL = OrderedDict()
    def on_connect(self):
        msg = "connected"
        path = './TASK_MANAGER.DB.PICKLE' 
        #if not os.path.exists(path): 
        #    with open(path, 'wb') as f:
        #        pass
        print(msg)
        
    def exposed_register_job(self, backup_parpath): 
        """
            returns: 
                status: run, completed, failed
        """
        pool = self.POOL
        job_id = hash(backup_parpath)
        if pool.has_key(job_id): 
            status= 'locked'
        else:
            status= 'open'
            self.POOL[job_id] = {'dir': backup_parpath,  }
            msg = 'job %s is added in POOL'%(job_id) 
            print self.POOL.keys()
            print msg 
        
        return status 
    
    def exposed_get_pool(self): 
        return self.POOL 
    
    def exposed_test(self, num):     # 对于服务端来说， 只有以"exposed_"打头的方法才能被客户端调用，所以要提供给客户端的方法都得加"exposed_" 
        return 1+num
    
    def exposed_load(self): 
        pass
    
    def exposed_save(self): 
        pass
    
    def on_disconnect(self):
        msg = 'disconnected'
        print(msg)
       
class TestIt(unittest.TestCase): 
    def setUp(self): 
        port = 0
        s = ThreadedServer(TaskManager, port=port, auto_register=False)
        self.s= s
        self.s.start()
       
        ssh = SshMachine("localhost", user='zhli')
        #conn = rpyc.classic.ssh_connect(ssh, port)
        conn = rpyc.ssh_connect(ssh, port)
        #cResult  = conn.root.test(11000)  # test是服务端的那个以"exposed_"开头的方法 
        self.conn = conn 
        self.ssh = ssh 
      
    def test_temp(self):
        from tempfile import mkdtemp
        for i in range(10): 
            pass 
    
    def tearDown(self): 
        print 'eeeeeeeeeeeee'
        self.conn.close()
        self.ssh.close()
        self.s.close()
        
        

if __name__ == '__main__' : 
    if 0: 
        if 0:
            TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
            unittest.main()
            
        else: 
            suite = unittest.TestSuite()
            add_list = [
                TestIt('test_temp'), 
            ]
            for a in add_list: 
                suite.addTest(a)

            unittest.TextTestRunner(verbosity=0).run(suite)
    
    #sr = ThreadedServer(TestRpyc, port=9999, auto_register=False)
    s = ThreadedServer(TaskManager, port=17012, auto_register=False)
    s.start()
    






