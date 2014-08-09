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
import time 
import datetime 

from collections import OrderedDict
from rpyc import Service
from rpyc.utils.server import ThreadedServer
import threading 
import rpyc
from plumbum import SshMachine

from merapy.context_util import rpyc_conn 
from merapy.utilities import load, save 

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
        #print(msg)
        self._conn._config.update(dict(
                    allow_all_attrs = True,
                    allow_pickle = True,
                    allow_getattr = True,
                    allow_setattr = True,
                    allow_delattr = True,
                    import_custom_exceptions = True,
                    instantiate_custom_exceptions = True,
                    instantiate_oldstyle_exceptions = True,
                ))        
        
    def exposed_register_job(self, backup_parpath, config=None): 
        """
            returns: 
                status: run, completed, failed
        """
        config = config if config is not None else {}
        pool = self.POOL
        job_id = hash(backup_parpath)
        if pool.has_key(job_id): 
            status = 'locked'
        else:
            status = 'open'
            self.POOL[job_id] = {'dir': backup_parpath,  }
            print  'job %s is added'%(str(job_id)[-6:], )
            
            pool['time_started'] = str(datetime.datetime.now()).split('.')[0]
            temp = ['hostname', 'pid']
            for i in temp: 
                self.POOL[i] = config.get(i)
        
        return status 
    
    def exposed_disregister_job(self, job_id): 
        self.POOL.pop(job_id)
        print 'job %s is disregisterd'%(str(job_id)[-6:], )
    
    def exposed_get_pool(self): 
        return self.POOL 
    
    def exposed_operate_pool(self, func): 
        print 'oooooperate'
        func(self.POOL)
    
    def exposed_test(self, num):     # 对于服务端来说， 只有以"exposed_"打头的方法才能被客户端调用，所以要提供给客户端的方法都得加"exposed_" 
        return 1+num
    
    def exposed_load(self): 
        pass
    
    def exposed_save(self): 
        pass
    
    def on_disconnect(self):
        msg = 'disconnected'
        #print(msg)

def get_pool(port=17012): 
    with rpyc_conn('local', 'service', port) as conn: 
        res = conn.root.get_pool() 
        res = rpyc.classic.obtain(res) 
        return res 
    
class TestIt(unittest.TestCase): 
    def setUp(self): 
        port = 17000
        server = ThreadedServer(TaskManager, port=port, auto_register=False)
        self.server = server 
        t = threading.Thread(target=self.server.start)
        t.start()
      
        ssh = SshMachine("localhost", user='zhli')
        #conn = rpyc.classic.ssh_connect(ssh, port)
        conn = rpyc.ssh_connect(ssh, 17000, config={"allow_pickle":True})
        self.conn = conn 
        self.ssh = ssh 
        
    def test_temp(self):
        from tempfile import mkdtemp
        for i in range(10): 
            self.conn.root.register_job(str(i)*10)
        print get_pool(17000)
        #print self.conn.root.get_pool()
    
    def tearDown(self): 
        self.conn.close()
        self.ssh.close()
        self.server.close()
        

if __name__ == '__main__' : 
    if 1: 
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
    
    if 0: 
        #sr = ThreadedServer(TestRpyc, port=9999, auto_register=False)
        s = ThreadedServer(TaskManager, port=17012, auto_register=False)
        s.start()
    






