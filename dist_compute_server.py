#!/usr/bin/env python        
# coding:utf-8

"""
    RPyC is a seamless library for integrating Python processes on many machines/processes
    to fully utilize computational resource, support communication between each hosts

"""
import sys 
import argparse 
import unittest 
import cPickle as pickle
import cloud
pickle_any = cloud.serialization.cloudpickle
import os
import time 
import datetime 
import pandas as pd 
from collections import OrderedDict 

from collections import OrderedDict
from rpyc import Service
from rpyc.utils.server import ThreadedServer
import threading 
import rpyc
from plumbum import SshMachine

from merapy.context_util import rpyc_conn , rpyc_conn_easy
from merapy.utilities import load, save 

path = os.path.abspath(__file__)
parpath = os.path.dirname(path)
DB_PATH = '/'.join([parpath, 'TASK_MANAGER.DB.PICKLE']) #'./TASK_MANAGER.DB.PICKLE' 

class TaskServer(Service):
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
        self.backup_path = path 
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
        
    def exposed_register_job(self, job_id, config=None, info=1): 
        """
            returns: 
                status: run, completed, failed
        """
        config = rpyc.classic.obtain(config) if config is not None else {}
        pool = self.POOL
        #job_id = hash(backup_parpath)
        if pool.has_key(job_id): 
            dir_status = 'locked'
        else:
            dir_status = 'open'
            
            job_info = OrderedDict()
            time_started = str(datetime.datetime.now()).split('.')[0]
            status = 'run'
            temp = [ 'status', 'time_started']
            dic = locals()
            for i in temp: 
                job_info[i] = dic[i]
            
            temp = ['backup_parpath', 'job_group_name', 'hostname', 'pid']
            for i in temp: 
                job_info[i] = config.get(i)
            pool[job_id] = job_info 
            if info>0: 
                print  'add job %s: %s'%(str(job_id)[-6:], job_info.items())
        return dir_status 
    
    def exposed_update_job(self, job_id, key, val): 
        self.POOL[job_id][key] = val 
    
    def exposed_disregister_job(self, job_id): 
        self.POOL.pop(job_id)
        print 'job %s is disregisterd'%(str(job_id)[-6:], )
    
    def exposed_get_pool(self): 
        return self.POOL 
    
    def exposed_set_pool(self, pool): 
        self.POOL = pool 
    
    def exposed_operate_pool(self, func): 
        print 'oooooperate'
        func(self.POOL)
    
    def exposed_test(self, num):     # 对于服务端来说， 只有以"exposed_"打头的方法才能被客户端调用，所以要提供给客户端的方法都得加"exposed_" 
        return 1+num
    
    def exposed_load(self): 
        pass
    
    def exposed_save(self): 
        pass

    def exposed_dump_pool(self): 
        msg = 'try tump pool to %s'%(self.backup_path[-30: ])
        save(self.POOL, self.backup_path)
        msg +=  '  ... done' 
        print msg 
    
    def exposed_recover(self): 
        msg = 'try recover pool from %s'%(self.backup_path[-30: ])
        temp = load(self.backup_path)
        self.__class__.POOL = temp 
        msg +=  '  ... done' 
        print self.POOL.keys() 
    
    def on_disconnect(self):
        msg = 'disconnected'
        #print(msg)
    

class TaskManager(object): 
    def __init__(self): 
        conn = rpyc_conn_easy('local', 'service', 17012)
        self.conn = conn 
        self.root = conn.root 
    
    def get_pool(self, as_df=1): 
        res = self.root.get_pool() 
        res = rpyc.classic.obtain(res) 
        if as_df: 
            res=  pd.DataFrame(res).T 
        return res 
    
    def clear_failed_jobs(self): 
        df = self.get_pool()
        #df = pd.DataFrame(pool).T 
        a=df[df['status']=='FAILED']
        for i in a.index:
            msg = 'clear %s'%(i, )
            print msg 
            self.root.disregister_job(i)
    
    def disregister_job(self, job_id_list): 
        for i in job_id_list: 
            msg = 'disregister %s'%(i, )
            self.root.disregister_job(i)
            msg += '  ... done' 
            print msg 
            
    def save(self): 
        pool = self.get_pool(as_df=0)
        save(pool, DB_PATH)
        print 'try to save pool to %s'%(DB_PATH), '   ... succeed'
    
    def close(self):
        self.conn.close()
        print  'conn closed'
        
    def __del__(self): 
        self.conn.close()
        

def get_pool(port=17012): 
    with rpyc_conn('local', 'service', port) as conn: 
        res = conn.root.get_pool() 
        res = rpyc.classic.obtain(res) 
        return res 

def clear_failed(): 
    pass 

class TestIt(unittest.TestCase): 
    def setUp(self): 
        port = 17000
        server = ThreadedServer(TaskServer, port=port, auto_register=False)
        self.server = server 
        t = threading.Thread(target=self.server.start)
        t.start()
      
        ssh = SshMachine("localhost", user='zhli')
        #conn = rpyc.classic.ssh_connect(ssh, port)
        conn = rpyc.ssh_connect(ssh, 17000, config={"allow_pickle":True})
        self.conn = conn 
        self.ssh = ssh 
        
    def test_temp(self):
        pass 
    
    def test_register_job(self):
        from tempfile import mkdtemp
        for i in range(5): 
            config = {'pid': i}
            self.conn.root.register_job(str(i), config=config)
        pool=get_pool(17000)
        df = pd.DataFrame(pool)
        print  df.T  
    
    def test_register_job_1(self):
        from vmps.run_ising.config import make_config 
        from vmps.main import Main 
        import numpy as np 
 
        hh = np.arange(0.5, 1.0, 0.1)

        dir_sur = ['a']
        alg_sur = ['', 'lam_guess', 'psi_guess']
        grid = Main.make_grid(hh, alg_sur[0:1], dir_sur)
        config_group = []

        for i in grid: 
            h, asur, dsur = i
            c = make_config(h, algorithm='idmrg', alg_surfix=asur, root1='', surfix = dsur)
            c.update(N0=4, NUM_OF_THREADS=1, energy_diff_tol=1e-10)
            c['which_increase_D'] = 'insert_sites'
            #c['schedule'] = [4, 8, 12]
            #c['schedule'] = range(2, 100, 2)
            c['schedule'] = [4, 10, 20, 40] 
            #c['backup_parpath'] = None
            c['num_iter_min'] = 20
            config_group.append(c)
            
        for c in config_group: 
            id = (hash(c['backup_parpath_local']), 'bbb', 'ccc')
            print  self.conn.root.register_job(id, c, info=0)

    def test_update_job(self): 
        #self.test_register_job()
        self.conn.root.register_job(1000,  )
        self.conn.root.update_job(1000, 'status', 'completed')
        pool=get_pool(17000)
        df = pd.DataFrame(pool)
        print  df.T  
        self.conn.root.update_job(1000, 'status', 'hello')
        pool=get_pool(17000)
        df = pd.DataFrame(pool)
        print  df.T  
    
    def tearDown(self): 
        self.conn.close()
        self.ssh.close()
        self.server.close()


if __name__ == '__main__' : 
    #if 1: 
    if len(sys.argv)>1: 
        parser = argparse.ArgumentParser(description='dont know')
        parser.add_argument('-s', '--start', action='store_true', default=False) 
        parser.add_argument('-p', '--port', type=int, default=17012)
        args = parser.parse_args()
        args = vars(args)
        port = args['port']
        if args['start']: 
            try: 
                s = ThreadedServer(TaskServer, port=port, auto_register=False)
                s.start()
                
            except KeyboardInterrupt: 
                msg = 'KeyboardInterrupt captured.'
                s.dump_pool()
                sys.exit( )
            except Exception as err: 
                print err 
        
    else: 
        #sr = ThreadedServer(TestRpyc, port=9999, auto_register=False)
        
        if 0:
            TestIt.test_temp=unittest.skip("skip test_temp")(TestIt.test_temp) 
            unittest.main()
            
        else: 
            suite = unittest.TestSuite()
            add_list = [
                #'test_temp', 
                #'test_register_job', 
                #'test_register_job_1', 
                'test_update_job', 
            ]
            for a in add_list: 
                suite.addTest(TestIt(a))

            unittest.TextTestRunner(verbosity=0).run(suite)
    






