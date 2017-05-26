#coding=utf8
"""
"""
import unittest 
import inspect 
import os 
import numpy as np
import collections
from tempfile import mkdtemp 
import cPickle as pickle 
import warnings 
import sys
import string
import random
import gzip 
import zlib 
try: 
    import cloud 
    pickle_any = cloud.serialization.cloudpickle
except ImportError as err: 
    warnings.warn(str(err))

#__all__=["getch"]

getch=sys.stdin.read

def repr_util(obj, keys):
    """
    deprecated
    repr obj inteligently
    """
    #dic = boj.__dict__
    #keys= dic.keys()
    
    for k in keys:
        temp = repr(self.__dict__[k])
        if k == 'data':
            temp = str(self.__dict__[k][:self.totDim].round(10)) 
        elif k == 'Dims':
            temp  = repr(self.__dict__[k][:self.rank]) 
        elif k == "Block_idx":
            temp = str(self.Block_idx[0, :self.nidx]) + "\n\t" + str(self.Block_idx[1, :self.nidx]) + "\n\t" + str(self.Block_idx[2, :self.nidx]) + "\n" 
            #temp = str(self.Block_idx[:, 0:self.nidx])
        elif k == "Addr_idx":
            temp = str(self.__dict__[k][:self.rank, :self.nidx]) 

        str1 += k+":\n\t" + temp +"\n" 

    res=str0+ str1 +str2
    return res
       
def repr_format(k, value_str):
    return  k+":\n\t" + value_str +"\n"  

x =r"""
    G3[1].size=8
    G3[1].nodes=-1
    G3[1].edges=-1
    G3[1].names[1]="Up" #Up_2
    G3[1].nodes[1]=4
    G3[1].edges[0:4,1]=[1,2,3,4]
    G3[1].names[2]="OO" #OO
    G3[1].nodes[2]=4
    G3[1].edges[0:4,2]=[5,6,7,8]
    G3[1].names[3]="ooo" #oo_4_5_6
    G3[1].nodes[3]=6
    G3[1].edges[0:6,3]=[9,3,4,10,11,12]
    G3[1].names[4]="U" #U_2
    G3[1].nodes[4]=4
    G3[1].edges[0:4,4]=[11,12,13,14]
    G3[1].names[5]="Vp" #Vp_1
    G3[1].nodes[5]=4
    G3[1].edges[0:4,5]=[5,15,9,1]
    G3[1].names[6]="Vp" #Vp_2
    G3[1].nodes[6]=4
    G3[1].edges[0:4,6]=[6,2,16,17]
    G3[1].names[7]="V" #V_1
    G3[1].nodes[7]=4
    G3[1].edges[0:4,7]=[15,10,13,7]
    G3[1].names[8]="V" #V_2
    G3[1].nodes[8]=4
    G3[1].edges[0:4,8]=[14,16,17,8]
#$ 6/8         2.01020001E+16
    order3[0:8, 1]=[4,3,1,5,7,2,6,8]
    
# Diagram for toplevel
    GTop.size = 4        
    GTop.nodes = -1
    GTop.edges = -1
    GTop.names[1] = "oo"
    GTop.nodes[1] = 4
    GTop.edges[0:4, 1] = [3,4,1,2]
    
    GTop.names[2] = "V"
    GTop.nodes[2] = 3
    GTop.edges[0:3, 2] = [1,2,7]
    
    GTop.names[3] = "Vp"
    GTop.nodes[3] = 3
    GTop.edges[0:3, 3] = [8,3,4]
    
    GTop.names[4] = "OO"
    GTop.nodes[4] = 2
    GTop.edges[0:2, 4] = [8,7]        
    order_top[0:4] = [1,2,4,3]    
"""


def get_local(l1, l2):
    res= ""
    for i in l2:
        #print i, l1[i], 
        if i in l1.keys():
            res+= repr(i) +" " +  str(l1[i])  + " " 
    return res

def replace_num(a1):
    import re
    """ regular expression """
    
    def repl1(matchobj):
        a = matchobj.group(1)
        b = eval(a)-1
        c = str(b)
        return  c +','
    patt1="(\d*),"
    #print re.findall(patt1, a1)
    a2=re.sub(patt1,repl1,a1)


    patt2 = "(\d*)\]"
    def repl2(matchobj):
        a = matchobj.group(1)
        b = eval(a)-1
        c = str(b)
        return  c +']'
    #print re.findall(patt2, a2)
    a3=re.sub(patt2,repl2,a2)
    #print a3

    patt3 = ":(\d*)"
    def repl3(matchobj):
        a = matchobj.group(1)
        b = eval(a) + 1
        c = str(b)
        return  ":" + c 
    #print re.findall(patt3, a3)
    a4=re.sub(patt3,repl3,a3)
    #print a4

    patt4 = "\[(-?\d*)\]\."
    def repl4(matchobj):
        a = matchobj.group(1)
        b = eval(a) + 1
        c = str(b)
        return  "[" + c  + "]."
    #print re.findall(patt4, a4)
    a5=re.sub(patt4,repl4,a4)
    print a5

    patt5 = "(\d)\]([=\n])"
    def repl5(matchobj):
        a = matchobj.group(1)
        a2 = matchobj.group(2)
        b = eval(a) - 1
        c = str(b)
        return   c + ']' + a2 
    #print re.findall(patt4, a4)
    #a6=re.sub(patt5,repl5,a5)
    #print a6
            
def random_str(size=6, head='', tail=''):
    chars=string.ascii_uppercase + string.digits 
    middle = ''.join(random.choice(chars) for _ in range(size))
    res=''.join([head, middle, tail])
    return res

def load(path, as_str=False, info=0):
    with open(path, 'rb') as f:
        s= f.read()
        head = s[:10]
        #print 'oooo', hex(ord(head[0])),  hex(ord(head[1]))
        #print 'hhhh', repr(head)
        #assert len(head)>2 
        if not len(head)>2: 
            raise EOFError('file is proberbly broken, %s'%(path, ))
        #head1 is different compress levels 
        if hex(ord(head[0]))=='0x78' and hex(ord(head[1])) in ['0x1', '0x9c', '0x5e', '0xda']: 
            msg = 'discompress file ... '
            s= zlib.decompress(s)
            msg += 'done' 
            if info>0: 
                print msg
        if as_str: 
            res= s
        else: 
            res= pickle.loads(s)
        #res= pickle_any.loads(s)
    return res

def save(obj, path=None, compress=False, compress_level=2, as_str=False, info=0): 
    """
        params:
            path: can be None if as_str=True
    """
    if not compress: 
        if not as_str: 
            out = open(path, "wb")
            #pickle.dump(obj, out)
            pickle_any.dump(obj, out)
            out.close()
        else: 
            return pickle_any.dumps(obj)
    else: 
        try: 
            s= pickle_any.dumps(obj, pickle.HIGHEST_PROTOCOL)
        except Exception as err: 
            print 'pickling error, diagonstic which value cant be dumped: '
            if isinstance(obj, dict): 
                temp = obj 
            else: 
                temp = dir(obj)
                temp = {t: getattr(obj, t) for t in temp if t[:2]!= '__'}
            print sorted(temp.keys())
            for k, v in sorted(temp.items()): 
                try: 
                    print k,   
                    pickle_any.dumps(v)
                    print '--pass'
                except Exception as e: 
                    print '--fail'
                    #print e
                    print str(e)[: 50]
                
            raise err 
        msg = 'compressing file ...'
        z = zlib.compress(s, compress_level)
        msg += 'done' 
        if info>0: 
            print msg 
        if not as_str: 
            with open(path, 'wb') as f:  
                #pickle.dump(obj, f)
                #s= pickle.dumps(obj, pickle.HIGHEST_PROTOCOL)
                f.write(z)
        else: 
            return z 

def print_vars(dic, var_name_list=None, head=None, sep=', ', key_val_sep='=', 
        show_header=0, return_str=False, fault_tol=True, round=5):
    """
        
    """
    if head is None: 
        lineno = inspect.currentframe().f_back.f_lineno
        frm = inspect.stack()[1]
        module_name = inspect.getmodule(frm[0]).__name__
        #inspect.getfile(inspect.currentframe())
        caller_name =  frm[3] #inspect.stack()[1][3] 
        #head = '\n---print at %s line %s---\n'%(module_name, lineno)
        head = '\n[%5s LINE%d]  '%(module_name, lineno, )
    dic['np'] = np  #inject np to the namespace so that it can be evaled 
    if var_name_list is None: 
        var_name_list = sorted(dic)
    if show_header: 
        msg = ''.join(['-'*10, str(var_name_list), '-'*10])
        print msg 
    def get_val(x): 
        #if not '.' in str(x): 
        if not isinstance(x, str): 
            res = dic[x]
        else: 
            #a, b = x.split('.')
            #res = dic[a].__getattribute__(b)
            try: 
                res= eval(x, dic)
            except NameError as err: 
                res= 'not defined'
            #except AttributeError as err: 
            except Exception as err: 
                if fault_tol:
                    res='error: ' +  str(err)
                else:
                    raise  
                
        #if isinstance(res, str):   # only print a str 
        #    return  str(res)
        #else: 
        #    return  str(x) + key_val_sep + str(res)
        if isinstance(res, np.ndarray) and res.dtype==float:
            res= res.round(round)
        res_str = str(res)
        if isinstance(res, np.ndarray) and res.ndim>1:
            res_str = '\n\t' + res_str.replace('\n', '\n\t')
        return  str(x) + key_val_sep + res_str
            
    
    res = head
    res  +=  sep.join([ get_val(i)  for i in  var_name_list])
    if not return_str:         
        print  res
    else: 
        return res 


def dict_to_object(dic): 
    class OBJECT:
        def __init__(self, **entries): 
            self.__dict__.update(entries)
        def __repr__(self): 
            #return '<%s>' % str('\n '.join('%s : %s' % (k, repr(v)) for (k, v) in self.__dict__.iteritems()))             
            return ', '.join(self.__dict__.keys())
    return OBJECT(**dic)


def send_email(msg): 
    """
    issue: not secure to save password in a file
        find solution refer to 
        https://www.google.com/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8#q=python%20email%20encrypt%20password
        http://stackoverflow.com/questions/157938/hiding-a-password-in-a-python-script
    """
    import smtplib
    import email.mime.text
    # my test mail
    mail_username = 'zhihuali@mail.ustc.edu.cn'
    s= ''
    #issue: this is risky to put it in a file, fix this!
    with open('./p.txt', 'r') as f:
        s=f.readline()[:-1]
    mail_password = s
    from_addr = mail_username
    to_addrs=('chiwhalee@gmail.com')

    # HOST & PORT
    #HOST = 'smtp.gmail.com'
    HOST = 'smtp.ustc.edu.cn'
    PORT = 25

    # Create SMTP Object
    smtp = smtplib.SMTP()
    print 'connecting ...'

    # show the debug log
    smtp.set_debuglevel(1)

    # connet
    try:
        print smtp.connect(HOST,PORT)
    except:
        print 'CONNECT ERROR ****'
    # gmail uses ssl
    #smtp.starttls()
    # login with username & password
    try:
        print 'loginning ...'
        smtp.login(mail_username,mail_password)
    except:
        print 'LOGIN ERROR ****'
    
    # fill content with MIMEText's object 
    the_email = email.mime.text.MIMEText('send from merapy')
    the_email['From'] = from_addr
    the_email['To'] = ';'.join(to_addrs)
    #the_email['Subject']='hello , today is a special day'
    the_email['Subject'] = msg.get('subject', 'no subject')
    print the_email.as_string()
    
    smtp.sendmail(from_addr,to_addrs,the_email.as_string())
    smtp.quit()


class OrderedSet(collections.MutableSet):
    """
        taken from 
            http://stackoverflow.com/questions/1653970/does-python-have-an-ordered-set
            http://code.activestate.com/recipes/576694/
    """
    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

class TestIt(unittest.TestCase): 
    def test_temp(self): 
        s = OrderedSet('sssssdfsdfdsf')
        print_vars(vars(),  ['s', 'type(s)', 'list(s)'])
    
    def test_save_load(self): 
        for i in [0, 1]: 
            z = {'sda': 3, 'sesdf': 1000000}
            path = '/tmp/'  +  random_str( )
            print  'path is ', path 
            save(z,  path, i)
            self.assertTrue(os.path.exists(path))
            b = load(path)
            print 'bbb', b 
            self.assertEqual(z, b)
            os.remove(path)
    
    def xtest_send_email(self):
        msg = {'a':1}
        send_email(msg)
           

if __name__ == "__main__": 
    

    if 0: 
        unittest.main()
    else: 
        suite = unittest.TestSuite()
        
        add_list = [
        #'test_temp', 
        'xtest_send_email', 
        #'test_save_load', 
          
        ]
        for a in add_list: 
            suite.addTest(TestIt(a))

        unittest.TextTestRunner().run(suite)





