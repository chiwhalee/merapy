#coding=utf8
"""
"""
import numpy as np
import sys
import string
import random

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


#replace_graph(x)
if __name__ == "__main__": 
    replace_num(x)




