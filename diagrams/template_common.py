#!/usr/bin/env python
#coding=utf8

"""
this file is a copy of Ternary.py

#连结结构不变，变的只是位置
#the keys are names of vertices,  the tuples in the values are infect edges labeled by (node, kth leg), 
note that the values are ordered.
I_n are just  place holders, they are to be replaced by hamiltonian ops. 
"""
if 1:
    U_tmpl = """ 
         "U_%(U)s":[("I_%(I1)s", 2), ("I_%(I2)s", 2), ("V_%(V1)s", 3), ("V_%(V2)s",1)], 
        "Up_%(U)s":[("Vp_%(V1)s", 4), ("Vp_%(V2)s",2), ("I_%(I1)s", 1), ("I_%(I2)s", 1)], 
        "I_%(I1)s":[("Up_%(U)s", 3), ("U_%(U)s", 1)],
        "I_%(I2)s":[("Up_%(U)s", 4), ("U_%(U)s", 2)], 
    """

    V_tmpl = """
         "V_%(V)s":[("U_%(U1)s", 4), ("I_%(I)s", 2), ("U_%(U2)s", 3), ("O_%(O)s", 2)],
        "Vp_%(V)s":[("O_%(O)s", 1),("Up_%(U1)s",2), ("I_%(I)s", 1), ("Up_%(U2)s",1)],
        "I_%(I)s":[("Vp_%(V)s", 3), ("V_%(V)s", 2)],
        "O_%(O)s":[("Vp_%(V)s", 1), ("V_%(V)s",4)], 
    """


ooN_tmpl = """{
    "%(name)s_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("%(name)s_%(1)s_%(2)s", 1), ("%(name)s_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("%(name)s_%(1)s_%(2)s", 2), ("%(name)s_%(1)s_%(2)s", 4)]
}    
"""

#one site op.
oo1_tmpl =  """{
    "oo_%(1)s":[("I_%(1)s",1), ("I_%(1)s", 2)],
    "I_%(1)s":[("oo_%(1)s", 1), ("oo_%(1)s", 2)]
    }
"""


#two site op.
oo2_tmpl =  """{
    "oo_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("oo_%(1)s_%(2)s", 1), ("oo_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("oo_%(1)s_%(2)s", 2), ("oo_%(1)s_%(2)s", 4)]
    }
"""

#3 site op.
oo3_tmpl =  """{
    "ooo_%(1)s_%(2)s_%(3)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(3)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2), ("I_%(3)s", 2)],
    "I_%(1)s":[("ooo_%(1)s_%(2)s_%(3)s", 1), ("ooo_%(1)s_%(2)s_%(3)s", 4)],
    "I_%(2)s":[("ooo_%(1)s_%(2)s_%(3)s", 2), ("ooo_%(1)s_%(2)s_%(3)s", 5)],
    "I_%(3)s":[("ooo_%(1)s_%(2)s_%(3)s", 3), ("ooo_%(1)s_%(2)s_%(3)s", 6)]
    }
"""

# where is this defined from?
G_3_top = {
    "V_0":[("U_1", 4), ("I_1", 2), ("U_0", 3), ("O_1", 2)],
    "V_1":[("U_0", 4), ("I_4", 2), ("U_1", 3), ("O_2", 2)],
    "U_0":[("I_2", 2), ("I_3", 2), ("V_0", 3), ("V_1", 1)],
    "U_1":[("I_5", 2), ("I_0", 2), ("V_1", 3), ("V_0", 1)],
    
    "Vp_0":[("O_1", 1), ("Up_1", 2), ("I_1", 1), ("Up_0", 1)],
    "Vp_1":[("O_2", 1), ("Up_0", 2), ("I_4", 1), ("Up_1", 1)],
    "Up_0":[("Vp_0", 4), ("Vp_1", 2), ("I_2", 1), ("I_3", 1)],
    "Up_1":[("Vp_1", 4), ("Vp_0", 2), ("I_5", 1), ("I_0", 1)],
    
    "I_0":[("Up_1", 4), ("U_1", 2)],
    "I_1":[("Vp_0", 3), ("V_0", 2)],
    "I_2":[("Up_0", 3), ("U_0", 1)],
    "I_3":[("Up_0", 4), ("U_0", 2)],
    "I_4":[("Vp_1", 3), ("V_1", 2)],
    "I_5":[("Up_1", 3), ("U_1", 1)],
    
    "O_1":[("Vp_0", 1), ("V_0", 4)],
    "O_2":[("Vp_1", 1), ("V_1", 4)]
    }



def make_template():
    k = 4
    for i in xrange(1,k):
        j=i*3
        #算子的相对位置
        D={"U":i, "V1":i-1, "V2":i, "I1":j-1, "I2":j}
        print U_tmpl % D
        
    for i in xrange(0,k):
        j=i*3+1
        D={"V":i, "I":j, "U1":i, "U2":i+1, "O":i}
        print V_tmpl % D


if __name__ == "__main__":
    #make_template()
    from merapy.diagrams.diagram import Diagram
    i = 1
    oo = oo2_tmpl % {"1":i, "2":i+1}
    oo = eval(oo)
    oo = Diagram(oo)
    oo.plot()
    
    



