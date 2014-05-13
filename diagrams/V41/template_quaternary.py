#!/usr/bin/env python
#coding=utf8
import pprint 
"""
V are of (4, 1) type tensors

"""



U_tmpl = """ "U_%(U)s":[("I_%(I1)s", 2), ("I_%(I2)s", 2), ("V_%(V1)s", 4), ("V_%(V2)s",1)], 
"Up_%(U)s":[("Vp_%(V1)s", 5), ("Vp_%(V2)s",2), ("I_%(I1)s", 1), ("I_%(I2)s", 1)], 
"I_%(I1)s":[("Up_%(U)s", 3), ("U_%(U)s", 1)],
"I_%(I2)s":[("Up_%(U)s", 4), ("U_%(U)s", 2)],
"""
V_tmpl = """ "V_%(V)s":[("U_%(U1)s", 4), ("I_%(I1)s", 2), ("I_%(I2)s", 2), ("U_%(U2)s", 3), ("O_%(O)s", 2)],
"Vp_%(V)s":[("O_%(O)s", 1), ("Up_%(U1)s", 2), ("I_%(I1)s", 1), ("I_%(I2)s", 1), ("Up_%(U2)s", 1)],
"I_%(I1)s":[("Vp_%(V)s", 3), ("V_%(V)s", 2)],
"I_%(I2)s":[("Vp_%(V)s", 4), ("V_%(V)s", 3)],
"O_%(O)s":[("Vp_%(V)s",1), ("V_%(V)s", 5)],
"""

oo2_tmpl =  """{
    "oo_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("oo_%(1)s_%(2)s", 1), ("oo_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("oo_%(1)s_%(2)s", 2), ("oo_%(1)s_%(2)s", 4)]
    }
"""

oo3_tmpl =  """{
    "ooo_%(1)s_%(2)s_%(3)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(3)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2), ("I_%(3)s", 2)],
    "I_%(1)s":[("ooo_%(1)s_%(2)s_%(3)s", 1), ("ooo_%(1)s_%(2)s_%(3)s", 4)],
    "I_%(2)s":[("ooo_%(1)s_%(2)s_%(3)s", 2), ("ooo_%(1)s_%(2)s_%(3)s", 5)],
    "I_%(3)s":[("ooo_%(1)s_%(2)s_%(3)s", 3), ("ooo_%(1)s_%(2)s_%(3)s", 6)]
    }
"""

ooN_tmpl = """{
    "%(name)s_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("%(name)s_%(1)s_%(2)s", 1), ("%(name)s_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("%(name)s_%(1)s_%(2)s", 2), ("%(name)s_%(1)s_%(2)s", 4)]
}    
"""

#G4 with 2 blocks, made by make_template(2), with little modifications
G4 = {
     "U_2":[("I_4", 2), ("I_5", 2), ("V_1", 4), ("V_2",1)], 
    "Up_2":[("Vp_1", 5), ("Vp_2",2), ("I_4", 1), ("I_5", 1)], 
    "I_4":[("Up_2", 3), ("U_2", 1)],
    "I_5":[("Up_2", 4), ("U_2", 2)],

     "V_1":[("I_1", 2), ("I_2", 2), ("I_3", 2), ("U_2", 3), ("O_1", 2)],
    "Vp_1":[("O_1", 1), ("I_1", 1), ("I_2", 1), ("I_3", 1), ("Up_2", 1)],
    "I_1":[("Vp_1", 2), ("V_1", 1)],
    "I_2":[("Vp_1", 3), ("V_1", 2)],
    "I_3":[("Vp_1", 4), ("V_1", 3)],
    "O_1":[("Vp_1",1), ("V_1", 5)],

     "V_2":[("U_2", 4), ("I_6", 2), ("I_7", 2), ("I_8", 2), ("O_2", 2)],
    "Vp_2":[("O_2", 1), ("Up_2", 2), ("I_6", 1), ("I_7", 1), ("I_8", 1)],
    "I_6":[("Vp_2", 3), ("V_2", 2)],
    "I_7":[("Vp_2", 4), ("V_2", 3)],
    "I_8":[("Vp_2", 5), ("V_2", 4)],
    "O_2":[("Vp_2",1), ("V_2", 5)]
}


#G4 with 7 blocks, made by make_template(7), with little modifications
G4_7block= {
     "U_2":[("I_4", 2), ("I_5", 2), ("V_1", 4), ("V_2",1)], 
    "Up_2":[("Vp_1", 5), ("Vp_2",2), ("I_4", 1), ("I_5", 1)], 
    "I_4":[("Up_2", 3), ("U_2", 1)],
    "I_5":[("Up_2", 4), ("U_2", 2)],

     "U_3":[("I_8", 2), ("I_9", 2), ("V_2", 4), ("V_3",1)], 
    "Up_3":[("Vp_2", 5), ("Vp_3",2), ("I_8", 1), ("I_9", 1)], 
    "I_8":[("Up_3", 3), ("U_3", 1)],
    "I_9":[("Up_3", 4), ("U_3", 2)],

     "U_4":[("I_12", 2), ("I_13", 2), ("V_3", 4), ("V_4",1)], 
    "Up_4":[("Vp_3", 5), ("Vp_4",2), ("I_12", 1), ("I_13", 1)], 
    "I_12":[("Up_4", 3), ("U_4", 1)],
    "I_13":[("Up_4", 4), ("U_4", 2)],

     "U_5":[("I_16", 2), ("I_17", 2), ("V_4", 4), ("V_5",1)], 
    "Up_5":[("Vp_4", 5), ("Vp_5",2), ("I_16", 1), ("I_17", 1)], 
    "I_16":[("Up_5", 3), ("U_5", 1)],
    "I_17":[("Up_5", 4), ("U_5", 2)],

     "U_6":[("I_20", 2), ("I_21", 2), ("V_5", 4), ("V_6",1)], 
    "Up_6":[("Vp_5", 5), ("Vp_6",2), ("I_20", 1), ("I_21", 1)], 
    "I_20":[("Up_6", 3), ("U_6", 1)],
    "I_21":[("Up_6", 4), ("U_6", 2)],

     "U_7":[("I_24", 2), ("I_25", 2), ("V_6", 4), ("V_7",1)], 
    "Up_7":[("Vp_6", 5), ("Vp_7",2), ("I_24", 1), ("I_25", 1)], 
    "I_24":[("Up_7", 3), ("U_7", 1)],
    "I_25":[("Up_7", 4), ("U_7", 2)],

     #"V_1":[("U_1", 4), ("I_2", 2), ("I_3", 2), ("U_2", 3), ("O_1", 2)],
     "V_1":[("I_1", 2), ("I_2", 2), ("I_3", 2), ("U_2", 3), ("O_1", 2)],
    #"Vp_1":[("O_1", 1), ("Up_1", 2), ("I_2", 1), ("I_3", 1), ("Up_2", 1)],
    "Vp_1":[("O_1", 1), ("I_1", 1), ("I_2", 1), ("I_3", 1), ("Up_2", 1)],
    #add I_1
    "I_1":[("Vp_1", 2), ("V_1", 1)],

    "I_2":[("Vp_1", 3), ("V_1", 2)],
    "I_3":[("Vp_1", 4), ("V_1", 3)],
    "O_1":[("Vp_1",1), ("V_1", 5)],

     "V_2":[("U_2", 4), ("I_6", 2), ("I_7", 2), ("U_3", 3), ("O_2", 2)],
    "Vp_2":[("O_2", 1), ("Up_2", 2), ("I_6", 1), ("I_7", 1), ("Up_3", 1)],
    "I_6":[("Vp_2", 3), ("V_2", 2)],
    "I_7":[("Vp_2", 4), ("V_2", 3)],
    "O_2":[("Vp_2",1), ("V_2", 5)],

     "V_3":[("U_3", 4), ("I_10", 2), ("I_11", 2), ("U_4", 3), ("O_3", 2)],
    "Vp_3":[("O_3", 1), ("Up_3", 2), ("I_10", 1), ("I_11", 1), ("Up_4", 1)],
    "I_10":[("Vp_3", 3), ("V_3", 2)],
    "I_11":[("Vp_3", 4), ("V_3", 3)],
    "O_3":[("Vp_3",1), ("V_3", 5)],

     "V_4":[("U_4", 4), ("I_14", 2), ("I_15", 2), ("U_5", 3), ("O_4", 2)],
    "Vp_4":[("O_4", 1), ("Up_4", 2), ("I_14", 1), ("I_15", 1), ("Up_5", 1)],
    "I_14":[("Vp_4", 3), ("V_4", 2)],
    "I_15":[("Vp_4", 4), ("V_4", 3)],
    "O_4":[("Vp_4",1), ("V_4", 5)],

     "V_5":[("U_5", 4), ("I_18", 2), ("I_19", 2), ("U_6", 3), ("O_5", 2)],
    "Vp_5":[("O_5", 1), ("Up_5", 2), ("I_18", 1), ("I_19", 1), ("Up_6", 1)],
    "I_18":[("Vp_5", 3), ("V_5", 2)],
    "I_19":[("Vp_5", 4), ("V_5", 3)],
    "O_5":[("Vp_5",1), ("V_5", 5)],

     "V_6":[("U_6", 4), ("I_22", 2), ("I_23", 2), ("U_7", 3), ("O_6", 2)],
    "Vp_6":[("O_6", 1), ("Up_6", 2), ("I_22", 1), ("I_23", 1), ("Up_7", 1)],
    "I_22":[("Vp_6", 3), ("V_6", 2)],
    "I_23":[("Vp_6", 4), ("V_6", 3)],
    "O_6":[("Vp_6",1), ("V_6", 5)],

     #"V_7":[("U_7", 4), ("I_26", 2), ("I_27", 2), ("U_8", 3), ("O_7", 2)],
     "V_7":[("U_7", 4), ("I_26", 2), ("I_27", 2), ("I_28", 2), ("O_7", 2)],
    #"Vp_7":[("O_7", 1), ("Up_7", 2), ("I_26", 1), ("I_27", 1), ("Up_8", 1)],
    "Vp_7":[("O_7", 1), ("Up_7", 2), ("I_26", 1), ("I_27", 1), ("I_28", 1)],
    #add I_28
    "I_28":[("Vp_7", 5), ("V_7", 4)],

    "I_26":[("Vp_7", 3), ("V_7", 2)],
    "I_27":[("Vp_7", 4), ("V_7", 3)],
    "O_7":[("Vp_7",1), ("V_7", 5)],
}
    
def make_template(num_of_V):
    """
        joint parts of U, V to G4raw.
        note that this func wont give exactly G4raw, rather one need to deal with boundary legs by hand!
    """
    num_of_U = num_of_V-1
    res= "{"
    
    for i in range(2, 2 + num_of_U ):
        res+= "\n" 
        j = 4*(i-1)+1
        D={"U":i, "V1":i-1, "V2":i, "I1":j-1, "I2":j}
        res +=  U_tmpl % D
    
    for i in range(1, 1 + num_of_V):
        res+= "\n" 
        j = 4*(i-1)+1
        D={"V":i, "U1":i, "U2":i+1, "O":i, "I1":j+1, "I2":j+2,"I3":j+3}
        res +=  V_tmpl % D


    res+= "}" 
    
    #res = eval(res)
    return res
   
if __name__ == "__main__":
    from merapy.diagrams.diagram import compare_template
    
    res = make_template(2)
    res= eval(res)
    
    compare_template(G4, res)
    


