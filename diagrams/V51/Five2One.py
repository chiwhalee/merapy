#!/usr/bin/env python
#from Diagram import Diagram

U_tmpl = """ "U_%(U)s":[("I_%(I1)s", 2), ("I_%(I2)s", 2), ("V_%(V1)s", 5), ("V_%(V2)s",1)], 
"Up_%(U)s":[("Vp_%(V1)s", 6), ("Vp_%(V2)s",2), ("I_%(I1)s", 1), ("I_%(I2)s", 1)], 
"I_%(I1)s":[("Up_%(U)s", 3), ("U_%(U)s", 1)],
"I_%(I2)s":[("Up_%(U)s", 4), ("U_%(U)s", 2)],
"""
V_tmpl = """ "V_%(V)s":[("U_%(U1)s", 4), ("I_%(I1)s", 2), ("I_%(I2)s", 2), ("I_%(I3)s", 2), ("U_%(U2)s", 3), ("O_%(O)s", 2)],
"Vp_%(V)s":[("O_%(O)s", 1), ("Up_%(U1)s", 2), ("I_%(I1)s", 1), ("I_%(I2)s", 1), ("I_%(I3)s", 1), ("Up_%(U2)s", 1)],
"I_%(I1)s":[("Vp_%(V)s", 3), ("V_%(V)s", 2)],
"I_%(I2)s":[("Vp_%(V)s", 4), ("V_%(V)s", 3)],
"I_%(I3)s":[("Vp_%(V)s", 5), ("V_%(V)s", 4)],
"O_%(O)s":[("Vp_%(V)s",1), ("V_%(V)s", 6)],
"""

G5 = {
"U_2":[("I_5", 2), ("I_6", 2), ("V_1", 5), ("V_2",1)], 
"Up_2":[("Vp_1", 6), ("Vp_2",2), ("I_5", 1), ("I_6", 1)], 
"I_5":[("Up_2", 3), ("U_2", 1)],
"I_6":[("Up_2", 4), ("U_2", 2)],

 "V_1":[("I_1", 2), ("I_2", 2), ("I_3", 2), ("I_4", 2), ("U_2", 3), ("O_1", 2)],
"Vp_1":[("O_1", 1), ("I_1", 1), ("I_2", 1), ("I_3", 1), ("I_4", 1), ("Up_2", 1)],
"I_1":[("Vp_1", 2), ("V_1", 1)],
"I_2":[("Vp_1", 3), ("V_1", 2)],
"I_3":[("Vp_1", 4), ("V_1", 3)],
"I_4":[("Vp_1", 5), ("V_1", 4)],
"O_1":[("Vp_1",1), ("V_1", 6)],

 "V_2":[("U_2", 4), ("I_7", 2), ("I_8", 2), ("I_9", 2), ("I_10", 2), ("O_2", 2)],
"Vp_2":[("O_2", 1), ("Up_2", 2), ("I_7", 1), ("I_8", 1), ("I_9", 1), ("I_10", 1)],
"I_7":[("Vp_2", 3), ("V_2", 2)],
"I_8":[("Vp_2", 4), ("V_2", 3)],
"I_9":[("Vp_2", 5), ("V_2", 4)],
"I_10":[("Vp_2", 6), ("V_2", 5)],
"O_2":[("Vp_2",1), ("V_2", 6)],    
}

oo2_tmpl =  """{
    "oo_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("oo_%(1)s_%(2)s", 1), ("oo_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("oo_%(1)s_%(2)s", 2), ("oo_%(1)s_%(2)s", 4)]
    }
"""

def main():
    i = 2
    j = 5*(i-1)+1
    D={"U":i, "V1":i-1, "V2":i, "I1":j-1, "I2":j}
    print U_tmpl % D
    
    i = 1
    j = 5*(i-1)+1
    D={"V":i, "U1":i, "U2":i+1, "O":i, "I1":j+1, "I2":j+2,"I3":j+3}
    print V_tmpl % D

    i = 2
    j = 5*(i-1)+1
    D={"V":i, "U1":i, "U2":i+1, "O":i, "I1":j+1, "I2":j+2,"I3":j+3}
    print V_tmpl % D
    
    
if __name__ == "__main__":
    main()
        
