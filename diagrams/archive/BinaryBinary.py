#!/usr/bin/env python

U_tmpl = """ "U_%(U)s":[("I_%(I1)s", 2), ("I_%(I2)s", 2), ("W_%(W1)s", 2), ("W_%(W2)s",1)], 
"Up_%(U)s":[("Wp_%(W1)s", 3), ("Wp_%(W2)s",2), ("I_%(I1)s", 1), ("I_%(I2)s", 1)], 
"I_%(I1)s":[("Up_%(U)s", 3), ("U_%(U)s", 1)],
"I_%(I2)s":[("Up_%(U)s", 4), ("U_%(U)s", 2)],
"""
W_tmpl = """ "W_%(W)s":[("U_%(U1)s", 4), ("U_%(U2)s", 3), ("V_%(V)s", %(leg1)s)], 
"Wp_%(W)s":[("Vp_%(V)s",%(leg2)s), ("Up_%(U1)s",2), ("Up_%(U2)s",1)],
"""
W1_tmpl = """" "W_%(W)s":[("I_%(I)s", 2), ("U_%(U)s", 3), ("V_%(V)s", 2)], 
"Wp_%(W)s":[("Vp_%(V)s",3), ("I_%(I)s",1), ("Up_%(U)s",1)],
"I_%(I)s":[("Wp_%(W)s", 2), ("W_%(W)s",1)],
"""
W2_tmpl = """" "W_%(W)s":[("U_%(U)s", 4), ("I_%(I)s", 2), ("V_%(V)s", 1)], 
"Wp_%(W)s":[("Vp_%(V)s",2), ("Up_%(U)s", 2), ("I_%(I)s",1)],
"I_%(I)s":[("Wp_%(W)s", 3), ("W_%(W)s",2)],
"""
V_tmpl = """ "V_%(V)s":[("W_%(W1)s", 3), ("W_%(W2)s", 3), ("O_%(O)s", 2)],
"Vp_%(V)s":[("O_%(O)s", 1), ("Wp_%(W1)s", 1), ("Wp_%(W2)s", 1)],
"O_%(O)s":[("Vp_%(V)s",1), ("V_%(V)s", 3)],
"""

G_BB={"U_0":[("I_0", 2), ("I_1", 2), ("W_0", 2), ("W_1",1)], 
"Up_0":[("Wp_0", 3), ("Wp_1",2), ("I_0", 1), ("I_1", 1)], 
"I_0":[("Up_0", 3), ("U_0", 1)],
"I_1":[("Up_0", 4), ("U_0", 2)],

 "U_1":[("I_2", 2), ("I_3", 2), ("W_1", 2), ("W_2",1)], 
"Up_1":[("Wp_1", 3), ("Wp_2",2), ("I_2", 1), ("I_3", 1)], 
"I_2":[("Up_1", 3), ("U_1", 1)],
"I_3":[("Up_1", 4), ("U_1", 2)],

 "U_2":[("I_4", 2), ("I_5", 2), ("W_2", 2), ("W_3",1)], 
"Up_2":[("Wp_2", 3), ("Wp_3",2), ("I_4", 1), ("I_5", 1)], 
"I_4":[("Up_2", 3), ("U_2", 1)],
"I_5":[("Up_2", 4), ("U_2", 2)],

 "U_3":[("I_6", 2), ("I_7", 2), ("W_3", 2), ("W_4",1)], 
"Up_3":[("Wp_3", 3), ("Wp_4",2), ("I_6", 1), ("I_7", 1)], 
"I_6":[("Up_3", 3), ("U_3", 1)],
"I_7":[("Up_3", 4), ("U_3", 2)],

 "U_4":[("I_8", 2), ("I_9", 2), ("W_4", 2), ("W_5",1)], 
"Up_4":[("Wp_4", 3), ("Wp_5",2), ("I_8", 1), ("I_9", 1)], 
"I_8":[("Up_4", 3), ("U_4", 1)],
"I_9":[("Up_4", 4), ("U_4", 2)],

 "U_5":[("I_10", 2), ("I_11", 2), ("W_5", 2), ("W_6",1)], 
"Up_5":[("Wp_5", 3), ("Wp_6",2), ("I_10", 1), ("I_11", 1)], 
"I_10":[("Up_5", 3), ("U_5", 1)],
"I_11":[("Up_5", 4), ("U_5", 2)],

 "U_6":[("I_12", 2), ("I_13", 2), ("W_6", 2), ("W_7",1)], 
"Up_6":[("Wp_6", 3), ("Wp_7",2), ("I_12", 1), ("I_13", 1)], 
"I_12":[("Up_6", 3), ("U_6", 1)],
"I_13":[("Up_6", 4), ("U_6", 2)],

 "U_7":[("I_14", 2), ("I_15", 2), ("W_7", 2), ("W_8",1)], 
"Up_7":[("Wp_7", 3), ("Wp_8",2), ("I_14", 1), ("I_15", 1)], 
"I_14":[("Up_7", 3), ("U_7", 1)],
"I_15":[("Up_7", 4), ("U_7", 2)],

 "U_8":[("I_16", 2), ("I_17", 2), ("W_8", 2), ("W_9",1)], 
"Up_8":[("Wp_8", 3), ("Wp_9",2), ("I_16", 1), ("I_17", 1)], 
"I_16":[("Up_8", 3), ("U_8", 1)],
"I_17":[("Up_8", 4), ("U_8", 2)],

 "W_0":[("I_M1", 2), ("U_0", 3), ("V_0", 1)], 
"Wp_0":[("Vp_0", 2), ("I_M1",1), ("Up_0",1)],
"I_M1":[("Wp_0", 2), ("W_0", 1)],
   
 "W_1":[("U_0", 4), ("U_1", 3), ("V_0", 2)], 
"Wp_1":[("Vp_0",3), ("Up_0",2), ("Up_1",1)],

 "W_2":[("U_1", 4), ("U_2", 3), ("V_1", 1)], 
"Wp_2":[("Vp_1",2), ("Up_1",2), ("Up_2",1)],

 "W_3":[("U_2", 4), ("U_3", 3), ("V_1", 2)], 
"Wp_3":[("Vp_1",3), ("Up_2",2), ("Up_3",1)],

 "W_4":[("U_3", 4), ("U_4", 3), ("V_2", 1)], 
"Wp_4":[("Vp_2",2), ("Up_3",2), ("Up_4",1)],

 "W_5":[("U_4", 4), ("U_5", 3), ("V_2", 2)], 
"Wp_5":[("Vp_2",3), ("Up_4",2), ("Up_5",1)],

 "W_6":[("U_5", 4), ("U_6", 3), ("V_3", 1)], 
"Wp_6":[("Vp_3",2), ("Up_5",2), ("Up_6",1)],

 "W_7":[("U_6", 4), ("U_7", 3), ("V_3", 2)], 
"Wp_7":[("Vp_3",3), ("Up_6",2), ("Up_7",1)],

 "W_8":[("U_7", 4), ("U_8", 3), ("V_4", 1)], 
"Wp_8":[("Vp_4",2), ("Up_7",2), ("Up_8",1)],

 "W_9":[("U_8", 4), ("I_20", 2), ("V_4", 2)], 
"Wp_9":[("Vp_4",3), ("Up_8",2), ("I_20",1)],
"I_20":[("Wp_9",3), ("W_9", 2)],

 "V_0":[("W_0", 3), ("W_1", 3), ("O_0", 2)],
"Vp_0":[("O_0", 1), ("Wp_0", 1), ("Wp_1", 1)],
"O_0":[("Vp_0",1), ("V_0", 3)],

 "V_1":[("W_2", 3), ("W_3", 3), ("O_1", 2)],
"Vp_1":[("O_1", 1), ("Wp_2", 1), ("Wp_3", 1)],
"O_1":[("Vp_1",1), ("V_1", 3)],

 "V_2":[("W_4", 3), ("W_5", 3), ("O_2", 2)],
"Vp_2":[("O_2", 1), ("Wp_4", 1), ("Wp_5", 1)],
"O_2":[("Vp_2",1), ("V_2", 3)],

 "V_3":[("W_6", 3), ("W_7", 3), ("O_3", 2)],
"Vp_3":[("O_3", 1), ("Wp_6", 1), ("Wp_7", 1)],
"O_3":[("Vp_3",1), ("V_3", 3)],

 "V_4":[("W_8", 3), ("W_9", 3), ("O_4", 2)],
"Vp_4":[("O_4", 1), ("Wp_8", 1), ("Wp_9", 1)],
"O_4":[("Vp_4",1), ("V_4", 3)]
}

def main():
    for i in xrange(0, 20, 2): 
        D={"U":i/2, "I1":i, "I2":i+1, "W1":i/2, "W2":i/2+1}
        print U_tmpl%D

    print
    for i in xrange(0, 10): 
        l,j=divmod(i,2)
        D={"W":i, "U1":i-1, "U2":i, "V":i/2, "leg1":j+1,"leg2":j+2}    
        print W_tmpl%D
    
    print
    for i in xrange(0, 10, 2): 
        l,j=divmod(i,2)
        D={"V":i/2, "W1":i, "W2":i+1, "O":i/2}
        print V_tmpl%D

def Four2One():
    for i in xrange(0, 5): 
        D={"U":i, "I1":4*i+1, "I2":4*i+2, "W1":2*i, "W2":2*i+1}
        print U_tmpl%D

    print
    for i in xrange(0, 10, 2): 
        D={"W":i, "U":i/2, "I":2*i, "V":i/2}
        print W1_tmpl%D
#        l,j=divmod(i,2)
#         D={"W":i, "U1":i-1, "U2":i, "V":i/2, "leg1":j+1,"leg2":j+2}    
#         print W_tmpl%D
    
#    print
#     for i in xrange(0, 10, 2): 
#         l,j=divmod(i,2)
#         D={"V":i/2, "W1":i, "W2":i+1, "O":i/2}
#         print V_tmpl%D
    
if __name__ == "__main__":
    Four2One()
#    main()
