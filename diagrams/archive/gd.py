#!/usr/bin/python
from FourToOne import *
from graph import *
import sys

D_0 = {}; D_8 = {}; D_M8={}
D_oo1 = {}; D_oo2 = {}; D_oo3 = {}
G_oo1 = {}; G_oo2 = {}; G_oo3 = {}

for i in xrange(-3,5):
    if (i<0): k="M"+str(abs(i))
    else: k = str(abs(i))
    D_0[k] = i
    D_8[k] = i+8
    D_M8[k] = i-8
    D_oo3["1"] = i; D_oo3["2"] = i+1; D_oo3["3"] = i+2
    G_oo1[i] = eval(oo1_tmpl % D_oo3)
    G_oo2[i] = eval(oo2_tmpl % D_oo3)
    G_oo3[i] = eval(oo3_tmpl % D_oo3)
    
G_0 = eval(G_tmpl % D_0)
G_8 = eval(G_tmpl % D_8)
G_M8 = eval(G_tmpl % D_M8)

G = Add_Graph(G_M8, Add_Graph(G_0, G_8))

for i in xrange(-3,5):
    for j in xrange(i+2,5): 
        og = Add_Graph(G_oo2[i], G_oo2[j])        
        gg = Simplify_Graph(Add_Graph(og, G))
        
        print_graph(gg)
        print 
        
# print G_0
# print G_8
# print G_M8

# print G_oo1[0]
# print G_oo2[0]
# print G_oo3[0]
