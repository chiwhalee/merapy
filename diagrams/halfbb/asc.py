#!/usr/bin/python
from FourToOne import *
from template_halfbb import G_HBB
from merapy.diagrams.graph import *
import sys

def Simplify_Graph_HBB(g):
    ng = Simplify_Graph_I(g)
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W1(ng)    
    ng = Simplify_Graph_W2(ng)    
    ng = Simplify_Graph_V(ng)    
    ng = Simplify_Graph_O(ng)
    return ng

G=G_HBB
OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}

for i in xrange(7,12):
    oo = eval(oo2_tmpl % {"1":i, "2":i+1})
    gg = Simplify_Graph_HBB(Add_Graph(oo, G))
    check_graph(gg)
    weight= "   weight2(%d)=1.d0/4.d0"    
    if i==7:
        gg = Combine_OO(gg,OO)
        weight= "   weight2(%d)=1.d0/8.d0"
    if i==11:
        gg = Combine_OO(OO,gg)
        weight= "   weight2(%d)=1.d0/8.d0"
    check_graph(gg)
#    print_graph(gg)
    ii = i-7
    print weight % (ii)
    Output_Fortran(ii, "G2", "order2", gg)


