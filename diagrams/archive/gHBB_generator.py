#!/usr/bin/env python

from FourToOne import *
from merapy.diagrams.graph import *
from HalfBB import G_HBB
import sys

def Simplify_Graph_HBB(g):
    ng = Simplify_Graph_I(g)
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W1(ng)    
    ng = Simplify_Graph_W2(ng)    
    ng = Simplify_Graph_V(ng)    
    ng = Simplify_Graph_O(ng)
    return ng

print '''subroutine init_mera_graph()
   use mGraphics
   use mMERA
   implicit none
'''

G=G_HBB
OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}

for i in xrange(7,12):
    oo = eval(oo2_tmpl % {"1":i, "2":i+1})
    print Add_Graph(oo, G).has_key("OO")
    gg = Simplify_Graph_HBB(Add_Graph(oo, G))
    check_graph(gg)
    weight= "   weight2(%d)=1.d0/4.d0"    
    print gg.keys()
    exit()
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

for i in xrange(7,11):
    oo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
    gg = Simplify_Graph_HBB(Add_Graph(oo, G))
    weight="   weight3(%d)=1.d0/4.d0"
#     if i==7 or i==11:
#         weight="   weight3(%d)=1.d0/8.d0"
    ii = i-7
    print weight % (ii)
    Output_Fortran(ii, "G3", "order3", gg)

for i in xrange(7,11):
    for j in xrange(i+2,12):
        gg = G_HBB
        oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i, "2":i+1})
        oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":j, "2":j+1})
        gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
        ii = i-7
        jj = j-7
        weight="   weight_2_2(%d,%d)=1.d0/4.d0"
        print weight % (ii,jj)
        Output_Fortran(ii, "G_2_2", "order_2_2", gg, jj=jj)
    
for i in xrange(7,11):
    for j in xrange(12,16):
        gg = G_HBB
        oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i, "2":i+1})
        oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":j, "2":j+1})
        gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
        ii = i-7
        jj = j-12
        weight="   weight_2_3(%d,%d)=1.d0/4.d0"
        print weight % (ii,jj)
        Output_Fortran(ii, "G_2_3", "order_2_3", gg, jj=jj)
#        Output_Fortran(ii, "G3", "order3", gg)

print "end subroutine init_mera_graph"
        
