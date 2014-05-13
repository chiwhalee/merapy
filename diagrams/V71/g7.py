#!/usr/bin/env python
from merapy.diagrams.diagram import Diagram
from Seven2One import G7, oo2_tmpl
import sys

def Simplify(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify

OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}
OO = Diagram(OO)

G = Diagram(G7)
G.check()
print '''subroutine init_mera_graph()
   use mGraphics
   use mMera
   implicit none
'''

for i in xrange(4,11):
    oo = eval(oo2_tmpl % {"1":i, "2":i+1})
    oo = Diagram(oo)
    gg = G+oo

    ii = i-3
    print ""
    print "   weight2(%d)=1.d0/7.d0" % ii
    ng = gg.Simplify()
    if ii==1 or ii==2:
        ng = ng.Combine_Node(OO, 'OO', 1, 1)
    if ii==7 or ii==6:
        ng = OO.Combine_Node(ng, 'OO', 1, 1)
        
    ng.toFortran("G2", "order2", ii)

print "end subroutine init_mera_graph"
    
