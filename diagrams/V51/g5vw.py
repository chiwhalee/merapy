#!/usr/bin/env python
from Five2OneVW import G5, oo2_tmpl, oo3_tmpl, ooN_tmpl
from Diagram import Diagram, graphs_info
import sys

def reflect(i):
    return 6-i

def reflect_g(gs, g):
    return (reflect(g[1]), reflect(g[0]))

graphs_info.reflect_g = reflect_g

def Simplify(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("W1", "W1p")
    ng = ng.Simplify_VNodeG("W2", "W2p")
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify

G=Diagram(G5)
OO = {"O":[("I_8", 1), ("I_8", 2)],
      "I_8":[("O", 1), ("O", 2)]}
OO = Diagram(OO)

print '''subroutine init_mera_graph()
   use mGraphics
   use mMera
   implicit none
'''

for i in xrange(8, 13):
    oo = eval(oo2_tmpl % {"1":i, "2":i+1})
    oo = Diagram(oo)
    gg = G+oo
    ng = gg.Simplify()
    ii = i-7

    if ng.has_key('O'):
        ng = ng.Combine_Node(OO, 'O', 1, 1)
        ng = ng.Node_Rename('O', 'OO')
        ng.show()
        print 
#    print "   weight_2(%d)=1.d0/5.d0" % (ii)
#    ng.toFortran("G_2", "order_2", ii)
#    print 

sys.exit(0)

for i in xrange(8, 13):
    oo = eval(oo3_tmpl % {"1":i, "2":i+1,"3":i+3})
    oo = Diagram(oo)
    gg = G+oo
    ng = gg.Simplify()
    ii = i-7
    if ng.has_key('O'):
        ng.show()
    print "   weight_3(%d)=1.d0/5.d0" % (ii)
    ng.toFortran("G_3", "order_2", ii)
    print 
    
GP_2_2 = graphs_info()
for i in xrange(8, 13):
    for j in xrange(i+2, 13):
        oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
        oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
        oo1 = Diagram(oo1)
        oo2 = Diagram(oo2)

        ii = i-7; jj = j-7; gid=(ii,jj)
        GP_2_2.Add(gid)
        
        gg = G+oo1
        gg = gg+oo2
        ng = gg.Simplify()
        print "   weight_2_2(%d,%d)=1.d0/5.d0" % (ii,jj)
        ng.toFortran("G_2_2", "order_2_2", ii,jj)
        print 
    
GP_2_2.Show("GP_2_2")

GP_2_3 = graphs_info()
for i in xrange(8, 13):
    for j in xrange(max(i+2,13), 18):
        oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
        oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
        oo1 = Diagram(oo1)
        oo2 = Diagram(oo2)

        ii = i-7; jj = j-12; gid=(ii,jj)
        GP_2_3.Add(gid)
        
        gg = G+oo1
        gg = gg+oo2
        ng = gg.Simplify()
        print "   weight_2_3(%d,%d)=1.d0/5.d0" % (ii,jj)
        ng.toFortran("G_2_3", "order_2_3", ii,jj)
        print 
    
GP_2_3.Show("GP_2_3")

print 'end subroutine init_mera_graph'
