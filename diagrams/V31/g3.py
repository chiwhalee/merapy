#!/usr/bin/env python
#coding=utf8
from Ternary import G3, G3_new, oo2_tmpl, oo3_tmpl, ooN_tmpl
#from merapy.diagrams import Diagram, graphs_info
from merapy.diagrams.diagram import Diagram
import pprint

def reflect(i):
    return 4-i

def reflect_g(gs, g):
    return (reflect(g[1]), reflect(g[0]))


graphs_info.reflect_g = reflect_g

def Simplify_Graph3(g):
    """
        here amounts to calss inhereiance, reload the method simplify for ternary
    """
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify_Graph3
#G3 is bare data structure of a graph
#这样可以使得G具有加法等运算

G = Diagram(G3)
if 0:
    i = 4        
    oo = eval(oo2_tmpl % {"1":i, "2":i+1})
    oo = Diagram(oo)
    gg = G+oo
    gg = gg.Simplify()
    print gg.keys()


if 1:
    print '''subroutine init_mera_graph()
       use mGraphics
       use mMera
       implicit none
    '''

    print "   Tot_GN_GP_2=3"
    print "   Ref_GN_GP_2=2"
    for i in xrange(4,7):
        #generating L, C, R in sequal; 4 is at the centre leg of V tensor
        #replacing I_i, I_{i + 1} with oo
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()
        ii = i-3
        print "   weight2(%d)=1.d0/3.d0" % ii
        print "   GP_2(%d)= %d" % (ii, reflect(ii))
        ng.toFortran("G2", "order2", ii)    
        print 
        

    print "   Tot_GN_GP_3=3"
    print "   Ref_GN_GP_3=2"
    for i in xrange(5,8):
        #q: why not range(4, 7)? --- more symmetric?
        ooo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        ooo = Diagram(ooo)
        gg = G+ooo
        ng = gg.Simplify()
        ii = i-4
        print "   weight3(%d)=1.d0/3.d0" % ii
        print "   GP_3(%d)= %d" % (ii, reflect(ii))    
        ng.toFortran("G3", "order3", ii)
        print 


    GP_2_2 = graphs_info()
    for i in xrange(4,5):
        for j in xrange(6,7):
            oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
            oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
            oo1 = Diagram(oo1)
            oo2 = Diagram(oo2)
            gg = G+oo1
            gg = gg+oo2
            ng = gg.Simplify()
            ii = i-3
            jj = j-3

            gid = (ii,jj)
            GP_2_2.Add(gid)
            print "   weight_2_2(%d,%d)=1.d0/3.d0" % (ii,jj)
            ng.toFortran("G_2_2", "order_2_2", ii,jj)
            print 

    GP_2_2.Show("GP_2_2")

    GP_2_3 = graphs_info()    
    for i in xrange(4,7):
        for j in xrange(max(i+2,7), 10):
            oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
            oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
            oo1 = Diagram(oo1)
            oo2 = Diagram(oo2)
            gg = G+oo1
            gg = gg+oo2
            ng = gg.Simplify()
            ii = i-3; jj = j-6

            gid = (ii,jj)
            GP_2_3.Add(gid)
            print "   weight_2_3(%d,%d)=1.d0/3.d0" % (ii,jj)
            ng.toFortran("G_2_3", "order_2_3", ii,jj)
            print 


    GP_2_3.Show("GP_2_3")
    print 'end subroutine init_mera_graph'



#pprint.pprint(G)
#print G3


