#!/usr/bin/python

from FourToOne import *
from graph import *
import sys

#from BinaryBinary import G_BB
from Ternary import G_3, G_3_top

G=G_3

oo = eval(oo3_tmpl % {"1":3, "2":4, "3":5})
oo_1 = eval(oo2_tmpl % {"1":2, "2":3})
oo_2 = eval(oo2_tmpl % {"1":6, "2":7})
G4 = G_3

print '''subroutine init_mera_graph()
   use mGraphics
   use mMERA
   implicit none
!$$    type(Graphics):: G3(-1:1), G2(-1:1), G_2_2(-3:3, -3:3) 
!$$    integer order3(MaxNode,-1:1), order2(MaxNode,-1:1), order_2_2(MaxNode,-3:3, -3:3)
'''

for i in xrange(3,6):
    oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G4)
    gg=Simplify_Graph(gg)    
    Output_Fortran(i-4, "G3", "order3", gg)
    
for i in xrange(3,6):
    oo = eval(oo2_Ntmpl % {"name":"oo", "1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G4)
    gg=Simplify_Graph(gg)
    Output_Fortran(i-4, "G2", "order2", gg)
    
for i in xrange(2,7):
    for l in xrange(i+2, 8):
        oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i-1, "2":i, "3":i+1})
        oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":l-1, "2":l, "3":l+1})
        
        gg = Add_Graph(oo1, Add_Graph(oo2, G4))
        gg=Simplify_Graph(gg)
        ge,gm = GC2GE(gg)
        gm_nodes = {}
    
        j = 0    
        ks = gm.keys()
        for k in ks:
            j += 1
            name=k.split("_",1)[0]
            gm_nodes[k] = j
        leg,comp, path, cost=get_path_cost(gg)
        order = []
        for k in path:
            order.append(gm_nodes[k])
    
        print_graph(gg, "!!$ ")
        ii = i-5; ll= l-5
        print "   G_2_2(%d,%d)%%nNode=%d" % (ii,ll,len(ks))
        print "   G_2_2(%d,%d)%%Nodes=-1" % (ii,ll)
        print "   G_2_2(%d,%d)%%Edges=-1" % (ii,ll)
        
        j = 0
        for k in ks:
            j=j+1
            name=k.split("_",1)[0]        
            edges = ",".join([str(x) for x in gm[k]])
            leng = len(gm[k])
            if (name == "OO" and leng == 6): name = "OOO"
            if (name == "oo" and leng == 6): name = "ooo"
            print """   G_2_2(%d,%d)%%Names(%d)="%s" """ % (ii,ll,j, name)
            print "   G_2_2(%d,%d)%%Nodes(%d)=%d" % (ii,ll,j,leng)
            print "   G_2_2(%d,%d)%%Edges(1:%d,%d)=(/%s/)" % (ii,ll,leng,j,edges)
#            print "G_2_2(%d,%d)%%Edges(%d)=(/%s/)" % (ii,ll,j,edges)
        print "!$$ %d/%d  %G" % (leg,comp, cost)
        print "   order_2_2(1:%d, %d, %d)=(/%s/)" % (len(ks), ii, ll,
                                           ",".join([str(x) for x in order]))
        print


for i in xrange(1,4):
    oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G_3_top)
    gg=Simplify_Graph(gg)
    Output_Fortran(i-2, "G_3_top", "order_3_top", gg)
    
for i in xrange(2,5):
    oo = eval(oo2_tmpl % {"1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G_3_top)
    gg=Simplify_Graph(gg)
    Output_Fortran(i-3, "G_2_top", "order_2_top", gg)

for i in xrange(1,3):
    if (i==1):
        oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":1, "2":2})
        oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":3, "2":4})    
    else:
        oo1 = eval(oo2_Ntmpl % {"name":"oo2", "1":1, "2":2})
        oo2 = eval(oo2_Ntmpl % {"name":"oo1", "1":3, "2":4})    
        
    gg = Add_Graph(oo2, Add_Graph(oo1, G_3_top))
    gg=Simplify_Graph(gg)    
    Output_Fortran(i, "G_2_2_top", "order_2_2_top", gg)
        
print "end subroutine init_mera_graph"
