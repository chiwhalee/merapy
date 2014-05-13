#!/usr/bin/env python
#coding=utf8
from Ternary import G3, G3_new, oo2_tmpl, oo3_tmpl, ooN_tmpl
from Diagram import Diagram, graphs_info
import pprint

def reflect(i):
    return 4-i

def reflect_g(gs, g):
    return (reflect(g[1]), reflect(g[0]))


graphs_info.reflect_g = reflect_g

def Simplify_Graph3(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify_Graph3

G = Diagram(G3)

i = 4
oo = eval(oo2_tmpl % {"1":i, "2":i+1})
#oo = oo2_tmpl % {"1":i, "2":i+1}
print "ooooo \n", oo


def add(self, other):
    d = Diagram()
    k1s = set(self.keys())
    k2s = set(other.keys())
    keys_comm = k1s & k2s
    new_keys = k1s ^ k2s

    print "ccc", keys_comm #, "\n", k2s-k1s, "\n", k1s-k2s

    for k1 in (k2s-k1s):
        nv = []
        for v in other[k1]:
            #v is a tuple, node linked to k1
            print "v", v, 
            if (v[0] not in keys_comm):
                nv.append(v)
            else:
                nv.append(self[v[0]][v[1]-1])
        d[k1] = nv
        print "kkk A", k1, "nv", nv 
    for k1 in (k1s-k2s):
        nv = []
        for v in self[k1]:
            if (v[0] not in keys_comm):
                nv.append(v)
            else:
                nv.append(other[v[0]][v[1]-1])
        d[k1] = nv
        print "kkk B", k1, nv
    return d

#goo=add(G, oo)
goo = G + oo
print "before Simplify "
k0=set(goo.keys())
#goo=goo.Simplify_VNodeG("V", "Vp")

if 1:
    goo=goo.Simplify_I()
    print "afer Simplify I\t", 
    k1=set(goo.keys())
    print k0-k1

    goo=goo.Simplify_U()
    print "afer Simplify U\t", 
    k2=set(goo.keys())
    print k1-k2

    goo=goo.Simplify_V()
    print "afer Simplify V\t", 
    goo.update()
    k3=set(goo.keys())
    print k2-k3

    goo=goo.Simplify_O()
    print "afer Simplify O\t", 
    goo.update()
    k4=set(goo.keys())
    print k3-k4

#goo=goo.Simplify()
print "FINALLY ", 
kf = goo.keys()
print kf

goo.show()
goo.show_edges()

if 0:
    print "   Tot_GN_GP_2=3"
    print "   Ref_GN_GP_2=2"
    for i in xrange(4,7):

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
    #this merely generating one graph o2o2->o2
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
            #these are all 2*2 to 3 site mappings
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



#pprint.pprint(G)
#print G3


