#!/usr/bin/python

from FourToOne import *
from HalfBB import G_HBB
from graph import *
import sys

from BinaryBinary import G_BB
from Ternary import G_3, G_3_top

def Simplify_Graph_BB(g):
    ng = Simplify_Graph_I(g)
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W(ng)    
    ng = Simplify_Graph_V(ng)    
    ng = Simplify_Graph_O(ng)
    return ng
def Simplify_Graph_HBB(g):
    ng = Simplify_Graph_I(g)
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W1(ng)    
    ng = Simplify_Graph_W2(ng)    
    ng = Simplify_Graph_V(ng)    
    ng = Simplify_Graph_O(ng)
    return ng

gg = G_HBB
i=7; j=12
# oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i, "2":i+1})
# oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":j, "2":j+1})
# gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
# Output_Fortran(i-7, "G_2", "order_2", gg)

i=7
oo = eval(oo2_tmpl % {"1":i, "2":i+1})
OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}

# sys.exit(0)
gg = Simplify_Graph_HBB(Add_Graph(oo, gg))
gg = Combine_OO(gg,OO)
#print_graph(gg)
#gg= Simplify_Graph_OO(gg)
# print_graph(gg)
# check_graph(gg)
Output_Fortran(i-7, "G_2", "order_2", gg)

sys.exit(0)
print_graph(gg, prefix="!!$ ")        
print
#        Output_Fortran(i-8, "G_2_2", "order_2_2", gg, jj=j-12)        

# for i in xrange(8,20):
#     gg = G_HBB
#     oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
#     gg = Simplify_Graph_HBB(Add_Graph(oo, gg))
#     print_graph(gg)
#     leg,comp,path,cost=get_path_cost(gg)#,pr=True)
#     print leg, comp, path, cost
#     print 
#     print    
    
# sys.exit(0)
for i in xrange(8,11):
    for j in xrange(12,15):
        gg = G_HBB
        oo1 = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo2 = eval(oo2_tmpl % {"1":j, "2":j+1})
        gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
        print_graph(gg)
        leg,comp,path,cost=get_path_cost(gg)#,pr=True)
        print leg, comp, path, cost
        print 
        print
sys.exit(0)

i=4
oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
#oo = eval(oo2_tmpl % {"1":i-1, "2":i, "3":i+1})
gg = Simplify_Graph(Add_Graph(oo, G_3))

print_graph(gg)
get_path_cost(gg, pr=True)
sys.exit(0)

for i in xrange(2,5):
    oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
    oo = eval(oo2_tmpl % {"1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G_3_top)
    gg=Simplify_Graph(gg)
    print_graph(gg)
    print
    continue

    ge,gm = GC2GE(gg)
    gm_nodes = {}
    
    j = 0    
    ks = gm.keys()
    for k in ks:
        j += 1
        name=k.split("_",1)[0]
        gm_nodes[k] = j
    leg,comp, path=get_path_cost(gg)
    order = []
    for k in path:
        order.append(gm_nodes[k])
    
    ii = i-1
    print "   G_3_top(%d)%%nNode=%d" % (ii,len(ks))
    print "   G_3_top(%d)%%Nodes=-1" % (ii)
    print "   G_3_top(%d)%%Edges=-1" % (ii)
    j = 0
    for k in ks:
        j=j+1
        name=k.split("_",1)[0]        
        print """   G_3_top(%d)%%Names(%d)="%s" """ % (ii,j, name)
        edges = ",".join([str(x) for x in gm[k]])
        leng = len(gm[k])
        print "   G_3_top(%d)%%Nodes(%d)=%d" % (ii,j,leng)
        print "   G_3_top(%d)%%Edges(1:%d,%d)=(/%s/)" % (ii,leng,j,edges)
    print "!$$ %d/%d" % (leg,comp)
    print "   order_3_top(1:%d, %d)=(/%s/)" % (len(ks), ii,
                                               ",".join([str(x) for x in order]))

sys.exit(0)

G=G_3_top

oo = eval(oo3_tmpl % {"1":3, "2":4, "3":5})
oo_1 = eval(oo2_tmpl % {"1":1, "2":2})
oo_2 = eval(oo2_tmpl % {"1":6, "2":7})

oo = eval(oo3_tmpl % {"1":0, "2":1, "3":2})
gg=Add_Graph(oo, G)
check_graph(gg)
gg=Simplify_Graph(gg)
print_graph(gg)
sys.exit(0)

for i in xrange(3,6):
    oo = eval(oo3_tmpl % {"1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G4)
    gg=Simplify_Graph(gg)
    ge,gm = GC2GE(gg)
    gm_nodes = {}
    
    j = 0    
    ks = gm.keys()
    for k in ks:
        j += 1
        name=k.split("_",1)[0]
        gm_nodes[k] = j
    leg,comp, path=get_path_cost(gg)
    order = []
    for k in path:
        order.append(gm_nodes[k])
    
    ii = i-4
    print "G3(%d)%%nNode=%d" % (ii,len(ks))
    j = 0
    for k in ks:
        j=j+1
        name=k.split("_",1)[0]        
        print """G3(%d)%%Names(%d)="%s" """ % (ii,j, name)
        edges = ",".join([str(x) for x in gm[k]])
        print "G3(%d)%%Edges(%d)=(/%s/)" % (ii,j,edges)
    print "!$$ %d/%d" % (leg,comp)
    print "order3(%d)=(/%s/)" % (ii, ",".join([str(x) for x in order]))
    print
    
for i in xrange(3,6):
    oo = eval(oo2_Ntmpl % {"name":"oo", "1":i-1, "2":i, "3":i+1})
    gg = Add_Graph(oo, G4)
    gg=Simplify_Graph(gg)
    ge,gm = GC2GE(gg)
    gm_nodes = {}
    
    j = 0    
    ks = gm.keys()
    for k in ks:
        j += 1
        name=k.split("_",1)[0]
        gm_nodes[k] = j
    leg,comp, path=get_path_cost(gg)
    order = []
    for k in path:
        order.append(gm_nodes[k])
    
    ii = i-4
    print "G2(%d)%%nNode=%d" % (ii,len(ks))
    j = 0
    for k in ks:
        j=j+1
        name=k.split("_",1)[0]        
        print """G2(%d)%%Names(%d)="%s" """ % (ii,j, name)
        edges = ",".join([str(x) for x in gm[k]])
        print "G2(%d)%%Edges(%d)=(/%s/)" % (ii,j,edges)
    print "!$$ %d/%d" % (leg,comp)
    print "order2(%d)=(/%s/)" % (ii, ",".join([str(x) for x in order]))
    print
    
for i in xrange(1,7):
    for l in xrange(i+2, 7):
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
        leg,comp, path=get_path_cost(gg)
        order = []
        for k in path:
            order.append(gm_nodes[k])
    
        ii = i-4; ll= l-4
        print "G_2_2(%d,%d)%%nNode=%d" % (ii,ll,len(ks))
        j = 0
        for k in ks:
            j=j+1
            name=k.split("_",1)[0]        
            print """G_2_2(%d,%d)%%Names(%d)="%s" """ % (ii,ll,j, name)
            edges = ",".join([str(x) for x in gm[k]])
            print "G_2_2(%d,%d)%%Edges(%d)=(/%s/)" % (ii,ll,j,edges)
        print "!$$ %d/%d" % (leg,comp)
        print "order_2_2(%d,%d)=(/%s/)" % (ii,ll,",".join([str(x) for x in order]))
        print
        
    
sys.exit(0)

ge,gm=GC2GE(gg)
for k in gm.keys():
    print k, gm[k]
g = gnode(gm)
max_leg = 1000; max_comp = 1000    

ks = g.keys()
ks.sort()
for k in ks:
    g[k].sort()    

for i in xrange(0,len(ks)):
    for j in xrange(i+1, len(ks)):
        start = ks[i]; end = ks[j]
        paths=find_all_path(g,start,end)
        for p in paths:
            leg,comp=contract_gmera(gm, p)#, True)
            if leg < 0:
                print 'error contracting'
                sys.exit(-1)
            if (comp == max_comp):
                if (leg < max_leg):
                    max_leg = leg
            elif (comp < max_comp):
                max_leg = leg
                max_comp = comp
            
            if (leg<=max_leg and comp<=max_comp):
                print leg, comp, p
                
#             if ((comp<=max_comp)):# or (leg<=max_leg)):
#                 max_comp = comp
#                 if (leg<=max_leg):
#                     max_leg = leg
#                 print leg,comp, p


sys.exit(0)

G0 = eval(G_tmpl % {"M3":-3, "M2":-2, "M1":-1, "0":0, "1":1, "2":2, "3":3, "4":4})
G8 = eval(G_tmpl % {"M3":5, "M2":6, "M1":7, "0":8, "1":9, "2":10, "3":11, "4":12})
# print_graph(G0)
# print
# print_graph(G8)
# print 

gg = Add_Graph(G0, G8)
oo_4_5 = eval(oo2_tmpl % {"1":4, "2":5})
oo_1_2 = eval(oo2_tmpl % {"1":1, "2":2})
#print oo3_tmpl % {"1":-3, "2":-2, "3":-1}
oo_M3_M2_M1 = eval(oo3_tmpl % {"1":-3, "2":-2, "3":-1})

gg = Add_Graph(oo_1_2, Add_Graph(oo_4_5, gg))
gg = Simplify_Graph(gg)
check_graph(gg)
print_graph(gg)
sys.exit(0)
g_M3_M2={
    "oo_M3_M2":[("I_M3", 1), ("I_M2", 1), ("I_M3",2), ("I_M2", 2)],
    "I_M3":[("oo_M3_M2",1), ("oo_M3_M2",3)],
    "I_M2":[("oo_M3_M2",2), ("oo_M3_M2",4)]
    }

g_M1_0 = {
    "oo_M1_0":[("I_M1", 1), ("I_0", 1), ("I_M1",2), ("I_0", 2)],
    "I_M1":[("oo_M1_0",1), ("oo_M1_0",3)],
    "I_0":[("oo_M1_0",2), ("oo_M1_0",4)]
    }


check_graph(g_M3_M2)
check_graph(g_M1_0)
#sys.exit()
gg = Simplify_Graph(G)#Add_Graph(g_M1_0,Add_Graph(g_M3_M2,G)))

g_5_6 = {
    "oo_5_6":[("I_5", 1), ("I_6", 1), ("I_5",2), ("I_6", 2)],
    "I_5":[("oo_5_6",1), ("oo_5_6",3)],
    "I_6":[("oo_5_6",2), ("oo_5_6",4)]
    }

gg = Add_Graph(g_5_6, gg)
#print_graph(gg)

gs = Break_Graph(gg)
for g in gs:
    print_graph(g)
    print 
# keys = gg.keys()
# keys.sort()
# for k in keys:
#     print k, gg[k]

check_graph(gg)
    
sys.exit(0)
