#coding=utf8

"""
define graphs for long range heisenberg model

"""

from merapy.diagrams.template_ternary import G3, G3_new, oo2_tmpl, oo3_tmpl, ooN_tmpl
from merapy.diagrams.Diagram import Diagram, graphs_info
import pprint

def Simplify_Graph3(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify_Graph3

G = Diagram(G3)


def G_2_2_gen():
    """
    generating h2->h2 ascending graph
    """
    ng = {}
    i = 4
    #for i in xrange(4,7):
    for a in ["L", "C", "R"]:   
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng[a] = gg.Simplify()#.connections
        i = i + 1
        #ng[a].name = a

    return ng

def G_3_3_gen():
    """
    generating h3->h3 ascending graph
    """
    ng = {}
    pos= ["R"]   
    #for i in xrange(5,8):
    
    i = 6
    for a in pos:
        ooo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        ooo = Diagram(ooo)
        gg = G+ooo
        ng[a] = gg.Simplify() #.connections
        #ii = i-4
        #ng[a].name = a
        i = i + 1
    return ng

def G_3_2_gen():
    """
    generating h3->h2 ascending graph
    2 graphs
    """
    ng = {}
    pos= ["C", "L"]   
    #for i in xrange(5,8):
    #for i in xrange(4,5):
    i = 4
    for a in pos:
        ooo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        ooo = Diagram(ooo)
        gg = G+ooo
        ng[a] = gg.Simplify()
        ii = i-4
        i = i + 1
    return ng

def G_22_2_gen():
    """
    o2*o2->o2 
    only 1
    """
    #for i in xrange(4,5):
    #    for j in xrange(6,7):
    #map for siis-> ss=h2 
    ng = {}
    for (i, j) in [(4, 6)]:
        oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
        oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
        oo1 = Diagram(oo1)
        oo2 = Diagram(oo2)
        gg = G+oo1
        gg = gg+oo2
        ng["C"] = gg.Simplify() #.connections
        #ng["C"].name = "C"
    return ng

def G_22_3_gen():
    """
    8 graphs in all
    (si, is)->H3
    """
    ng = {}
    for i in xrange(4,7):
        for j in xrange(max(i+2,7), 10):
            oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
            oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
            oo1 = Diagram(oo1)
            oo2 = Diagram(oo2)
            gg = G+oo1
            gg = gg+oo2
            
            ii = i-3; jj = j-6

            gid = (ii,jj)
            ng[gid] = gg.Simplify() #.connections
            #GP_2_3.Add(gid)
    return ng



if 0:
    def G_22_2_gen():
        pass

    def G_22_22_1_gen():
        pass
    def G_22_22_2_gen():
        pass
    def G_22_22_3_gen():
        pass
    def G_22_22_4_gen():
        pass



#pprint.pprint(G_2_2, indent=4)

if __name__ == "__main__" :
    pass
    G_2_2 = G_2_2_gen()
    G_3_3 = G_3_3_gen()
    G_3_2 = G_3_2_gen()
    
    G_22_2 = G_22_2_gen()
    G_22_3 = G_22_3_gen()
    
    #pprint.pprint(G_2_2["L"], indent=4)
    if 0:
        G_2_2["L"].show()
        print G_2_2["L"].edges
        print 
        print G_2_2["L"].nodes()
    for g in G_2_2:    
        pass
        #G_2_2["L"].toGraphics("G_2_2", "order2", "L")    
    for g in G_3_3:
        pass
        #G_3_3[g].toGraphics("G_3_3", "order_3_3", g)    

    for g in G_3_2:
        pass
        #G_3_2[g].toGraphics("G_3_2", "order_3_2", g)    

    for g in G_22_2:
        pass
        #G_22_2[g].toGraphics("G_22_2", "order_22_2", g)    

    for g in G_22_3:
        pass
        G_22_3[g].toGraphics("G_22_3", "order_22_3", g)    

    
    #print G_2_2["L"].show_edges()
    #pprint.pprint(G_2_2["L"], indent=2)
    #print G_2_2["L"]  #.connections
    
    #pprint.pprint(G_3_3)
    #pprint.pprint(G_22_2, indent=3)
    #pprint.pprint(G_22_3)








