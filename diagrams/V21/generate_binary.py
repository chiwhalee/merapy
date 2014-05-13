#!/usr/bin/env python
from pprint import pprint as pp
from merapy.diagrams.diagram import Diagram
from  template_binary  import G2, oo2_tmpl, oo3_tmpl
from merapy.diagrams.template_common import oo1_tmpl

from merapy.diagrams.graph import Combine_OO

def Simplify(ng, info=0):
    ng = ng.Simplify_I(info)
    ng = ng.Simplify_U(info)
    ng = ng.Simplify_VNodeG("V", "Vp", info)
    ng = ng.Simplify_O(info)
    ng.update()
    return ng

Diagram.Simplify = Simplify


def G_2_2_gen(G):
    """
    generating h2->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(1, 4):
    for i in [2]:
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        ng.toGraphics("G_2_2", "order2", i-2, weight=0.5, calc_order=True)    

def G_2_3_gen(G):
    """
    generating h2->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(1, 4):
    for i in [3]:
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        
        ng.toGraphics("G_2_2", "order2", i-3, weight=0.5, calc_order=True)    
    print "there is an issue in simplify_O. one needs to manualy modify 'O' to 'OO' in the generated graph"

def G_3_3_gen(G):
    """
        generating h2->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(1, 4):
    for i in [2, 3]:
        oo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        ng.toGraphics("G_3_3", "order_3_2", i-2, weight=0.5, calc_order=True)    



def generate_graph(G, inn, out, start, period):
    inn_map = {2:oo2_tmpl, 3:oo3_tmpl, 22:ooN_tmpl}
    #eval_pos = lambda i: {"1":i, "2":i+1, "3":i+2}
    end = start + period - 1
    o_tmpl = inn_map(inn)
    
    G = Diagram(G)
    G.check()
    for i in xrange(start, start + period):
        g_o = eval(o_tmpl % {"1":i, "2":i+1, "3":i+2})
        g_o = Diagram(g_o)

        gg = G  + g_o
        ng = gg.Simplify()#.connections
        if i==2:
            ng = ng.Combine_Node(OO, 'OO', 1, 1)
        if i==6:
            ng = O.Combine_Node(ng, 'O', 1, 1)
        
        ng.toGraphics("G_2_2", "order2", i, calc_order=False)    




if __name__ == "__main__":
    pass
    #G_2_2_gen(G2)
    #G_2_3_gen(G2)
    G_3_3_gen(G2)





