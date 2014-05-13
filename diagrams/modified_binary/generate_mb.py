#!/usr/bin/env python
from merapy.diagrams.diagram import Diagram
from  template_mb import G_mb, oo2_tmpl, oo3_tmpl
from merapy.diagrams.template_common import oo1_tmpl

from merapy.diagrams.graph import Combine_OO
from pprint import pprint as pp
def Simplify(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify


#G = Diagram(G4)
#G.check()
if 1:
    OO = {"OO":[("I_8", 1), ("I_8", 2)],
          "I_8":[("OO", 1), ("OO", 2)]}
    OO = Diagram(OO)
    
    O = {"O":[("I_8", 1), ("I_8", 2)],
          "I_8":[("O", 1), ("O", 2)]}
    O = Diagram(O)



def G_2_2_gen(G):
    """
    generating h2->h2 ascending graph
    """
    G = Diagram(G)
    #pp(G.edges)
    #pp(G.edge_info)
    #exit()
    G.check()
    #ng = {}
    for i in xrange(0, 4):
    #for i in [1]:
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        #ng[a].name = a
        #print ng.connections
        
        ng.toGraphics("G_2_2", "order2", i, weight=1/4., calc_order=True)    
    print "there is an issue in simplify_O. one needs to manualy modify 'O' to 'OO' in the generated graph"

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
    from template_mb import *
    G_2_2_gen(G_mb)
    #G_3_2_gen(G4)
    #G_22_2_gen(G4_7block)
    #G_22_3_gen(G4_7block)
    #G_1_1_gen(G4)





