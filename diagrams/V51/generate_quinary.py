#!/usr/bin/env python
from merapy.diagrams.diagram import Diagram
from  template_quinary import G5, G5_2block
from merapy.diagrams.template_common import oo1_tmpl, oo2_tmpl, oo3_tmpl

from merapy.diagrams.graph import Combine_OO

def Simplify(g):
    ng = g.Simplify_I()
    #exit()
    ng = ng.Simplify_U()
    ng = ng.Simplify_VNodeG("V", "Vp")
    ng = ng.Simplify_O()
    ng.update()
    return ng

Diagram.Simplify = Simplify


#G = Diagram(G4)
#G.check()
if 1:
    OO = {"OO":[("I", 1), ("I", 2)],
          "I":[("OO", 1), ("OO", 2)]}
    OO = Diagram(OO)
    
    O = {"O":[("I", 1), ("I", 2)],
          "I":[("O", 1), ("O", 2)]}
    O = Diagram(O)

def G_1_1_gen(G):
    """
        this is only used to calc scaling dim
        so just use 1 out of 4 equivlent pos of o1

    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(2, 7):
    for i in [2]:
        oo = eval(oo1_tmpl % {"1":i})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        #ng[a].name = a
        if 1:
            pass
            #if i==2:
            #    ng = ng.Combine_Node(O, 'O', 1, 1)
            #if i==6:
            #    ng = O.Combine_Node(ng, 'O', 1, 1)
        else:
            pass
            #if i==2:
            #    ng = Combine_OO(ng, OO)
            #if i==6:
            #    ng = Combine_OO(OO, ng)
        
        ng.toGraphics("G_1_1", "order2", i, weight=1.0)    


def G_2_2_gen(G):
    """
    generating h2->h2 ascending graph
    """
    G = Diagram(G)
    G.check()

    #ng = {}
    for i in xrange(2, 7):
    #for i in [2]:
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()#.connections
        #ng[a].name = a
        if 1:     
            if i in [2, 3]:
                ng = ng.Combine_Node(O, 'O', 1, 1)
            if i==3 and False:
                ng = O.Combine_Node(ng, 'O', 1, 1)
        if 0:
            if i==2:
                ng = Combine_OO(ng, OO)
            if i==6:
                ng = Combine_OO(OO, ng)

        name = 'G_2_2'
        #ng.plot(path=name + "-%d"%i + '.png', remove_O=False, layout='neato', len=2.5) 
        ng.toGraphics(name, "order2", i-2, weight=1./5)    
    print "there is an issue in simplify_O. one needs to manualy modify 'O' to 'OO' in the generated graph"

def G_3_2_gen(G):
    """
        generating h3->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    for i in xrange(2, 6):
        ooo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        ooo = Diagram(ooo)

        gg = G  + ooo
        ng = gg.Simplify()#.connections
        #if i==2:
        #    ng = ng.Combine_Node(OO, 'OO', 1, 1)
        #if i==6:
        #    ng = O.Combine_Node(ng, 'O', 1, 1)
        
        #ng.plot(path=str(i) + '.png', remove_O=True) 
        ng.toGraphics("G_3_2", "order_3_2", i, weight=0.25)    

def G_22_2_gen(G):
    """
        generating h3->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(2, 6):
    #    for j in xrange(2, 6):
    for i, j in [(2, 4), (3, 5), (4, 6), 
                (2, 5), (3, 6), 
                (2, 6), 
                (2, 10)]:
        oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
        oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
        oo1 = Diagram(oo1)
        oo2 = Diagram(oo2)
        gg = G+oo1
        gg = gg+oo2
        
        ng = gg.Simplify()#.connections
        if i%2 == 4 or j%2 == 4 :
            ng = ng.Combine_Node(OO, 'OO', 1, 1)
        if 0:
            if i==2:
                ng = ng.Combine_Node(OO, 'OO', 1, 1)
            if i==6:
                ng = O.Combine_Node(ng, 'O', 1, 1)
        
        #ng.plot(path=str(i) + '.png', remove_O=True) 
        ng.toGraphics("G_22_2", "order_22_2", (i-2, j-2+1), weight=1./4)    

def G_22_3_gen(G):
    """
        generating h3->h2 ascending graph
    """
    G = Diagram(G)
    G.check()
    #ng = {}
    #for i in xrange(2, 6):
    #    for j in xrange(2, 6):
    #for i, j in [(2, 4), (3, 5), (4, 6), 
    #            (2, 5), (3, 6), 
    #            (2, 6)]:
    
    for ii in [0, 1, 2, 3]:
        for jj in range(max(ii + 3, 6), 10):
            
            if ii == 0 and jj == 9:  #it is special that 0, 9 would yeild 'OO', so move it to G_22_2_gen
                continue

            i = ii + 2
            j = jj + 2 - 1
            
            oo1 = eval(ooN_tmpl % {"name":"oo1", "1":i, "2":i+1, "3":i+2})
            oo2 = eval(ooN_tmpl % {"name":"oo2", "1":j, "2":j+1, "3":j+2})
            oo1 = Diagram(oo1)
            oo2 = Diagram(oo2)
            gg = G+oo1
            gg = gg+oo2
            
            ng = gg.Simplify()#.connections
            #if i%2 == 4 or j%2 == 4 :
            #    ng = ng.Combine_Node(OO, 'OO', 1, 1)
            if 0:
                if i==2:
                    ng = ng.Combine_Node(OO, 'OO', 1, 1)
                if i==6:
                    ng = O.Combine_Node(ng, 'O', 1, 1)
            
            ng.toGraphics("G_22_3", "order_22_3", (ii, jj), weight=1./4, calc_order=True)    

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

def G_2_2_gen_old():
    for i in xrange(2,7):
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        ng = gg.Simplify()
        #the following deals with boundary iregular legs/nodes
        if i==2:
            ng = ng.Combine_Node(OO, 'O', 1, 1)
        if i==6:
            ng = OO.Combine_Node(ng, 'O', 1, 1)
            
        print ""
        ng.toFortran("G2", "order2", i)
    if 0:        
        print ng.connections
        ng.check()
        max_leg, max_comp, optimal_path, max_cost = ng.find_optimal_path()
        print max_leg, max_comp, max_cost
        print optimal_path


if __name__ == "__main__":
    pass
    G_2_2_gen(G5_2block)
    #G_3_2_gen(G4)
    #G_22_2_gen(G4_7block)
    #G_22_3_gen(G4_7block)
    #G_1_1_gen(G4)





