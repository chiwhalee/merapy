#coding=utf8

"""
define graphs for long range heisenberg model

"""

from merapy.diagrams.template_halfbb import G_HBB, oo1_tmpl, oo2_tmpl, oo2_Ntmpl, oo3_tmpl#G3, G3_new, oo2_tmpl, oo3_tmpl, ooN_tmpl
from merapy.diagrams.diagram import Diagram #, graphs_info
import pprint

if 0:
    def Simplify_Graph3(g):
        ng = g.Simplify_I()
        ng = ng.Simplify_U()
        ng = ng.Simplify_VNodeG("V", "Vp")
        ng = ng.Simplify_O()
        ng.update()
        return ng
    Diagram.Simplify = Simplify_Graph3

def Simplify(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_W1()
    ng = ng.Simplify_W2()
    ng = ng.Simplify_V()
    ng = ng.Simplify_O()
    return ng

Diagram.Simplify = Simplify

G = Diagram(G_HBB)
OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}

OO = Diagram(OO)


def gen_G_2_2():
    ng = {}
    for i in xrange(7, 12):
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        gg = G+oo
        gg = gg.Simplify()
        
        #print gg.keys()
        #gg["OO"] = gg["O"]
        #gg = gg + OO
        #gg.Combine_Node( OO, node, n1, n2)

        if 1:
            if i==7:
                gg = gg.Combine_Node(OO, 'OO', 1, 1)
                #weight= "   weight2(%d)=1.d0/8.d0"
            if i==11:
                gg = OO.Combine_Node(gg, 'OO', 1, 1)
                #weight= "   weight2(%d)=1.d0/8.d0"
        gid = i-7
        ng[gid] = gg
        #ng.show()
    for g in ng:
        ng[g].toGraphics("G_2_2", "order_2_2", g)

#gen_G_2_2()

def gen_G_3_2():
    ng = {}
    for i in xrange(7, 11):  #11):
        oo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        gg = G + oo
        gg = gg.Simplify()
        
        
        gid = i-7
        #weight="   weight3(%d)=1.d0/4.d0"
        #print weight % (gid)
        
        ng[gid] = gg
    
    for g in ng:
        ng[g].toGraphics("G_3_2", "order_3", g)

gen_G_3_2()



if 0:
#g2
    for i in xrange(7,12):
        oo = eval(oo2_tmpl % {"1":i, "2":i+1})
        gg = Simplify_Graph_HBB(Add_Graph(oo, G))
        check_graph(gg)
        weight= "   weight2(%d)=1.d0/4.d0"    
        if i==7:
            gg = Combine_OO(gg,OO)
            weight= "   weight2(%d)=1.d0/8.d0"
        if i==11:
            gg = Combine_OO(OO,gg)
            weight= "   weight2(%d)=1.d0/8.d0"
        check_graph(gg)
    #    print_graph(gg)
        ii = i-7
        print weight % (ii)
        Output_Fortran(ii, "G2", "order2", gg)
#g3
    for i in xrange(7,11):
        oo = eval(oo3_tmpl % {"1":i, "2":i+1, "3":i+2})
        gg = Simplify_Graph_HBB(Add_Graph(oo, G))
        weight="   weight3(%d)=1.d0/4.d0"
    #     if i==7 or i==11:
    #         weight="   weight3(%d)=1.d0/8.d0"
        ii = i-7
        print weight % (ii)
        Output_Fortran(ii, "G3", "order3", gg)

    for i in xrange(7,11):
        for j in xrange(i+2,12):
            gg = G_HBB
            oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i, "2":i+1})
            oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":j, "2":j+1})
            gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
            ii = i-7
            jj = j-7
            weight="   weight_2_2(%d,%d)=1.d0/4.d0"
            print weight % (ii,jj)
            Output_Fortran(ii, "G_2_2", "order_2_2", gg, jj=jj)
        
    for i in xrange(7,11):
        for j in xrange(12,16):
            gg = G_HBB
            oo1 = eval(oo2_Ntmpl % {"name":"oo1", "1":i, "2":i+1})
            oo2 = eval(oo2_Ntmpl % {"name":"oo2", "1":j, "2":j+1})
            gg = Simplify_Graph_HBB(Add_Graph(oo2,Add_Graph(oo1, gg)))
            ii = i-7
            jj = j-12
            weight="   weight_2_3(%d,%d)=1.d0/4.d0"
            print weight % (ii,jj)
            Output_Fortran(ii, "G_2_3", "order_2_3", gg, jj=jj)
    #        Output_Fortran(ii, "G3", "order3", gg)



if 0:
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






#pprint.pprint(G_2_2, indent=4)

if __name__ == "__main__" :
    if 0:
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

        







