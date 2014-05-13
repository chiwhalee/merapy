#!/usr/bin/env python

from FourToOne import *
from graph import *
from HalfBB import G_HBB
import sys

def Simplify_Graph_HBB(g):
    ng = Simplify_Graph_I(g)
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W1(ng)    
    ng = Simplify_Graph_W2(ng)    
    ng = Simplify_Graph_V(ng)    
    ng = Simplify_Graph_O(ng)
    return ng

def find_pos(a,b):
    pos = []
    for ib in b:
        for i in xrange(0,len(a)):
            if ib == a[i]: pos.append(i)
    return pos
    
base = 100.
special_node = ["W", "V"]
def leg_cost_weight(name, SS):
    if not name[0] in special_node:
        cost = math.pow(base, comp)
        return cost
    
    sslist = list(SS)
    pos = find_pos(gmera[path[0]], sslist)
    print path[0], pos
    
def contract_gmera_chi(gmera, path, pr= False):
    if pr:
        print path
    S1 = set(gmera[path[0]])
    S2 = set(gmera[path[1]])
    SS = S1&S2
    
    St = (S1|S2)-SS
    max_leg = len(St)
    comp = len(St)+len(SS) 
    max_comp = comp
    base = 100.
    cost = math.pow(base, comp)

    if pr:
        print S1, S2
        print str(path[1])+',', len(St), len(St)+len(SS), St
        
    for p in path[2:]:
        S1 = St
        S2 = set(gmera[p])
        SS = S1&S2
        St = (S1|S2)-SS        
        comp = len(St)+len(SS)
        max_leg = max(max_leg, len(St))
        max_comp = max(max_comp, comp)
        cost = cost + math.pow(base, comp)
        if pr:
            print str(p)+',', len(St), len(St)+len(SS), St
    if len(St) != 0:
        print St
        return -1,-1
    return max_leg, max_comp, cost

def get_path_cost_chi(gx, pr=False):
    ge,gm=GC2GE(gx)
    graph = gnode(gm)
    ks = graph.keys()
    ks.sort()
    for k in ks:
        graph[k].sort()    
        
    max_leg = 1000; max_comp = 10000
    max_cost = 1e200
    for i in xrange(0,len(ks)):
        for j in xrange(i+1, len(ks)):
            start = ks[i]; end = ks[j]
            paths=find_all_path_pool(graph,start,end)
            for p in paths:
                leg,comp,cost=contract_gmera_chi(gm, p, pr)
                if leg < 0:
                    print 'error contracting'
                    sys.exit(-1)
                if (cost < max_cost):
                    max_leg = leg
                    max_comp = comp
                    path = p
                    max_cost = cost
                    if pr: print leg, comp, p, cost
                if (cost == max_cost):
                    if (leg < max_leg):
                        max_leg = leg
                        path = p
                        if pr: print leg, comp, p, cost
    return max_leg, max_comp, path, max_cost

G=G_HBB
OO = {"OO":[("I_8", 1), ("I_8", 2)],
      "I_8":[("OO", 1), ("OO", 2)]}

i=8
oo = eval(oo2_tmpl % {"1":i, "2":i+1})
gg = Simplify_Graph_HBB(Add_Graph(oo, G))
check_graph(gg)

max_leg, max_comp, path, max_cost=get_path_cost_chi(gg)
print max_leg, max_comp, path, max_cost
sys.exit(0)
