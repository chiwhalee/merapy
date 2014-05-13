#!/usr/bin/env python
from sets import Set
import sys

# g_x={
#      1:[(8,2), (4,3), (12,1)],
#      2:[(4,4), (5,3), (12,2)],
#      3:[(5,4), (9,3), (12,3)],
#      4:[(11,4), (11,5), (1,2), (2,1)],
#      5:[(11,6), (7,4), (2,2), (3,1)],
#      6:[(8,3),  (9,2), (11,1), (11,2)],
#      7:[(9,3),  (10,2), (11,3), (5,2)],
#      8:[(12,4), (1,1), (6,1)],
#      9:[(12,5), (6,2), (7,1)],
#     10:[(12,6), (7,2), (3,2)],
#     11:[(6,3), (6,4), (7,3), (4,1), (4,2), (5,1)],
#     12:[(1,3), (2,3), (3,3), (8,1), (9,1), (10,1)]
#     }

# #Binary 
# g_mera={
#     1:[1,2,3],
#     2:[4,5,6],
#     3:[7,8,9],
#     4:[10,11,2,4],
#     5:[12,13,5,7],
#     6:[14,15,16,17],
#     7:[18,19,20,13],
#     8:[21,1,14],
#     9:[22,15,18],
#     10:[23,19,8],
#     11:[16,17,20,10,11,12],
#     12:[3,6,9, 21,22,23]
#     }

# # 4 to 1
# g_mera = {
#     1:[1,2,3,4,  5],
#     2:[6,7,8,9,  10],
#     3:[11,12,13,14,  15],
#     4:[16,17,  4,6],
#     5:[18,19,  9,11],
#     6:[20,21,  16,22],
#     7:[23,24,  25,19],
#     8:[26,     1,2,3,20],
#     9:[27,     21,28,29,23],
#    10:[30,     24,12,13,14],
#    11:[22,28,  17,7],
#    12:[29,25,  8,18],
#    13:[5,10,15,   26,27,30]
#     }

# g_mera = {
#     1:[1,2,3,4,  5],
#     2:[6,7,8,9,  10],
#     3:[11,12,13,14,  15],
#     4:[16,17,   4,6],
#     5:[18,19,   9,11],
#     6:[20,21,   22,17],
#     7:[23,24,   18,25],
#     8:[26,      1,2,27,20],
#     9:[28,      21,7,8,23],
#     10:[29,     24,30,13,14],
#     11:[27,22,  3,16],
#     12:[25,30,  19,12],
#     13:[5,10,15,  26,28,29]
#     }

# #Ternary 
# g_mera = {
#     1:[1,2,3,      4],
#     2:[5,6,7,      8],
#     3:[9,10,11,    12],
#     4:[13,14,      3,5],
#     5:[15,16,      7,9],
#     6:[17,18,      19,14],
#     7:[20,21,      15,22],
#     8:[23,         1,24,17],
#     9:[25,         18,6,20],
#     10:[26,        21,27,11],
#     11:[24,19,     2,13],
#     12:[22,27,     16,10],
#     13:[4,8,12,    23,25,26]
#     }


# g_mera = {
#     1:[1,2,3,      4],
#     2:[5,6,7,      8],
#     3:[9,10,11,    12],
#     4:[13,14,      3,5],
#     5:[15,16,      7,9],
#     6:[17,18,      13,19],
#     7:[20,21,      22,23],
#     8:[24,         1,2,17],
#     9:[25,         18,26,20],
#     10:[27,        21,10,11],
#     11:[19,26,     14,6],
#     12:[22,23,     15,16],
#     13:[4,8,12,    24,25,27]
#     }


def find_all_path(g, start, end, path=[]):
    path = path+[start]
    if (start==end):
        if len(path) == len(g.keys()):
            return [path]
        return []
    if not g.has_key(start):
        return []
    paths = []
    for node in g[start]:
        if node not in path:
            newpaths = find_all_path(g, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths

def gnode(gmera):
    g = {}
    for k in gmera.keys():
        if not g.has_key(k): 
            g[k] = []
        for leg in gmera[k]:
            for kp in gmera.keys():
                if kp != k:
                    if (leg in gmera[kp]) and (not kp in g[k]):
                        g[k].append(kp)
    return g

def contract_gmera(gmera, path, pr= False):
    if pr:
        print path
    S1 = Set(gmera[path[0]])
    S2 = Set(gmera[path[1]])
    SS = S1&S2
    St = (S1|S2)-SS
    max_leg = len(St)
    max_comp = len(St)+len(SS) 
    if pr:
        print S1, S2
        print str(path[1])+',', len(St), len(St)+len(SS), St
        
    for p in path[2:]:
        S1 = St
        S2 = Set(gmera[p])
        SS = S1&S2
        St = (S1|S2)-SS        
        max_leg = max(max_leg, len(St))
        max_comp = max(max_comp, len(St)+len(SS))
        if pr:
            print str(p)+',', len(St), len(St)+len(SS), St
    if len(St) != 0:
        print St
        return -1,-1
    return max_leg, max_comp

def GC2GE(gc):
    ge = {}
    gm = {}
    ne = 0
    for node in gc.keys():
        i = 0
        gm[node] = []
        for edge in gc[node]:
            i+=1
            if ge.has_key((node,i)):
                print "error"
            if not ge.has_key(edge):
                ne += 1                
                ge[(node,i)] = [edge, ne]
                gm[node].append(ne)
            else:
                ge[(node,i)] = [edge, ge[edge][1]]
                gm[node].append(ge[edge][1])

    return ge,gm


# g_x = {
#     1:[(6,2), (4,3), (11,1)],
#     2:[(4,4), (10,3), (11,2)],
#     3:[(10,4), (8,3), (11,3)],
#     4:[(9,3), (9,4), (1,2), (2,1)],
#     5:[(6,3), (7,2), (9,1), (9,2)],
#     6:[(11,4), (1,1), (5,1)],
#     7:[(11,5), (5,2), (10,1)],
#     8:[(11,6), (10,2), (3,2)],
#     9:[(5,3), (5,4), (4,1), (4,2)],
#     10:[(7,3), (8,2), (2,2), (3,1)],
#     11:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
#     }

def add_graph(gx1, gx2):
    k1s = gx1.keys()
    k1s.sort()

    gx3 = {}
    k2s = gx2.keys()
    k2s.sort()
    for k2 in k2s:
        gx3[k2+100] = []
        for v in gx2[k2]:
            gx3[k2+100].append((v[0]+100, v[1]))
#        print k2+100, gx3[k2+100]

#    print
    k3s = gx3.keys()
    k3s.sort()
        
    k1_link = k1s[-1]
    k3_link = k3s[-2]

    gx4 = {}
    for k1 in k1s:
        if k1 == k1_link:
            continue
        else:
            gx4[k1] = []
            for n,e in gx1[k1]:
                if n != k1_link:
                    gx4[k1].append((n,e))
                else:
                    nn,ee = gx3[k3_link][e-1]
                    gx4[k1].append((nn,ee))
#            print k1, gx4[k1]
    
    for k3 in k3s:
        if k3 == k3_link:
            continue
        else:
            gx4[k3] = []
            for n,e in gx3[k3]:
                if n != k3_link:
                    gx4[k3].append((n,e))
                else:
                    nn, ee = gx1[k1_link][e-1]
                    gx4[k3].append((nn,ee))
#            print k3, gx4[k3]

    return gx4

def get_path_cost(gx):
    ge,gm=GC2GE(gx)
    for k in gm.keys():
        print k, gm[k]
    graph = gnode(gm)

    ks = graph.keys()
    ks.sort()
    for k in ks:
        graph[k].sort()    

    print    
    max_leg = 1000; max_comp = 1000    
    for i in xrange(0,len(ks)):
        for j in xrange(i+1, len(ks)):
            start = ks[i]; end = ks[j]
            paths=find_all_path(graph,start,end)
            for p in paths:
                leg,comp=contract_gmera(gm, p)#, True)
                if leg < 0:
                    print 'error contracting'
                    sys.exit(-1)
                if ((comp<=max_comp)):# or (leg<=max_leg)):
                    max_comp = comp
                    if (leg<=max_leg):
                        max_leg = leg
                    print leg,comp, p

                        
from Four2One_graph import *

gx = gx1_3[4]

for i in xrange(-3,5):
    gx = gx1_2[i]
    get_path_cost(gx)    
    gx = gx1_3[i]
    get_path_cost(gx)
#     for j in xrange(-3,5):
#         gx = add_graph(gx1_3[i], gx1_2[j])
#         get_path_cost(gx)
sys.exit(0)
                       
ge,gm=GC2GE(gx)
for k in gm.keys():
    print k, gm[k]
graph = gnode(gm)

ks = graph.keys()
ks.sort()
for k in ks:
    graph[k].sort()    

print    
#path = [8, 3, 102, 2, 9, 4, 5, 7, 104, 106, 101, 1, 6, 103]
#contract_gmera(gm, path, pr=True)

#sys.exit(0)
# path=[8, 1, 4, 6, 11, 9, 13, 2, 12, 7, 5, 3, 10]
# contract_gmera(g_mera, path, pr=True)
# sys.exit(0)

# for k in ks:
#     print k, graph[k]
    
# sys.exit(0)
# print
# print

max_leg = 1000; max_comp = 1000    
for i in xrange(0,len(ks)):
    for j in xrange(i+1, len(ks)):
        start = ks[i]; end = ks[j]
        paths=find_all_path(graph,start,end)
        for p in paths:
            leg,comp=contract_gmera(gm, p)#, True)
            if leg < 0:
                print 'error contracting'
                sys.exit(-1)
            if ((comp<=max_comp)):# or (leg<=max_leg)):
                max_comp = comp
                if (leg<=max_leg):
                    max_leg = leg
                print leg,comp, p
