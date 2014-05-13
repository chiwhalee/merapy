#!/usr/bin/env python
#from sets import Set
import sys, math

def find_all_path_pool(g, start, end, pool=set([]), path=[]):
    path = path+[start]
    if start==end:
        if len(path) == len(g.keys()):
            return [path]
        return []
    if not g.has_key(start):
        return []
    paths = []
    path_nodes = set(path)
    new_nodes = set(g[start])-path_nodes
    pool = pool | new_nodes
    if len(pool)==0:
        return []
    while len(pool)>0:
        node = pool.pop()
        newpaths = find_all_path_pool(g, node, end, pool, path)
        for newpath in newpaths:
            paths.append(newpath)
    return paths
        
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
    S1 = set(gmera[path[0]])
    S2 = set(gmera[path[1]])
    SS = S1&S2
    St = (S1|S2)-SS
    max_leg = len(St)
    max_comp = len(St)+len(SS) 
    base = 100.
    cost = 0.
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

def check_graph(g):
    res = True
    for k in g.keys():
        e = 0
        for v in g[k]:
            e += 1
            nk,ne = g[v[0]][v[1]-1]
            if k != nk :
                res = False
                print 'graph not consistent', k, e, nk, ne
    return res

def Add_Graph(gx1, gx2):
    k1s = set(gx1.keys())
    k2s = set(gx2.keys())
    keys_comm = k1s & k2s    
    new_keys = k1s ^ k2s
    
    ng = {}
    for k1 in (k2s-k1s):
        nv = []
        for v in gx2[k1]:
            if (v[0] not in keys_comm) :
                nv.append(v)
            else:
                nv.append(gx1[v[0]][v[1]-1])
        ng[k1] = nv
        
    for k1 in (k1s-k2s):
        nv = []
        for v in gx1[k1]:
            if (v[0] not in keys_comm):
                nv.append(v)
            else:
                nv.append(gx2[v[0]][v[1]-1])
        ng[k1] = nv
        
    return ng

def Simplify_Graph(g):
    ng = Simplify_Graph_I(g)
#    print_graph(ng)    
    ng = Simplify_Graph_U(ng)    
    ng = Simplify_Graph_W3(ng)    
    ng = Simplify_Graph_O(ng)

    return ng

def Combine_OO(g1,g2):
    ng = {}
    o1 = g1['OO']
    o2 = g2['OO']

    leg1 = [0]*100
    leg2 = [0]*100
    len1 = len(o1)/2
    len2 = len(o2)/2    
    for i in xrange(1, len1+1):
        leg1[i] = i
        leg1[i+len1] = i+len1+len2

    for i in xrange(1, len2+1):
        leg2[i] = i+len1
        leg2[i+len2] = i+len1*2+len2

    noo = []
    for i in xrange(0, len1):
        noo.append(o1[i])
    for i in xrange(0, len2):
        noo.append(o2[i])
    for i in xrange(len1, len1+len1):
        noo.append(o1[i])
    for i in xrange(len2, len2+len2):
        noo.append(o2[i])
        
    ng['OO'] = noo
    for k in g1.keys():
        if k == "OO": continue
        noo = []
        oo = g1[k]
        for kk in oo:
            if kk[0] == 'OO': 
                nleg = leg1[kk[1]]
                noo.append((kk[0],nleg))
            else:
                noo.append(kk)                
        ng[k] = noo
        
    for k in g2.keys():
        if k == "OO": continue
        noo = []
        oo = g2[k]
        for kk in oo:
            if kk[0] == 'OO': 
                nleg = leg2[kk[1]]
                noo.append((kk[0],nleg))
            else:
                noo.append(kk)                
        ng[k] = noo
    
    return ng
#     ng = g
#     I_V = set([])
#     nng = {}
#     for k in ng.keys():
#         if k == "OO":
#             ks = "OO"
#             s = "0"
#         else:
#             ks,s = k.split("_", 1)
#         if (ks == "OO"):
#             s = True
#             for v in ng[k]:
#                 if v[0] != k:
#                     s = False
#                     break
#             if not s:
#                 nng[k] = ng[k]
#         else:
#             nng[k] = ng[k]
#     nk = []
#     for k in nng.keys():
#         if k == "OO":
#             ks = "OO"
#             e = "0"
#         else:
#             ks,e = k.split("_",1)
#         if ks == "OO":
#             nk.append(k)
#     nk.sort(cmp)
#     i = 0
#     l = len(nk)
#     gv = {}
#     ov = [(0,0)]*2*l
#     for k in nk:
#         i+= 1
#         nv = [("OO", i,), ("OO", i+l)]
#         gv[k] = nv
#         ov[i-1] = (k, 1)
#         ov[i+l-1] = (k,2)
#     gv["OO"] = ov
#    
#     return Add_Graph(gv, nng)

def Simplify_Graph_O(g):
#    ng = Simplify_Graph_W(g)
#    ng = Simplify_Graph_W4(g)
    ng = g
    I_V = set([])
    nng = {}
    for k in ng.keys():
        ks,s = k.split("_", 1)
        if (ks == "O"):
            s = True
            for v in ng[k]:
                if v[0] != k:
                    s = False
                    break
            if not s:
                nng[k] = ng[k]
        else:
            nng[k] = ng[k]
    nk = []
    for k in nng.keys():
        ks,e = k.split("_",1)
        if ks == "O":
            nk.append(k)
    nk.sort(cmp)
    i = 0
    l = len(nk)
    gv = {}
    ov = [(0,0)]*2*l
    for k in nk:
        i+= 1
        nv = [("OO", i,), ("OO", i+l)]
        gv[k] = nv
        ov[i-1] = (k, 1)
        ov[i+l-1] = (k,2)
    gv["OO"] = ov
    
    return Add_Graph(gv, nng)

def cmp(k1,k2):
    k,e1 = k1.split("_",1)
    e1 = int(e1)
    k,e2 = k2.split("_",1)
    e2 = int(e2)
    return e1-e2

def Simplify_Graph_W3(g):
#    ng = Simplify_Graph_U(g)
    ng = g
    I_V = set([])
    
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "W"):
            kp = "Wp_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            k3,e3 = ng[k][2]
            if (k1==kp) and (e1==2) and (k2==kp) and \
                    (e2==3) and (k3==kp) and (e3==4):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "Wp_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (kp, 4), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2), (k, 3)],
            Is:[(kp,1), (k,4)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng

def Simplify_Graph_W4(g):
#    ng = Simplify_Graph_U(g)
    ng = g
    I_V = set([])
    
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "W"):
            kp = "Wp_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            k3,e3 = ng[k][2]
            k4,e4 = ng[k][3]
            if (k1==kp) and (e1==2) and (k2==kp) and \
                    (e2==3) and (k3==kp) and (e3==4) and \
                    (k4==kp) and (e4==5):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "Wp_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (kp, 4), (kp, 5), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2), (k, 3), (k,4)],
            Is:[(kp,1), (k,5)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng

def Simplify_Graph_W1(g):
    ng = g
    I_V = set([])
    
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "W1"):
            kp = "W1p_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            if (k1==kp) and (e1==2) and (k2==kp) and (e2==3):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "W1p_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2)],
            Is:[(kp,1), (k,3)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng    

def Simplify_Graph_W2(g):
    ng = g
    I_V = set([])
    
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "W2"):
            kp = "W2p_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            if (k1==kp) and (e1==2) and (k2==kp) and (e2==3):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "W2p_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2)],
            Is:[(kp,1), (k,3)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng    

def Simplify_Graph_W(g):
#    ng = Simplify_Graph_U(g)
    ng = g
    I_V = set([])
    
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "W"):
            kp = "Wp_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            if (k1==kp) and (e1==2) and (k2==kp) and (e2==3):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "Wp_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2)],
            Is:[(kp,1), (k,3)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng    
    
def Simplify_Graph_U(g):
#    ng = Simplify_Graph_V(g)
    ng = g
    
    I_V = set([])
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "U"):
            kp = "Up_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            if (k1==kp) and (e1==3) and (k2==kp) and (e2==4):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "Up_"+s
        I1s = "I1_"+s
        I2s = "I2_"+s
        gv = {
            k:[(kp,3), (kp,4), (I1s, 2), (I2s,2)],
            kp:[(I1s,1), (I2s,1), (k,1), (k,2)],
            I1s:[(kp,1), (k,3)],
            I2s:[(kp,2), (k,4)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng
    
    
def Simplify_Graph_V(g):
#    ng = Simplify_Graph_I(g)
    ng = g
    
    I_V = set([])
    for k in ng.keys():
        ks,s = k.split("_",1)
        if (ks == "V"):
            kp = "Vp_"+s
            k1,e1 = ng[k][0]
            k2,e2 = ng[k][1]
            if (k1==kp) and (e1==2) and (k2==kp) and (e2==3):
                I_V.add(k)
    
    for k in I_V:
        ks,s = k.split("_",1)
        kp = "Vp_"+s
        Is = "I_"+s
        gv = {
            k:[(kp,2), (kp,3), (Is, 2)],
            kp:[(Is,1), (k,1), (k,2)],
            Is:[(kp,1), (k,3)]
            }
        ng = Add_Graph(gv, ng)
    ng = Simplify_Graph_I(ng)
    return ng

def print_graph(g, prefix=""):
    keys = g.keys()
    keys.sort()
    for k in keys:
        print prefix, k, g[k]
        
def Simplify_Graph_I(g):
    keys = set(g.keys())
    
    I_keys = set([])
    for k in keys:
        if k[0] == "I":
            I_keys.add(k)

    ng = {}
    for k in keys - I_keys:
        nv = []
        for v in g[k]:
            if (v[0] not in I_keys):
                nv.append(v)
            else:
                e = v[1]-1
                if (e==0): e = 1
                else: e = 0
                nv.append(g[v[0]][e])
        ng[k] = nv
        
    return ng

def Break_Graph(g):
    keys = set(g.keys())
    key_class = []
    while (len(keys)>0):
        kp = set([])
        kc = set([])
        k = keys.pop()
        kp.add(k)
        while (len(kp) > 0):
            k = kp.pop()
            kc.add(k)
            for v in g[k]:
                if (v[0] not in kc):
                    kp.add(v[0])
                    keys.discard(v[0])
        key_class.append(kc)

    gs = []
    for kc in key_class:
        ng = {}
        for k in kc:
            ng[k] = g[k]
        gs.append(ng)
    return gs
#     for kc in key_class:
#         print kc
#         print 
#         for v in g[k]:
#             kc.add(v[0])
#             keys.discard(v[0])

def ChangeNodeName(g):
    ng = {}
    if g.has_key('OO') and len(g['OO']) == 6:
        for k in g.keys():
            nvs = []
            for v in g[k]:                
                if v[0] == "OO":
                    nv = ""
            
    
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
    return gx4

def get_path_cost(gx, pr=False):
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
                leg,comp,cost=contract_gmera(gm, p)#, True)
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
                        
#                 if (comp == max_comp):
#                     if (leg < max_leg):
#                         max_leg = leg
#                         path = p
#                 elif (comp < max_comp):
#                     max_leg = leg
#                     max_comp = comp
#                     path = p
            
#                 if (leg<=max_leg and comp<=max_comp):
#                     if pr: print leg, comp, p, cost
                    
    return max_leg, max_comp, path, max_cost

def Output_Fortran(ii, Gname, Oname, gg, jj=-1):
    ge,gm = GC2GE(gg)
    gm_nodes = {}
    
    j = 0    
    ks = gm.keys()
    for k in ks:
        j += 1
        name=k.split("_",1)[0]
        gm_nodes[k] = j
    order = []
    leg,comp, path, cost=get_path_cost(gg)#,pr=True)
    for k in path:
        order.append(gm_nodes[k])
    
    print_graph(gg, "!!$ ")
    if jj < 0:
        print "   %s(%d)%%nNode=%d" % (Gname, ii,len(ks))
        print "   %s(%d)%%Nodes=-1" % (Gname, ii)
        print "   %s(%d)%%Edges=-1" % (Gname, ii)
    else:
        print "   %s(%d,%d)%%nNode=%d" % (Gname, ii,jj,len(ks))
        print "   %s(%d,%d)%%Nodes=-1" % (Gname, ii,jj)
        print "   %s(%d,%d)%%Edges=-1" % (Gname, ii,jj)
        
    j = 0
    for k in ks:
        j=j+1
        name=k.split("_",1)[0]        
        edges = ",".join([str(x) for x in gm[k]])
        leng = len(gm[k])
        if (name == "OO" and leng == 6): name = "OOO"
        if (name == "OO" and leng == 2): name = "O"
        if (name == "oo" and leng == 6): name = "ooo"
        if jj < 0:
            print '''   %s(%d)%%Names(%d)="%s" !%s''' % (Gname, ii,j, name, k) 
            print "   %s(%d)%%Nodes(%d)=%d" % (Gname, ii,j,leng)
            print "   %s(%d)%%Edges(1:%d,%d)=(/%s/)" % (Gname, ii,leng,j,edges)
        else:
            print '''   %s(%d,%d)%%Names(%d)="%s" !%s''' % (Gname, ii,jj,j, name, k) 
            print "   %s(%d,%d)%%Nodes(%d)=%d" % (Gname, ii,jj,j,leng)
            print "   %s(%d,%d)%%Edges(1:%d,%d)=(/%s/)" % (Gname, ii,jj,leng,j,edges) 
    print "!$$ %d/%d %20.14G" % (leg,comp, cost)
    if jj<0:
        print "   %s(1:%d, %d)=(/%s/)" % (Oname, len(ks), ii, 
                                          ",".join([str(x) for x in order]))
    else:
        print "   %s(1:%d, %d,%d)=(/%s/)" % (Oname, len(ks), ii, jj,
                                          ",".join([str(x) for x in order]))
        
    print


oo1_tmpl = """{
    "oo_%(1)s":[("I_%(1)s",1), ("I_%(1)s", 2)],
    "I_%(1)s":[("oo_%(1)s", 1), ("oo_%(1)s", 2)]
}"""

oo2_tmpl =  """{
    "oo_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("oo_%(1)s_%(2)s", 1), ("oo_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("oo_%(1)s_%(2)s", 2), ("oo_%(1)s_%(2)s", 4)]
    }
"""

oo2_Ntmpl =  """{
    "%(name)s_%(1)s_%(2)s":[("I_%(1)s",1), ("I_%(2)s", 1), ("I_%(1)s", 2), ("I_%(2)s", 2)],
    "I_%(1)s":[("%(name)s_%(1)s_%(2)s", 1), ("%(name)s_%(1)s_%(2)s", 3)],
    "I_%(2)s":[("%(name)s_%(1)s_%(2)s", 2), ("%(name)s_%(1)s_%(2)s", 4)]
    }
"""

oo3_tmpl =  """{
    "oo_%(1)s_%(2)s_%(3)s":[("I_%(1)s",1), ("I_%(2)s",1), ("I_%(3)s",1), 
                            ("I_%(1)s",2), ("I_%(2)s",2), ("I_%(3)s",2)],
    "I_%(1)s":[("oo_%(1)s_%(2)s_%(3)s",1), ("oo_%(1)s_%(2)s_%(3)s",4)],
    "I_%(2)s":[("oo_%(1)s_%(2)s_%(3)s",2), ("oo_%(1)s_%(2)s_%(3)s",5)],
    "I_%(3)s":[("oo_%(1)s_%(2)s_%(3)s",3), ("oo_%(1)s_%(2)s_%(3)s",6)]
    }
"""
