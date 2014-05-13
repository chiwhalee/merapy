#!/usr/bin/env python
#coding=utf8
"""


"""
import sys 
import os
import math
import pprint 
import networkx as nx
import pygraphviz as pgz
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint as pp


def cmp(k1,k2):
    k,e1 = k1.split("_",1)
    e1 = int(e1)
    k,e2 = k2.split("_",1)
    e2 = int(e2)
    return e1-e2

class Diagram(dict):
    def __init__(self, old_graph={}):
        """
            self is inited by __setitem__, update, update_key
        """
        dict.__init__(self)
        self.edges = {}
        self.connections = {}
        self.edge_info = {}
        self.NumOfEdges = 0
        self.Nodes = {}
        self.name = ""
        
        for k in old_graph.keys():
            self.__setitem__(k, old_graph[k])
        
    def __repr__(self):
        r = ""
        for k,v in self.connections.items():
            r = r+str(k)+": "+str(v)+"\r"
        return r
    
    def __repr__new(self):
        print "rrrr"
        #return self.connections.__repr__()
        return pprint.saferepr(self.connections)
    
    def __str__new(self):
        return self.connections.__str__()

    def __setitem__(self, key, value):
        self.connections[key] = value
        #print "updating ", key 
        self.update_key(key)
        #print "after update", self.edges[key]
        
    def __len__(self):
        return self.connections.__len__()

    def __getitem__(self, key):
        return self.connections[key]

    def __add__(self, other):
        d = Diagram()
        k1s = set(self.keys())
        k2s = set(other.keys())
        
        #keys in common will be abbondoned
        keys_comm = k1s & k2s
        new_keys = k1s ^ k2s

        for k1 in (k2s-k1s):
            nv = []
            for v in other[k1]:
                #v is a tuple, leg linked to k1
                if (v[0] not in keys_comm):
                    nv.append(v)
                else:
                    #用self[v[0]][v[1]-1] 替代原来的leg
                    nv.append(self[v[0]][v[1]-1])
            d[k1] = nv
        for k1 in (k1s-k2s):
            nv = []
            for v in self[k1]:
                if (v[0] not in keys_comm):
                    nv.append(v)
                else:
                    #用other[v[0]][v[1]-1] 替代原来的leg
                    nv.append(other[v[0]][v[1]-1])
            d[k1] = nv
        return d
    
    @staticmethod
    def node_property(node_name):
        #print 'nnn', node_name
        try:
            name, xpos = node_name.split('_')
        except :
            count = node_name.count('_')
            if count >= 2:
               name, xpos = node_name.split('_')[:2]
            elif count == 0:
                #issue: for 'O' xpos is not defined
                name = node_name
                xpos= '-1'

            else:
                raise ValueError(node_name)

        xpos = eval(xpos)*100
        #xpos = float(xpos)
        type_map = {'U':'unitary', 'V':'disentangler', 'o':'H', 'O':'rho', 'I':'I'}
        dual = False
        try:
            if node_name[1] == 'p':
                dual = True
        except: pass
        type = type_map[name[0]]
        ypos = 200.
        if name in ['U', 'V']: ypos = 100.
        if name in ['Up', 'Vp']: ypos = 300. 
        pos= (xpos, ypos)
        

        if 0:
            color_map = {"V":"yellow", "U":"green", "o":"red", "O":"red"}
            color_map['Vp'] = color_map['V']
            color_map['Up'] = color_map['U']
            shape_map = {'U':'square', 'Up':'square', 'V':'triangle', 'Vp':'invtriangle', 'I':'square'}
            #shape_map = {'U':'square', 'Up':'square', 'V':'hexagon', 'Vp':'hexagon', 'I':'square'}
            for i in G:
                node_ppt= Diagram.node_property(i)
                #x = eval(i.split('_')[1])
                #name = i[0]
                name = i.split('_')[0]
                node = G.node[i]
                #node['pos'] = nodep['pos']
                if color_map.has_key(name):
                    node["color"] = color_map[name] 
                if shape_map.has_key(name):
                    node["shape"] = shape_map[name] 
                node.update(node_ppt)
        propt = {
                'U':    {'shape':'square',      'color':'green'}, 
                'Up':   {'shape':'square',      'color':'green'}, 
                'V':    {'shape':'triangle',    'color':'yellow'}, 
                'Vp':   {'shape':'invtriangle', 'color':'yellow'}, 
                'o':    {'shape':'circle',      'color':'red'}, 
                'O':    {'shape':'circle',      'color':'red'}, 
                'I':    {'shape':'square',      'color':'black'}, 
                }
        propt['ooo'] = propt['oo']  = propt['o']
        propt['OOO'] = propt['OO']  = propt['O']
        #return vars()
        #res= {'name', ''}
        #return vars()
        temp = ['name', 'type', 'dual', 'pos']
        #temp = ['name', 'type', 'dual', ]
        variables= vars()
        res= {i:variables[i] for i in temp}
        res.update(propt[name])
        return res
        

    def has_key(self, key):
        return self.connections.has_key(key)
    
    def update(self):
        """
        更新edge, nodes 等信息
        """
        self.NumOfEdges = 0
        self.edges = {}
        nn = 1
        for key in self.connections.keys():
            self.update_key(key)
            self.Nodes[key] = nn
            nn = nn+1
        
    def update_key(self,key):
        """
            this meth mainly calc self.edges[key] using the info of self.connections[key]
        """

        #key,value looks like this:  'I_1', [('Vp_0', 2), ('V_0', 1)]
        value = self.connections[key]
        edges = []
        #print "eee", key, value
        for i in xrange(0, len(value)):
            
            ii = i+1
            v = value[i]
            node = v[0]
            edge = v[1]   #edge meas node's edge's leg
            
            try: 
                weight = v[2]
            except: 
                weight = 1
                k=key.split("_")[0]
                if (k in ["W1", "W2"]) and (ii == 3):
                    weight = 2
                if (k in ["W1p", "W2p"]) and (ii == 1):
                    weight = 2
                if (k == "V") and (ii in [1, 2]):
                    weight = 2
                if (k == "Vp") and (ii in [2, 3]):
                    weight = 2
            es = self.get_edges(node)
            #if node[0] == "V": 
            #print "\tiii", "key", node, "#leg", edge, "edges", es
            if len(es) == 0:
                self.NumOfEdges += 1
                en = self.NumOfEdges
                self.edge_info[en]=[weight,(node,edge)]
            else:
                en = es[edge-1]
                self.edge_info[en].append((node,edge))
            edges.append((en,weight))
        self.edges[key] = edges


    def show(self, prefix=""):
        
        #print "definition of diagram: " + self.name
        keys = self.connections.keys()
        keys.sort()
        for k in keys:
            v = self.connections[k]
            print prefix+k+": "+str(v)

    if 1: #edges and connections
        def show_edges(self):
            keys = self.edges.keys()
            keys.sort()
            for k in keys:
                edges = self.edges[k]
                estr = k+": "
                for e in edges:
                    estr = estr+ " " +str(e[0])
                print estr
                
        def show_edge_info(self):
            keys = self.edge_info.keys()
            keys.sort()
            for k in keys:
                print k, self.edge_info[k]
                
        def keys(self):
            return self.connections.keys()
        
        def nodes(self):
            return self.keys()
        
        def get_edges(self, node):
            """
                return a list such as: [1, 3, 12, 13]
            """
            try:
                return [e[0] for e in self.edges[node]]
            except:
                return []
        
        def get_all_connections(self, node):
            try:
                return self.connections[node]
            except:
                return []
            
        def get_connection(self, node, edge):
            try:
                v = self.connections[node][edge-1]
                return v[0],v[1]
            except:
                return "empty", -1000
        
    def weight(self, node, edge):
        return self.edges[node][edge][1]    
    
    def check(self):
        res = ""
        for k in self.connections.keys():
            e = 0
            for v in self.connections[k]:
                e += 1
                nk, ne = self.get_connection(v[0], v[1])
                if k != nk :
                    res +=  'graph not consistent ' + str((k, e)) + \
                           " is " + str((v[0], v[1])) + " but the latter is " + \
                           str((nk, ne)) + '\n'
        return res

    def node_rename(self, oldnode, newnode):
        ng = Diagram({})
        for k in self.keys():
            v = self[k]
            nv = []
            for n,e in v:
                if n==oldnode:
                    nv.append((newnode,e,))
                else:
                    nv.append((n,e,))
            if k==oldnode:
                ng[newnode]=nv
            else:
                ng[k]=nv
        ng.update()   #this line added by lzh
        return ng

    if 1:# Simplify_I, etc  ==================================================                                                  
        def Simplify_I(self, info=0):
            """
            remove all I_n
            Category theoretically that is identifying the identity tensor with a direct link
            """
            keys = set(self.connections.keys())
            I_keys = set([])
            nd = Diagram()
            for k in keys:
                # k is a string
                if k[0] == "I":
                    I_keys.add(k)

            for k in keys - I_keys:
                nv = []
                for v in self.connections[k]:
                    if (v[0] not in I_keys):
                        nv.append(v)
                    else:
                        #v[1] ranges from 1, 2 representing two legs of I
                        e = v[1]-1
                        if (e==0): e = 1
                        else: e = 0
                        nv.append(self.connections[v[0]][e])
                nd[k] = nv
            if info>0:
                print "\n after simplify ", "I "*20
                pp(nd.edges)
            return nd


        def Simplify_U(self, info=0):
            """
            identifying UU^+  with I\otimes I and then remove Is
            """
            ng = Diagram(self)
            #I_V stores all U that will be removed
            I_V = set([])
            for k in ng.keys():
                ks,s = k.split("_",1)
                if (ks == "U"):
                    kp = "Up_"+s
                    k1,e1 = ng.get_connection(k,1)
                    k2,e2 = ng.get_connection(k,2)                    
                    #when U and Up connect directly
                    if (k1==kp) and (e1==3) and (k2==kp) and (e2==4):
                        I_V.add(k)
            
            #print "IIII"
            #pp(ng.edges)
            #print I_V
            
            #I_V mabe somthing like  set(['U_2'])
            #replace U, Up with I, and finally remove I
            for k in I_V:
                ks,s = k.split("_",1)

                kp = "Up_"+s
                
                #left and right leg of U
                I1s = "I1_"+s  
                I2s = "I2_"+s
                gv = {
                    k:[(kp,3), (kp,4), (I1s, 2), (I2s,2)],
                    kp:[(I1s,1), (I2s,1), (k,1), (k,2)],
                    I1s:[(kp,1), (k,3)],
                    I2s:[(kp,2), (k,4)]
                    }
                gv = Diagram(gv)
                ng = ng+gv
            #print "eee"
            #pp(ng.edges); exit()
            ng = ng.Simplify_I()
            
            if info>0:
                print "\n after simplify ", "U "*20
                pp(ng.edges)
            
            return ng

        def Simplify_VNodeG(self, name1, name2, info=0):
            ng = self
            I_V = set([])
            for k in ng.keys():
                ks,s = k.split("_",1)
                if (ks == name1):
                    kp = name2+"_"+s
                    match = True
                    cons = ng.connections[k]
                    for en in xrange(0, len(cons)-1):
                        kc,ec = cons[en]
    #                    print k, kp, kc, en, ec
                        if not ((ec==en+2) and (kc==kp)):
                            match = False
                    if match:
                        I_V.add(k)
            
            for k in I_V:
                ks,s = k.split("_",1)
                kp = name2+"_"+s
                Is = "I_"+s
                kv = []
                kvp = [(Is,1,)]
                nleg = len(ng.connections[k])
                for i in xrange(1, nleg):
                    kv.append((kp,i+1,))
                    kvp.append((k,i,))
                kv.append((Is,2,))

                gv = {
                    k:kv,
                    kp:kvp,
                    Is:[(kp,1), (k,nleg)]
                    }
                gv = Diagram(gv)            
                ng = ng + gv
            ng = ng.Simplify_I()
            
            if info>0:
                print "\n after simplify ", "V "*20
                pp(ng.edges)

            return ng
            
        def Simplify_W1(self):
            """ for halfBB diagrams """
            return self.Simplify_VNodeG("W1", "W1p")
        
        def Simplify_W2(self):
            return self.Simplify_VNodeG("W2", "W2p")

        def Simplify_V(self):
            return self.Simplify_VNodeG("V", "Vp")

   
        def Simplify_O(self, info=0):
            """
            remove Os, or replace with OO, OOO,etc
            this implicitly comprimise an ascending process, i.e. oo is ascending to OO,  or likewise
            """
            #self.show()
            ng = self
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
            
            nk = []  #num of key that is "O"
            for k in nng.keys():
                ks,e = k.split("_",1)
                if ks == "O":
                    nk.append(k)
            
            nk.sort(cmp)
            i = 0
            l = len(nk)
            gv = {}
            ov = [(0,0)]*2*l

            
            Oname = "OO"
            if l==3:
                Oname="OOO"
            elif l==1:
                #Oname="O"  this line cause bug, modified by lzh to following, but in future and in some curcumstance, it should be change back
                #Oname = "OO"
                Oname = "O"

            for k in nk:
                i+= 1
                nv = [(Oname, i,), (Oname, i+l)]
                gv[k] = nv
                ov[i-1] = (k, 1)
                ov[i+l-1] = (k,2)
            gv[Oname] = ov
            gv = Diagram(gv)
            nng = gv+nng
            nng.update()
            if info>0:
                print "\n after simplify ", "O "*20
                pp(nng.edges)

            return nng

    if 1: #calc cost of paths
        def find_all_path(self, start, end, pool=set([]), path=[]):
            """
                here uses recursion to find all paths between two distinct nodes start and end, the
                paths should pass all nodes
                e.g. one path is somthing like this:
                ['OOO', 'V_0', 'V_1', 'U_2', 'Up_2', 'V_2', 'Vp_2', 'Vp_1', 'Vp_0', 'Up_1', 'ooo_2_3_4', 'U_1']
                if start='OOO' and end='U_1'
            
            """
            path = path+[start]
            if start==end:  #recursion return condition
                if len(path) == len(self.keys()):
                    return [path]
                return []
            if not self.connections.has_key(start):
                return []
            paths = []  #start from an empty
            path_nodes = set(path)

            cnodes = [] #nodes connected to start node
            for node in self.connections[start]:
                cnodes.append(node[0])            
            new_nodes = set(cnodes)-path_nodes
            pool = pool | new_nodes
            if len(pool)==0:
                return []
            while len(pool)>0:
                node = pool.pop()
                newpaths = self.find_all_path(node, end, pool, path)
                for newpath in newpaths:
                    paths.append(newpath)
            return paths
        
        def get_path_cost(self, path, pr=False, allow_open_leg=0):
            """
            returns:
                max_leg: num of open legs of two nodes which would contract into one node
                max_comp: comp may means complete, num of all legs of two nodes
                cost: cumulating all the costs of contractiong two nodes. this in effect count the number of max_comp
                
            """
            S1 = set(self.get_edges(path[0]))
            S2 = set(self.get_edges(path[1]))
            
            SS = S1&S2
            S_open = (S1|S2)-SS
            max_leg = len(S_open)   
            comp = len(S_open)+len(SS)
            max_comp = comp     
            
            weight = 1
            for e in S1|S2:
                weight = weight*self.edge_info[e][0]
    #        weight = 1
            
            base = 100. # base can be other values
            cost = math.pow(base, comp)*weight
            if pr:
                print S1, S2
                print str(path[1])+',', len(S_open), len(S_open)+len(SS), S_open
            
            for p in path[2:]:
                S1 = S_open  #note this line
                S2 = set(self.get_edges(p))
                SS = S1&S2
                S_open = (S1|S2)-SS        
                comp = len(S_open)+len(SS)
                max_leg = max(max_leg, len(S_open))
                max_comp = max(max_comp, comp)
                weight = 1
                for e in S1|S2:
                    weight = weight*self.edge_info[e][0]            

                cost = cost + math.pow(base, comp)*weight
                if pr:
                    print str(p)+',', len(S_open), len(S_open)+len(SS), S_open
            if len(S_open) != 0 and not allow_open_leg:  # for a closed graph, there should be no open leg in the end!!
                print "SSSS_open", S_open
                return -1,-1,-1
            return max_leg, max_comp, cost        

        def find_optimal_path(self, pr=False, allow_open_leg=0):
            keys = self.keys()
            keys.sort()
            
            max_leg = 1000; max_comp = 10000
            max_cost = 1e200
            for i in xrange(0,len(keys)):
                for j in xrange(i+1, len(keys)):  #transverse all starts and ends
                    start = keys[i]; end = keys[j]
                    paths = self.find_all_path(start, end)
                    for path in paths:  # transverse all paths for given start and end
                        leg, comp, cost = self.get_path_cost(path, pr=False, allow_open_leg=allow_open_leg)
                        #print "ppp", path; exit()
                        if (cost < max_cost):
                            max_leg = leg
                            max_comp = comp
                            optimal_path = path
                            max_cost = cost
                            if pr: print leg, comp, path, cost
                        if (cost == max_cost):
                            if (leg < max_leg):
                                max_leg = leg
                                optimal_path = path
                                if pr: print leg, comp, path, cost
                        
            return max_leg, max_comp, optimal_path, max_cost

    def Combine_Node(self, d, node, n1, n2):
        """
        d is another Diagram
        """

        ng = {}
        for k in self.connections.keys():
            ng[k] = self.connections[k]
        for k in d.connections.keys():
            ng[k] = d.connections[k]
        o1 = self.connections[node]
        o2 = d.connections[node]
        noo = []
        nn = 1
        for i in xrange(0,n1):
            noo.append(o1[i])
            ng[o1[i][0]][o1[i][1]-1] = (node, nn)
            nn = nn+1
            
        for i in xrange(0,n2):
            noo.append(o2[i])
            ng[o2[i][0]][o2[i][1]-1] = (node, nn)
            nn = nn+1

        for i in xrange(n1,len(o1)):
            noo.append(o1[i])
            ng[o1[i][0]][o1[i][1]-1] = (node, nn)
            nn = nn+1

        if i in xrange(n2,len(o2)):
            noo.append(o2[i])                    
            ng[o2[i][0]][o2[i][1]-1] = (node, nn)
            nn = nn+1
            
        ng[node] = noo
        ng = Diagram(ng)
        ng.update()
        return ng
        
    def toFortran(self, Gname, Oname, i1, i2=-10000):
        self.show(prefix="!!$ ")
        ks = self.connections.keys()
        if i2<-1000:
            print "   %s(%d)%%nNode=%d" % (Gname, i1, len(ks))
            print "   %s(%d)%%Nodes=-1" % (Gname, i1)
            print "   %s(%d)%%Edges=-1" % (Gname, i1)
        else:
            print "   %s(%d,%d)%%nNode=%d" % (Gname, i1,i2, len(ks))
            print "   %s(%d,%d)%%Nodes=-1" % (Gname, i1,i2)
            print "   %s(%d,%d)%%Edges=-1" % (Gname, i1,i2)
        for k in ks:
            j = self.Nodes[k]
            name = k.split("_",1)[0]
            es = self.get_edges(k)
            leng = len(es)
            es = ",".join([str(x) for x in es])
	    if i2<-1000:
                print '''   %s(%d)%%Names(%d)="%s" !%s''' %(Gname, i1, j, name, k)
                print '''   %s(%d)%%Nodes(%d)=%d''' % (Gname, i1, j, leng)
                print '''   %s(%d)%%Edges(1:%d,%d) = (/%s/)''' % (Gname, i1, leng, j, es)
	    else:
                print '''   %s(%d,%d)%%Names(%d)="%s" !%s''' %(Gname, i1, i2, j, name, k)
                print '''   %s(%d,%d)%%Nodes(%d)=%d''' % (Gname, i1, i2, j, leng)
                print '''   %s(%d,%d)%%Edges(1:%d,%d) = (/%s/)''' % (Gname, i1, i2, leng, j, es)

        order = []
        leg,comp, path, cost= self.find_optimal_path()
        for k in path:
            order.append(self.Nodes[k])
        print "!$$ %d/%d %20.14G" % (leg,comp, cost)            
        if i2 < -1000:
            print "   %s(1:%d, %d)=(/%s/)" % (Oname, len(ks), i1, 
                                          ",".join([str(x) for x in order]))
        else:
            print "   %s(1:%d, %d,%d)=(/%s/)" % (Oname, len(ks), i1, i2, 
                                          ",".join([str(x) for x in order]))
    
    def toGraphics(self, Gname, Oname, i1, i2=-10000, weight=None, calc_order=True):
        from merapy.graphics import Graphics

        self.show(prefix="# ")
        ks = self.connections.keys()
        print "%s[%s]=Graphics()" %(Gname, i1)
        print "%s[%s].graph_name='%s[%s]'" %(Gname, i1, Gname, i1)
        if i2<-1000:
            print "%s[%s].size=%d" % (Gname, i1, len(ks))
            print "%s[%s].nodes[:]=-1" % (Gname, i1)
            print "%s[%s].edges[:, :]=-1" % (Gname, i1)
        else:
            print "%s[%d,%d].size=%d" % (Gname, i1,i2, len(ks))
            print "%s[%d,%d].nodes[:]=-1" % (Gname, i1,i2)
            print "%s[%d,%d].edges[:, :]=-1" % (Gname, i1,i2)
        for k in ks:
            #j = self.Nodes[k]
            j = self.Nodes[k]-1
            name = k.split("_",1)[0]
            es = self.get_edges(k)
            es= [i-1 for i in es]

            leng = len(es)
            es = ",".join([str(x) for x in es])
            if i2<-1000:
                    print '''%s[%s].names[%d]="%s"  #%s''' %(Gname, i1, j, name, k)
                    print '''%s[%s].nodes[%d]=%d''' % (Gname, i1, j, leng)
                    print '''%s[%s].edges[0:%d,%d] = [%s]\n''' % (Gname, i1, leng, j, es)
            else:
                    print '''%s[%d,%d].names[%d]="%s"  #%s''' %(Gname, i1, i2, j, name, k)
                    print '''%s[%d,%d].nodes[%d]=%d''' % (Gname, i1, i2, j, leng)
                    print '''%s[%d,%d].edges[0:%d,%d] = [%s]\n''' % (Gname, i1, i2, leng, j, es)
        
        
        if calc_order:
            order = []
            leg, comp, path, cost= self.find_optimal_path()
            
            for k in path:
                order.append(self.Nodes[k])

            order = [i-1 for i in order]
            print "# %d/%d %20.14G" % (leg,comp, cost)            
            if weight is not None:
                print '''%s[%s].weight = %1.15f  ''' %(Gname, i1, weight)
            if i2 < -1000:
                #print "%s[%s][0:%d]=[%s]" % (Oname, i1, len(ks), 
                #                              ",".join([str(x) for x in order]))
                print '''%s[%s].contract_order=[%s]  \n''' %(Gname, i1, ",".join([str(x) for x in order]))
                
            else:
                print "%s[0:%d, %d,%d]=[%s]" % (Oname, len(ks), i1, i2, 
                                              ",".join([str(x) for x in order]))
    
    def toGraphics2(self, Gname, Oname, i1, i2=-10000, calc_order=True):
        from merapy.graphics import Graphics
        
        G_dic = {}
        G = Graphics()

        self.show(prefix="# ")
        ks = self.connections.keys()

        G_dic[i1].graph_name = Gname + "['%s']"%i1
        #print "%s[%s].graph_name='%s[%s]'" %(Gname, i1, Gname, i1)
        if i2<-1000:
            print "%s[%s].size=%d" % (Gname, i1, len(ks))
            print "%s[%s].nodes[:]=-1" % (Gname, i1)
            print "%s[%s].edges[:, :]=-1" % (Gname, i1)
        else:
            print "%s[%d,%d].size=%d" % (Gname, i1,i2, len(ks))
            print "%s[%d,%d].nodes[:]=-1" % (Gname, i1,i2)
            print "%s[%d,%d].edges[:, :]=-1" % (Gname, i1,i2)
        for k in ks:
            #j = self.Nodes[k]
            j = self.Nodes[k]-1
            name = k.split("_",1)[0]
            es = self.get_edges(k)
            es= [i-1 for i in es]

            leng = len(es)
            es = ",".join([str(x) for x in es])
            if i2<-1000:
                    print '''%s[%s].names[%d]="%s"  #%s''' %(Gname, i1, j, name, k)
                    print '''%s[%s].nodes[%d]=%d''' % (Gname, i1, j, leng)
                    print '''%s[%s].edges[0:%d,%d] = [%s]''' % (Gname, i1, leng, j, es)
            else:
                    print '''%s[%d,%d].names[%d]="%s"  #%s''' %(Gname, i1, i2, j, name, k)
                    print '''%s[%d,%d].nodes[%d]=%d''' % (Gname, i1, i2, j, leng)
                    print '''%s[%d,%d].edges[0:%d,%d] = [%s]''' % (Gname, i1, i2, leng, j, es)
        
        if calc_order:
            order = []
            leg, comp, path, cost= self.find_optimal_path()
            
            for k in path:
                order.append(self.Nodes[k])

            order = [i-1 for i in order]
            print "# %d/%d %20.14G" % (leg,comp, cost)            
            if i2 < -1000:
                print "%s[%s][0:%d]=[%s]" % (Oname, i1, len(ks), 
                                              ",".join([str(x) for x in order]))
            else:
                print "%s[0:%d, %d,%d]=[%s]" % (Oname, len(ks), i1, i2, 
                                              ",".join([str(x) for x in order]))
    
    def to_networkx(self, directed=False, remove_O=False):
        """
        ops 'O', 'OO', etc are omitted!
        """
        #G=nx.MultiDiGraph()
        #print 'eeee'; exit()
        if not directed:
            G=nx.MultiGraph()
            for k, v in self.connections.iteritems():
                leg0 = 1
                for k1, leg1 in v:
                    e = (k, k1, (k, k1, leg0, leg1))
                    e_reverse = (k, k1, (k1, k, leg1, leg0))
                    if not G.has_edge(*e_reverse):
                        G.add_edge(*e)
                    leg0 += 1 
        else:
            raise NotImplemented
            G=nx.MultiDiGraph()
            for k, v in self.connections.iteritems():
                leg0 = 1
                for k1, leg1 in v:
                    e = (k, k1, (k, k1, leg0, leg1))
                    e_reverse = (k, k1, (k1, k, leg1, leg0))
                    if not G.has_edge(*e_reverse):
                        G.add_edge(*e)
                    leg0 += 1 
            #print G.edges()
        
        if remove_O:
            print "ops 'O', 'OO', etc are omitted!"
            nodes= G.nodes()
            nodes= [i for i in nodes if i[0] not in ['O']  ]
            G = G.subgraph(nodes)

        for i in G:
            node = G.node[i]
            node_ppt= Diagram.node_property(i)
            node.update(node_ppt)

        return G
    

    #def plot(self, path=None, remove_O=True, graph_attr={}, node_attr={}, edge_attr={}):
    def plot(self, path=None, remove_O=True, **attr):
        G = self.to_networkx(remove_O=remove_O)
        nodes=G.nodes()

        if remove_O:
            nodes= [i for i in nodes if i[0] not in ['O']  ]
            G = G.subgraph(nodes)
        
        g = nx.to_agraph(G)

        node_class= {}
        gsub = {}
        if 1:
            temp = ['V_', 'U_', 'I', 'o', 'Up', 'Vp', ]
            #temp = ['V_',  'I', 'o', 'Vp', ]
            #temp.reverse()
            #rank = [1, 2, 3, 4, 5]
            for i in temp:
                node_class[i] = [j for j in nodes if j[0:len(i)] == i]
                gsub[i] = g.add_subgraph(node_class[i])
                rank = 'same'
                if i  == 'Vp': rank = 'max'
                if i  == 'V_': rank = 'min' 
                #if i  == 'U_': rank = None
                gsub[i].graph_attr['rank'] = rank
                #gsub[i].node_attr['group'] = i + '_group'
                if i == 'I': 
                    pass
                    #gsub[i].node_attr['group'] = 'I_group'


                #gsub[i].graph_attr['rank'] =  temp.index(i) + 1# 'same'

        if 0:
            properties = [('disentangler', False), ('unitary', False), ('disentangler', True), ('unitary', True)]
            for ppt in properties:
                node_class[ppt] = [j for j in nodes if g[j]['type'] == ppt[0]  and g[j]['dual'] == ppt[1]]
                gsub[i] = g.add_subgraph(node_class[i])
                gsub[i].graph_attr['rank'] = 'same'
        
        graph_attr_default = {
                'dim':2, 
                #'size': (1, 1)    #not working
                #'page':     (2, 2)  #not working
                }
        graph_attr_adj = {
                'layout':   'neato', #neato, dot, twopi, circo, fdp, nop
                'splines':  'curved',    #'line', 
                'ranksep':  0.5, #(0.3, 'equally'),   #equally is not usefull
                'nodesep':  1.0,   #not work
                'landscape':    False, 
                'rankdir':  'TB',   #'LR'
                #'pin':      True
                #'sep':      1, 
                #'minlen':   3,   #dot only
                #neato only
                #'mode':     'hier',   'levelsgap': 1.0,  #neato only
                }


        node_attr_default = {}
        node_attr_adj = {
                #'style':    'invis', 
                #'ordering': 'in',      #default is ""   not work
                }
        
        if 1:
            edge_attr_default = {}
            edge_attr_adj = {
                #'style':    'invis'

                #neato only
                'len':      1.0,    #this works!
                #dot only
                'minlen':   1, 
                    }
        
        #for a in [g.graph_attr, g.node_attr, g.edge_attr]:
        variables= vars()
        for a in ['graph_attr', 'node_attr', 'edge_attr']:
            variables[a + '_default'].update( variables[a + '_adj'])
            for k in attr:
                if variables[a + '_default'].has_key(k):
                    variables[a + '_default'][k] = attr[k]
                    print "%s k is updated in %s to be %s"%(k, a, attr[k])
            g.__dict__[a].update(variables[a + '_default'])

        #g.graph_attr.update(graph_attr_default)
        #g.graph_attr.update(graph_attr_adj)

        #g.node_attr.update(node_attr_adj)
        #g.edge_attr.update(edge_attr_adj)
        
        #g.layout('dot')   #neato, dot, twopi, circo, fdp, nop
        g.layout()
        path = path if path is not None else 'test.png'
        print g.string()
        g.draw(path)
        
        os.popen("display " + path)


def compare_template(g1, g2):
    g_1m2 = set(g1.keys())-set(g2.keys())
    print "only in g1",  g_1m2
    for i in g_1m2:
        print i, g1[i]
    
    g_2m1 = set(g2.keys())-set(g1.keys())
    print "\nonly in g2", g_2m1
    for i in g_2m1:
        print i,  g2[i]
    
    print "\nin both but different"
    g12 = set(g1.keys())&set(g2.keys())
    for i in g12:
        s1 = set(g1[i])
        s2 = set(g2[i])
        eq = s1 == s2 
        #print i, eq
        if not eq:
            print i
            print  "\t", (s1-s2),  s2-s1


class graphs_info:
    def __init__(self):
        self.cur_gn = 1
        self.cur_grn= -1
        self.Tot_GN = 0
        self.Tot_GN_Refl = 0
        self.D={}
        
    def Add(self, gid):
        if not self.D.has_key(gid):
            self.D[gid] = self.cur_gn
            self.cur_gn +=1
            self.Tot_GN += 1
            self.Tot_GN_Refl +=1
            
        grid = self.reflect_g(gid)
        if not self.D.has_key(grid):
            self.D[grid] = self.cur_grn
            self.cur_grn -= 1            
            self.Tot_GN += 1

    def Get_GN(self, gid):
        gn = self.D[gid]
        if gn<0: gn = self.Tot_GN+gn+1
        return gn
    
    def Show(self, name):
        ks = self.D.keys()
        ks.sort()
        print "   TotGN_%s=%d" % (name, self.Tot_GN)
        print "   RefGN_%s=%d" % (name, self.Tot_GN_Refl)
        
        for k in ks:
            gn = self.Get_GN(k)
            grn = self.Get_GN(self.reflect_g(k))            
            print "   %s(:, %d)=(/%d,%d,  %d/)" % \
                  (name, gn, k[0],k[1], grn)
        print
 

def Simplify(g):
    ng = g.Simplify_I()
    ng = ng.Simplify_U()
    ng = ng.Simplify_W1()
    ng = ng.Simplify_W2()
    ng = ng.Simplify_V()
    ng = ng.Simplify_O()
    return ng

def main():
    Diagram.Simplify = Simplify

    #import HalfBB
    import template_halfbb as HalfBB
    G=Diagram(HalfBB.G_HBB)
    for i in xrange(7, 8):
        oo = eval(HalfBB.oo2_tmpl % {"1":i, "2":i+1})
        oo = Diagram(oo)
        
        gg = G+oo
        ng = gg.Simplify()
        ng.show()
        
        max_leg, max_comp, optimal_path, max_cost = ng.find_optimal_path()
        print max_leg, max_comp, max_cost
        print optimal_path

class Test_diagram():
    def __init__(self):
        #from merapy.diagrams.V31.template_ternary import  G3 as G, V_tmpl
        from merapy.diagrams.V21.template_binary import  G2 as G, V_tmpl
        if 0:
            i = 1
            j=i*3+1
            D={"V":i, "I":j, "U1":i, "U2":i+1, "O":i}
            V = eval("{" + V_tmpl%D  + "}")

            G = Diagram(V)
        if 1:
            G = Diagram(G)
        self.graph = G
        #print G.edges
        #G.show_edges()
        #G.show_edge_info()
        #e1=G.get_edges(node="V_0")
        #print e1
    def to_networkx(self):
        self.graph.to_networkx()



if __name__=="__main__":
    td = Test_diagram()

    #td.to_networkx()
    td.graph.plot()
    #print td.graph.connections
    #print td.graph.edges
    
    if 0:
        Diagram.Simplify = Simplify

        #import HalfBB
        import Ternary

        #G=Diagram(HalfBB.G_HBB)
        G = Diagram(Ternary.G3)
        #print G.edges
        #G.show_edges()
        #G.show_edge_info()
        e1=G.get_edges(node="V_0")
        print e1




