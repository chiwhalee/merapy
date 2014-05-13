#!/usr/bin/env python

def reflect(i):
    return 4-i

def reflect_g(g):
    return (reflect(g[1]), reflect(g[0]))

D={}
gn=1; grn=9
for j in xrange(3, 0, -1):
    for i in xrange(1, 4):
        g = (i,j)
        gr = reflect_g(g)
        if D.has_key(g) or D.has_key(gr):
            continue
        if (gr == g):
            D[g] = gn
            gn = gn+1
        else:
            D[g] = gn
            D[gr] = grn
            gn = gn+1
            grn = grn-1


k = D.keys()
k.sort()
for g in k:
    p = D[g]
    gr = reflect_g(g)
    pr = D[gr]
    print "    GP_2_3(:,%d)=(/%d,%d, %d/)" % (p, g[0], g[1], pr)

