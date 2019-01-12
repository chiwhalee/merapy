#!/usr/bin/env python

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import *
def reflect(i):
    return 4-i

def reflect_g(g):
    return (reflect(g[1]), reflect(g[0]))

D={}
gn=1; grn=9
for j in range(3, 0, -1):
    for i in range(1, 4):
        g = (i,j)
        gr = reflect_g(g)
        if g in D or gr in D:
            continue
        if (gr == g):
            D[g] = gn
            gn = gn+1
        else:
            D[g] = gn
            D[gr] = grn
            gn = gn+1
            grn = grn-1


k = list(D.keys())
k.sort()
for g in k:
    p = D[g]
    gr = reflect_g(g)
    pr = D[gr]
    print("    GP_2_3(:,%d)=(/%d,%d, %d/)" % (p, g[0], g[1], pr))

