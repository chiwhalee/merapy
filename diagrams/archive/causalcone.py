#!/usr/bin/env python

m = {
    -3:0, -2:0, -1:(0,1), 0:(0,1), 
     1:(0,1), 2:(0,1), 3:1, 4:1,
     (-4,-3):(-1,0), (-3,-2):0, (-2,-1):(0,1), (-1,0):(0,1),
     (0,1):(0,1), (1,2):(0,1), (2,3):(0,1), (3,4):1, (4,5):(1,2)     
     }

m = {
    -3:(0,1), -2:(0,1), -1:(0,1), 0:(0,1), 
     1:(0,1), 2:(0,1), 3:(0,1), 4:(0,1),
     (-4,-3):(-1,0), (-3,-2):(0,1), (-2,-1):(0,1), (-1,0):(0,1),
     (0,1):(0,1), (1,2):(0,1), (2,3):(0,1), (3,4):(0,1), (4,5):(1,2)    
    }
# [8*l-3, 8*l]->[2*l,2*l+1]

import types
def isTuple(l):
    return type(l) == types.TupleType

def Asc(i):
    if isTuple(i):
        l,j = divmod(i[0]+3, 8)
        j = (j-3,j-2)
        if isTuple(m[j]):
            return (2*l+m[j][0], 2*l+m[j][1])
        else:
            return 2*l+m[j]
#        return (i[0]/4, i[1]/4+1)
    else:
        l,j = divmod(i+3,8)
        j = j-3
        if isTuple(m[j]):
            return (2*l+m[j][0], 2*l+m[j][1])
        else:
            return 2*l+m[j]
    
x=64
for i in xrange(-x+1,x+1):
    i1 = Asc(i); i2 = Asc(i1); i3 = Asc(i2)
    l1,j = divmod(i+3,4); l2,j = divmod(l1+3,4)
    l3,j = divmod(l2+3,4)
    print i, '|', l1,l2,l3, '|', i1, i2, i3
#     l1,j = divmod(l+3,4)
#     l2,j = divmod(l1+3,4)
#     print "%5d  %5d  %5d  %5d" % (i,l,l1,l2)
    
#     l,j = divmod(i+3,8)
#     j=j-3
#     print i, l, j, M[0][j]    
