#!/usr/bin/env python2.7

from cg_topol import *
from numpy import *
from num_deriv import num_deriv

# manually create a term
#term = PolyBond("test", [(0,1)])
#term = PolyAngle("test", [(0,1,2), (1,2,3)]) # try 2 angles at once
#term = PolyTorsion("test", [(0,1,2,3)])
term = LJPair("test", [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])

random.seed(1)
# create a config.
#x = zeros((4,3))
#x[:,2] = arange(4)*1.2
#x[1,0] = -0.2
x = random.random((4,3))*10
#print x

# grab the derivative (shape = x.shape + (params,))
D = term.design(x, 1)[1]
#print term.design(x)
#print D

# shape = x.shape + (params,)
D2 = num_deriv(lambda x: term.design(x)[0], x)

print abs(D2 - D) > 1e-5
print abs(D - D2).max()
