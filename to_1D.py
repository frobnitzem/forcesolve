from numpy import *
from bspline import spline_func, Bspline
from pspline import pspline_match

def kpt(M):
    k = arange(M)
    k[(M+1)/2:] -= M
    return k

# Treat M as a bunch of samples of a distance-dependent function
# on a 3D lattice.
def to_1D(M, n=None, order=4, deriv=2):
    if n == None:
        n = array(M.shape).max()

    r = (kpt(M.shape[0])**2)[:,newaxis,newaxis] \
      + (kpt(M.shape[1])**2)[newaxis,:,newaxis] \
      + (kpt(M.shape[2])**2)[newaxis,newaxis,:]
    r = sqrt(r)

    x0 = 0.0
    L = r.max()
    h = L/n

    n += order-1 # Add external points.
    shift = order-1 # and shift left

    f = spline_func((x0,h,n),Bspline(order),False,shift,(x0,x0+L), "grid")
    
    m = pspline_match(f, deriv, 2, 1.0e-8)
    m.append(reshape(r, -1), reshape(M, -1))
    m.alpha = 0.0 # turn off smoothing
    m.maximize(maxiter=1)
    #m.dimensionality()
    #m.sample(samples,skip,toss)
    f.c = m.theta
    return m
