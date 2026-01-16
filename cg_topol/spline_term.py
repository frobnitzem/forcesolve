from numpy import reshape
import numpy.random as rand
from .bspline import spline_func

# Default parameter count, constraints, prior, and pri_sample for 
class SplineTerm:
    def __init__(self, name, bspl, n, x0, L, m):
        order = 2 # derivative of spline for smoothness
        self.hyp_params = 1 # hyper-parameters

        if x0 == None: # periodic
            self.params = n
            self.f = spline_func((-0.5*L, L/n, self.params), bspl, \
                                   True, bspl.order*0.5, (-0.5*L,0.5*L), name)
        else:
            self.params = n+bspl.order-1
            self.f = spline_func((x0, L/n, self.params), bspl, \
                                   False, bspl.order-1, (x0,x0+L), name)
        self.name = name
        self.constraints = [self.f.integral(0, m)[0]]
        self.ineqs = []

        scale = (m + 1.0)/abs(self.f.rng[1]**(m + 1.0) \
                              - self.f.rng[0]**(m + 1.0))
        self.prior = [ (0, self.f.quadratic_integral(order, m)*scale) ]
        self.pri_rank = [self.params - order]
        self.ind = [0]

    # Returns "order+1" arrays giving up to the order-th derivative of f(x).
    # X can be any shape, thus each return array is x.shape+(n,)
    def spline(self, x, order=0):
        if order == 0:
            spl = self.f.spline(reshape(x,x.size), None)
            return reshape(spl, x.shape+(spl.shape[-1],))
        else:
            spl = self.f.spline(reshape(x,x.size), order)
            return reshape( spl, \
                    (spl.shape[0],) + x.shape + (spl.shape[-1],) )

    def write(self, pre, c):
        self.f.c = c
        self.f.write_spl(pre + self.name+".espl")
    
    #def read(self, base, name):
    #    self.name = name
    #    self.f = read_spline_func(base+"%s.espl"%name, name)
