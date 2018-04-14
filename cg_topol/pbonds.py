from poly_term import PolyTerm
from bspline import Bspline
from concat_term import FFconcat
from numpy import *
from bonds import bond, dbond

# Single bond term type shared by all edges in the list
# e.g. C-C, or C-H
class PolyBond(PolyTerm):
    def __init__(self, name, edges):
        # internal vars
        PolyTerm.__init__(self, name, 2)
	self.ineqs = [array([0.,  0.0, 1.0]),
		      array([0., -1.0, 0.0])]
        self.edges = edges

    def energy(self, c, x):
        A = bond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
        return sum(self.f.y(c, A), -1)

    def force(self, c, x):
        b, db = dbond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
        u, du = self.f.y(c, b, 1)

        en = sum(u, -1)
        F = zeros(x.shape)
        for k, (i,j) in enumerate(self.edges):
            F[...,i,:] += du[k]*db[...,k,:]
            F[...,j,:] -= du[k]*db[...,k,:]

        return en, F

    def design(self, x, order=0):
        if order == 0:
            b = bond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
            spl = self.spline(b, order)
            return sum(spl, -2)
        elif order == 1:
            Ad = zeros(x.shape + (self.params,))
            b, db = dbond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
            spl, dspl = self.spline(b, order)
            for k, (i,j) in enumerate(self.edges):
                Ad[...,i,:,:] -= db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
                Ad[...,j,:,:] += db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
            A = sum(spl, -2)
            return A, Ad
        raise RuntimeError, "Error! >1 energy derivative not "\
                              "supported."

# UB is the same as a bond, but requires a different naming scheme.
class PolyUB(PolyBond):
    def __init__(self, name, angles):
        PolyTerm.__init__(self, name, 2)
        self.edges = set([(i, k) for i,j,k in angles])
	self.ineqs = [array([0.,  0.0, 1.0]),
		      array([0., -1.0, 0.0])]

