from numpy import *
from .spline_term import SplineTerm
from .bspline import Bspline
from .concat_term import FFconcat

# Single bond term type shared by all edges in the list
# e.g. C-C, or C-H
class SplineBond(SplineTerm):
    def __init__(self, name, edges):
        # internal vars
        SplineTerm.__init__(self, name, Bspline(4), 40, 0.0, 4.0, 2)
        self.edges = edges

    def energy(self, c, x):
        self.f.c = c
        A = bond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
        return sum(self.f.y(A), -1)

    def force(self, c, x):
        self.f.c = c
        b, db = dbond(array([x[...,j,:]-x[...,i,:] for i,j in self.edges]))
        u, du = self.f.y(b, 1)

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
        raise NotImplementedError("Error! >1 energy derivative not "\
                              "supported.")

def bond(x):
	if len(x.shape) > 2:
		trns = list(range(len(x.shape)))
		trns = trns[1:-1]+trns[:1]+trns[-1:]
		x = transpose(x, trns)
	return sqrt(sum(x*x,-1))

def dbond(x):
	if len(x.shape) > 2:
		trns = list(range(len(x.shape)))
		trns = trns[1:-1]+trns[:1]+trns[-1:]
		x = transpose(x, trns)
	# x=j-i => dr(j-i)/dj = x/r
	r = sqrt(sum(x*x,-1))
	x /= r[...,newaxis]
	return r, x

# Compute force design matrix
# Array Float (blk, N, 3) -> Array Float (blk, N, 3)

# Compute force
# Array Float M -> Array Float (N, 3) -> Array Float (N, 3)

# Create an FFTerm from all bonds in PDB
# PDB -> (name : String -> edges : Set (Int,Int) -> FFTerm) -> FFTerm
def bond_terms(pdb, mkterm, edge=None):
    if edge is None:
        edge = pdb.edge
    bond_index = {}
    for i,j in edge:
        ti = pdb.names[i][2]
        tj = pdb.names[j][2]
        if ti > tj: # "Strike that; reverse it." -- Willy Wonka.
            ti, tj = (tj, ti)
            i, j = (j, i)

        name = "%s_%s"%(ti, tj)
        if not bond_index.has_key(name):
            bond_index[name] = []
        bond_index[name].append((i,j))

    print("%d Bonds"%sum(map(len, bond_index.values())))
    terms = [mkterm(n,l) for n,l in bond_index.iteritems()]
    return FFconcat(terms)

