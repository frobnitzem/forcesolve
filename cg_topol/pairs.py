# 2-body pairwise interactions

# This file is part of ForceSolve, Copyright (C) 2008 David M. Rogers.
#
#   ForceSolve is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ForceSolve is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ForceSolve (i.e. frc_solve/COPYING).
#   If not, contact the author(s) immediately \at/ wantye \/ gmail.com or
#   http://forceSolve.sourceforge.net/. And see http://www.gnu.org/licenses/
#   for a copy of the GNU GPL.


# David M. Rogers
# Dr. Thomas Beck Lab
# University of Cincinnati
# This work was supported by a DOE CSGF.

from spline_term import SplineTerm
from edge import modprod
from bonds import bond, dbond
from concat_term import FFconcat
from numpy import *

# Adds the "pair" forcefield term into the list of atomic interactions.
class SplinePair(SplineTerm):
    def __init__(self, name, edges, L=None):
        # internal vars
        SplineTerm.__init__(self, "pair_" + name, Bspline(6), 170, 0.0, 17.0, 2)
        self.edges = edges
        if L != None:
            self.periodic = True
            self.box_L = L
        else:
            self.periodic = False

# Returns the energy of configuration x and/or the conserved force on each atom.
    def energy(self, c, x):
        self.f.c = c
        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)
        A = bond(delta)
        return sum(info.f.y(A), -1)
	
    def force(self, c, x):
        self.f.c = c

        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)
        b, db = dbond(delta)
        u, du = self.f.y(b, 1)

        en = sum(u, -1)
        F = zeros(x.shape)
        for k,(i,j) in enumerate(plist):
            F[...,i,:] += du[k]*db[...,k,:]
            F[...,j,:] -= du[k]*db[...,k,:]
        return en,F

    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design_pair(self, x, order=0):
        A = []
        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)

        if order == 0:
            return sum(self.spline(bond(delta), order),-2)
        elif order == 1:
            Ad = zeros(x.shape + (info.f.n,))
            b,db = dbond(delta)
            spl, dspl = info.spline(b, order)
            for k,(i,j) in enumerate(plist):
                Ad[...,i,:,:] -= db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
                Ad[...,j,:,:] += db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
            A = sum(spl,-2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                  "supported."

# set join
mconcat = lambda m: reduce(lambda x,y: x|y, m, set())
# extend neighbors
extend = lambda pdb, x: x | mconcat(pdb.conn[b] for b in x)

def orderset(a, x):
    s = set()
    for b in x:
	if a > b:
	    s.add((b,a))
	else:
	    s.add((a,b))
    return s

# n = 4 => count 1,4 pairs
def pair_terms(pdb, mkterm, n=4):
        assert n >= 2, "Can't count self-pairs."
        xpair = [set([a]) for a in range(pdb.atoms)]
	for i in range(n-2): # Extend table by 1 bond.
	    for a in range(pdb.atoms):
		xpair[a] = extend(pdb, xpair[a])
	xpair = mconcat([orderset(a,x) for a,x in enumerate(xpair)])
	pair = []
	for i in range(pdb.atoms-1):
	    pair += [(i,j) for j in range(i+1,pdb.atoms)]
	pdb.pair = set(pair)-xpair

	pair_index = {}
	for i,j in pdb.pair:
		ti = pdb.names[i][2]
		tj = pdb.names[j][2]
		if ti > tj:
                    ti, tj = (tj, ti)
                    i, j = (j, i)
		name = "%d+%s_%s"%(n,ti,tj)
                if not pair_index.has_key(name):
                    pair_index[name] = []
                pair_index[name].append((i,j))
        terms = [mkterm(name, l, pdb.L) \
                        for name, l in pair_index.iteritems()]
        print "%d types of >= 1,%d pairs."%(len(terms), n)
        return FFconcat(terms)

# count only 1,n pairs
# still add to pdb.pair
def pair_n_terms(pdb, mkterm, n=4, pair_n=None):
        assert n >= 2, "Need at least 2 atoms to make a pair!"
        if pair_n is None:
            xpair = [set([a]) for a in range(pdb.atoms)]

            for i in range(n-2): # Extend table by 1 bond.
                for a in range(pdb.atoms):
                    xpair[a] = extend(pdb, xpair[a])

            pair_n = [extend(pdb, x) - x for x in xpair]
            pair_n = mconcat([orderset(a, x) for a,x in enumerate(pair_n)])

        pdb.pair |= set(pair_n)

	pair_index = {}
	for i,j in pair_n:
		ti = pdb.names[i][2]
		tj = pdb.names[j][2]
		if ti > tj:
                    ti, tj = (tj, ti)
                    i, j = (j, i)
		name = "1,%d_%s_%s"%(n,ti,tj)
                if not pair_index.has_key(name):
                    pair_index[name] = []
                pair_index[name].append((i,j))
        terms = [mkterm(name, l, pdb.L) \
                        for name, l in pair_index.iteritems()]
        print "%d types of special 1,%d pairs"%(len(terms),n)
        return FFconcat(terms)

