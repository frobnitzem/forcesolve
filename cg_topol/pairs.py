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

from bspline import Bspline
from spline_term import SplineTerm
from edge import modprod, srt2
from bonds import bond, dbond
from concat_term import FFconcat
from numpy import *

# Generate the actual pairs to consider.
def ex_gen(edges, excl):
    if excl is None:
        for i,j in edges:
            yield i,j
    else:
        if excl[0][0] == excl[1][0]: # same list
           for ii,i in enumerate(excl[0][:-1]):
               for j in excl[1][ii+1:]:
                   if (i,j) in edges:
                       continue
                   yield i,j
        else:
            for i in excl[0]:
                for j in excl[1]:
                    if srt2(i,j) in edges:
                        continue
                    yield i,j

def pbc_delta(x, edges, excl, L):
    assert False, "This routine needs to be fixed (output of >1 delta-r per pair won't currently work)."
    delta = []
    xv = reshape(x, (-1,x.shape[-2],3))
    N = x.shape[-2]
    for i,j in ex_gen(edges, excl):
	for p in x[:,j] - x[:,i]:
	    for z in lat_pts(p, self.R2, self.box_L):
		delta.append(-z[1])
    return array(delta)

# TODO:  This should be implemented to return a stream of chunks (generator).
#        That way, it could eventually be coupled to a consumer
#        without using wasting memory in the midterm.
#      However, memory use would only be O(N) for pbc
#      if a distance-test were done here.
def calc_delta(x, edges, excl, L):
    if L is not None:
        return pbc_delta(x, edges, L)
    if excl is None:
        N = len(edges)
    elif excl[0][0] == excl[1][0]: # same list
        N = len(excl[0])*(len(excl[0])-1)/2 - len(edges)
    else:
        N = len(excl[0])*len(excl[1]) - len(edges)
    delta = zeros((N,) + x.shape[:-2] + (3,))
    k = -1
    for k,(i,j) in enumerate(ex_gen(edges, excl)):
        delta[k] = x[...,j,:]-x[...,i,:]
    assert N == k+1, "Incorrectly formed Exclusion term!"
    return delta

# Adds the "pair" forcefield term into the list of atomic interactions.
# The default (excl == None) is to treat 'edges' as pairs to *include*.
# If excl == ([i], [j]), then 'edges' is the *exclude* list instead,
# and all pair interactions between atoms in list [i] and [j] are
# computed unless they are in 'edges'.
# If present, L should be a 3x3 lower-diagonal array of translation row-vectors.
class SplinePair(SplineTerm):
    def __init__(self, name, edges, L=None, excl=None):
        # internal vars
        SplineTerm.__init__(self, name, Bspline(6), 100, 0.0, 11.0, 2)
        self.edges = set( srt2(*e) for e in edges if e[0] != e[1] )
        self.excl = excl
        self.L = L

# Returns the energy of configuration x and/or the conserved force on each atom.
    def energy(self, c, x):
        self.f.c = c
        delta = calc_delta(x, self.edges, self.excl, self.L)
        A = bond(delta)
        return sum(self.f.y(A), -1)
	
    def force(self, c, x):
        self.f.c = c
        delta = calc_delta(x, self.edges, self.excl, self.L)
        b, db = dbond(delta)
        u, du = self.f.y(b, 1)

        en = sum(u, -1)
        F = zeros(x.shape)
        for k,(i,j) in enumerate(ex_gen(self.edges, self.excl)):
            delta[k] = x[...,j,:]-x[...,i,:]
            F[...,i,:] += du[k]*db[...,k,:]
            F[...,j,:] -= du[k]*db[...,k,:]
        return en,F

    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design_pair(self, x, order=0):
        A = []
        delta = calc_delta(x, self.edges, self.excl, self.L)

        if order == 0:
            return sum(self.spline(bond(delta), order),-2)
        elif order == 1:
            Ad = zeros(x.shape + (self.f.n,))
            b,db = dbond(delta)
            spl, dspl = self.spline(b, order)
            for k,(i,j) in enumerate(ex_gen(self.edges, self.excl)):
                Ad[...,i,:,:] -= db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
                Ad[...,j,:,:] += db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
            A = sum(spl,-2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                  "supported."

