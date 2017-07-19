# 4-body improper interactions

# This file is part of ForceSolve, Copyright (C) 2008-2015 David M. Rogers.
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

from edge import modprod
from poly_term import PolyTerm
from bspline import Bspline
from concat_term import FFconcat
from numpy import *
import numpy.linalg as la
from ptorsions import torsionc, dtorsionc
from torsions import cross_product

# Adds the "pimprop" forcefield term into the list of atomic interactions.
class PolyImprop(PolyTerm):
    def __init__(self, name, tors):
        self.tors = tors

        self.hyp_params = 0
        self.params = 1
        self.name = "pimprop_" + name
        self.constraints = []
	self.ineqs = [array([1.0])]
        self.prior = []
        self.pri_rank = []
        self.ind = []
	
    def energy(self, c, x):
        A = atorsion(array([[x[...,i,:]-x[...,j,:],\
                             x[...,k,:]-x[...,j,:],\
                             x[...,l,:]-x[...,k,:]] \
                                 for i,j,k,l in self.tors]))
        return self.c[0]*sum(A*A, -1)
    
    def force(self, c, x):
        t, dt = datorsion(array([[x[...,i,:]-x[...,j,:],\
                                  x[...,k,:]-x[...,j,:],\
                                  x[...,l,:]-x[...,k,:]] \
                                  for i,j,k,l in self.tors]))
        u, du = self.c[0]*t*t, self.c[0]*2.0*t
        en = sum(u, -1)
        F = zeros(x.shape)
        for z,(i,j,k,l) in enumerate(alist):
            F[...,i,:] -= du[z]*dt[...,z,0,:]
            F[...,j,:] += du[z]*(dt[...,z,1,:] + dt[...,z,0,:])
            F[...,k,:] -= du[z]*(dt[...,z,1,:] - dt[...,z,2,:])
            F[...,l,:] -= du[z]*dt[...,z,2,:]
        return en, F
    
    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design(self, x, order=0):
        A = []
        if order == 0:
            tor = atorsion(array([[x[...,i,:]-x[...,j,:],\
                                   x[...,k,:]-x[...,j,:],\
                                   x[...,l,:]-x[...,k,:]] \
                                     for i,j,k,l in self.tors]))
            return sum(self.spline(tor, order),-2)
        elif order == 1:
            Ad = zeros(x.shape+(self.params,))
            t, dt = datorsion(array([[x[...,i,:]-x[...,j,:],\
                                      x[...,k,:]-x[...,j,:],\
                                      x[...,l,:]-x[...,k,:]] \
                                         for i,j,k,l in self.tors]))
            spl, dspl = self.spline(t, order)
            for z,(i,j,k,l) in enumerate(self.tors):
                Ad[...,i,:,:] += dt[...,z,0,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
                Ad[...,j,:,:] -= (dt[...,z,1,:,newaxis]+dt[...,z,0,:,newaxis]) \
                                        * dspl[...,z,newaxis,:]
                Ad[...,k,:,:] += (dt[...,z,1,:,newaxis]-dt[...,z,2,:,newaxis]) \
                                        * dspl[...,z,newaxis,:]
                Ad[...,l,:,:] += dt[...,z,2,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
            A = sum(spl, -2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                  "supported."

    # Used to construct vectors which multiply parameters.
    # If x is a N-dim vector, the return value is an (nd+1)xNxP matrix
    def spline(self, x, nd=0):
	if nd == None:
	    return x[..., newaxis]**2
        Mp = x[newaxis,...]**(2 - arange(nd+1).reshape([nd+1]+[1]*len(x.shape)))
        # d^n/dx^n [x^a] = x^{a-n} prod_{i=0}^{n-1} a - i
        #                = x^{a-n} prod_{j=a+1-n}^a j
	if nd >= 1:
	    Mp[1:] *= 2.0
	    if nd >= 3:
		Mp[3:] = 0.0
        return Mp[...,newaxis] # P = 1

    def write(self, pre, c, mode='w'):
        name = pre + self.name + ".impr"
        out = open(name, mode)
        out.write("#IMPR %s %e\n"%(self.name,c[0]))
        out.close()

# arccos(cosine of torsion) -- is in [0,pi]
def atorsion(x):
    return arccos(torsionc(x))

# derivative of cosine of torsion
def datorsion(x):
    x, dx = dtorsionc(x)
    return arccos(x), -dx/sqrt(1.0-x*x)[...,newaxis,newaxis]

# List out all improper terms (atoms with 3 bonds)
def improper_terms(pdb, mkterm):
    tors = set()
    for i in range(pdb.atoms):
	if len(pdb.conn[i]) == 3:
	    tors.add((i,) + tuple(pdb.conn[i]))

    imp_index = {}
    for i,j,k,l in tors:
        ti = pdb.names[i][2]
        tj = pdb.names[j][2]
        tk = pdb.names[k][2]
        tl = pdb.names[l][2]
        if tk > tl: # Bubble sort
	    tk, tl = tl, tk
	    k, l   = l, k
	if tj > tk:
	    tj, tk = tk, tj
	    j, k   = k, j
	    if tk > tl: # Bubble sort
		tk, tl = tl, tk
		k,   l = l, k
	if tk == tl: # according to the CHARMM docs, we swap again!
            tj, tk, tl = tk, tl, tj
            j, k, l    = k, l, j
        name = "%s-%s-%s-%s"%(ti,tj,tk,tl)
        if not imp_index.has_key(name):
            imp_index[name] = []
        imp_index[name].append((i,j,k,l))

    #pdb.impr = imp_index

    print "%d Impropers"%sum(map(len, imp_index.values()))
    terms = [mkterm(n,l) for n,l in imp_index.iteritems()]
    return FFconcat(terms)

