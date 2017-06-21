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

from poly_term import PolyTerm
from bonds import bond, dbond
from concat_term import FFconcat
from numpy import *
from ewsum import lat_pts

# Adds the "pair" forcefield term into the list of atomic interactions.
# If present, L should be a 3x3 lower-diagonal array of translation row-vectors.
class LJPair(PolyTerm):
    def __init__(self, name, edges, L=None, R2=11.0**2):
        # internal vars
        PolyTerm.__init__(self, "ljpair_" + name, 2)
	#self.ineqs = [array([0.,  0.0, 1.0]), 0.0,
	#	      array([0.,  0.0,-1.0]), 100.0,
	#	      array([0., -1.0, 0.0]), 0.0,
	#	      array([0.,  1.0, 0.0]), -100.0]
	self.ineqs = [array([0.,  0.0, 1.0]),
		      array([0., -1.0, 0.0]),
		      array([0.,  1.0, 1.0])]
        self.edges = edges
	self.R2 = R2
        if L != None:
            self.periodic = True
            self.box_L = L
        else:
            self.periodic = False

    def deltas(self, x): # get pair distance array
        if self.periodic: # Wrap for periodicity.
	    delta = []
	    xv = reshape(x, (-1,x.shape[-2],3))
	    for i,j in self.edges:
		for p in x[:,j] - x[:,i]:
		    for z in lat_pts(delta, self.R2, self.box_L):
			delta.append(-z[1])
	    return array(delta)
	else:
	    return array([x[...,j,:]-x[...,i,:] for i,j in self.edges])

# Returns the energy of configuration x and/or the conserved force on each atom.
    def energy(self, c, x):
	delta = self.deltas(x)
        b = bond(delta)**-6
        return sum(info.f.y(c, b), -1)
	
    def force(self, c, x):
        delta = self.deltas(x)
        r, db = dbond(delta) # r, dr/dx
        u, du = self.f.y(c, r**-6, 1) # f(r^-6), df/du (u = r^-6)
        # du/dr = -6 r^-7
        du *= -6*r**-7

        en = sum(u, -1)
        de = zeros(x.shape)
        for k,(i,j) in enumerate(self.edges):
            de[...,i,:] += du[k]*db[...,k,:] # df/du du/dr dr/dx
            de[...,j,:] -= du[k]*db[...,k,:]
        return en, de

    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design(self, x, order=0):
        delta = self.deltas(x)

        if order == 0:
            return sum(self.spline(bond(delta)**-6, order), -2)
        elif order == 1:
            Ad = zeros(x.shape + (self.params,))
            r, db = dbond(delta) # r, dr/dx
            spl, dspl = self.spline(r**-6, order) # D(r^-6), dD(r^-6)/du
            dspl *= -6*r[...,newaxis]**-7 # du/dr
            for k,(i,j) in enumerate(self.edges): # dD/dx
                Ad[...,i,:,:] -= db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
                Ad[...,j,:,:] += db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
            A = sum(spl, -2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                  "supported."

