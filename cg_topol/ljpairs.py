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

# Adds the "pair" forcefield term into the list of atomic interactions.
class LJPair(PolyTerm):
    def __init__(self, name, edges, L=None):
        # internal vars
        PolyTerm.__init__(self, "ljpair_" + name, 2)
        self.edges = edges
        if L != None:
            self.periodic = True
            self.box_L = L
        else:
            self.periodic = False

# Returns the energy of configuration x and/or the conserved force on each atom.
    def energy(self, c, x):
        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)
        b = bond(delta)**-6
        return sum(info.f.y(c, b), -1)
	
    def force(self, c, x):
        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)
        r, db = dbond(delta) # r, dr/dx
        u, du = self.f.y(c, r**-6, 1) # f(r^-6), df/du (u = r^-6)
        # du/dr = -6 r^-7
        du *= -6*r**-7

        en = sum(u, -1)
        F = zeros(x.shape)
        for k,(i,j) in enumerate(plist):
            F[...,i,:] += du[k]*db[...,k,:] # df/du du/dr dr/dx
            F[...,j,:] -= du[k]*db[...,k,:]
        return en,F

    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design_pair(self, x, order=0):
        delta = array([x[...,j,:]-x[...,i,:] for i,j in self.edges])
        if self.periodic: # Wrap for periodicity.
            delta -= self.box_L * floor(delta/self.box_L+0.5)

        if order == 0:
            return sum(self.spline(bond(delta)**-6, order), -2)
        elif order == 1:
            Ad = zeros(x.shape + (info.f.n,))
            r, db = dbond(delta) # r, dr/dx
            spl, dspl = info.spline(r**-6, order) # D(r^-6), dD(r^-6)/du
            dspl *= -6*r**-7 # du/dr
            for k,(i,j) in enumerate(plist): # dD/dx
                Ad[...,i,:,:] -= db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
                Ad[...,j,:,:] += db[...,k,:,newaxis] \
                                * dspl[...,k,newaxis,:]
            A = sum(spl, -2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                  "supported."

