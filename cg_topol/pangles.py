# 3-body angle interaction code.

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

# David M. Rogers
# Dr. Thomas Beck Lab
# University of Cincinnati
# This work was supported by a DOE CSGF.

from spline_term import SplineTerm
from bspline import Bspline
from edge import modprod
from concat_term import FFconcat
from numpy import *

# Adds the "angle" forcefield term into the list of atomic interactions.
class PolyAngle(PolyTerm):
    def __init__(self, name, angs):
        PolyTerm.__init__(self, "pangle_" + name, 2)
        self.angs = angs

    def energy(self, c, x):
        A = angle(array([[x[...,i,:]-x[...,j,:],\
                          x[...,k,:]-x[...,j,:]] for i,j,k in self.angs]))
        return sum(self.f.y(c, A), -1)
    
    def force(self, c, x):
        a, da = dangle(array([[x[...,i,:]-x[...,j,:],\
                               x[...,k,:]-x[...,j,:]] for i,j,k in self.angs]))
        u, du = self.f.y(c, a, 1)

        en = sum(u, -1)
        F = zeros(x.shape)
        for z,(i,j,k) in enumerate(self.angs):
            F[...,i,:] -= du[z]*da[...,z,0,:]
            F[...,j,:] += du[z]*sum(da[...,z,:,:],-2)
            F[...,k,:] -= du[z]*da[...,z,1,:]
        return en, F

    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    def design(self, x, order=0):
        if order == 0:
            a = angle(array([[x[...,i,:]-x[...,j,:],\
                              x[...,k,:]-x[...,j,:]] for i,j,k in self.angs]))
            spl = self.spline(a, order)
            A = sum(spl,-2) # Sum over all atoms in ea. structure.
            return A
        elif order == 1:
            Ad = []
            Ad = zeros(x.shape+(self.f.n,), float)
            a,da = dangle(array([ [x[...,i,:]-x[...,j,:],\
                                   x[...,k,:]-x[...,j,:]] \
                                 for i,j,k in self.angs]))
            spl, dspl = self.spline(a, order)
            for z,(i,j,k) in enumerate(self.angs):
                Ad[...,i,:,:] += da[...,z,0,:,newaxis] \
                                * dspl[...,z,newaxis,:]
                Ad[...,j,:,:] -= sum(da[...,z,:,:],-2)[...,newaxis] \
                                * dspl[...,z,newaxis,:]
                Ad[...,k,:,:] += da[...,z,1,:,newaxis] \
                                * dspl[...,z,newaxis,:]
            A = sum(spl,-2)
            return A, Ad
        else:
            raise RuntimeError, "Error! >1 energy derivative not "\
                                "supported."
def angle(x):
    if len(x.shape) > 3:
            trns = range(len(x.shape))
            trns = trns[2:-1]+trns[:2]+trns[-1:]
            x = transpose(x, trns)
    cth = sum(x[...,0,:]*x[...,1,:],-1) \
      / sqrt(sum(x[...,0,:]*x[...,0,:],-1)*sum(x[...,1,:]*x[...,1,:],-1))
    return arccos(cth)

def dangle(x):
    if len(x.shape) > 3:
        trns = range(len(x.shape))
        trns = trns[2:-1]+trns[:2]+trns[-1:]
        x = transpose(x, trns)
    rijk_inv = 1.0/sqrt(sum(x*x,-1)) # x.shape by 2
    x *= rijk_inv[..., newaxis] # Normalize x.
    cth = sum(x[...,0,:]*x[...,1,:],-1) # x.shape vector of cosines.
    A = cth[...,newaxis,newaxis]*x # Future derivatives.
    A[...,0,:] = x[...,1,:]-A[...,0,:]
    A[...,1,:] = x[...,0,:]-A[...,1,:]
    A *= rijk_inv[..., newaxis] # d(cth) / di and dk
    return arccos(cth), -A/sqrt(1.0-cth*cth)

