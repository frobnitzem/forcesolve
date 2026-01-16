# 4-body torsional interactions

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

from numpy import *
import numpy.linalg as la
from .edge import modprod
from .poly_term import PolyTerm
from .bspline import Bspline
from .concat_term import FFconcat
from .torsions import cross_product

# [cos(n phi)] = V . [(cos phi)^n]
#V = array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.],
#           [ 0.,  1.,  0.,  0.,  0.,  0.,  0.],
#           [-1.,  0.,  2.,  0.,  0.,  0.,  0.],
#           [ 0., -3.,  0.,  4.,  0.,  0.,  0.],
#           [ 1.,  0., -8.,  0.,  8.,  0.,  0.],
#          [ 0.,  5.,  0.,-20.,  0., 16.,  0.],
#          [-1.,  0., 18.,  0.,-48.,  0., 32.]])
W = array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.,  0.,  0.,  0.],
           [ 2.,  0.,  1.,  0.,  0.,  0.,  0.],
           [ 0.,  3.,  0.,  1.,  0.,  0.,  0.],
           [ 6.,  0.,  4.,  0.,  1.,  0.,  0.],
           [ 0., 10.,  0.,  5.,  0.,  1.,  0.],
           [20.,  0., 15.,  0.,  6.,  0.,  1.]])
W[:,1:] *= 2.0
W /= (2**arange(len(W)))[:,newaxis]

# Generate constraints for all but cos(n phi + d)
# since d = 0 or 180, it only changes the sign, so
# the constraints are unaffected.
def gen_constrain_n(n,m):
    nparam = len(n)
    if nparam == 0:
        return identity(m)
    un = [i for i in range(m) if i not in n]
    return transpose(W)[un] # gets a bunch of rows

# Adds the "tor" forcefield term into the list of atomic interactions.
class PolyTorsion(PolyTerm):
    def __init__(self, name, tors, constrain_n=False):
        PolyTerm.__init__(self, name, 6)
        if constrain_n != False and constrain_n(name) is not None:
            self.constraints = gen_constrain_n(constrain_n(name),7)
        self.tors = tors
        
    def energy(self, c, x):
        A = torsionc(array([[x[...,i,:]-x[...,j,:],\
                             x[...,k,:]-x[...,j,:],\
                             x[...,l,:]-x[...,k,:]] \
                                 for i,j,k,l in self.tors]))
        return sum(self.f.y(c, A), -1)
    
    def force(self, c, x):
        t, dt = dtorsionc(array([[x[...,i,:]-x[...,j,:],\
                                  x[...,k,:]-x[...,j,:],\
                                  x[...,l,:]-x[...,k,:]] \
                                  for i,j,k,l in self.tors]))
        u, du = self.f.y(c, t, 1)
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
            tor = torsionc(array([[x[...,i,:]-x[...,j,:],\
                                   x[...,k,:]-x[...,j,:],\
                                   x[...,l,:]-x[...,k,:]] \
                                     for i,j,k,l in self.tors]))
            return sum(self.spline(tor, order),-2)
        elif order == 1:
            Ad = zeros(x.shape+(self.params,))
            t, dt = dtorsionc(array([[x[...,i,:]-x[...,j,:],\
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
            raise NotImplementedError("Error! >1 energy derivative not "\
                                  "supported.")

# cosine of torsion
def torsionc(x):
    trsp = list(range(len(x.shape)))
    # Operate on last 2 dim.s (atom and xyz) + (..., number)
    trsp = trsp[1:2] + trsp[-1:] + trsp[2:-1] + trsp[:1]
    x = transpose(x, trsp)

    A0 = cross_product(x[0],x[1])
    A1 = cross_product(x[2],x[1])
    return sum(A0*A1,0)/sqrt(sum(A0*A0,0)*sum(A1*A1,0))

# derivative of cosine of torsion
def dtorsionc(x):
    trsp = list(range(len(x.shape))) # (tor serial, tor atom, ..., xyz)
    trsp = trsp[1:2] + trsp[-1:] + trsp[2:-1] + trsp[:1]
    x = transpose(x, trsp) # (tor atom, xyz, ..., tor serial)

    d = zeros(x.shape)
    
    A0 = cross_product(x[0],x[1]) # -(a ^ b)^*
    A1 = cross_product(x[2],x[1]) # -(c ^ b)^*
    x0 = sqrt(sum(A0*A0,0))
    x1 = sqrt(sum(A1*A1,0))
    A0 /= x0
    A1 /= x1

    cphi = sum(A0*A1, 0)
    # a x (b x c) = - a . (b ^ c) = b a.c - c a.b
    d[0] = cross_product(x[1], A1 - cphi*A0)/x0
    d[2] = cross_product(x[1], A0 - cphi*A1)/x1
    d[1] = -cross_product(x[0], A1 - cphi*A0)/x0 \
           -cross_product(x[2], A0 - cphi*A1)/x1

    trsp = list(range(2, len(x.shape))) + [0,1] # Move (atom,xyz) to last 2 dim.s
    return cphi, transpose(d, trsp)

