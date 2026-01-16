# 4-body torsional interactions

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

from numpy import array, zeros, transpose, arctan2
from .edge import modprod
from .spline_term import SplineTerm
from .bspline import Bspline
from .concat_term import FFconcat

# Adds the "tor" forcefield term into the list of atomic interactions.
class SplineTorsion(SplineTerm):
    def __init__(self, name, tors):
        SplineTerm.__init__(self, name, Bspline(4), 180, None, 2*pi, 0)
        self.tors = tors
        
    def energy(self, c, x):
        self.f.c = c
        A = torsion(array([[x[...,i,:]-x[...,j,:],\
                            x[...,j,:]-x[...,k,:],\
                            x[...,l,:]-x[...,k,:]] \
                                for i,j,k,l in self.tors]))
        return sum(self.f.y(A), -1)
    
    def force(self, c, x):
        self.f.c = c
        t, dt = dtorsion(array([[x[...,i,:]-x[...,j,:],\
                                 x[...,j,:]-x[...,k,:],\
                                 x[...,l,:]-x[...,k,:]] \
                                 for i,j,k,l in self.tors]))
        u, du = self.f.y(t, 1)
        en = sum(u, -1)
        F = zeros(x.shape)
        for z,(i,j,k,l) in enumerate(alist):
            F[...,i,:] -= du[z]*dt[...,z,0,:]
            F[...,l,:] -= du[z]*dt[...,z,1,:]
            F[...,j,:] += du[z]*dt[...,z,2,:]
            F[...,k,:] -= du[z]*dt[...,z,3,:]
        return en, F
    
    # The design array multiplies the spline coefficients to produce the
    # total bonded energy/force.
    # TODO: re-check j - k for central bond???
    def design_tor(self, x, order=0):
        A = []
        if order == 0:
            tor = torsion(array([[x[...,i,:]-x[...,j,:],\
                                  x[...,j,:]-x[...,k,:],\
                                  x[...,l,:]-x[...,k,:]] \
                                    for i,j,k,l in self.tors]))
            return sum(self.spline(tor, order),-2)
        elif order == 1:
            Ad = zeros(x.shape+(self.f.n,))
            t, dt = dtorsion(array([[x[...,i,:]-x[...,j,:],\
                                     x[...,j,:]-x[...,k,:],\
                                     x[...,l,:]-x[...,k,:]] \
                                        for i,j,k,l in alist]))
            spl, dspl = self.spline(t, order)
            for z,(i,j,k,l) in enumerate(self.tors):
                Ad[...,i,:,:] += dt[...,z,0,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
                Ad[...,l,:,:] += dt[...,z,1,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
                Ad[...,j,:,:] -= dt[...,z,2,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
                Ad[...,k,:,:] += dt[...,z,3,:,newaxis] \
                                        * dspl[...,z,newaxis,:]
            A = sum(spl, -2)
            return A, Ad
        else:
            raise NotImplementedError("Error! >1 energy derivative not "\
                                  "supported.")

def torsion(x):
        trsp = range(len(x.shape))
        # Operate on last 2 dim.s (atom and xyz) + (..., number)
        trsp = trsp[1:2] + trsp[-1:] + trsp[2:-1] + trsp[:1]
        x = transpose(x, trsp)
        
        x[1] = -x[1]/sqrt(sum(x[1]*x[1],0))
        x[0] = x[0] - x[1]*sum(x[0]*x[1],0)
        x[2] = x[2] - x[1]*sum(x[2]*x[1],0)
        proj = sum(x[2]*cross_product(x[1],x[0]),0)
        return arctan2(proj,sum(x[0]*x[2],0))
        
def dtorsion(x):
        trsp = range(len(x.shape)) # (tor serial, tor atom, ..., xyz)
        trsp = trsp[1:2] + trsp[-1:] + trsp[2:-1] + trsp[:1]
        x = transpose(x, trsp) # (tor atom, xyz, ..., tor serial)
        
        A = zeros((4,3) + x.shape[2:])
        
        A[1:4] = x
        x1 = sqrt(sum(x[1]*x[1],0))
        x2 = sum(x[0]*x[1],0)
        x3 = sum(x[2]*x[1],0)
        A[0] = cross_product(A[1],A[2])
        A[1] = cross_product(A[3],A[2])
        v = sum(A[0]*A[0],0)
        w = sum(A[1]*A[1],0)
        A[2] = A[1]*(x3/(w*x1))[newaxis,...]
        
        t = x2/(v*x1)
        x2 = sum(A[0]*A[1],0)
        x3 = sum(A[0]*A[3],0)*x1
        A[3] = A[0]*t[newaxis,...]
        A[2] -= A[3]
        A[0] *= (-x1/v)[newaxis,...]
        A[1] *= (x1/w)[newaxis,...]
        #t = sqrt(v*w)
        #x2 /= t # cos(phi)
        #x3 /= t # sin(phi)
        phi = -arctan2(x3,x2)
        A[3] = A[2]-A[1]
        A[2] += A[0]
        
        trsp = range(2, len(x.shape)) + [0,1] # Move (atom,xyz) to last 2 dim.s
        return phi, transpose(A, trsp)
        
def cross_product(x,y):
        return array( [ x[1]*y[2] - x[2]*y[1], \
                        x[2]*y[0] - x[0]*y[2], \
                        x[0]*y[1] - x[1]*y[0] ] )

def torsion_terms(pdb, mkterm, tors=None):
    if tors is None:
        tors = set()
        for j,k in pdb.edge:
            tors |= set([(i,j,k,l) for i,l in modprod(pdb.conn[j]-set([k]), \
                                                      pdb.conn[k]-set([j])) \
                                           if i != l])

    tor_index = {}
    for i,j,k,l in tors:
        ti = pdb.names[i][2]
        tj = pdb.names[j][2]
        tk = pdb.names[k][2]
        tl = pdb.names[l][2]
        # "Strike that; reverse it." -- Willy Wonka.
        if tj > tk or (tj == tk and ti > tl):
            ti, tj, tk, tl = (tl, tk, tj, ti)
            i, j, k, l = (l, k, j, i)
        name = "%s_%s_%s_%s"%(ti,tj,tk,tl)
        if not tor_index.has_key(name):
            tor_index[name] = []
        tor_index[name].append((i,j,k,l))

    pdb.tors = tors

    print("%d Torsions"%sum(map(len, tor_index.values())))
    terms = [mkterm(n,l) for n,l in tor_index.iteritems()]
    return FFconcat(terms)
