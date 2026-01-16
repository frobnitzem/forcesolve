# Langevin dynamics class using B-spline topology files.

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
# 5/12/2008
# This work was supported by a DOE CSGF.

import os
from numpy import *
import numpy.random as rand
from cg_topol.ucgrad.array_io import *

class md_integrate:
        def __init__(self, x0, v0, topol, dt, kT):
                self.topol = topol
                self.topol.read_all_splines()
                assert topol.can_do_md()
                
                self.dt = dt
                self.kT = kT
                self.build_type_index()
                self.read_iz(topol.const_dir)
                
                self.calc_factors()
                self.set_x(x0)
                if v0 != None:
                        self.set_v(v0)
                else:
                        self.v = rand.standard_normal(self.x.shape) \
                                * self.fac[:,newaxis]
                        self.v -= sum(self.v*self.mass[:,newaxis], 0) \
                                / sum(self.mass)
                
        def set_x(self, x): # Read and non-dimensionalize.
                self.x = x.copy()
        def set_v(self, v):
                self.v = v.copy()
        
        def do_md(self, ns, nout, out):
                for ext in [".pdb",".x",".v",".f",".en"]:
                        of = open(out+ext, 'w')
                        of.truncate()
                        of.close()
                
                self.log_state(out, 1)
                for i in range(ns):
                        self.integrate(nout)
                        self.log_state(out, i+2)
        
        def log_state(self, out, n=1):
                self.topol.pdb.x = self.x*.5291772083
                self.topol.pdb.write(out+".pdb", 'a', n)
                write_matrix(out+".x", self.x, 'a')
                write_matrix(out+".v", self.v, 'a')
                
                E, F = self.topol.calc_force(self.x)
                write_matrix(out+".en", array([[n,E]]), 'a')
                write_matrix(out+".f", F, 'a')
        
        def integrate(self, ns=1):
                assert ns > 0
                
                self.first_half_step()
                for i in range(1,ns):
                        self.full_step()
                self.last_half_step()
        
        def first_half_step(self):
                E, F = self.topol.calc_force(self.x)
                F *= self.dt*0.5/self.mass[:,newaxis]
                R = rand.standard_normal(self.v.shape)*self.lam[:,newaxis]
                #R = 0.0
                self.v = self.vdamp[:,newaxis]*(self.v + F+R)
                self.lR = R
                self.x += self.xdamp[:,newaxis]*self.v
        
        def last_half_step(self):
                E, F = self.topol.calc_force(self.x)
                F *= self.dt*0.5/self.mass[:,newaxis]
                self.v = self.vdamp[:,newaxis]*self.v + F+self.lR
                
        def full_step(self):
                E, F = self.topol.calc_force(self.x)
                F *= self.dt/self.mass[:,newaxis]
                R = rand.standard_normal(self.v.shape)*self.lam[:,newaxis]
                #R = 0.0
                self.v = self.vdamp[:,newaxis]*(self.vdamp[:,newaxis]*self.v\
                                + F + self.lR + R)
                self.lR = R
                self.x += self.xdamp[:,newaxis]*self.v
        
        # Build an index to atoms by designated atom type.
        def build_type_index(self):
                self.type_names = list(set( [aname[2] for \
                                aname in self.topol.names] ))
                self.types = len(self.type_names)
                #print(self.type_names)
                
                type_num = {} # Dictionary for quick reference.
                self.tmass = ones(self.types)
                for i,t in enumerate(self.type_names):
                        type_num[t] = i
                        if not self.topol.tmass.has_key(t):
                                print("Warning! No mass has been input for "\
                                        "atom type %s: assuming mass=1."%(t))
                        else:
                                self.tmass[i] = self.topol.tmass[t]
                #print(type_num)
                
                self.type_index = []
                for aname in self.topol.names:
                        tn = type_num[aname[2]]
                        self.type_index.append(tn)
        
        def calc_factors(self): # Factor giving units of x
                self.mass = zeros(self.topol.atoms)
                self.z = zeros(self.topol.atoms)
                for i,t in enumerate(self.type_index):
                        self.z[i] = self.zt[t]
                        self.mass[i] = self.tmass[t]
                self.fac = sqrt(self.kT/self.mass) # sqrt(kT/m)
                #self.ifac = 1.0/self.fac
                
                # Magnitude of random force.
                self.lam = sqrt(0.5*self.z)*self.fac
                self.vdamp = exp(-0.5*self.z) # Integration damping factors.
                self.xdamp = self.dt*(1.0-exp(-self.z))/(self.z*self.vdamp)
                #print self.z, self.mass, self.lam, self.vdamp, self.xdamp
        
        def read_iz(self, base):
                #self.zt = zeros(len(self.type_names))+1.e-12*self.dt
                #return
                try:
                    z = open(os.path.join(base,"v.out"))
                except IOError:
                    self.zt = zeros(len(self.type_names))+1.e-12*self.dt
                    print("Warning! Setting damping const. to 1.0e-12 !")
                    return
                
                self.zt = zeros(len(self.type_names))
                marked = zeros(self.zt.shape, int)
                for line in z.xreadlines():
                        tok = line.split()
                        if not tok or tok[0][0] == '#':
                                continue
                        try:
                                i = self.type_names.index(tok[0])
                        except ValueError:
                                continue
                        self.zt[i] = float(tok[1])
                        marked[i] = 1
                assert all(marked), "Didn't find z-records for all atom types."
