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

from cg_topol import *
from bspline import *

# Adds the "tor" forcefield term into the list of atomic interactions.
class cg_topol:
	__metaclass__=ExtendInplace
	def init_tor(self):
		self.tor = set()
	
	def append_tor(self, line):
		if line[0] == "all":
			self.add_all_tor()
			return
		else:
			raise RuntimeError, "Error! Non-all bond table not "\
					"supported."
		pass
	def add_all_tor(self):
		for j,k in self.edge:
		    self.tor |= \
		     set([(i,j,k,l) for i,l in modprod(self.conn[j]-set([k]), \
				self.conn[k]-set([j])) if i != l])
		   
	def add_tor(self, i, j, k, l):
		if j < k:
			self.tor.add((i,j,k,l))
		else:
			self.tor.add((l,k,j,i))
	
	def finalize_tor(self):
		print "Finalizing torsions..."
		
		spl = Bspline(4) # Default spline.
		
		tor = [] # List of atoms involved.
		tor_info = [] # List of bond_type_calc classes.
		tor_index = {}
		for i,j,k,l in self.tor:
			ti = self.names[i][2]
			tj = self.names[j][2]
			tk = self.names[k][2]
			tl = self.names[l][2]
			if tj < tk or (tj == tk and ti<=tl):
				name = "%s-%s-%s-%s"%(ti,tj,tk,tl)
				try:
					t = tor_index[name]
				except KeyError:
					t = len(tor)
					tor_index[name] = t
					tor.append([])
					tor_info.append(tor_type_calc(\
							name, spl))
				tor[t].append((i,j,k,l))
			else: # "Strike that; reverse it." -- Willy Wonka.
				name = "%s-%s-%s-%s"%(tl,tk,tj,ti)
				try:
					t = tor_index[name]
				except KeyError:
					t = len(tor)
					tor_index[name] = t
					tor.append([])
					tor_info.append(tor_type_calc(\
							name, spl))
				tor[t].append((l,k,j,i))
		
		self.tor = tor
		self.tor_info = tor_info
	
	def prior_tor(self):
		return [info.prior for info in self.tor_info]
	def constrain_tor(self):
		return [info.constraint for info in self.tor_info]
	
	def read_tor(self, base):
		for info in self.tor_info:
			info.read(base)
	def write_tor(self, base):
		for info in self.tor_info:
			info.write(base)
	
	def calc_tor(self, x):
		tor = []
		for alist,info in zip(self.tor, self.tor_info):
		  tor.append(torsion(array([[x[...,i,:]-x[...,j,:],\
					x[...,j,:]-x[...,k,:],\
					x[...,l,:]-x[...,k,:]] \
						for i,j,k,l in alist])))
		return tor
	
	# The design array multiplies the spline coefficients to produce the
	# total bonded energy/force.
	def design_tor(self, x, order=0):
		A = []
		if order == 0:
		  for alist,info in zip(self.tor, self.tor_info):
		    tor = torsion(array([[x[...,i,:]-x[...,j,:],\
					x[...,j,:]-x[...,k,:],\
					x[...,l,:]-x[...,k,:]] \
						for i,j,k,l in alist]))
		    A.append(sum(info.spline(tor, order),-2))
		  return A
		elif order == 1:
		  Ad = []
		  for alist,info in zip(self.tor, self.tor_info):
		    ad = zeros(x.shape+(info.f.n,), float)
		    t, dt = dtorsion(array([[x[...,i,:]-x[...,j,:],\
					x[...,j,:]-x[...,k,:],\
					x[...,l,:]-x[...,k,:]] \
						for i,j,k,l in alist]))
		    spl, dspl = info.spline(t, order)
		    for z in range(len(alist)):
			i,j,k,l = alist[z]
			ad[...,i,:,:] += dt[...,z,0,:,newaxis] \
						* dspl[...,z,newaxis,:]
			ad[...,l,:,:] += dt[...,z,1,:,newaxis] \
						* dspl[...,z,newaxis,:]
			ad[...,j,:,:] -= dt[...,z,2,:,newaxis] \
						* dspl[...,z,newaxis,:]
			ad[...,k,:,:] += dt[...,z,3,:,newaxis] \
						* dspl[...,z,newaxis,:]
		    A.append(sum(spl,-2))
		    Ad.append(ad)
		  return A, Ad
		else:
		  raise RuntimeError, "Error! >1 energy derivative not "\
					"supported."

	def commit_tor(self):
		for info in self.tor_info:
			info.commit()
		
# Returns the energy of configuration x and/or the conserved force on each atom.
	def tor_energy(self, x):
		A = self.calc_tor(x)
		en = 0.0
		for a, info in zip(A, self.tor_info):
			en += sum(info.f.ay(a),-1)
		return en
	
	def tor_force(self, x):
		en = 0.0
		F = zeros(x.shape, float)
		for alist,info in zip(self.tor, self.tor_info):
		    t, dt = dtorsion(array([[x[...,i,:]-x[...,j,:],\
					x[...,j,:]-x[...,k,:],\
					x[...,l,:]-x[...,k,:]] \
						for i,j,k,l in alist]))
		    u, du = info.f.ay(t, 1)
		    en += sum(u,-1)
		    for z in range(len(alist)):
			i,j,k,l = alist[z]
			F[...,i,:] -= du[z]*dt[...,z,0,:]
			F[...,l,:] -= du[z]*dt[...,z,1,:]
			F[...,j,:] += du[z]*dt[...,z,2,:]
			F[...,k,:] -= du[z]*dt[...,z,3,:]
		return en,F
	
def tor_type_calc(name, bspl):
	n = 180
	id = "tor_"+name
	f = spline_func((-pi, 2.0*pi/n, n), bspl, True, \
				bspl.order*0.5, (-pi, pi), id)
	m = 0
	return type_calc(id, f, f.integral(0, m), m)
	
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
	
	A = zeros((4,3) + x.shape[2:], float)
	
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

