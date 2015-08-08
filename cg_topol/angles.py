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

from cg_topol import *
from bspline import *

# Adds the "angle" forcefield term into the list of atomic interactions.
class cg_topol:
	__metaclass__=ExtendInplace
	def init_angle(self):
		self.angle = set()
	
	def append_angle(self, line):
		if line[0] == "all":
		    for j in range(self.atoms):
			angles = [(i,j,k) for i,k in modprod(self.conn[j], \
					self.conn[j]) if not i>=k]
			self.angle |= set(angles)
		    return
		else:
		    raise RuntimeError, "Error! Non-all angle table not "\
					"supported."
		pass
	def add_angle(self, i, j, k):
		if i < k:
			self.angle.add((i,j,k))
		else:
			self.angle.add((k,j,i))
	
	def finalize_angle(self):
		print "Finalizing angles..."
		bspl = Bspline(4) # Default bspline.
		
		angle = [] # List of atoms involved.
		angle_info = [] # List of bond_type_calc classes.
		angle_index = {}
		for i,j,k in self.angle:
			ti = self.names[i][2]
			tj = self.names[j][2]
			tk = self.names[k][2]
			if ti <= tk:
				name = "%s-%s-%s"%(ti,tj,tk)
				try:
					t = angle_index[name]
				except KeyError:
					t = len(angle)
					angle_index[name] = t
					angle.append([])
					angle_info.append(angle_type_calc(\
							name, bspl))
				angle[t].append((i,j,k))
			else: # "Strike that; reverse it." -- Willy Wonka.
				name = "%s-%s-%s"%(tk,tj,ti)
				try:
					t = angle_index[name]
				except KeyError:
					t = len(angle)
					angle_index[name] = t
					angle.append([])
					angle_info.append(angle_type_calc(\
							name, bspl))
				angle[t].append((k,j,i))
		
		self.angle = angle
		self.angle_info = angle_info
	
	def prior_angle(self):
		return [info.prior for info in self.angle_info]
	def constrain_angle(self):
		return [info.constraint for info in self.angle_info]
	
	def read_angle(self, base):
		for info in self.angle_info:
			info.read(base)
	def write_angle(self, base):
		for info in self.angle_info:
			info.write(base)
	
	def calc_angle(self, x):
		ang = []
		for alist,info in zip(self.angle, self.angle_info):
			ang.append( anglec(array([[x[...,i,:]-x[...,j,:],\
				x[...,k,:]-x[...,j,:]] for i,j,k in alist])) )
		return ang
	
	# The design array multiplies the spline coefficients to produce the
	# total bonded energy/force.
	def design_angle(self, x, order=0):
		A = []
		if order == 0:
		  for alist,info in zip(self.angle, self.angle_info):
		    a = anglec(array([[x[...,i,:]-x[...,j,:],\
			x[...,k,:]-x[...,j,:]] for i,j,k in alist]))
		    spl = info.spline(a, order)
		    A.append(sum(spl,-2))# Sum over all atoms in ea. structure.
		  return A
		elif order == 1:
		  Ad = []
		  for alist,info in zip(self.angle, self.angle_info):
		    ad = zeros(x.shape+(info.f.n,), float)
		    a,da = danglec(array([[x[...,i,:]-x[...,j,:],\
			x[...,k,:]-x[...,j,:]] for i,j,k in alist]))
		    spl, dspl = info.spline(a, order)
		    for z in range(len(alist)):
			i,j,k = alist[z]
			ad[...,i,:,:] += da[...,z,0,:,newaxis] \
					* dspl[...,z,newaxis,:]
			ad[...,j,:,:] -= sum(da[...,z,:,:],-2)[...,newaxis] \
					* dspl[...,z,newaxis,:]
			ad[...,k,:,:] += da[...,z,1,:,newaxis] \
					* dspl[...,z,newaxis,:]
		    A.append(sum(spl,-2))
		    Ad.append(ad)
		  return A,Ad
		else:
		  raise RuntimeError, "Error! >1 energy derivative not "\
					"supported."

	def commit_angle(self):
		for info in self.angle_info:
			info.commit()
		
# Returns the energy of configuration x and/or the conserved force on each atom.
	def angle_energy(self, x): # Requires commit.
		en = 0.0
		A = self.calc_angle(x)
		for a, info in zip(A, self.angle_info):
			en += sum(info.f.ay(a),-1)
		return en
	
	def angle_force(self, x):
		en = 0.0
		F = zeros(x.shape, float)
		for alist,info in zip(self.angle, self.angle_info):
			a, da = danglec(array([[x[...,i,:]-x[...,j,:],\
				x[...,k,:]-x[...,j,:]] for i,j,k in alist]))
			u, du = info.f.ay(a, 1)
			en += sum(u,-1)
			for z in range(len(alist)):
				i,j,k = alist[z]
				F[...,i,:] -= du[z]*da[...,z,0,:]
				F[...,j,:] += du[z]*sum(da[...,z,:,:],-2)
				F[...,k,:] -= du[z]*da[...,z,1,:]
		return en,F
	
def angle_type_calc(name, bspl):
	n = 90+bspl.order-1
	id = "angle_"+name
	f = spline_func((-1.0, 2.0/(n-(bspl.order-1)), n), bspl, \
		False, bspl.order-1, (-1.0, 1.0), id)
	m = 0
	return type_calc(id, f, f.integral(0, m), m)
	
def anglec(x):
	if len(x.shape) > 3:
		trns = range(len(x.shape))
		trns = trns[2:-1]+trns[:2]+trns[-1:]
		x = transpose(x, trns)
	return sum(x[...,0,:]*x[...,1,:],-1) \
	  / sqrt(sum(x[...,0,:]*x[...,0,:],-1)*sum(x[...,1,:]*x[...,1,:],-1))
	
def danglec(x):
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
	return cth, A

