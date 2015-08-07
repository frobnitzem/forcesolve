# 2-body bonded interactions

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

# Represents bonded energy as b-splined f(r^2), r in (0,h*n).

# Adds the "bond" forcefield term into the list of atomic interactions.
class cg_topol:
	__metaclass__=ExtendInplace
	def init_bond(self):
		self.bond = set()
	
	def append_bond(self, line):
		if line[0] == "all":
			self.bond |= self.edge
			return
		else:
			raise RuntimeError, "Error! Non-all bond table not "\
					"supported."
		#rname = line[0]
		#anames = parse_aname(line[1:])
		
		pass
	def add_bond(self, i, j):
		if i < j:
			self.bond.add((i,j))
		else:
			self.bond.add((j,i))
	
	def finalize_bond(self):
		print "Finalizing bonds..."
		
		bspl = Bspline(4) # Default spline order.
		bond = [] # List of atoms involved.
		bond_info = [] # List of bond_type_calc classes.
		bond_index = {}
		for i,j in self.bond:
			ti = self.names[i][2]
			tj = self.names[j][2]
			if ti <= tj:
				name = "%s-%s"%(ti,tj)
				try:
					t = bond_index[name]
				except KeyError:
					t = len(bond)
					bond_index[name] = t
					bond.append([])
					bond_info.append(bond_type_calc(\
							name,bspl))
				bond[t].append((i,j))
			else: # "Strike that; reverse it." -- Willy Wonka.
				name = "%s-%s"%(tj,ti)
				try:
					t = bond_index[name]
				except KeyError:
					t = len(bond)
					bond_index[name] = t
					bond.append([])
					bond_info.append(bond_type_calc(\
							name,bspl))
				bond[t].append((j,i))
		
		self.bond = bond
		self.bond_info = bond_info
	
	def prior_bond(self):
		return [info.prior for info in self.bond_info]
	def constrain_bond(self):
		return [info.constraint for info in self.bond_info]
	
	def read_bond(self, base):
		for info in self.bond_info:
			info.read(base)
	def write_bond(self, base):
		for info in self.bond_info:
			info.write(base)
	
	def calc_bond(self, x):
		bondt = []
		for blist, info in zip(self.bond, self.bond_info):
			bondt.append(bond(array([x[...,j,:]-x[...,i,:] \
						for i,j in blist])))
		return bondt
	
	# The design array multiplies the spline coefficients to produce the
	# total bonded energy/force.
	def design_bond(self, x, order=0):
		A = []
		if order == 0:
		  for blist,info in zip(self.bond, self.bond_info):
		    b = bond(array([x[...,j,:]-x[...,i,:] for i,j in blist]))
		    spl = info.spline(b, order)
		    A.append(sum(spl,-2)) # Sum over all such interactions.
		  return A
		elif order == 1:
		  Ad = [] # "They are still puzzled, Mr. Bond." -- Dr. No
		  for blist,info in zip(self.bond, self.bond_info):
		    ad = zeros(x.shape+(info.f.n,), float)
		    b, db = dbond(array([x[...,j,:]-x[...,i,:] \
						for i,j in blist]))
		    spl, dspl = info.spline(b, order)
		    for k in range(len(blist)):
			i,j = blist[k]
			ad[...,i,:,:] -= db[...,k,:,newaxis] \
					* dspl[...,k,newaxis,:]
			ad[...,j,:,:] += db[...,k,:,newaxis] \
					* dspl[...,k,newaxis,:]
		    A.append(sum(spl,-2))
		    Ad.append(ad)
		  return A, Ad
		else:
		  raise RuntimeError, "Error! >1 energy derivative not "\
					"supported."

	def commit_bond(self):
		for info in self.bond_info:
			info.commit()
		
# Returns the energy of configuration x and/or the conserved force on each atom.
	def bond_energy(self, x):
		en = 0.0
		A = self.calc_bond(x)
		for a, info in zip(A, self.bond_info):
			en += sum(info.f.ay(a),-1)
		return en
	
	def bond_force(self, x):
		en = 0.0
		F = zeros(x.shape, float)
		for blist,info in zip(self.bond, self.bond_info):
			b, db = dbond(array([x[...,j,:]-x[...,i,:] \
						for i,j in blist]))
			u, du = info.f.ay(b, 1)
			en += sum(u,-1)
			for k in range(len(blist)):
				i,j = blist[k]
				F[...,i,:] += du[k]*db[...,k,:]
				F[...,j,:] -= du[k]*db[...,k,:]
		return en,F
	
def bond(x):
	if len(x.shape) > 2:
		trns = range(len(x.shape))
		trns = trns[1:-1]+trns[:1]+trns[-1:]
		x = transpose(x, trns)
	return sqrt(sum(x*x,-1))

def dbond(x):
	if len(x.shape) > 2:
		trns = range(len(x.shape))
		trns = trns[1:-1]+trns[:1]+trns[-1:]
		x = transpose(x, trns)
	# x=j-i => dr(j-i)/dj = x/r
	r = sqrt(sum(x*x,-1))
	x /= r[...,newaxis]
	return r, x

def bond_type_calc(name, bspl):
	n = 40+bspl.order-1
	x0 = 0.0
	L = 4.0
	#x0 = 1.5
	#L = 7.5
	id = "bond_"+name
	f = spline_func((x0,L/(n-(bspl.order-1)),n), bspl, \
		False, bspl.order-1, (x0,x0+L), id)
	m = 2
	return type_calc(id, f, f.integral(0, m), m)
