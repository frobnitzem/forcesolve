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

from cg_topol import *
from bspline import *
from bonds import bond, dbond

# Adds the "pair" forcefield term into the list of atomic interactions.
class cg_topol:
	__metaclass__=ExtendInplace
	def init_pair(self):
		self.pair = set()
		if not hasattr(self, "periodic"): # Don't know append_ order.
			self.periodic = False
	
	def append_pair(self, line):
		if line[0][0:3] == "all":
			if len(line[0]) <= 3 or line[0][3] not in "0123456789":
				print "PAIR all? statement lacks number, "\
					"ignoring..."
			skip = int(line[0][3:])-1
			self.add_all_pairs(skip)
			return
		else:
			raise RuntimeError, "Error! Non-all pair table not "\
					"supported."
		pass
	
	def append_periodic(self, line):
		if len(line) < 3:
			raise InputError, "Error parsing PERIODIC line, "\
					"should have 3 box lengths!"
		self.periodic = True
		self.box_L = array([float(line[0]), \
				float(line[1]), float(line[2])])
	
	def add_all_pairs(self, skip):
		xpair = []
		for a in range(self.atoms):
			xpair.append(set([a]))
		for i in range(1,skip): # Extend table by 1 bond.
		    for a in range(self.atoms):
			xpair[a] |= reduce(lambda x,y: x|y, \
					[self.conn[b] for b in xpair[a]])
		self.xpair = xpair
		xpair = reduce(lambda x,y: x|y, \
				[modprod([a],x) for a,x in enumerate(xpair)])
		pair = []
		for i in range(self.atoms-1):
		    pair += [(i,j) for j in range(i+1,self.atoms)]
		self.pair = set(pair)-xpair
	def add_pair(self, i, j):
		if i < j:
			self.pair.add((i,j))
		else:
			self.pair.add((j,i))
	
	def finalize_pair(self):
		print "Finalizing pairs..."
		
		bspl = Bspline(6) # Default bspline.
		
		pair = [] # List of atoms involved.
		pair_info = [] # List of bond_type_calc classes.
		pair_index = {}
		for i,j in self.pair:
			ti = self.names[i][2]
			tj = self.names[j][2]
			if ti <= tj:
				name = "%s-%s"%(ti,tj)
				try:
					t = pair_index[name]
				except KeyError:
					t = len(pair)
					pair_index[name] = t
					pair.append([])
					pair_info.append(pair_type_calc(\
						name, bspl))
				pair[t].append((i,j))
			else: # "Strike that; reverse it." -- Willy Wonka.
				name = "%s-%s"%(tj,ti)
				try:
					t = pair_index[name]
				except KeyError:
					t = len(pair)
					pair_index[name] = t
					pair.append([])
					pair_info.append(pair_type_calc(\
						name, bspl))
				pair[t].append((j,i))
		
		self.pair = pair
		self.pair_info = pair_info
	
	def prior_pair(self):
		return [info.prior for info in self.pair_info]
	def constrain_pair(self):
		return [info.constraint for info in self.pair_info]
	
	def read_pair(self, base):
		for info in self.pair_info:
			info.read(base)
			info.f.verb = False
	def write_pair(self, base):
		for info in self.pair_info:
			info.write(base)
	
	def calc_pair(self, x):
		pair = []
		for plist, info in zip(self.pair, self.pair_info):
		    delta = array([x[...,j,:]-x[...,i,:] \
				for i,j in plist])
		    if self.periodic: # Wrap for periodicity.
			delta -= self.box_L * floor(delta/self.box_L+0.5)
		    pair.append(bond(delta))
		return pair
	
	# The design array multiplies the spline coefficients to produce the
	# total bonded energy/force.
	def design_pair(self, x, order=0):
		A = []
		if order == 0:
		  for plist,info in zip(self.pair, self.pair_info):
		    delta = array([x[...,j,:]-x[...,i,:] \
				for i,j in plist])
		    if self.periodic: # Wrap for periodicity.
			delta -= self.box_L * floor(delta/self.box_L+0.5)
		    A.append(sum(info.spline(bond(delta), order),-2))
		  return A
		elif order == 1:
		  Ad = []
		  for plist,info in zip(self.pair, self.pair_info):
		    ad = zeros(x.shape+(info.f.n,), float)
		    delta = array([x[...,j,:]-x[...,i,:] \
				for i,j in plist])
		    if self.periodic: # Wrap for periodicity.
			delta -= self.box_L * floor(delta/self.box_L+0.5)
		    b,db = dbond(delta)
		    spl, dspl = info.spline(b, order)
		    for k in range(len(plist)):
			i,j = plist[k]
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

	def commit_pair(self):
		for info in self.pair_info:
			info.commit()
	
# Returns the energy of configuration x and/or the conserved force on each atom.
	def pair_energy(self, x):
		en = 0.0
		A = self.calc_pair(x)
		for a, info in zip(A, self.pair_info):
			en += sum(info.f.ay(a),-1)
		return en
	
	def pair_force(self, x):
		en = 0.0
		F = zeros(x.shape, float)
		for plist,info in zip(self.pair, self.pair_info):
		    delta = array([x[...,j,:]-x[...,i,:] \
				for i,j in plist])
		    if self.periodic: # Wrap for periodicity.
			delta -= self.box_L * floor(delta/self.box_L+0.5)
		    b,db = dbond(delta)
		    u, du = info.f.ay(b, 1)
		    en += sum(u,-1)
		    for k in range(len(plist)):
			i,j = plist[k]
			F[...,i,:] += du[k]*db[...,k,:]
			F[...,j,:] -= du[k]*db[...,k,:]
		return en,F
		
def pair_type_calc(name, bspl):
	#n = 100
	#L = 8.4
	n = 170
	L = 17.0
	id = "pair_"+name
	f = spline_func((0.,L/n,n), bspl, False, bspl.order-1, (0.,inf), \
				id, False) # Ignore out-of-range errors.
	return type_calc(id, f, None, 2)
