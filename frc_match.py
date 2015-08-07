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
# 12/10/2007
# This work was supported by a DOE CSGF.

from numpy import *
import numpy.linalg as la
import numpy.random as rand
from cg_topol.ucgrad import write_matrix, write_list

# cg_topol uses integrated first, second and third derivatives for prior which,
# when multiplied by the corresponding values in alpha, serve to push the
# energy function toward a flat, linear, or quadratic shape -- respectively.

class frc_match:
	def __init__(self, topol, dt, kT, \
			E0=1.0e-8, calpha=1.0e+3, logalpha=None):
		# No inputs defined yet...
		self.topol = topol
		assert topol.can_force_match()
		
		self.logalpha = logalpha
		self.build_type_index()
		types = self.types
		params = self.topol.parameters
		
		self.dt = dt
		self.kT = kT
		
		self.D = zeros((types,params), float)
		self.D2 = zeros((types,params,params), float)
		self.DF = zeros((types,params), float)
		self.F2 = zeros((types), float)
		self.S = 0
		
		self.E0 = E0
		self.calpha = calpha
		self.prior, self.prior_rank, self.prng = self.topol.calc_prior()
		for i in range(len(self.prior)):
		    if self.prior_rank[i] != len(self.prior[i]):
			val, A = la.eigh(self.prior[i])
			for j in range(len(self.prior[i])-self.prior_rank[i]):
				val[j] = self.E0*self.E0*0.01
			self.prior[i] = dot(A*val[newaxis,:], transpose(A))
			#A = A[:,-self.prior_rank[i]:]
			#val = val[-self.prior_rank[i]:]
			#self.prior[i] = dot(A*val[newaxis,:], transpose(A))
		self.alpha = ones(len(self.prior))
		
		constr = self.topol.calc_constraints()
		# Normalize.
		self.orthonormalize_constraints(constr)
		
		self.theta = self.theta_from_topol() # Parameters.
		self.z = ones(types, float)
		
		# Sampling accumulators.
		self.sum_v = []
		self.sum_alpha = []
		self.sum_theta = []
		self.sum_dmu2 = []
		self.sum_df2 = []
		self.df2 = zeros(len(self.topol.nt))
		self.samples = 0
	
	# Build an index to atoms by designated atom type.
	def build_type_index(self):
		self.type_names = list(set( [aname[2] for \
                                aname in self.topol.names] ))
		self.types = len(self.type_names)
		#print self.type_names
		
		type_num = {} # Dictionary for quick reference.
		mass = ones(self.types)
		for i,t in enumerate(self.type_names):
			type_num[t] = i
			if not self.topol.tmass.has_key(t):
				print "Warning! No mass has been input for "\
					"atom type %s: assuming mass=1."%(t)
				pass
			else:
				mass[i] = self.topol.tmass[t]
		#print type_num
		#print mass
		
		nt = [0]*self.types
		
		self.type_index = []
		for aname in self.topol.names:
			tn = type_num[aname[2]]
			nt[tn] += 1
			self.type_index.append(tn)
		self.mass = zeros(self.topol.atoms)
		for i,t in enumerate(self.type_index):
			self.mass[i] = mass[t]
		self.nt = array(nt)
		
	# Estimate the actual dimensionality of iC.
	def dimensionality(self, sigma_tol=1.0e-5):
		tol = sigma_tol*sigma_tol
		itol = 1.0/tol
		iC = self.calc_iC()
		#print self.constraints
		w, A = la.eigh(iC)
		free_dim = len([l for l in w if l <= tol])
		fixed_dim = len([l for l in w if l >= itol])
		print "Free directions = %d, Fixed directions = %d"%(\
					free_dim, fixed_dim)
		if free_dim > 0:
		    for i in range(len(w)):
			if w[i] > tol:
				continue
			print "\tValue %e == "%(w[i]) +\
			  str(dot(dot(iC,A[:,i]), A[:,i]))
			w[i] = 10.0*itol
			#print "\tVector:"
			#print A[:,i]
		return w, A

# Append a set of data points to the present frc_match object.
# Takes as input data from a single specified atom type which may contain
# differing remote interaction types
	def append(self, x,f):
		chunk = 2
		
		if self.samples > 0:
			raise runtimeError, "Error! Cannot add more data "\
				"points to mature frc_match object."
		if x.shape != f.shape:
			print x.shape, f.shape
			raise runtimeError, "Error! x and f trajectory shapes "\
				"differ!"
		if x.shape[1] != self.topol.atoms:
			raise InputError, "Error! Number of atoms in "\
				"trajectory does not match topology!"
		
		# Multiply by ugly constants here.
		Xfac = self.dt*sqrt(self.kT/self.mass) # Non-dimensionalize.
		f *= (self.dt/sqrt(self.mass*self.kT))[newaxis,:,newaxis]
		
		print "Appending %d samples..."%(len(x))
		self.S += len(x)
		# Operate on "chunk" structures at once.
		for i in range(0,len(x)-chunk+1,chunk):
		    D = -1.0*self.topol.calc_design(x[i:i+chunk],1) \
			* Xfac[newaxis,:,newaxis,newaxis] # Factor cancels 1/dx
		    self.D += self.type_sum(sum(sum(D,-2),0))
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:i+chunk])
		    self.F2 += self.type_sum(sum(sum(\
			f[i:i+chunk]*f[i:i+chunk],-1),0))
		i = len(x)%chunk
		if i != 0:
		    i = len(x)-i
		    D = -1.0*self.topol.calc_design(x[i:],1) \
			* Xfac[newaxis,:,newaxis,newaxis] # Factor cancels 1/dx
		    self.D += self.type_sum(sum(sum(D,-2),0))
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:])
		    self.F2 += self.type_sum(sum(sum(f[i:]*f[i:],-1),0))
		
		# 1 structure at a time.
		#for xi, fi in zip(x,f):# Design matrices are for energy deriv.s
		#	D = -1.0*self.topol.calc_design(xi,1) # -ize -> F.
		#	self.D += self.type_sum(sum(D, 1))
		#	self.D2 += self.type_sum(sum(\
		#		D[:,:,:,newaxis]*D[:,:,newaxis,:],1))
		#	self.DF += self.type_sum(sum(D*fi[:,:,newaxis],1))
		#	self.F2 += self.type_sum(sum(fi*fi,1))
	
	# Given an (atoms by ?) array, reduce to a (types by ?) array
	# by summation.
	def type_sum(self, x):
		xt = zeros((self.types,)+x.shape[1:], float)
		for i,t in enumerate(self.type_index):
			xt[t] += x[i]
		return xt
	# Special type_sum for accumulating D2 matrices.
	def type_sum_D2(self, DS):
		for D in DS:
		  for i,t in enumerate(self.type_index):
		    for j in range(3):
			self.D2[t] += D[i,j,:,newaxis]*D[i,j,newaxis,:]
	# Special type_sum for accumulating DF matrices.
	def type_sum_DF(self, DS, FS):
		for D, F in zip(DS, FS):
		  for i,t in enumerate(self.type_index):
		    for j in range(3):
			self.DF[t] += D[i,j]*F[i,j,newaxis]
	
	# Make constraints orthonormal.
	def orthonormalize_constraints(self, constraints):
		if constraints == None:
			self.constraints = zeros((1,self.topol.parameters), \
							float)
			return
		
		self.constraints = orthonormalize(constraints)
	
	# Find maximum likelihood estimate.
	def maximize(self, tol=1.0e-5, maxiter=1000):
		if self.S < 1:
			raise ProgramError, "Error! no data has been collected!"
		
		print "Maximizing posterior PDF..."
		az, ibz, aa, iba = self.calc_za_ab()
		# const. - log(P)
		llp = dot(ibz,self.z) - dot(az-1,log(self.z)) \
			+ dot(iba,self.alpha) - dot(aa-1,log(self.alpha))
		
		delta = tol + 1.0
	        iter = 0
		while delta > tol and iter < maxiter:
		    iter += 1
		    iC = self.calc_iC()
		    theta = dot(self.z, self.DF)
		    theta = la.solve(iC, theta)
		    self.theta = theta
		    
		    az, ibz, aa, iba = self.calc_za_ab()
		    self.z = (az-1.0)/ibz
		    self.alpha = (aa-1.0)/iba
		    
		    lp = dot(ibz,self.z) - dot(az-1,log(self.z)) \
			+ dot(iba,self.alpha) - dot(aa-1,log(self.alpha))
		    delta = llp-lp
		    llp = lp
		    print "  Iteration %d, delta = %e"%(iter,delta)
		
		theta, dmu2, df2 = self.calc_theta_stats()
		df2 = array(df2)/( self.S*3.0*sum(self.nt) )
		self.df2 = df2
	
	def calc_penalty(self):
		pen = zeros(len(self.prior))
		i = 0
		for P, r in zip(self.prior, self.prng):
			pen[i] = dot(dot(P, self.theta[r[0]:r[1]]), \
					self.theta[r[0]:r[1]])
			i += 1
		return pen*(pen > 0.0) # Forces negative pen -> 0
	
	# Note: const(D)-log(P) = dot(bz,self.z) - dot(az-1,log(self.z))\ 
	#                + dot(ba,self.alpha) - dot(aa-1,log(self.alpha))
	def calc_za_ab(self):
		self.ft2 = dot(dot(self.D2,self.theta), self.theta)
		bz = (self.ft2+self.F2) - 2*dot(self.DF, self.theta)
		bz = 0.5*(bz*(bz > 0.0) + self.E0)
		
		ba = 0.5*(self.calc_penalty()+self.E0)
		return 1.5*self.nt*self.S, bz, 0.5*self.prior_rank, ba
		
	def calc_iC(self):
		iC = sum(self.z[:,newaxis,newaxis]*self.D2, 0)
		# Add in constraints.
		iC += self.calpha*dot(transpose(self.constraints), \
					self.constraints)
		# Add in prior info.
		for a, P, r in zip(self.alpha, self.prior, self.prng):
			iC[r[0]:r[1],r[0]:r[1]] += a*P
		return iC
	
	# Generate conditional samples.
	def update_sample(self, n=0, logalpha=None):
	    for i in xrange(n):
		#print "    Updating."
		iC = self.calc_iC()
		b = zeros((self.topol.parameters,2), float)
		b[:,0] = rand.standard_normal(self.topol.parameters) # Sample
		b[:,1] = dot(self.z, self.DF) # Mean
		try:
			L = la.cholesky(iC)
			b[:,1:] = forward_subst(L, b[:,1:])
			b = back_subst(transpose(L), b)
		except la.linalg.LinAlgError:
			w, A = self.dimensionality()
			print "\nForce design matrix is degenerate!"
			hC = dot(A, transpose(A)/sqrt(w)[:,newaxis])
			b = dot(hC, b)
			b[:,1] = dot(hC, b[:,1])
			#import sys
			#sys.exit()
		self.mean = b[:,1]
		self.theta = b[:,0]+b[:,1]
		
		az, ibz, aa, iba = self.calc_za_ab()
		self.z = array([rand.gamma(ai,1.0) for ai in az])/ibz
		self.alpha = array([rand.gamma(ai,1.0) for ai in aa])/iba
		if logalpha != None:
		    write_matrix(logalpha, reshape(self.alpha, (1,-1)), 'a')
		
	def calc_theta_stats(self):
		iC = self.calc_iC()
		b = dot(self.z, self.DF)
		try:
			#theta = la.solve(iC, b)
			C = la.inv(iC)
		except la.linalg.LinAlgError:
			w, A = self.dimensionality()
			C = dot(A, transpose(A)/w[:,newaxis])
		theta = dot(C, b)
		#return theta, [trace(la.solve(iC, D2)) for D2 in self.D2]
		fv = []
		i = 0
		for ti in self.topol.nt:
			ip = i+ti
			D2t = sum(self.D2[:,i:ip,i:ip],0)
			fv.append(trace(dot(C[i:ip,i:ip],D2t)))
			i = ip
		return theta, [trace(dot(C, D2)) for D2 in self.D2], fv
	
	def sample(self, samples, skip=100, toss=10):
		if self.S < 1:
			raise ProgramError, "Error! no data has been collected!"
		if skip > 0:
			print "Doing sampling burn-in..."
		for i in xrange(skip):
			self.update_sample(toss)
		print "Collecting %d samples..."%samples
		if self.logalpha != None:
			out = open(self.logalpha, 'w')
			out.truncate()
			out.close()
		for i in xrange(samples):
			self.update_sample(toss, self.logalpha)
			theta, dmu2, df2 = self.calc_theta_stats()
			self.sum_theta.append(theta)
			self.sum_dmu2.append(dmu2)
			self.sum_df2.append(df2)
			self.sum_v.append(1.0/self.z)
			self.sum_alpha.append(self.alpha)
			self.samples += 1
		self.posterior_estimate()
	
	# Best estimates from posterior distribution.
	def posterior_estimate(self):
		if self.samples < 1:
			raise ProgramError, "Error! No samples collected!"
		self.theta = sum(self.sum_theta,0)/self.samples
		dmu2 = sum(array(self.sum_dmu2), 0)
		df2 = sum(array(self.sum_df2), 0)
		for dmu in array(self.sum_theta)-self.theta:
			dmu2 += dot(dot(self.D2, dmu),dmu)
			f2 = []
			i = 0
			for ti in self.topol.nt:
				ip = i+ti
				D2t = sum(self.D2[:,i:ip,i:ip],0)
				f2.append(dot(dot(D2t, dmu[i:ip]), dmu[i:ip]))
				i = ip
			df2 += array(f2)
		dmu2 /= self.samples*self.S*3.0*self.nt
		df2 /= self.samples*self.S*3.0*sum(self.nt)
		v = sum(self.sum_v,0)/self.samples + dmu2
		self.df2 = df2
		#print "LINPROB = %e"%(float(sum(array(self.sum_E)<v*1e-3))\
		#	/ self.samples)
		#print "ERFAC = %e"%(self.ab0*2.0/v) # product of z and E0
		self.z = 1.0/v
	
	def show_index(self):
		params = 0
		print "\t     name        params  number present"
		for k in self.topol.ff_terms:
			t = getattr(self.topol, "%s"%k)
			if len(t) < 1:
				continue
			
			print "%s:"%k
			tinfo = getattr(self.topol, "%s_info"%k)
			for tlist, info in zip(t,tinfo):
				print "\t%-16s %5d %5d"%(info.type_name, \
					info.f.n, len(tlist))
				params += info.f.n
		assert(params == self.topol.parameters)
		print "%d total parameters."%params
		print "%d total independent linear constraints."%\
				len(self.constraints)
		print
	
	def theta_from_topol(self):
		theta = zeros(self.topol.parameters, float)
		i = 0
		for k in self.topol.ff_terms:
			tinfo = getattr(self.topol, "%s_info"%k)
			for info in tinfo:
				np = len(info.f.c)
				theta[i:i+np] = info.f.c
				i += np
		assert(i == len(theta))
		
		return theta
		
	def theta_to_topol(self, theta):
		# Put constants into topol.
		i = 0
		for k in self.topol.ff_terms:
			tinfo = getattr(self.topol, "%s_info"%k)
			for info in tinfo:
				np = len(info.f.c)
				info.f.c = theta[i:i+np]
				i += np
		assert(i == len(theta))
	
	def write_out(self, name):
		self.theta_to_topol(self.theta*self.kT) # Dimensionalize
		# Use topol's own write methods.
		self.topol.write_all_splines(name)
		
		#write_matrix(name+"alpha.out", reshape(self.alpha,(1,-1)), 'a')
		
		if self.samples > 0:
			avg_v = sum(self.sum_v,0)/self.samples
			s_v = sum((array(self.sum_v)-avg_v)**2,0)
			s_v = sqrt(s_v/self.samples)
		else:
			avg_v = 1./self.z
			s_v = zeros(len(self.z))
		lam = open(name+"v.out", 'w')
		lam.write("#type\tv\t<v>\tsigma_v\n")
		for t,l2,avg,sd in zip(self.type_names, 1./self.z, \
							avg_v, s_v):
			lam.write("%s\t%e\t%e\t%e\n"%(t,l2,avg,sd))
		lam.close()
		
		df = open(name+"df.out", 'w')
		i = 0
		df.write("#type\t<stdev>\n")
		for k in self.topol.ff_terms:
			t = getattr(self.topol, "%s"%k)
			if len(t) < 1:
				continue
			tinfo = getattr(self.topol, "%s_info"%k)
			for tlist, info in zip(t,tinfo):
				df.write( "%-16s %e\n"%(info.type_name, \
						sqrt(self.df2[i])) )
				i += 1
		assert i == len(self.df2)
		df.close()

# Operates on row space of A
def orthonormalize(A, tol=1.0e-10):
	B = []
	for i in A/sqrt(sum(A*A)):
		for j in B: # Remove all projections.
			i -= j*dot(i,j)
		m = sum(dot(i,i))
		if m > tol: # Linearly independent constraint.
			B.append(i/sqrt(m))
	return array(B)

# Solve Ax=b, where A is lower diagonal.
def forward_subst(A, b):
	x = b.copy()
	for i in range(len(x)):
		x[i] = (x[i]-sum(A[i,:i,newaxis]*x[:i,...],0))/A[i,i]
	return x

# Solve Ax=b, where A is upper diagonal.
def back_subst(A, b):
	x = b.copy()
	for i in reversed(range(len(x))):
		x[i] = (x[i]-sum(A[i,i+1:,newaxis]*x[i+1:,...],0))/A[i,i]
	return x
