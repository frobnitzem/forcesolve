# General B-spline sampling class, effectively frc_match for 1D problems.

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

from numpy import *
import numpy.linalg as la
import numpy.random as rand
from .ucgrad import write_matrix, write_list, stats

# Looking back, the only confusing thing about this implementation
# is that I call the paper's spline tension (lambda) as alpha in
# this code.  The paper's ratio (alpha) is nowhere calculated in here.

class pspline_match:
	def __init__(self, func, deriv=1, m=0, aa0=1, ab0=.1, logalpha=None, rewt=False):
		# No inputs defined yet...
		self.f = func
		
		params = self.f.n
		self.logalpha = logalpha
		
		self.D = zeros(params)
		self.D2 = zeros((params,params))
		self.DF = zeros(params)
		self.F2 = 0.0
		self.S = 0
		
		self.aa0 = aa0 # Prior for theta's prior hyperparameters.
		self.ab0 = ab0
		self.rewt = rewt
		
		w, A = la.eigh( self.f.quadratic_integral(deriv,m) \
			/ (self.f.rng[1]**(m+1)-self.f.rng[0]**(m+1)) )
		#w, A = la.eigh(self.f.quadratic_integral(deriv,m) \
		#		* self.f.h**(2*deriv-1))
		#w, A = la.eigh(self.f.quadratic_integral(deriv,m) \
		#		* self.f.h**(2*deriv) * (m+1.) \
		#	/ (self.f.rng[1]**(m+1)-self.f.rng[0]**(m+1))
		for i in range(len(w)): # Flush to numerical zero.
			if w[i] < 1.0e-10:
				#w[i] = 0.0
				w[i] = 1.0e-10
		self.prior = dot(A, w[:,newaxis]*transpose(A))
		self.prior_rank = params-deriv
		
		self.theta = self.f.c
		self.z = 1.0
		self.alpha = 1.0
		
		# Sampling accumulators.
		self.sum_v = []
		self.sum_alpha = []
		self.sum_theta = []
		self.sum_dmu2 = []
		self.sum_wt = []
		self.sum_E = []
		self.samples = 0
	
	# Estimate the actual dimensionality of iC.
	def dimensionality(self, sigma_tol=1.0e-5):
		tol = sigma_tol*sigma_tol
		itol = 1.0/tol
		iC = self.calc_iC()
		w, A = la.eigh(iC)
		free_dim = len([l for l in w if l <= tol])
		fixed_dim = len([l for l in w if l >= itol])
		print "Alpha = %e, z = %e"%(self.alpha, self.z)
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

# Append a set of data points to the present matching object.
# Takes as input x,f(x) samples.
	def append(self, x, fin):
		chunk = 100
		
		if x.shape != fin.shape:
			print x.shape, f.shape
			raise runtimeError, "Error! x and f(x) shapes "\
				"differ!"
		
		# Standardizing the responses:
		#self.scale = sqrt(sum(fin*fin,0)/len(fin))
		self.scale = 1.0
		#print "Dividing responses by "+str(self.scale)
		f = fin/self.scale
		
		print "Appending %d samples..."%(len(x))
		self.S += len(x)
		# Operate on "chunk" structures at once.
		for i in range(0,len(x)-chunk+1,chunk):
		    D = self.f.spline(x[i:i+chunk])
		    self.D += sum(D,0)
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:i+chunk])
		    self.F2 += sum(f[i:i+chunk]*f[i:i+chunk],0)
		i = len(x)%chunk
		if i != 0:
		    i = len(x)-i
		    D = self.f.spline(x[i:])
		    self.D += sum(D,0)
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:])
		    self.F2 += sum(f[i:]*f[i:],0)
		
	# Special type_sum for accumulating D2 matrix.
	def type_sum_D2(self, DS):
		for D in DS:
			self.D2 += D[:,newaxis]*D[newaxis,:]
	# Special type_sum for accumulating DF vector.
	def type_sum_DF(self, D, F):
		self.DF += dot(F, D)
	
	# Find maximum likelihood estimate.
	def maximize(self, tol=array([1.0e-5,1.0e-5]), maxiter=1000):
		if self.S < 1:
			raise ProgramError, "Error! no data has been collected!"
		
		print "Maximizing posterior PDF..."
		tol2 = tol*tol
		delta = tol2 + 1.0
	        iter = 0
		while any(delta > tol2) and iter < maxiter:
		    iter += 1
		    iC = self.calc_iC()
		    theta = self.z*self.DF
		    theta = la.solve(iC, theta)
		    delta[0] = sum((self.theta-theta)**2)/len(theta)
		    self.theta = theta
		    
		    az, bz, aa, ba = self.calc_za_abab()
		    
		    z = (az-1)/bz
		    #z = min((az-1)/bz, 1.0e+18) # Cap the max. z.
		    delta[1] = (z-self.z)**2
		    self.z = z
		    self.alpha = (aa-1)/ba
		    
		    print "  Iteration %d, delta = "%iter + str(sqrt(delta))
	
	def calc_za_abab(self):
		ba = dot(dot(self.prior,self.theta),self.theta)
		ba = 0.5*max(ba, 0.0)+self.ab0
		
		self.ft2 = dot(dot(self.D2, self.theta), self.theta)
		bz = 0.5*(self.ft2+self.F2) - dot(self.DF, self.theta) + 1e-8
		return 0.5*self.S, bz, 0.5*self.prior_rank+self.aa0, ba
	
	def calc_iC(self):
		return self.z*self.D2 + self.alpha*self.prior
	
	# Generate conditional samples.
	def update_sample(self, n=0, logalpha=None):
	    for i in range(n):
		#print "    Updating."
		iC = self.calc_iC()
		b = zeros((self.f.n,2))
		b[:,0] = rand.standard_normal(self.f.n) # Sample
		b[:,1] = self.z*self.DF # Mean
		try:
			L = la.cholesky(iC)
			b[:,1:] = forward_subst(L, b[:,1:])
			b = back_subst(transpose(L), b)
		except la.linalg.LinAlgError:
			w, A = self.dimensionality()
			print "\nForce design matrix is degenerate!"
			#self.sum_dmu2.append(1.0e+10)
			hC = dot(A, transpose(A)/sqrt(w)[:,newaxis])
			b = dot(hC, b)
			b[:,1] = dot(hC, b[:,1])
			#import sys
			#sys.exit()
		#self.mean = b[:,1]
		self.theta = b[:,0]+b[:,1]
		
		az, ibz, aa, iba = self.calc_za_abab()
		
		self.z = rand.gamma(az, 1.0)/ibz
		self.alpha = rand.gamma(aa, 1.0)/iba
		
		if logalpha != None:
		    write_matrix(logalpha, array([[self.alpha, self.z]]),'a')
	
	def calc_theta_stats(self):
		iC = self.calc_iC()
		theta = la.solve(iC, self.z*self.DF)
		return theta, trace(la.solve(iC, self.D2))
	def sample(self, samples, skip=100, toss=10):
		if self.S < 1:
			raise ProgramError, "Error! no data has been collected!"
		if skip > 0:
			print "Doing sampling burn-in..."
		for i in range(skip):
			self.update_sample(toss, None)
		print "Collecting %d samples..."%samples
		if self.logalpha != None:
			out = open(self.logalpha, 'w')
			out.truncate()
			out.close()
		for i in range(samples):
			self.update_sample(toss, self.logalpha)
			theta, dmu2 = self.calc_theta_stats()
			self.sum_theta.append(theta)
			self.sum_dmu2.append(dmu2)
			self.sum_v.append(1./self.z)
			self.sum_alpha.append(self.alpha)
			self.sum_E.append(dot(dot(self.prior,theta),theta))
			self.sum_wt.append(self.ab0*self.alpha)
			self.samples += 1
		if self.rewt:
			self.rewt_posterior_estimate()
		else:
			self.posterior_estimate()
	
	# Best estimates from posterior distribution by the criterion
	# of minimum information loss.
	def posterior_estimate(self):
		if self.samples < 1:
			raise ProgramError, "Error! No samples collected!"
		self.theta = sum(self.sum_theta,0)/self.samples
		dmu2 = sum(self.sum_dmu2)
		for dmu in array(self.sum_theta)-self.theta:
			dmu2 += dot(dot(self.D2, dmu), dmu)
		dmu2 /= self.samples*self.S
		v = sum(self.sum_v)/self.samples + dmu2
		print "LINPROB = %e"%(float(sum(array(self.sum_E)<v*1e-3))\
					/ self.samples)
		print "ERFAC = %e"%(self.ab0*2.0/v) # product of z and E0
		self.z = 1.0/v
		self.alpha = sum(self.sum_alpha)/self.samples
	
	def rewt_posterior_estimate(self):
		if self.samples < 1:
			raise ProgramError, "Error! No samples collected!"
		wt = exp(self.sum_wt - max(self.sum_wt)) # <= 1.0 ea.
		fac = 1./sum(wt)
		
		self.theta = sum(self.sum_theta*wt[:,newaxis],0)*fac
		dmu2 = sum(self.sum_dmu2*wt)
		for dmu, w in zip(array(self.sum_theta)-self.theta, wt):
			dmu2 += dot(dot(self.D2, dmu), dmu)*w
		dmu2 /= self.S
		v = (sum(self.sum_v*wt) + dmu2)*fac
		print "LINPROB = %e"%(sum((array(self.sum_E)<v*1e-3)*wt)*fac)
		print "ERFAC = %e"%(self.ab0*2.0/v) # product of iz and E0
		self.z = 1.0/v
	
	def show_index(self):
		print "%d total parameters."%self.f.n
		print
	
	def write_out(self, name):
		#print "Re-scaling resulting spline parameters by %e"%self.scale
		self.f.c = self.theta*self.scale
		self.f.write_spl(name+".espl")
		
		if self.rewt:
			wt = exp(self.sum_wt - max(self.sum_wt)) # <= 1.0 ea.
			fac = 1./sum(wt)
		
		if self.samples > 0:
			if self.rewt:
			  avg_la = sum(log(self.sum_alpha)*wt)*fac
			  s_la = sqrt(sum(wt*(log(self.sum_alpha)-avg_la)**2)*fac)
			else:
			  m, v, me, ve = stats(log(self.sum_alpha))
			  avg_la = m
			  s_la = sqrt(v)
		else:
			avg_la = log(self.alpha)
			s_la = 0.0
		
		if self.samples > 0:
			if self.rewt:
			  avg_v = sum(self.sum_v*wt)*fac
			  s_v = sqrt(sum(wt*(self.sum_v-avg_v)**2)*fac)
			else:
			  m, v, me, ve = stats(array(self.sum_v))
			  avg_v = m
			  s_v = sqrt(v)
		else:
			avg_v = 1.0/self.z
			s_v = 0.0
		avg_v *= self.scale**2
		s_v *= self.scale**2
		lam = open(name+".zla", 'w')
		lam.write("#type\tval\t<val>\terr\n")
		lam.write("v   %e\t%e\t%e\n"%(self.scale**2/self.z,avg_v,s_v))
		lam.write("la  %e\t%e\t%e\n"%(log(self.alpha), avg_la, s_la))
		lam.close()

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
