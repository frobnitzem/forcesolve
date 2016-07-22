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
from cg_topol import write_topol, show_index

try:
    from quad_prog import quad_prog
except ImportError:
    quad_prog = None

# cg_topol uses integrated first, second and third derivatives for prior which,
# when multiplied by the corresponding values in alpha, serve to push the
# energy function toward a flat, linear, or quadratic shape -- respectively.

class frc_match:
	def __init__(self, topol, pdb, dt, kT, \
			E0=1.0e-8, calpha=1.0e+3, logalpha=None):
		# No inputs defined yet...
		self.topol = topol
                self.pdb = pdb
                self.ind = topol.ind + [topol.params]
		#assert topol.can_force_match()
		
		self.logalpha = logalpha
		self.build_type_index()
		types = self.types
		params = self.topol.params
		
		self.dt = dt
		self.kT = kT
		
		self.D = zeros((types,params))
		self.D2 = zeros((types,params,params))
		self.DF = zeros((types,params))
		self.F2 = zeros((types))
		self.S = 0
		
		self.E0 = E0
		self.calpha = calpha
		#for i in range(len(self.topol.prior)):
		#    if self.prior_rank[i] != len(self.prior[i]):
	#		val, A = la.eigh(self.prior[i])
	#		for j in range(len(self.prior[i])-self.prior_rank[i]):
	#			val[j] = self.E0*self.E0*0.01
	#		self.prior[i] = dot(A*val[newaxis,:], transpose(A))
	#		#A = A[:,-self.prior_rank[i]:]
	#		#val = val[-self.prior_rank[i]:]
	#		#self.prior[i] = dot(A*val[newaxis,:], transpose(A))
		self.alpha = ones(self.topol.hyp_params)
		
		# Normalize.
		self.orthonormalize_constraints(array(self.topol.constraints))
                self.ineqs = array(self.topol.ineqs)
		
		#self.theta = self.theta_from_topol() # Parameters.
                self.theta = zeros(self.topol.params)
		self.z = ones(types)
		
		# Sampling accumulators.
		self.sum_v = []
		self.sum_alpha = []
		self.sum_theta = []
		self.sum_dmu2 = []
		self.sum_df2 = []
		self.df2 = zeros(len(topol.ind))
		self.samples = 0
	
	# Build an index to atoms by designated atom type.
	def build_type_index(self):
                tmass = dict([(aname[2],m) for aname,m in
                                zip(self.pdb.names,
                                    self.pdb.mass)])
		self.type_names = tmass.keys()
		self.types = len(self.type_names)
		#print self.type_names
		
                # Dictionary for quick reference.
		type_num = dict([(t,i) for i,t in enumerate(self.type_names)])
		mass = array([tmass[n] for n in self.type_names])
		
		nt = [0]*self.types
		
		self.type_index = []
		for aname in self.pdb.names:
			tn = type_num[aname[2]]
			nt[tn] += 1
			self.type_index.append(tn)
		self.mass = self.pdb.mass
		self.nt = array(nt)
                #print self.nt, mass, self.type_index

	# Estimate the actual dimensionality of iC.
	def dimensionality(self, sigma_tol=1.0e-5):
		tol = sigma_tol*sigma_tol
		itol = 1.0/tol
		iC = self.calc_iC()
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
	def append(self, x,f):
		chunk = 100
		
		if self.samples > 0:
			raise ValueError, "Error! Cannot add more data "\
				"points to mature frc_match object."
		if x.shape != f.shape:
			print x.shape, f.shape
			raise ValueError, "Error! x and f trajectory shapes "\
				"differ!"
		if x.shape[-2] != self.pdb.atoms:
			raise ValueError, "Error! Number of atoms in "\
				"trajectory does not match topology!"
                if x.shape[-1] != 3:
			raise ValueError, "Error! last dim should be crd xyz!"
		
		# Multiply by ugly constants here.
		Xfac = self.dt*sqrt(self.kT/self.mass) # Non-dimensionalize.
		f *= (self.dt/sqrt(self.mass*self.kT))[newaxis,:,newaxis]
		
		print "Appending %d samples..."%(len(x))
		self.S += len(x)
		# Operate on "chunk" structures at once.
		for i in range(0,len(x)-chunk+1,chunk):
		    D = -1.0*self.topol.design(x[i:i+chunk],1)[1] \
			* Xfac[newaxis,:,newaxis,newaxis] # Factor cancels 1/dx
		    self.D += self.type_sum(sum(sum(D,-2),0))
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:i+chunk])
		    self.F2 += self.type_sum(sum(sum(\
			f[i:i+chunk]*f[i:i+chunk],-1),0))
		i = len(x)%chunk
		if i != 0:
		    i = len(x)-i
		    D = -1.0*self.topol.design(x[i:],1)[1] \
			* Xfac[newaxis,:,newaxis,newaxis] # Factor cancels 1/dx
		    self.D += self.type_sum(sum(sum(D,-2),0))
		    self.type_sum_D2(D)
		    self.type_sum_DF(D, f[i:])
		    self.F2 += self.type_sum(sum(sum(f[i:]*f[i:],-1),0))
		
		# 1 structure at a time.
		#for xi, fi in zip(x,f):# Design matrices are for energy deriv.s
		#	D = -1.0*self.topol.design(xi,1)[1] # -ize -> F.
		#	self.D += self.type_sum(sum(D, 1))
		#	self.D2 += self.type_sum(sum(\
		#		D[:,:,:,newaxis]*D[:,:,newaxis,:],1))
		#	self.DF += self.type_sum(sum(D*fi[:,:,newaxis],1))
		#	self.F2 += self.type_sum(sum(fi*fi,1))
	
	# Given an (atoms by ?) array, reduce to a (types by ?) array
	# by summation.
	def type_sum(self, x):
		xt = zeros((self.types,)+x.shape[1:])
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
                self.constraints = zeros((1,self.topol.params))
            else:
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
		    rhs = dot(self.z, self.DF) # weight residual by atom type

                    if len(self.ineqs) > 0 and quad_prog != None:
			print "        Using constrained solve."
                        # solve with inequality constraints
                        if abs(self.constraints).max() < 1e-10:
			  theta = quad_prog(iC, -rhs, \
                                        G = -self.ineqs, \
                                        h = zeros(len(self.ineqs)))
                        else:
			  theta = quad_prog(iC, -rhs, \
                                        G = -self.ineqs, \
                                        h = zeros(len(self.ineqs)), \
					A = self.constraints,
					b = zeros(len(self.constraints)))
		    else:
                        # TODO: hack by setting violated constraints to zero
			theta = la.solve(iC, rhs)
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
		df2 = array(df2)/( self.S*3.0*self.pdb.atoms )
		self.df2 = df2
	
	def calc_penalty(self):
            pen = zeros(self.topol.hyp_params)
            for i,P in self.topol.prior:
                r0 = self.ind[i]
                r1 = self.ind[i+1]
                pen[i] = dot(dot(P, self.theta[r0:r1]), self.theta[r0:r1])
            return pen*(pen > 0.0) # Forces negative pen -> 0
	
	# Note: const(D)-log(P) = dot(bz,self.z) - dot(az-1,log(self.z))\ 
	#                + dot(ba,self.alpha) - dot(aa-1,log(self.alpha))
	def calc_za_ab(self):
            self.ft2 = dot(dot(self.D2,self.theta), self.theta)
            bz = (self.ft2+self.F2) - 2*dot(self.DF, self.theta)
            bz = 0.5*(bz*(bz > 0.0) + self.E0)
            
            ba = 0.5*(self.calc_penalty()+self.E0)
            return 1.5*self.nt*self.S, bz, 0.5*array(self.topol.pri_rank), ba
		
	def calc_iC(self):
		iC = sum(self.z[:,newaxis,newaxis]*self.D2, 0)
		# Add in constraints.
		iC += self.calpha*dot(transpose(self.constraints), \
					self.constraints)
		# Add in prior info.
                for a, (i,P) in zip(self.alpha, self.topol.prior):
                    r0 = self.ind[i]
                    r1 = self.ind[i+1]
                    iC[r0:r1, r0:r1] += a * P
		return iC
	
	# Generate conditional samples.
	def update_sample(self, n=0, logalpha=None):
	    for i in xrange(n):
		#print "    Updating."
		iC = self.calc_iC()
		b = zeros((self.topol.params,2))
		b[:,0] = rand.standard_normal(self.topol.params) # Sample
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
                for k,i in enumerate(self.ind[:-1]):
			ip = self.ind[k+1]
			D2t = sum(self.D2[:,i:ip,i:ip],0)
			fv.append(trace(dot(C[i:ip,i:ip],D2t)))
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
                        for k,i in enumerate(self.ind[:-1]):
				ip = self.ind[k+1]
				D2t = sum(self.D2[:,i:ip,i:ip],0)
				f2.append(dot(dot(D2t, dmu[i:ip]), dmu[i:ip]))
			df2 += array(f2)
		dmu2 /= self.samples*self.S*3.0*self.nt
		df2 /= self.samples*self.S*3.0*sum(self.nt)
		v = sum(self.sum_v,0)/self.samples + dmu2
		self.df2 = df2
		#print "LINPROB = %e"%(float(sum(array(self.sum_E)<v*1e-3))\
		#	/ self.samples)
		#print "ERFAC = %e"%(self.ab0*2.0/v) # product of z and E0
		self.z = 1.0/v

        def write_out(self, name):
            # Dimensionalize and use topol's own write methods.
            write_topol(self.topol, name, self.theta*self.kT)

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
            
            # cheating...
            def get_names(t):
                if hasattr(t, "terms"):
                    return reduce(lambda a,b: a+get_names(b), t.terms, [])
                return [t.name]
            tname = get_names(self.topol)
            df = open(name+"df.out", 'w')
            df.write("#type\t<stdev>\n")
            for df2, (i,P) in zip(self.df2, self.topol.prior):
                df.write( "%-16s %e\n"%(tname[i], sqrt(df2)) )
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
