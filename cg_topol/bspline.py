# Bspline master classes

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

##################### General Support Routines #############################

def eval_poly1(x, c):
	n = len(c)
	y = c[n-1]
	for i in range(n-2,-1,-1):
		y = c[i]+y*x
	return y
def eval_poly2(x, c):
	n = len(c)
	
	dy = 0.0
	y = c[n-1]
	for i in range(n-2, -1, -1):
		dy = y + dy*x
		y = c[i] + y*x
	
	return y, dy

# N-dimensional derivative form. -- TODO More efficient
def eval_poly(x, c, nd=None):
	if nd == None:
		return eval_poly1(x, c)
	elif nd == 0: # see TODO.
		return reshape(eval_poly1(x, c), (len(x),1))
	elif nd == 1: # see TODO.
		return transpose(array(eval_poly2(x, c)))
	print "Warning! Using eval_poly is unstable!"
	k = array([sum(Darr(i,len(c)),0) for i in range(nd+1)])
	z = ones(x.shape)
	print c,z,k
	s = c[0]*z[:,newaxis]*k[newaxis,:,0]
	for i in range(1, len(c)):
		z *= x
		s += c[i]*z[:,newaxis]*k[newaxis,:,i]
	ix = 1.0/x
	for j in range(1,nd+1):
		s[:,j:] *= ix[:,newaxis]
	return s

# Add A to B at position starting at i,j of A
def wrap_radd_matrix(B, A, i, j):
		# Put the two starting indices into the proper range for A.
		i -= int(A.shape[0]*floor(float(i)/A.shape[0]))
		j -= int(A.shape[1]*floor(float(j)/A.shape[1]))
		ei = i+B.shape[0]
		ej = j+B.shape[1]
		#print i,ei,j,ej
		if ei <= A.shape[0]:
		    if ej <= A.shape[1]:
			B += A[i:ei,j:ej]
		    else:
			B[:,:A.shape[1]-j] += A[i:ei,j:]
			B[:,A.shape[1]-j:] += A[i:ei,:ej-A.shape[1]]
		else:
		    if ej <= A.shape[1]:
			B[:A.shape[0]-i,:] += A[i:,j:ej]
			B[A.shape[0]-i:,:] += A[:ei-A.shape[0],j:ej]
		    else:
			B[:A.shape[0]-i,:A.shape[1]-j] += A[i:,j:]
			B[:A.shape[0]-i,A.shape[1]-j:] += A[i:,:ej-A.shape[1]]
			B[A.shape[0]-i:,:A.shape[1]-j] += A[:ei-A.shape[0],j:]
			B[A.shape[0]-i:,A.shape[1]-j:] += \
				A[:ei-A.shape[0],:ej-A.shape[1]]
# Add B to A at position starting at i,j of A
def wrap_add_matrix(A, B, i, j):
		# Put the two starting indices into the proper range for A.
		i -= int(A.shape[0]*floor(float(i)/A.shape[0]))
		j -= int(A.shape[1]*floor(float(j)/A.shape[1]))
		ei = i+B.shape[0]
		ej = j+B.shape[1]
		#print i,ei,j,ej
		if ei <= A.shape[0]:
		    if ej <= A.shape[1]:
			A[i:ei,j:ej] += B
		    else:
			A[i:ei,j:] += B[:,:A.shape[1]-j]
			A[i:ei,:ej-A.shape[1]] += B[:,A.shape[1]-j:]
		else:
		    if ej <= A.shape[1]:
			A[i:,j:ej] += B[:A.shape[0]-i,:]
			A[:ei-A.shape[0],j:ej] += B[A.shape[0]-i:,:]
		    else:
			A[i:,j:] += B[:A.shape[0]-i,:A.shape[1]-j]
			A[i:,:ej-A.shape[1]] += B[:A.shape[0]-i,A.shape[1]-j:]
			A[:ei-A.shape[0],j:] += B[A.shape[0]-i:,:A.shape[1]-j]
			A[:ei-A.shape[0],:ej-A.shape[1]] += \
					B[A.shape[0]-i:,A.shape[1]-j:]

# Add A to B at position starting at i,j of A
def intersect_radd_matrix(B, A, i, j):
		A0 = [i,j]
		A1 = [i+B.shape[0],j+B.shape[1]]
		B0 = [0,0]
		B1 = [B.shape[0],B.shape[1]]
		# Out of bounds!
		if A1[0] <= 0 or A1[1] <= 0 \
				or i >= A.shape[0] or j >= A.shape[1]:
			print "Warning! intersect_add_matrix called"\
				" without any intersection!"
			return
		if A0[0] < 0:
			B0[0] -= A0[0]
			A0[0] = 0
		if A0[1] < 0:
			B0[1] -= A0[1]
			A0[1] = 0
		if A1[0] > A.shape[0]:
			B1[0] -= A1[0]-A.shape[0]
			A1[0] = A.shape[0]
		if A1[1] > A.shape[1]:
			B1[1] -= A1[1]-A.shape[1]
			A1[1] = A.shape[1]
		
		B[B0[0]:B1[0], B0[1]:B1[1]] += A[A0[0]:A1[0], A0[1]:A1[1]]
# Add B to A at position starting at i,j of A
def intersect_add_matrix(A, B, i, j):
		A0 = [i,j]
		A1 = [i+B.shape[0],j+B.shape[1]]
		B0 = [0,0]
		B1 = [B.shape[0],B.shape[1]]
		# Out of bounds!
		if A1[0] <= 0 or A1[1] <= 0 \
				or i >= A.shape[0] or j >= A.shape[1]:
			print "Warning! intersect_add_matrix called"\
				" without any intersection!"
			return
		if A0[0] < 0:
			B0[0] -= A0[0]
			A0[0] = 0
		if A0[1] < 0:
			B0[1] -= A0[1]
			A0[1] = 0
		if A1[0] > A.shape[0]:
			B1[0] -= A1[0]-A.shape[0]
			A1[0] = A.shape[0]
		if A1[1] > A.shape[1]:
			B1[1] -= A1[1]-A.shape[1]
			A1[1] = A.shape[1]
		
		A[A0[0]:A1[0], A0[1]:A1[1]] += B[B0[0]:B1[0], B0[1]:B1[1]]

######################### Integration support routines ####################
def factorial(n):
	if n <= 1.5:
		return 1
	
	return n*factorial(n-1)

# Derivative array for nth order derivative. d(n)P/dx**n = [x^i]^T*D(n)*C
def Darr(n, order):
	D = zeros((order, order))
	i = 0
	D[i,n] = factorial(n)
	for j in range(n+1,order):
		i += 1
		D[i,j] = D[i-1,j-1]*j/(j-n)
	return D

# Binomial formula coefficients for powers of x in (s*x+t)**I
# B_IJ = s**J t**(I-J) (I choose J)
def Barr(t, order, s=1.0):
	B = identity(order)
	for i in range(1, order):
		B[i,i] = B[i-1,i-1]*s
		B[i,0] = B[i-1,0]*t
		ip = i
		for j in range(1,order-i):
			ip = ip+1
			B[ip,j] = B[ip-1,j-1]*s*ip/j
	return B

######################## Integration Routines ##############################
def sum_elem(i, j, k):
	return i+j+k

# Returns a matrix corresponding to the integral from a to b of
# x**(i+j+k)
def poly_integral3(shape, a, b):
	I = zeros(shape[0:3]) # powers of b
	A = zeros(shape[0:3]) # powers of a
	I[0,0,0] = b
	A[0,0,0] = a
	for i in range(1, shape[0]):
		I[i,0,0] = I[i-1,0,0]*b
		A[i,0,0] = A[i-1,0,0]*a
	for j in range(1, shape[1]):
		I[:,j,0] = I[:,j-1,0]*b
		A[:,j,0] = A[:,j-1,0]*a
	for k in range(1, shape[2]):
		I[:,:,k] = I[:,:,k-1]*b
		A[:,:,k] = A[:,:,k-1]*a
	I -= A
	del A
	
	return I / (fromfunction(sum_elem, shape[0:3])+1.0)

# Returns a matrix corresponding to the integral from a to b of
# x**(i+j)*(s*x+t)**m
def poly_integral(order, a, b, m=0, t=0, s=1.0):
	if m == 0: # Problem with rank of matrix returned from Barr...
		return poly_integral3((1,order,order), a, b)[0]
	else:
		return dot(Barr(t, m+1, s)[m], transpose(\
			poly_integral3((m+1, order, order), a, b), (1,0,2)))

def calc_carr(n):
	C = zeros((n/2+1, n))
	
	c = zeros(n-1)
	C[0,n-1] = 1.0
	for i in range(1, len(C)):
		C[i,n-1] = C[i-1,n-1]*(i-n-1)/i
	
	rng = arange(len(C))
	for j in range(n-1, 0, -1):
		C[:,j-1] = C[:,j]*rng*j/(j-n)
	for i in range(1, len(C)):
		C[i] += C[i-1]
	C /= factorial(n-1)
	
	return C

################### Top-Level Spline Class Interfaces ###############
# Bspline functions for a given spline order.
class Bspline:
	def __init__(self, order):
		self.order = order
		self.odd = order%2
		
		self.C = calc_carr(order)
		self.M = self.calc_polycoef()
	
# This powerful routine constructs the matrix which, when left-multiplied by
# spline coefficients at knots int(u)-(order-1), ..., int(u)
# and right-multiplied by [x^i] gives the value of the spline at point x.
	def calc_polycoef(self):
		M = zeros((self.order, self.order))
		for i in range(self.order/2):
			M[i] = dot(self.C[i], Barr(i+1, self.order, -1.0))
			M[self.order-1-i] = dot(self.C[i], Barr(i, self.order))
		if self.odd:
			i = self.order/2
			M[i] = dot(self.C[i], Barr(i+1, self.order, -1.0))
		
		return M

# This routine assumes that x=ceil(s)-s is in [0,1) and calculates all
# (spline order) of the shifted spline coefficients.  These should be placed
# in positions ceil(s)-order to ceil(s) and wrapped if the spline is periodic.
# If nd is defined, it calculates the Pxnd+1 array (for derivatives 0,...,nd)
	def calc_splcoef(self, x, nd=None):
		if nd == None:
		    a = zeros((len(x), self.order))
		    for i in range(self.order/2):
			a[:,i] = eval_poly(i+x, self.C[i])
			a[:,self.order-1-i] = eval_poly(i+1-x, self.C[i])
		    if self.odd:
			i = self.order/2
			a[:,i] = eval_poly(i+x, self.C[i])
		else:
		    a = zeros((len(x), self.order, nd+1))
		    for i in range(self.order/2):
			a[:,i] = eval_poly(i+x, self.C[i], nd)
			a[:,self.order-1-i] = eval_poly(i+1-x, self.C[i], nd)
		    if self.odd:
			i = self.order/2
			a[:,i] = eval_poly(i+x, self.C[i], nd)
		    for d in range(1,nd+1,2): # Negate odd dimensions
						#  (chain rule)
			a[:,:self.order/2+self.order%2,d] *= -1.0
			
		
		return a

# Calculates the matrix corresponding to the integral of the nth derivative 
# of polynomial 1 times nth derivative of polynomial 2 times (s*x+t)**m
# -- Useful for calculating the quadratic integral matrix.
	def dpoly_integral(self, a, b, n=0, m=0, t=0, s=1.0):
		if n >= self.order:
			raise RuntimeError, "%d order spline is only %d times "\
				"differentiable!"%(self.order,self.order-1)
		D = Darr(n, self.order)
		
		I = dot(dot(transpose(D), \
				poly_integral(self.order, a, b, m, t, s)), D)
		
		return I

# A Cardinal Bsplined function. Bin specification is (x0,h,number of knots)
# Periodic Bsplines have L=n*h
# Aperiodic should have shift=order-1, L=(n-shift)*h, since this gives
# full order polynomials on each interval.
# Spl is a spline class.
class spline_func:
	def __init__(self, bin, spl, periodic=False, shift=0, rng=(-inf,inf), id="", verb=True):
		self.spl = spl
		self.id = id
		self.x0 = bin[0]
		self.h = bin[1] # Height.
		self.ih = 1.0/bin[1]
		self.n = bin[2] # Number of points.
		self.c = zeros(self.n) # Spline constants.
		if self.n < spl.order:
			raise "InputError","Error! Number of spline points is "\
				"less than spline order!"
		
		self.periodic = periodic # Periodicity flag.
		# Shift: order/2 for periodic and order-1 o.w.
		self.shift = shift
		
		# Must be intersected with actual nonzero range for non-periodic
		if not periodic: # splines to function correctly.
			assert rng[0] < rng[1]
			if (rng[0]-self.x0)*self.ih+shift > self.spl.order-1:
			    raise "InputError", "Error! Lower range of spline "\
					"creates meaningless coeff.s"
			if (rng[1]-self.x0)*self.ih+shift < self.n-1:
			    raise "InputError", "Error! Upper range of spline "\
					"creates meaningless coeff.s"
			self.rng = [0,0]
			self.rng[0] = max([rng[0], self.x0-self.h*shift])
			self.rng[1] = min([rng[1], \
				self.x0+self.h*(self.n-1+self.spl.order-shift)])
		else: # Indicates unique interval range.
			self.rng = [self.x0, self.x0+self.h*self.n]
		self.verb = verb # Verbose feedback?
	
	# Pre-evaluate all spline intervals to get polynomial coeff.s
	# eliminating storage of a large matrix and a multiplication
	# step upon spline function calculation.
	def commit(self):
		pcoef = []
		lpcoef = []
		r = self.spl.order
		cm = reshape(self.c, (1,len(self.c)))
		
		if self.periodic: # Use add_matrix for convenience.
			add_matrix = wrap_radd_matrix
		else:
			add_matrix = intersect_radd_matrix
		urng = ((self.rng[0]-self.x0)*self.ih+self.shift,\
			(self.rng[1]-self.x0)*self.ih+self.shift)
		
		# First partial-interval?
		if (floor(urng[0]+0.5)-urng[0])**2 > 1.0e-20:
			u0 = int(urng[0])+1 # Since we know urng0 is non-int
			c = zeros((1,r))
			add_matrix(c, cm, 0, u0-r)
			pcoef.append(dot(c,self.spl.M))
			self.pcoef_shift = -(u0-1)
		else: # Begins on integer.
			u0 = int(urng[0]+0.5)
			self.pcoef_shift = -u0
		
		# Last partial-interval?
		if (floor(urng[1]+0.5)-urng[1])**2 > 1.0e-20:
			u1 = int(urng[1])
			c = zeros((1,r))
			add_matrix(c, cm, 0, u1+1-r)
			lpcoef.append(dot(c,self.spl.M))
		else:
			u1 = int(urng[1]+0.5)
		
		# Everything in-between.
		for i in range(u0,u1):
			c = zeros((1,r))
			add_matrix(c, cm, 0, i+1-r)
			pcoef.append(dot(c,self.spl.M))
		
		self.pcoef = reshape(array(pcoef + lpcoef), (-1,r))
	
	def ay(self, x, nd=None):
		# Calculate range for polynomial coeff. lookup.
		if self.periodic: # Wrap out-of range points back in range.
			L = self.rng[1]-self.rng[0]
			x -= L*floor((x-self.rng[0])/L)
		u = (x-self.x0)*self.ih+self.shift
		if nd != None:
			y = zeros((nd+1,len(x)))
		else:
			y = zeros(len(x))
		for i,w in enumerate(u):
			if x[i] < self.rng[0] or x[i] >= self.rng[1]:
			    if self.verb:
				print "Warning! Out-of range point in %s: %f"%(\
						self.id, x[i])
			    continue
			y[...,i] = eval_poly(w-floor(w), \
				self.pcoef[int(w)+self.pcoef_shift], nd)
		if nd != None: # Multiply by appropriate Jacobian.
			for i in range(1,len(y)):
				y[i:] *= self.ih
		return y
	
	# Calculates deriviatives 0,...,n of the function.
	def y(self, x, nd=None):
		return dot(self.spline(x, nd), self.c)
	
	# Used to construct vectors which multiply parameters.
	# If x is a N-dim vector, the return value is an nd+1xNxP matrix
	def spline(self, x, nd=None):
		if(len(x.shape) != 1):
			raise ValueError, "Spline passed != 1 dim. array!"
		
		s = (x-self.x0)*self.ih+self.shift
		u = ceil(s) # floaing point
		
		# Calculate all nonzero spline coeff.s up to derivative nd.
		if nd != None:
			M = self.spl.calc_splcoef(u-s, nd)
			Mp = zeros((len(x), self.n, nd+1))
		else:
			M = self.spl.calc_splcoef(u-s)
			Mp = zeros((len(x), self.n))
		if self.periodic: # Wrap out-of range points back in range.
			u -= self.n*floor((u-self.spl.order)/self.n)
			u = u.astype(int)
			ui = [(u, i) for i,u in enumerate(u)]
			ui.sort()
			n = 0
			while n < len(ui) and ui[n][0] <= self.n:
				j,i = ui[n]
				Mp[i,j-self.spl.order:j] = M[i]
				n += 1
			while n < len(ui):
				j,i = ui[n]
				k = self.n-j+self.spl.order
				Mp[i,j-self.spl.order:] = M[i,:k]
				Mp[i,:j-self.n] = M[i,k:]
				n += 1
		else: # Discard out-of range points.
			u = u.astype(int)
			ui = [(u, i) for i,u in enumerate(u)]
			ui.sort()
			n = 0
			while n < len(ui) and x[ui[n][1]] < self.rng[0]:
				n += 1
			if self.verb and n > 0:
				print "Warning! %d points are smaller "\
					"than min of spline range for "\
					"type id %s"%(n, self.id)
			while n < len(ui) and ui[n][0] < self.spl.order:
				j,i = ui[n]
				Mp[i,:j] = M[i,self.spl.order-j:]
				n += 1
			while n < len(ui) and ui[n][0] < self.n:
				j,i = ui[n]
				Mp[i,j-self.spl.order:j] = M[i]
				n += 1
			while n < len(ui) and \
				    x[ui[n][1]] <= self.rng[1]:
				j,i = ui[n]
				k = self.n-j+self.spl.order
				Mp[i,j-self.spl.order:] = M[i,:k]
				n += 1
			if self.verb and n != len(ui):
				print "Warning! %d points are larger than max "\
					"of spline range for type id %s"%(\
						len(ui)-n, self.id)
				#print ui[n:]
		
		if nd == None:
			return Mp
		# Make derivative index 1st dimension.
		Mp = transpose(Mp, (2,0,1))
		
		# Multiply by appropriate powers of du/dx
		ih = self.ih
		for i in range(1,nd+1):
			Mp[i] *= ih
			ih *= self.ih
		return Mp
	
	# Returns the vector to be multiplied by the spline
	# coefficients, corresponding to integration over y^{(n)}(x) * x**m.
	def spl_integral(self, x, m=0):
		if(len(x.shape) != 1):
			raise ValueError, "Spline passed != 1 dim. array!"
		
		s = (x-self.x0)*self.ih+self.shift
		u = ceil(s) # floaing point
		
		Mi = self.spl.calc_splcoef(u-s)
		Mp = zeros((len(x), self.n))
		
		n = 0
		M = self.spl.M
		r = self.spl.order
		
		Q = zeros((1, self.n))
		if self.periodic:
			add_matrix = wrap_add_matrix
			u -= self.n*floor((u-self.spl.order)/self.n)
		else:
			add_matrix = intersect_add_matrix
		u = u.astype(int)
		ui = [(u, i) for i,u in enumerate(u)]
		ui.sort()
		urng = ((self.rng[0]-self.x0)*self.ih+self.shift,\
			(self.rng[1]-self.x0)*self.ih+self.shift)
		
		k = 0
		while k < len(ui) and s[ui[k][1]] < urng[0]:
			k += 1
		if self.verb and k > 0:
			print "Warning! %d points are smaller "\
				"than min of spline range for "\
				"type id %s"%(k, self.id)
		
		# First partial-interval? - ends on u0
		if (floor(urng[0]+0.5)-urng[0])**2 > 1.0e-20:
			u0 = int(urng[0])+1 # Since we know urng0 is non-int
			while k < len(ui) and ui[k][0] == u0:
			  j,i = ui[k]
			  I = self.spl.dpoly_integral(urng[0]-floor(urng[0]),\
				s[i]-floor(urng[0]),\
				n, m, self.x0*self.ih+u0-1-self.shift)[n:n+1]
			  I = dot(I, transpose(Mi[i]))
			  #Mp[i] = Q.copy()
			  add_matrix(Mp[i], I, 0, j-r)
			  
			I = self.spl.dpoly_integral(urng[0]-floor(urng[0]), 1,\
				n, m, self.x0*self.ih+u0-1-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. array, put it into Q.
			add_matrix(Q, I, 0, u0-r)
			#print "Range: " + str((urng[0],u0))
		else: # Begins on integer.
			u0 = int(urng[0]+0.5)
		
		if (floor(urng[1]+0.5)-urng[1])**2 > 1.0e-20:
			u1 = int(urng[1])
		else:
			u1 = int(urng[1]+0.5)
		
		# Everything in-between. - start at i
		for i in range(u0,u1):
			while k < len(ui) and ui[k][0] == i+1:
			  j,l = ui[k]
			  I = self.spl.dpoly_integral(0, s[l]-i,
				n, m, self.x0*self.ih+i-self.shift)[n:n+1]
			  I = dot(I, transpose(Mi[l]))
			  Mp[l] = Q.copy()
			  add_matrix(Mp[l], I, 0, i+1-r)
			I = self.spl.dpoly_integral(0, 1, n, m, \
					self.x0*self.ih+i-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, 0, i+1-r)
			#print "Range: " + str((i,i+1))
		
		# Last partial-interval? - starts at u1
		if (floor(urng[1]+0.5)-urng[1])**2 > 1.0e-20:
			while k < len(ui) and ui[k][0] == u1+1 and s[ui[k][1]] < urng[1]:
			  j,l = ui[k]
			  I = self.spl.dpoly_integral(0, s[l]-u1,
				n, m, self.x0*self.ih+u1-self.shift)[n:n+1]
			  I = dot(I, transpose(Mi[l]))
			  Mp[l] = Q.copy()
			  add_matrix(Mp[l], I, 0, u1+1-r)
			I = self.spl.dpoly_integral(0, urng[1]-u1, n, m, \
				self.x0*self.ih+u1-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, 0, u1+1-r)
			#print "Range: " + str((u1,urng[1]))
		
		# Must multiply by scale factor, since u=(x-x0)/h
		return Mp*self.h**(m+1-n)
	
	# Returns the vector to be multiplied by the spline
	# coefficients, corresponding to integration over y^{(n)}(x) * x**m.
	def integral(self, n=0, m=0):
		if n != 0:
			raise RuntimeError, "Error! nonzero values of n not yet supported."
		M = self.spl.M
		r = self.spl.order
		
		Q = zeros((1, self.n))
		if self.periodic:
			add_matrix = wrap_add_matrix
		else:
			add_matrix = intersect_add_matrix
		urng = ((self.rng[0]-self.x0)*self.ih+self.shift,\
			(self.rng[1]-self.x0)*self.ih+self.shift)
		
		# First partial-interval?
		if (floor(urng[0]+0.5)-urng[0])**2 > 1.0e-20:
			u0 = int(urng[0])+1 # Since we know urng0 is non-int
			I = self.spl.dpoly_integral(urng[0]-floor(urng[0]), 1,\
				n, m, self.x0*self.ih+u0-1-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. array, put it into Q.
			add_matrix(Q, I, 0, u0-r)
			#print "Range: " + str((urng[0],u0))
		else: # Begins on integer.
			u0 = int(urng[0]+0.5)
		
		# Last partial-interval?
		if (floor(urng[1]+0.5)-urng[1])**2 > 1.0e-20:
			u1 = int(urng[1])
			I = self.spl.dpoly_integral(0, urng[1]-u1, n, m, \
				self.x0*self.ih+u1-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, 0, u1+1-r)
			#print "Range: " + str((u1,urng[1]))
		else:
			u1 = int(urng[1]+0.5)
		
		# Everything in-between.
		for i in range(u0,u1):
			I = self.spl.dpoly_integral(0, 1, n, m, \
					self.x0*self.ih+i-self.shift)[n:n+1]
			I = dot(I, transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, 0, i+1-r)
			#print "Range: " + str((i,i+1))
		
		# Must multiply by scale factor, since u=(x-x0)/h
		return Q*self.h**(m+1-n)
	
	# Returns the quadratic form to be multiplied by the spline
	# coefficients, corresponding to integration over y^{(n)}(x)^2 * x**m.
	def quadratic_integral(self, n=0, m=0):
		M = self.spl.M
		r = self.spl.order
		
		Q = zeros((self.n, self.n))
		if self.periodic:
			add_matrix = wrap_add_matrix
		else:
			add_matrix = intersect_add_matrix
		urng = ((self.rng[0]-self.x0)*self.ih+self.shift,\
			(self.rng[1]-self.x0)*self.ih+self.shift)
		
		# First partial-interval?
		if (floor(urng[0]+0.5)-urng[0])**2 > 1.0e-20:
			u0 = int(urng[0])+1 # Since we know urng0 is non-int
			I = self.spl.dpoly_integral(urng[0]-floor(urng[0]), 1,\
				n, m, self.x0*self.ih+u0-1-self.shift)
			I = dot(dot(M, I), transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, u0-r, u0-r)
			#print "Range: " + str((urng[0],u0))
		else: # Begins on integer.
			u0 = int(urng[0]+0.5)
		
		# Last partial-interval?
		if (floor(urng[1]+0.5)-urng[1])**2 > 1.0e-20:
			u1 = int(urng[1])
			I = self.spl.dpoly_integral(0, urng[1]-u1, n, m, \
				self.x0*self.ih+u1-self.shift)
			I = dot(dot(M, I), transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, u1+1-r, u1+1-r)
			#print "Range: " + str((u1,urng[1]))
		else:
			u1 = int(urng[1]+0.5)
		
		# Everything in-between.
		for i in range(u0,u1):
			I = self.spl.dpoly_integral(0, 1, n, m, \
						self.x0*self.ih+i-self.shift)
			I = dot(dot(M, I), transpose(M))
			# Now I is the coeff. matrix, put it into Q.
			add_matrix(Q, I, i+1-r, i+1-r)
			#print "Range: " + str((i,i+1))
		
		# Must multiply by scale factor, since u=(x-x0)/h
		return Q*self.h**(m+1-2*n)
	
	def write_spl(self, name, mode='w'):
		out = open(name, mode)
		if len(self.id) < 1:
			id = "file"
		else:
			id = self.id.replace(' ', '-').replace(\
				'\t','_').replace('\n', '')
		hdr = "#SPLINE %s order %d "%(id,self.spl.order) # File header.
		if self.periodic:
			hdr += "periodic. "
		else:
			hdr += "aperiodic. "
		hdr += "x0= %f h= %f n= %d s= %f a= %f b= %f\n"%(\
		  self.x0, self.h, self.n, self.shift, self.rng[0], self.rng[1])
		out.write(hdr)
		
		# x(u+r/2) = x0+h*(u+r/2-shift) -- maxima of hat functions.
		x = self.x0+self.h*(self.spl.order*0.5-self.shift)
		for c in self.c:
			out.write("%8.3f %12e\n"%(x, c))
			x += self.h
	        out.close()

def read_spline_func(file, name=None):
	digits = "01234567890.-+"
	
	lines = open(file).readlines()
	info = lines[0][1:].split()
	#print info
	if info[0] != "SPLINE" or len(info) != 17:
		raise InputError, "Error! Input file is not SPLINE type."
	if name == None:
		id = info[1]
	else:
		id = name
	
	order = int(info[3])
	periodic = not "aperiodic" in info[4]
	x0 = float(info[6])
	h = float(info[8])
	n = int(info[10])
	s = float(info[12])
	a = float(info[14])
	b = float(info[16])
	f = spline_func((x0,h,n), Bspline(order), periodic, s, (a,b), id)
	i = 0 # Parse out the spline coeff.s
	for line in lines[1:]:
		tok = line.split()
		if len(tok) > 1:
		   if(tok[0][0] in digits and tok[1][0] in digits):
		      if i>=n:
			print "Warning! Extra data in SPLINE file "\
				"ignored!!"
			break
		      f.c[i] = float(tok[1])
		      i += 1
	if i != n:
		print "Warning! Not all spline coefficients read from %s!"%file
	return f
