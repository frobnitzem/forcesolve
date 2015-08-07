# This file is part of ucgrad, Copyright 2008 David M. Rogers.
#
#   ucgrad is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ucgrad is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ucgrad (i.e. frc_solve/COPYING).
#   If not, contact the author(s) immediately \at/ wantye \/ gmail.com or
#   http://forceSolve.sourceforge.net/. And see http://www.gnu.org/licenses/
#   for a copy of the GNU GPL.

from numpy import *
import numpy.linalg as la

def arg_amin(list):
	min = fabs(list[0][1])
	I = 0
	for i in range(len(list)):
		if(list[i][1] < min):
			I = i
			min = fabs(list[i][1])
	return I

def to_tuple_list(x):
	list = []
	for i in x:
		list.append( (i[0], i[1]) )
	return list

# numerical integration using the triangle rule
def Nintegrate(x, i0=-1, y0=0.0):
	if(type(x) == type(array([0]))):
		x = to_tuple_list(x)
	x.sort()	# assume list of tuples
	
	i0 = int(i0)
	I = zeros(len(x), Float)
	
	# minimize numerical error by starting small
	k = arg_amin(x)
	
	for i in range(k+1, len(x)):
		I[i] = I[i-1] + 0.5*(x[i][0]-x[i-1][0])*(x[i][1]+x[i-1][1])
	
	for i in range(k, 0, -1):
		I[i-1] = I[i] - 0.5*(x[i][0]-x[i-1][0])*(x[i][1]+x[i-1][1])
	
	fix = 0.0
	if(i0 >= 0):
		fix = y0 - I[i0]
	return I + fix

# Solve a system of overdetermined equations,
# returning the sol.n, and its covariance matrix.
# Note: if debugging, make sure U and V are not switched, as this
# may have changed when switching from numarray->numpy!
def solve_system(A, b, tol=1.0e-10):
	U,w,V = la.svd(A, 0) # Don't need full matrices, just main ones.
	
	# divide each row by w (i.e. each column by w_i)
	#V.transpose()
	#V = V/w #-- replaced with following for stability
	low_val = max(w)*tol # too low -- not good!
	out = []
	for i,v in enumerate(w):
		if(v < low_val):
			#print "Parameter %d is unimportant -- %f"%(i+1, v)
			out.append(i)
			V[i] *= 0.0
		else:
			V[i] /= v
	V.transpose()
	
	# Solution
	x = dot(V, dot(transpose(U), b))
	# Covariances
	C = dot(transpose(V), V)
	# Residual
	#res = sum((dot(A, x)-b)**2)
	
	return x, C, out

# Returns the mean, variance, and the 
# "asymptotic normality assumption" error (standard deviation) in the same.
def stats(x):
	#reduce(multiply, a.shape, 1) #total size
	n = float(x.shape[0])
	m = sum(x,0)/n
	v = x-m
	v = sum(v*v,0)/n
	#me = sqrt(v/n))
	me = sqrt(v/(n-1))
	ve = v*2/sqrt(n-1)
	return m,v,me,ve
