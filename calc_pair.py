#!/usr/bin/python

# This code assumes the system has only one type of interaction -- pairwise energy between two atoms of the same type.  It fits this energy for each frame (using least-squares fitting) to give an indication of what the distribution of empirical forces being sent to the FM code is.

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
# 1/17/2008
# This work was supported by a DOE CSGF.

import sys
from cg_topol.ucgrad import *
from cg_topol import *

def main(argv):
	if len(argv) < 4:
		print "Usage: %s <traj base> <topology> <out>"%(argv[0])
		return 1
	
	t = cg_topol(argv[2])
	x = read_matrix(argv[1]+".x")
	x = reshape(x, (-1,t.atoms,3))
	f = read_matrix(argv[1]+".f")
	f = reshape(f, (-1,t.atoms,3))
	
	plist = t.pair[0]
	info = t.pair_info[0]
	delta = array([x[...,j,:]-x[...,i,:] for i,j in plist])
	delta -= t.box_L * floor(delta/t.box_L+0.5)
	b,db = dbond(delta)
	
	S = len(b) # Number of samples
	N = t.atoms
	F = zeros((S,N,2))
	
	for i in range(S): # Calculate least-squares matches to force coeff.s
		calc_coef(b[i],db[i],f[i], plist, N, F[i])
	write_array(argv[3], F)
	
	return 0

def calc_coef(b, db, f, plist, N, F): # Assuming b is (N_p) and db is (N_p, 3)
	Q = N # Must be <= 3N
	y, z = quantile(b, Q, 16.0) # All quantiles with distances < max.
	# Quantial medians:
	F[:,0] = array([y[int((i+0.5)*len(y)/N+0.49999)][0] for i in range(N)])
	
	M = zeros((N,3,Q))
	q = 0
	for d,k in y:
		if z[q+1] < d:
			q += 1
		i,j = plist[k]
		M[i,:,q] -= db[k]
		M[j,:,q] += db[k]
	M = reshape(M, (3*N,Q))
	U,w,V = la.svd(M,0,1)# Quick SVD-based solution of linear-least squares.
	low_val = max(w)*1.0e-8
	n = 0
	for i,v in enumerate(w):
		if(v < low_val):
			V[i] *= 0.0
			n += 1
		else:
			V[i] /= v
	if n:
		print "%d singular values zeroed."%n
	V = transpose(V)
	F[:,1] = dot(V, dot(transpose(U), reshape(f,-1)))

def quantile(x, N, m):
	z = zeros(N+1)
	y = [(xi,i) for i,xi in enumerate(x) if xi < m]
	y.sort()
	
	h = float(len(y))/N
	for i in range(N):
		z[i] = y[int(i*h+0.5)][0]
	z[N] = y[-1][0]+abs(y[-1][0])*.000001
	return y, z

if __name__=="__main__":
	main(sys.argv)
