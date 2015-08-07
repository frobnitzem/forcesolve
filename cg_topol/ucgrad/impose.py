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

from math import sqrt
from numpy import *
import numpy.linalg as la

TOL = 1.0e-6
MAXCYC = 25

class impose:
	def __init__(self, in_coords, weights, g_align):
		self.nstructs = len(in_coords)
		weights /= sum(weights)		# normalize weights
		if self.nstructs < 2:
			raise ValueError, "Cannot calculate superposition on" \
						"single/no structure."
		self.natoms = len(weights)
		msd_structs = self.nstructs*(self.nstructs-1)/2;
		if self.natoms != in_coords.shape[1]:
			raise ValueError, "Input weights do not match " \
						"structure set."
		
		# Center all structures.
		coords = self.center_all(in_coords, weights)
		print "Aligning %d structures using %d atoms."%(self.nstructs,\
								self.natoms)
		if(not g_align):
			self.simple_align(coords, weights)
			return
		elif(g_align < 0):
			self.align_one(coords, weights, -g_align-1)
			return
		
			# start calculation
		print "  Calculating initial msd error..."
		self.calcE0(coords, weights)
		E0 = sum(self.E0)/msd_structs
		print "    r<msd> = %f A"%(sqrt(E0))
		
		print "  Calculating P matrices..."
		self.calcP(coords, weights)
		print "  Doing initial rotation to first structure..."
		self.initp()
		E = E0 + sum(self.calcE())/msd_structs
		print "    Initial rotation, r<msd> = %f"%(sqrt(E))
		
		lE = 0.0
		cyc = 0		# convergence criteria
		
		while( abs(E - lE) > TOL and cyc < MAXCYC ):
			self.iterate()
			
			lE = E		# old error sum
			E = sum(self.E, 0)/msd_structs
			E = E + E0	# new error sum
			cyc += 1
			print "\tCycle %d, r<msd> = %f"%( cyc, sqrt(E) )
		if(abs(E - lE) < TOL):
			print "Fit converged in %d cycles!"%(cyc)
		else:
			print "Reached cycle quota, insert more quarters..."
	
	def center_all(self, coords, weights):
		self.cm = sum(weights[newaxis,:,newaxis]*coords, 1)
		#print "Centered from:"
		#print self.cm
		return coords-self.cm[:,newaxis,:]
	
	def getCM(self, s):
		return self.cm[s]
	
	def getR(self, s):
		p = self.p[s]
		
		R = array([ \
			[p[0]**2 - p[1]**2 - p[2]**2 + p[3]**2, \
			    2.0*(p[0]*p[1] - p[2]*p[3]), \
			    2.0*(p[0]*p[2] + p[1]*p[3])], \
			[2.0*(p[0]*p[1] + p[2]*p[3]), \
			    -p[0]**2 + p[1]**2 - p[2]**2 + p[3]**2, \
			    2.0*(p[1]*p[2] - p[0]*p[3])], \
			[2.0*(p[0]*p[2] - p[1]*p[3]), \
			    2.0*(p[1]*p[2] + p[0]*p[3]), \
			    -p[0]**2 - p[1]**2 + p[2]**2 + p[3]**2] \
			  ], float)
		
		return R
	
	def calcPI(self, A):
		PI = zeros((4,4), float)
		
		for B in range(self.nstructs):
			if(A == B):
				continue
			PI += self.getP(A,B)
		
		return PI
	
	def calcE(self):	# this is really (Err-E0)/2
		E = zeros(self.nstructs, float)
		
		for A in range(self.nstructs):
			E[A] = -1.0 * dot( \
				dot(self.p[A], self.calcPI(A)), \
				self.p[A])
	#		for B in range(self.nstructs()):
	#			E[A] += max(self.E0[A,B])
		
		self.E = E
		return E
	
	def calcE0(self, coords, weights):
		E0 = []
		
		for A in range(self.nstructs):
			E0.append([])
			for B in range(self.nstructs):
				if(A >= B):
					E0[A].append(0.0)
					continue
				E0[A].append( do_msd(weights, \
					coords[A], coords[B]) )
		self.E0 = E0
		
	def calcP(self, coords, weights):
		P = zeros((self.nstructs,self.nstructs, 4,4), float)
		
					# upper triangle
		for A in range(self.nstructs):
		  for B in range(self.nstructs):
		    if(A == B):
			pass
		    elif(A > B):		# Bermuda triangle
			P[A,B] = P[B,A]		# copy and negate
			P[A,B,3,:] *= -1.0	# last column and row
			P[A,B,:,3] *= -1.0
		    else:
			P[A,B] = self.makeP(coords[A],coords[B],weights)
		
		self.P = P
		
	def makeP(self,A,B,weights):
		P = zeros((4,4), float)
					# store M in P[0:3,0:3]
		for a in range(self.natoms):	# calculate upper diagonal
		  for i in range(3):
		    for j in range(3):
			P[i,j] += weights[a] * A[a,i] * B[a,j]
		
		V = [	P[1,2] - P[2,1], \
			P[2,0] - P[0,2], \
			P[0,1] - P[1,0] ]
		
		P = P + transpose(P) - 2*identity(4,float)*trace(P)
		P[3,3] = 0.0
		
		for i in range(3):
			P[3,i] = V[i]
			P[i,3] = V[i]
		
		return P
	
	def simple_align(self, coords, weights): # align to first structure
		p = [ array([0.0, 0.0, 0.0, 1.0], float) ]
		
		#print self.P[1,0]
		#print la.eigenvectors(self.P[1,0])
		
		print "    n, initial, final RMSD from zeroth structure:"
		for A in range(1,self.nstructs):
			PI = self.makeP(coords[A],coords[0],weights)
			p.append( getvec(PI) )
			E0 = do_msd(weights, coords[A], coords[0])
			E = -1.0 * dot(dot(p[A], PI), p[A])
			print "%5d %8.3f %8.3f"%(A, E0, E0+E)
		
		self.p = array(p)
		#print self.p
	
	
	def align_one(self, coords, weights, A):
		"""Align one structure to the rest of 'em."""
		PI = zeros((4,4), float)
		E = 0.0
		for B in range(self.nstructs):
			if B == A: # A must be >= 0
				continue
			PI += self.makeP(coords[A],coords[B],weights)
			E += do_msd(weights, coords[A], coords[B])
		p = getvec(PI)
		E -= 2*dot(dot(p, PI), p)
		E /= self.nstructs-1
		
		self.p = zeros((self.nstructs,4), float)
		self.p[:,0] = 1.0
		self.p[A] = p
		print "r<msd> for %d = %f\n"%(A+1,sqrt(E))
	
	def initp(self): # initially, align to structures before you
		p = [ array([0.0, 0.0, 0.0, 1.0], float) ]
		
		#print self.P[1,0]
		#print getvec(self.P[1,0])
		
		for A in range(1,self.nstructs):
			#PI = zeros((4,4), float)
			
			self.p = array(p)
			#for B in range(A):
			#	PI += self.getP(A,B)
			#p.append( getvec(PI) )
			p.append( getvec(self.getP(A, 0)) )
		#print 
		
		self.p = array(p)
		#print self.p
	
	def iterate(self):	# calculate PI, and set p to top eigenvec
		for A in range(self.nstructs):
			PI = self.calcPI(A)
			self.p[A] = getvec(PI)
			
			# update E[A]
			self.E[A] = -1.0 * dot( \
				dot(self.p[A], PI), \
				self.p[A])
		#print "\t    p[0] = " + str(self.p[0])
	
	def getP(self,A,B):
		p=array([[self.p[B,3], -self.p[B,2], self.p[B,1], self.p[B,0]],\
			 [self.p[B,2], self.p[B,3], -self.p[B,0], self.p[B,1]],\
			 [-self.p[B,1], self.p[B,0], self.p[B,3], self.p[B,2]],\
		       [-self.p[B,0], -self.p[B,1], -self.p[B,2], self.p[B,3]]\
				])
		P = dot(p, self.P[A,B])
		P = dot(P, transpose(p))
		
		return P

def do_msd(weights, A, B):
	atoms = len(weights)
	
	if(len(A) != atoms or len(B) != atoms):
		print "Error! cannot do RMSD calculation without some respect"
		return -1.0
	msd = 0.0
	
	for a in range(atoms):
		for i in range(3):
			msd += weights[a]*(A[a,i] - B[a,i])**2
	
	# weights should sum to one
	return msd # / sum(weights, 0)

def getvec(M):
	if(len(M.shape) != 2 or M.shape[0] != M.shape[1]):
		print "Invalid matrix for eigenvalue calculation!"
		return(zeros(M.shape[0], float))
	
	val, vec = la.eig(M)
	#print M
	#print val
	#print vec
	val = val.real
	
	top = 0
	for i in range(1, len(val)):
		if(val[i] > val[top]):
			top = i
	return vec[:,top].real

