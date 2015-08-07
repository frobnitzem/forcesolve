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

# misc. functions

from numpy import *

def angle(x, y):
	x /= sqrt(sum(x**2))
	y /= sqrt(sum(y**2))
	th = x[0]*y[0]+x[1]*y[1]+x[2]*y[2]
	th = arccos(th)
	return th*180.0/pi

def getpt(coord):
	pt = []
	tok = coord.split()
	for i in range(3):
		pt.append(float(tok[i]))
	return array(pt)

# computes torsion assuming vectors are connected just so:
#	 -y--->
#	  ^ (right thumb)
#	  |
#	<-x-
#   i.e. the angle required to rotate x through to y
def get_tors(x, y, n):
	n /= sqrt(sum(n**2))
	x = x - dot(x, n)*n
	y = y - dot(y, n)*n
	
	proj = y[0] * (n[1]*x[2]-n[2]*x[1]) \
		+ y[1] * (n[2]*x[0]-n[0]*x[2]) \
		+ y[2] * (n[0]*x[1]-n[1]*x[0])
	
	return arctan2(proj, dot(x, y))

def get_th_phi(x, y, n):
	r_n = sqrt(dot(n, n))
	r_x = sqrt(dot(x, x))
	r_y = sqrt(dot(y, y))
	n /= r_n
	
	th1 = arccos(dot(x, n)/r_x)*180.0/pi
	th2 = arccos(dot(y, n)/r_y)*180.0/pi
	x = x - dot(x, n)*n
	y = y - dot(y, n)*n
	
	proj = y[0] * (n[1]*x[2]-n[2]*x[1]) \
		+ y[1] * (n[2]*x[0]-n[0]*x[2]) \
		+ y[2] * (n[0]*x[1]-n[1]*x[0])
	
	phi = arctan2(proj, dot(x, y))
	return th1, th2, phi

def get_ic(x, connect=[]):
	if(len(x.shape) != 2 or x.shape[1] != 3):
	    raise ValueError, "Coordinates have wrong shape: " + str(x.shape)
	
	if(connect == []): # linear chain (default behavior)
		connect = range(-1, x.shape[0]-1)
	
	if(len(connect) != x.shape[0]):
		raise ValueError, "Connection list has wrong length: %d out of"\
				" %d atoms"%(len(connect), x.shape[0])
	
	# Coordinates + 3 universal orientation points
	y = resize(x, (x.shape[0]+3, 3))
	y[-1] = array([0., 0.,  0.]) # closest to an actual atom
	y[-2] = array([0., 0., -1.])
	y[-3] = array([1., 0., -1.]) # first to get dropped off list
					# as we travel along the chain
	#parent = [p-1 for p in connect] # Convert to 0-offset list.
	parent = connect[:]
	parent.append(-3) # -2 has parent -3 and nobody asks about -3's mama.
	parent.append(-2) # -1 has parent -2
	
	n = x.shape[0]
	ic = zeros((n, 3), float)
	
	for a,c3 in enumerate(x):
		p2 = parent[a]
		p1 = parent[p2]
		p0 = parent[p1]
		c0 = y[p0] # Path to atom is 0-1-2-3 with 3
		c1 = y[p1] # as atom "a"
		c2 = y[p2]
		
		# bond length
		d = c3 - c2
		r = sqrt(sum(d*d))
		if r > 1.0e-20:
			d /= r
		
		# angle
		n = c2-c1
		nr = sum(n*n)
		
		if(nr > 1.0e-20):
			n /= sqrt(nr)
		th = arccos(-dot(n, d))
		
		# torsion
		phi = get_tors(c0-c1, d, n)
		
		ic[a,0] = r
		ic[a,1] = th
		ic[a,2] = phi
	
	return ic

# Rotates x with the unique rotation
# that transforms b0 to lie along the x-axis,
# and b1 is on the x,z plane with positive z-value.
def orientate(x, b0, b1):
	R = build_rot_to(b0, array([1.,0.,0.]))
	y = dot(R, b1) # New b1
	y[0] = 0. # Subtract projection onto x-axis
	ym = sum(y*y)
	if(ym < 1.0e-20):
		print "Warning! Second vector is colinear with first!"
		print "Unable to satisfy second orientation criterion uniquely."
	else:
		R = dot(build_rot_to(y,array([0.,0.,1.])), R)
	return dot(x, transpose(R))

# returns the cosine of the torsion angle, assuming n is normalized
def get_ctor(x, y, n):
	x = x - n*dot(x, n)
	y = y - n*dot(y, n)
	
	x = x/sqrt(dot(x, x))
	y = y/sqrt(dot(y, y))
	
	return dot(x, y)

# calculate the coordinate frame transformation matrix
def ic_matrix(r, th, phi):
	cth = cos(th)
	sth = sin(th)
	cphi = cos(phi)
	sphi = sin(phi)
	
	T = array([ [ -cth*cphi, -sphi,  sth*cphi,  r*sth*cphi ], \
		    [ -cth*sphi,  cphi,  sth*sphi,  r*sth*sphi ], \
		    [ -sth     ,   0.0, -cth     , -r*cth      ], \
		    [       0.0,   0.0,       0.0,         1.0 ] ], float)
	
	
	return T

# return x cross y in 3D
def cross_prod(x, y):
	return array( [ x[1]*y[2] - x[2]*y[1], \
			x[2]*y[0] - x[0]*y[2], \
			x[0]*y[1] - x[1]*y[0] ] )

# get the angle between the two vectors
def get_angle(x, y):
	x = dot(x, y)/sqrt(dot(x,x)*dot(y,y))
	return arccos(x)

# build a rotation matrix that will rotate vector x onto the objective
def build_rot_to(cur, obj):
	x = array(cur)
	y = array(obj)
	n = cross_prod(x, y)
	m = sum(n*n)
	if(m > 1.0e-20):
		th = get_angle(x, y)
		n /= sqrt(m)
		return build_rotation(n, th)
	else: # colinear... (obj = a*cur) -- Note: <obj,obj> = <cur,obj>*a
		if(dot(x,y) < 0.0): # a < 0.0 -- need to flip
			return array([[-1.,0.,0.],[0.,0.,-1.],[0.,-1.,0.]])
			# invert through center and swap y and z coord.s
			# effectively rotating about the x-axis
			# -- necessary to preserve chirality
		else:
			return identity(3, float)

# build rotation matrix to rotate about an arbitrary vector
def build_rotation(n, theta):
	m = sum(n*n)
	if(m < 1.0-1.0e-10 or m > 1.0+1.0e-10):
		n /= sqrt(m)
	
	trans = zeros((3,3), float)
	
	s = sin(theta)
	c = cos(theta)
	trans[0,0] = c + n[0]*n[0]*(1.0-c)	# Rodrigues' Rotation Formula
	trans[0,1] = n[0]*n[1]*(1.0-c) - n[2]*s
	trans[0,2] = n[1]*s + n[0]*n[2]*(1.0-c)
	
	trans[1,0] = n[2]*s + n[0]*n[1]*(1.0-c)
	trans[1,1] = c + n[1]*n[1]*(1.0-c)
	trans[1,2] = -n[0]*s + n[1]*n[2]*(1.0-c)
	
	trans[2,0] = -n[1]*s + n[0]*n[2]*(1.0-c)
	trans[2,1] = n[0]*s + n[1]*n[2]*(1.0-c)
	trans[2,2] = c + n[2]*n[2]*(1.0-c)
	
	return trans
