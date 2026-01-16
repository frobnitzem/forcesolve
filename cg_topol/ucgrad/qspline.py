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

# Quartic spline classes to do 1D and 2D quartic spline interpolation.
# The quartic splines developed allow fitting a function and its derivative
# at every point and require solution of only one unknown.

# Author: David Rogers
# Sponsored by a Department of Energy CSGF
# University of Cincinnati, February, 2007

from numpy import *
from .array_io import *

# FIXME: Currently assumes input list is sorted in ascending x-order
#        and that all x-s are distinct.

# Dat contains at least 2 columns for x and y values
# the third is optional and supplies the first and last y' or y'' value.
class qspline:
	def __init__(self, dat, a="min"):
		if(len(dat.shape) != 2 or dat.shape[1] < 3):
		    raise runtimeError, "Data must contain at least 3 columns!"
		
		if(dat.shape[1] > 2): # A and B ignored.
			self.spl = spline_table(dat,mode, dat[0,2],dat[-1,2])
		else:
			self.spl = spline_table(dat, mode, a, b)
		self.n = self.spl.shape[0]
	
	def find_x(self, x):
		l = 0
		u = self.n-1
		
		# Find relevant section of table
		while(u-l > 1):
			t = (u+l)/2
			#print "l %d  u %d  t %d"%(l,u,t)
			if(self.spl[t,0] > x):
				u = t
			else:
				l = t
		#print "Using section: " + str(self.spl[l:u+1])
		return l,u
		
	# Note that this routine will still return values for
	# points outside its data range using the cubic
	# equation coefficients from the end-ranges!
	def y(self, x):
		l, u = self.find_x(x)
		
		h = self.spl[u,0]-self.spl[l,0]
		h26 = h*h/6.0
		a = (self.spl[u,0]-x)/h
		b = 1.0-a
		ya = self.spl[l,1]+(a*a-1.0)*self.spl[l,2]*h26
		ya *= a
		yb = self.spl[u,1]+(b*b-1.0)*self.spl[u,2]*h26
		yb *= b
		return ya+yb
	
	def yp(self, x): # quadratic interpolation
		l, u = self.find_x(x)
		
		h = self.spl[u,0]-self.spl[l,0]
		a = (self.spl[u,0]-x)/h
		b = 1.0-a
		
		a = h*(3*a*a-1.0)/6.0
		b = h*(3*b*b-1.0)/6.0
		
		yp = (self.spl[u,1]-self.spl[l,1])/h
		yp += b*self.spl[u,2]-a*self.spl[l,2]
		
		return yp
	
	def ypp(self, x): # linear interpolation
		l, u = self.find_x(x)
		
		h = self.spl[u,0]-self.spl[l,0]
		a = (self.spl[u,0]-x)/h
		
		return a*self.spl[l,2]+(1.0-a)*self.spl[u,2]

# Cspline class for equally-spaced meshes.
# Data contains 1 or 2 columns -- for y and optionally y' or y''
class cspline_h:
	def __init__(self, data, x0, h, mode=2, a=0., b=0.):
		dat = zeros((data.shape[0],2), float)
		dat[:,0] = arange(x0, h*(data.shape[0]-0.5)+x0, h)
		if(len(data.shape) > 1): # test for broken >2 cond.?
			dat[:,1] = data[:,0]
		else:
			dat[:,1] = data
		
		if(len(data.shape) > 1 and data.shape[1] > 1): # A and B ignored
			self.spl = spline_table(dat,mode, data[0,1],data[-1,1])
		else:
			self.spl = spline_table(dat, mode, a, b)
		self.n = self.spl.shape[0]
		self.x0 = x0
		self.h = h
		self.h26 = h*h/6.0
	
	def find_x(self, x):
		# Find relevant section of table
		x -= self.x0
		x /= self.h
		l = int(x)
		return l, l+1
	
	# Note that this routine will not return 0.0 for
	# points outside its data range.
	def y(self, x):
		l, u = self.find_x(x)
		
		if(l<0 or u>=self.n):
			return 0.0
		
		a = u-x
		b = 1.0-a
		ya = self.spl[l,1]+(a*a-1.0)*self.spl[l,2]*self.h26
		ya *= a
		yb = self.spl[u,1]+(b*b-1.0)*self.spl[u,2]*self.h26
		yb *= b
		return ya+yb
	
	def yp(self, x): # quadratic interpolation
		l, u = self.find_x(x)
		
		if(l<0 or u>=self.n):
			return 0.0
		
		h = self.spl[u,0]-self.spl[l,0]
		a = (self.spl[u,0]-x)/h
		b = 1.0-a
		
		a = h*(3*a*a-1.0)/6.0
		b = h*(3*b*b-1.0)/6.0
		
		yp = (self.spl[u,1]-self.spl[l,1])/h
		yp += b*self.spl[u,2]-a*self.spl[l,2]
		
		return yp
	
	def ypp(self, x): # linear interpolation
		l, u = self.find_x(x)
		
		if(l<0 or u>=self.n):
			return 0.0
		
		h = self.spl[u,0]-self.spl[l,0]
		
		h = self.spl[u,0]-self.spl[l,0]
		a = (self.spl[u,0]-x)/h
		
		return a*self.spl[l,2]+(1.0-a)*self.spl[u,2]


# Read an input x,y list and return an appropriate cspline object.
# Intelligent enough to figure out if you want cspline or cspline_h
def cspline_file(filename, x_col=0, y_col=1, mode=2, e_col=-1, a=0., b=0.):
	data = read_matrix(filename)
	n = data.shape[0]
	if(data.shape[1] < 2):
		raise runtimeError, "Error! Data file must have >1 col."
	
	dat = zeros((data.shape[0],3), float)
	dat[:,0] = data[:,x_col]
	dat[:,1] = data[:,y_col]
	
	h = dat[1,0]-dat[0,0]
	test_h = dat[:,0] - arange(dat[0,0], (n-0.5)*h+dat[0,0], h)
	# Chi-Squared-test for equally-spacedness
	#print sum(test_h*test_h)/sum(dat[:,0]*dat[:,0])
	if(sum(test_h*test_h)/sum(dat[:,0]*dat[:,0]) < 1.0e-10):
	    print "Creating cspline_h object for \"%s\"."%(filename)
	    if(e_col >= 0):
		dat[:,2] = data[:,e_col]
		return cspline_h(dat[:,1:3], dat[0,0], h, mode)
	    else:
		return cspline_h(dat[:,1], dat[0,0], h, mode, a, b)
	else:
	    print "Creating cspline object for \"%s\"."%(filename)
	    if(e_col >= 0):
		dat[:,2] = data[:,e_col]
		return cspline(dat, mode)
	    else:
		return cspline(dat[:,0:2], mode, a, b)
	
# Create a spline table for the given x,y data with a and b
# as the "mode"nd derivatives at the endpoints.
def spline_table(dat, mode, a, b):
	if(mode != 1 and mode != 2):
		raise runtimeError, "Unacceptable value entered for mode."
	
	# Table of position, value and second derivative
	table = zeros((dat.shape[0],3), float)
	table[:,0:2] = dat[:,0:2] # Copy x and y values.
	
	x = table[:,0]
	y = table[:,1]
	y2 = table[:,2]
	u = zeros(dat.shape[0], float)
	
	if(mode == 2):
		y2[0]  = u[0]  = a
		y2[-1] = u[-1] = b
	else:
		y2[0] = -0.5
		u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-a)
		
		y2[-1] = 0.5
		u[-1] = (3.0/(x[-1]-x[-2]))*(b-(y[-1]-y[-2])/(x[-1]-x[-2]))
	
	for i in range(1, dat.shape[0]-1):
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
		p = sig*y2[i-1]+2.0
		y2[i] = (sig-1.0)/p
		u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p
	
	y2[-1] = (u[-1]-y2[-1]*u[-2])/(y2[-1]*y2[-2]+1.0)
	for i in range(dat.shape[0]-2, -1, -1):
		y2[i] = y2[i]*y2[i+1]+u[i]
	
	return table

def plot_spline(spl, h):
	x = arange(spl.spl[0,0], spl.spl[-1,0]+0.5*h, h)
	y = zeros((len(x),3), float)
	for i,v in enumerate(x):
		y[i,0] = spl.y(v)
		y[i,1] = spl.yp(v)
		y[i,2] = spl.ypp(v)
	return x, y
