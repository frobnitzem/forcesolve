#!/usr/bin/env python

# Wish your measurement data could be turned into a function? Look no further.
# This program takes column-formatted data (1st col. is x, "col"-th col. is y)
# and performs P-spline inference to create a B-spline coefficient file.

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
# 3/20/2008
# This work was sponsored by a DOE CSGF.

import sys
from cg_topol.ucgrad import *
from cg_topol.bspline import *
from cg_topol.pspline import *

UsageInfo = "Usage: %s [options] <data file> <out>\n"\
	"Options:\n"\
	"\t-h        Help me!\n"\
	"\t-p        Periodic spline (default = nonperiodic).\n"\
	"\t-r [4]    Spline order\n"\
	"\t-n [100]  Number of (internal) spline knots\n"\
	"\t-c [1]    Column of input to use as y-value\n"\
	"\t-d [2]    Order of derivative penalty\n"\
	"\t-m [0]    Power of x to be used in volume element\n"\
	"\t-s [2000] Number of MCMC samples to collect\n"\
	"\t-f [10]   Frequency of collecting parameter statistics\n"\
	"\t-t [10]    Number of statistics to toss at start\n"

def main(argv):
	periodic = False
	r = 4
	n = 100
	x0 = None
	L = None
	shift = None
	samples = 2000
	skip = 10
	toss = 10
	deriv = 2
	rho = 0
	col = 1
	while len(argv) > 1 and argv[1][0] == '-':
		if argv[1] == '-h':
			print UsageInfo % argv[0]
			return 0
		elif argv[1] == '-p':
			periodic = True
			del argv[1]
		elif argv[1] == '-r':
			r = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-n':
			n = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-c':
			col = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-d':
			deriv = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-m':
			rho = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-s':
			samples = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-f':
			skip = int(argv[2])
			del argv[1:3]
		elif argv[1] == '-t':
			toss = int(argv[2])
			del argv[1:3]
		elif argv[1] == '--':
			del argv[1]
			break
		else:
			print UsageInfo % argv[0]
			print "Unrecognized option: " + argv[1]
			return 1
	
	if len(argv) != 3:
		print UsageInfo % argv[0]
		return 1
	
	D = read_matrix(argv[1])
	assert len(D.shape) == 2, "Input is not a 2D matrix!"
	if x0 == None: # Default range is input range.
		x0 = D[0,0]
	if L == None:
		L = D[-1,0]-x0
	
	h = L/n
	if not periodic:
		n += r-1 # Add external points.
	if shift == None:
		if periodic:
			shift = 0.5*r
		else:
			shift = r-1
	f = spline_func((x0,h,n),Bspline(r),periodic,shift,(x0,x0+L), argv[1])
	
	m = pspline_match(f, deriv, rho, 1.0e-8)
	m.append(D[:,0], D[:,col])
	m.dimensionality()
	m.sample(samples,skip,toss)
	f.c = m.theta
	m.write_out(argv[2])
	
	# To do something with the collected samples, call this function.
	#calc_functional(m, argv[2]+".func")

# Code your own functional analysis here.
def calc_functional(m, out):
	# m.f - original spline function object
	# m.sum_theta -- list of theta, variance, and tension sample values.
	# m.sum_v
	# m.sum_alpha
	F = functional(m.sum_theta) # you have to write your own functional...
		# here I assume the output is (samples x what-have you)
	mean = sum( F, 0 )/m.samples
	var = sum( (F-mean)**2, 0 )/m.samples
	
	# move the data into a matrix, M
	write_matrix(out, M) 
		
if __name__=="__main__":
	main(sys.argv)
