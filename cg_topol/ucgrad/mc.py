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
import numpy.random as rand

# calculates the acceptance probability for a given value of dev
def calc_acc(sys, dev, beta, n=500):
	acc = 0
	i = 0
	
	en = sys.energy()
	x = sys.get_x()
	for i in xrange(n):
		sys.kick(dev)
		n_en = sys.energy()
		if(n_en < en or exp(beta*(en-n_en)) > rand.random()):
			en = n_en
			x = sys.get_x()
			acc += 1
		else:
			sys.set_x(x)
	
	return float(acc)/float(n)

def find_dev(sys, beta, dev=0.01, ratio=0.5, delta = 0.1, TOL = 0.01):
	maxit = 1000
	minit = 10
	i = 1
	
	acc = calc_acc(sys, dev, beta)
	err = abs(acc-ratio)
	while((err > TOL and i < maxit) or i < minit):
		print "%d: dev = %f, acc = %f, err = %f, delta = %f"%(i, \
				dev, acc, acc-ratio, delta)
		dev *= exp(delta*(acc-ratio))
		
		l_err = err
		acc = calc_acc(sys, dev, beta)
		err = abs(acc-ratio)
		if(err > l_err):
			delta *= 0.8
		else:
			delta *= 1.2
		i += 1
	print "%d: dev = %f, acc = %f, err = %f, delta = %f"%(i, \
			dev, acc, acc-ratio, delta)
	
	return dev

def do_mc(sys, beta, sys_func, n, skip=1000, dev=0.01, toss=100):
	acc = 0
	i = 0
	maxit = skip*(n+toss)
	
	en = sys.energy()
	x = sys.get_x()
	for i in range(n+toss):
	    for j in range(skip):
		sys.kick(dev)
		n_en = sys.energy()
		if(n_en < en or exp(beta*(en-n_en)) > rand.random()):
			en = n_en
			x = sys.get_x()
			acc += 1
		else:
			sys.set_x(x)
	    if(i >= toss):
		sys_func(sys)
	
	return float(acc)/float(maxit)

