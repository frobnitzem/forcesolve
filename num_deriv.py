# Numerical differentiation function for use in
# testing derivative implementations for new forcefield terms.

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
# 12/5/2007
# This work was supported by a DOE CSGF.

from numpy import *

def num_deriv(f, x):
	ih = 2.0**10
	h = 0.5/ih
	
	oshape = x.shape
	x = reshape(x, x.size)
	dx = zeros(x.size, float)
	x0 = x.copy()
	x1 = x.copy()
	for i in range(x.size):
		x0[i] = x[i]-h
		x1[i] = x[i]+h
		dx[i] = (f(reshape(x1,oshape))-f(reshape(x0,oshape)))*ih
		x0[i] = x[i]
		x1[i] = x[i]
	x = reshape(x, oshape)
	return reshape(dx, oshape)
