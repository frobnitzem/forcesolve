#!/usr/bin/env python

# Simple code to plot the function specified by an input spline file
# (along with its derivative in the second column).

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
# 12/11/2007
# This work was supported by a DOE CSGF.


import sys

from cg_topol.ucgrad.array_io import *
from numpy import *
from cg_topol.bspline import *

def main(argv):
        if len(argv) < 3:
                print "Usage: %s <in> <out>"%(argv[0])
                return 1
        
        f = read_spline_func(argv[1])
        x = arange(f.rng[0], f.rng[1], f.h*0.1)
        y = f.y(x,1)
        write_matrix(argv[2], transpose(array([x,y[0],y[1]])))

if __name__ == "__main__":
        main(sys.argv)
