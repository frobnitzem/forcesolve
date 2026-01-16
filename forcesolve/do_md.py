#!/usr/bin/python

# This code is a wrapper for md.py, a set of functions implementing Skeel and Izaguirre's Langevin dynamics method (Mol. Phys. 2002, 100(24):3885-91).  The pairwise interactions are implemented in an N^2 fashion, so don't expect too much.

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
# 5/13/2008
# This work was supported by a DOE CSGF.

import sys
from cg_topol import *
from cg_topol.ucgrad.array_io import *
from cg_topol.ucgrad.parser import *
from .md import *

UsageInfo =  "Usage: do_md.py [options] [files]\n" \
"Files:\n" \
"\t-p  [topology]    Molecule topology file\n"\
"\t-o  [output]      Prefix for all output files\n"\
"Options:\n" \
"\t-x  [coordinates] Initial coordinate file\n"\
"\t-T  [298.15]      Langevin thermostat temperature.\n"\
"\t-v  [velocities]  Initial velocity file\n"\
"\t-sx [1.0]         Scale coordinates by this amount on input.\n"\
"\t-sv [1.0]         Scale velocities by this amount on input.\n"\
"\t-dt [1.0]         Time-step between input frames (dx = v*dt).\n"\
"\t-ns 100           Number of steps to run.\n"\
"\t-nt 1             Number of time-steps per step.\n"\

RequiredFlags = ['p', 'o']
AcceptedFlags = ['x', 'T', 'v', 'sx', 'sv', 'dt', 'ns', 'nt']+RequiredFlags

def main(argv):
	if(len(argv) < 5):
		print(UsageInfo)
		return 1
	kT = 1.380603e-23/4.3597482e-18*300 # Hartree
	
	# Parse arguments.
	args, flags = parseflags(argv, UsageInfo, \
			AcceptedFlags, RequiredFlags)
	#print args, flags
	if flags < 0:
		return 1
	
	if flags.has_key('ns'):
		ns = int(flags['ns'][0])
	else:
		ns = 100
	if flags.has_key('nt'):
		nt = int(flags['nt'][0])
	else:
		nt = 1
	if flags.has_key('dt'):
		dt = float(flags['dt'][0])
	else:
		dt = 1.0
	if flags.has_key('T'):
		temp = float(flags['T'][0])
		print("Setting thermostat temperature to %f K"%(temp))
		kT = 1.380603e-23/4.3597482e-18*temp # Hartree
	
	# read parameter file
	topol = cg_topol(flags['p'][0])
	
	if flags.has_key('x'):
		x0 = read_array(flags['x'][0], (topol.atoms, 3))
	else:
		x0 = topol.pdb.x
	if flags.has_key('v'):
		v0 = read_array(flags['v'][0], (topol.atoms, 3))
	else:
		v0 = None
	
	md = md_integrate(x0, v0, topol, dt, kT)
	md.do_md(ns, nt, flags['o'][0])

def run():
	main(sys.argv)
