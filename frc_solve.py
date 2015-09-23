#!/usr/bin/python

# Main program for frc_match.py force matching methods.
# Please convert your data to atomic units
# -- Hartree, Bohr, e- mass, atomic time units --
# before using to avoid any unpleasantness.
# Although you could use -sx and -kT I haven't tested this yet...

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
# 12/10/2007
# This work was supported by a DOE CSGF.

import sys
from cg_topol.ucgrad.array_io import *
from cg_topol.ucgrad.parser import *
from frc_match import *
from cg_topol import cg_topol, show_index

UsageInfo =  "Usage: frc_solve.py [options] [files]\n" \
"Files:\n" \
"\t-p  [topology]    Molecule topology file\n"\
"\t-x  [coordinates] Input coordinate sample/trajectory file\n"\
"\t-o  [output]      Prefix for all output files\n"\
"Options:\n" \
"\t-f  [forces]      Input forces sample file\n"\
"\t-sx [1.0]         Scale coordinates by this amount on input.\n"\
"\t-sv [1.0]         Scale velocities by this amount on input.\n"\
"\t-sf [1.0]         Scale forces by this amount on input.\n"\
"\t-dt [1.0]         Time-step between input frames (dx = v*dt).\n"\
"\t-kT [9.5e-4]      Boltzmann constant times temperature.\n"\
"\t-ns 1000          Number of samples to collect.\n"\
"\t-mle None         Prefix for writing initial maximum liklihood estimate.\n"

RequiredFlags = ['p', 'o', 'x']
AcceptedFlags = ['f', 'sx', 'sv', 'sf', 'dt', 'ns', 'mle', 'kT']\
				+RequiredFlags

def main(argv):
	if(len(argv) < 7):
		print UsageInfo
		return 1
	
	# Parse arguments.
	args, flags = parseflags(argv, UsageInfo, \
			AcceptedFlags, RequiredFlags)
	#print args, flags
	if flags < 0:
		return 1
	
	if flags.has_key('ns'):
		samples = int(flags['ns'][0])
	else:
		samples = 5000
	
	# read parameter file
	topol, pdb = cg_topol(flags['p'][0])
	
	dt = 1.0
	if flags.has_key('dt'):
		dt = float(flags['dt'][0])
		print "dt = %f"%dt
	else:
		print "Assuming dt = %f"%dt
	
	kT = 1.380603e-23/4.3597482e-18*300.0 # Hartree
	if flags.has_key('kT'):
		kT = float(flags['kT'][0])
		print "kT = %f"%kT
	else:
		print "Assuming kT = %f"%kT
	
	forces = frc_match(topol, pdb, dt, kT)
	
	# Classify and combine input into types
	print "Interaction type index:"
	show_index(forces.topol)
	if(forces.topol.params-len(forces.constraints) == 0):
		print
		print Usage
		print "Yow! No fit-able parameters present."
		print "Are you sure your input is complete and correctly "\
				"formatted?"
		return -1
	
	for xf, ff in zip(flags['x'],flags['f']):
		print "Appending %s %s"%(xf,ff)
		x = read_matrix(xf)
		x = reshape(x, (len(x)/pdb.atoms, pdb.atoms, 3))
		f = read_array(ff, x.shape)
		if flags.has_key('sx'):
			scale = float(flags['sx'][0])
			print "Scaling positions by %f"%(scale)
			x *= scale
		if flags.has_key('sf'):
			scale = float(flags['sf'][0])
			print "Scaling forces by %f"%(scale)
			f *= scale
		forces.append(x, f) # Force matching.
	
	forces.dimensionality()

        # Just do the right thing (TM).
        if not flags.has_key('mle') \
            and (forces.topol.hyp_params == 0 or samples == 0):
                flags['mle'] = flags['o']

        if flags.has_key('mle'):
		print "Finding maximum posterior estimate..."
		forces.maximize()
		print "Writing mle to prefix \"%s\"..."%(\
					flags['mle'][0])
		forces.write_out(flags['mle'][0])

        if forces.topol.hyp_params == 0 or samples == 0:
            return 0
	
	print "Sampling spline coefficients..."
	forces.sample(samples, 500, 5)
	
	print "Writing output to prefix \"%s\"..."%(flags['o'][0])
	forces.write_out(flags['o'][0])
	return 0

if __name__ == "__main__":
        main(sys.argv)

