#!/usr/bin/python

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

# David Rogers
# Dr. Thomas Beck Lab
# October 9th, 2006
# Sponsored in part by the Department of Energy CSGF
# and another part by a DoD/Army grant.

import sys
from cg_topol.ucgrad.prot_files import *
from cg_topol.ucgrad.parser import *
from cg_topol.ucgrad.array_io import write_matrix

AcceptedFlags = ['d', 'i', 'o', 's', 'm']
RequiredFlags = ['d', 'i', 'o']
UsageInfo="Usage: coarsen.py -d <site definition> -i <pdb in> -o <pdb out>\n" \
	"\tOptions:\n" \
	"\t\t-s [1.0]                Scale factor to apply to output\n"\
	"\t\t-m <out.mat>            Matrix for coarsening raw coord.s\n"

nowrap = 1 # Set to 0 if you want to wrap the edges of residues.

def main(argv):
	args, flags = parseflags(argv, UsageInfo, AcceptedFlags, RequiredFlags)
	if(flags == -1):
		return -1
	
	if('s' in flags.keys()):
		scale = float(flags['s'][0])
	else:
		scale = 1.0
	
	j = 0
	for d,i,o in zip(flags['d'], flags['i'], flags['o']):
		cg = coarsen(d, o, scale)
		structloop(i, cg.addstruct)
		if 'm' in flags.keys():
			if j >= len(flags['m']):
				print "Warning! Unable to write matrix for "\
				"%s %s -> %s (not enough output names)."%(d,i,o)
			else:
				write_matrix(flags['m'][j], cg.M)
		j += 1
	
	return 0

class coarsen:
	def __init__(self, def_name, out_name, scale=1.0):
		self.out_name = out_name
		out = open(out_name, 'w')
		out.truncate()
		out.close()
		
		print "Creating coarsening from definitions in %s"%(def_name)
		self.clist = parse_cgdef(def_name, scale)
		#print self.clist
	
	def addstruct(self, pdb, warn=0):
		cpdb = protein()
		M = []
		n = 0
		i = 0
		for res in pdb.res: # number residue atoms consecutively.
			res.atom_zero = i
			i += len(res.names)
		for rn,res in enumerate(pdb.res):
		    if res.name in self.clist.keys():
			#print "Attempting to add residue %s %d using spec:"%(\
			#		res.name, rn+1)
			#print self.clist[res.name]
			for cres,rspec in self.clist[res.name].iteritems():
			    n += 1
			    ra, rm = calcres(pdb, rn, rspec, cres, n)
			    cpdb.append(ra)
			    M += rm
		    elif(warn):
			print "Warning! Fine residue %s %d does not map to "\
				"any coarse coordinates."%(res.name, res.num)
				
		cpdb.write(self.out_name, 'a')
		del cpdb
		self.M = array(M)

# Given a dictonary connecting coarse atom names with coarsening spec.s,
# calculate all atom positions and return a residue.
def calcres(pdb, rn, spec, name, serial=1):
	#print "Adding coarse residue %s %d using spec:"%(name, serial)
	#print spec
	
	M = []
	res = residue([], name, serial, pdb.res[rn].chain_id)
	for aname, aspec in spec.iteritems():
		if(nowrap):
		  if(sometrue(rn+array(aspec[1]) < zeros(len(aspec[1]))) or \
		     sometrue(rn+array(aspec[1]) >= len(pdb.res) * \
						ones(len(aspec[1])))):
			continue # ignore end-effects
		#print "Calculating coarse atom %s from spec:"%(aname)
		#print aspec
		x, m = calc_cg_pos(pdb, rn, aspec)
		if sum(m) > 1.0e-10:
			res.append(aname, x)
			M.append(m)
	return res, M

def calc_cg_pos(pdb, rn, aspec):
	x = zeros(3)
	m = zeros(len(pdb.x))
	sc = 0.0
	for aname,ro,c in zip(*aspec):
		sc += c
		ri = (rn+ro+len(pdb.res))%len(pdb.res)
		ai = pdb.res[ri].get_atom_index(aname)-pdb.res[ri].start_atom
		if ai < 0: # One of the atoms is missing.
			print "Abandoning residue %s %d: missing %d \"%s\""%(\
					pdb.res[rn].name, rn+1, ro, aname)
			return x, -1.0
		ai += pdb.res[ri].atom_zero
		m[ai] += c
		#print "Atom %s at "%(aname) + str(pdb.x[ai])
		x += c*pdb.x[ai]
	#print "Returning: " + str(x)
	return x, m

# Parses site definition file into dictionary data object:
# index = fine residue name
# value = {coarse_residue_name:coarse_residue_spec, 
#	coarse_residue_name:coarse_residue_spec, ... }
# coarse_residue_spec = dictionary:
# index = atom_name
# value = [ [fine atom names], [fine atom res_offsets], [fine_atom_wts] ]
def parse_cgdef(filename, scale=1.0):
	file = open(filename)
	natom = 0
	cpatom = 0
	cdef = {}
	n = 0
	for line in file.xreadlines():
		n += 1
		chk = line.strip()
		if not chk or chk[0] == '#':
			continue
		if natom == 0 and cpatom == 0: # New residue def.
		    if(line[0:4].upper() == "SITE"):
			ca,cr,natom,fr,w = parse_site(line, scale, n)
			catom = [ca] # Possible list of def.ns
			cres = [cr]
			fres = [fr]
			wt = [w]
		    elif(line[0:4].upper() == "COPY"):
			cr,cpatom,fr,w = parse_copy(line, scale, n)
			cres = [cr]
			fres = [fr]
			wt = [w]
		    else:
			raise RuntimeError, "Syntax error on line %d: "\
				"expected SITE or COPY keyword\n"%(n)
		    catom_spec = [[], [], []]
		# Deal with multiple definitions using same atom list
		# and different residue names and/or coarse atom names.
		elif(cpatom == 0 and line[0:4].upper() == "SITE"):
			if(len(catom_spec[0]) > 0):
			    raise RuntimeError,"Error on line %d! New SITE def"\
				"inition before end of last atom list!"%(n)
			ca,cr,nat,fr,w = parse_site(line, scale, n)
			if nat != natom:
			  raise RuntimeError, "Error on line %d! New SITE defi"\
				"tion is atom # incompatible with previous."%(n)
			catom.append(ca) # Possible list of def.ns
			cres.append(cr)
			fres.append(fr)
			wt.append(w)
		elif(natom == 0 and line[0:4].upper() == "COPY"):
			cr,cat,fr,w = parse_copy(line, scale, n)
			if cat != cpatom:
			  raise RuntimeError, "Error on line %d! New COPY defi"\
			    "nition is atom # incompatible with previous."%(n)
			cres.append(cr)
			fres.append(fr)
			wt.append(w)
		elif cpatom:
			if(len(line) < 8):
			  raise RuntimeError, "Error on line %d! Copy line "\
				"must contain 2 4-char. atom names."%(n)
			catom = line[0:4]
			catom_spec[0] = [ line[4:8] ]
			catom_spec[1] = [ 0 ] # offset assumed to always be 0
			catom_spec[2] = [ 1.0 ]
			cpatom -= 1
			j = 0
			for fr,cr,w in zip(fres,cres,wt):
			    catom_spec[2][0] = w
			    # Add to spec. list.
			    if(fr not in cdef.keys()):
				cdef[fr] = {cr:{}}
			    cdef[fr][cr][catom] = catom_spec[:]
			    j += 1
		else:
			catom_spec[0].append(line[0:4])
			tok = line[4:].split()
			catom_spec[1].append(int(tok[0]))
			if(len(tok) >= 2 and tok[1][0] != '#'):
				catom_spec[2].append(float(tok[1]))
			else:
				catom_spec[2].append(1.0)
			natom -= 1
			if(natom == 0): # Time to add atom into spec.
				norm = 0.0 # normalize weights
				for a in catom_spec[2]:
					norm += abs(a)
				#print "Line %d: Appending def: "%(n)
				#print catom_spec
				#print "After dividing by %f"%(norm)
				# Normalize
				for i in range(len(catom_spec[2])):
					catom_spec[2][i] *= 1.0/norm
				wt.append(1.0) # initial wt = 1.0
				j = 0 # copy counter
				for fr,cr,ca in zip(fres,cres,catom):
				    # Scale by individual weight.
				    for i in range(len(catom_spec[2])):
					catom_spec[2][i] *= wt[j]/wt[j-1]
				    # Add to spec. list.
				    if(fr not in cdef.keys()):
					cdef[fr] = {cr:{}}
				    cdef[fr][cr][ca] = catom_spec[:]
				    j += 1
				#print catom_spec
	file.close()
	if(cpatom != 0):
		raise RuntimeError, "Error! Line %d: Residue %s from %s "\
				"COPY is incomplete!"%(n, cres, fres)
	if(natom != 0):
		raise RuntimeError, "Error! Line %d: Residue %s %s from %s "\
				"SPEC is incomplete!"%(n, cres, catom, fres)
	
	return cdef

def parse_site(line, scale, n):
	if len(line) < 9:
		raise RuntimeError, "Error parsing SITE statement, line %d!"%(n)
	tok = line[9:].split()
	if len(tok) < 4:
		raise RuntimeError, "Error parsing SITE statement, line %d!"%(n)
	if(tok[2].upper() != "FROM"):
		raise RuntimeError, "Syntax error on line %d: "\
		    "expected FROM keyword in SITE def.\n"%(n)
	
	catom = line[5:9]
	cres = tok[0]
	natom = int(tok[1])
	fres = tok[3]
	if(len(tok) > 5 and tok[4].upper() == "WT"):
		wt = float(tok[5])*scale
	else:
		wt = scale
	return catom, cres, natom, fres, wt

def parse_copy(line, scale, n):
	tok = line[5:].split()
	if len(tok) < 4:
		raise RuntimeError, "Error parsing COPY statement, line %d!"%(n)
	if(tok[2].upper() != "FROM"):
		raise RuntimeError, "Syntax error on line %d: "\
		    "expected FROM keyword in COPY def.\n"%(n)
	cres = tok[0]
	natom = int(tok[1])
	fres = tok[3]
	if(len(tok) > 5 and tok[4].upper() == "WT"):
		wt = float(tok[5])*scale
	else:
		wt = scale
	return cres, natom, fres, wt

if __name__=="__main__":
	main(sys.argv)
