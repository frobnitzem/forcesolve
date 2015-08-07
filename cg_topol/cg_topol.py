# Coarse-grained topology creation, MD interface, and I/O routines.

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
# 12/7/2007
# This work was supported by a DOE CSGF.

from numpy import *
from ucgrad.parser import *
from ucgrad.prot_files import *
from bspline import *

RequiredParams = ["EXPDB", "CONST_DIR"]

class cg_topol:
	def __init__(self, pfile):
		params = parseparam(pfile, RequiredParams)
		if params == -1:
			raise InputError, "Invald topology file."
		
		# Parse type indices into self.tindex.
		self.tindex = {}
		self.tmass = {}
		self.const_dir = params['CONST_DIR'][0][0]
		del params['CONST_DIR']
		
		# Parse example PDB.
		self.pdb = read_pdb(params['EXPDB'][0][0], 1)[0]
		del params['EXPDB']
		self.atoms = self.pdb.atoms
		# Computes the atom index at the start of each residue.
		res_start = [ len(res.names) for res in self.pdb.res ]
		res_start = reduce(lambda x, y: x+[x[-1]+y], \
						res_start, [0])
		for r, s in zip(self.pdb.res, res_start):
			r.atom_zero = s
		
		self.init_edge()
		# Add the following important keys in sequence.
		for k in ['EDGE', 'REDGE']:
		    app_method = "append_"+k.lower()
		    if not hasattr(self, app_method):
			raise RuntimeError, "No %s method found!"%app_method
		    app_method = getattr(self, app_method)
		    
		    if params.has_key(k):
			for line in params[k]:
				app_method(line)
			del params[k]
		print "%d total edges initialized."%(len(self.edge))
		self.build_conn()
		
		# Keywords add to forcefield term sets.
		self.ff_terms = []
		for k in params:
		    try: app_method = getattr(self, "append_"+k.lower())
		    except AttributeError:
			print "Warning! Unknown keyword in input file: \""+k+\
					"\" ignoring..."
			continue
		    
		    # Initialize term.
		    try:
			init_method = getattr(self, "init_%s"%k.lower())
		    except AttributeError:
			#print "Warning, no init method for ff term %s"%(k)
			init_method = None
		    if init_method != None:
			init_method()
		
		    # Add all such terms.
		    for line in params[k]:
			app_method(line)
		    try: # Print a nice summary.
			print "%d total %ss initialized."%(\
				len(getattr(self, k.lower())), k.lower())
		    except AttributeError:
			pass
		    except TypeError:
			pass
		    
		    if hasattr(self, "finalize_%s"%k.lower()):
		    	self.ff_terms.append(k.lower())
		
		self.build_name_index()
		# Finalize all terms.
		for k in self.ff_terms:
		    fin_method = getattr(self, "finalize_%s"%k)
		    fin_method()
	
	# Method to determine if the present function has all the bits and
	# pieces necessary to do MD.
	def can_do_md(self):
		for i in self.ff_terms:
		    if not hasattr(self, "%s_force"%i):
			return False
		self.commit_ff()
		return True
	# Method to determine if the present function has all the bits and
	# pieces necessary to do force matching.
	def can_force_match(self):
		for i in self.ff_terms:
		    if not hasattr(self, "design_%s"%i):
			return False
		    if not hasattr(self, "constrain_%s"%i):
			return False
		self.count_parameters()
		return self.parameters > 0
	
	def count_parameters(self):
		nt = []
		for k in self.ff_terms:
			design = getattr(self, "design_%s"%k)
			dmat = design(self.pdb.x)
			nt += [d.shape[-1] for d in dmat]
		self.types = len(nt)
		self.nt = array(nt, int)
		self.parameters = sum(self.nt)
	
	def write_all_splines(self, base):
		print "Writing parameters to %s..."%base
		for i in self.ff_terms:
		    try:
			write_method = getattr(self, "write_%s"%i)
		    except AttributeError:
			if hasattr(self, "read_%s"%i):
			  print "Warning, no write method for ff term %s"%(i)
			continue
		    write_method(base)
	def read_all_splines(self):
		print "Reading parameters from directory %s"%self.const_dir
		base = os.path.join(self.const_dir, "")
		for i in self.ff_terms:
		    try:
			read_method = getattr(self, "read_%s"%i)
		    except AttributeError:
			if hasattr(self, "write_%s"%i):
			  print "Warning, no read method for ff term %s"%(i)
			continue
		    read_method(base)
	
	# Computes the (residue, atom) name tuple from an atom index.
	def lookup_rname(self, i):
		for r in pdb.res:
		    if r.atom_zero <= i:
			i -= r.atom_zero
			return (r.name, r.names[i].split()[0])
		else:
			raise ValueError, "Invalid atom passed to lookup_rname."
	
	# Builds the internal connection table from a set of edges (self.edge).
	def build_conn(self):
		self.conn = [] # Connection table.
		for i in range(self.atoms):
			self.conn.append(set())
		for i,j in self.edge:
			self.conn[i].add(j)
			self.conn[j].add(i)
	
	def build_name_index(self):
		self.name_con = {} # (rname, aname) -> (res,atom) number index.
		self.names = [] # i->(rname,aname,type)
		r = 0
		for i in range(self.atoms):
			k = i-self.pdb.res[r].atom_zero
			while k >= len(self.pdb.res[r].names):
				r += 1
				k = i-self.pdb.res[r].atom_zero
			rname = (self.pdb.res[r].name, \
					self.pdb.res[r].names[k].split()[0])
			if not self.name_con.has_key(rname):
				self.name_con[rname] = []
			self.name_con[rname].append((r,i))
			try: # Checks atom typing as well.
				t = self.tindex[rname]
			except KeyError:
				t = "%s:%s"%(rname[0],rname[1])
				print "Warning! Auto-typing %s."%t
				self.tindex[rname] = t
			self.names.append((rname[0],rname[1],t))
	
	def init_edge(self):
		self.edge = set()
	
	def add_edge(self, i, j):
		if i < j:
			self.edge.add((i,j))
		else:
			self.edge.add((j,i))
	
	def add_all_redge(self, res_table, oi, ai, oj, aj):
		mino = min(oi,oj)
		maxo = max(oi,oj)
		for i in reversed(range(len(res_table))):
			if res_table[i]+mino < 0 \
			    or res_table[i]+maxo >= len(self.pdb.res):
				continue
			cchain = self.pdb.res[res_table[i]].chain_id
			
			ir = self.pdb.res[res_table[i]+oi]
			jr = self.pdb.res[res_table[i]+oj]
			if cchain != ir.chain_id or cchain != jr.chain_id:
				continue
			try:
				iat = [a.split()[0] for a in ir.names].index(ai)
			except ValueError:
				print "Warning! From atom %s not present in "\
				   "residue %s %d offset %d from %s"%(ai, \
				    ir.name, res_table[i]+oi+1, oi, \
				    self.pdb.res[res_table[i]].name)
				del res_table[i]
				continue
			try:
				jat = [a.split()[0] for a in jr.names].index(aj)
			except ValueError:
				# Common for variable-composition residues...
				continue
			self.add_edge(ir.atom_zero+iat, jr.atom_zero+jat)
		
	def append_edge(self, line):
		i = int(line[0])-1
		if i < 0 or i > self.atoms:
			raise inputError, "Error EDGE line contains out-of"\
				"-range from atom number %d.\n"%(i+1)
		for to in line[1:]:
			j = int(to)
			if j < 0 or j > self.atoms:
			    raise inputError, "Error EDGE line contains out-of"\
				"-range to atom number %d->%d.\n"%(i+1, j+1)
			self.append_edge(i, j)
	
	def append_redge(self, line):
		rname = line[0]
		anames = parse_aname(line[1:])
		res_table = [ r for r in range(len(self.pdb.res)) \
				if self.pdb.res[r].name == rname ]
		#print anames
		oi = anames[0][0]
		ai = anames[0][1]
		for oj, aj in anames[1:]:
			self.add_all_redge(res_table, oi,ai, oj,aj)
	
# Commit constants to memory for faster computation.
	def commit_ff(self):
		for i in self.ff_terms:
		    try:
			commit_method = getattr(self, "commit_%s"%i)
		    except AttributeError:
			continue
		    commit_method()
# Calculate energy for given configuration.
	def calc_energy(self, x):
		en = 0.0
		for i in self.ff_terms:
		    try:
			en_method = getattr(self, "%s_energy"%i)
		    except AttributeError:
			continue
		    en += en_method(x)
		return en
	
# Calculate energy and forces for given configuration.
	def calc_force(self, x):
		en = 0.0
		F = zeros((self.atoms,3), float)
		for i in self.ff_terms:
		    try:
			F_method = getattr(self, "%s_force"%i)
		    except AttributeError:
			continue
		    en_i,F_i = F_method(x)
		    en += en_i
		    F += F_i
		return en,F
	
	# Append all linear constraints into a single uber-constraint.
	def calc_constraints(self):
		if not hasattr(self, "parameters"):
			raise ProgramError, "Error! can_force_match() must "\
					"be called first!"
		
		constr = []
		i = 0
		t = 0
		for k in self.ff_terms:
		    tconst = getattr(self, "constrain_%s"%k)()
		    for const in tconst:
		        np = self.nt[t]
			if const != None:
		          for ci in const:
			    c = zeros(self.parameters, float)
			    c[i:i+np] = ci
			    constr.append(c)
			i += np
			t += 1
		assert(i == self.parameters)
		
		return array(constr)
	
	# Produces an uber-matrix from the quadratic prior penalty functions for
	# each interaction type.
	def calc_prior(self):
		if not hasattr(self, "parameters"):
			raise ProgramError, "Error! can_force_match() must "\
					"be called first!"
		prior = []
		prank = [] # Corresponding rank of matrix.
		prng = [] # Corresponding range in theta.
		i = 0
		t = 0
		for k in self.ff_terms:
		    info = getattr(self, "%s_info"%k)
		    for Ki in info:
			np = self.nt[t]
			prior += Ki.prior
			prank += Ki.prank
			prng += [[i,i+np]]*len(Ki.prior)
			i += np
			t += 1
		assert len(prior) == len(prng)
		
		return prior, array(prank), prng
	
	# Append all design matrices into a single uber-design matrix.
	def calc_design(self, x, order=0):
		if not hasattr(self, "parameters"):
			raise ProgramError, "Error! can_force_match() must "\
					"be called first!"
		if order == 0:
			D = zeros(x.shape[:-2]+(1,self.parameters), float)
		elif order == 1:
			D = zeros(x.shape+(self.parameters,), float)
		else:
		    raise ProgramError, "Error order %d not supported."%order
		i = 0
		for k in self.ff_terms:
		    design = getattr(self, "design_%s"%k)(x, order)
		    if order == 1:
				design = design[1]
		    for d in design:
			np = d.shape[-1]
			D[...,i:i+np] = d
			i += np
		assert(i == self.parameters)
		
		return D

# Container class for bspline-related calculations.
# Defines several important variables:
# type_name - name of the bond type used here
# bspl - the underlying bspline object
# (x0,h,n), periodic - spline spacing, number, and periodic flag.
# constr - the function returning the bspline coeff. constraints
# m - the power of x to include in the prior array.
# c - spline constants
# constraint - ? x n array defining null space of c
class type_calc:
	def __init__(self, type_name, f, constr, m=0):
		self.type_name = type_name
		self.f = f
		self.constraint = constr
		self.m = m
		self.calc_attributes()
	
	# Read in a new spline for ourselves.
	def read(self, base):
		self.f = read_spline_func(base+"%s.espl"%self.type_name, \
								self.type_name)
		self.calc_attributes()
	
	def commit(self):
		self.f.commit()
	
	def calc_attributes(self):
		# Normalize by range to level the playing field for
		# different types of energy contributions.
		scale = self.f.rng[1]**(self.m+1)-self.f.rng[0]**(self.m+1)
		scale = (self.m+1.0)/scale
		order = [2]
		#print "Term %s: %e"%(self.type_name,scale)
		#self.prior = [ self.f.quadratic_integral(1,self.m)*scale,
		#		self.f.quadratic_integral(2,self.m)*scale,
		#		self.f.quadratic_integral(3,self.m)*scale ]
		self.prior = [self.f.quadratic_integral(n,self.m)*scale \
					for n in order]
		self.prank = [self.f.n-n for n in order]
		#print self.prior[0]
	
	def write(self, base):
		self.f.write_spl(base+self.type_name+".espl")
	
	# Returns "order+1" arrays giving up to the order-th derivative of f(x).
	# X can be any shape, thus each return array is x.shape+(n,)
	def spline(self, x, order=0):
		if order == 0:
			spl = self.f.spline(reshape(x,x.size), None)
			return reshape(spl, x.shape+(spl.shape[-1],))
		else:
			spl = self.f.spline(reshape(x,x.size), order)
			return reshape( spl, \
				(spl.shape[0],) + x.shape + (spl.shape[-1],) )

# Parses a list of names and optional offsets into a list
# of the form [(off, name), (off, name), ...]
def parse_aname(tok):
	resname = []
	i = 0
	off = 0 # Default.
	while i < len(tok):
		if tok[i][0] in "+-":
			off = int(tok[i]) # Read an optional offset.
			i += 1
		else:
			off = 0
		resname.append((off, tok[i]))
		i += 1
	
	return resname

# Modular product set of input sets.
# Very useful for enumerating those pesky angles/torsions...
def modprod(*a):
	b = [ [(i,) for i in a[0]] ]
	for d in range(1, len(a)):
	    b.append([])
	    for j in a[d]:
	        b[d] += [i+(j,) for i in b[d-1]]
	c = set(b[-1])
	del b
	return c

# Daniel Brodie's Class Extension Example from aspn.activestate.com 2005/04/29.
class ExtendInplace(type):
    def __new__(self, name, bases, dict):
        prevclass = globals()[name]
        del dict['__module__']
        del dict['__metaclass__']

        # We can't use prevclass.__dict__.update since __dict__
        # isn't a real dict
        for k,v in dict.iteritems():
            setattr(prevclass, k, v)
        return prevclass
