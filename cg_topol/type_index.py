# properties by atom type indexing services

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

# Copyright 2008 David M. Rogers.  See COPYING for copyright information.
# Dr. Thomas Beck Lab
# University of Cincinnati
# This work was supported by a DOE CSGF.

from cg_topol import *

# Adds the tindex (type index) by defining addition functions for ra-t and t-ra.
class cg_topol:
	__metaclass__=ExtendInplace
	def append_rat(self, line):
		self.add_rat((line[0],line[1]), line[2])
	def append_tra(self, line):
		t = line[0]
		for i in range(1,len(line),2):
			self.add_rat((line[i],line[i+1]), t)
	def append_tmass(self, line): # Mass for a given type.
		self.tmass[line[0]] = 1836.1527*float(line[1]) # Add a mass to the index.
	
	def add_rat(self, ra, t):
		if self.tindex.has_key(ra):
		 if self.tindex[ra] != t:
		   raise InputError, "Error! ra pair " + str(ra) + " was " +\
		    str(self.tindex[ra])+", but is redefined as type "+t+"."
		else:
			self.tindex[ra] = t

