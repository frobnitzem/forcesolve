#!/bin/sh

# Configuration coarse-graining tutorial.
# David M. Rogers
# Dr. Thomas Beck Lab
# University of Cincinnati
# 12/10/2008
# This work was sponsored by a DOE CSGF and an NSF grant.

# Commands and description for the coarsen.py tutorial test system.

# Please run these commands manually so you know what you are doing.
# also, 
echo "This is not a conventional shell script." && exit 1

# First, you will need a DNA pdb -- use one from the Nucleic Acids Database:
wget http://ndbserver.rutgers.edu/ftp/NDB/models/mdna001/mdna001.pdb

# Next, you have to define what atomic sites go to what beads, and 
# with what weighting (e.g. center of mass, etc.)
# Have a look at the pdb and the dna_site.cpar files to understand
# how to put this together.
vi dna_site.cpar

# Now you are ready to coarsen.
../../coarsen.py -d dna_site.cdef -i mdna001.pdb -o mdna.cg.pdb -m mdna.cg.mat

# Which generates a CG representation of every MODEL
# in mdna001.pdb.
# You can compare the results with mine in mdna001.cg.pdb

# Note the "Abandoning residue T 13" message due to a missing atom
# you told it to use in the CG site.  This is a hint that we didn't
# bother to correctly fix the 5' and 3' ends of the DNA!  Do this
# before you actually use this system, OK?

# It also generates a matrix meant to multiply atomic config.s by
# to get coarse coords.  This is the preferred choice for generating
# input (coord.s and forces) to frc_solve.py.
