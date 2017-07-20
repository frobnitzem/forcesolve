#!/usr/bin/env python

import sys, os
dn = os.path.dirname
sys.path.append(dn(dn(os.path.abspath(__file__))))

from psf import read_psf
from cg_topol import *
from numpy import array, sqrt, sum, argmax, newaxis, abs, dot, pi, arange
from read_prm import PRM

def main(argv):
    assert len(argv) == 3, "Usage: %s <param dir> <out.prm>"%argv[0]
    prm = prm_of_param(argv[1])
    prm.write(argv[2])

# returns [(n, K, delta)]
# where the list contains only parameters with |K| > tol
def charmm_tor(c, tol=1e-4):
    # W . [cos(n phi)] = [(cos phi)^n]
    W = array(
          [[ 1.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.,  0.,  0.,  0.],
	   [ 2.,  0.,  1.,  0.,  0.,  0.,  0.],
	   [ 0.,  3.,  0.,  1.,  0.,  0.,  0.],
	   [ 6.,  0.,  4.,  0.,  1.,  0.,  0.],
           [ 0., 10.,  0.,  5.,  0.,  1.,  0.],
           [20.,  0., 15.,  0.,  6.,  0.,  1.]])
    W[:,1:] *= 2.0
    W /= (2**arange(len(W)))[:,newaxis]

    # list of K values
    #K = dot(W.transpose(), c)
    K = dot(c, W)

    out = []
    for n,k in enumerate(K):
	if abs(k) < tol:
	    continue
	if k < 0.0: # cos(n phi + 180) = -cos(n phi)
	    out.append((n, -k, 180.0))
	else:
	    out.append((n, k, 0.0))
    return out

def get_term(name):
    id, c = read_poly_term(name)
    return tuple(id[1].split("-")), c

def prm_of_param(path):
    # Does not output an atoms section.
    #atoms = set(psf.t)
    #atoms = dict([(n, (i+1, psf.m[psf.t.index(n)])) \
#		  for i,n in enumerate(sorted(atoms))])
    atoms = {}
    bonds = {}
    angles = {}
    dihedrals = {}
    impropers = {}
    nonbonded = {}
    to_deg = 180./pi
    for fname in os.listdir(path):
	tp = fname.split("_")[0]
	name = os.path.join(path, fname)
	if tp == "pbond":
	    id, c = get_term(name)
	    bonds[id] = (c[2], -0.5*c[1]/c[2])
	elif tp == "pangle":
	    id, c = get_term(name)
	    if angles.has_key(id):
		r = angles[id][2:]
	    else:
		r = ()
	    angles[id] = (c[2], -to_deg*0.5*c[1]/c[2]) + r
	elif tp == "pub":
	    id, c = get_term(name)
	    if angles.has_key(id):
		r = angles[id]
	    else:
		r = 0.0, 100.0
	    angles[id] = r + (c[2], -0.5*c[1]/c[2])
	elif tp == "ptor":
	    id, c = get_term(name)
	    dihedrals[id] = charmm_tor(c)
	elif tp == "ljpair":
	    # Handles 2 cases:
	      # pair_terms name = "%d+%s-%s"
	      # pair_n_terms name = "1,%d-%s-%s"
	    id, c = get_term(name)
	    c6, c12 = c[1:]
	    eps, R0 = 0.25*c6*c6/c12, (-2*c12/c6)**(1/6.0)
	    if "+" in id[0]: # fix "4+type_type" notation
		num = id[0].split('+')
		id = (num[1],) + id[1:]
		if num[0] == '4':
		    nonbonded[id] = [eps, R0]
		elif num[0] == '5':
		    if nonbonded.has_key(id):
			nonbonded[id][0:2] = [eps, R0]
		    else:
			nonbonded[id] = [eps, R0, nan, nan]
		else:
		    raise ValueError, "LJ pairs are not 4+ or 5+?"
	    else:
		if id[0][0:3] != "1,4":
		    raise ValueError, "Special LJ pairs not 1,4?"
		id = (id[0][3:],) + id[1:]
		if nonbonded.has_key(id):
		    nonbonded[id][2:] = eps, R0
		else:
		    nonbonded[id] = [nan, nan, eps, R0]
	elif tp == "pimprop":
	    # cg_topol/pimprop.py:117
	    # just writes the line, "#IMPR <name> <K>"
	    line = open(name).read()[0].split()
	    id = line[1].split("_")[1].split("-")
	    impropers[tuple(id)] = float(line[2])
    return PRM(atoms, bonds, angles, dihedrals, impropers, nonbonded)

if __name__=="__main__":
    main(sys.argv)

