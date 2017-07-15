#!/usr/bin/env python
# Driver program for processing CHARMM-specific
# information into PDB and cg_topol data structures
# used by frc_solve.

import sys, os, argparse
dn = os.path.dirname
sys.path.append(dn(dn(os.path.abspath(__file__))))

from numpy import load, newaxis, sum, zeros
from psf import read_psf
from frc_match import *
from cg_topol import *
from read_prm import read_prm

from ewsum import ES_seed, ES_frc, dES_frc

# Re-format psf data into PDB : {
#     conn : [Set Int], -- complete connection table,
#     edge : Set (Int, Int), -- unique edges,
#     x : Array Float (n, 3), -- example / reference config.
#     names : [(name : String, res : String, atype : String)],
#     mass : Array Float n,
#     atoms : Int, -- number of atoms
#     L : None | Array Float (3,3),
# }
def pdb_of_psf(psf):
    names = [ (psf.names[i], psf.res[i], psf.t[i]) \
                for i in range(psf.N) ]
    return PDB(names, psf.m, zeros((psf.N,3)),
               set([(i,j) for i,j in psf.bonds]))

# Note: LJPair (as called here) currently uses an LJ-cutoff of 11 Ang.
def topol_of_pdb(pdb, dihedrals, UB=True, LJ=True):
    # Create constraints to fix n (looked up from dihedrals).
    #   If no n, no torsion!
    def constrain_n(type):
	t = tuple(type.split("-"))
	if dihedrals.has_key(t):
	    return [u[0] for u in dihedrals[t]]
	print "Warning: no dihedral n found for type:", type
	print("Assuming n = 3!")
	return [3]
    def mk_polytors(*a, **b):
	return PolyTorsion(*a, constrain_n=constrain_n, **b)

    terms = []
    terms.append(bond_terms(pdb, PolyBond))
    terms.append(angle_terms(pdb, PolyAngle))
    if UB:
	terms.append(angle_terms(pdb, PolyUB))
    terms.append(torsion_terms(pdb, mk_polytors))
    terms.append(improper_terms(pdb, PolyImprop))
    if LJ == "14":
	# Achtung! Because it creates pdb.pair, this must be called first.
	terms.append(pair_terms(pdb, LJPair, n=5))
	terms.append(pair_n_terms(pdb, LJPair, n=4)) # 1-4 appends to pdb.pair
    elif LJ:
	terms.append(pair_terms(pdb, LJPair, n=4))
    else:
	pair_terms(pdb, LJPair, n=4) # still need to create pairlist
    return FFconcat(terms)

# Creates MQ by giving 1 DOF to ea. unique charge/atom type in the input.
# Then removes one by setting sum = 0.
def wrassle_es(pairs, q, t):
    mask = []
    for i in range(len(q)-1):
	for j in range(i+1, len(q)):
	    if (i,j) not in pairs:
		mask.append((i,j))

    # Create an extended type bundling charge and type name.
    ext_type = [(qi,ti) for qi,ti in zip(q,t)]
    un = list(set(ext_type))
    mult = array([ext_type.count(u) for u in un])
    MQ = zeros((len(q),len(un)-1))
    for i in range(len(q)):
	j = un.index(ext_type[i])
	if j == len(un)-1: # short straw
	    MQ[i] = -mult[:-1]/float(mult[-1])
	else:
	    MQ[i,j] = 1.0

    print "Using %d minus 1 charge group types:\n"%len(mult)
    print "\n".join("%s %f (%d)"%(ti,qi,m) \
		    for (qi,ti),m in zip(un, mult))
    #print "MQ = " + str(MQ)

    un = array([i for i,j in un[:-1]]) # strip type info.

    return un, mask, MQ

def save_charges(psf, q, out):
    open(os.path.join(out, "charges"), 'w').write(
	  "# atom type old_chg new_chg\n" \
        + "".join("%d %s %f %f\n"%(i+1, psf.t[i], psf.q[i], q[i]) \
    	  for i in range(len(q))) )
    open(os.path.join(out, "charges.str"), 'w').write(
	 "".join("scalar charge set %f sele bynu %d end ! %s %f\n" % (\
		    q[i], i+1, psf.t[i], psf.q[i]) \
	 for i in range(len(q))) )

def LJ_frc(pairs, tors, x, coef):
    print coef
    print len(tors)
    f = zeros(x.shape)
    LJ14=set()
    for i,j,k,l in tors:
        if i<l:
           LJ14.add((i,l))
        else:
           LJ14.add((l,i))
    for i,j in pairs:
        r = x[:,i] - x[:,j]
        drij = (sum(r*r, -1)**(-0.5))
        if (i,j) in LJ14:
           eps = (coef[i][2]*coef[j][2])**(0.5)
           rmin = coef[i][3] + coef[j][3]
        else:
           eps = (coef[i][0]*coef[j][0])**(0.5)
           rmin = coef[i][1] + coef[j][1]
        sij = rmin*drij
        r *= (12*eps*(rmin**(-2))*(sij**(14) - sij**(8)))[...,newaxis]
        f[:,i] += r
        f[:,j] -= r
    return f

def main(argv):
    parser = argparse.ArgumentParser(description='Match Forces Using PSF.')
    parser.add_argument('psf', metavar='sys.psf', type=str,
			help='System PSF file.')
    parser.add_argument('xf', metavar='xf.npy', type=str,
			help='Coordinate and force data.')
    parser.add_argument('out', metavar='out_dir', type=str,
			help='Output directory name.')
    parser.add_argument('--prm', type=str,
		        help='Read dihedral multiplicities from prm.')
    parser.add_argument('--box', metavar=('Lx', 'Ly', 'Lz'),
			type=float, nargs=3, default=None,
		        help='Box diagonal lengths.')
    parser.add_argument('--tilt', metavar=('Lxy', 'Lxz', 'Lyz'),
			type=float, nargs=3, default=None,
		        help='Box tilt factors.')
    # Toggles.
    parser.add_argument('--noUB', dest='UB', action='store_false',
		        help='Don\'t Fit Urey-Bradley Terms.')
    parser.add_argument('--noLJ', dest='LJ', action='store_false',
		        help='Don\'t Fit LJ parameters at all.')
    parser.add_argument('--14', dest='LJ', action='store_const',
			const="14",
		        help='Fit LJ and Special 14 LJ parameters.')
    parser.add_argument('--chg', action='store_true',
		        help='Fit charges.')
    args = parser.parse_args()
    print(args)
    #exit(0)

    if args.prm != None:
	prm = read_prm(args.prm)
	dih = prm.dihedrals
    else:
	prm = None
	dih = {}
    out = args.out

    # Read input data.
    psf = read_psf(args.psf)
    pdb = pdb_of_psf(psf)
    xf = load(args.xf)
    if args.box != None:
	pdb.L = diag(args.box)
	if args.tilt != None:
	    pdb.L[1,0] = args.tilt[0]
	    pdb.L[2,0] = args.tilt[1]
	    pdb.L[2,1] = args.tilt[2]
    else:
	pdb.L = None

    # Create topol and FM object.
    topol = topol_of_pdb(pdb, dih, args.UB, args.LJ)
    if args.LJ == False:
	if prm == None:
	    raise LookupError, "PRM required for subtracting LJ (using --noLJ)"
	xf[:,1] -= LJ_frc(pdb.pair, pdb.tors, xf[:,0], \
			  [prm.nonbonded[t] for t in psf.t])

    q, mask, MQ = wrassle_es(pdb.pair, psf.q, psf.t)

    forces = frc_match(topol, pdb, 1.0, 1.0, do_nonlin=args.chg)
    forces.add_nonlin("es", q, ES_seed, ES_frc, dES_frc, (mask, pdb.pair, MQ, pdb.L))
    show_index(forces.topol)

    # Do work.
    forces.append(xf[:,0], xf[:,1])
    forces.dimensionality() # double-checks well-formedness
    forces.maximize()
    forces.write_out(out)
    if args.chg:
	q = dot(MQ, forces.nonlin["es"][0])
	print "charges = " + str(q)
	save_charges(psf, q, out)

if __name__=="__main__":
    main(sys.argv)

