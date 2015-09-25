#!/usr/bin/env python2.7

import sys
from ucgrad import read_matrix
from numpy import load, newaxis, sum, zeros
from frc_match import *
from cg_topol import *

def mk_topol():
    terms = [ PolyBond("OH", [(0,1), (0,2)]), \
              PolyAngle("HOH", [(1,0,2)]) ]
    return FFconcat(terms)

class PDB:
    def __init__(self, names, mass):
        self.names = [["UNK", 0, n] for n in names]
        self.mass = array(mass)
        self.atoms = len(mass)

def main(argv):
    # Set up output locations.
    out = "mle/"

    # Read input data.
    x = read_matrix("test.x").reshape((-1,3,3))
    f = read_matrix("test.f").reshape((-1,3,3))

    # Create topol and FM object.
    topol = mk_topol()
    pdb = PDB(["O", "H", "H"], [16., 1., 1.])

    forces = frc_match(topol, pdb, 1.0, 1.0)
    show_index(forces.topol)

    # Do work.
    forces.append(x, f)
    forces.dimensionality() # double-checks well-formedness
    forces.maximize()
    forces.write_out(out)

if __name__=="__main__":
    main(sys.argv)

