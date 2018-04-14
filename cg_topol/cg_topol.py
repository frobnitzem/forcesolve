from pdb import build_pdb
from ucgrad import parseparam

from bonds import SplineBond
from torsions import SplineTorsion
from angles import SplineAngle
from pairs import SplinePair
from concat_term import FFconcat

from poly_term import read_poly_term
from pbonds import PolyBond, PolyUB
from ptorsions import PolyTorsion
from pangles import PolyAngle
from ljpairs import LJPair
from pimprop import PolyImprop

from read_term import read_terms

import os

RequiredParams = ["EXPDB", "TERMS"]
def cg_topol(pfile):
    params = parseparam(pfile, RequiredParams)
    if params == -1:
        raise InputError, "Invald topology file."

    pdb = build_pdb(params)
    top = read_terms(pdb, params['TERMS'][0][0])

    return top, pdb

######## These lines are a bit of a hack ########
# since they assume more than the published CgTopol API...
def show_index(topol):
    params = 0
    print "\t     name        params"
    def show_term(t):
        if hasattr(t, "terms"):
            map(show_term, t.terms)
            return

        print "%-16s: %5d"%(t.name, t.params)
    show_term(topol)
    print "%d total parameters." % topol.params
    print "%d total constraints." % len(topol.constraints)
    print "%d total inequalities." % len(topol.ineqs)
    print

# Traverse a tree of topology terms.
def write_topol_r(t, pre, c):
    if hasattr(t, "terms"):
        i = 0
        for ti in t.terms:
            ip = i + ti.params
            write_topol_r(ti, pre, c[i:ip])
            i = ip
    else:
        t.write(pre, c)

def write_topol(t, pre, c):
    ensure_dir(pre)
    write_topol_r(t, os.path.join(pre, ""), c)

def ensure_dir(f):
    d = f #os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

