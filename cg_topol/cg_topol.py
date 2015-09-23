from pdb import build_pdb
from ucgrad import parseparam

from bonds import bond_terms, SplineBond
from torsions import torsion_terms, SplineTorsion
from angles import angle_terms, SplineAngle
from pairs import pair_terms, SplinePair
from concat_term import FFconcat

from pbonds import PolyBond
from ptorsions import PolyTorsion
from pangles import PolyAngle
from ljpairs import LJPair

RequiredParams = ["EXPDB", "CONST_DIR"]
def cg_topol(pfile):
    params = parseparam(pfile, RequiredParams)
    if params == -1:
        raise InputError, "Invald topology file."

    pdb = build_pdb(params)
    terms = []
    if params.has_key("BOND"):
        terms.append(bond_terms(pdb, SplineBond))
    if params.has_key("ANGLE"):
        terms.append(angle_terms(pdb, SplineAngle))
    if params.has_key("TORSION"):
        terms.append(torsion_terms(pdb, SplineTorsion))
    if params.has_key("PAIR"):
        terms.append(pair_terms(pdb, SplinePair))

    return FFconcat(terms), pdb

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
    print

def write_topol(t, pre, c):
    if hasattr(t, "terms"):
        i = 0
        for ti in t.terms:
            ip = i + ti.params
            write_topol(ti, pre, c[i:ip])
            i = ip
    else:
        t.write(pre, c)

