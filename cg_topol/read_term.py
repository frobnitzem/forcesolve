# Read term.py files
# TODO: - add Ewald term

#import imp
from edge import modprod, srt2
from concat_term import FFconcat

# Basic term types
from bonds import SplineBond
from torsions import SplineTorsion
from angles import SplineAngle
from pairs import SplinePair

from pbonds import PolyBond, PolyUB
from ptorsions import PolyTorsion
from pangles import PolyAngle
from ljpairs import LJPair

from pimprop import improper_terms, PolyImprop

# (gen or filter) and (gen or filter)
# = (gen and gen) or (gen and filter) or (filter and gen) or (filter and filter)
def mkAnd(a, b):
    def gen(pdb):
        if a.gen != None:
            if b.gen != None:
                # have to explicitly create one of the sets
                s1 = set(a.gen(pdb))
                for i in b.gen(pdb):
                    if i in s1:
                        yield i
                    s1.remove(i)
            else:
                s1 = a.gen(pdb)
            if b.test != None:
                f = b.test(pdb)
                for i,t in s1:
                    if f(i,t):
                        yield i,t
        if b.gen != None and a.test != None:
            f = a.test(pdb)
            for i,t in b.gen(pdb):
                if f(i,t):
                    yield i,t
    if a.gen == None and b.gen == None:
        gen = None

    def test(pdb):
        if a.test != None and b.test != None:
            f = a.test(pdb)
            g = b.test(pdb)
        return lambda i,t: f(i,t) and g(i,t)
    if a.test == None or b.test == None:
        test = None

    return Selection(gen, test)

# (gen or filter) or (gen or filter)
# = (gen or gen) or (filter or filter)
def mkOr(a, b):
    def gen(pdb):
        if a.gen != None:
            for i in a.gen(pdb):
                yield i
        if b.gen != None:
            for i in b.gen(pdb):
                yield i
    if a.gen == None and b.gen == None:
        gen = None
    
    def test(pdb):
        if a.test != None and b.test != None:
            f = a.test(pdb)
            g = b.test(pdb)
        return lambda i,t: f(i,t) or g(i,t)
    if a.test == None and b.test == None:
        test = None
    elif a.test == None:
        test = b.test
    elif b.test == None:
        test = a.test
    return Selection(gen, test)

# not (gen or filter) = not gen and not filter
def mkNot(a, pdb):
    def test(pdb):
        if a.gen != None: # have to explicitly create...
            s1 = set(a.gen(pdb))
        else:
            s1 = set()
        if a.test != None:
            g = a.test(pdb)
        else: g = lambda i,t: False
        def f(i,t):
            if i in s1:
                return False
            return not g(i,t)
    return Selection(None, test)

""" Lazy selection code -- implements an "or" of a generator
    or a filter.
"""
# This is meant to be eventually used using set( Selection(...).run(pdb) ).
# That's beacause 'run' is a generator that may return duplicates.
class Selection:
    def __init__(self, gen, test=None):
        self.gen = gen
        self.test = test
    def __and__(a, b):
        return mkAnd(a, b)
    def __mul__(a, b):
        return mkAnd(a, b)
    def __or__(a, b):
        return mkOr(a, b)
    def __add__(a, b):
        return mkOr(a, b)
    def __div__(a,b):
        return mkAnd(a, mkNot(b))
    def __sub__(a,b):
        return mkAnd(a, mkNot(b))
    def run(self, pdb):
        if self.test != None:
            raise RuntimeError, "Improper selection."
        # Run generator on the pdb.
        for i in self.gen(pdb):
            yield i

# Filters return only True/False
def IdFn(u):
    return Selection(None, lambda pdb: (lambda ij,tij: u(*ij)))

def TFn(u):
    return Selection(None, lambda pdb: (lambda ij,tij: u(*tij)))

def Id(*ij):
    def test(n,t): # as a filter
        assert len(n) == len(ij), "Atom number mismatch"
        for i,j in zip(n, ij):
            if i != j and j != None:
                return False
        return True
    if ij.count(None) > 0: # 'Id' contains a wildcard
        return Selection(None, lambda pdb: test)
    def gen(pdb): # generate a literal, singlet list
        yield lookup(map(pred, ij), pdb)
    return Selection(gen, None)

def Type(*tij):
    def test(*t):
        assert len(t) == len(tij), "Atom number mismatch"
        for i,j in zip(t, tij):
            if i != j and j != None:
                return False
        return True
    return TFn(test)

# generate all 1,2,3,4 pairs
def tors(pdb):
    tors = set()
    for j,k in pdb.edge:
        tors |= set([(i,j,k,l) for i,l in modprod(pdb.conn[j]-set([k]), \
                                                  pdb.conn[k]-set([j])) \
                                           if i != l])
    return tors

def oop(pdb):
    tors = set()
    for i in range(pdb.atoms):
        if len(pdb.conn[i]) == 3:
            tors.add((i,) + tuple(pdb.conn[i]))
    return tors

def pred(x):
    return x-1
def succ(x):
    return x+1
# return pairs of 1-based index list, type list
def lookup(ij, pdb):
    ijk = tuple(map(succ, ij))
    typ = tuple( map(lambda i: pdb.names[i][2], ij) )
    return order(ijk, typ)

def Conn(*ijk):
    ijk = tuple(ijk)
    M = max(ijk)
    assert min(ijk) == 1 and M <= 4, "Invalid connectivity: %s" % str(ijk)
    def gen(pdb):
        if M == 2:
            v = pdb.edge
        elif M == 3:
            v = []
            for j in range(pdb.atoms):
                v += [(i,j,k) for i,k in modprod(pdb.conn[j], \
                                pdb.conn[j]) if i < k]
        elif M == 4:
            v = tors(pdb)
        for i in v:
            ij = [i[u-1] for u in ijk] # user's permutation (1,..,4)
            yield lookup(ij, pdb)
    return Selection(gen, None)

# Generate all out-of-plane centers.
def OOP():
    def gen(pdb):
      for i,j,k,l in oop(pdb):
        ti = pdb.names[i][2]
        tj = pdb.names[j][2]
        tk = pdb.names[k][2]
        tl = pdb.names[l][2]
        if tk > tl: # Bubble sort
            tk, tl = tl, tk
            k, l   = l, k
        if tj > tk:
            tj, tk = tk, tj
            j, k   = k, j
            if tk > tl: # Bubble sort
                tk, tl = tl, tk
                k,   l = l, k
        if tk == tl: # according to the CHARMM docs, we swap again!
            tj, tk, tl = tk, tl, tj
            j, k, l    = k, l, j

        yield (i+1,j+1,k+1,l+1), (ti,tj,tk,tl)
    return Selection(gen, None)

# canonically order term (alpha by type)
def order(i, t):
    if len(i) == 2 and t[0] > t[1]:
        i = i[1],i[0]
        t = t[1],t[0]
    elif len(i) == 3 and t[0] > t[2]:
        i = i[2],i[1],i[0]
        t = t[2],t[1],t[0]
    elif len(i) == 4 and ( t[1] > t[2] or \
                           t[1] == t[2] and t[0] > t[3] ):
        i = i[3],i[2],i[1],i[0]
        t = t[3],t[2],t[1],t[0]
    return i, t

# Symbolic term -- knows how to get from a PDB to an actual term.
class Term:
    def __init__(self, name, cl, gen):
        self.name = name # e.g. "bond"
        self.cl = cl # e.g. PolyBond
        self.gen = gen # e.g. Conn(1,2)

    # Group all terms by types of the atoms involved.
    def build(self, pdb):
        index = {}
        for ijk, typ in self.gen.run(pdb):
            name = "_".join(typ)
            if not index.has_key(name):
                index[name] = set()
            index[name].add(tuple(map(pred, ijk)))
        print( "%d/%d unique %s terms"%( len(index), \
                sum(map(len, index.values())), self.name) )
        terms = [self.cl(self.name+"_"+n,l) for n,l in index.iteritems()]
        return FFconcat(terms)

class PairTerm:
    def __init__(self, name, cl, gen):
        self.name = name # e.g. "pair_4+"
        self.cl = cl # e.g. SplinePair
        self.gen = gen # e.g. Conn(1,2)

    # Group all terms by types of the atoms involved.
    def build(self, pdb):
        unt = list( set( name[2] for name in pdb.names ) )
        unt.sort()
        ant = dict((t,[]) for t in unt) # atoms with type 't'
        aex = {} # exclusions with type 'ti,tj'
        for i,name in enumerate(pdb.names):
            ant[name[2]].append(i)
        for i in range(len(unt)):
            for j in range(i, len(unt)):
                aex[(unt[i], unt[j])] = set()

        for ij, typ in self.gen.run(pdb):
            aex[typ].add( srt2(ij[0]-1, ij[1]-1) )

        index = {}
        count = 0
        for i,t in enumerate(unt):
            for j in range(i, len(unt)): # all pairs of types
                u  = unt[j]
                ex = aex[(t,u)]
                aa = [k for k in ant[t] \
                        if not all(srt2(k,z) in ex for z in ant[u] \
                                                    if z != k)]
                bb = [k for k in ant[u] \
                        if not all(srt2(z,k) in ex for z in ant[t] \
                                                    if z != k)]
                nex = set(aa + bb)
                for i,j in ex:
                    if i not in nex or j not in nex:
                        ex.remove((i,j))
                if len(aa) > 0 and len(bb) > 0:
                    index["%s_%s"%(t,u)] = ex, aa, bb
                    if i == j:
                        count += len(aa)*(len(aa)-1)/2
                    else:
                        count += len(aa)*len(bb)
                    count -= len(ex)
        print( "%d/%d unique %s terms"%(
                len(index), count, self.name) )
        terms = [ self.cl(self.name+"_"+n, l[0], excl=l[1:]) \
                    for n,l in index.iteritems() ]

        return FFconcat(terms)

# Read and interpret the term file.
def read_terms(pdb, name):
    #mod = imp.load_source("TermFile", name)
    #if not hasattr(mod, 'terms'):
    #    raise KeyError, "Term File: %s does not define 'terms'"%(name)
    #terms = mod.terms

    exec(open(name).read())
    if not locals().has_key('terms'):
        raise KeyError, "Term File: %s does not define 'terms'"%(name)

    return FFconcat( [t.build(pdb) for t in terms] )

# 2-atom graph: 1-2 (a)
# 3-atom graph: 1-2-3 (a) | 1-2-3..1(b) [a is a subgraph of b]
#
# 4-atom graphs: six possibilities
#
# 1-2-3-4 (a)
#
# 3   2 (b)
#  \ / 
#   1
#   |
#   4
#
# a is a subgraph of c-f.
# 
# 3-2 (c)
# | |
# 4-1
#      
# 1-2-3 (d)
#   |/
#   4
#
#   2 (e)
#  / \
# 3---1
#  \ /
#   4
#
# tetrahedral (f)
#   3-2
#   |X|
#   4-1

if __name__=="__main__":
    from pdb import PDB
    pdb = PDB([ ("1",1,'H3'),
                ("2",1,'H3'),
                ("3",1,'H3'),
                ("4",1,'C3'),
                ("5",1,'H2'),
                ("6",1,'H2'),
                ("7",1,'C2'),
                ("8",1,'H2'),
                ("9",1,'H2'),
                ("10",1,'C2'),
                ("11",1,'H3'),
                ("12",1,'H3'),
                ("13",1,'H3'),
                ("14",1,'C3') ],
               [], [], edge=[(0,3), (1,3), (2,3), (3,6),
                             (4,6), (5,6), (6,9),
                             (7,9), (8,9), (9,13),
                             (10,13), (11,13), (12,13)])
    #def f(i,j):
    #    return i == 4 or j == 4
    #sel = Conn(1,2) * Type('C', 'H') | Conn(1,3) * IdFn(f)
    #sel = OOP()
    #sel = Conn(1,2,3)
    #sel = Conn(1,3)
    #sel = ( Conn(1,2,3,4) | (Conn(1,3,2,1) & Id(1,None,None,1)) ) \
    #            & Type(None, 'C', 'C', None)
    #sel = Conn(1,2) | Conn(1,3) | Conn(1,4)
    #print set(sel.run(pdb))
    top = read_terms(pdb, '../share/vanilla.py')
    from cg_topol import show_index
    show_index(top)

