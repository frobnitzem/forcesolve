Term Files
==========

ForceSolve understands a molecular mechanical forcefield
as a list of additive terms, `t`.
Each term has a vector, :math:`\theta_t` of `N` parameters,
and a function, :math:`U_t`, which outputs `N` numbers.
To keep things simple, the function only takes
the positions of 2, 3, or 4 atoms.  Which atoms
are actually used is given by a permutation, :math:`\pi`.
Terms are used multiple times, so each term has a corresponding
set, :math:`\Gamma_t` of atom permutations.

All of this can be summed up by the equations,

.. math::

  U(x; \theta) = \sum_{t} \sum_{\pi \in \Gamma_t} \theta_t \cdot U_t(\pi \circ x)

  \Gamma_t = \{ (i,j,\ldots), i = 1, 2, \ldots \}
  
  U_t(x_i, x_j, \ldots) \in \mathbb R^N

  \theta_t \in \mathbb R^N

Specifying a forcefield requires listing which atoms (permutations)
go with each :math:`U_t`.  This is the purpose of a term file.
A term file is an ordinary python file whose job is to create
a list of Term objects::

  terms = [
    Term( "bond",     PolyBond,    Conn(1,2) ),
    Term( "angle",    PolyAngle,   Conn(1,2,3) ),
    Term( "torsion",  PolyTorsion, Conn(1,2,3,4) ),
    Term( "improper", PolyImprop,  OOP() ),
    PairTerm("ljpair_4+",  LJPair, Conn(1,2) | Conn(1,3))
  ]

The example term file above illustrates the basic idea.
Each `Term` is created by assigning an arbitrary (but unique) name,
selecting the term type, and providing an atom selection.
Atom selections are documented in the next section.
There is a special `PairTerm` object that takes an
*exclusion* selection list instead of an inclusion selection.
For `PairTerm`-s, every (non-self) pair of atoms is included except
those that match the selection.
The list of implemented term types includes:

============= ===== ===========================================
Type          Atoms
============= ===== ===========================================
SplineBond     2     spline-based bond
SplineAngle    3     spline-based angle
SplineTorsion  4     spline-based torsion
SplinePair     2     spline-based nonbonded
PolyBond       2     quadratic bond term
PolyAngle      3     quadratic angle term
PolyUB         3     quadratic Urey-Bradley
PolyImprop     4     quadratic improper dihedral
PolyTorsion    4     six-parameter Fourier dihedral expansion
LJPair         2     2-parameter LJ term
============= ===== ===========================================

Note by the example that the `Pair` terms are special
in that they list excluded atoms rather than included
ones.  By default, `LJPair` and `SplinePair` compute
pairwise interactions between all atoms *except* those
pairs it has been told about.  Self-interactions (i.e. 1,1)
are always excluded.

Selections
----------

The simplest selection is `Conn`, which
uses the molecule's bond connectivity to generate
all valid sequences of connected atoms.
The arguments to `Conn` must start with 1 (the central atom).
Selections can be combined with `and` and `or`
in order to filter out undesired terms and/or
add extra terms.  A list of selection primitives follows::

  Conn(1,...)

Select connected chains of atoms, counting from 1.
Up to 4 neighbors are supported.  Re-ordering is
allowed, for example `(1,3)` and `(3,1)` both give
an equivalent selection of pairs of end-atoms
involved in an angle.::

  a & b (aka a * b)
  a | b (aka a + b)
  a - b (aka a / b)

Selections can be combined using `&` (for set intersection)
and `|` for (set union).  The operators `*` and
`+` are convenient synonyms for `&` and `|`.  However,
python will not allow using the words `and` / `or` / `not`.
For set exclusion, use `/` or `-` instead.
Be careful of operator precedence (Python's precedence is used).::

  Id(i,j,k,...)

Atom ID selection (sequential atom numbers start at 1).
A value of "None" indicates a wildcard match.::

  Type(ti,tj,...)

Type selection (all atoms of a given type).
A value of "None" indicates a wildcard match.::

  TFn(test)

A general type-based test function, `test(ti,tj,...)`.
Like filter, it should return True for all keepers.::

  IdFn(test)

A general atom id-based test function, `test(i,j,...)`
Like filter, it should return True for all keepers.

Each atom may be involved in any number of terms, and no attempt is
made to remove possible duplication.  This means that it
is possible to accidentally add a 1,4 atom pair both
to a special 1,4 pair term and to a "generic" all-pairs term.
Only judicious use of selections will prevent this -- you have been warned.

For writing `TFn` or `IdFn` tests, it's important to know that
atoms within a term are always ordered "canonically"
based on their type.  For a 2-center interaction, `ti <= tj`
(where the comparison is alphabetic).  For a 3-center one,
`ti <= tk`.  For a 4-center one, `tj <= tk`.  For out-of-plane
interactions, the out-of-plane atom is first and the rest
are sorted except in one case.  If two of the non-central
out-of-plane atoms have equal types, they immediately follow
`i` as `j` and `k`.
All of this means that `IdFn` tests could see either
ordering of the atoms id-s, while `TFn` will always see a nice, ordered
list of types.

