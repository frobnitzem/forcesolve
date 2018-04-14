Running ForceSolve
==================

At present, the simplest way to run ForceSolve is
to use the https://github.com/frobnitzem/chemparam package,
or the webserver at https://predictivestatmech.org:300/.

There are also tutorials within the project source
that describe using ForceSolve directly.

Before running forcesolve, you need a set of system configurations
(atomic coordinates) and forces, as well as an idea of what
forcefield terms you want to fit.

To run forcesolve directly, use::

  frc_solve.py -p topol.ff -x coords.xyz -f frc.xyz -o out/

The topology is an input file, described below.

The files `coords.xyz` and `frc.xyz` should be
text files containing arrays of coordinates (ordered
as usual by frame, then atom number, then xyz).
The input file should be in atomic units (Bohr and Hartree/Bohr).
The program's help documentation also describes how
to change the units by setting a scale converting
your file to the required atomic units.

Topology file
-------------

The topology needs to define all bonds, and the mapping
from atom/residue name to atom type.
It consists of a set of statements, one per line.
Can also include other topology files,
using `include` statements.
Statement meanings are given by example below.

* `TRA hw SOL HW1 SOL HW2 SOL HW3`

  - `TRA` names all atoms with a given type (parameter 1) by giving residue-atom name pairs.

* `RAT SOL OW ow`

  - `RAT` defines a single residue-name pair that maps to a single type.

* `TMASS hw 1`

  - `TMASS` defines the mass of an atom by type.  This is used for mass-weighting the residual.

* `REDGE SOL OW HW1 HW2`

  - `REDGE` defines edges by giving a residue-atom pair (central atom) and then a list of other atoms within the residue that are also bonded.
  - Far-atoms can optionally be prefixed by a `-1` as in `REDGE SOL OW -1 OW`,
    to indicate the central atom (`OW` here) is bonded to the far atom
    in the preceding residue.

* `file EXPDB wat.pdb`

  - A single `file EXPDB` argument is required to give an example PDB file.

* `file TERMS term.py`

  - A single `TERMS` file is required to compute atom selections.
    Its format is described in the section on :doc:`term files`.

* `PERIODIC 35.18859 35.18859 35.18859`

  - The `PERIODIC` directive lets you specify a cubic periodic box size.

Run Process
------------------------

Internally, the code goes through the following steps.

1. Load an input structure.
2. Create a list of edges and atom types from the parameter file.
3. Read the term file and execute it to generate lists of bonds, angles,
   torsions, etc. from the molecule graph.
4. Read the coordinate and force data.
5. Pre-compute linear regression coefficient matrices.
6. Find maximum likelihood parameters and possible carry out MC sampling.
7. Store all parameters (possibly averages) to an output directory.

Running MD follows the same procedure up to step 3, but then
reads a previous output parameter directory
and does time-stepping.

Extending Terms
------------------------

Extending the default list of interaction types
can be done by writing a new code following
the examples in `cg_topol/pbond.py` or `cg_topol/ptorsions.py`.
Send a pull-request if you have tested that
its numerical derivative is correct (using `run_deriv.py`).

