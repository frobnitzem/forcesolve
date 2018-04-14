#!/bin/sh

# Water Hamiltonian inference through force-matching tutorial.
# David M. Rogers
# Dr. Thomas Beck Lab
# University of Cincinnati
# 12/7/2008
# This work was sponsored by a DOE CSGF.

# Commands and description for the water tutorial test system.

# Please run these commands manually so you know what you are doing.
# also, 
echo "This is not a conventional shell script." && exit 1

# 1. First, you will need to collect some sample data. 
# ====================================================
grompp -c start.gro -f run.mdp -o run
mdrun -s run -o run -g run -e run -c run
trjconv -f run.trr -o run.pdb -s run.tpr
# Clean up the output files a little, gromacs makes lots of them.

# Create the .x and .f files using numerical differencing (which adds some error, but its easy to show and tell.)
python <<EOF
# I suggest adding the directory with ucgrad to PYTHONPATH so this will work
from ucgrad import *

pdb, x = read_pdb("run.pdb")
x.shape

# Change x -> cg x here if you are getting rid of degrees of freedom.
# x = dot(identity(x.shape[1]), x) # Dot operates on second to last dimension.
# x = transpose(x, (1,0,2)) # Don't ask me why.

x /= Bohr # Convert to atomic units.
write_array("end.x", x[-1]) # Save final frame for later.

v = x[1:]-x[:-1]
f = v[1:]-v[:-1]

f.shape
x = x[1:-1] # fix to same times as numerical f

# Remember, atomic units.
m = resize(array([16.,1.,1.]), 648)*Pmass
dt = dt=2*pi*Hartree/Planck*1e-15 # 1 fs = 41.341414985923947 atomic time units

f *= m[newaxis,:,newaxis]/dt**2 # f = m a

write_array("run.x", x)
write_array("run.f", f)

EOF
cd .. # Back to main tutorial directory.

# 2. Force Matching
# =================

# First, make sure that your CG topology is correct.
vi topol/cg.ff # Remember to check box volume too!
vi topol/types.index
vi topol/wat.pdb # Layout of atoms matches .x
# It turns out this is the easy part!
# The trailing foward slashes after output dir.s are important.
../../frc_solve.py -p topol/cg.ff -x gmx/run.x -f gmx/run.f -mle mle/ -o out/ -dt 41.341414985923947 -ns 10 >fm.log &
# And take a nap.  Note that 10 samples is way too few,
# except for show and tell.

# Note that the output directories contain ENERGIES, which are the integrated
# (matched) forces plus an arbitrary constant.  This makes them nicer
# to work with.

# 3. Analyzing results.
# =====================
# Right away, we can plot the fitted functions (although the B-spl. constants aren't exactly the function, you get the idea, see misc. tutorial for a better one):
gnuplot -persist out.gplot

# Now, you can't actually simulate using do_md until you give your topology
# energy function parameters:
cp out/*.espl param/*.espl
mkdir sim
# Now, we can do an MD run to get properties to compare.
../../do_md.py -p param/cg.ff -x end.x -o sim/run -nt 10 -ns 5000 -dt 41.341414985923947

