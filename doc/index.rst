.. ForceSolve documentation master file, created by
   sphinx-quickstart on Thu Apr 12 14:59:40 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ForceSolve User Guide
======================================

Contents:

.. toctree::
   :maxdepth: 2

   running forcesolve
   term files
   dynamics
   analyzing results

Tutorial
--------

ForceSolve is a numerical code for performing
Bayesian inference on molecular mechanics parameters
from data on molecular forces.

What it does:

  * Lists terms based on the molecule's bond list.
  * Creates topologies capable of computing forces (when given parameters).
  * Computes Ewald electrostatic sums.
  * Computes the log-likelihood of a given parameter set (when given force data).
  * Finds parameters maximizing or sampling the log-likelihood.
  * Can even fit charges (nonlinear terms).
  * Works with general B-spline functions.
  * Simulates MD trajectories (slow).

What it does not do:

  * Run MD or compute forces without any provided parameters.
  * Compute parameters without any provided (position,force) data.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

