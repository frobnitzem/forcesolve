# See maggotroot.blogspot.com/2013/11/constrained-linear-least-squares-in.html
# for examples using sparse matrices.
# Unfortunately, the conversions to/from cxvopt.matrix didn't work for me.

import numpy as np
from cvxopt import solvers, matrix #, spmatrix, mul
#import itertools
#from scipy import sparse
import numpy.linalg as la

def mat(x):
    if x is None:
	return None
    if len(x.shape) == 2: # gotta love fortran
	return matrix(np.transpose(x).tolist(), x.shape, 'd')
    else: # the column vector
	return matrix(x.tolist(), (len(x), 1), 'd')

def arr(m):
    if m == None:
	return None
    s = m.size
    if s[0] == 1: # vectors
	m = m.T
	s = m.size
    if s[1] == 1:
	x = np.zeros(s[0])
	for i in range(s[0]):
	    x[i] = m[i,0]
	return x
    x = np.zeros(s)
    for i in range(s[0]):
	for j in range(s[1]):
	    x[i,j] = m[i,j]
    return x

# Solve min x^T P x / 2 + q^T x
#   subject to G x <= h
#   and A x = b
# documentation at:
# http://cvxopt.org/userguide/coneprog.html#quadratic-programming
def quad_prog(P, q, G=None,h=None, A=None,b=None, x0=None):
    assert len(P) == len(q) and len(P) == P.shape[1]
    P = 0.5*(P + np.transpose(P))
    if G is None: # lsqlin will throw a fit
	if A == None: # Totally unconstrained solve.
	    return la.solve(P, -q)
	n = len(A) # number of constraints
	U, s, V = la.svd(A)
	if np.any(np.abs(s) < 1e-8):
	    print "Error: constraints are degenerate!"
	    return np.zeros(len(q))
	ib = np.dot(b, np.dot(U/s, V[:n]))
	y = la.solve(np.dot(V[n:], np.dot(P, np.transpose(V[n:]))),
		     -np.dot(V[n:], q + np.dot(P, ib)))
	return np.dot(y, V[n:]) + ib

    sol = solvers.qp(mat(P),mat(q), mat(G),mat(h), mat(A),mat(b), mat(x0))
    if sol['status'] != "optimal":
	print "Unable to minimize."
	return np.zeros(len(q))
    return arr(sol['x'])

