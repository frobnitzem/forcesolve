# Generate a list of configurations and forces for the
# lowly water molecule

import sys
from ucgrad import *

# intra-molecular ff for water, kb, ka, 
# E = kb (r - r0)^2/2 + ka (th - th0)^2
#   = kb/2 r^2 - r0 kb r + kb r0^2/2
kb = 100.0
r0 = 1.1
ka = 200.0
a0 = 109.0*pi/180.0

bnd = [0.5*r0*r0, -r0*kb, 0.5*kb]
ang = [0.5*a0*a0, -a0*ka, 0.5*ka]

# v = [r0, r1, ang]
# Molecule painting function (sets the geometry).
def mk_x(v):
    s = sin(0.5*v[2])
    c = cos(0.5*v[2])
    cx = (v[0]-v[1])*s/18.0 # shift left by x-center of mass
    cz = (v[0]+v[1])*c/18.0 # shift down by z-center of mass
    x = array([[       -cx, 0.,       -cz],
               [ v[0]*s-cx, 0., v[0]*c-cz],
               [-v[1]*s-cx, 0., v[1]*c-cz]])
    return x

def dE(v):
    F = zeros((3,3))
    s = sin(0.5*v[2])
    c = cos(0.5*v[2])
    rh0 = array([ s, 0., c])
    rh1 = array([-s, 0., c])
    dEdr0 = 2*v[0]*bnd[2] + bnd[1]
    dEdr1 = 2*v[1]*bnd[2] + bnd[1]
    # z = (r0 . r1) / (|r0| |r1|)
    # d(arccos(z))/dz = -1/sqrt(1-z^2)
    dEdca = (2*v[2]*ang[2] + ang[1])*(-1./sin(v[2]))

    #dEdr0 = 0
    #dEdr1 = 0
    #dEdca = -1./sin(v[2]) # validating angle code
    # da / dx0 = d(arccos(r0 . r1)) / dx0 = darccos(z)/dz dr0/dx0
    z = cos(v[2])
    d0 = dEdr0 * rh0
    d0 += dEdca*(rh1 - z*rh0)/v[0]
    d1 = dEdr1 * rh1
    d1 += dEdca*(rh0 - z*rh1)/v[1]
    return array([-d0-d1, d0, d1])

def main(argv):
    S = 100
    x = zeros((S, 3, 3))
    F = zeros((S, 3, 3))
    for i in xrange(S):
	v = rand.standard_normal(3)*array([0.1, 0.1, 10*pi/180.0]) + array([r0,r0,a0])
	x[i] = mk_x(v)
	F[i] = -dE(v)
    write_array("test.x", x)
    write_array("test.f", F)
    return 0

# validating angle code.
def dangle(x):
	rijk_inv = 1.0/sqrt(sum(x*x,-1)) # x.shape by 2
	x *= rijk_inv[..., newaxis] # Normalize x.
	cth = sum(x[...,0,:]*x[...,1,:],-1) # x.shape vector of cosines.
	A = cth[...,newaxis,newaxis]*x # Future derivatives.
	A[...,0,:] = x[...,1,:]-A[...,0,:]
	A[...,1,:] = x[...,0,:]-A[...,1,:]
	A *= rijk_inv[..., newaxis] # d(cth) / di and dk
	return arccos(cth), -A/sqrt(1.0-cth*cth)

def test():
    v = rand.standard_normal(3)*array([0.1, 0.1, 10*pi/180.0]) \
	    + array([r0,r0,a0])
    x = mk_x(v)
    F = -dE(v)
    a, da = dangle(reshape(x[1:]-x[0], (1,2,3)))
    print a-a0, v[2]-a0
    print F
    print da

if __name__ == "__main__":
    main(sys.argv)

