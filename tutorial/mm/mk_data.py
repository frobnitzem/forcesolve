# Generate a list of configurations and forces for the
# lowly water molecule

import sys
from ucgrad import *

# intra-molecular ff for water, kb, ka, 
# E = kb (r - r0)^2/2 + ka (th - th0)^2
#   = kb/2 r^2 - r0 kb r + kb r0^2/2
kb = 100.0
r0 = 1.0
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
    # da / dx0 = d(arccos(r0 . r1)) / dx0 = darccos(z)/dz dr0/dx0
    d0 = dEdr0 * rh0
    d0 += dEdca*(rh1 / v[0] - cos(v[2]) * rh0 / v[0]**2)
    d1 = dEdr1 * rh1
    d1 += dEdca*(rh0 / v[1] - cos(v[2]) * rh1 / v[1]**2)
    return array([-d0-d1, d0, d1])

# dE/ds = sum dE/dx dx/ds
def prj_de(v, x, de):
    s = sin(0.5*v[2])
    c = cos(0.5*v[2])
    cxa = s/18.0; cxb = -s/18.0; cxc =  (v[0]-v[1])*c/18.0
    cza = c/18.0; czb =  c/18.0; czc = -(v[0]+v[1])*s/18.0
    dx = array([[[       -cxa, 0.,       -cza],
                 [      s-cxa, 0.,      c-cza],
                 [       -cxa, 0.,       -cza]],
                [[       -cxb, 0.,       -czb],
                 [       -cxb, 0.,       -czb],
                 [     -s-cxb, 0.,      c-czb]],
                [[       -cxc, 0.,       -czc],
                 [ v[0]*c-cxc, 0.,-v[0]*s-czc],
                 [-v[1]*c-cxc, 0.,-v[1]*s-czc]]])
    dx[2] *= 0.5 # chain rule
    return tensordot(dx, de, axes=[(1,2), (0,1)])

def main(argv):
    S = 100
    x = zeros((S, 3, 3))
    F = zeros((S, 3, 3))
    for i in xrange(S):
	v = rand.random(3)*array([0.1, 0.1, 10*pi/180.0]) + array([r0,r0,a0])
	x[i] = mk_x(v)
	F[i] = -dE(v)
    write_array("test.x", x)
    write_array("test.f", F)
    return 0

if __name__ == "__main__":
    main(sys.argv)
