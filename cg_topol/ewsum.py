from numpy import *
from scipy.special import erf, erfc

sign = lambda x: 1-2*(x<0)
def seq(z):
    i0 = int(round(z))
    i = 0
    if i0 < z: # rounded down.
        while True:
            yield i0+i
            i = -i+(i<=0) # 1, -1, 2, -2, ...
    else:
        while True:
            yield i0+i
            i = -i-(i>=0) # -1, 1, -2, 2, ...

# Enumerate all lattice vectors, m^T L, within R2 of x0
# Returns (m, m^T L - x0, |m^T L - x0|^2)
def lat_pts(x0, L, R2):
    pt = zeros(3)
    m = zeros(3)
    for k in seq(x0[2]/L[2,2]):
        m[2] = k
        pt[2] = k*L[2,2] - x0[2]
        Cz = R2 - pt[2]**2
        if Cz < 0.0:
            break
        pt1   = k*L[2,1] - x0[1]
        for j in seq(-pt1/L[1,1]):
            m[1] = j
            pt[1] = j*L[1,1] + pt1
            Cy = Cz - pt[1]**2
            if Cy < 0.0:
                break
            pt2 = k*L[2,0] + j*L[1,0] - x0[0]
            for i in seq(-pt2/L[0,0]):
                m[0] = i
                pt[0] = i*L[0,0] + pt2
                Cx = Cy - pt[0]**2
                if Cx < 0.0:
                    break
                yield m, pt, R2 - Cx

# Newton-Rhapson for monotonic functions only!
def newton(f, df, x0):
    x = x0
    err = f(x)
    while abs(err) > 1e-10:
        x -= err/df(x)
        err = f(x)
    return x

def Eta_M2(V, tol=1e-8):
    C = -log(tol)
    Eta = sqrt(pi*V**-0.33333)
    fac = (pi/Eta)**2
    M2 = newton(lambda x: log(x) + fac*x - C,
                lambda x: 1./x + fac, C/fac)
    return Eta, M2

def rev2D(x):
    return array([[x[2,2], x[2,1], x[2,0]],
                  [x[1,2], x[1,1], x[1,0]],
                  [x[0,2], x[0,1], x[0,0]]])

def rev(x):
    return array([x[2], x[1], x[0]])

# Ewald sum for 2-particle interactions.
# Here, rij = ri - rj, and f = -dE(rij)/dri
# As it's a pairwise energy, the Ewald correction constants are not included.
# count = 1 => full interaction, 0 => only LR interaction (subtract 1/r)
def esum(rij, count, L, iL, Eta, M2, recip_only=False):
    Eta2 = Eta*Eta
    fac = Eta*2/pi**0.5
    iV = iL[0,0]*iL[1,1]*iL[2,2]
    e = zeros(len(rij))
    f = zeros(rij.shape)

    rij -= L[2]*floor(rij[...,2]/L[2,2]+0.5)[...,newaxis]
    rij -= L[1]*floor(rij[...,1]/L[1,1]+0.5)[...,newaxis]
    rij -= L[0]*floor(rij[...,0]/L[0,0]+0.5)[...,newaxis]
    #sij = dot(rij, iL)

    # iL is lower-diagonal
    # and its columns represent reciprocal lattice vectors.
    g = lat_pts(zeros(3), rev2D(iL.transpose()), M2)
    g.next() # manual origin contribution
    for _, m, m2 in g:
        m = rev(m)
        t = 2*pi*sum(rij*m, -1)
        s = sin(t)
        c = cos(t)
        phi = exp(-m2*pi*pi/Eta2)/m2
        e += c*phi
        f += (s*phi)[:,newaxis]*m

    e *= iV/pi
    f *= 2.*iV

    if recip_only:
        return e, f

    g = lat_pts(zeros(3), L, M2*iV**-0.66666)
    _, dr, d2 = g.next() # nearest image subject to mask
    d2 = sum(rij*rij, -1)
    dist = sqrt(d2)
    # count = 1 => full count, 0 => mask out
    h = (count - erf(Eta*dist))/dist
    e += h
    f += ((fac*exp(-Eta2*d2) + h)/d2)[...,newaxis] * rij
    for _, dr, d2 in g:
        q = rij - dr
        d2 = sum(q*q, -1)
        dist = sqrt(d2)
        e += erfc(Eta*dist)/dist
        f += ((fac*exp(-Eta2*d2) + erfc(Eta*dist)/dist)/d2)[...,newaxis] * q

    return e, f

# generate a plot of pot. and derivative vs. separation distance
def test():
    L = array([[10., 0., 0.],
               [0., 10., 0.],
               [0., 0., 10.]])
    L = array([[12., 0., 0.],
               [1., 11., 0.],
               [0., 0., 10.]])
    iL = la.inv(L)

    d = arange(200)*0.1-10.05
    dx = array([0., 0., 1.])*d[:,newaxis]

    Eta, M2 = Eta_M2(prod(diag(L)))
    #print Eta, M2
    Eta = 0.5
    M2 = 0.5**2
    e, f = esum(dx, 1.0, L, iL, Eta, M2, True)

    M = zeros((200,5))
    M[:,0] = d
    M[:,1] = e
    M[:,2:] = f
    write_matrix("ew_pot.dat", M)

# Calculate ES residual and Jacobian
# mask is a dictionary mapping ordered (i,j) pairs to scaling factors
# for the nearest-neighbor 1/r (Coulomb) energy interaction.
# All other (non-self) pairwise interactions are included with a scale of 1.
#
# cfac = 332.0716 # kcal/mol * Ang / esu^2
# E_ij = mask[(i,j)] * cfac * q_i q_j / r_{ij}
# q_atoms = dot(MQ, q), MQ.shape = (N, Nchg)
def ES_seed(x, mask, MQ, L, cfac = 332.0716):
    if L == None: # non-periodic
	return inf_ES_seed(x, mask, MQ, cfac)
    iL = la.inv(L)
    Eta, M2 = Eta_M2(prod(diag(L)))

    P = MQ.shape[-1]
    fs = zeros(x.shape + (P,P))
    for i in range(x.shape[1]-1):
      for j in range(i+1, x.shape[1]):
        try:
            c = mask[(i,j)]
        except KeyError:
            c = 1.0
        e, fq = esum(x[:,i,:]-x[:,j,:], c, L, iL, Eta, M2)
	S = MQ[i,:,newaxis]*MQ[j,newaxis,:] \
	  + MQ[j,:,newaxis]*MQ[i,newaxis,:]
	fs[:,i] += fq[...,newaxis,newaxis]*S
	fs[:,j] -= fq[...,newaxis,newaxis]*S

    return fs*cfac

# Non-periodic version.
def inf_ES_seed(x, mask, MQ, cfac = 332.0716):
    P = MQ.shape[-1]
    fs = zeros(x.shape + (P,P))

    for i in range(x.shape[1]-1):
      for j in range(i+1, x.shape[1]):
        try:
            c = mask[(i,j)]
            if c <= 0.0:
                continue
        except KeyError:
            c = 1.0
	r = x[:,i] - x[:,j]
	r *= (sum(r*r, -1)**-1.5)[...,newaxis]
	S = MQ[i,:,newaxis]*MQ[j,newaxis,:] \
	  + MQ[j,:,newaxis]*MQ[i,newaxis,:]
        S *= c
	fs[:,i] += r[...,newaxis,newaxis]*S
	fs[:,j] -= r[...,newaxis,newaxis]*S

    return fs*cfac

def ES_frc(q, fs, *args):
    return 0.5*dot(dot(fs, q), q)

def dES_frc(q, fs, *args):
    J = dot(fs, q)
    #f = 0.5*dot(J, q)
    #return f, J
    return J

if __name__=="__main__":
    test()
