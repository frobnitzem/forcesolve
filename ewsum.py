from ucgrad import *
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

# Ewald sum for 2-particle interactions.
# here, rij = ri - rj, and f = -dE(rij)/dri
def esum(rij, nbr, L, iL, Eta, M2):
    Eta2 = Eta*Eta
    fac = Eta*2/pi**0.5
    iV = iL[0,0]*iL[1,1]*iL[2,2]
    e = zeros(len(ri))
    f = zeros(ri.shape)

    rij -= L[2]*floor(rij[...,2]/L[2,2]+0.5)[...,newaxis]
    rij -= L[1]*floor(rij[...,1]/L[1,1]+0.5)[...,newaxis]
    rij -= L[0]*floor(rij[...,0]/L[0,0]+0.5)[...,newaxis]
    sij = dot(rij, iL)

    g = lat_pts(zeros(3), identity(3), M2)
    g.next() # manual origin contribution
    for m, _, m2 in g:
        s = sin(2*pi*sum(sij*m, -1))
        c = sqrt(1.0 - s*s)
        phi = exp(-m2*pi*pi/Eta2)/m2
        e += c*phi
        f += s*phi*m

    e *= iV/pi
    f = 2.*iV * dot(f, transpose(iL))

    g = lat_pts(zeros(3), L, M2*iV**-0.66666)
    _, dr, d2 = g.next() # nearest image depends on nbr flag
    d2 = sum(rij*rij, -1)
    dist = sqrt(d2)
    if nbr:
        e += -erf(Eta*dist)/dist
        f += ((fac*exp(-Eta2*d2) -  erf(Eta*dist)/dist)/d2)[...,newaxis] * rij
    else:
        e += erfc(Eta*dist)/dist
        f += ((fac*exp(-Eta2*d2) + erfc(Eta*dist)/dist)/d2)[...,newaxis] * rij
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
    iL = la.inv(L)

    d = arange(200)*0.1-10.05
    dx = array([0., 0., 1.])*d[:,newaxis]

    Eta, M2 = Eta_M2(prod(diag(L)))
    e, f = esum(dx, False, L, iL, Eta, M2)

    M = zeros((200,5))
    M[:,0] = d
    M[:,1] = e
    M[:,2:] = f
    write_matrix("ew_pot.dat", M)

# Calculate ES residual and Jacobian
cfac = 332.0716 # kcal/mol * Ang / esu^2
# E = cfac * q_i q_j / r_{ij}
# q_atoms = dot(MQ, q), MQ.shape = (N, Nchg)
def ES_frc(mask, pairs, x, MQ, q, L):
    iL = la.inv(L)
    Eta, M2 = Eta_M2(prod(diag(L)))
    qat = dot(MQ, q)

    f = zeros(x.shape)
    J = zeros(x.shape + (len(q),))
    for i,j in mask:
        e, fq = esum(x[:,i,:]-x[:,j,:], True, L, iL, Eta, M2)
        f[:,i] += fq*qat[i]*qat[j]
        f[:,j] -= fq*qat[i]*qat[j]
        J[:,i] += fq[...,newaxis]*qat[j]*MQ[i]
                + fq[...,newaxis]*qat[i]*MQ[j]
        J[:,j] -= fq[...,newaxis]*qat[j]*MQ[i]
                - fq[...,newaxis]*qat[i]*MQ[j]

    for i,j in pairs:
        e, fq = esum(x[:,i,:]-x[:,j,:], False, L, iL, Eta, M2)
        f[:,i] += fq*qat[i]*qat[j]
        f[:,j] -= fq*qat[i]*qat[j]
        J[:,i] += fq[...,newaxis]*qat[j]*MQ[i]
                + fq[...,newaxis]*qat[i]*MQ[j]
        J[:,j] -= fq[...,newaxis]*qat[j]*MQ[i]
                - fq[...,newaxis]*qat[i]*MQ[j]
    return f*cfac, J*cfac

