#!/usr/bin/env python
# Parse a directory containing ".cor" and ".de" files
# into ".npy" binary matrix data used by "match.py".
# This negates "de" to store forces.
# See line 35 for storage format (xf).

import sys, os
from ucgrad import *
from psf import *

def main(argv):
    assert len(argv) == 4, "Usage: %s <psf> <file dir> <out.npy>"%argv[0]
    psf = read_psf(argv[1])
    xf = read_crds(argv[2], psf.N)
    save(argv[3], xf)

def read_crds(path, N):
    files = os.listdir(path)
    x = {}
    de = {}
    for name in files:
        if name[-4:] == ".cor" or name[-3:] == ".de":
            if name[-4:] == ".cor":
                off = 4
                u = x
            else:
                off = 1
                u = de
            i = int(name.split('.')[-2])
            u[i] = read_vecs(os.path.join(path,name), N, off)

    n = list(set(x.keys()) & set(de.keys())) # set intersection
    n.sort() # sorted frame numbers
    S = len(n)
    xf = zeros((S, 2, N, 3))
    for i,j in enumerate(n):
        xf[i,0] = x[j]
        xf[i,1] = -de[j]
    return xf

# col indicates the starting (whitespace-delim.) column
# where x,y,z coordinates start
def read_vecs(name, N, col):
    with open(name, 'r') as f:
        lines = []
        for line in f.readlines():
            t = line.split()
            if len(t) >= col+3:
                lines.append(t[col:col+3])
    x = [map(float, u) for u in lines[-N:]]
    return array(x)

if __name__=="__main__":
    main(sys.argv)
