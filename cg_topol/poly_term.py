from numpy import *

class PolyTerm:
    def __init__(self, name, n):
        self.hyp_params = 0
        self.params = n+1
        self.name = name
        # first coeff = 0
        self.constraints = [zeros(self.params)]
        self.constraints[0][0] = 1.0
        self.ineqs = []
        self.prior = []
        self.pri_rank = []
        self.ind = []

    # Calculates deriviatives 0,...,n of the function.
    def y(self, c, x, nd=None):
        return eval_poly(x, c, nd)

    # Used to construct vectors which multiply parameters.
    # If x is a N-dim vector, the return value is an (nd+1)xNxP matrix
    def spline(self, x, nd=0):
        if nd == None:
            return x[..., newaxis]**arange(self.params)

        Mp = x[...,newaxis,newaxis]**(arange(self.params)[:,newaxis] - arange(nd+1)[newaxis,:])
        # d^n/dx^n [x^a] = x^{a-n} prod_{i=0}^{n-1} a - i
        #                = x^{a-n} prod_{j=a+1-n}^a j
        for i in range(nd):
            Mp[..., :, i+1:] *= (arange(self.params)-i)[:,newaxis]
        return transpose(Mp, [len(Mp.shape)-1,] + range(len(Mp.shape)-1))

    def write(self, pre, c, mode='w'):
        #name = self.name.replace(' ', '-').replace('\t', '_'\
        #                    ).replace('\n', '')
        name = pre + self.name + ".poly"
        out = open(name, mode)
        hdr = "#POLY %s order %d\n"%(self.name,self.params-1) # File header.
        out.write(hdr)
        for i,c in enumerate(c):
            out.write("%d %12e\n"%(i, c))
        out.close()

def read_poly_term(file):
	digits = "01234567890.-+"
	
	lines = open(file).readlines()
	info = lines[0][1:].split()
	#print info
	if info[0] != "POLY" or len(info) != 4:
		raise InputError, "Error! Input file is not POLY type."

	id = info[1].split("_")
	
	n = int(info[3])
	c = zeros(n+1)
	i = 0
	for line in lines[1:]:
		tok = line.split()
		if len(tok) > 1:
		   if(tok[0][0] in digits and tok[1][0] in digits):
		      if i > n:
			print "Warning! Extra data in POLY file "\
				"ignored!!"
			break
		      c[i] = float(tok[1])
		      i += 1
	if i != n+1:
		print "Warning! Not all spline coefficients read from %s!"%file
	return id, c
