from numpy import *

class poly_func:
        def __init__(self, n, peridic=False, rng=[-100,100], id=""):
            self.n = n+1
            self.id = id
            self.rng = [-100, 100]
            self.verb = False
            self.c = zeros(n+1)

        def commit(self):
            pass

	#TODO: calc. coeff. matrix for integral y^{(n)}(x)^2 * x**m
	def quadratic_integral(self, n, m):
	    return zeros((self.n, self.n))

        def ay(self, x, nd=None):
            return self.y(x, nd)

	# Calculates deriviatives 0,...,n of the function.
        def y(self, x, nd=None):
            return eval_poly(x, self.spl.c, nd)

	# Used to construct vectors which multiply parameters.
	# If x is a N-dim vector, the return value is an (nd+1)xNxP matrix
        def spline(self, x, nd=0):
            if nd == None:
                return x[..., newaxis]**arange(self.n)

            Mp = x[...,newaxis,newaxis]**(arange(self.n)[:,newaxis] - arange(nd+1)[newaxis,:])
            # d^n/dx^n [x^a] = x^{a-n} prod_{i=0}^{n-1} a - i
            #                = x^{a-n} prod_{j=a+1-n}^a j
            for i in range(nd):
                Mp[..., :, i+1:] *= (arange(self.n)-i)[:,newaxis]
            return transpose(Mp, [len(Mp.shape)-1,] + range(len(Mp.shape)-1))

        def write_spl(self, name, mode='w'):
            out = open(name, mode)
            if len(self.id) < 1:
                id = "file"
            else:
                id = self.id.replace(' ', '-').replace('\t', '_'\
                                ).replace('\n', '')
		hdr = "#POLY %s order %d\n"%(id,self.n-1) # File header.
                out.write(hdr)
                for i,c in enumerate(self.c):
                    out.write("%d %12e\n"%(i, c))
                out.close()

