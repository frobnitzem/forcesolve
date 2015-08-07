class poly_func:
        def __init__(self, n, id=""):
            self.n = n+1
            self.id = id
            self.rng = [-inf, inf]
            self.verb = False
            self.c = zeros(n+1)

        def commit(self):
            pass

        def ay(self, x, nd=None):
            return self.y(x, nd)

	# Calculates deriviatives 0,...,n of the function.
        def y(self, x, nd=None):
            return eval_poly(x, self.spl.c, nd)

	# Used to construct vectors which multiply parameters.
	# If x is a N-dim vector, the return value is an nd+1xNxP matrix
        def spline(self, x, nd=None):
            if nd != None:
                Mp = x[:,newaxis,newaxis]**(arange(self.n)[newaxis,:,newaxis] - arange(nd+1)[newaxis,newaxis,:])
                # d^n/dx^n [x^a] = x^{a-n} prod_{i=0}^{n-1} a - i
                #                = x^{a-n} prod_{j=a+1-n}^a j
                for i in range(nd):
                    Mp[:,:,i+1:] *= arange(self.n)-i
            else:
                Mp = x[:,newaxis]**(arange(self.n)[newaxis,:])
            return Mp

        def write_spl(self, name, mode='w'):
            out = open(name, mode)
            if len(self.id < 1):
                id = "file"
            else:
                id = self.id.replace(' ', '-').replace('\t', '_'\
                                ).replace('\n', '')
		hdr = "#POLY %s order %d\n"%(id,self.n-1) # File header.
                out.write(hdr)
                for i,c in enumerate(self.c):
                    out.write("%d %12e\n"%(i, c))
                out.close()

