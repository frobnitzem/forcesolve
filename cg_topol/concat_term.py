from numpy import array, concatenate, zeros

# The central object in the fitting is the FFTerm module,
# which is any class implementing the following API.
#
# These objects are monoids, so can be concatenated to
# make a combined FF using FFconcat.
#
# If no terms are present, the monoid is not technically correct,
# since it returns 0.0 as the force.
#
# Inputs:
#   pdb : PDB,
#   terms : FFTerm, -- parameterized interaction list
# Output:
#   CgTopol : {
#     name   : String, -- term name
#     params : Int, -- number of FF parameters (count_parameters)
#     hyp_params : Int, -- number of hyper-parameters
#     prior  : [(Int, Array Float (d,d))], -- term # and penalty matrix
#     ind    : [Int], -- start indices for ea. term coeff. set
#     pri_rank : [Int] -- len=hyp_params list of alpha for beta-dist'n
#     constraints : [Array Float (params,)],
#     ineqs    : [Array Float (params,)],
#
#     energy : Array Float (params,) -> Array Float (N,3) -> Float,
#     force  : Array Float (params,) -> Array Float (N,3) -> Array Float (N,3),
#     design : Array Float (N, 3) -> order : Int ->
#               {|} Array Float (1,params), order == 0
#               {|} Array Float (N,3,params), order == 1
#               {|} Error
#          -- energy / force design matrix
#   }
# concatenate FFTerms
class FFconcat:
    def __init__(self, terms):
        self.terms = terms
        self.params = sum(t.params for t in terms)
        self.hyp_params = sum(t.hyp_params for t in terms)
        self.pri_rank = reduce(lambda a,b: a+b.pri_rank, self.terms, [])
        
        # Manage starting indices
        self.ind = []
        i = 0
        for t in self.terms:
            self.ind += [k + i for k in t.ind]
            i += t.params

        # Manage term indices
        self.prior = []
        i = 0
        for t in self.terms:
            self.prior += [(k+i, P) for k,P in t.prior]
            i += len(t.ind)

        # Pad constraints with zeros to increase vectors to length = params
        i = 0
        self.constraints = []
        for t in self.terms:
            self.constraints += map(lambda u: pad(u, i, self.params),
                                     t.constraints)
            i += t.params

        # cut and paste for ineqs
        i = 0
        self.ineqs = []
        for t in self.terms:
            self.ineqs += map(lambda u: pad(u, i, self.params), t.ineqs)
            i += t.params

    # Misc. reductions
    def energy(self, c, x):
        E = 0.0
        i = 0
        for t in self.terms:
            E += t.energy(c[i:i+t.params], x)
            i += t.params
        return E
    def force(self, c, x):
        F = 0.0
        i = 0
        for t in self.terms:
            F = F + t.force(c[i:i+t.params], x)
            i += t.params
        return F

    def design(self, x, order=0):
        if order == 0:
            return concatenate([t.design(x, order) \
                                 for t in self.terms], axis=-1)
        else:
            r = [[] for i in range(order+1)]
            for t in self.terms:
                u = t.design(x, order)
                if not isinstance(u[0], list):
                    for i in range(order+1):
                        r[i].append(u[i])
            if len(self.terms) > 0:
                for i in range(order+1):
                    r[i] = concatenate(r[i], axis=-1)
            return r

def pad(x, st, sz):
    u = zeros(sz)
    u[st:st+len(x)] = x
    return u
