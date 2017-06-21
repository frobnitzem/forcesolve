# Scan prm file for hand-coded information, including:
#
# 1. CHARMM's decision on multiplicity for each torsion
# 2. (future) list of atom types requiring impropers
# 3. (one possible future) nonbonded cutoff scheme in use
""" # example:
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
"""
# shift between CTONnb and CTOFnb (10 to 12)
# setting NBXMod 5 means skipping 1-2, 1-3 interactions and treating 1-4 interactions specially
# e14fac is a 1-4 scaling for both ES and LJ (can usually be assumed 1.0)

class PRM:
    def __init__(self, atoms={}, bonds={},
                       angles={}, dihedrals={},
                       impropers={}, nonbonded={}):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.impropers = impropers
        self.nonbonded = nonbonded

    def write(self, name):
	f = open(name, 'w')
	f.write("\nATOMS\n")
	writeem(write_atom, f, self.atoms)
	f.write("\nBONDS\n")
	writeem(write_bond, f, self.bonds)
	f.write("\nANGLES\n")
	writeem(write_angle, f, self.angles)
	f.write("\nDIHEDRALS\n")
	writeem(write_dihedral, f, self.dihedrals)
	f.write("\nIMPROPERS\n")
	writeem(write_improper, f, self.impropers)
        f.write("\n" \
	 +"NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -\n"\
         +"cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n")
	f.write("\nNBFIX\n")
	writeem(write_nb, f, self.nonbonded)
	f.write(""" 
HBOND CUTHB 0.5  ! If you want to do hbond analysis (only), then use
                 ! READ PARAM APPEND CARD
                 ! to append hbond parameters from the file: par_hbond.inp

END
""")

def write_atom(k,v): # (n, mass)
    return "MASS %5d %-6s %9.5f"%(v[0], k, v[1])

def write_bond(k, v): # K, r0
    return "%-6s %-6s  %6.2f   %8.4f"%(k+v)

def write_angle(k, v): # K, theta0, K_UB, r0_ub (deg)
    if len(v) == 2:
        return "%-6s %-6s %-6s  %6.2f   %7.2f"%(k+v)
    return "%-6s %-6s %-6s  %6.2f   %7.2f %6.2f %9.5f"%(k+v)

def write_dihedral(k, vl): # [(n, n, delta (deg))]
    return "\n".join("%-6s %-6s %-6s %-6s %10.5f %d   %6.2f"%(k + (v[1],v[0],v[2]))
			for v in vl)

def write_improper(k, v): # K
    return "%-6s %-6s %-6s %-6s  %10.4f  0     0.00"%(k + (v,))

def write_nb(k, v): # -eps R0
    if isinstance(k, str):
	s = "%-6s  0.0   %10.4f     %0.4f"%(k, -v[0],v[1])
    s = "%-6s %-6s  %10.4f     %0.4f"%(k[0], k[1], -v[0],v[1])
    if len(v) >= 4: # Write special 1,4 2nd.
	s += "  %10.4f     %0.4f"%(-v[2],v[3])
    return s

def writeem(f1, f, kv):
    for key in sorted(kv.keys()):
	f.write(f1(key, kv[key]) + "\n")

# [[String]] -> {[(String,String,String,String)] : (n=Int,k=Float,delta=Float) }
def read_dihedrals(words):
    p = {}
    for line in words:
        if len(line) < 7:
            continue
	# textual ordering
        if line[1] > line[2] \
	   or line[1] == line[2] and line[0] > line[3]:
            line[0], line[3] = line[3], line[0]
            line[1], line[2] = line[2], line[1]
        if not p.has_key(tuple(line[:4])):
	    p[tuple(line[:4])] = []
        p[tuple(line[:4])].append((int(line[5]), float(line[4]), float(line[6])))
    return p

# [[String]] -> {(String,String,String,String) : Float }
# these are harmonic constraint energies for a torsion
def read_impropers(words):
    p = {}
    for line in words:
        if len(line) < 7:
            continue
        #line[1:4] = sorted(line[1:4]) # textual sort on 'other 3'
        line[1:3] = sorted(line[1:3]) # central atoms are exchangeable
        p[tuple(line[:4])] = float(line[4])
    return p

# TODO: do something useful with first 2 lines of nonbonded sec.
def read_nonbonded(words):
    # skip first 2 lines
    p = {}
    for line in words[2:]:
        if len(line) < 4:
            continue
	if len(line) >= 7:
	    p[line[0]] = (-float(line[2]), float(line[3]),
			  -float(line[5]), float(line[6]))
	else:
	    p[line[0]] = (-float(line[2]), float(line[3]),
			  -float(line[2]), float(line[3]))
    return p

# String -> IO(PRM)
def read_prm(name):
    # what / how to parse
    parse = { 'dihedrals': read_dihedrals,
              'impropers': read_impropers,
              'nonbonded': read_nonbonded,
            }
    # parse output
    out = {}

    sec = None
    words = [] # current section parse

    def process(out, sec, words):
        if sec != None: # process last good parse
            out[sec] = parse[sec](words)
        return None, []

    with open(name, 'r') as f:
        for line in f.xreadlines():
            line = line.split('!', 1)[0] # remove comments
            tok = line.split()
            if len(tok) < 1:
                continue
            if len(tok) == 1: # new sections are id-ed by a single token
                sec, words = process(out, sec, words)
                if parse.has_key(tok[0].lower()): # begin new parse
                    sec = tok[0].lower()
                continue
            if tok[0] == "NONBONDED": # Special case, keep first line
                sec, words = process(out, sec, words)
                if parse.has_key(tok[0].lower()): # begin new parse
                    sec = tok[0].lower()
                else:
		    continue
            if sec != None: # store all lines in a parsed section
                words.append(tok)
    process(out, sec, words)

    if not out.has_key('dihedrals'):
        print "Dihedrals not found in prm!"

    return PRM(**out)

