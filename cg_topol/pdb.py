from ucgrad import read_pdb, array
from edge import *

#   PDB : {
#     conn : [Set Int], -- complete connection table,
#     edge : Set (Int,Int), -- unique edges,
#     x : Array Float n 3, -- example / reference config.
#     names : [(name : String, res : String, atype : String)],
#     mass : Array Float n,
#     atoms : Int, -- number of atoms
#   }
class PDB:
    def __init__(self, names, m, x, edge):
        self.atoms = len(names)
        self.names = names
        self.mass = m
        self.x = x
        self.edge = edge
        self.conn = build_conn(self.atoms, edge)

# Inputs:
#   params -- parsed parameter file
# Output:
#   PDB
def build_pdb(params):
    name = params['EXPDB'][0][0]
    pdb = read_pdb(name, 1)[0]

    # Compute the atom index at the start of each residue.
    res_start = reduce(lambda x, y: x+[x[-1]+y], \
                       [ len(res.names) for res in pdb.res ], [0])
    for r, s in zip(pdb.res, res_start): # and store in pdb
        r.atom_zero = s

    edge = build_edges(pdb, params)
    tindex, tmass = build_tindex(params)
    names = build_name_index(pdb, tindex)
    mass = array([tmass[n[2]] for n in names])
    return PDB(names, mass, pdb.x, edge)

############## Parse EDGE, REDGE, RAT, TRA, and TMASS ##############
def build_edges(pdb, params):
    edge = set()
    if params.has_key("EDGE"):
        for line in params["EDGE"]:
            append_edge(edge, pdb, line)
    if params.has_key("REDGE"):
        for line in params["REDGE"]:
            append_redge(edge, pdb, line)
    return edge

def build_tindex(params):
    tindex = {}
    tmass  = {}
    if params.has_key("RAT"):
        for line in params["RAT"]:
            add_rat(tindex, (line[0],line[1]), line[2])
    if params.has_key("TRA"):
        for line in params["TRA"]:
            t = line[0]
            for i in range(1, len(line), 2):
                add_rat(tindex, (line[i],line[i+1]), t)
    if params.has_key("TMASS"):
        for line in params["TMASS"]:
            tmass[line[0]] = 1836.1527*float(line[1])
    print tmass
    return tindex, tmass

def add_rat(tindex, ra, t):
    if tindex.has_key(ra):
        if tindex[ra] != t:
            raise InputError, "Error! ra pair " + str(ra) + " was " +\
              str(self.tindex[ra])+", but is redefined as type "+t+"."
    else:
        tindex[ra] = t

######################################################

def build_name_index(pdb, tindex):
	names = [] # i -> (rname,aname,type)
	r = 0
	for i in range(pdb.atoms):
		k = i - pdb.res[r].atom_zero
		while k >= len(pdb.res[r].names):
			r += 1
			k = i - pdb.res[r].atom_zero
		rname = (pdb.res[r].name, \
				pdb.res[r].names[k].split()[0])
		try: # Checks atom typing as well.
			t = tindex[rname]
		except KeyError:
			t = "%s:%s"%(rname[0],rname[1])
			print "Warning! Auto-typing %s."%t
			tindex[rname] = t
		names.append((rname[0],rname[1],t))
	return names

# Builds the internal connection table from a set of edges (self.edge).
def build_conn(atoms, edge):
    conn = [] # Connection table.
    for i in range(atoms):
        conn.append(set())
    for i,j in edge:
        conn[i].add(j)
        conn[j].add(i)
    return conn

