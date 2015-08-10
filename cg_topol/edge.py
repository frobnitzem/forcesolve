# Edge creation logic
def add_edge(edge, i, j):
    if i < j:
        edge.add((i,j))
    else:
        edge.add((j,i))
	
def add_all_redge(edge, pdb, res_table, oi, ai, oj, aj):
    mino = min(oi,oj)
    maxo = max(oi,oj)
    for i in reversed(range(len(res_table))):
        if res_table[i]+mino < 0 \
            or res_table[i]+maxo >= len(pdb.res):
                continue
        cchain = pdb.res[res_table[i]].chain_id
                    
        ir = pdb.res[res_table[i]+oi]
        jr = pdb.res[res_table[i]+oj]
        if cchain != ir.chain_id or cchain != jr.chain_id:
            continue
        try:
            iat = [a.split()[0] for a in ir.names].index(ai)
        except ValueError:
            print "Warning! From atom %s not present in "\
               "residue %s %d offset %d from %s"%(ai, \
                ir.name, res_table[i]+oi+1, oi, \
                pdb.res[res_table[i]].name)
            del res_table[i]
            continue
        try:
            jat = [a.split()[0] for a in jr.names].index(aj)
        except ValueError:
        # Common for variable-composition residues...
            continue
        add_edge(edge, ir.atom_zero+iat, jr.atom_zero+jat)
    
def append_edge(edge, pdb, line):
    i = int(line[0])-1
    if i < 0 or i > pdb.atoms:
        raise InputError, "Error EDGE line contains out-of"\
                "-range from atom number %d.\n"%(i+1)
    for to in line[1:]:
        j = int(to)
        if j < 0 or j > pdb.atoms:
            raise InputError, "Error EDGE line contains out-of"\
                "-range to atom number %d->%d.\n"%(i+1, j+1)
        add_edge(edge, i, j)

def append_redge(edge, pdb, line):
    rname = line[0]
    anames = parse_aname(line[1:])
    res_table = [ r for r in range(len(pdb.res)) \
                    if pdb.res[r].name == rname ]
    #print anames
    oi = anames[0][0]
    ai = anames[0][1]
    for oj, aj in anames[1:]:
        add_all_redge(edge, pdb, res_table, oi,ai, oj,aj)

# Parses a list of names and optional offsets into a list
# of the form [(off, name), (off, name), ...]
def parse_aname(tok):
    resname = []
    i = 0
    off = 0 # Default.
    while i < len(tok):
        if tok[i][0] in "+-":
            off = int(tok[i]) # Read an optional offset.
            i += 1
        else:
            off = 0
        resname.append((off, tok[i]))
        i += 1

    return resname

# Modular product set of input sets.
# Very useful for enumerating those pesky angles/torsions...
def modprod(*a):
    b = [ [(i,) for i in a[0]] ]
    for d in range(1, len(a)):
        b.append([])
        for j in a[d]:
            b[d] += [i+(j,) for i in b[d-1]]
    c = set(b[-1])
    del b
    return c

