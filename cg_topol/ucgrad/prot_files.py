# This file is part of ucgrad, Copyright 2008 David M. Rogers.
#
#   ucgrad is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ucgrad is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ucgrad (i.e. frc_solve/COPYING).
#   If not, contact the author(s) immediately \at/ wantye \/ gmail.com or
#   http://forceSolve.sourceforge.net/. And see http://www.gnu.org/licenses/
#   for a copy of the GNU GPL.

from numpy import *
from .array_io import *

mass_lookup = { "N":14.00674,        "C":12.011,        "O":15.9994, \
                "S":32.06,        "H":1.00794,        "P":30.973761, \
                "1H":1.00794,        "2H":1.00794,        "3H":1.00794, \
                "CL":35.453 }

# Read in a pdb file with multiple copies
# and return a pdb object and a multiconformer coordinate array.
def read_pdb(filename, max_structs=0, paranoia=0):
        mypdb = 0
        coords = []
        
        file = open(filename)
        pdb = []
        
        for line in file.xreadlines():
            if(len(line) < 5):
                pass
            elif(line[0:5] == "MODEL"):
                pass
            elif(line[0:6] == "ENDMDL"):
                if(mypdb == 0): # add pdb
                        mypdb = protein(pdb)
                        coords.append(mypdb.x)
                else:
                    if paranoia == 0:
                        coords.append(coords[0].copy())
                        get_coords(coords[-1], pdb, len(coords[0]))
                    else:
                        prot = protein(pdb)
                        if paranoia > 1 and cmp_pdb_seq(prot, mypdb, 1):
                            print("Discrepancy found between model 1 "\
                                        "and %d"%(len(coords)+1))
                            raise RuntimeError("Error! Input PDB contains"\
                                        " models of differing types.")
                        coords[-1] = prot.x.copy()
                        del prot
                
                pdb = []
                if(len(coords) == max_structs):
                        break
            else:
                pdb.append(line)
        
        if(len(pdb) != 0): # add last pdb
                if(mypdb == 0):
                        mypdb = protein(pdb)
                        coords.append(mypdb.x)
                else:
                        prot = protein(pdb)
                        if(cmp_pdb_seq(prot, mypdb)):
                                raise RuntimeError("Error! Input PDB contains"\
                                        " models of differing types.")
                        coords.append(prot.x)
                
        file.close()
        return mypdb, array(coords)

# Runs a specified function for each MODEL,
# interpreted as its own protein object
def protloop(filename, eafunc):
        file = open(filename)
        pdb = []
        inmodel = 0
        n = 0
        
        for line in file.xreadlines():
            if(len(line) < 5):
                pass
            elif(line[0:5] == "MODEL"):        # handle multiple models
                if(inmodel):
                        raise ValueError("Encountered a MODEL " \
                                        "inside of a MODEL")
                inmodel = 1
            elif(line[0:6] == "ENDMDL"):
                if(not inmodel):
                        raise ValueError("Encountered ENDMDL " \
                                        "outside of a MODEL")
                eafunc(protein(pdb))
                n += 1
                if(n % 1000 == 0):
                        print("Added structure %d"%(n))
                inmodel = 0
                del pdb
                pdb = []
            else:
                pdb.append(line)
        
        if(len(pdb) != 0):
                eafunc(protein(pdb))
                n += 1
                print("Added final structure, %d"%(n))
                del pdb
                
        file.close()
        return n

# Runs a specified fuction for each MODEL, where the protein object is created
# only for the first MODEL and all others just read coordinates from the pdb.
def structloop(filename, eafunc):
        file = open(filename)
        pdb = []
        inmodel = 0
        n = 0
        atoms = 0
        prot = protein()
        
        for line in file.xreadlines():
            if(len(line) < 5):
                pass
            elif(line[0:5] == "MODEL"):        # handle multiple models
                if(inmodel):
                        raise ValueError("Encountered a MODEL " \
                                        "inside of a MODEL")
                inmodel = 1
            elif(line[0:6] == "ENDMDL"):
                if(not inmodel):
                        raise ValueError("Encountered ENDMDL " \
                                        "outside of a MODEL")
                if n == 0:
                        prot = protein(pdb)
                        atoms = len(prot.x)
                else:
                        get_coords(prot.x, pdb, atoms)
                
                eafunc(prot)
                n += 1
                if(n % 1000 == 0):
                        print("Added structure %d"%(n))
                inmodel = 0
                del pdb
                pdb = []
            else:
                pdb.append(line)
        
        if(len(pdb) != 0):
                eafunc(protein(pdb))
                n += 1
                print("Added final structure, %d"%(n))
                del pdb
                
        file.close()
        return n

# Runs a specified fuction for each MODEL.
# This one ignores pdb information completely and just reads coordinates.
def coordloop(filename, eafunc):
        file = open(filename)
        pdb = []
        inmodel = 0
        n = 0
        
        for line in file.xreadlines():
            if(len(line) < 3):
                pass
            elif(line[0:5] == "MODEL"):        # handle multiple models
                if(inmodel):
                        raise ValueError("Encountered a MODEL " \
                                        "inside of a MODEL")
                inmodel = 1
            #elif(line[0:6] == "ENDMDL"):
            elif(line[0:3] == "END"):
                #if(not inmodel):
                #        raise ValueError("Encountered ENDMDL " \
                #                        "outside of a MODEL")
                
                eafunc(scam_coords(pdb))
                n += 1
                if(n % 1000 == 0):
                        print("Added structure %d"%(n))
                inmodel = 0
                del pdb
                pdb = []
            else:
                pdb.append(line)
        
        if(len(pdb) != 0):
                eafunc(scam_coords(pdb))
                n += 1
                print("Added final structure, %d"%(n))
                del pdb
                
        file.close()
        return n

# Checks for same-ness in protein obj.s a and b
# Ret True if different and False if same
def cmp_pdb_seq(a, b, verb=0):
        if a.seq != b.seq:
                if verb:
                        print("Sequences differ.")
                return True
        if a.start_res != b.start_res:
                if verb:
                        print("Starting residue numbers differ.")
                return True
        if a.start_atom != b.start_atom:
                if verb:
                        print("Starting atom numbers differ.")
                return True
        
        # shouldn't be necessary, but do it anyway
        if len(a.res) != len(b.res):
                return True
        # compare atom names and ordering within residues
        for i in range(len(a.res)):
                ra = a.res[i]
                rb = a.res[i]
                if ra.start_atom != rb.start_atom:
                        return True
                if ra.name != rb.name:
                        return True
                if ra.num != rb.num:
                        return True
                if ra.chain_id != rb.chain_id:
                        return True
                if ra.names != rb.names: # the big one
                        if verb:
                                print("Atoms in residue %s %d differ."%( \
                                                ra.name, a.start_res+i))
                        return True
        
        return False

class residue:
        def __init__(self, lines, name, num, chain_id=' '):
                self.name = name
                self.num = num
                self.chain_id = chain_id
                self.coords = []
                self.names = []
                self.mass = []
                self.ob = []
                if len(lines) > 0:
                    try:
                        self.start_atom = int(lines[0][6:11], 10)
                    except ValueError:
                        try: # Hey.
                            self.start_atom = int(lines[0][6:11], 16)
                        except ValueError:
                            self.start_atom = 1
                else:
                        self.start_atom = 1
                
                for line in lines:
                        self.names.append( line[11:30] )
                        if(len(line) >= 78 and \
                                line[76:78].strip() in mass_lookup.keys()):
                            self.mass.append(mass_lookup[line[76:78].strip()])
                        elif(line[12:14].strip() in mass_lookup.keys()):
                            self.mass.append(mass_lookup[line[12:14].strip()])
                        else:
                            self.mass.append(0.0)
                            print("Warning! PDB line conatains unknown atom:")
                            print("\t" + line)
                        if(len(line) >= 60):
                                self.ob.append([float(line[54:60])])
                        else:
                                self.ob.append([1.0])
                        if(len(line) >= 66):
                                self.ob[-1].append(float(line[60:66]))
                        else:
                                self.ob[-1].append(1.0)
                        
                                # add coordinates
                        pt = array([ float(line[30:38]), float(line[38:46]), \
                                float(line[46:54]) ])
                        
                        self.coords.append(pt)
                self.coords = array(self.coords, float)
                self.mass = array(self.mass, float)
                self.ob = array(self.ob, float)
        
        def __doc__(self):
                return "Residue object contains the following variables:\n" \
                        "  name:       residue name\n"\
                        "  num:        residue number\n"\
                        "  chain_id:   chain identifier character\n"\
                        "  start_atom: starting atom from pdb file\n"\
                        "  names:      list of names, formatted as in: "\
                                "'  N   ARG A   1    '\n"\
                        "  coords:     last read coordinate array for atoms\n"\
                        "  mass:       inferred masses from atom names\n"\
                        "  bfac:       beta factors (1.0 if not present)\n"\
                        "and the following methods:\n"\
                        "  append(name, coord=[0,0,0], mass=1) appends atom\n"\
                        "  get_atom_index(name) returns the index of a given "\
                                "atom name\n"\
                        "  get_atom(name) returns coordinate of a given atom\n"
        
        def __str__(self):
                str = "Residue %s %d on chain %c with %d atoms, "\
                        "numbered starting from %d."%(self.name, self.num, \
                        self.chain_id, len(self.names), self.start_atom)
                return str
        
        def num(self):
                return self.num
        def name(self):
                return self.name
        def chain_id(self):
                return self.chain_id
        def names(self):
                return self.names
        
        # Append an atom
        def append(self, name, coord=[0.,0.,0.], mass=1.0, ob=[1.,1.]):
                self.names.append(" %4s %3s %c%4d    "%(name, self.name,\
                                self.chain_id, self.num))
                self.coords = resize(self.coords, (len(self.names), 3))
                self.coords[-1] = array(coord, float)
                self.mass = resize(self.mass, (len(self.names),))
                self.mass[-1] = mass
                self.ob = resize(self.ob, (len(self.names), 2))
                self.ob[-1] = array(ob)
                #print "Appended atom %s, coord, coords = "%(name) + str(coord)
                #print self.coords
        
        def __del__(self):
                del self.coords
                del self.names
                del self.name
                del self.num
                del self.chain_id
        
        def get_pt(self, n):
                return self.coords[n]
        
        # these routines match the atom name, optimally 4 characters
        def get_atom_index(self, name):
                ncmp = len(name)+1
                for i in range(len(self.names)):
                        if(len(self.names[i][1:]) < ncmp):
                                print("Warning! atom name '" + \
                                    self.names[i][1:] + "' shorter than query!")
                                continue
                        if(self.names[i][1:ncmp] == name):
                                return i+self.start_atom
                return -1
        
        def get_atom(self, name):
                ncmp = len(name)+1
                for i in range(len(self.names)):
                        if(len(self.names[i][1:]) < ncmp):
                                print("Warning! atom name '" + \
                                    self.names[i][1:] + "' shorter than query!")
                                continue
                        if(self.names[i][1:ncmp] == name):
                                return self.coords[i]
                print("Warning! Atom \"%s\" not found in residue %s %d %c!\n"%(\
                        name, self.name, self.num, self.chain_id))
                return zeros(3, float)

# Simple protein structure reading class, ignores MODEL/ENDMDL tags
#   if your pdb contains those, use structloop() or read_pdb() instead.
class protein:
    def __init__(self, lines=[]):
        # initial conditions for reading a frame
        self.seq = []
        self.res = []
        self.start_res = 0
        self.start_atom = 1
        self.comments = []
        
        comment = 1
        cresn = -1
        cchain = -1
        res_lines = []
        for line in lines:
            if(not line):
                continue
            
            if(comment):
                if(line[0:4] == "ATOM" or line[0:6] == "HETATM"):
                        comment = 0
                else:
                        self.comments.append(line)
                        continue
            
                # First, check record
            type = line.split()[0]
            if(not type):
                continue
            if(type not in ["ATOM", "HETATM"] ):
                if(type not in [ "TER" ] ):
                        break
                continue
            if(len(line) < 54):
                print("Error importing %s record with " \
                        "length %d!"%(type, len(line)))
                continue
                
                # Record is genuine, import it
            resn = int(line[22:26])
            chain = line[21:22]
            if(resn == cresn and chain == cchain):
                res_lines.append(line)
            else:
                if(cresn == -1):        # first residue
                        # TODO change to make chain-specific
                        self.start_res = resn
                        self.start_atom = int(line[6:11])
                        
                        name = line[16:20].strip()
                        res_lines.append(line)
                        cresn = resn
                        cchain = chain
                        continue
                self.res.append(residue(res_lines, name, cresn, cchain))
                self.seq.append(name)
                
                cresn = resn
                cchain = chain
                name = line[16:20].strip()
                
                del res_lines
                res_lines = [ line ]
        
        if(len(res_lines) > 0):
                self.seq.append(name)
                self.res.append(residue(res_lines, name, cresn, cchain))
                del res_lines
        
        self.atoms = 0
        newat = 0
        self.x = array([], float)
        self.mass = array([], float)
        self.ob = array([], float)
        for res in self.res:
            newat += len(res.coords)
            self.x = resize(self.x, (newat, 3))
            self.x[self.atoms:] = res.coords
            self.mass = resize(self.mass, (newat,))
            self.mass[self.atoms:] = res.mass
            self.ob = resize(self.ob, (newat, 2))
            self.ob[self.atoms:] = res.ob
            self.atoms = newat
        i = 0
        j = 0
        for res in self.res:
            j = len(res.coords)
            res.coords = self.x[i:j]
            res.mass = self.mass[i:j]
            res.ob = self.ob[i:j]
            i = j
        
    # Appends a residue to the protein
    def append(self, res, num=-1, chain=-1):
        if(num != -1):
                res.num = num
        if(chain != -1):
                res.chain_id = chain
        if(len(self.res) == 0):
                self.start_atom = res.start_atom
                self.start_res = res.num
        
        for r in self.res: # Avoid massive copying and swapping.
            r.coords = None
            r.mass = None
            r.ob = None
        
        self.res.append(res)
        self.seq.append(res.name)
        
        # Add coordinates
        self.atoms = self.x.shape[0] + res.coords.shape[0]
        self.x = resize(self.x, (self.atoms, 3))
        self.x[-res.coords.shape[0]:] = res.coords
        self.mass = resize(self.mass, (self.atoms,))
        self.mass[-res.coords.shape[0]:] = res.mass
        self.ob = resize(self.ob, (self.atoms,2))
        self.ob[-res.coords.shape[0]:] = res.ob
        
        # Renumber all residue-based information.
        i = 0
        j = 0
        for res in self.res:
            j = len(res.names)
            res.coords = self.x[i:j]
            res.mass = self.mass[i:j]
            res.ob = self.ob[i:j]
            i = j

    # consecutively renumbers residues
    def renumber_residues(self, n=1):
        self.start_res = n
        for r in self.res:
                r.num = n
                n += 1

    def x(self):
        return self.x
    def __del__(self):
        del self.seq
        del self.res
        del self.start_res
        
    def res(self):
        return self.res
 
    def write(self, name, mode='w', model=1):
        if type(name) == type("string"):
                file = open(name, mode)
        else:
                file = name # Assume we can write to it, then.
        
        if(mode == 'w'):
            for line in self.comments:
                file.write(line)
        
        if(max(fabs(self.ob[:,0]) < 1000.0)):
                fmt_ob = "%6.2f"
        else:
                fmt_ob = "%6d"
        if(max(fabs(self.ob[:,1]) < 1000.0)):
                fmt_ob += "%6.2f\n"
        else:
                fmt_ob += "%6d\n"
            
        if(model):
                file.write("MODEL %5d\n"%(model))
        
        i = 0
        for res in self.res:
            for j,n in enumerate(res.names):
                if i > 99999:
                        fmt = "ATOM  %5x%s%8.3f%8.3f%8.3f"+fmt_ob
                else:
                        fmt = "ATOM  %5d%s%8.3f%8.3f%8.3f"+fmt_ob
                file.write(fmt%( \
                  self.start_atom+i, n, self.x[i,0], self.x[i,1], self.x[i,2],
                        self.ob[i,0], self.ob[i,1] ))
                i += 1
        if(model):
                file.write("TER\nENDMDL\n")
        if file != name:
                file.close()

# Read the coordinates from a pdb, ignoring the atom name information.
def scam_coords(lines):
        # Grep out atom info.
        lines = [ line for line in lines if line[0:6] in ["ATOM  ", "HETATM"] ]
        x = zeros((len(lines), 3), float)
        
        for n, line in enumerate(lines):
                x[n,0] = float(line[30:38])
                x[n,1] = float(line[38:46])
                x[n,2] = float(line[46:54])
        
        return x

def get_coords(x, lines, atoms):
        # Grep out atom info.
        lines = [ line for line in lines if line[0:6] in ["ATOM  ", "HETATM"] ]
        
        if len(lines) != atoms:
                raise ValueError("Error processing pdb: "\
                    "number of atoms is %d instead of %d."%(len(lines),atoms))
        
        for n, line in enumerate(lines):
                x[n,0] = float(line[30:38])
                x[n,1] = float(line[38:46])
                x[n,2] = float(line[46:54])
        
        return

# Reads Amber restrt files, returning x,v,L.
def read_restrt(name):
        f = open(name)
        head = f.readline()
        head += f.readline()
        crd = read_list(f)
        L = crd[-6:]
        crd = reshape(crd[:-6], (-1,3))
        
        N = len(crd)/2
        
        f.close()
        return crd[:N], crd[N:], L
