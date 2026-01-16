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

def read_array(name, shape):
        list = read_list(name)
        pi = 1
        for n in shape:
                pi *= n
        if(len(list) != pi):
                if(len(list) > pi):
                 print("Warning! Data file \"%s\" has %d values and " \
                  "is incompatible with shape: "%(name, len(list)) + str(shape))
                 print("Using last %d values to construct array."%(pi))
                 list = list[-pi:] # take end of file
                else:
                 raise RuntimeError("Data file \"%s\" has %d values and " \
                  "is incompatible with shape: "%(name, len(list)) + str(shape))
        
        return reshape(list, shape)

# matrix data, find
def read_matrix(name, sep=" ,\t\n"):
        import re
        digits = "01234567890.-+"
        
        isfile = type(name) == file
        if isfile:
                ifile = name
        else:
                ifile = open(name)
        list = []
        
        sep = re.compile("["+re.escape(sep)+"]*")
        
        cols = -1
        i = 0
        for line in ifile.xreadlines():
                j = 0
                for tok in sep.split(line):
                        if not tok: # Ignore blank generated from line's \n.
                                continue
                        if tok[0] not in digits:
                                break
                        list.append(float(tok))
                        j += 1
                if(j > 0):
                    if(cols == -1):
                        cols = j
                    elif(j != cols):
                      print("Warning! line %d contains %d columns, "\
                          "terminating read and ignoring last line."%(i+1, j))
                      #raise ValueError, "Error, line %d contains wrong " \
                        #        "number of columns"%(i+1)
                      list = list[:-j]
                      break
                    i += 1
        if(cols == -1):
                raise ValueError("No data!")
        
        if not isfile:
                ifile.close()
        return reshape(array(list), (i,cols))

def read_list(name):
        digits = "01234567890.-+"
        
        isfile = type(name) == file
        if isfile:
                ifile = name
        else:
                ifile = open(name)
        list = []
        
        for line in ifile.xreadlines():
                for tok in line.split():
                        if(tok[0] not in digits):        # stop line parsing
                                break                        # at first non-numeric
                        list.append(float(tok))                # value
        
        if not isfile:
                ifile.close()
        return array(list)

def write_graph(name, list, *dims):
        if len(list.shape) == 1:
                write_plot(name, list, *dims)
        elif len(list.shape) == 2:
                write_surface(name, list, *dims)
        else:
            raise RuntimeError("Error! %dD graphs not supported."%(len(dims)))

def write_plot(name, list, dim=(0., 1.)):
        out = open(name, 'w')
        
        if(len(list.shape) != 1):
                raise ValueError("Error: List should be 1D!")
        
        for i in range(len(list)):
                out.write("%10f %e\n"%(dim[0]+(i+0.5)*dim[1], list[i]))
        out.write('\n')
        
        out.close()

def write_surface(name, mat, idim=(0., 1.), jdim=(0., 1.)):
        out = open(name, 'w')
        
        if(len(mat.shape) != 2):
                raise ValueError("Error: Matrix should be 2D!")
        
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                out.write("%8.3f %8.3f %e\n"%(idim[0]+(i+0.5)*idim[1], \
                        jdim[0]+(j+0.5)*jdim[1], mat[i,j]))
            out.write('\n')
        out.write('\n')
        
        out.close()

def flatten(d, fmt=" %d"):
        d = reshape(d, d.size)
        s = ""
        for i in d:
                s += fmt%i
        return s

def write_array(name, arr, mode='w'):
        if(len(arr.shape) < 3):
                if(len(arr.shape) == 1):
                    write_list(name, arr, mode)
                elif(len(arr.shape) == 2):
                    write_matrix(name, arr, mode)
                else:
                    print("Warning! zero-size array not written to %s"%(name))
                return 0
        
        # loop over highest dimension and recursively call myself
        write_array(name, arr[0], mode) # possibly the first call to me
        ofile = open(name, 'a')
        ofile.write("\n")  # append a newline to give dimension to the data
        ofile.close()
        for i in range(1, arr.shape[0]):
                write_array(name, arr[i], 'a')
                ofile = open(name, 'a')
                ofile.write("\n")
                ofile.close()

def write_list(name, list, mode='w'):
        out = open(name, mode)
        
        if(len(list.shape) != 1):
                raise ValueError("Error: List should be 1D!")
        
        for i in list:
                out.write(str(i) + '\n')
        
        out.close()

def write_matrix(name, mat, mode='w'):
        out = open(name, mode)
        
        if(len(mat.shape) != 2):
                raise ValueError("Error: Matrix should be 2D!")
        
        for i in range(mat.shape[0]):
            line = ""
            for j in range(mat.shape[1]):
                line += " " + str(mat[i,j])
            out.write(line + '\n')
            del line
        
        out.close()

# really just here for me to remember HOW-TO
def read_binarray(ifile, shape="linear", type='float'):
        if shape == "linear":
                return fromfile(ifile, type)
        else:
                return reshape(fromfile(ifile, type), shape)

def write_binarray(ofile, x, type='float'):
        x.astype(type).tofile(ofile)
