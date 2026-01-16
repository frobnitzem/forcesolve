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

import os, shlex
# parser for command line options and parameter files

# parse all options into name:value pairs
def parseflags(blobs, UsageInfo, AcceptedFlags=[], RequiredFlags=[]):
        args = []
        flags = {}
        
        inflag = 0
        for blob in blobs:
            if(inflag == 0):
                if(blob[0] == '-'):
                    if(blob[1:] in AcceptedFlags):
                        inflag = blob[1:]
                    elif(blob[1:].lower() == 'help'):
                        print(UsageInfo)
                        sys.exit(0)
                    else:
                        flags[blob[1:]] = 1
                else:
                        args.append(blob)
            else:
                if(inflag in flags.keys()):
                        flags[inflag].append(blob)
                else:
                        flags[inflag] = [blob]
                inflag = 0
        
        for i in RequiredFlags:
                if(i not in flags.keys()):
                        print("Error! %s not found in command\n" \
                                "options: "%(i) + str(flags))
                        print(UsageInfo)
                        return args, -1
        
        return args, flags

def addition(*a):
        if len(a) < 1:
                return a
        if type(a[0]) == type({}):
                return merge_dicts(addition, *a)
        
        #print("Adding: " + str(a))
        r = a[0]
        for i in a[1:]:
                r += i
        #print "\t-> " + str(r)
        return r

def merge_dicts(merge_fn, *ind):
        #print "Merging: " + str(ind)
        r = ind[0].copy()
        
        for d in ind[1:]:
            for k in d:
                if r.has_key(k):
                    try:
                        r[k] = merge_fn(r[k], d[k])
                    except KeyError:
                        raise KeyError("Error encountered adding key " + str(k))
                else:
                        r[k] = d[k]
        #print "\t-> " + str(r)
        return r

def get_fullpath(name, cwd=None):
        if os.path.isabs(name): # Already done.
                return name
        
        if cwd == None: # Default cwd.
                return os.path.abspath(name)
        else: # Given cwd.
                return os.path.abspath(os.path.join(cwd,name))

# Get all parameters from the file listed below.
# Comment must be a single, delimiting comment string (e.g. '#' or '//').
def scam_params(name, cwd=None, comment='#', included=set()):
        params = {}
        #print "Called with file \"%s\"."%name
        #print included
        
        name = get_fullpath(name, cwd)
        if name in included: # We have already been included!
                print("Warning! >1 reference to %s detected, ignoring."%name)
                return {}
        included.add(name) # Add myself to ''in'' list.
        print("Parsing file \"%s\"."%name)
        
        cwd = os.path.split(name)[0] # Change effective working directory
                                        # to that of the file we are reading.
        file = open(name)
        for line in file.xreadlines():
                tok = shlex.split(line.split(comment)[0])
                if not tok:
                        continue
                
                # file keyword -> change all relative to absolute pathnames.
                if tok[0] == "file":
                        del tok[0]
                        for i in range(1,len(tok)):
                                tok[i] = get_fullpath(tok[i], cwd)
                
                if tok[0] == "include": # Recursively include files.
                        inc_params,included = scam_params(tok[1], \
                                        cwd, comment, included)
                        params = addition(params, inc_params)
                        continue
                
                if tok[0] not in params:        # encountered new parameter
                        params[tok[0]] = []
                
                # One list per line the parameter name is mentioned.
                params[tok[0]].append(tok[1:])
        file.close()
        
        return params, included

def parseparam(name, RequiredParams):
        params, included = scam_params(name, None, '#', set())
        
        for i in RequiredParams:
                if not params.has_key(i):
                  print("Error! %s not found in parameter file(s):"%i)
                  for j in included:
                    print("\t\"" + j + "\"")
                  return -1
        
        return params
