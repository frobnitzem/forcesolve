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

# Creates a counter object to collect histogram data.
# bin_spec is a list of tuples: [spec, spec, ...].
# spec = (start, step, N)
class counter:
	def __init__(self, bin_spec):
		self.dim = len(bin_spec)
		self.start = []
		self.step = []
		self.n = []	# holds shape of bin array
		self.min = []
		self.max = []
		self.out = 0	# number of out-of bounds entries
		self.N = 0
		l = 1
		for i in range(self.dim):
			self.start.append(bin_spec[i][0])
			self.step.append(bin_spec[i][1])
			self.n.append(int(bin_spec[i][2]))
			l *= self.n[-1]
			#print "start = %f, step = %f, n = %d"%(self.start[i],\
			#		self.step[i], self.n[i])
		print "len = %d"%(l)
		if(l > 1.0e+7):
			raise ValueError, "Error! Bin array is more than 20 Mb!"
		self.bin = zeros(self.n)
	
	# sums all bins not in marginal distribution list, 
	# whose numbering starts from 0
	def marginaldist(self, bins):
		# useful?
		bins.sort()
		
		# construct list of dimensions to sum over
		sumit = []
		for i in range(self.dim):
			if(i not in bins):
				sumit.append(i)
		
		print "Constructing marginal distribution of %d " \
			"variables by summing over the other %d.\n"%( \
				len(bins), len(sumit))
		
		bin = self.bin
		for i in range(len(sumit)):
			index = sumit[i]-i
			bin = sum(bin, index)
		
		return bin
		
	def show_status(self):
		print "%d data points taken, %d out of range.\n"%(self.N, \
					self.out)
		
	def append(self, x):
		if(len(x) != self.dim):
			raise RuntimeError, "Bad call to counter.collect()"
		
		# find correct bin
		piece = self.bin
		for d in range(0, self.dim-1):
			j = int((x[d]-self.start[d])/self.step[d])
			if(j < 0 or j >= self.n[d]):
				#print "Point out of range: " + str(x)
				#print "%f < x[%d] < %f"%(self.start[d], d, \
				#	self.start[d]+self.n[d]*self.step[d])
				if self.out%100 == 0:
					print "Collected %d points out " \
						"of range."%(self.out)
				self.out += 1
				return
			piece = piece[j]
		
		d = self.dim-1
		j = int((x[d]-self.start[d])/self.step[d])
		if(j < 0 or j >= self.n[d]):
			#print "Point out of range: " + str(x)
			#print "%f < x[%d] < %f"%(self.start[d], d, \
			#	self.start[d]+self.n[d]*self.step[d])
			self.out += 1
			return
		
		piece[j] += 1
		self.N += 1
	
	def freq(self):
		return self.bin.astype('float')/self.N

# Same as above, except that points can be optionally weighted.
class wt_counter:
	def __init__(self, bin_spec):
		self.dim = len(bin_spec)
		self.start = []
		self.step = []
		self.n = []	# holds shape of bin array
		self.min = []
		self.max = []
		self.out = 0	# number of out-of bounds entries
		self.N = 0.
		l = 1
		for i in range(self.dim):
			self.start.append(bin_spec[i][0])
			self.step.append(bin_spec[i][1])
			self.n.append(int(bin_spec[i][2]))
			l *= self.n[-1]
			#print "start = %f, step = %f, n = %d"%(self.start[i],\
			#		self.step[i], self.n[i])
		print "len = %d"%(l)
		if(l > 1.0e+7):
			raise ValueError, "Error! Bin array is more than 40 Mb!"
		self.bin = zeros(self.n, float)
	
	# sums all bins not in marginal distribution list, 
	# whose numbering starts from 0
	def marginaldist(self, bins):
		if(len(bins) < 0):
			raise RuntimeError, "No Bins selected!"
		# useful?
		bins.sort()
		
		# construct list of dimensions to sum over
		sumit = []
		for i in range(self.dim):
			if(i not in bins):
				sumit.append(i)
		
		print "Constructing marginal distribution of %d " \
			"variables by summing over the other %d.\n"%( \
				len(bins), len(sumit))
		
		bin = self.bin
		for i in range(len(sumit)):
			index = sumit[i]-i
			bin = sum(bin, index)
		
		return bin
		
	def show_status(self):
		print "%d data points taken, %d out of range.\n"%(self.N, \
					self.out)
		
	def append(self, x, wt=1.0):
		if(len(x) != self.dim):
			raise RuntimeError, "Bad call to counter.append()"
		
		#print "Adding " + str(x) + " with wt = %e."%(wt)
		# find correct bin
		piece = self.bin
		for d in range(0, self.dim-1):
			j = int((x[d]-self.start[d])/self.step[d])
			if(j < 0 or j >= self.n[d]):
				#print "Point out of range: " + str(x)
				#print "\t%f < x[%d] < %f"%(self.start[d], d, \
				#	self.start[d]+self.n[d]*self.step[d])
				self.out += 1
				if self.out%100 == 0:
					print "Collected %d points out " \
						"of range."%(self.out)
				return
			piece = piece[j]
		
		d = self.dim-1
		j = int((x[d]-self.start[d])/self.step[d])
		if(j < 0 or j >= self.n[d]):
			#print "Point out of range: " + str(x)
			#print "\t%f < x[%d] < %f"%(self.start[d], d, \
			#	self.start[d]+self.n[d]*self.step[d])
			self.out += 1
			return
		
		piece[j] += wt
		self.N += wt
	
	def freq(self):
		return self.bin/self.N

def fixperiod(v, period):
	if(v < 0.0):
		v += period * (1+int(v/period))
	elif(v >= period):
		v -= period * int(v/period)
	
	return v
	
