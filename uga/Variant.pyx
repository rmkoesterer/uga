## Copyright (c) 2015 Ryan Koesterer GNU General Public License v3
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import Geno
cimport Geno
cimport numpy as np
cimport cython
import math
import logging

module_logger = logging.getLogger("Variant")

cdef class Ref(object):
	def __cinit__(self):
		logger = logging.getLogger("Variant.Ref.__cinit__")
		logger.debug("initialize variant reference")
		self.db = {}

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef load(self, np.ndarray v):
		logger = logging.getLogger("Variant.Ref.load")
		logger.debug("load")
		for row in v:
			self.db[row['___uid___']] = {}
			self.db[row['___uid___']]['chr'] = row['chr']
			self.db[row['___uid___']]['pos'] = row['pos']
			self.db[row['___uid___']]['a1'] = row['a1']
			self.db[row['___uid___']]['a2'] = row['a2']
			self.db[row['___uid___']]['id'] = row['id']
			self.db[row['___uid___']]['id_unique'] = row['id_unique']

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef update(self, np.ndarray v):
		logger = logging.getLogger("Variant.Ref.update")
		logger.debug("update")
		for row in v:
			if not row['___uid___'] in self.db:
				self.db[row['___uid___']] = {}
				self.db[row['___uid___']]['chr'] = row['chr']
				self.db[row['___uid___']]['pos'] = row['pos']
				self.db[row['___uid___']]['a1'] = row['a1']
				self.db[row['___uid___']]['a2'] = row['a2']
				self.db[row['___uid___']]['id'] = row['id']
				self.db[row['___uid___']]['id_unique'] = row['id_unique']

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef complement(allele):
	cdef str x = allele
	if x != "NA":
		letters = list(x)
		comp = []
		for l in letters:
			if l == 'T': 
				c = 'A'
			elif l == 'A':
				c = 'T'
			elif l == 'G': 
				c = 'C'
			elif l == 'C':
				c = 'G'
			elif l == '0':
				c = '0'
			elif l == ',':
				c = ','
			elif l == 'NA':
				c = 'NA'
			elif l == '-':
				c = '-'
			elif l == 'I':
				c = 'D'
			elif l == 'D':
				c = 'I'
			elif l in ['1','2','3','4','5','6','7','8','9','0']:
				c = l
			else:
				c = 'X'
			comp.append(c)
	else:
		comp = ['NA']
	return ''.join(comp)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef get_universal_variant_id(chr_py,pos_py,a1_py,a2_py,delim_py):
	cdef str chr = chr_py
	cdef str pos = pos_py
	cdef str a1 = a1_py[0:1000]
	cdef str a2 = a2_py[0:1000]
	cdef str delim = delim_py
	analogs = [chr + delim + pos + delim + a1 + delim + a2, 
					chr + delim + pos + delim + complement(a1) + delim + complement(a2)]
	if a2 != 'NA':	
		analogs = analogs + [chr + delim + pos + delim + complement(a2) + delim + complement(a1),
							chr + delim + pos + delim + a2 + delim + a1, 
							chr + delim + pos + delim + a1 + delim + 'NA', 
							chr + delim + pos + delim + complement(a1) + delim + 'NA', 
							chr + delim + pos + delim + complement(a2) + delim + 'NA', 
							chr + delim + pos + delim + a2 + delim + 'NA']
	return "_".join(sorted(list(set(analogs))))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double calc_callrate(np.ndarray[np.float64_t, ndim=1] x):
	cdef double xlen = len(x)
	cdef unsigned int ylen = len(x[~np.isnan(x)])
	if ylen > 0:
		return ylen/xlen if xlen > 0 else float('nan')
	else:
		return 0.0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double calc_freq(np.ndarray[np.float64_t, ndim=1] x):
	x = x[~np.isnan(x)]
	cdef unsigned int n = 2 * len(x)
	cdef double count = x.sum()
	return count / n if n > 0 else float('nan')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double calc_mac(np.ndarray[np.float64_t, ndim=1] x):
	cdef double n
	x = x[~np.isnan(x)]
	if len(x) > 0:
		n = min(2.0 * len(x) - np.sum(x),np.sum(x))
	else:
		n = float('nan')
	return n

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double calc_freqx(np.ndarray[np.float64_t, ndim=1] male, np.ndarray[np.float64_t, ndim=1] female):
	male = male[~np.isnan(male)]
	female = female[~np.isnan(female)]
	cdef unsigned int n = len(male) + 2 * len(female)
	cdef double count = (male.sum()/2) + female.sum()
	return count / n if len(male) > 0 and len(female) > 0 else float('nan')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double calc_rsq(np.ndarray[np.float64_t, ndim=1] x):
	x = x[~np.isnan(x)]
	cdef double rsq = float('nan')
	if len(x) > 0 and not set(x).issubset(set([0,1,2])):
		rsq = x.var(ddof=1) / (2 * (x.mean() / 2) * (1 - (x.mean() / 2))) if x.mean() != 0 and x.mean() != 2 else 0.0
		return rsq if rsq >= 0 and rsq <= 1 else 1 / rsq
	else:
		return rsq

def calc_hwe(np.ndarray[np.float64_t, ndim=1] x):
	cdef int i,h
	cdef unsigned int x_hets
	cdef unsigned int x_hom1
	cdef unsigned int x_hom2
	cdef unsigned int x_homc
	cdef unsigned int x_homr
	cdef unsigned int rare_copies
	cdef unsigned int genotypes
	cdef unsigned int rare
	cdef unsigned int common
	x = x[~np.isnan(x)]
	if len(x) > 0 and set(x).issubset(set([0,1,2])):
		x_hets = np.sum(x == 1)
		x_hom1 = np.sum(x == 2)
		x_hom2 = np.sum(x == 0)
		if x_hets >= 0 and x_hom1 >= 0 and x_hom2 >= 0:
			rare   = 2*min(x_hom1,x_hom2)+x_hets
			common = 2*max(x_hom1,x_hom2)+x_hets
			if not rare:
				return 1.0
			hets = rare*common/(rare+common)
			if rare%2 != hets%2:
				hets += 1
			hom_r = (rare-hets)/2
			hom_c = (common-hets)/2
			prx = [0]*(rare/2+1)
			prx[hets/2] = 1.0
			for i,h in enumerate(xrange(hets,1,-2)):
				prx[h/2-1] = prx[h/2]*h*(h-1) / (4*(hom_r+i+1)*(hom_c+i+1))
			for i,h in enumerate(xrange(hets,rare-1,2)):
				prx[h/2+1] = prx[h/2]*4*(hom_r-i)*(hom_c-i) / ((h+1)*(h+2))
			p_x = prx[x_hets/2]
			p_hwe = sum(p for p in prx if p <= p_x)/sum(prx)
			return p_hwe
		else:
			return float('nan')
	else:
		return float('nan')