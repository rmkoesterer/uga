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

from Model import *
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
from libc.string cimport strtok

@cython.boundscheck(False)
@cython.wraparound(False)
def GetDelimiterCy(delimiter):
	cdef str d = delimiter
	if d == 'tab':
		d = '\t'
	elif d == 'space':
		d = ' '
	else:
		d = ','
	return d

@cython.boundscheck(False)
@cython.wraparound(False)
def GetRowCallsCy(np.ndarray row):
	cdef unsigned int i = row[8].split(':').index('GT')
	cdef unsigned int j = 0
	cdef str het1 = '1/0'
	cdef str het2 = '0/1'
	cdef str hom1 = '0/0'
	cdef str hom2 = '1/1'
	for j in xrange(len(row[9:])):
		row[j+9] = str(strtok(str(row[j+9]),':'))
		row[j+9] = float('NaN') if not str(row[j+9]) in [hom1,het2,hom2,het1] else row[j+9]
		row[j+9] = '2' if str(row[j+9]) == hom1 else row[j+9]
		row[j+9] = '1' if str(row[j+9]) == het1 or str(row[j+9]) == het2 else row[j+9]
		row[j+9] = '0' if str(row[j+9]) == hom2 else row[j+9]
	return row

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcCallrateCy(np.ndarray x):
	cdef double xlen = len(x)
	cdef unsigned int ylen = len(x[~np.isnan(x)])
	if ylen > 0:
		return ylen/xlen
	else:
		return 0.0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcFreqCy(np.ndarray x):
	x = x[~np.isnan(x)]
	cdef unsigned int n = 2 * len(x)
	cdef double count = x.sum()
	return count / n if len(x) > 0 else float('nan')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcFreq23Cy(np.ndarray male, np.ndarray female):
	male = male[~np.isnan(male)]
	female = female[~np.isnan(female)]
	cdef unsigned int n = len(male) + 2 * len(female)
	cdef double count = (male.sum()/2) + female.sum()
	return count / n if len(male) > 0 and len(female) > 0 else float('nan')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcRsqCy(np.ndarray x):
	x = x[~np.isnan(x)]
	cdef double rsq = float('nan')
	if len(x) > 0 and not set(x).issubset(set([0,1,2])):
		rsq = x.var(ddof=1) / (2 * (x.mean() / 2) * (1 - (x.mean() / 2))) if x.mean() != 0 and x.mean() != 2 else 0.0
		return 1 / rsq if rsq > 1 else rsq
	else:
		return rsq

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcHWECy(np.ndarray x):
	cdef unsigned int obs_hets
	cdef unsigned int obs_hom1
	cdef unsigned int obs_hom2
	cdef unsigned int obs_homc
	cdef unsigned int obs_homr
	cdef unsigned int rare_copies
	cdef unsigned int genotypes
	x = x[~np.isnan(x)]
	if len(x) > 0 and set(x).issubset(set([0,1,2])):
		obs_hets = np.sum(x == 1)
		obs_hom1 = np.sum(x == 2)
		obs_hom2 = np.sum(x == 0)
		if not (obs_hets < 0 and obs_hom1 < 0 and obs_hom2 < 0):
			rare   = 2*min(obs_hom1,obs_hom2)+obs_hets
			common = 2*max(obs_hom1,obs_hom2)+obs_hets
			if not rare:
				return 1.0
			hets = rare*common/(rare+common)
			if rare%2 != hets%2:
				hets += 1
			hom_r = (rare-hets)/2
			hom_c = (common-hets)/2
			probs = [0]*(rare/2+1)
			probs[hets/2] = 1.0
			for i,h in enumerate(xrange(hets,1,-2)):
				probs[h/2-1] = probs[h/2]*h*(h-1) / (4*(hom_r+i+1)*(hom_c+i+1))
			for i,h in enumerate(xrange(hets,rare-1,2)):
				probs[h/2+1] = probs[h/2]*4*(hom_r-i)*(hom_c-i) / ((h+1)*(h+2))
			p_obs = probs[obs_hets/2]
			p_hwe = sum(p for p in probs if p <= p_obs)/sum(probs)
			return p_hwe
		else:
			return float('nan')
	else:
		return float('nan')

@cython.boundscheck(False)
@cython.wraparound(False)
def GenerateFilterCodeCy(np.float row_callrate, np.float row_freq, np.float row_freq_unrel, np.float row_rsq, np.float row_rsq_unrel, np.float row_hwe, np.float row_hwe_unrel, np.float miss_thresh=None, np.float maf_thresh=None, np.float maxmaf_thresh=None, np.float rsq_thresh=None, np.float hwe_thresh=None, no_mono=True):
	cdef unsigned int filter = 0
	if (not miss_thresh is None and not np.isnan(row_callrate) and row_callrate < miss_thresh) or (not np.isnan(row_callrate) and row_callrate == 0) or np.isnan(row_callrate):
		filter += 1000
	if not np.isnan(row_freq): 
		if no_mono and (row_freq == 0 or row_freq == 1):
			filter += 100
		else:
			if ((	not maf_thresh is None
				and 
					(		row_freq < maf_thresh
						 or row_freq > 1-maf_thresh
					)
				and 
					(		row_freq_unrel < maf_thresh
						 or row_freq_unrel > 1-maf_thresh
					)
			   ) 
			  or
			   (	not maxmaf_thresh is None
				and 
					(	   
						(		row_freq >= maxmaf_thresh
							and row_freq <= 1-maxmaf_thresh
						)
					)
				and
					(
						(		row_freq_unrel >= maxmaf_thresh
							and row_freq_unrel <= 1-maxmaf_thresh
						)
						or row_freq_unrel == 0
						or row_freq_unrel == 1
					)
			   )):
				filter += 100
	if not rsq_thresh is None and not np.isnan(row_rsq) and (row_rsq < rsq_thresh and row_rsq_unrel < rsq_thresh):
		filter += 10
	if not hwe_thresh is None and not np.isnan(row_hwe) and (row_hwe < hwe_thresh and row_hwe_unrel < hwe_thresh):
		filter += 1
	return filter

@cython.boundscheck(False)
@cython.wraparound(False)
def ComplementCy(allele):
	cdef str x = allele
	if x != "NA":
		letters = list(x)
	else:
		letters = "NA"
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
	return ''.join(comp)

@cython.boundscheck(False)
@cython.wraparound(False)
def ListCompatibleMarkersMetaCy(chr_py,pos_py,a1_py,a2_py,delim_py):
	cdef str chr = chr_py
	cdef str pos = pos_py
	cdef str a1 = a1_py
	cdef str a2 = a2_py
	cdef str delim = delim_py
	cdef str r1, r2, r3, r4
	r1 = chr + delim + pos + delim + a1 + delim + a2
	r2 = chr + delim + pos + delim + ComplementCy(a1) + delim + ComplementCy(a2)
	r3 = chr + delim + pos + delim + ComplementCy(a2) + delim + ComplementCy(a1)
	r4 = chr + delim + pos + delim + a2 + delim + a1
	return "_".join(sorted([r1,r2,r3,r4]))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipEffectCy(refa1_py, refa2_py, a1_py, a2_py, effect_py):
	cdef str refa1 = refa1_py
	cdef str refa2 = refa2_py
	cdef str a1 = a1_py
	cdef str a2 = a2_py
	cdef float effect = effect_py
	if (refa1 + refa2 == ComplementCy(a2) + ComplementCy(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * effect
	else:
		return effect

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipFreqCy(refa1_py, refa2_py, a1_py, a2_py, freq_py):
	cdef str refa1 = refa1_py
	cdef str refa2 = refa2_py
	cdef str a1 = a1_py
	cdef str a2 = a2_py
	cdef float freq = freq_py
	if (refa1 + refa2 == ComplementCy(a2) + ComplementCy(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 - freq
	else:
		return freq

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipORCy(refa1_py, refa2_py, a1_py, a2_py, o_r_py):
	cdef str refa1 = refa1_py
	cdef str refa2 = refa2_py
	cdef str a1 = a1_py
	cdef str a2 = a2_py
	cdef float o_r = o_r_py
	if (refa1 + refa2 == ComplementCy(a2) + ComplementCy(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 / o_r
	else:
		return o_r

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipZCy(refa1_py, refa2_py, a1_py, a2_py, z_py):
	cdef str refa1 = refa1_py
	cdef str refa2 = refa2_py
	cdef str a1 = a1_py
	cdef str a2 = a2_py
	cdef float z = z_py
	if (refa1 + refa2 == ComplementCy(a2) + ComplementCy(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * z
	else:
		return z
