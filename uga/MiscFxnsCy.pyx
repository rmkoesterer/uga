from Model import *
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
from libc.string cimport strtok

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
def GetDelimiterCy(str delimiter):
	if delimiter == 'tab':
		delimiter = '\t'
	elif delimiter == 'space':
		delimiter = ' '
	else:
		delimiter = ','
	return delimiter

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
	return count / n if len(x) > 0 else float('NaN')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcFreq23Cy(np.ndarray male, np.ndarray female):
	male = male[~np.isnan(male)]
	female = female[~np.isnan(female)]
	cdef unsigned int n = len(male) + 2 * len(female)
	cdef double count = (male.sum()/2) + female.sum()
	return count / n if len(male) > 0 and len(female) > 0 else float('NaN')

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CalcRsqCy(np.ndarray x):
	x = x[~np.isnan(x)]
	cdef double rsq = x.var(ddof=1) / (2 * (x.mean() / 2) * (1 - (x.mean() / 2)))
	return 1 / rsq if rsq > 1 else rsq

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
	if list(set(x)) == [0,1,2]:
		obs_hets = np.sum(x == 1)
		obs_hom1 = np.sum(x == 2)
		obs_hom2 = np.sum(x == 0)
		if not (obs_hets == 0 and obs_hom1 == 0 and obs_hom2 == 0):
			obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
			obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

			rare_copies = 2 * obs_homr + obs_hets
			genotypes   = obs_hets + obs_homc + obs_homr

			het_probs = [0.0] * (rare_copies + 1)

			#start at midpoint
			mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)

			#check to ensure that midpoint and rare alleles have same parity
			if (rare_copies & 1) ^ (mid & 1):
				mid += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			het_probs[mid] = 1.0
			sum = float(het_probs[mid])

			for curr_hets in xrange(mid, 1, -2):
				het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
				sum += het_probs[curr_hets - 2]
				# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
				curr_homr += 1
				curr_homc += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			for curr_hets in xrange(mid, rare_copies - 1, 2):
				het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
				sum += het_probs[curr_hets + 2]
				#add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
				curr_homr -= 1
				curr_homc -= 1

			for i in xrange(0, rare_copies + 1):
				het_probs[i] /= sum

			#alternate p-value calculation for p_hi/p_lo
			p_hi = float(het_probs[obs_hets])
			for i in xrange(obs_hets, rare_copies+1):
				p_hi += het_probs[i]

			p_lo = float(het_probs[obs_hets])
			for i in xrange(obs_hets-1, -1, -1):
				p_lo += het_probs[i]

			p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

			p_hwe = 0.0
			#  p-value calculation for p_hwe
			for i in xrange(0, rare_copies + 1):
				if het_probs[i] > het_probs[obs_hets]:
					continue
				p_hwe += het_probs[i]

			p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

			return p_hwe
		else:
			return float('NaN')
	else:
		return float('NaN')

@cython.boundscheck(False)
@cython.wraparound(False)
def GenerateFilterCodeCy(np.float row_callrate, np.float row_freq, np.float row_freq_unrel, np.float row_rsq, np.float row_rsq_unrel, np.float row_hwe, np.float row_hwe_unrel, np.float miss_thresh, np.float freq_thresh, np.float max_freq_thresh, np.float rsq_thresh, np.float hwe_thresh, no_mono=True):
	cdef unsigned int filter = 0
	if (not miss_thresh is None and not np.isnan(row_callrate) and row_callrate < miss_thresh) or (not np.isnan(row_callrate) and row_callrate == 0) or np.isnan(row_callrate):
		filter += 1000
	if not np.isnan(row_freq): 
		if no_mono and (row_freq == 0 or row_freq == 1):
			filter += 100
		else:
			if ((	not freq_thresh is None
				and 
					(		row_freq < freq_thresh
						 or row_freq > 1-freq_thresh
					)
				and 
					(		row_freq_unrel < freq_thresh
						 or row_freq_unrel > 1-freq_thresh
					)
			   ) 
			  or
			   (	not max_freq_thresh is None
				and 
					(	   
						(		row_freq >= max_freq_thresh
							and row_freq <= 1-max_freq_thresh
						)
					)
				and
					(
						(		row_freq_unrel >= max_freq_thresh
							and row_freq_unrel <= 1-max_freq_thresh
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