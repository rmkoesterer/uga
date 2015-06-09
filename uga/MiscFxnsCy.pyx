from Model import *
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
from libc.string cimport strtok

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
	cdef double rsq = x.var(ddof=1) / (2 * (x.mean() / 2) * (1 - (x.mean() / 2))) if x.var(ddof=1) != 0 else 0
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
	if set(x).issubset(set([0,1,2])):
		obs_hets = np.sum(x == 1)
		obs_hom1 = np.sum(x == 2)
		obs_hom2 = np.sum(x == 0)
		if not (obs_hets < 0 and obs_hom1 < 0 and obs_hom2 < 0):
			rare   = 2*min(obs_hom1,obs_hom2)+obs_hets
			common = 2*max(obs_hom1,obs_hom2)+obs_hets
			if not rare:
				return 1.0

			# expected heterogygotes
			hets = rare*common/(rare+common)

			# check for non-matching parity of rare and hets
			if rare%2 != hets%2:
				hets += 1

			# expected rare and common homozygotes
			hom_r = (rare-hets)/2
			hom_c = (common-hets)/2

			# initialize heterozygote probability vector
			probs = [0]*(rare/2+1)

			# Set P(expected hets)=1
			probs[hets/2] = 1.0

			# fill in relative probabilities for less than the expected hets
			for i,h in enumerate(xrange(hets,1,-2)):
				probs[h/2-1] = probs[h/2]*h*(h-1) / (4*(hom_r+i+1)*(hom_c+i+1))

			# fill in relative probabilities for greater than the expected hets
			for i,h in enumerate(xrange(hets,rare-1,2)):
				probs[h/2+1] = probs[h/2]*4*(hom_r-i)*(hom_c-i) / ((h+1)*(h+2))

			# compute the pvalue by summing the probabilities <= to that of the
			# observed number of heterogygotes and normalize by the total
			p_obs = probs[obs_hets/2]
			p_hwe = sum(p for p in probs if p <= p_obs)/sum(probs)
			return p_hwe
		else:
			return float('NaN')
	else:
		return float('NaN')

@cython.boundscheck(False)
@cython.wraparound(False)
def GenerateFilterCodeCy(np.float row_callrate, np.float row_freq, np.float row_freq_unrel, np.float row_rsq, np.float row_rsq_unrel, np.float row_hwe, np.float row_hwe_unrel, np.float miss_thresh=None, np.float maf_thresh=None, np.float maf_max_thresh=None, np.float rsq_thresh=None, np.float hwe_thresh=None, no_mono=True):
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
			   (	not maf_max_thresh is None
				and 
					(	   
						(		row_freq >= maf_max_thresh
							and row_freq <= 1-maf_max_thresh
						)
					)
				and
					(
						(		row_freq_unrel >= maf_max_thresh
							and row_freq_unrel <= 1-maf_max_thresh
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