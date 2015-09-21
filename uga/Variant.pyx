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

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

import numpy as np
import pandas as pd
cimport numpy as np
cimport cython
import rpy2.robjects as ro
import pandas.rpy.common as py2r
from rpy2.rinterface import RRuntimeError
import math

cdef double CalcCallrate(np.ndarray[np.float64_t, ndim=1] x):
	cdef double xlen = len(x)
	cdef unsigned int ylen = len(x[~np.isnan(x)])
	if ylen > 0:
		return ylen/xlen
	else:
		return 0.0

cdef double CalcFreq(np.ndarray[np.float64_t, ndim=1] x):
	x = x[~np.isnan(x)]
	cdef unsigned int n = 2 * len(x)
	cdef double count = x.sum()
	return count / n if len(x) > 0 else float('nan')

cdef double CalcMAC(np.ndarray[np.float64_t, ndim=1] x):
	cdef double n
	x = x[~np.isnan(x)]
	if len(x) > 0 and set(x).issubset(set([0,1,2])):
		n = min(2.0 * len(x[x == 2]) + len(x[x == 1]),2.0 * len(x[x == 0]) + len(x[x == 1]))
	else:
		n = float('nan')
	return n

cdef double CalcFreqX(np.ndarray[np.float64_t, ndim=1] male, np.ndarray[np.float64_t, ndim=1] female):
	male = male[~np.isnan(male)]
	female = female[~np.isnan(female)]
	cdef unsigned int n = len(male) + 2 * len(female)
	cdef double count = (male.sum()/2) + female.sum()
	return count / n if len(male) > 0 and len(female) > 0 else float('nan')

cdef double CalcRsq(np.ndarray[np.float64_t, ndim=1] x):
	x = x[~np.isnan(x)]
	cdef double rsq = float('nan')
	if len(x) > 0 and not set(x).issubset(set([0,1,2])):
		rsq = x.var(ddof=1) / (2 * (x.mean() / 2) * (1 - (x.mean() / 2))) if x.mean() != 0 and x.mean() != 2 else 0.0
		return 1 / rsq if rsq > 1 else rsq
	else:
		return rsq

def CalcHWE(np.ndarray[np.float64_t, ndim=1] x):
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
		if not (x_hets < 0 and x_hom1 < 0 and x_hom2 < 0):
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

cpdef int GenerateFilterCode(np.float row_callrate, np.float row_freq, np.float row_mac, np.float row_rsq, np.float row_hwe, np.float row_freq_case=None, 
						np.float row_freq_ctrl=None, np.float miss_thresh=None, np.float maf_thresh=None, np.float maxmaf_thresh=None, np.float mac_thresh=None, 
							np.float rsq_thresh=None, np.float hwe_thresh=None, np.float hwe_maf_thresh=None, no_mono=True):
	cdef unsigned int filter = 0
	if (not miss_thresh is None and not np.isnan(row_callrate) and row_callrate < miss_thresh) or (not np.isnan(row_callrate) and row_callrate == 0) or np.isnan(row_callrate):
		filter += 10000
	if not np.isnan(row_freq): 
		if no_mono and ((row_freq == 0 or row_freq == 1) or ((row_freq_case is not None and not np.isnan(row_freq_case) and (row_freq_case == 0 or row_freq_case == 1)) or (row_freq_ctrl is not None and not np.isnan(row_freq_ctrl) and (row_freq_ctrl == 0 or row_freq_ctrl == 1)))):
			filter += 1000
		else:
			if ((	not maf_thresh is None
				and 
					(		row_freq < maf_thresh
						 or row_freq > 1-maf_thresh
					)
			   ) 
			  or
			   (	not maxmaf_thresh is None
				and    
					(		row_freq >= maxmaf_thresh
						and row_freq <= 1-maxmaf_thresh
					)
			   )):
				filter += 1000
	if not mac_thresh is None and not np.isnan(row_mac) and (row_mac < mac_thresh):
		filter += 100
	if not rsq_thresh is None and not np.isnan(row_rsq) and (row_rsq < rsq_thresh):
		filter += 10
	if not hwe_thresh is None and not hwe_maf_thresh is None and not np.isnan(row_hwe) and not np.isnan(row_freq) and ((row_freq <= 0.5 and row_freq > hwe_maf_thresh and row_hwe < hwe_thresh) or (row_freq > 0.5 and 1-row_freq > hwe_maf_thresh and row_hwe < hwe_thresh)):
		filter += 1
	return filter

def CalcVarStatsBinomial(marker_info, model_df, chr):
	model_df_nondup = model_df[model_df['___uga_nondup'].isin([1])]
	marker_info['callrate']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcCallrate(np.array(col).astype(np.float)), raw=True)
	if chr != 23:
		marker_info['freq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcFreq(np.array(col).astype(np.float)), raw=True)
		marker_info['freq.ctrl']=model_df_nondup[model_df_nondup['___uga_ctrl'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcFreq(np.array(col).astype(np.float)), raw=True)
		marker_info['freq.case']=model_df_nondup[model_df_nondup['___uga_case'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcFreq(np.array(col).astype(np.float)), raw=True)
		marker_info['rsq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcRsq(np.array(col).astype(np.float)), raw=True)
		marker_info['hwe']=model_df_nondup[model_df_nondup['___uga_ctrl'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcHWE(np.array(col).astype(np.float)), raw=True)
		marker_info['mac']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcMAC(np.array(col).astype(np.float)), raw=True)
	else:
		male_idx = model_df_nondup[model_df_nondup['___uga_male'].isin([1])].index.values
		female_idx = model_df_nondup[model_df_nondup['___uga_female'].isin([1])].index.values
		marker_info['freq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcFreqX(np.array(col[male_idx]).astype(np.float), np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['freq.ctrl']=model_df_nondup[model_df_nondup['___uga_ctrl'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcFreqX(np.array(col[male_idx]).astype(np.float), np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['freq.case']=model_df_nondup[model_df_nondup['___uga_case'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcFreqX(np.array(col[male_idx]).astype(np.float), np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['rsq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcRsq(np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['hwe']=model_df_nondup[model_df_nondup['___uga_ctrl'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcHWE(np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['mac']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcMAC(np.array(col).astype(np.float)), raw=True)
	return marker_info

def CalcVarStats(marker_info, model_df, chr):
	model_df_nondup = model_df[model_df['___uga_nondup'].isin([1])]
	marker_info['callrate']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcCallrate(np.array(col).astype(np.float)), raw=True)
	if chr != 23:
		marker_info['freq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcFreq(np.array(col).astype(np.float)), raw=True)
		marker_info['rsq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcRsq(np.array(col).astype(np.float)), raw=True)
		marker_info['hwe']=model_df_nondup[model_df_nondup['___uga_founder'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcHWE(np.array(col).astype(np.float)), raw=True)
		marker_info['mac']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcMAC(np.array(col).astype(np.float)), raw=True)
	else:
		male_idx = model_df_nondup[model_df_nondup['___uga_male'].isin([1])].index.values
		female_idx = model_df_nondup[model_df_nondup['___uga_female'].isin([1])].index.values
		marker_info['freq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcFreqX(np.array(col[male_idx]).astype(np.float), np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['rsq']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcRsq(np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['hwe']=model_df_nondup[model_df_nondup['___uga_founder'].isin([1])][marker_info['marker_unique']].apply(lambda col: CalcHWE(np.array(col[female_idx]).astype(np.float)), raw=True)
		marker_info['mac']=model_df_nondup[marker_info['marker_unique']].apply(lambda col: CalcMAC(np.array(col).astype(np.float)), raw=True)
	return marker_info

#cpdef FlipMinor(np.ndarray row, flip_a1, flip_a2, flip_freq, flip_effect, flip_or, flip_z):
#	cdef unsigned int i
#	a1 = row[flip_a1]
#	a2 = row[flip_a2]
#	row[flip_a1] = a2
#	row[flip_a2] = a1
#	for i in flip_freq:
#		row[i] = 1 - row[i]
#	for i in flip_effect:
#		row[i] = -1 * row[i]
#	for i in flip_or:
#		row[i] = 1 / row[i]
#	for i in flip_z:
#		row[i] = -1 * row[i]
#	return row

cpdef np.ndarray[np.float64_t, ndim=2] Bglm(np.ndarray[np.float64_t, ndim=2] a, formula, np.ndarray focus, np.ndarray markers):
	cdef unsigned int i = 0
	cdef unsigned int k = 0
	for k in xrange(len(markers)):
		cmd = "".join(['glm(',formula.replace('marker',markers[k]),',data=model_df,family=binomial())'])
		ro.globalenv['model_out'] = ro.r(cmd)
		if ro.r('model_out$converged')[0] and not ro.r('model_out$boundary')[0]:
			ro.globalenv['coef'] = ro.r('coef(summary(model_out))')
			i = 0
			for x in [z.replace('marker',markers[k]) for z in focus]:
				a[k,0 + i*5] = ro.r('coef["' + x + '","Estimate"]')[0]
				a[k,1 + i*5] = ro.r('coef["' + x + '","Std. Error"]')[0]
				a[k,2 + i*5] = math.exp(a[k,0 + i*5]) if not a[k,0 + i*5] > 709.782712893384 and not a[k,0 + i*5] < -709.782712893384 else float('nan')
				a[k,3 + i*5] = ro.r('coef["' + x + '","z value"]')[0]
				a[k,4 + i*5] = ro.r('coef["' + x + '","Pr(>|z|)"]')[0]
				i = i + 1
			a[k,0 + i*5] = len(ro.r('model_df$' + markers[k] + '[! is.na(model_df$' + markers[k] + ')]'))
			a[k,1 + i*5] = 1
		else:
			a[k,1 + i*5] = -2
	return a
"""
cpdef np.ndarray[np.float64_t, ndim=2] Bssmeta(np.ndarray[np.float64_t, ndim=2] a, formula, np.ndarray focus, np.ndarray markers, np.ndarray model_cols):
	ro.globalenv['markers'] = ro.StrVector(markers)
	ro.globalenv['model_cols'] = ro.StrVector(model_cols)
	ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': markers, 'gene': 'NA'}), strings_as_factors=False)
	ro.globalenv['z'] = ro.r('model_df[,names(model_df) %in% markers]')
	ro.globalenv['pheno'] = ro.r('model_df[model_cols]')
	cmd = 'prepScores(Z=z,formula=' + formula + ',SNPInfo=snp_info,data=pheno,family=binomial())'
	ro.globalenv['ps'] = ro.r(cmd)
	cmd = 'singlesnpMeta(ps,SNPInfo=snp_info)'
	try:
		ro.globalenv['result'] = ro.r('try(' + cmd + ',silent=TRUE)')
	except RRuntimeError:
		a[:,7] = -3
	else:
		a[:,0] = ro.r('result$beta')
		a[:,1] = ro.r('result$se')
		a[:,2] = np.exp(ro.r('result$beta'))
		a[:,3] = np.divide(ro.r('result$beta'), ro.r('result$se'))
		a[:,4] = ro.r('result$p')
		a[:,5] = ro.r('result$nmiss')
		a[:,6] = ro.r('result$ntotal')
		a[:,7] = 1
	return a
"""