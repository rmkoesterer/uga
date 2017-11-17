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

import pandas as pd
import numpy as np
import readline
from re import split as re_split
import sys
import Process
cimport numpy as np
cimport cython
import Geno
import Variant
cimport Variant
import Fxns
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.rinterface import RRuntimeError
from __version__ import version
import time
import math
import scipy.stats as scipy
import numpy.lib.recfunctions as recfxns
import logging
pandas2ri.activate()
ro.r('options(warn=-1)')
ro.r('options(stringsAsFactors=FALSE)')
ro.r('options(na.action=na.omit)')
ro.r('suppressMessages(library(R.utils))')

module_logger = logging.getLogger("Model")

cdef class Model(object):
	cdef public unsigned int case_code, ctrl_code, tbx_start, tbx_end, \
								male, female, nobs, nunique, nfounders, nlongitudinal, \
								nfamilies, nunrelated, ncases, nctrls, nmales, nfemales
	cdef public np.ndarray cases_idx, ctrls_idx, male_cases_idx, male_ctrls_idx, female_cases_idx, female_ctrls_idx, \
							model_cols, results_header, calc_hwe_idx, \
							variant_stats, results, unique_idx, founders_idx, founders_ctrls_idx, male_idx, female_idx, \
							geno_cases_idx, geno_ctrls_idx, geno_male_cases_idx, geno_male_ctrls_idx, geno_female_cases_idx, geno_female_ctrls_idx, \
							geno_calc_hwe_idx, geno_unique_idx, geno_male_idx, geno_female_idx
	cdef public bytes fxn, format, pheno, variants_file, type, samples_file, drop_file, keep_file, \
						iid, fid, matid, patid, sex, sep, a1, a2
	cdef public str metadata, metadata_cc, family, formula, focus, dep_var, interact, random_effects, covars
	cdef public object pheno_df, variants, out, results_dtypes, pedigree, drop, keep
	cdef public bint all_founders, reverse, reml, kr
	def __cinit__(self, fxn, format, variants_file, pheno, type, iid, fid, 
					case_code = None, ctrl_code = None, all_founders = False, dep_var = None, covars = None, interact = None, random_effects = None, reml = False, kr = False, reverse = False, samples_file = None, 
					drop_file = None, keep_file = None, matid = None, patid = None, sex = None, male = 1, female = 2, sep = 'tab', **kwargs):
		super(Model, self).__init__(**kwargs)
		logger = logging.getLogger("Model.Model.__cinit__")
		logger.debug("initialize model")
		self.out = None
		self.fxn = fxn
		self.dep_var = dep_var
		self.covars = covars
		self.interact = interact
		self.random_effects = random_effects
		self.reml = reml
		self.kr = kr
		self.reverse = reverse
		self.format = format
		self.variants_file = variants_file
		self.samples_file = samples_file
		self.drop_file = drop_file
		self.keep_file = keep_file
		self.pheno = pheno
		self.type = type
		self.iid = iid
		self.fid = fid
		self.case_code = case_code
		self.ctrl_code = ctrl_code
		self.matid = matid
		self.patid = patid
		self.sex = sex
		self.male = male
		self.female = female
		self.sep = sep
		self.all_founders = all_founders
		self.family = 'gaussian'

		# parse formula into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		# snv				snv					placeholder for snv (will be replaced with each snv during iteration and can be used in an interaction)

		self.model_cols = np.array([self.dep_var])
		self.model_cols = np.append(self.model_cols,np.array([x.replace('factor(','').replace(')','') for x in self.covars.split('+')])) if self.covars else self.model_cols
		self.model_cols = np.append(self.model_cols,np.array(self.interact.replace('factor(','').replace(')',''))) if self.interact else self.model_cols
		self.model_cols = np.append(self.model_cols,np.array(self.random_effects.split('+'))) if self.random_effects else self.model_cols
		self.model_cols = np.unique(self.model_cols)

		self.male_idx = np.array([])
		self.female_idx = np.array([])
		self.cases_idx = np.array([])
		self.ctrls_idx = np.array([])
		self.male_cases_idx = np.array([])
		self.male_ctrls_idx = np.array([])
		self.female_cases_idx = np.array([])
		self.female_ctrls_idx = np.array([])
		self.founders_ctrls_idx = np.array([])
		self.calc_hwe_idx = np.array([])

		self.geno_male_idx = np.array([])
		self.geno_female_idx = np.array([])
		self.geno_cases_idx = np.array([])
		self.geno_ctrls_idx = np.array([])
		self.geno_male_cases_idx = np.array([])
		self.geno_male_ctrls_idx = np.array([])
		self.geno_female_cases_idx = np.array([])
		self.geno_female_ctrls_idx = np.array([])
		self.geno_calc_hwe_idx = np.array([])

		try:
			self.variants = getattr(Geno,self.format.capitalize())(self.variants_file, self.samples_file)
		except Process.Error as err:
			raise Process.Error(err.msg)
		else:
			print "extracting model fields from pheno file and reducing to complete observations ..."
			p_names = (self.fid,self.iid) + tuple(x for x in [self.matid, self.patid] if x is not None)
			p_names = p_names + tuple(x for x in self.model_cols if x not in [self.fid,self.iid,self.matid,self.patid])
			p_names = p_names + (self.sex,) if self.sex not in self.model_cols else p_names
			p_dtypes = ('|S100', '|S100') + tuple('|S100' for x in [self.matid, self.patid] if x is not None)
			p_dtypes = p_dtypes + tuple('f8' for x in self.model_cols if x not in [self.fid,self.iid,self.matid,self.patid])
			p_dtypes = p_dtypes + ('f8',) if self.sex not in self.model_cols else p_dtypes
			dtypes = dict(zip(p_names, p_dtypes))
			try:
				self.pheno_df = np.genfromtxt(fname=self.pheno, delimiter=Fxns.get_delimiter(self.sep), dtype=dtypes.values(), names=True, usecols=dtypes.keys())
			except:
				raise Process.Error("unable to load phenotype file " + self.pheno + " with columns " + ', '.join(p_names) + "; " + str(sys.exc_info()[0]))
			for x in self.model_cols:
				if x in self.pheno_df.dtype.names:
					print "   model column %s found" % (x)
				else:
					raise Process.Error("column " + x + " not found in phenotype file " + self.pheno)
			if self.fid in self.pheno_df.dtype.names:
				print "   fid column %s found" % self.fid
			else:
				raise Process.Error("column " + self.fid + " not found in phenotype file " + self.pheno)
			if self.iid in self.pheno_df.dtype.names:
				print "   iid column %s found" % self.iid
			else:
				raise Process.Error("column " + self.iid + " not found in phenotype file " + self.pheno)
			if self.matid is not None:
				if self.matid in self.pheno_df.dtype.names:
					print "   matid column %s found" % self.matid
				else:
					raise Process.Error("column " + self.matid + " not found in phenotype file " + self.pheno)
			if self.patid is not None:
				if self.patid in self.pheno_df.dtype.names:
					print "   patid column %s found" % self.patid
				else:
					raise Process.Error("column " + self.patid + " not found in phenotype file " + self.pheno)
			if self.sex in self.pheno_df.dtype.names:
				print "   sex column %s found" % self.sex
			else:
				raise Process.Error("column " + self.sex + " not found in phenotype file " + self.pheno)
			if self.matid is not None and self.patid is not None and self.sex is not None and set([self.iid,self.fid,self.matid,self.patid,self.sex]) <= set(self.pheno_df.dtype.names):
				print "   extracting pedigree"
				self.pedigree = pd.DataFrame(self.pheno_df[[self.iid,self.fid,self.matid,self.patid,self.sex]])
				self.pedigree.loc[(self.pedigree[self.sex] != self.male) & (self.pedigree[self.sex] != self.female),self.sex]=-999
				self.pedigree[self.sex].replace({self.male: 1, self.female: 2, -999: 3})
				self.pedigree[self.matid].replace({'0': 'NA'})
				self.pedigree[self.patid].replace({'0': 'NA'})
			for x in [y for y in dtypes if dtypes[y] == 'f8']:
				self.pheno_df = self.pheno_df[~np.isnan(self.pheno_df[x])]
			for x in [y for y in dtypes if dtypes[y] == '|S100']:
				self.pheno_df = self.pheno_df[~(self.pheno_df[x] == 'NA')]
			self.pheno_df = self.pheno_df[np.in1d(self.pheno_df[self.iid],np.intersect1d(self.pheno_df[self.iid],self.variants.samples))]
			print "phenotype file and data file contain " + str(self.pheno_df.shape[0]) + " common samples"
			if self.drop_file is not None:
				try:
					self.drop=np.genfromtxt(fname=self.drop_file, dtype='object')
				except:
					raise Process.Error("unable to load sample drop file " + self.drop_file)
				print "dropping " + str(len([a for a in np.in1d(self.pheno_df[self.iid],self.drop) if a])) + " samples from file " + self.drop_file
				self.pheno_df = self.pheno_df[np.in1d(self.pheno_df[self.iid],self.drop,invert=True)]
			if self.keep_file is not None:
				try:
					self.keep=np.genfromtxt(fname=self.keep_file, dtype='object')
				except:
					raise Process.Error("unable to load sample keep file " + self.keep_file)
				print "keeping " + str(len([a for a in np.in1d(self.pheno_df[self.iid],self.keep) if a])) + " samples from file " + self.keep_file
				self.pheno_df = self.pheno_df[np.in1d(self.pheno_df[self.iid],self.keep)]
			if self.pheno_df.shape[0] > 0:
				iids_unique, iids_counts = np.unique(self.pheno_df[self.iid], return_counts=True)
				fids_unique, fids_counts = np.unique(self.pheno_df[self.fid], return_counts=True)
				self.unique_idx = np.in1d(self.pheno_df[self.iid],iids_unique)
				if self.all_founders or self.matid is None or self.patid is None:
					self.founders_idx = self.unique_idx.copy()
				else:
					self.founders_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.unique_idx][(self.pheno_df[self.unique_idx][self.matid] == '0') & (self.pheno_df[self.unique_idx][self.patid] == '0')][self.iid])
				self.calc_hwe_idx = self.founders_idx.copy()
				self.nobs = self.pheno_df.shape[0]
				self.nunique = len(iids_unique)
				self.nlongitudinal = len(iids_counts[iids_counts != 1])
				self.nfamilies = len(fids_counts[fids_counts > 1])
				self.nunrelated = len(fids_counts[fids_counts == 1])
				print "data summary ..."
				print "   " + str(self.nobs) + " total observations"
				print "   " + str(self.nunique) + " unique samples"
				print "   " + str(self.nlongitudinal) + " multiple observation samples"
				print "   " + str(self.nfamilies) + " families"
				print "   " + str(self.nunrelated) + " unrelated samples"
				if self.matid is not None and self.patid is not None:
					self.nfounders = self.pheno_df[self.unique_idx][(self.pheno_df[self.unique_idx][self.matid] == '0') & (self.pheno_df[self.unique_idx][self.patid] == '0')].shape[0]
					print "   " + str(self.nfounders) + " founders"
				else:
					self.nfounders = self.nunique
					print "   " + str(self.nfounders) + " founders"
				if self.male is not None:
					self.nmales = self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.sex] == self.male].shape[0]
					self.male_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.sex] == self.male][self.iid])
					print "   " + str(self.nmales) + " male"
				if self.female is not None:
					self.nfemales = self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.sex] == self.female].shape[0]
					self.female_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.sex] == self.female][self.iid])
					print "   " + str(self.nfemales) + " female"
				if len(np.unique(self.pheno_df[self.dep_var])) == 2:
					self.family = 'binomial'
					self.ncases = self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.dep_var] == self.case_code].shape[0]
					self.cases_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.dep_var] == self.case_code][self.iid])
					self.male_cases_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.male_idx][self.pheno_df[self.male_idx][self.dep_var] == self.case_code][self.iid])
					self.female_cases_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.female_idx][self.pheno_df[self.female_idx][self.dep_var] == self.case_code][self.iid])
					print "   " + str(self.ncases) + " cases"
					self.nctrls = self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.dep_var] == self.ctrl_code].shape[0]
					self.ctrls_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.unique_idx][self.pheno_df[self.unique_idx][self.dep_var] == self.ctrl_code][self.iid])
					self.founders_ctrls_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.founders_idx][self.pheno_df[self.founders_idx][self.dep_var] == self.ctrl_code][self.iid])
					self.calc_hwe_idx = self.founders_ctrls_idx.copy()
					self.male_ctrls_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.male_idx][self.pheno_df[self.male_idx][self.dep_var] == self.ctrl_code][self.iid])
					self.female_ctrls_idx = np.in1d(self.pheno_df[self.iid],self.pheno_df[self.female_idx][self.pheno_df[self.female_idx][self.dep_var] == self.ctrl_code][self.iid])
					print "   " + str(self.nctrls) + " controls"
					self.pheno_df[self.dep_var][self.cases_idx] = 1
					self.pheno_df[self.dep_var][self.ctrls_idx] = 0

				if len(np.unique(self.pheno_df[self.dep_var])) == 1:
					raise Process.Error("phenotype consists of only a single unique value")
			else:
				raise Process.Error("phenotype file and data file contain no common samples")

		self.metadata = '## source: uga' + version + '\n' + \
						'## date: ' + time.strftime("%Y%m%d") + '\n' + \
						'## method: ' + self.fxn + '\n' + \
						'## unique samples: ' + str(self.nunique) + '\n' + \
						'## multiple observation samples: ' + str(self.nlongitudinal) + '\n' + \
						'## families: ' + str(self.nfamilies) + '\n' + \
						'## unrelated samples: ' + str(self.nunrelated) + '\n' + \
						'## founders: ' + str(self.nfounders) + '\n' + \
						'## males: ' + str(self.nmales) + '\n' + \
						'## females: ' + str(self.nfemales)

		self.metadata_cc = '## cases: ' + str(self.ncases) + '\n' + \
							'## controls: ' + str(self.nctrls)

	def get_region(self, region, group_id = None):
		logger = logging.getLogger("Model.Model.get_region")
		logger.debug("get_region " + region)
		try:
			self.variants.get_region(region, group_id)
		except:
			print "region " + region + " returned empty"
			raise

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_variant_stats(self):
		logger = logging.getLogger("Model.Model.calc_variant_stats")
		logger.debug("calc_variant_stats")
		if self.family == "binomial":
			self.variant_stats = np.zeros((self.variants.info.shape[0],1), dtype=[('filter','uint32'),('mac','f8'),('callrate','f8'),('freq','f8'),('freq.case','f8'),('freq.ctrl','f8'),('rsq','f8'),('hwe','f8'),('n','f8')])
		else:
			self.variant_stats = np.zeros((self.variants.info.shape[0],1), dtype=[('filter','uint32'),('mac','f8'),('callrate','f8'),('freq','f8'),('rsq','f8'),('hwe','f8'),('n','f8')])
		cdef unsigned int i
		self.geno_unique_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.unique_idx])
		self.geno_calc_hwe_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.calc_hwe_idx])
		self.geno_male_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.male_idx])
		self.geno_female_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.female_idx])
		if self.family == "binomial":
			self.geno_cases_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.cases_idx])
			self.geno_ctrls_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.ctrls_idx])
			self.geno_male_cases_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.male_cases_idx])
			self.geno_male_ctrls_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.male_ctrls_idx])
			self.geno_female_cases_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.female_cases_idx])
			self.geno_female_ctrls_idx = np.in1d(self.variants.data[:,0],self.pheno_df[self.iid][self.female_ctrls_idx])
		if self.variants.chr == 23:
			for i in xrange(self.variants.info.shape[0]):
				self.variant_stats['mac'][i] = Variant.calc_mac(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['callrate'][i] = Variant.calc_callrate(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['freq'][i] = Variant.calc_freqX(male=self.variants.data[self.geno_male_idx,i+1].astype('float64'), female=self.variants.data[self.geno_female_idx,i+1].astype('float64')) if len(self.geno_male_idx) > 0 and len(self.geno_female_idx) > 0 else np.nan
				if self.family == "binomial":
					self.variant_stats['freq.case'][i] = Variant.calc_freqX(male=self.variants.data[self.geno_male_cases_idx,i+1].astype('float64'), female=self.variants.data[self.geno_female_cases_idx,i+1].astype('float64')) if len(self.geno_male_cases_idx) > 0 and len(self.geno_female_cases_idx) > 0 else np.nan
					self.variant_stats['freq.ctrl'][i] = Variant.calc_freqX(male=self.variants.data[self.geno_male_ctrls_idx,i+1].astype('float64'), female=self.variants.data[self.geno_female_ctrls_idx,i+1].astype('float64')) if len(self.geno_male_ctrls_idx) > 0 and len(self.geno_female_ctrls_idx) > 0 else np.nan
				self.variant_stats['rsq'][i] = Variant.calc_rsq(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['hwe'][i] = Variant.calc_hwe(self.variants.data[self.geno_calc_hwe_idx,i+1].astype('float64')) if len(self.geno_calc_hwe_idx) > 0 else np.nan
				self.variant_stats['n'][i] = round(self.variant_stats['callrate'][i] * self.nunique)
		else:
			for i in xrange(self.variants.info.shape[0]):
				self.variant_stats['mac'][i] = Variant.calc_mac(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['callrate'][i] = Variant.calc_callrate(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['freq'][i] = Variant.calc_freq(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				if self.family == "binomial":
					self.variant_stats['freq.case'][i] = Variant.calc_freq(self.variants.data[self.geno_cases_idx,i+1].astype('float64')) if len(self.geno_cases_idx) > 0 else np.nan
					self.variant_stats['freq.ctrl'][i] = Variant.calc_freq(self.variants.data[self.geno_ctrls_idx,i+1].astype('float64')) if len(self.geno_ctrls_idx) > 0 else np.nan
				self.variant_stats['rsq'][i] = Variant.calc_rsq(self.variants.data[self.geno_unique_idx,i+1].astype('float64'))
				self.variant_stats['hwe'][i] = Variant.calc_hwe(self.variants.data[self.geno_calc_hwe_idx,i+1].astype('float64')) if len(self.geno_calc_hwe_idx) > 0 else np.nan
				self.variant_stats['n'][i] = round(self.variant_stats['callrate'][i] * self.nunique)

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef filter(self, np.float miss_thresh=None, np.float maf_thresh=None, np.float maxmaf_thresh=None, np.float mac_thresh=None, 
								np.float rsq_thresh=None, np.float hwe_thresh=None, np.float hwe_maf_thresh=None, allow_mono=False):
		logger = logging.getLogger("Model.Model.filter")
		logger.debug("filter")
		cdef unsigned int i
		for i in xrange(self.variant_stats.shape[0]):
			if (not miss_thresh is None and not np.isnan(self.variant_stats['callrate'][i]) and self.variant_stats['callrate'][i] < miss_thresh) or (not np.isnan(self.variant_stats['callrate'][i]) and self.variant_stats['callrate'][i] == 0) or np.isnan(self.variant_stats['callrate'][i]):
				self.variant_stats['filter'][i] += 10000
			if not np.isnan(self.variant_stats['freq'][i]): 
				if not allow_mono and ((self.variant_stats['freq'][i] == 0 or self.variant_stats['freq'][i] == 1) or 
								(self.family == "binomial" and 
									((self.variant_stats['freq.case'][i] is not None and not np.isnan(self.variant_stats['freq.case'][i]) and 
									(self.variant_stats['freq.case'][i] == 0 or self.variant_stats['freq.case'][i] == 1)) or 
									(self.variant_stats['freq.ctrl'][i] is not None and not np.isnan(self.variant_stats['freq.ctrl'][i]) and 
									(self.variant_stats['freq.ctrl'][i] == 0 or self.variant_stats['freq.ctrl'][i] == 1))))):
					self.variant_stats['filter'][i] += 1000
				else:
					if ((	not maf_thresh is None
						and 
							(		self.variant_stats['freq'][i] < maf_thresh
								or self.variant_stats['freq'][i] > 1-maf_thresh
							)
					) 
					or
					(	not maxmaf_thresh is None
						and    
							(		self.variant_stats['freq'][i] >= maxmaf_thresh
								and self.variant_stats['freq'][i] <= 1-maxmaf_thresh
							)
					)):
						self.variant_stats['filter'][i] += 1000
			if not mac_thresh is None and not np.isnan(self.variant_stats['mac'][i]) and (self.variant_stats['mac'][i] < mac_thresh):
				self.variant_stats['filter'][i] += 100
			if not rsq_thresh is None and not np.isnan(self.variant_stats['rsq'][i]) and (self.variant_stats['rsq'][i] < rsq_thresh):
				self.variant_stats['filter'][i] += 10
			if not hwe_thresh is None and not hwe_maf_thresh is None and not np.isnan(self.variant_stats['hwe'][i]) and not np.isnan(self.variant_stats['freq'][i]) and ((self.variant_stats['freq'][i] <= 0.5 and self.variant_stats['freq'][i] > hwe_maf_thresh and self.variant_stats['hwe'][i] < hwe_thresh) or (self.variant_stats['freq'][i] > 0.5 and 1-self.variant_stats['freq'][i] > hwe_maf_thresh and self.variant_stats['hwe'][i] < hwe_thresh)):
				self.variant_stats['filter'][i] += 1

cdef class SnvModel(Model):
	cdef public str metadata_snv, metadata_snv_cc
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.SnvModel.__cinit__")
		logger.debug("initialize SnvModel")
		super(SnvModel, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 1
		if self.family == "binomial":
			self.results_header = np.array(['chr','pos','id','a1','a2','filter','callrate','rsq','hwe','n','mac','freq','freq.case','freq.ctrl'])
		else:
			self.results_header = np.array(['chr','pos','id','a1','a2','filter','callrate','rsq','hwe','n','mac','freq'])

		self.metadata_snv = '## chr: chromosome' + '\n' + \
							'## pos: chromosomal position' + '\n' + \
							'## id: snv name' + '\n' + \
							'## a1: reference (coded) allele used for stat calculations' + '\n' + \
							'## a2: alternate (non-coded) allele' + '\n' + \
							'## filter: filter code (+1 = failed hwe, +10 = failed rsq, +100 = failed mac, +1000 = failed maf, +10000 = failed miss' + '\n' + \
							'## callrate: callrate' + '\n' + \
							'## rsq: imputation info metric (calculated only for variants exhibiting probabilistic genotypes, otherwise coded as NA)' + '\n' + \
							'## hwe: Hardy Weinberg p-value (calculated only if dosages are in [0,1,2], otherwise coded as NA; calculated in founders only if pedigree information available, otherwise calculated in all samples)' + '\n' + \
							'## mac: minor allele count' + '\n' + \
							'## n: samples analyzed' + '\n' + \
							'## freq: reference (coded) allele frequency'

		self.metadata_snv_cc = '## freq.case: reference (coded) allele frequency in cases' + '\n' + \
								'## freq.ctrl: reference (coded) allele frequency in controls'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Model.SnvModel.get_snvs")
		logger.debug("get_snvs")
		try:
			self.variants.get_snvs(buffer)
		except:
			logger.debug("found 0 variants for " + str(self.variants.data.shape[0]) + " samples")
			raise
		else:
			try:
				self.variants.data = self.variants.data[np.in1d(self.variants.data[:,0],np.intersect1d(self.pheno_df[self.iid],self.variants.data[:,0]))]
			except:
				raise Process.Error("ped file and data file contain no common samples")
			else:
				logger.debug("found " + str(self.variants.data.shape[1]) + " variants for " + str(self.variants.data.shape[0]) + " samples")
				self.calc_variant_stats()

cdef class SnvgroupModel(Model):
	cdef public str snvgroup_map, metadata_gene
	cdef public unsigned int cmac
	def __cinit__(self, snvgroup_map = None, cmac = None, **kwargs):
		logger = logging.getLogger("Model.SnvgroupModel.__cinit__")
		logger.debug("initialize SnvgroupModel")
		self.snvgroup_map = snvgroup_map
		self.cmac = cmac if cmac is not None else 1
		super(SnvgroupModel, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 2
		self.results_header = np.array(['chr','start','end','id','total','passed'])

		if self.snvgroup_map is not None:
			try:
				self.variants.load_snvgroup_map(self.snvgroup_map)
			except Process.Error as err:
				raise Process.Error(err.msg)

		self.metadata_gene = '## chr: chromosome' + '\n' + \
								'## start: start chromosomal position' + '\n' + \
								'## end: end chromosomal position' + '\n' + \
								'## id: snv group name (ie. gene name)'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvgroup(self, int buffer, group_id):
		logger = logging.getLogger("Model.SnvgroupModel.get_snvgroup")
		logger.debug("get_snvgroup")
		try:
			self.variants.get_snvgroup(buffer, group_id)
		except:
			raise
		else:
			try:
				self.variants.data = self.variants.data[np.in1d(self.variants.data[:,0],np.intersect1d(self.pheno_df[self.iid],self.variants.data[:,0]))]
			except:
				raise Process.Error("ped file and data file contain no common samples")
			else:
				self.calc_variant_stats()

cdef class Score(SnvModel):
	cdef public bint adjust_kinship
	def __cinit__(self, adjust_kinship = False, **kwargs):
		logger = logging.getLogger("Model.Score.__cinit__")
		logger.debug("initialize Score model")
		self.adjust_kinship = adjust_kinship
		super(Score, self).__init__(**kwargs)
		print "setting score test family option to " + self.family

		self.formula = self.dep_var + '~' + self.covars if self.covars is not None else self.dep_var + '~1'
		print "formula: " + self.formula

		self.results_dtypes=[('err','f8'),('nmiss','f8'),('ntotal','f8'),('effect','f8'),('stderr','f8'),('or','f8'),('p','f8')]
		if self.family == 'binomial':
			self.results_header = np.append(self.results_header,np.array(['err','nmiss','ntotal','effect','stderr','or','p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','nmiss','ntotal','effect','stderr','p']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')
		if self.adjust_kinship:
			print "loading R package kinship2"
			ro.r('suppressMessages(library(kinship2))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + '## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + '## adjust kinship: True' if self.adjust_kinship else self.metadata + '\n' + '## adjust kinship: False'
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed singlesnpMeta)' + '\n' + \
						'## nmiss: number of missing samples' + '\n' + \
						'## ntotal: total number of included samples' + '\n' + \
						'## effect: effect size' + '\n' + \
						'## stderr: standard error'
		self.metadata = self.metadata + '\n' + '## *.or: odds ratio (exp(effect), not provided by seqMeta)' if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + '## *.p: p-value' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=self.results_dtypes)
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['id_unique'][passed]))
			ro.globalenv['snp_info'] = pd.DataFrame({'Name': list(self.variants.info['id_unique'][passed]), 'gene': 'NA'})
			ro.r('snp_info$Name<-as.character(snp_info$Name)')
			ro.r('snp_info$gene<-as.character(snp_info$gene)')
			ro.r('z<-data.matrix(model_df[,names(model_df) %in% variants])')
			if len(passed) == 1:
				ro.r('colnames(z)<-"' + self.variants.info['id_unique'][passed][0] + '"')
			if self.adjust_kinship:
				ro.globalenv['ped'] = self.pedigree
				ro.globalenv['kins'] = ro.r("kinship(pedigree(famid=ped$" + self.fid + ",id=ped$" + self.iid + ",dadid=ped$" + self.patid + ",momid=ped$" + self.matid + ",sex=ped$" + self.sex + ",missid='NA'))")
				cmd = "prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "',kins=kins,sparse=FALSE)"
			else:
				cmd = "prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "')"
			try:
				ro.globalenv['ps'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][passed] = 2
				raise Process.Error(rerr.message)
			cmd = 'singlesnpMeta(ps,SNPInfo=snp_info)'
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][passed] = 3
				raise Process.Error(rerr.message)
			else:
				ro.r('result$err<-0')
				ro.r('result$err[! is.finite(result$beta) | ! is.finite(result$se) | ! is.finite(result$p) | result$p == 0]<-1')
				ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$ntotal[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$beta[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$se[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
				self.results['err'][passed] = np.array(ro.r('result$err'))[:,None]
				self.results['nmiss'][passed] = np.array(ro.r('result$nmiss'))[:,None]
				self.results['ntotal'][passed] = np.array(ro.r('result$ntotal'))[:,None]
				self.results['effect'][passed] = np.array(ro.r('result$beta'))[:,None]
				self.results['stderr'][passed] = np.array(ro.r('result$se'))[:,None]
				self.results['or'][passed] = np.array(np.exp(ro.r('result$beta')))[:,None]
				self.results['p'][passed] = np.array(ro.r('result$p'))[:,None]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Gee(SnvModel):
	cdef public str corstr
	def __cinit__(self, corstr = None, **kwargs):
		logger = logging.getLogger("Model.Gee.__cinit__")
		logger.debug("initialize Gee model")
		self.corstr = corstr if corstr is not None else 'exchangeable'
		super(Gee, self).__init__(**kwargs)
		print "setting gee test family option to " + self.family

		if self.reverse:
			if self.interact is not None:
				self.formula = '___snv___~' + self.dep_var + '*' + self.interact
				self.focus = self.dep_var + ':' + self.interact
			else:
				self.formula = '___snv___~' + self.dep_var
				self.focus = self.dep_var
			self.family = 'gaussian'
		else:
			if self.interact is not None:
				self.formula = self.dep_var + '~___snv___*' + self.interact
				self.focus = '___snv___:' + self.interact
			else:
				self.formula = self.dep_var + '~___snv___'
				self.focus = '___snv___'
		self.formula = self.formula + '+' + self.covars if self.covars is not None else self.formula
		print "formula: " + self.formula

		self.results_dtypes = [('err','f8'),('effect','f8'),('stderr','f8'),('or','f8'),('wald','f8'),('p','f8')]
		if self.family == 'binomial':
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','or','wald','p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','wald','p']))

		print "loading R package geepack"
		ro.r('suppressMessages(library(geepack))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + '## formula: ' + self.formula
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						'## corstr: ' + self.corstr + '\n' + \
						self.metadata_snv + '\n' + \
						'## err: error code (0: no error, 1: geeglm error reported, 2: infinite value or zero p-value detected, 3: geeglm failed)' + '\n' + \
						'## effect: effect size' + '\n' + \
						'## stderr: standard error'
		self.metadata = self.metadata + '\n' + '## or: odds ratio (included only if binomial family)' if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## wald: wald-statistic' + '\n' + \
						'## p: p-value' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=self.results_dtypes)
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			ro.r('model_df$' + self.fid + '<-as.factor(model_df$' + self.fid + ')')
			ro.r('model_df<-model_df[order(model_df$' + self.fid + '),]')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['id_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [self.fid] + [vu]))
				cmd = 'geeglm(' + self.formula.replace('___snv___',vu) + ',id=' + self.fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + self.family + ',corstr="' + self.corstr + '")'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					print vu + ": " + rerr.message
				else:
					ro.r('result<-summary(result)')
					vu = self.focus.replace('___snv___',vu)
					if ro.r('result$error')[0] == 0:
						ro.r('err<-0')
						ro.r('err[! is.finite(result$coefficients["' + vu + '",1]) | ! is.finite(result$coefficients["' + vu + '",2]) | ! is.finite(result$coefficients["' + vu + '",4]) | result$coefficients["' + vu + '",4] == 0]<-2')
						self.results['err'][v] = np.array(ro.r('err'))[:,None]
						self.results['effect'][v] = np.array(ro.r('result$coefficients["' + vu + '",1]'))[:,None]
						self.results['stderr'][v] = np.array(ro.r('result$coefficients["' + vu + '",2]'))[:,None]
						self.results['or'][v] = np.array(np.exp(ro.r('result$coefficients["' + vu + '",1]')))[:,None]
						self.results['wald'][v] = np.array(ro.r('result$coefficients["' + vu + '",3]'))[:,None]
						self.results['p'][v] = np.array(ro.r('result$coefficients["' + vu + '",4]'))[:,None]
					else:
						self.results['err'][v] = 1
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Glm(SnvModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Glm.__cinit__")
		logger.debug("initialize Glm model")
		super(Glm, self).__init__(**kwargs)
		print "setting glm test family option to " + self.family

		if self.reverse:
			if self.interact is not None:
				self.formula = '___snv___~' + self.dep_var + '*' + self.interact
				self.focus = self.dep_var + ':' + self.interact
			else:
				self.formula = '___snv___~' + self.dep_var
				self.focus = self.dep_var
			self.family = 'gaussian'
		else:
			if self.interact is not None:
				self.formula = self.dep_var + '~___snv___*' + self.interact
				self.focus = '___snv___:' + self.interact
			else:
				self.formula = self.dep_var + '~___snv___'
				self.focus = '___snv___'
		self.formula = self.formula + '+' + self.covars if self.covars is not None else self.formula
		print "formula: " + self.formula

		self.results_dtypes = [('err','f8'),('effect','f8'),('stderr','f8'),('or','f8'),('z','f8'),('p','f8')]
		if self.family == 'binomial':
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','or','z','p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','z','p']))

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + '## formula: ' + self.formula
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## err: error code (0: no error, 1: missing values returned by glm, 2: infinite value or zero p-value detected, 3: glm failed)' + '\n' + \
						'## effect: effect size' + '\n' + \
						'## stderr: standard error'
		self.metadata = self.metadata + '\n' + '## or: odds ratio (included only if binomial family)' if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## z: z-statistic' + '\n' + \
						'## p: p-value' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=self.results_dtypes)
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['id_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [vu]))
				cmd = 'glm(' + self.formula.replace('___snv___',vu) + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + self.family + ')'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					print vu + ": " + rerr.message
				else:
					vu = self.focus.replace('___snv___',vu)
					ro.r('result<-summary(result)')
					ro.r('err<-0')
					ro.r('err[result$coefficients["' + vu + '",1] == "NA"]<-1')
					ro.r('err[! is.finite(result$coefficients["' + vu + '",1]) | ! is.finite(result$coefficients["' + vu + '",2]) | ! is.finite(result$coefficients["' + vu + '",4]) | result$coefficients["' + vu + '",4] == 0]<-2')
					self.results['err'][v] = np.array(ro.r('err'))[:,None]
					self.results['effect'][v] = np.array(ro.r('result$coefficients["' + vu + '",1]'))[:,None]
					self.results['stderr'][v] = np.array(ro.r('result$coefficients["' + vu + '",2]'))[:,None]
					self.results['or'][v] = np.array(np.exp(ro.r('result$coefficients["' + vu + '",1]')))[:,None]
					self.results['z'][v] = np.array(ro.r('result$coefficients["' + vu + '",3]'))[:,None]
					self.results['p'][v] = np.array(ro.r('result$coefficients["' + vu + '",4]'))[:,None]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Lm(SnvModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Lm.__cinit__")
		logger.debug("initialize Lm model")
		super(Lm, self).__init__(**kwargs)

		if self.reverse:
			if self.interact is not None:
				self.formula = '___snv___~' + self.dep_var + '*' + self.interact
				self.focus = self.dep_var + ':' + self.interact
			else:
				self.formula = '___snv___~' + self.dep_var
				self.focus = self.dep_var
			self.family = 'gaussian'
		else:
			if self.interact is not None:
				self.formula = self.dep_var + '~___snv___*' + self.interact
				self.focus = '___snv___:' + self.interact
			else:
				self.formula = self.dep_var + '~___snv___'
				self.focus = '___snv___'
		self.formula = self.formula + '+' + self.covars if self.covars is not None else self.formula
		print "formula: " + self.formula

		self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','t','p']))
		self.results_dtypes = [('err','f8'),('effect','f8'),('stderr','f8'),('t','f8'),('p','f8')]

		self.metadata = self.metadata + '\n' + \
						'## formula: ' + self.formula + '\n' + \
						self.metadata_snv + '\n' + \
						'## err: error code (0: no error, 1: missing values returned by lm, 2: infinite value or zero p-value detected, 3: lm failed)' + '\n' + \
						'## effect: effect size' + '\n' + \
						'## stderr: standard error' + '\n' + \
						'## t: t-statistic' + '\n' + \
						'## p: p-value' + '\n#'
						
	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=self.results_dtypes)
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['id_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [vu]))
				cmd = 'lm(' + self.formula.replace('___snv___',vu) + ',data=na.omit(model_df[,names(model_df) %in% cols]))'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					print vu + ": " + rerr.message
				else:
					vu = self.focus.replace('___snv___',vu)
					ro.r('result<-summary(result)')
					ro.r('err<-0')
					ro.r('err[result$coefficients["' + vu + '",1] == "NA"]<-1')
					ro.r('err[! is.finite(result$coefficients["' + vu + '",1]) | ! is.finite(result$coefficients["' + vu + '",2]) | ! is.finite(result$coefficients["' + vu + '",4]) | result$coefficients["' + vu + '",4] == 0]<-2')
					self.results['err'][v] = np.array(ro.r('err'))[:,None]
					self.results['effect'][v] = np.array(ro.r('result$coefficients["' + vu + '",1]'))[:,None]
					self.results['stderr'][v] = np.array(ro.r('result$coefficients["' + vu + '",2]'))[:,None]
					self.results['t'][v] = np.array(ro.r('result$coefficients["' + vu + '",3]'))[:,None]
					self.results['p'][v] = np.array(ro.r('result$coefficients["' + vu + '",4]'))[:,None]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Lmer(SnvModel):
	def __cinit__(self, corstr = None, **kwargs):
		logger = logging.getLogger("Model.Lmer.__cinit__")
		logger.debug("initialize Lmer model")
		super(Lmer, self).__init__(**kwargs)
		print "lmer test will treat outcome as quantitative"

		if self.reverse:
			if self.interact is not None:
				self.formula = '___snv___~' + self.dep_var + '*' + self.interact
				self.focus = self.dep_var + ':' + self.interact
			else:
				self.formula = '___snv___~' + self.dep_var
				self.focus = self.dep_var
		else:
			if self.interact is not None:
				self.formula = self.dep_var + '~___snv___*' + self.interact
				self.focus = '___snv___:' + self.interact
			else:
				self.formula = self.dep_var + '~___snv___'
				self.focus = '___snv___'
		self.formula = self.formula + '+' + self.covars if self.covars is not None else self.formula

		if self.random_effects is not None:
			for reff in self.random_effects.split("+"):
				print "adding random effect " + reff
				self.formula = self.formula + '+(1|' + reff + ')' 
		print "formula: " + self.formula

		self.results_dtypes = [('effect','f8'),('stderr','f8'),('df','f8'),('t','f8'),('p','f8')]
		self.results_header = np.append(self.results_header,np.array(['effect','stderr','df','t','p']))
		if self.kr:
			self.results_dtypes = self.results_dtypes + [('df_kr','f8'),('f_kr','f8'),('p_kr','f8')]
			self.results_header = np.append(self.results_header,np.array(['df_kr','f_kr','p_kr']))

		print "loading R package lmerTest and pbkrtest and their dependencies"
		ro.r('suppressMessages(library(lmerTest))')
		ro.r('suppressMessages(library(pbkrtest))')

		print "setting R contrasts to Type III sum of squares"
		ro.r('options(contrasts = c("contr.sum","contr.poly"))')

		self.metadata = self.metadata + '\n' + '## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## effect: effect size' + '\n' + \
						'## stderr: standard error' + '\n' + \
						'## df: satterthwaite estimated degrees of freedom' + '\n' + \
						'## t: t statistic' + '\n' + \
						'## p: p value'
		if self.kr:
			self.metadata = self.metadata + '\n' + \
						'## df_kr: kenward-roger estimated degrees of freedom' + '\n' + \
						'## f_kr: kenward-roger f statistic' + '\n' + \
						'## p_kr: kenward-roger p-value'
		self.metadata = self.metadata + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=self.results_dtypes)
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			ro.r('model_df$' + self.fid + '<-as.factor(model_df$' + self.fid + ')')
			ro.r('model_df<-model_df[order(model_df$' + self.fid + '),]')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['id_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [self.fid] + [vu]))
				cmd = 'lmer(' + self.formula.replace('___snv___',vu) + ',data=na.omit(model_df[,names(model_df) %in% cols]),REML=' + str(self.reml).upper() + ')'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					print vu + ": " + rerr.message
				else:
					ro.r('result_base<-summary(result)')
					vu = self.focus.replace('___snv___',vu)
					self.results['effect'][v] = np.array(ro.r('result_base$coefficients["' + vu + '",1]'))[:,None]
					self.results['stderr'][v] = np.array(ro.r('result_base$coefficients["' + vu + '",2]'))[:,None]
					self.results['df'][v] = np.array(ro.r('result_base$coefficients["' + vu + '",3]'))[:,None]
					self.results['t'][v] = np.array(ro.r('result_base$coefficients["' + vu + '",4]'))[:,None]
					self.results['p'][v] = np.array(ro.r('result_base$coefficients["' + vu + '",5]'))[:,None]

					if self.kr:
						ro.r('result_kr<-anova(result, ddf="Kenward-Roger")')
						self.results['df_kr'][v] = np.array(ro.r('result_kr["' + vu + '","DenDF"]'))[:,None]
						self.results['f_kr'][v] = np.array(ro.r('result_kr["' + vu + '","F.value"]'))[:,None]
						self.results['p_kr'][v] = np.array(ro.r('result_kr["' + vu + '","Pr(>F)"]'))[:,None]

		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Skat(SnvgroupModel):
	cdef public str skat_wts, skat_method, mafrange
	cdef unsigned int timeout
	cdef public bint adjust_kinship
	def __cinit__(self, skat_wts = None, skat_method = None, mafrange = None, timeout = None, adjust_kinship = False, **kwargs):
		logger = logging.getLogger("Model.Skat.__cinit__")
		logger.debug("initialize Skat model")
		self.skat_wts = skat_wts if skat_wts is not None else 'function(maf){dbeta(maf,1,25)}'
		self.skat_method = skat_method if skat_method is not None else 'saddlepoint'
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.timeout = timeout if timeout is not None else 3600
		self.adjust_kinship = adjust_kinship
		super(Skat, self).__init__(**kwargs)
		print "setting skat test family option to " + self.family

		if self.covars is not None:
			self.formula = self.dep_var + '~' + self.covars
		else:
			self.formula = self.dep_var + '~1'
		print "formula: " + self.formula

		self.results_header = np.append(self.results_header,np.array(['total','passed','cmac','err','nmiss','nsnps','cmaf','p','q']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')
		if self.adjust_kinship:
			print "loading R package kinship2"
			ro.r('suppressMessages(library(kinship2))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		
		self.metadata = self.metadata + '\n' + \
						'## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + '## adjust kinship: True' if self.adjust_kinship else self.metadata + '\n' + '## adjust kinship: False'
		self.metadata = self.metadata + '\n' + '## skat wts: ' + self.skat_wts + '\n' + \
						'## skat method: ' + self.skat_method + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## total: group total snv count' + '\n' + \
						'## passed: group snv count that passed filters' + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatMeta timed out' + '\n' + \
						'## nmiss: number of missing genotypes in group' + '\n' + \
						'## nsnps: number of snps in the group' + '\n' + \
						'## cmaf: cumulative minor allele frequency' + '\n' + \
						'## p: group p-value' + '\n' + \
						'## q: skat q-statistic' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('total','uint32'),('passed','uint32'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnps','f8'),('cmaf','f8'),('p','f8'),('q','f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.group_id
		self.results['total'][0] = self.variants.info.shape[0]
		self.results['passed'][0] = len(passed)
		self.results['nmiss'][0] = np.nan
		self.results['nsnps'][0] = np.nan
		self.results['cmaf'][0] = np.nan
		self.results['p'][0] = np.nan
		self.results['q'][0] = np.nan
		if len(passed) > 1:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['id_unique'][passed]))
			ro.globalenv['snp_info'] = pd.DataFrame({'Name': list(self.variants.info['id_unique'][passed]), 'gene': 'NA'})
			ro.r('snp_info$Name<-as.character(snp_info$Name)')
			ro.r('snp_info$gene<-as.character(snp_info$gene)')
			ro.r('z<-data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['cmac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['cmac'][0] >= self.cmac:
				if len(passed) == 1:
					ro.r('colnames(z)<-"' + self.variants.info['id_unique'][passed][0] + '"')
				if self.adjust_kinship:
					ro.globalenv['ped'] = self.pedigree
					ro.globalenv['kins'] = ro.r("kinship(pedigree(famid=ped$" + self.fid + ",id=ped$" + self.iid + ",dadid=ped$" + self.patid + ",momid=ped$" + self.matid + ",sex=ped$" + self.sex + ",missid='NA'))")
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "',kins=kins,sparse=FALSE), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				else:
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					print rerr.message
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = "tryCatch(expr = { evalWithTimeout(skatMeta(ps,SNPInfo=snp_info,wts=" + self.skat_wts + ",method='" + self.skat_method + "',mafRange=" + self.mafrange + "), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						print rerr.message
					else:
						if ro.r('class(result) == "data.frame"')[0]:
							ro.r('result$err<-0')
							ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
							ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$nsnps[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$q[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$cmaf[! is.na(result$err) & result$err == 1]<-NA')
							self.results['err'][0] = np.array(ro.r('result$err'))[:,None]
							self.results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
							self.results['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
							self.results['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
							self.results['p'][0] = np.array(ro.r('result$p'))[:,None]
							self.results['q'][0] = np.array(ro.r('result$Q'))[:,None]
						else:
							self.results['err'][0] = 8
				else:
					self.results['err'][0] = 7
			else:
				self.results['err'][0] = 6
		else:
			self.results['err'][0] = 5
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['total','passed','cmac','err','nmiss','nsnps','cmaf','p','q']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + 'total','uint32'),(tag + 'passed','uint32'),(tag + '.cmac','f8'),(tag + '.err','f8'),(tag + '.nmiss','f8'),(tag + '.nsnps','f8'),(tag + '.cmaf','f8'),(tag + '.p','f8'),(tag + '.q','f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')
			ro.globalenv[tag + '_cmac'] = self.out[tag + '.cmac'][0]

cdef class Skato(SnvgroupModel):
	cdef public str skat_wts, burden_wts, skat_method, mafrange, skato_rho
	cdef unsigned int timeout
	cdef public bint adjust_kinship
	def __cinit__(self, skat_wts = None, burden_wts = None, skat_method = None, skato_rho = None, mafrange = None, timeout = None, adjust_kinship = False, **kwargs):
		logger = logging.getLogger("Model.Skato.__cinit__")
		logger.debug("initialize Skato model")
		self.skat_wts = skat_wts if skat_wts is not None else 'function(maf){dbeta(maf,1,25)}'
		self.burden_wts = burden_wts if burden_wts is not None else 'function(maf){maf < 0.01}'
		self.skat_method = skat_method if skat_method is not None else 'saddlepoint'
		self.skato_rho = skato_rho if skato_rho is not None else 'seq(0,1,0.1)'
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.timeout = timeout if timeout is not None else 3600
		self.adjust_kinship = adjust_kinship
		super(Skato, self).__init__(**kwargs)
		print "setting skat-o test family option to " + self.family

		if self.covars is not None:
			self.formula = self.dep_var + '~' + self.covars
		else:
			self.formula = self.dep_var + '~1'
		print "formula: " + self.formula

		self.results_header = np.append(self.results_header,np.array(['total','passed','cmac','err','nmiss','nsnps','cmaf','p','pmin','rho']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')
		if self.adjust_kinship:
			print "loading R package kinship2"
			ro.r('suppressMessages(library(kinship2))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + '## adjust kinship: True' if self.adjust_kinship else self.metadata + '\n' + '## adjust kinship: False'
		self.metadata = self.metadata + '\n' + '## skat wts: ' + self.skat_wts + '\n' + \
						'## burden wts: ' + self.burden_wts + '\n' + \
						'## skat method: ' + self.skat_method + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## total: group total snv count' + '\n' + \
						'## passed: group snv count that passed filters' + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatOMeta, 4: skatOMeta errflag > 0, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatOMeta timed out' + '\n' + \
						'## nmiss: number of missing genotypes in group' + '\n' + \
						'## nsnps: number of snps in the group' + '\n' + \
						'## cmaf: cumulative minor allele frequency' + '\n' + \
						'## p: group p-value' + '\n' + \
						'## pmin: minimum snp p-value' + '\n' + \
						'## rho: skato rho parameter' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('total','uint32'),('passed','uint32'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnps','f8'),('cmaf','f8'),('p','f8'),('pmin','f8'),('rho','f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.group_id
		self.results['total'][0] = self.variants.info.shape[0]
		self.results['passed'][0] = len(passed)
		self.results['nmiss'][0] = np.nan
		self.results['nsnps'][0] = np.nan
		self.results['cmaf'][0] = np.nan
		self.results['pmin'][0] = np.nan
		self.results['p'][0] = np.nan
		self.results['rho'][0] = np.nan
		if len(passed) > 1:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['id_unique'][passed]))
			ro.globalenv['snp_info'] = pd.DataFrame({'Name': list(self.variants.info['id_unique'][passed]), 'gene': 'NA'})
			ro.r('snp_info$Name<-as.character(snp_info$Name)')
			ro.r('snp_info$gene<-as.character(snp_info$gene)')
			ro.r('z<-data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['cmac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['cmac'][0] >= self.cmac:
				if passed == 1:
					ro.r('colnames(z)<-"' + self.variants.info['id_unique'][passed][0] + '"')
				if self.adjust_kinship:
					ro.globalenv['ped'] = self.pedigree
					ro.globalenv['kins'] = ro.r("kinship(pedigree(famid=ped$" + self.fid + ",id=ped$" + self.iid + ",dadid=ped$" + self.patid + ",momid=ped$" + self.matid + ",sex=ped$" + self.sex + ",missid='NA'))")
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "',kins=kins,sparse=FALSE), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				else:
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					print rerr.message
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = "tryCatch(expr = { evalWithTimeout(skatOMeta(ps,SNPInfo=snp_info,rho=" + self.skato_rho + ",skat.wts=" + self.skat_wts + ",burden.wts=" + self.burden_wts + ",method='" + self.skat_method + "',mafRange=" + self.mafrange + "), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						print rerr.message
					else:
						if ro.r('class(result) == "data.frame"')[0]:
							ro.r('result$err<-0')
							ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
							ro.r('result$err[! is.finite(result$errflag) | result$errflag > 0]<-4')
							ro.r('result$nmiss[! is.na(result$err) & result$err > 0]<-NA')
							ro.r('result$nsnps[! is.na(result$err) & result$err > 0]<-NA')
							ro.r('result$p[! is.na(result$err) & result$err > 0]<-NA')
							ro.r('result$pmin[! is.na(result$err) & result$err > 0]<-NA')
							ro.r('result$rho[! is.na(result$err) & result$err > 0]<-NA')
							ro.r('result$cmaf[! is.na(result$err) & result$err > 0]<-NA')
							self.results['err'][0] = np.array(ro.r('result$err'))[:,None]
							self.results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
							self.results['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
							self.results['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
							self.results['pmin'][0] = np.array(ro.r('result$pmin'))[:,None]
							self.results['p'][0] = np.array(ro.r('result$p'))[:,None]
							self.results['rho'][0] = np.array(ro.r('result$rho'))[:,None]
						else:
							self.results['err'][0] = 8
				else:
					self.results['err'][0] = 7
			else:
				self.results['err'][0] = 6
		else:
			self.results['err'][0] = 5
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['total','passed','cmac','err','nmiss','nsnps','cmaf','p','pmin','rho']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + '.total','uint32'),(tag + '.passed','uint32'),(tag + '.cmac','f8'),(tag + '.err','f8'),(tag + '.nmiss','f8'),(tag + '.nsnps','f8'),(tag + '.cmaf','f8'),(tag + '.p','f8'),(tag + '.pmin','f8'),(tag + '.rho','f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')
			ro.globalenv[tag + '_cmac'] = self.out[tag + '.cmac'][0]

cdef class Burden(SnvgroupModel):
	cdef public str mafrange, burden_wts
	cdef unsigned int timeout
	cdef public bint adjust_kinship
	def __cinit__(self, mafrange = None, burden_wts = None, timeout = None, adjust_kinship = False, **kwargs):
		logger = logging.getLogger("Model.Burden.__cinit__")
		logger.debug("initialize Burden model")
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.burden_wts = burden_wts if burden_wts is not None else '1'
		self.timeout = timeout if timeout is not None else 3600
		self.adjust_kinship = adjust_kinship
		super(Burden, self).__init__(**kwargs)
		print "setting burden test family option to " + self.family

		if self.covars is not None:
			self.formula = self.dep_var + '~' + self.covars
		else:
			self.formula = self.dep_var + '~1'
		print "formula: " + self.formula

		self.results_header = np.append(self.results_header,np.array(['total','passed','cmac','err','nmiss','nsnpsTotal','nsnpsUsed','cmafTotal','cmafUsed','beta','se','p']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')
		if self.adjust_kinship:
			print "loading R package kinship2"
			ro.r('suppressMessages(library(kinship2))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + '## adjust kinship: True' if self.adjust_kinship else self.metadata + '\n' + '## adjust kinship: False'
		self.metadata = self.metadata + '\n' + '## burden wts: ' + self.burden_wts + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## total: group total snv count' + '\n' + \
						'## passed: group snv count that passed filters' + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed burdenMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: burdenMeta timed out' + '\n' + \
						'## nmiss: number of missing snps' + '\n' + \
						'## nsnpsTotal: number of snps in group' + '\n' + \
						'## nsnpsUsed: number of snps used in analysis' + '\n' + \
						'## cmafTotal: cumulative minor allele frequency in group' + '\n' + \
						'## cmafUsed: cumulative minor allele frequency for snps used in analysis' + '\n' + \
						'## beta: coefficient for effect of genotype' + '\n' + \
						'## se: standard error for effect of genotype' + '\n' + \
						'## p: p value for burden test' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('total','uint32'),('passed','uint32'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnpsTotal','f8'),('nsnpsUsed','f8'),('cmafTotal','f8'),('cmafUsed','f8'),('beta','f8'),('se','f8'),('p','f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.group_id
		self.results['total'][0] = self.variants.info.shape[0]
		self.results['passed'][0] = len(passed)
		self.results['nmiss'][0] = np.nan
		self.results['nsnpsTotal'][0] = np.nan
		self.results['nsnpsUsed'][0] = np.nan
		self.results['cmafTotal'][0] = np.nan
		self.results['cmafUsed'][0] = np.nan
		self.results['beta'][0] = np.nan
		self.results['se'][0] = np.nan
		self.results['p'][0] = np.nan
		if len(passed) > 1:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['id_unique'][passed])
			ro.globalenv['model_df'] = pheno_df.merge(variants_df, on=self.iid, how='left')
			for col in list(self.model_cols) + list(self.variants.info['id_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['id_unique'][passed]))
			ro.globalenv['snp_info'] = pd.DataFrame({'Name': list(self.variants.info['id_unique'][passed]), 'gene': 'NA'})
			ro.r('snp_info$Name<-as.character(snp_info$Name)')
			ro.r('snp_info$gene<-as.character(snp_info$gene)')
			ro.r('z<-data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['cmac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['cmac'][0] >= self.cmac:
				if len(passed) == 1:
					ro.r('colnames(z)<-"' + self.variants.info['id_unique'][passed][0] + '"')
				if self.adjust_kinship:
					ro.globalenv['ped'] = self.pedigree
					ro.globalenv['kins'] = ro.r("kinship(pedigree(famid=ped$" + self.fid + ",id=ped$" + self.iid + ",dadid=ped$" + self.patid + ",momid=ped$" + self.matid + ",sex=ped$" + self.sex + ",missid='NA'))")
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "',kins=kins,sparse=FALSE), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				else:
					cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					print rerr.message
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = 'tryCatch(expr = { evalWithTimeout(burdenMeta(ps,SNPInfo=snp_info,mafRange=' + self.mafrange + ',wts=' + self.burden_wts + '), timeout=' + str(self.timeout) + ') }, error = function(e) return(1))'
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						print rerr.message
					else:
						if ro.r('class(result) == "data.frame"')[0]:
							ro.r('result$err<-0')
							ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
							ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$nsnpsTotal[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$nsnpsUsed[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$cmafTotal[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$cmafUsed[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$beta[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$se[! is.na(result$err) & result$err == 1]<-NA')
							ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
							self.results['err'][0] = np.array(ro.r('result$err'))[:,None]
							self.results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
							self.results['nsnpsTotal'][0] = np.array(ro.r('result$nsnpsTotal'))[:,None]
							self.results['nsnpsUsed'][0] = np.array(ro.r('result$nsnpsUsed'))[:,None]
							self.results['cmafTotal'][0] = np.array(ro.r('result$cmafTotal'))[:,None]
							self.results['cmafUsed'][0] = np.array(ro.r('result$cmafUsed'))[:,None]
							self.results['beta'][0] = np.array(ro.r('result$beta'))[:,None]
							self.results['se'][0] = np.array(ro.r('result$se'))[:,None]
							self.results['p'][0] = np.array(ro.r('result$p'))[:,None]
						else:
							self.results['err'][0] = 8
				else:
					self.results['err'][0] = 7
			else:
				self.results['err'][0] = 6
		else:
			self.results['err'][0] = 5
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['total','passed','cmac','err','nmiss','nsnpsTotal','nsnpsUsed','cmafTotal','cmafUsed','beta','se','p']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + '.total','uint32'),(tag + '.passed','uint32'),(tag + '.cmac','f8'),(tag + '.err','f8'),(tag + '.nmiss','f8'),(tag + '.nsnpsTotal','f8'),(tag + '.nsnpsUsed','f8'),(tag + '.cmafTotal','f8'),(tag + '.cmafUsed','f8'),(tag + '.beta','f8'),(tag + '.se','f8'),(tag + '.p','f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')
			ro.globalenv[tag + '_cmac'] = self.out[tag + '.cmac'][0]

cdef class Neff(SnvgroupModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Neff.__cinit__")
		logger.debug("initialize Neff model")
		super(Neff, self).__init__(**kwargs)
		print "setting neff test family option to " + self.family

		if self.covars is not None:
			self.formula = self.dep_var + '~' + self.covars
		else:
			self.formula = self.dep_var + '~1'
		print "formula: " + self.formula

		self.results_header = np.append(self.results_header,np.array(['total','passed','eff']))

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## formula: ' + self.formula
		self.metadata = self.metadata + '\n' + \
						self.metadata_gene + '\n' + \
						'## total: group total snv count' + '\n' + \
						'## passed: group snv count that passed filters' + '\n' + \
						'## eff: number of effective tests in group' + '\n#'

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno_df,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('total','uint32'),('passed','uint32'),('eff','f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.group_id
		self.results['total'][0] = self.variants.info.shape[0]
		self.results['passed'][0] = len(passed)
		self.results['eff'][0] = np.nan
		if len(passed) > 0:
			if len(passed) > 1:
				variants_df = pd.DataFrame(self.variants.data[:,passed_data],dtype='float64')
				variants_df.dropna(inplace=True)
				markers_cor = variants_df.corr()
				markers_cor_eigvals = np.linalg.eigvalsh(markers_cor)
				markers_cor_eigvalsnew = [x if x > 0 else 0 for x in markers_cor_eigvals]
				self.results['eff'][0] = sum([1 if x else 0 for x in markers_cor_eigvals>1]+(markers_cor_eigvalsnew-np.floor(markers_cor_eigvalsnew)))
			else:
				self.results['eff'][0] = 1.0
		else:
			self.results['eff'][0] = 0.0
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')

cdef class Meta(object):
	cdef public unsigned int tbx_start, tbx_end
	cdef public str tag, meta, metadata
	cdef public np.ndarray results_header
	cdef public object out
	def __cinit__(self, tag, meta, **kwargs):
		super(Meta, self).__init__(**kwargs)
		logger = logging.getLogger("Model.Meta.__cinit__")
		logger.debug("initialize meta")
		self.tag = tag
		self.meta = meta
		self.metadata = '## source: uga' + version + '\n' + \
						'## date: ' + time.strftime("%Y%m%d") + '\n' + \
						'## method: meta\n' + \
						'## formula: ' + self.meta

cdef class SnvMeta(Meta):
	cdef public str type
	cdef public object df
	cdef public np.ndarray meta_incl
	def __cinit__(self, type = "sample_size", **kwargs):
		logger = logging.getLogger("Model.SnvMeta.__cinit__")
		logger.debug("initialize SnvMeta")
		self.type = type
		super(SnvMeta, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 1
		self.meta_incl = np.array(self.meta.split('+'))

		if self.type == 'stderr':
			self.results_header = np.array(['chr','pos','id','a1','a2','effect','stderr','or','z','p','dir','n','hetq','hetdf','heti2','hetp'])
		else:
			self.results_header = np.array(['chr','pos','id','a1','a2','z','p','dir','n'])

		self.metadata = self.metadata + '\n' + \
						'## chr: chromosome' + '\n' + \
						'## pos: chromosomal position' + '\n' + \
						'## id: snv name' + '\n' + \
						'## a1: reference (coded) allele used for stat calculations' + '\n' + \
						'## a2: alternate (non-coded) allele' + '\n'
		self.metadata = self.metadata + '## effect: effect size' + '\n' + '## stderr: standard error' + '\n' + '## or: odds ratio' + '\n' if self.type == "stderr" else self.metadata
		self.metadata = self.metadata + \
						'## z: z score' + '\n' + \
						'## p: p-value' + '\n' + \
						"## dir: direction of effect of each cohort (same order as formula; 'x' indicates cohort excluded from meta)" + '\n' + \
						'## n: sample size'
		self.metadata = self.metadata + '\n' + '## hetq: heterogeneity chi-squared statistic' + '\n' + '## hetdf: heterogeneity degrees of freedom' + '\n' + '## heti2: heterogeneity inconsistency statistic (I-squared)' + '\n' + '## hetp: heterogeneity p-value\n#' if self.type == "stderr" else self.metadata + '\n#'

	def calc_meta(self, df):
		# df: a dataframe with tagged input data using meta_incl tags
		header = list(df.columns)
		if self.type == 'sample_size':
			df['meta.dir'] = ''
			for tag in self.meta_incl:
				filter_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.filter')][0]
				N_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.p')][0]
				Eff_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				df['meta.' + tag + '.dir'] = df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				df['meta.' + tag + '.zi'] = df.apply(lambda x: (-1 * scipy.norm.ppf(1 - (x[P_idx]/2)) if x[Eff_idx] < 0 else scipy.norm.ppf(1 - (x[P_idx]/2))) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.' + tag + '.n'] = df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.' + tag + '.wi'] = df.apply(lambda x: math.sqrt(x[tag + '.n']) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.' + tag + '.ziwi'] = df.apply(lambda x: x['meta.' + tag + '.zi'] * x['meta.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.dir'] = df['meta.dir'] + df['meta.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith('meta') and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith('meta') and s.endswith('.wi')]
			ZiWi_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith('meta') and s.endswith('.ziwi')]
			df['meta.n'] = df.apply(lambda x: x[N_idx_all].sum() if len(x['meta.dir'].replace('x','')) > 0 else float('nan'),axis=1)
			df['meta.z'] = df.apply(lambda x: x[ZiWi_idx_all].sum()/math.sqrt(x[N_idx_all].sum()) if len(x['meta.dir'].replace('x','')) > 0 else float('nan'), axis=1)
			df['meta.stderr'] = df.apply(lambda x: float('nan'), axis=1)
			df['meta.effect'] = df.apply(lambda x: float('nan'), axis=1)
			df['meta.or'] = df.apply(lambda x: float('nan'), axis=1)
		else:
			df['meta.dir'] = ''
			for tag in self.meta_incl:
				filter_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.filter')][0]
				N_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.p')][0]
				Eff_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				StdErr_idx=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tag) and s.endswith('.stderr')][0]
				df['meta.' + tag + '.dir'] = df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				df['meta.' + tag + '.n'] = df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.' + tag + '.wi'] = df.apply(lambda x: 1/(x[StdErr_idx]**2) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.' + tag + '.biwi'] = df.apply(lambda x: x[Eff_idx] * x['meta.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				df['meta.dir'] = df['meta.dir'] + df['meta.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tuple('meta.' + x for x in self.meta_incl)) and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tuple('meta.' + x for x in self.meta_incl)) and s.endswith('.wi')]
			BiWi_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tuple('meta.' + x for x in self.meta_incl)) and s.endswith('.biwi')]
			df['meta.n'] = df.apply(lambda x: x[N_idx_all].sum() if len(x['meta.dir'].replace('x','')) > 0 else x[N_idx] if len(x['meta.dir'].replace('x','')) == 1 else  float('nan'),axis=1)
			df['meta.stderr'] = df.apply(lambda x: math.sqrt(1/(x[Wi_idx_all].sum())) if len(x['meta.dir'].replace('x','')) > 0 else float('nan'), axis=1)
			df['meta.effect'] = df.apply(lambda x: (x[BiWi_idx_all].sum())/(x[Wi_idx_all].sum()) if len(x['meta.dir'].replace('x','')) > 0 else float('nan'), axis=1)
			df['meta.or'] = df.apply(lambda x: math.exp(x['meta.effect']) if len(x['meta.dir'].replace('x','')) > 0 and len([k + '.or' in df for k in self.meta_incl]) > 0 and not x['meta.effect'] > 709.782712893384 and not x['meta.effect'] < -709.782712893384 else float('nan'), axis=1)
			df['meta.z'] = df.apply(lambda x: x['meta.effect']/x['meta.stderr'] if len(x['meta.dir'].replace('x','')) > 0 else float('nan'), axis=1)
			for tag in self.meta_incl:
				df['meta.' + tag + '.hetq_pre'] = df.apply(lambda x: x['meta.' + tag + '.wi'] * (x[tag + '.effect'] - x['meta.effect'])**2 if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
			Hetq_pre_idx_all=[i for i, s in enumerate(list(df.columns.values)) if s.startswith(tuple('meta.' + x for x in self.meta_incl)) and s.endswith('.hetq_pre')]
			df['meta.hetq'] = df.apply(lambda x: x[Hetq_pre_idx_all].sum() if len(x['meta.dir'].replace('x','')) > 1 else float('nan'),axis=1)
			df['meta.hetdf'] = df.apply(lambda x: len([a for a in x['meta.dir'] if a != 'x'])-1 if len(x['meta.dir'].replace('x','')) > 1 else float('nan'),axis=1)
			df['meta.heti2'] = df.apply(lambda x: ((x['meta.hetq']-x['meta.hetdf'])/x['meta.hetq'])*100 if len(x['meta.dir'].replace('x','')) > 1 and x['meta.hetq'] != 0 else float('nan'),axis=1)
			df['meta.hetp'] = df.apply(lambda x: 1-scipy.chi2.cdf(x['meta.hetq'], x['meta.hetdf']) if len(x['meta.dir'].replace('x','')) > 1 else float('nan'),axis=1)
		df['meta.p'] = df.apply(lambda x: 2 * scipy.norm.cdf(-1 * abs(float(x['meta.z']))) if len(x['meta.dir'].replace('x','')) > 0 else float('nan'), axis=1)
		df['meta.dir'] = df.apply(lambda x: x['meta.dir'] if not math.isnan(x['meta.p']) else float('nan'), axis=1)
		df = df[['chr','pos','id','a1','a2'] + [x for x in df.columns if 'meta.' in x]]
		df.columns = [x.replace('meta.','') for x in df.columns]
		self.out = df[self.results_header]

cdef class SkatMeta(Meta):
	cdef public object df
	cdef public np.ndarray meta_incl
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.SkatMeta.__cinit__")
		logger.debug("initialize SkatMeta")
		super(SkatMeta, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 2
		self.meta_incl = np.array(self.meta.split('+'))
		self.results_header = np.array(['chr','start','end','id','incl','cmac','err','nmiss','nsnps','cmaf','p','q'])
		self.metadata = self.metadata + '\n' + \
						'## chr: chromosome' + '\n' + \
						'## start: start chromosomal position' + '\n' + \
						'## end: end chromosomal position' + '\n' + \
						'## id: snv group name (ie. gene name)' + '\n' + \
						"## incl: cohort inclusion ('+' = included, 'x' = excluded)" + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatMeta timed out' + '\n' + \
						'## nmiss: number of missing genotypes in group' + '\n' + \
						'## nsnps: number of snps in the group' + '\n' + \
						'## cmaf: cumulative minor allele frequency' + '\n' + \
						'## p: group p-value' + '\n' + \
						'## q: skat q-statistic' + '\n#'

	def calc_meta(self, chr, start, end, id, Skat obj, incl):
		# obj: the first model object from individual analyses
		# incl: a list of individual model tags to include in meta
		results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('incl','|S100'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnps','f8'),('cmaf','f8'),('p','f8'),('q','f8')])
		results['chr'][0] = chr
		results['start'][0] = start
		results['end'][0] = end
		results['id'][0] = id
		for m in incl:
			if m == incl[0]:
				ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
			else:
				ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
		results['incl'][0] = ''.join(['+' if a in incl else 'x' for a in self.meta_incl])
		results['cmac'][0] = ro.r('sum(' + ",".join([x + "_cmac" for x in self.meta_incl if x in incl]) + ')')[0]
		if str(results['incl'][0]).count('+') > 1:
			cmd = "tryCatch(expr = { evalWithTimeout(skatMeta(" + ",".join([x + "_ps" for x in self.meta_incl if x in incl]) + ",SNPInfo=snp_info_meta,wts=" + obj.skat_wts + ",method='" + obj.skat_method + "',mafRange=" + obj.mafrange + "), timeout=" + str(obj.timeout) + ") }, error = function(e) return(1))"
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				results['err'][0] = 3
				print rerr.message
			else:
				if ro.r('class(result) == "data.frame"')[0]:
					ro.r('result$err<-0')
					ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
					ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$nsnps[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$q[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$cmaf[! is.na(result$err) & result$err == 1]<-NA')
					results['err'][0] = np.array(ro.r('result$err'))[:,None]
					results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
					results['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
					results['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
					results['p'][0] = np.array(ro.r('result$p'))[:,None]
					results['q'][0] = np.array(ro.r('result$Q'))[:,None]
				else:
					results['err'][0] = 8
		else:
			results['err'][0] = 5
			results['nmiss'][0] = np.nan
			results['nsnps'][0] = np.nan
			results['cmaf'][0] = np.nan
			results['p'][0] = np.nan
			results['q'][0] = np.nan
		self.out = pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')[self.results_header]

cdef class SkatoMeta(Meta):
	cdef public object df
	cdef public np.ndarray meta_incl
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.SkatoMeta.__cinit__")
		logger.debug("initialize SkatoMeta")
		super(SkatoMeta, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 2
		self.meta_incl = np.array(self.meta.split('+'))
		self.results_header = np.array(['chr','start','end','id','incl','cmac','err','nmiss','nsnps','cmaf','p','pmin','rho'])
		self.metadata = self.metadata + '\n' + \
						'## chr: chromosome' + '\n' + \
						'## start: start chromosomal position' + '\n' + \
						'## end: end chromosomal position' + '\n' + \
						'## id: snv group name (ie. gene name)' + '\n' + \
						"## incl: cohort inclusion ('+' = included, 'x' = excluded)" + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatMeta timed out' + '\n' + \
						'## nmiss: number of missing genotypes in group' + '\n' + \
						'## nsnps: number of snps in the group' + '\n' + \
						'## cmaf: cumulative minor allele frequency' + '\n' + \
						'## p: group p-value' + '\n' + \
						'## pmin: minimum snp p-value' + '\n' + \
						'## rho: skato rho parameter' + '\n#'

	def calc_meta(self, chr, start, end, id, Skato obj, incl):
		# obj: the first model object from individual analyses
		# incl: a list of individual model tags to include in meta
		results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('incl','|S100'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnps','f8'),('cmaf','f8'),('p','f8'),('pmin','f8'),('rho','f8')])
		results['chr'][0] = chr
		results['start'][0] = start
		results['end'][0] = end
		results['id'][0] = id
		for m in incl:
			if m == incl[0]:
				ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
			else:
				ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
		results['incl'][0] = ''.join(['+' if a in incl else 'x' for a in self.meta_incl])
		results['cmac'][0] = ro.r('sum(' + ",".join([x + "_cmac" for x in self.meta_incl if x in incl]) + ')')[0]
		if str(results['incl'][0]).count('+') > 1:
			cmd = "tryCatch(expr = { evalWithTimeout(skatOMeta(" + ",".join([x + "_ps" for x in self.meta_incl if x in incl]) + ",SNPInfo=snp_info_meta,rho=" + obj.skato_rho + ",skat.wts=" + obj.skat_wts + ",burden.wts=" + obj.burden_wts + ",method='" + obj.skat_method + "',mafRange=" + obj.mafrange + "), timeout=" + str(obj.timeout) + ") }, error = function(e) return(1))"
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				results['err'][0] = 3
				print rerr.message
			else:
				if ro.r('class(result) == "data.frame"')[0]:
					ro.r('result$err<-0')
					ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
					ro.r('result$err[! is.finite(result$errflag) | result$errflag > 0]<-4')
					ro.r('result$nmiss[! is.na(result$err) & result$err > 0]<-NA')
					ro.r('result$nsnps[! is.na(result$err) & result$err > 0]<-NA')
					ro.r('result$p[! is.na(result$err) & result$err > 0]<-NA')
					ro.r('result$pmin[! is.na(result$err) & result$err > 0]<-NA')
					ro.r('result$rho[! is.na(result$err) & result$err > 0]<-NA')
					ro.r('result$cmaf[! is.na(result$err) & result$err > 0]<-NA')
					results['err'][0] = np.array(ro.r('result$err'))[:,None]
					results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
					results['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
					results['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
					results['pmin'][0] = np.array(ro.r('result$pmin'))[:,None]
					results['p'][0] = np.array(ro.r('result$p'))[:,None]
					results['rho'][0] = np.array(ro.r('result$rho'))[:,None]
				else:
					results['err'][0] = 8
		else:
			results['err'][0] = 5
			results['nmiss'][0] = np.nan
			results['nsnps'][0] = np.nan
			results['cmaf'][0] = np.nan
			results['pmin'][0] = np.nan
			results['p'][0] = np.nan
			results['rho'][0] = np.nan
		self.out = pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')[self.results_header]

cdef class BurdenMeta(Meta):
	cdef public object df
	cdef public np.ndarray meta_incl
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.BurdenMeta.__cinit__")
		logger.debug("initialize BurdenMeta")
		super(BurdenMeta, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 2
		self.meta_incl = np.array(self.meta.split('+'))

		self.results_header = np.array(['chr','start','end','id','incl','cmac','err','nmiss','nsnpsTotal','nsnpsUsed','cmafTotal','cmafUsed','beta','se','p'])

		self.metadata = self.metadata + '\n' + \
						'## chr: chromosome' + '\n' + \
						'## start: start chromosomal position' + '\n' + \
						'## end: end chromosomal position' + '\n' + \
						'## id: snv group name (ie. gene name)' + '\n' + \
						"## incl: cohort inclusion ('+' = included, 'x' = excluded)" + '\n' + \
						'## cmac: group minor allele count' + '\n' + \
						'## err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed burdenMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: burdenMeta timed out' + '\n' + \
						'## nmiss: number of missing snps' + '\n' + \
						'## nsnpsTotal: number of snps in group' + '\n' + \
						'## nsnpsUsed: number of snps used in analysis' + '\n' + \
						'## cmafTotal: cumulative minor allele frequency in group' + '\n' + \
						'## cmafUsed: cumulative minor allele frequency for snps used in analysis' + '\n' + \
						'## beta: coefficient for effect of genotype' + '\n' + \
						'## se: standard error for effect of genotype' + '\n' + \
						'## p: p value for burden test' + '\n#'

	def calc_meta(self, chr, start, end, id, Burden obj, incl):
		# obj: the first model object from individual analyses
		# incl: a list of individual model tags to include in meta
		results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('incl','|S100'),('cmac','f8'),('err','f8'),('nmiss','f8'),('nsnpsTotal','f8'),('nsnpsUsed','f8'),('cmafTotal','f8'),('cmafUsed','f8'),('beta','f8'),('se','f8'),('p','f8')])
		results['chr'][0] = chr
		results['start'][0] = start
		results['end'][0] = end
		results['id'][0] = id
		for m in incl:
			if m == incl[0]:
				ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
			else:
				ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
		results['incl'][0] = ''.join(['+' if a in incl else 'x' for a in self.meta_incl])
		results['cmac'][0] = ro.r('sum(' + ",".join([x + "_cmac" for x in self.meta_incl if x in incl]) + ')')[0]
		if str(results['incl'][0]).count('+') > 1:
			cmd = 'tryCatch(expr = { evalWithTimeout(burdenMeta(' + ",".join([x + "_ps" for x in self.meta_incl if x in incl]) + ',SNPInfo=snp_info_meta,mafRange=' + obj.mafrange + ',wts=' + obj.burden_wts + '), timeout=' + str(obj.timeout) + ') }, error = function(e) return(1))'
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				results['err'][0] = 3
				print rerr.message
			else:
				if ro.r('class(result) == "data.frame"')[0]:
					ro.r('result$err<-0')
					ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
					ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$nsnpsTotal[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$nsnpsUsed[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$cmafTotal[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$cmafUsed[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$beta[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$se[! is.na(result$err) & result$err == 1]<-NA')
					ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
					results['err'][0] = np.array(ro.r('result$err'))[:,None]
					results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
					results['nsnpsTotal'][0] = np.array(ro.r('result$nsnpsTotal'))[:,None]
					results['nsnpsUsed'][0] = np.array(ro.r('result$nsnpsUsed'))[:,None]
					results['cmafTotal'][0] = np.array(ro.r('result$cmafTotal'))[:,None]
					results['cmafUsed'][0] = np.array(ro.r('result$cmafUsed'))[:,None]
					results['beta'][0] = np.array(ro.r('result$beta'))[:,None]
					results['se'][0] = np.array(ro.r('result$se'))[:,None]
					results['p'][0] = np.array(ro.r('result$p'))[:,None]
				else:
					results['err'][0] = 8
		else:
			results['err'][0] = 5
			results['nmiss'][0] = np.nan
			results['nsnpsTotal'][0] = np.nan
			results['nsnpsUsed'][0] = np.nan
			results['cmafTotal'][0] = np.nan
			results['cmafUsed'][0] = np.nan
			results['beta'][0] = np.nan
			results['se'][0] = np.nan
			results['p'][0] = np.nan
		self.out = pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')[self.results_header]
