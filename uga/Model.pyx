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
import Process
cimport numpy as np
cimport cython
import Geno
cimport Variant
import Variant
import Fxns
import rpy2.robjects as ro
import pandas.rpy.common as py2r
from rpy2.rinterface import RRuntimeError
from __version__ import version
import time
import numpy.lib.recfunctions as recfxns
import logging
ro.r('options(warn=1)')
ro.r('options(na.action=na.omit)')
pd.options.mode.chained_assignment = None
ro.r('suppressMessages(library(R.utils))')

module_logger = logging.getLogger("Model")

cdef class Model(object):
	cdef public unsigned int case_code, ctrl_code, tbx_start, tbx_end, \
								male, female, nobs, nunique, nfounders, nlongitudinal, \
								nfamilies, nunrelated, ncases, nctrls, nmales, nfemales
	cdef public np.ndarray cases_idx, ctrls_idx, male_cases_idx, male_ctrls_idx, female_cases_idx, female_ctrls_idx, \
							focus, model_cols, results_header, results_dtypes, calc_hwe_idx, \
							variant_stats, results, unique_idx, founders_idx, founders_ctrls_idx, male_idx, female_idx
	cdef public bytes fxn, formula, format, pheno_file, variants_file, type, samples_file, \
						iid, fid, matid, patid, sex, pheno_sep, a1, a2
	cdef public str metadata, metadata_cc, dep_var, family
	cdef public dict fields
	cdef public object pheno, variants, out
	cdef public bint all_founders
	def __cinit__(self, fxn, formula, format, variants_file, pheno_file, type, iid, fid, 
					case_code = None, ctrl_code = None, all_founders = False, samples_file = None, 
					matid = None, patid = None, sex = None, male = 1, female = 2, pheno_sep = 'tab', **kwargs):
		super(Model, self).__init__(**kwargs)
		logger = logging.getLogger("Model.Model.__cinit__")
		logger.debug("initialize model")
		self.out = None
		self.fxn = fxn
		self.formula = formula
		self.format = format
		self.variants_file = variants_file
		self.samples_file = samples_file
		self.pheno_file = pheno_file
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
		self.pheno_sep = pheno_sep
		self.all_founders = all_founders
		self.fields = {}
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
		for x in [a for a in list(set([b for b in re_split('~|\+|-|\*|:|factor|\(|\)',self.formula) if b not in ['1','0']])) if a != '']:
			mtype = "dependent" if x in re_split('factor|\(|\)',re_split('~',self.formula)[0]) else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				self.fields[x] = {'class': 'factor', 'type': mtype, 'dtype': '>f8'}
			else:
				self.fields[x] = {'class': 'numeric', 'type': mtype, 'dtype': '>f8'}

		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'snv'])))

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
		try:
			self.variants = getattr(Geno,self.format.capitalize())(self.variants_file, self.samples_file)
		except Process.Error as err:
			raise Process.Error(err.msg)
		else:
			print "extracting model fields from pheno file and reducing to complete observations ..."
			p_names = (self.fid,self.iid) + tuple(x for x in [self.matid, self.patid] if x is not None)
			p_names = p_names + tuple(x for x in self.fields if x not in [self.fid,self.iid,self.matid,self.patid])
			p_names = p_names + (self.sex,) if self.sex is not None and self.sex not in self.fields else p_names
			p_dtypes = ('|S100', '|S100') + tuple('|S100' for x in [self.matid, self.patid] if x is not None)
			p_dtypes = p_dtypes + tuple(self.fields[x]['dtype'] for x in self.fields if x not in [self.fid,self.iid,self.matid,self.patid])
			p_dtypes = p_dtypes + ('>f8',) if self.sex is not None and self.sex not in self.fields else p_dtypes
			dtypes = dict(zip(p_names, p_dtypes))
			try:
				self.pheno = np.genfromtxt(fname=self.pheno_file, delimiter=Fxns.get_delimiter(self.pheno_sep), dtype=p_dtypes, names=True, usecols=p_names)
			except:
				raise Process.Error("unable to load phenotype file " + self.pheno_file + " with columns " + ', '.join(p_names))
			for x in [y for y in dtypes if dtypes[y] == '>f8']:
				self.pheno = self.pheno[~np.isnan(self.pheno[x])]
			for x in [y for y in dtypes if dtypes[y] == '|S100']:
				self.pheno = self.pheno[~(self.pheno[x] == 'NA')]
			for x in self.fields:
				if x in self.pheno.dtype.names:
					print "   %s variable %s found" % (self.fields[x]['type'], x)
				elif x == 'snv':
					print "   %s variable %s skipped" % (self.fields[x]['type'], x)
				else:
					raise Process.Error("column " + x + " not found in phenotype file " + self.pheno_file)
			if self.fid in self.pheno.dtype.names:
				print "   fid column %s found" % self.fid
			else:
				raise Process.Error("column " + self.fid + " not found in phenotype file " + self.pheno_file)
			if self.iid in self.pheno.dtype.names:
				print "   iid column %s found" % self.iid
			else:
				raise Process.Error("column " + self.iid + " not found in phenotype file " + self.pheno_file)
			if self.matid is not None:
				if self.matid in self.pheno.dtype.names:
					print "   matid column %s found" % self.matid
				else:
					raise Process.Error("column " + self.matid + " not found in phenotype file " + self.pheno_file)
			if self.patid is not None:
				if self.patid in self.pheno.dtype.names:
					print "   patid column %s found" % self.patid
				else:
					raise Process.Error("column " + self.patid + " not found in phenotype file " + self.pheno_file)
			if self.sex is not None:
				if self.sex in self.pheno.dtype.names:
					print "   sex column %s found" % self.sex
				else:
					raise Process.Error("column " + self.sex + " not found in phenotype file " + self.pheno_file)
			self.pheno = self.pheno[np.in1d(self.pheno[self.iid],np.intersect1d(self.pheno[self.iid],self.variants.samples))]
			if self.pheno.shape[0] > 0:
				iids_unique, iids_counts = np.unique(self.pheno[self.iid], return_counts=True)
				fids_unique, fids_counts = np.unique(self.pheno[self.fid], return_counts=True)
				self.unique_idx = np.in1d(self.pheno[self.iid],iids_unique)
				if self.all_founders or self.matid is None or self.patid is None:
					self.founders_idx = self.unique_idx.copy()
				else:
					self.founders_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][(self.pheno[self.unique_idx][self.matid] == '0') & (self.pheno[self.unique_idx][self.patid] == '0')][self.iid])
				self.calc_hwe_idx = self.founders_idx.copy()
				self.nobs = self.pheno.shape[0]
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
					self.nfounders = self.pheno[self.unique_idx][(self.pheno[self.unique_idx][self.matid] == '0') & (self.pheno[self.unique_idx][self.patid] == '0')].shape[0]
					print "   " + str(self.nfounders) + " founders"
				else:
					self.nfounders = self.nunique
					print "   " + str(self.nfounders) + " founders"
				if self.sex is not None:
					if self.male is not None:
						self.nmales = self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.male].shape[0]
						self.male_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.male][self.iid])
						print "   " + str(self.nmales) + " male"
					if self.female is not None:
						self.nfemales = self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.female].shape[0]
						self.female_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.female][self.iid])
						print "   " + str(self.nfemales) + " female"
				if len([v for v in self.fields if self.fields[v]['type'] == 'dependent']) == 1:
					self.dep_var = [v for v in self.fields if self.fields[v]['type'] == 'dependent'][0]
					if len(np.unique(self.pheno[self.dep_var])) == 2:
						self.family = 'binomial'
						self.ctrl_code = min(self.pheno[self.dep_var]) if not self.ctrl_code else self.ctrl_code
						self.case_code = max(self.pheno[self.dep_var]) if not self.case_code else self.case_code
						self.ncases = self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.dep_var] == self.case_code].shape[0]
						self.cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.dep_var] == self.case_code][self.iid])
						if self.sex is not None:
							self.male_cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.male_idx][self.pheno[self.male_idx][self.dep_var] == self.case_code][self.iid])
							self.female_cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.female_idx][self.pheno[self.female_idx][self.dep_var] == self.case_code][self.iid])
						print "   " + str(self.ncases) + " cases"
						self.nctrls = self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.dep_var] == self.ctrl_code].shape[0]
						self.ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.dep_var] == self.ctrl_code][self.iid])
						self.founders_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.founders_idx][self.pheno[self.founders_idx][self.dep_var] == self.ctrl_code][self.iid])
						self.calc_hwe_idx = self.founders_ctrls_idx.copy()
						if self.sex is not None:
							self.male_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.male_idx][self.pheno[self.male_idx][self.dep_var] == self.ctrl_code][self.iid])
							self.female_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.female_idx][self.pheno[self.female_idx][self.dep_var] == self.ctrl_code][self.iid])
						self.pheno[self.dep_var][self.cases_idx] = 1
						self.pheno[self.dep_var][self.ctrls_idx] = 0
						print "   " + str(self.nctrls) + " controls"
			else:
				raise Process.Error("phenotype file and data file contain no common samples")

		self.metadata = '## source: uga' + version + '\n' + \
						'## date: ' + time.strftime("%Y%m%d") + '\n' + \
						'## method: ' + self.fxn + '\n' + \
						'## formula: ' + self.formula + '\n' + \
						'## total observations: ' + str(self.nobs) + '\n' + \
						'## unique samples: ' + str(self.nunique) + '\n' + \
						'## multiple observation samples: ' + str(self.nlongitudinal) + '\n' + \
						'## families: ' + str(self.nfamilies) + '\n' + \
						'## unrelated samples: ' + str(self.nunrelated) + '\n' + \
						'## founders: ' + str(self.nfounders) + '\n' + \
						'## males: ' + str(self.nmales) + '\n' + \
						'## females: ' + str(self.nfemales)

		self.metadata_cc = '## cases: ' + str(self.ncases) + '\n' + \
							'## controls: ' + str(self.nctrls)

	def print_fields(self):
		print "model fields ..."
		print "      {0:>{1}}".format("field", len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format("type", len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format("class", len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))
		print "      {0:>{1}}".format("-----", len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format("----", len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format("-----", len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))
		for k in self.fields:
			print "      {0:>{1}}".format(str(k), len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format(str(self.fields[k]['type']), len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format(str(self.fields[k]['class']), len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))

	def get_region(self, region, id = None):
		logger = logging.getLogger("Model.Model.get_region")
		logger.debug("get_region " + region)
		try:
			self.variants.get_region(region, id)
		except:
			raise

	cpdef calc_variant_stats(self):
		logger = logging.getLogger("Model.Model.calc_variant_stats")
		logger.debug("calc_variant_stats")
		self.variant_stats = np.zeros((self.variants.info.shape[0],1), dtype=[('filter','uint32'),('mac','>f8'),('callrate','>f8'),('freq','>f8'),('freq.case','>f8'),('freq.ctrl','>f8'),('rsq','>f8'),('hwe','>f8')])
		cdef unsigned int i
		if self.variants.chr == 23:
			for i in xrange(self.variants.info.shape[0]):
				self.variant_stats['mac'][i] = Variant.calc_mac(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['callrate'][i] = Variant.calc_callrate(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['freq'][i] = Variant.calc_freqX(male=self.variants.data[self.male_idx,i+1].astype('float64'), female=self.variants.data[self.female_idx,i+1].astype('float64')) if len(self.male_idx) > 0 and len(self.female_idx) > 0 else np.nan
				self.variant_stats['freq.case'][i] = Variant.calc_freqX(male=self.variants.data[self.male_cases_idx,i+1].astype('float64'), female=self.variants.data[self.female_cases_idx,i+1].astype('float64')) if len(self.male_cases_idx) > 0 and len(self.female_cases_idx) > 0 else np.nan
				self.variant_stats['freq.ctrl'][i] = Variant.calc_freqX(male=self.variants.data[self.male_ctrls_idx,i+1].astype('float64'), female=self.variants.data[self.female_ctrls_idx,i+1].astype('float64')) if len(self.male_ctrls_idx) > 0 and len(self.female_ctrls_idx) > 0 else np.nan
				self.variant_stats['rsq'][i] = Variant.calc_rsq(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['hwe'][i] = Variant.calc_hwe(self.variants.data[self.calc_hwe_idx,i+1].astype('float64')) if len(self.calc_hwe_idx) > 0 else np.nan
		else:
			for i in xrange(self.variants.info.shape[0]):
				self.variant_stats['mac'][i] = Variant.calc_mac(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['callrate'][i] = Variant.calc_callrate(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['freq'][i] = Variant.calc_freq(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['freq.case'][i] = Variant.calc_freq(self.variants.data[self.cases_idx,i+1].astype('float64')) if len(self.cases_idx) > 0 else np.nan
				self.variant_stats['freq.ctrl'][i] = Variant.calc_freq(self.variants.data[self.ctrls_idx,i+1].astype('float64')) if len(self.ctrls_idx) > 0 else np.nan
				self.variant_stats['rsq'][i] = Variant.calc_rsq(self.variants.data[self.unique_idx,i+1].astype('float64'))
				self.variant_stats['hwe'][i] = Variant.calc_hwe(self.variants.data[self.calc_hwe_idx,i+1].astype('float64')) if len(self.calc_hwe_idx) > 0 else np.nan

	cpdef filter(self, np.float miss_thresh=None, np.float maf_thresh=None, np.float maxmaf_thresh=None, np.float mac_thresh=None, 
								np.float rsq_thresh=None, np.float hwe_thresh=None, np.float hwe_maf_thresh=None, no_mono=True):
		logger = logging.getLogger("Model.Model.filter")
		logger.debug("filter")
		cdef unsigned int i
		for i in xrange(self.variant_stats.shape[0]):
			if (not miss_thresh is None and not np.isnan(self.variant_stats['callrate'][i]) and self.variant_stats['callrate'][i] < miss_thresh) or (not np.isnan(self.variant_stats['callrate'][i]) and self.variant_stats['callrate'][i] == 0) or np.isnan(self.variant_stats['callrate'][i]):
				self.variant_stats['filter'][i] += 10000
			if not np.isnan(self.variant_stats['freq'][i]): 
				if no_mono and ((self.variant_stats['freq'][i] == 0 or self.variant_stats['freq'][i] == 1) or 
								((self.variant_stats['freq.case'][i] is not None and not np.isnan(self.variant_stats['freq.case'][i]) and 
									(self.variant_stats['freq.case'][i] == 0 or self.variant_stats['freq.case'][i] == 1)) or 
									(self.variant_stats['freq.ctrl'][i] is not None and not np.isnan(self.variant_stats['freq.ctrl'][i]) and 
									(self.variant_stats['freq.ctrl'][i] == 0 or self.variant_stats['freq.ctrl'][i] == 1)))):
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
		self.results_header = np.array(['chr','pos','snv','a1','a2','filter','callrate','rsq','hwe','mac','freq','freq.case','freq.ctrl'])

		self.metadata_snv = '## chr: chromosome' + '\n' + \
							'## pos: chromosomal position' + '\n' + \
							'## snv: snv name' + '\n' + \
							'## a1: reference (coded) allele used for stat calculations' + '\n' + \
							'## a2: alternate (non-coded) allele' + '\n' + \
							'## filter: filter code (+1 = failed hwe, +10 = failed rsq, +100 = failed mac, +1000 = failed maf, +10000 = failed miss' + '\n' + \
							'## callrate: callrate' + '\n' + \
							'## rsq: imputation info metric (calculated only for variants exhibiting probabilistic genotypes, otherwise coded as NA)' + '\n' + \
							'## hwe: Hardy Weinberg p-value (calculated only if dosages are in [0,1,2], otherwise coded as NA; calculated in founders only if pedigree information available, otherwise calculated in all samples)' + '\n' + \
							'## mac: minor allele count' + '\n' + \
							'## freq: reference (coded) allele frequency'

		self.metadata_snv_cc = '## freq.case: reference (coded) allele frequency in cases' + '\n' + \
								'## freq.ctrl: reference (coded) allele frequency in controls'

	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Model.SnvModel.get_snvs")
		logger.debug("get_snvs")
		try:
			self.variants.get_snvs(buffer)
		except:
			raise
		else:
			try:
				self.variants.data = self.variants.data[np.in1d(self.variants.data[:,0],np.intersect1d(self.pheno[self.iid],self.variants.data[:,0]))]
			except:
				raise Process.Error("phenotype file and data file contain no common samples")
			else:
				self.calc_variant_stats()

cdef class SnvgroupModel(Model):
	cdef public str snvgroup_map, metadata_gene
	cdef public unsigned int snvgroup_mac
	def __cinit__(self, snvgroup_map = None, snvgroup_mac = 1.0, **kwargs):
		logger = logging.getLogger("Model.SnvgroupModel.__cinit__")
		logger.debug("initialize SnvgroupModel")
		self.snvgroup_map = snvgroup_map
		self.snvgroup_mac = snvgroup_mac
		super(SnvgroupModel, self).__init__(**kwargs)
		self.tbx_start = 1
		self.tbx_end = 2
		self.results_header = np.array(['chr','start','end','id'])

		if self.snvgroup_map is not None:
			try:
				self.variants.load_snvgroup_map(self.snvgroup_map)
			except Process.Error as err:
				raise Process.Error(err.msg)

		self.metadata_gene = '## chr: chromosome' + '\n' + \
								'## start: start chromosomal position' + '\n' + \
								'## end: end chromosomal position' + '\n' + \
								'## id: snv group name (ie. gene name)'

	cpdef get_snvgroup(self, int buffer, id):
		logger = logging.getLogger("Model.SnvgroupModel.get_snvgroup")
		logger.debug("get_snvgroup")
		try:
			self.variants.get_snvgroup(buffer, id)
		except:
			raise
		else:
			try:
				self.variants.data = self.variants.data[np.in1d(self.variants.data[:,0],np.intersect1d(self.pheno[self.iid],self.variants.data[:,0]))]
			except:
				raise Process.Error("phenotype file and data file contain no common samples")
			else:
				self.calc_variant_stats()

cdef class Score(SnvModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Score.__cinit__")
		logger.debug("initialize Score model")
		super(Score, self).__init__(**kwargs)
		print "setting score test family option to " + self.family

		self.focus = np.array(['snv'])

		if len(self.focus) > 1:
			for x in self.focus:
				self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.nmiss',x + '.ntotal',x + '.effect',x + '.stderr',x + '.or',x + '.p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','nmiss','ntotal','effect','stderr','or','p']))
		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'snv'])))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed singlesnpMeta)' + '\n' + \
						'## *.nmiss: number of missing samples' + '\n' + \
						'## *.ntotal: total number of included samples' + '\n' + \
						'## *.effect: effect size' + '\n' + \
						'## *.stderr: standard error' + '\n' + \
						'## *.or: odds ratio (exp(effect), not provided by seqMeta)' + '\n' + \
						'## *.p: p-value' + '\n#'

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('nmiss','>f8'),('ntotal','>f8'),('effect','>f8'),('stderr','>f8'),('or','>f8'),('p','>f8')])
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['variant_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.variants.info['variant_unique'][passed]), 'gene': list(self.variants.info['variant_unique'][passed])}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% variants])')
			if len(passed) == 1:
				ro.r('colnames(z)<-"' + self.variants.info['variant_unique'][passed][0] + '"')
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
		for i in xrange(self.variants.info.shape[0]):
			if self.variant_stats['freq'][i] > 0.5:
				a1 = self.variants.info['a1'][i]
				a2 = self.variants.info['a2'][i]
				self.variants.info['a1'][i] = a2
				self.variants.info['a2'][i] = a1
				self.variant_stats['freq'][i] = 1 - self.variant_stats['freq'][i]
				self.variant_stats['freq.case'][i] = 1 - self.variant_stats['freq.case'][i]
				self.variant_stats['freq.ctrl'][i] = 1 - self.variant_stats['freq.ctrl'][i]
				self.results['effect'][i] = -1 * self.results['effect'][i]
				self.results['or'][i] = 1 / self.results['or'][i]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Gee(SnvModel):
	cdef public str corstr
	def __cinit__(self, corstr = None, **kwargs):
		logger = logging.getLogger("Model.Gee.__cinit__")
		logger.debug("initialize Gee model")
		self.corstr = corstr if corstr is not None else 'exchangeable'
		super(Gee, self).__init__(**kwargs)
		print "setting gee test family option to " + self.family

		# to get list of model variables from R, use
		# self.focus = np.array([x for x in list(ro.r('labels')(ro.r('terms')(ro.r('formula')(self.formula)))) if 'snv' in x])
		self.focus = np.array(['snv'])

		if len(self.focus) > 1:
			for x in self.focus:
				self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.effect',x + '.stderr',x + '.or',x + '.wald',x + '.p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','or','wald','p']))
		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'snv'])))

		print "loading R package geepack"
		ro.r('suppressMessages(library(geepack))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						'## corstr: ' + self.corstr + '\n' + \
						self.metadata_snv + '\n' + \
						'## *.err: error code (0: no error, 1: geeglm error reported, 2: infinite value or zero p-value detected, 3: geeglm failed)' + '\n' + \
						'## *.effect: effect size' + '\n' + \
						'## *.stderr: standard error' + '\n' + \
						'## *.or: odds ratio (included only if binomial family)' + '\n' + \
						'## *.wald: wald-statistic' + '\n' + \
						'## *.p: p-value' + '\n#'

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('effect','>f8'),('stderr','>f8'),('or','>f8'),('wald','>f8'),('p','>f8')])
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			ro.r('model_df$' + self.fid + '<-as.factor(model_df$' + self.fid + ')')
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['variant_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [self.fid] + [vu]))
				cmd = 'geeglm(' + self.formula + '+' + vu + ',id=' + self.fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + self.family + ',corstr="' + self.corstr + '")'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					raise Process.Error(rerr.message)
				else:
					ro.r('result<-summary(result)')
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
		for i in xrange(self.variants.info.shape[0]):
			if self.variant_stats['freq'][i] > 0.5:
				a1 = self.variants.info['a1'][i]
				a2 = self.variants.info['a2'][i]
				self.variants.info['a1'][i] = a2
				self.variants.info['a2'][i] = a1
				self.variant_stats['freq'][i] = 1 - self.variant_stats['freq'][i]
				self.variant_stats['freq.case'][i] = 1 - self.variant_stats['freq.case'][i]
				self.variant_stats['freq.ctrl'][i] = 1 - self.variant_stats['freq.ctrl'][i]
				self.results['effect'][i] = -1 * self.results['effect'][i]
				self.results['or'][i] = 1 / self.results['or'][i]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Glm(SnvModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Glm.__cinit__")
		logger.debug("initialize Glm model")
		super(Glm, self).__init__(**kwargs)
		print "setting glm test family option to " + self.family

		# to get list of model variables from R, use
		# self.focus = np.array([x for x in list(ro.r('labels')(ro.r('terms')(ro.r('formula')(self.formula)))) if 'snv' in x])
		self.focus = np.array(['snv'])

		if len(self.focus) > 1:
			for x in self.focus:
				self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.effect',x + '.stderr',x + '.or',x + '.z',x + '.p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','or','z','p']))
		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'snv'])))

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata_snv = self.metadata_snv + '\n' + self.metadata_snv_cc if self.family == 'binomial' else self.metadata_snv
		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## *.err: error code (0: no error, 1: missing values returned by glm, 2: infinite value or zero p-value detected, 3: glm failed)' + '\n' + \
						'## *.effect: effect size' + '\n' + \
						'## *.stderr: standard error' + '\n' + \
						'## *.or: odds ratio (included only if binomial family)' + '\n' + \
						'## *.z: z-statistic' + '\n' + \
						'## *.p: p-value' + '\n#'
						

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('effect','>f8'),('stderr','>f8'),('or','>f8'),('z','>f8'),('p','>f8')])
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['variant_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [vu]))
				cmd = 'glm(' + self.formula + '+' + vu + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + self.family + ')'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					raise Process.Error(rerr.message)
				else:
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
		for i in xrange(self.variants.info.shape[0]):
			if self.variant_stats['freq'][i] > 0.5:
				a1 = self.variants.info['a1'][i]
				a2 = self.variants.info['a2'][i]
				self.variants.info['a1'][i] = a2
				self.variants.info['a2'][i] = a1
				self.variant_stats['freq'][i] = 1 - self.variant_stats['freq'][i]
				self.variant_stats['freq.case'][i] = 1 - self.variant_stats['freq.case'][i]
				self.variant_stats['freq.ctrl'][i] = 1 - self.variant_stats['freq.ctrl'][i]
				self.results['effect'][i] = -1 * self.results['effect'][i]
				self.results['or'][i] = 1 / self.results['or'][i]
				self.results['z'][i] = -1 * self.results['z'][i]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Lm(SnvModel):
	def __cinit__(self, **kwargs):
		logger = logging.getLogger("Model.Lm.__cinit__")
		logger.debug("initialize Lm model")
		super(Lm, self).__init__(**kwargs)

		# to get list of model variables from R, use
		# self.focus = np.array([x for x in list(ro.r('labels')(ro.r('terms')(ro.r('formula')(self.formula)))) if 'snv' in x])
		self.focus = np.array(['snv'])

		if len(self.focus) > 1:
			for x in self.focus:
				self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.effect',x + '.stderr',x + '.t',x + '.p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','effect','stderr','t','p']))
		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'snv'])))

		self.metadata = self.metadata + '\n' + \
						self.metadata_snv + '\n' + \
						'## *.err: error code (0: no error, 1: missing values returned by lm, 2: infinite value or zero p-value detected, 3: lm failed)' + '\n' + \
						'## *.effect: effect size' + '\n' + \
						'## *.stderr: standard error' + '\n' + \
						'## *.t: t-statistic' + '\n' + \
						'## *.p: p-value' + '\n#'
						

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.variants.info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('effect','>f8'),('stderr','>f8'),('t','>f8'),('p','>f8')])
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			for v in passed:
				vu = self.variants.info['variant_unique'][v]
				ro.globalenv['cols'] = list(set([a for a in self.model_cols] + [vu]))
				cmd = 'lm(' + self.formula + '+' + vu + ',data=na.omit(model_df[,names(model_df) %in% cols]))'
				try:
					ro.globalenv['result'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][v] = 3
					raise Process.Error(rerr.message)
				else:
					ro.r('result<-summary(result)')
					ro.r('err<-0')
					ro.r('err[result$coefficients["' + vu + '",1] == "NA"]<-1')
					ro.r('err[! is.finite(result$coefficients["' + vu + '",1]) | ! is.finite(result$coefficients["' + vu + '",2]) | ! is.finite(result$coefficients["' + vu + '",4]) | result$coefficients["' + vu + '",4] == 0]<-2')
					self.results['err'][v] = np.array(ro.r('err'))[:,None]
					self.results['effect'][v] = np.array(ro.r('result$coefficients["' + vu + '",1]'))[:,None]
					self.results['stderr'][v] = np.array(ro.r('result$coefficients["' + vu + '",2]'))[:,None]
					self.results['t'][v] = np.array(ro.r('result$coefficients["' + vu + '",3]'))[:,None]
					self.results['p'][v] = np.array(ro.r('result$coefficients["' + vu + '",4]'))[:,None]
		for i in xrange(self.variants.info.shape[0]):
			if self.variant_stats['freq'][i] > 0.5:
				a1 = self.variants.info['a1'][i]
				a2 = self.variants.info['a2'][i]
				self.variants.info['a1'][i] = a2
				self.variants.info['a2'][i] = a1
				self.variant_stats['freq'][i] = 1 - self.variant_stats['freq'][i]
				self.variant_stats['freq.case'][i] = 1 - self.variant_stats['freq.case'][i]
				self.variant_stats['freq.ctrl'][i] = 1 - self.variant_stats['freq.ctrl'][i]
				self.results['effect'][i] = -1 * self.results['effect'][i]
				self.results['t'][i] = -1 * self.results['t'][i]
		self.out = pd.to_numeric(pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.variants.info,self.variant_stats),flatten=True),self.results),flatten=True), dtype='object'),errors='coerce')

cdef class Skat(SnvgroupModel):
	cdef public str skat_wts, skat_method, mafrange
	cdef unsigned int timeout
	def __cinit__(self, skat_wts = None, skat_method = None, mafrange = None, timeout = None, **kwargs):
		logger = logging.getLogger("Model.Skat.__cinit__")
		logger.debug("initialize Skat model")
		self.skat_wts = skat_wts if skat_wts is not None else 'function(maf){dbeta(maf,1,25)}'
		self.skat_method = skat_method if skat_method is not None else 'saddlepoint'
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.timeout = timeout if timeout is not None else 3600
		super(Skat, self).__init__(**kwargs)
		print "setting skat test family option to " + self.family

		self.results_header = np.append(self.results_header,np.array(['mac','err','nmiss','nsnps','cmaf','p','pmin','rho']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## skat wts: ' + self.skat_wts + '\n' + \
						'## skat method: ' + self.skat_method + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## *.mac: minor allele count for the group' + '\n' + \
						'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatMeta timed out' + '\n' + \
						'## *.nmiss: number of missing genotypes in group' + '\n' + \
						'## *.nsnps: number of snps in the group' + '\n' + \
						'## *.cmaf: cumulative minor allele frequency' + '\n' + \
						'## *.p: group p-value' + '\n' + \
						'## *.q: skat q-statistic' + '\n#'

	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('mac','>f8'),('err','>f8'),('nmiss','>f8'),('nsnps','>f8'),('cmaf','>f8'),('p','>f8'),('q','>f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.id
		self.results['nmiss'][0] = np.nan
		self.results['nsnps'][0] = np.nan
		self.results['cmaf'][0] = np.nan
		self.results['p'][0] = np.nan
		self.results['q'][0] = np.nan
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['variant_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.variants.info['variant_unique'][passed]), 'gene': 'NA'}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['mac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['mac'][0] >= self.snvgroup_mac:
				if len(passed) == 1:
					ro.r('colnames(z)<-"' + self.variants.info['variant_unique'][passed][0] + '"')
				cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					raise Process.Error(rerr.message)
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = "tryCatch(expr = { evalWithTimeout(skatMeta(ps,SNPInfo=snp_info,wts=" + self.skat_wts + ",method='" + self.skat_method + "',mafRange=" + self.mafrange + "), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						raise Process.Error(rerr.message)
					
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

	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['mac','err','nmiss','nsnps','cmaf','p','q']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + '.mac','>f8'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnps','>f8'),(tag + '.cmaf','>f8'),(tag + '.p','>f8'),(tag + '.q','>f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')

cdef class Skato(SnvgroupModel):
	cdef public str skat_wts, burden_wts, skat_method, mafrange, skato_rho
	cdef unsigned int timeout
	def __cinit__(self, skat_wts = None, burden_wts = None, skat_method = None, skato_rho = None, mafrange = None, timeout = None, **kwargs):
		logger = logging.getLogger("Model.Skato.__cinit__")
		logger.debug("initialize Skato model")
		self.skat_wts = skat_wts if skat_wts is not None else 'function(maf){dbeta(maf,1,25)}'
		self.burden_wts = burden_wts if burden_wts is not None else 'function(maf){maf < 0.01}'
		self.skat_method = skat_method if skat_method is not None else 'saddlepoint'
		self.skato_rho = skato_rho if skato_rho is not None else 'seq(0,1,0.1)'
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.timeout = timeout if timeout is not None else 3600
		super(Skato, self).__init__(**kwargs)
		print "setting skat-o test family option to " + self.family
		self.results_header = np.append(self.results_header,np.array(['mac','err','nmiss','nsnps','cmaf','p','pmin','rho']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## skat wts: ' + self.skat_wts + '\n' + \
						'## burden wts: ' + self.burden_wts + '\n' + \
						'## skat method: ' + self.skat_method + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## *.mac: minor allele count for the group' + '\n' + \
						'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatOMeta, 4: skatOMeta errflag > 0, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: skatOMeta timed out' + '\n' + \
						'## *.nmiss: number of missing genotypes in group' + '\n' + \
						'## *.nsnps: number of snps in the group' + '\n' + \
						'## *.cmaf: cumulative minor allele frequency' + '\n' + \
						'## *.p: group p-value' + '\n' + \
						'## *.pmin: minimum snp p-value' + '\n' + \
						'## *.rho: skato rho parameter' + '\n#'

	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('mac','>f8'),('err','>f8'),('nmiss','>f8'),('nsnps','>f8'),('cmaf','>f8'),('p','>f8'),('pmin','>f8'),('rho','>f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.id
		self.results['nmiss'][0] = np.nan
		self.results['nsnps'][0] = np.nan
		self.results['cmaf'][0] = np.nan
		self.results['pmin'][0] = np.nan
		self.results['p'][0] = np.nan
		self.results['rho'][0] = np.nan
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['variant_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.variants.info['variant_unique'][passed]), 'gene': 'NA'}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['mac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['mac'][0] >= self.snvgroup_mac:
				if len(passed) == 1:
					ro.r('colnames(z)<-"' + self.variants.info['variant_unique'][passed][0] + '"')
				cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					raise Process.Error(rerr.message)
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = "tryCatch(expr = { evalWithTimeout(skatOMeta(ps,SNPInfo=snp_info,rho=" + self.skato_rho + ",skat.wts=" + self.skat_wts + ",burden.wts=" + self.burden_wts + ",method='" + self.skat_method + "',mafRange=" + self.mafrange + "), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						raise Process.Error(rerr.message)
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

	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['mac','err','nmiss','nsnps','cmaf','p','pmin','rho']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + '.mac','>f8'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnps','>f8'),(tag + '.cmaf','>f8'),(tag + '.p','>f8'),(tag + '.pmin','>f8'),(tag + '.rho','>f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')

cdef class Burden(SnvgroupModel):
	cdef public str mafrange, burden_wts
	cdef unsigned int timeout
	def __cinit__(self, mafrange = None, burden_wts = None, timeout = None, **kwargs):
		logger = logging.getLogger("Model.Burden.__cinit__")
		logger.debug("initialize Burden model")
		self.mafrange = mafrange if mafrange is not None else 'c(0,0.5)'
		self.burden_wts = burden_wts if burden_wts is not None else '1'
		self.timeout = timeout if timeout is not None else 3600
		super(Burden, self).__init__(**kwargs)
		print "setting burden test family option to " + self.family

		self.results_header = np.append(self.results_header,np.array(['mac','err','nmiss','nsnpsTotal','nsnpsUsed','cmafTotal','cmafUsed','beta','se','p']))

		print "loading R package seqMeta"
		ro.r('suppressMessages(library(seqMeta))')

		self.metadata = self.metadata + '\n' + self.metadata_cc if self.family == 'binomial' else self.metadata
		self.metadata = self.metadata + '\n' + \
						'## burden wts: ' + self.burden_wts + '\n' + \
						'## maf range: ' + self.mafrange + '\n' + \
						'## timeout: ' + str(self.timeout) + '\n' + \
						self.metadata_gene + '\n' + \
						'## *.mac: minor allele count for the group' + '\n' + \
						'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed burdenMeta, 5: analysis skipped, 6: failed group minor allele count threshold), 7: prepScores2 timed out, 8: burdenMeta timed out' + '\n' + \
						'## *.nmiss: number of missing snps' + '\n' + \
						'## *.nsnpsTotal: number of snps in group' + '\n' + \
						'## *.nsnpsUsed: number of snps used in analysis' + '\n' + \
						'## *.cmafTotal: cumulative minor allele frequency in group' + '\n' + \
						'## *.cmafUsed: cumulative minor allele frequency for snps used in analysis' + '\n' + \
						'## *.beta: coefficient for effect of genotype' + '\n' + \
						'## *.se: standard error for effect of genotype' + '\n' + \
						'## *.p: p value for burden test' + '\n#'

	cpdef calc_model(self, meta = None):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.variant_stats['filter'] == 0)[0])
		passed_data = list(np.where(self.variant_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('mac','>f8'),('err','>f8'),('nmiss','>f8'),('nsnpsTotal','>f8'),('nsnpsUsed','>f8'),('cmafTotal','>f8'),('cmafUsed','>f8'),('beta','>f8'),('se','>f8'),('p','>f8')])
		self.results['chr'][0] = self.variants.chr
		self.results['start'][0] = self.variants.start
		self.results['end'][0] = self.variants.end
		self.results['id'][0] = self.variants.id
		self.results['nmiss'][0] = np.nan
		self.results['nsnpsTotal'][0] = np.nan
		self.results['nsnpsUsed'][0] = np.nan
		self.results['cmafTotal'][0] = np.nan
		self.results['cmafUsed'][0] = np.nan
		self.results['beta'][0] = np.nan
		self.results['se'][0] = np.nan
		self.results['p'][0] = np.nan
		if len(passed) > 0:
			variants_df = pd.DataFrame(self.variants.data[:,[0] + passed_data],dtype='object')
			variants_df.columns = [self.iid] + list(self.variants.info['variant_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(variants_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.variants.info['variant_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['variants'] = ro.StrVector(list(self.variants.info['variant_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.variants.info['variant_unique'][passed]), 'gene': 'NA'}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% variants])')
			self.results['mac'][0] = np.sum(self.variant_stats['mac'][passed])
			if self.results['mac'][0] >= self.snvgroup_mac:
				if len(passed) == 1:
					ro.r('colnames(z)<-"' + self.variants.info['variant_unique'][passed][0] + '"')
				cmd = "tryCatch(expr = { evalWithTimeout(prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='" + self.family + "'), timeout=" + str(self.timeout) + ") }, error = function(e) return(1))"
				try:
					ro.globalenv['ps'] = ro.r(cmd)
				except RRuntimeError as rerr:
					self.results['err'][0] = 2
					raise Process.Error(rerr.message)
				if ro.r('class(ps) == "seqMeta"')[0]:
					cmd = 'tryCatch(expr = { evalWithTimeout(burdenMeta(ps,SNPInfo=snp_info,mafRange=' + self.mafrange + ',wts=' + self.burden_wts + '), timeout=' + str(self.timeout) + ') }, error = function(e) return(1))'
					try:
						ro.globalenv['result'] = ro.r(cmd)
					except RRuntimeError as rerr:
						self.results['err'][0] = 3
						raise Process.Error(rerr.message)
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

	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['mac','err','nmiss','nsnpsTotal','nsnpsUsed','cmafTotal','cmafUsed','beta','se','p']]))
		self.results = self.results.view(dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),(tag + '.mac','>f8'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnpsTotal','>f8'),(tag + '.nsnpsUsed','>f8'),(tag + '.cmafTotal','>f8'),(tag + '.cmafUsed','>f8'),(tag + '.beta','>f8'),(tag + '.se','>f8'),(tag + '.p','>f8')])
		self.out = pd.to_numeric(pd.DataFrame(self.results.flatten(), dtype='object',index=[0]),errors='coerce')
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')
			ro.r(tag + '_snp_info<-snp_info')

def SkatMeta(Skat obj, tag, meta, meta_incl):
	for m in meta_incl:
		if m == meta_incl[0]:
			ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
		else:
			ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
	results = np.full((1,1), fill_value=np.nan, dtype=[(tag + '.incl','|S100'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnps','>f8'),(tag + '.cmaf','>f8'),(tag + '.p','>f8'),(tag + '.q','>f8')])
	results[tag + '.incl'][0] = ''.join(['+' if a in meta_incl else 'x' for a in meta.split('+')])
	if str(results[tag + '.incl'][0]).count('+') > 1:
		cmd = "tryCatch(expr = { evalWithTimeout(skatMeta(" + ",".join([x + "_ps" for x in meta.split('+') if x in meta_incl]) + ",SNPInfo=snp_info_meta,wts=" + obj.skat_wts + ",method='" + obj.skat_method + "',mafRange=" + obj.mafrange + "), timeout=" + str(obj.timeout) + ") }, error = function(e) return(1))"
		try:
			ro.globalenv['result'] = ro.r(cmd)
		except RRuntimeError as rerr:
			results[tag + '.err'][0] = 3
			raise Process.Error(rerr.message)
		else:
			if ro.r('class(result) == "data.frame"')[0]:
				ro.r('result$err<-0')
				ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
				ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$nsnps[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$q[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$cmaf[! is.na(result$err) & result$err == 1]<-NA')
				results[tag + '.err'][0] = np.array(ro.r('result$err'))[:,None]
				results[tag + '.nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
				results[tag + '.nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
				results[tag + '.cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
				results[tag + '.p'][0] = np.array(ro.r('result$p'))[:,None]
				results[tag + '.q'][0] = np.array(ro.r('result$Q'))[:,None]
			else:
				results[tag + '.err'][0] = 8
	else:
		results[tag + '.err'][0] = 5
		results[tag + '.nmiss'][0] = np.nan
		results[tag + '.nsnps'][0] = np.nan
		results[tag + '.cmaf'][0] = np.nan
		results[tag + '.p'][0] = np.nan
		results[tag + '.q'][0] = np.nan
	return pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')

def SkatoMeta(Skato obj, tag, meta, meta_incl):
	for m in meta_incl:
		if m == meta_incl[0]:
			ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
		else:
			ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
	results = np.full((1,1), fill_value=np.nan, dtype=[(tag + '.incl','|S100'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnps','>f8'),(tag + '.cmaf','>f8'),(tag + '.p','>f8'),(tag + '.pmin','>f8'),(tag + '.rho','>f8')])
	results[tag + '.incl'][0] = ''.join(['+' if a in meta_incl else 'x' for a in meta.split('+')])
	if str(results[tag + '.incl'][0]).count('+') > 1:
		cmd = "tryCatch(expr = { evalWithTimeout(skatOMeta(" + ",".join([x + "_ps" for x in meta.split('+') if x in meta_incl]) + ",SNPInfo=snp_info_meta,rho=" + obj.skato_rho + ",skat.wts=" + obj.skat_wts + ",burden.wts=" + obj.burden_wts + ",method='" + obj.skat_method + "',mafRange=" + obj.mafrange + "), timeout=" + str(obj.timeout) + ") }, error = function(e) return(1))"
		try:
			ro.globalenv['result'] = ro.r(cmd)
		except RRuntimeError as rerr:
			results[tag + '.err'][0] = 3
			raise Process.Error(rerr.message)
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
				results[tag + '.err'][0] = np.array(ro.r('result$err'))[:,None]
				results[tag + '.nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
				results[tag + '.nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
				results[tag + '.cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
				results[tag + '.pmin'][0] = np.array(ro.r('result$pmin'))[:,None]
				results[tag + '.p'][0] = np.array(ro.r('result$p'))[:,None]
				results[tag + '.rho'][0] = np.array(ro.r('result$rho'))[:,None]
			else:
				results[tag + '.err'][0] = 8
	else:
		results[tag + '.err'][0] = 5
		results[tag + '.nmiss'][0] = np.nan
		results[tag + '.nsnps'][0] = np.nan
		results[tag + '.cmaf'][0] = np.nan
		results[tag + '.pmin'][0] = np.nan
		results[tag + '.p'][0] = np.nan
		results[tag + '.rho'][0] = np.nan
	return pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')

def BurdenMeta(Burden obj, tag, meta, meta_incl):
	for m in meta_incl:
		if m == meta_incl[0]:
			ro.globalenv['snp_info_meta'] = ro.r(m + '_snp_info')
		else:
			ro.r('snp_info_meta<-merge(snp_info_meta,' + m + '_snp_info,all=TRUE)')
	results = np.full((1,1), fill_value=np.nan, dtype=[(tag + '.incl','|S100'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnpsTotal','>f8'),(tag + '.nsnpsUsed','>f8'),(tag + '.cmafTotal','>f8'),(tag + '.cmafUsed','>f8'),(tag + '.beta','>f8'),(tag + '.se','>f8'),(tag + '.p','>f8')])
	results[tag + '.incl'][0] = ''.join(['+' if a in meta_incl else 'x' for a in meta.split('+')])
	if str(results[tag + '.incl'][0]).count('+') > 1:
		cmd = 'tryCatch(expr = { evalWithTimeout(burdenMeta(' + ",".join([x + "_ps" for x in meta.split('+') if x in meta_incl]) + ',SNPInfo=snp_info_meta,mafRange=' + obj.mafrange + ',wts=' + obj.burden_wts + '), timeout=' + str(obj.timeout) + ') }, error = function(e) return(1))'
		try:
			ro.globalenv['result'] = ro.r(cmd)
		except RRuntimeError as rerr:
			results[tag + '.err'][0] = 3
			raise Process.Error(rerr.message)
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
				results[tag + '.err'][0] = np.array(ro.r('result$err'))[:,None]
				results[tag + '.nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
				results[tag + '.nsnpsTotal'][0] = np.array(ro.r('result$nsnpsTotal'))[:,None]
				results[tag + '.nsnpsUsed'][0] = np.array(ro.r('result$nsnpsUsed'))[:,None]
				results[tag + '.cmafTotal'][0] = np.array(ro.r('result$cmafTotal'))[:,None]
				results[tag + '.cmafUsed'][0] = np.array(ro.r('result$cmafUsed'))[:,None]
				results[tag + '.beta'][0] = np.array(ro.r('result$beta'))[:,None]
				results[tag + '.se'][0] = np.array(ro.r('result$se'))[:,None]
				results[tag + '.p'][0] = np.array(ro.r('result$p'))[:,None]
			else:
				results[tag + '.err'][0] = 8
	else:
		results[tag + '.err'][0] = 5
		results[tag + '.nmiss'][0] = np.nan
		results[tag + '.nsnpsTotal'][0] = np.nan
		results[tag + '.nsnpsUsed'][0] = np.nan
		results[tag + '.cmafTotal'][0] = np.nan
		results[tag + '.cmafUsed'][0] = np.nan
		results[tag + '.beta'][0] = np.nan
		results[tag + '.se'][0] = np.nan
		results[tag + '.p'][0] = np.nan
	return pd.to_numeric(pd.DataFrame(results.flatten(), dtype='object',index=[0]),errors='coerce')
