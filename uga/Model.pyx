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
from re import split as re_split
from Process import Error
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
ro.r('options(warn=1)')
ro.r('options(na.action=na.omit)')
pd.options.mode.chained_assignment = None

cdef class Model(object):
	cdef public unsigned int case_code, ctrl_code, tbx_start, tbx_end
	cdef public np.ndarray cases_idx, ctrls_idx, male_cases_idx, male_ctrls_idx, female_cases_idx, female_ctrls_idx
	cdef public np.ndarray focus, model_cols, results_header, results_dtypes, marker_stats, results, unique_idx, founders_idx, founders_ctrls_idx, males_idx, females_idx
	cdef public bytes formula, format, pheno_file, biodata_file, gene_map, type, iid, fid, matid, patid, sex, pheno_sep, results_header_metadata, a1, a2
	cdef public unsigned int male, female, nobs, nunique, nfounders, nlongitudinal, nfamilies, nunrelated, ncases, nctrls, nmales, nfemales
	cdef public dict fields
	cdef public object pheno, biodata, out
	cdef public bint all_founders
	def __cinit__(self, formula, format, biodata_file, pheno_file, type, iid, fid, gene_map = None, case_code = None, ctrl_code = None, all_founders = False, matid = None, patid = None, sex = None, male = 1, female = 2, pheno_sep = '\t'):
		self.formula = formula
		self.format = format
		self.biodata_file = biodata_file
		self.pheno_file = pheno_file
		self.gene_map = gene_map
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

	def load(self):
		try:
			self.biodata = getattr(Geno,self.format.capitalize())(self.biodata_file)
		except Error as err:
			raise err
		else:
			print "extracting model fields from pheno file and reducing to complete observations ..."
			p_names = (self.fid,self.iid) + tuple(x for x in [self.matid, self.patid] if x is not None)
			p_names = p_names + tuple(x for x in self.fields if x not in [self.fid,self.iid,self.matid,self.patid])
			p_names = p_names + tuple(self.sex) if self.sex is not None and self.sex not in self.fields else p_names
			p_dtypes = ('|S100', '|S100') + tuple('|S100' for x in [self.matid, self.patid] if x is not None)
			p_dtypes = p_dtypes + tuple(self.fields[x]['dtype'] for x in self.fields if x not in [self.fid,self.iid,self.matid,self.patid])
			p_dtypes = p_dtypes + tuple('>f8') if self.sex is not None and self.sex not in self.fields else p_dtypes
			dtypes = dict(zip(p_names, p_dtypes))
			try:
				self.pheno = np.genfromtxt(fname=self.pheno_file, delimiter=self.pheno_sep, dtype=p_dtypes, names=True, usecols=p_names)
			except:
				raise Error("unable to load phenotype file " + self.pheno_file + " with columns " + ', '.join(p_names))
			for x in [y for y in dtypes if dtypes[y] == '>f8']:
				self.pheno = self.pheno[~np.isnan(self.pheno[x])]
			for x in [y for y in dtypes if dtypes[y] == '|S100']:
				self.pheno = self.pheno[~(self.pheno[x] == 'NA')]
			for x in self.fields:
				if x in self.pheno.dtype.names:
					print "   %s variable %s found" % (self.fields[x]['type'], x)
				elif x in ['marker','marker1','marker2','marker.interact']:
					print "   %s variable %s skipped" % (self.fields[x]['type'], x)
				else:
					raise Error("column " + x + " not found in phenotype file " + self.pheno_file)
			if self.fid in self.pheno.dtype.names:
				print "   fid column %s found" % self.fid
			else:
				raise Error("column " + self.fid + " not found in phenotype file " + self.pheno_file)
			if self.iid in self.pheno.dtype.names:
				print "   iid column %s found" % self.iid
			else:
				raise Error("column " + self.iid + " not found in phenotype file " + self.pheno_file)
			if self.matid is not None:
				if self.matid in self.pheno.dtype.names:
					print "   matid column %s found" % self.matid
				else:
					raise Error("column " + self.matid + " not found in phenotype file " + self.pheno_file)
			if self.patid is not None:
				if self.patid in self.pheno.dtype.names:
					print "   patid column %s found" % self.patid
				else:
					raise Error("column " + self.patid + " not found in phenotype file " + self.pheno_file)
			if self.sex is not None:
				if self.sex in self.pheno.dtype.names:
					print "   sex column %s found" % self.sex
				else:
					raise Error("column " + self.sex + " not found in phenotype file " + self.pheno_file)
			try:
				self.pheno = self.pheno[np.in1d(self.pheno[self.iid],np.intersect1d(self.pheno[self.iid],self.biodata.samples))]
			except:
				raise Error("phenotype file and data file contain no common samples")
			iids_unique, iids_counts = np.unique(self.pheno[self.iid], return_counts=True)
			fids_unique, fids_counts = np.unique(self.pheno[self.fid], return_counts=True)
			self.unique_idx = np.in1d(self.pheno[self.iid],iids_unique)
			if self.all_founders or self.matid is None or self.patid is None:
				self.founders_idx = self.unique_idx.copy()
			else:
				self.founders_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][(self.pheno[self.unique_idx][self.matid] == '0') & (self.pheno[self.unique_idx][self.patid] == '0')][self.iid])
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
					self.males_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.male][self.iid])
					print "   " + str(self.nmales) + " male"
				if self.female is not None:
					self.nfemales = self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.female].shape[0]
					self.females_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][self.sex] == self.female][self.iid])
					print "   " + str(self.nfemales) + " female"
			if hasattr(self, 'case_code') and self.case_code is not None:
				dep_var = [v for v in self.fields if self.fields[v]['type'] == 'dependent'][0]
				self.ncases = self.pheno[self.unique_idx][self.pheno[self.unique_idx][dep_var] == self.case_code].shape[0]
				self.cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][dep_var] == self.case_code][self.iid])
				if self.sex is not None:
					self.male_cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.males_idx][self.pheno[self.males_idx][dep_var] == self.case_code][self.iid])
					self.female_cases_idx = np.in1d(self.pheno[self.iid],self.pheno[self.females_idx][self.pheno[self.females_idx][dep_var] == self.case_code][self.iid])
				print "   " + str(self.ncases) + " cases"
			if hasattr(self, 'ctrl_code') and self.ctrl_code is not None:
				dep_var = [v for v in self.fields if self.fields[v]['type'] == 'dependent'][0]
				self.nctrls = self.pheno[self.unique_idx][self.pheno[self.unique_idx][dep_var] == self.ctrl_code].shape[0]
				self.ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.unique_idx][self.pheno[self.unique_idx][dep_var] == self.ctrl_code][self.iid])
				self.founders_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.founders_idx][self.pheno[self.founders_idx][dep_var] == self.ctrl_code][self.iid])
				if self.sex is not None:
					self.male_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.males_idx][self.pheno[self.males_idx][dep_var] == self.ctrl_code][self.iid])
					self.female_ctrls_idx = np.in1d(self.pheno[self.iid],self.pheno[self.females_idx][self.pheno[self.females_idx][dep_var] == self.ctrl_code][self.iid])
				print "   " + str(self.nctrls) + " controls"

	def get_region(self, region, id = None):
		try:
			self.biodata.get_region(region, id)
		except:
			raise

	cpdef get_chunk(self, int buffer):
		try:
			self.biodata.get_chunk(buffer)
		except:
			raise
		else:
			try:
				self.biodata.marker_data = self.biodata.marker_data[np.in1d(self.biodata.marker_data[:,0],np.intersect1d(self.pheno[self.iid],self.biodata.marker_data[:,0]))]
			except:
				raise Error("phenotype file and data file contain no common samples")
			else:
				self.calc_marker_stats()

	cpdef get_gene(self, gene):
		try:
			self.biodata.get_gene(gene)
		except:
			raise
		else:
			try:
				self.biodata.marker_data = self.biodata.marker_data[np.in1d(self.biodata.marker_data[:,0],np.intersect1d(self.pheno[self.iid],self.biodata.marker_data[:,0]))]
			except:
				raise Error("phenotype file and data file contain no common samples")
			else:
				self.calc_marker_stats()

	def print_fields(self):
		print "model fields ..."
		print "      {0:>{1}}".format("field", len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format("type", len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format("class", len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))
		print "      {0:>{1}}".format("-----", len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format("----", len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format("-----", len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))
		for k in self.fields:
			print "      {0:>{1}}".format(str(k), len(max(["field"] + [key for key in self.fields.keys()],key=len))) + "   " + "{0:>{1}}".format(str(self.fields[k]['type']), len(max(["type"] + [self.fields[key]['type'] for key in self.fields],key=len))) + "   " + "{0:>{1}}".format(str(self.fields[k]['class']), len(max(["class"] + [self.fields[key]['class'] for key in self.fields],key=len)))

	cdef add_id_col(self):
		self.results_header = np.append(self.results_header,np.array(['id']))

	cpdef calc_marker_stats(self):
		self.marker_stats = np.zeros((self.biodata.marker_info.shape[0],1), dtype=[('filter','uint32'),('mac','uint32'),('callrate','>f8'),('freq','>f8'),('freq.case','>f8'),('freq.ctrl','>f8'),('rsq','>f8'),('hwe','>f8')])
		cdef unsigned int i
		if self.biodata.chr == 23:
			for i in xrange(self.biodata.marker_info.shape[0]):
				self.marker_stats['mac'][i] = Variant.CalcMAC(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['callrate'][i] = Variant.CalcCallrate(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['rsq'][i] = Variant.CalcRsq(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['freq'][i] = Variant.CalcFreqX(male=self.biodata.marker_data[self.male_idx,i+1].astype('float64'), female=self.biodata.marker_data[self.female_idx,i+1].astype('float64')) if len(self.male_idx) > 0 and len(self.female_idx) > 0 else np.nan
				self.marker_stats['freq.case'][i] = Variant.CalcFreqX(male=self.biodata.marker_data[self.male_cases_idx,i+1].astype('float64'), female=self.biodata.marker_data[self.female_cases_idx,i+1].astype('float64')) if len(self.male_cases_idx) > 0 and len(self.female_cases_idx) > 0 else np.nan
				self.marker_stats['freq.ctrl'][i] = Variant.CalcFreqX(male=self.biodata.marker_data[self.male_ctrls_idx,i+1].astype('float64'), female=self.biodata.marker_data[self.female_ctrls_idx,i+1].astype('float64')) if len(self.male_ctrls_idx) > 0 and len(self.female_ctrls_idx) > 0 else np.nan
				self.marker_stats['hwe'][i] = Variant.CalcHWE(self.biodata.marker_data[self.founders_ctrls_idx,i+1].astype('float64')) if len(self.founders_ctrls_idx) > 0 else np.nan
		else:
			for i in xrange(self.biodata.marker_info.shape[0]):
				self.marker_stats['mac'][i] = Variant.CalcMAC(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['callrate'][i] = Variant.CalcCallrate(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['rsq'][i] = Variant.CalcRsq(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['freq'][i] = Variant.CalcFreq(self.biodata.marker_data[self.unique_idx,i+1].astype('float64'))
				self.marker_stats['freq.case'][i] = Variant.CalcFreq(self.biodata.marker_data[self.cases_idx,i+1].astype('float64')) if len(self.cases_idx) > 0 else np.nan
				self.marker_stats['freq.ctrl'][i] = Variant.CalcFreq(self.biodata.marker_data[self.ctrls_idx,i+1].astype('float64')) if len(self.ctrls_idx) > 0 else np.nan
				self.marker_stats['hwe'][i] = Variant.CalcHWE(self.biodata.marker_data[self.founders_ctrls_idx,i+1].astype('float64')) if len(self.founders_ctrls_idx) > 0 else np.nan

	cpdef filter(self, np.float miss_thresh=None, np.float maf_thresh=None, np.float maxmaf_thresh=None, np.float mac_thresh=None, 
								np.float rsq_thresh=None, np.float hwe_thresh=None, np.float hwe_maf_thresh=None, no_mono=True):
		cdef unsigned int i
		for i in xrange(self.marker_stats.shape[0]):
			if (not miss_thresh is None and not np.isnan(self.marker_stats['callrate'][i]) and self.marker_stats['callrate'][i] < miss_thresh) or (not np.isnan(self.marker_stats['callrate'][i]) and self.marker_stats['callrate'][i] == 0) or np.isnan(self.marker_stats['callrate'][i]):
				self.marker_stats['filter'][i] += 10000
			if not np.isnan(self.marker_stats['freq'][i]): 
				if no_mono and ((self.marker_stats['freq'][i] == 0 or self.marker_stats['freq'][i] == 1) or ((self.marker_stats['freq.case'][i] is not None and not np.isnan(self.marker_stats['freq.case'][i]) and (self.marker_stats['freq.case'][i] == 0 or self.marker_stats['freq.case'][i] == 1)) or (self.marker_stats['freq.ctrl'][i] is not None and not np.isnan(self.marker_stats['freq.ctrl'][i]) and (self.marker_stats['freq.ctrl'][i] == 0 or self.marker_stats['freq.ctrl'][i] == 1)))):
					self.marker_stats['filter'][i] += 1000
				else:
					if ((	not maf_thresh is None
						and 
							(		self.marker_stats['freq'][i] < maf_thresh
								or self.marker_stats['freq'][i] > 1-maf_thresh
							)
					) 
					or
					(	not maxmaf_thresh is None
						and    
							(		self.marker_stats['freq'][i] >= maxmaf_thresh
								and self.marker_stats['freq'][i] <= 1-maxmaf_thresh
							)
					)):
						self.marker_stats['filter'][i] += 1000
			if not mac_thresh is None and not np.isnan(self.marker_stats['mac'][i]) and (self.marker_stats['mac'][i] < mac_thresh):
				self.marker_stats['filter'][i] += 100
			if not rsq_thresh is None and not np.isnan(self.marker_stats['rsq'][i]) and (self.marker_stats['rsq'][i] < rsq_thresh):
				self.marker_stats['filter'][i] += 10
			if not hwe_thresh is None and not hwe_maf_thresh is None and not np.isnan(self.marker_stats['hwe'][i]) and not np.isnan(self.marker_stats['freq'][i]) and ((self.marker_stats['freq'][i] <= 0.5 and self.marker_stats['freq'][i] > hwe_maf_thresh and self.marker_stats['hwe'][i] < hwe_thresh) or (self.marker_stats['freq'][i] > 0.5 and 1-self.marker_stats['freq'][i] > hwe_maf_thresh and self.marker_stats['hwe'][i] < hwe_thresh)):
				self.marker_stats['filter'][i] += 1

cdef class Bssmeta(Model):
	def __cinit__(self, formula, format, pheno_file, biodata_file, type, iid, fid, case_code=2, ctrl_code=1, gene_map = None, all_founders = False, matid = None, patid = None, sex = None, male = 1, female = 2, pheno_sep='\t'):
		super(Bssmeta, self).__init__(formula=formula, format=format, case_code=case_code, ctrl_code=ctrl_code, all_founders=all_founders, pheno_file=pheno_file, biodata_file=biodata_file, type=type, iid=iid, fid=fid, matid=matid, patid=patid, sex=sex, male=male, female=female, pheno_sep=pheno_sep)
		self.results_header = np.array(['chr','pos','marker','a1','a2','filter','callrate','mac','freq','freq.case','freq.ctrl','rsq','hwe'])
		self.out = None
		self.tbx_start = 1
		self.tbx_end = 1

		# set marker as the focus variable (which is only statistic provided by singlesnpMeta())
		# to get list of model variables from R, use
		# self.focus = np.array([x for x in list(ro.r('labels')(ro.r('terms')(ro.r('formula')(self.formula)))) if 'marker' in x])
		self.focus = np.array(['marker'])

		# parse formula into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		for x in [a for a in list(set([b for b in re_split('~|\+|-|\*|:|factor|\(|\)',self.formula) if b not in ['1','0']])) if a != '']:
			mtype = "dependent" if x in re_split('factor|\(|\)',re_split('~',self.formula)[0]) else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				self.fields[x] = {'class': 'factor', 'type': mtype, 'dtype': '>f8'}
			else:
				self.fields[x] = {'class': 'numeric', 'type': mtype, 'dtype': '>f8'}
		if len(self.focus) > 1:
			for x in self.focus:
				self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.nmiss',x + '.ntotal',x + '.effect',x + '.stderr',x + '.or',x + '.p']))
		else:
			self.results_header = np.append(self.results_header,np.array(['err','nmiss','ntotal','effect','stderr','or','p']))
		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'marker'])))

		from rpy2.robjects.packages import importr
		print "loading R package seqMeta"
		importr('seqMeta')

		try:
			self.load()
		except Error as err:
			raise err

		self.results_header_metadata = '## source: uga' + version + '\n' + \
									'## date: ' + time.strftime("%Y%m%d") + '\n' + \
									'## method: bssmeta' + '\n' + \
									'## formula: ' + self.formula + '\n' + \
									'## total observations: ' + str(self.nobs) + '\n' + \
									'## unique samples: ' + str(self.nunique) + '\n' + \
									'## multiple observation samples: ' + str(self.nlongitudinal) + '\n' + \
									'## families: ' + str(self.nfamilies) + '\n' + \
									'## unrelated samples: ' + str(self.nunrelated) + '\n' + \
									'## founders: ' + str(self.nfounders) + '\n' + \
									'## males: ' + str(self.nmales) + '\n' + \
									'## females: ' + str(self.nfemales) + '\n' + \
									'## cases: ' + str(self.ncases) + '\n' + \
									'## controls: ' + str(self.nctrls) + '\n' + \
									'## chr: chromosome' + '\n' + \
									'## pos: chromosomal position' + '\n' + \
									'## marker: variant name' + '\n' + \
									'## a1: reference (coded) allele used for stat calculations' + '\n' + \
									'## a2: alternate (non-coded) allele' + '\n' + \
									'## filter: filter code (+1 = failed hwe, +10 = failed rsq, +100 = failed mac, +1000 = failed maf, +10000 = failed miss' + '\n' + \
									'## callrate: callrate' + '\n' + \
									'## mac: minor allele count (calculated only if dosages are in [0,1,2], otherwise coded as NA)' + '\n' + \
									'## freq: reference (coded) allele frequency' + '\n' + \
									'## freq.case: reference (coded) allele frequency in cases' + '\n' + \
									'## freq.ctrl: reference (coded) allele frequency in controls' + '\n' + \
									'## rsq: imputation info metric (calculated only for variants exhibiting probabilistic genotypes, otherwise coded as NA)' + '\n' + \
									'## hwe: Hardy Weinberg p-value (calculated only if dosages are in [0,1,2], otherwise coded as NA; calculated in founders only if pedigree information available, otherwise calculated in all samples)' + '\n' + \
									'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed singlesnpMeta)' + '\n' + \
									'## *.nmiss: number of missing samples' + '\n' + \
									'## *.ntotal: total number of included samples' + '\n' + \
									'## *.effect: effect size' + '\n' + \
									'## *.stderr: standard error' + '\n' + \
									'## *.or: odds ratio (exp(effect), not provided by seqMeta)' + '\n' + \
									'## *.p: p-value' + '\n#'

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.marker_stats['filter'] == 0)[0])
		passed_marker_data = list(np.where(self.marker_stats['filter'] == 0)[0]+1)
		self.results = np.full((self.biodata.marker_info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('nmiss','>f8'),('ntotal','>f8'),('effect','>f8'),('stderr','>f8'),('or','>f8'),('p','>f8')])
		if len(passed) > 0:
			biodata_df = pd.DataFrame(self.biodata.marker_data[:,[0] + passed_marker_data],dtype='object')
			biodata_df.columns = [self.iid] + list(self.biodata.marker_info['marker_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(biodata_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.biodata.marker_info['marker_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['markers'] = ro.StrVector(list(self.biodata.marker_info['marker_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.biodata.marker_info['marker_unique'][passed]), 'gene': list(self.biodata.marker_info['marker_unique'][passed])}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% markers])')
			if len(passed) == 1:
				ro.r('colnames(z)<-"' + self.biodata.marker_info['marker_unique'][passed][0] + '"')
			cmd = "prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='binomial')"
			try:
				ro.globalenv['ps'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][passed] = 2
				raise Error(rerr.message)
			cmd = 'singlesnpMeta(ps,SNPInfo=snp_info)'
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][passed] = 3
				raise Error(rerr.message)
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
			for i in xrange(self.biodata.marker_info.shape[0]):
				if self.marker_stats['freq'][i] > 0.5:
					a1 = self.biodata.marker_info['a1'][i]
					a2 = self.biodata.marker_info['a2'][i]
					self.biodata.marker_info['a1'][i] = a2
					self.biodata.marker_info['a2'][i] = a1
					self.marker_stats['freq'][i] = 1 - self.marker_stats['freq'][i]
					self.marker_stats['freq.case'][i] = 1 - self.marker_stats['freq.case'][i]
					self.marker_stats['freq.ctrl'][i] = 1 - self.marker_stats['freq.ctrl'][i]
					self.results['effect'][i] = -1 * self.results['effect'][i]
					self.results['or'][i] = 1 / self.results['or'][i]
		self.out = pd.DataFrame(recfxns.merge_arrays((recfxns.merge_arrays((self.biodata.marker_info,self.marker_stats),flatten=True),self.results),flatten=True), dtype='object').convert_objects(convert_numeric=True)

cdef class Bskato(Model):
	def __cinit__(self, formula, format, pheno_file, biodata_file, gene_map, type, iid, fid, case_code=2, ctrl_code=1, all_founders=False, matid = None, patid = None, sex = None, male = 1, female = 2, pheno_sep='\t'):
		super(Bskato, self).__init__(formula=formula, format=format, case_code=case_code, ctrl_code=ctrl_code, all_founders=all_founders, gene_map=gene_map, pheno_file=pheno_file, biodata_file=biodata_file, type=type, iid=iid, fid=fid, matid=matid, patid=patid, sex=sex, male=male, female=female, pheno_sep=pheno_sep)
		self.results_header = np.array(['chr','start','end','id'])
		self.out = None
		self.tbx_start = 1
		self.tbx_end = 2

		# set marker as the focus variable (which is only statistic provided by singlesnpMeta())
		# to get list of model variables from R, use
		# self.focus = np.array([x for x in list(ro.r('labels')(ro.r('terms')(ro.r('formula')(self.formula)))) if 'marker' in x])
		#self.focus = np.array(['marker'])

		# parse formula into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		for x in [a for a in list(set([b for b in re_split('~|\+|-|\*|:|factor|\(|\)',self.formula) if b not in ['1','0']])) if a != '']:
			mtype = "dependent" if x in re_split('factor|\(|\)',re_split('~',self.formula)[0]) else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				self.fields[x] = {'class': 'factor', 'type': mtype, 'dtype': '>f8'}
			else:
				self.fields[x] = {'class': 'numeric', 'type': mtype, 'dtype': '>f8'}
		#if len(self.focus) > 1:
		#	for x in self.focus:
		#		self.results_header = np.append(self.results_header,np.array([x + '.err',x + '.nmiss',x + '.nsnps',x + '.cmaf',x + '.p',x + '.pmin',x + '.rho']))
		#else:
		self.results_header = np.append(self.results_header,np.array(['mac','err','nmiss','nsnps','cmaf','p','pmin','rho']))

		self.model_cols = np.array(list(set([a for a in self.fields.keys() if a != 'marker'])))

		from rpy2.robjects.packages import importr
		print "loading R package seqMeta"
		importr('seqMeta')

		try:
			self.load()
		except Error as err:
			raise err

		try:
			self.biodata.load_gene_map(self.gene_map)
		except Error as err:
			raise err

		self.results_header_metadata = '## source: uga' + version + '\n' + \
									'## date: ' + time.strftime("%Y%m%d") + '\n' + \
									'## method: bskato' + '\n' + \
									'## formula: ' + self.formula + '\n' + \
									'## total observations: ' + str(self.nobs) + '\n' + \
									'## unique samples: ' + str(self.nunique) + '\n' + \
									'## multiple observation samples: ' + str(self.nlongitudinal) + '\n' + \
									'## families: ' + str(self.nfamilies) + '\n' + \
									'## unrelated samples: ' + str(self.nunrelated) + '\n' + \
									'## founders: ' + str(self.nfounders) + '\n' + \
									'## males: ' + str(self.nmales) + '\n' + \
									'## females: ' + str(self.nfemales) + '\n' + \
									'## cases: ' + str(self.ncases) + '\n' + \
									'## controls: ' + str(self.nctrls) + '\n' + \
									'## chr: chromosome' + '\n' + \
									'## start: start chromosomal position' + '\n' + \
									'## end: end chromosomal position' + '\n' + \
									'## id: region id (ie. gene)' + '\n' + \
									'## *.mac: minor allele count for the gene' + '\n' + \
									'## *.err: error code (0: no error, 1: infinite value or zero p-value detected, 2: failed prepScores2, 3: failed skatOMeta)' + '\n' + \
									'## *.nmiss: number of missing genotypes in gene' + '\n' + \
									'## *.nsnps: number of snps in the gene' + '\n' + \
									'## *.cmaf: cumulative minor allele frequency' + '\n' + \
									'## *.p: gene p-value' + '\n' + \
									'## *.pmin: minimum snp p-value' + '\n' + \
									'## *.rho: skato rho parameter' + '\n#'

	cpdef calc_model(self):
		pheno_df = pd.DataFrame(self.pheno,dtype='object')
		passed = list(np.where(self.marker_stats['filter'] == 0)[0])
		passed_marker_data = list(np.where(self.marker_stats['filter'] == 0)[0]+1)
		self.results = np.full((1,1), fill_value=np.nan, dtype=[('chr','|S100'),('start','|S100'),('end','|S100'),('id','|S100'),('mac','>f8'),('err','>f8'),('nmiss','>f8'),('nsnps','>f8'),('cmaf','>f8'),('p','>f8'),('pmin','>f8'),('rho','>f8')])
		self.results['chr'][0] = self.biodata.chr
		self.results['start'][0] = self.biodata.start
		self.results['end'][0] = self.biodata.end
		self.results['id'][0] = self.biodata.id
		if len(passed) > 0:
			biodata_df = pd.DataFrame(self.biodata.marker_data[:,[0] + passed_marker_data],dtype='object')
			biodata_df.columns = [self.iid] + list(self.biodata.marker_info['marker_unique'][passed])
			ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(biodata_df, on=self.iid, how='left'), strings_as_factors=False)
			for col in list(self.model_cols) + list(self.biodata.marker_info['marker_unique'][passed]):
				ro.r('class(model_df$' + col + ')<-"numeric"')
			ro.globalenv['markers'] = ro.StrVector(list(self.biodata.marker_info['marker_unique'][passed]))
			ro.globalenv['model_cols'] = ro.StrVector(list(self.model_cols))
			ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': list(self.biodata.marker_info['marker_unique'][passed]), 'gene': 'NA'}), strings_as_factors=False)
			ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% markers])')
			self.results['mac'][0] = ro.r('sum(z)')
			if len(passed) == 1:
				ro.r('colnames(z)<-"' + self.biodata.marker_info['marker_unique'][passed][0] + '"')
			cmd = "prepScores2(Z=z,formula=" + self.formula + ",SNPInfo=snp_info,data=model_df,family='binomial')"
			try:
				ro.globalenv['ps'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][0] = 2
				raise Error(rerr.message)
			cmd = 'skatOMeta(ps,SNPInfo=snp_info,rho=seq(0,1,0.1))'
			try:
				ro.globalenv['result'] = ro.r(cmd)
			except RRuntimeError as rerr:
				self.results['err'][0] = 3
				raise Error(rerr.message)
			else:
				ro.r('result$err<-0')
				ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
				ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$nsnps[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$pmin[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$rho[! is.na(result$err) & result$err == 1]<-NA')
				ro.r('result$cmaf[! is.na(result$err) & result$err == 1]<-NA')
				self.results['err'][0] = np.array(ro.r('result$err'))[:,None]
				self.results['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
				self.results['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
				self.results['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
				self.results['pmin'][0] = np.array(ro.r('result$pmin'))[:,None]
				self.results['p'][0] = np.array(ro.r('result$p'))[:,None]
				self.results['rho'][0] = np.array(ro.r('result$rho'))[:,None]
		self.out = pd.DataFrame(self.results.flatten(), dtype='object',index=[0]).convert_objects(convert_numeric=True)

	cpdef tag_results(self, tag):
		self.results_header = np.append(np.array(['chr','start','end','id']),np.array([tag + '.' + x for x in ['mac','err','nmiss','nsnps','cmaf','p','pmin','rho']]))
		self.results = self.results.view(dtype=[('chr','|S100'),('start','|S100'),('end','|S100'),('id','|S100'),(tag + '.mac','>f8'),(tag + '.err','>f8'),(tag + '.nmiss','>f8'),(tag + '.nsnps','>f8'),(tag + '.cmaf','>f8'),(tag + '.p','>f8'),(tag + '.pmin','>f8'),(tag + '.rho','>f8')])
		if not np.isnan(self.results[tag + '.p'][0]):
			ro.r(tag + '_ps<-ps')

def align_biodata(model_order, models):
	biodata_db = {}
	for m in model_order:
		i = 0
		for row in models[m].biodata.marker_info:
			i += 1
			if not row['uid'][0] in biodata_db:
				biodata_db[row['uid']] = {}
				biodata_db[row['uid']]['marker'] = row['marker']
				biodata_db[row['uid']]['marker_unique'] = row['marker_unique']
				biodata_db[row['uid']]['a1'] = row['a1'][0]
				biodata_db[row['uid']]['a2'] = row['a2'][0]
			else:
				if (biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] == Geno.complement(row['a2'][0]) + Geno.complement(row['a1'][0]) or biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] == row['a2'] + row['a1'] or biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] == Geno.complement(row['a2'][0]) + 'NA' or biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] == row['a2'] + 'NA') and biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] != "AT" and biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] != "TA" and biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] != "GC" and biodata_db[row['uid']]['a1'] + biodata_db[row['uid']]['a2'] != "CG":
					models[m].biodata.marker_info['a1'][i] = biodata_db[row['uid']]['a1']
					models[m].biodata.marker_info['a2'][i] = biodata_db[row['uid']]['a2']
					models[m].biodata.marker_info['marker'][i] = biodata_db[row['uid']]['marker']
					models[m].biodata.marker_info['marker_unique'][i] = biodata_db[row['uid']]['marker_unique']
					models[m].biodata.marker_data[:,i+1] = 2.0 - models[m].biodata.marker_data[:,i+1].astype(float)



	###### UPDATE DATABASE AND FLIP MARKERS WHERE NECESSARY #####
	#if len(cfg['model_order']) > 1:
	#	##### GENERATE ANALOGS FOR ALL MARKERS IN CHUNK #####
	#	chunkdf['analogs'] = chunkdf.apply(lambda row: MiscFxnsCy.ListCompatibleMarkersCy(row['chr'],row['pos'],row['a1'],row['a2'],'><'),axis=1)
	#	if k == cfg['model_order'][0] or mdb is None:
	#		##### GENERATE INITIAL REFERENCE MARKER DATABASE #####
	#		mdb = pd.DataFrame({'chr': list(chunkdf['chr']),'pos': list(chunkdf['pos']),'marker': list(chunkdf['marker']),'a1': chunkdf['a1'].astype(str).replace('-','_'),'a2': chunkdf['a2'].astype(str).replace('-','_'), 
	#									'analogs': chunkdf.apply(lambda row: MiscFxnsCy.ListCompatibleMarkersCy(row['chr'],row['pos'],row['a1'],row['a2'],'><'),axis=1)})
	#	else:
	#		##### FLIP CHUNK MARKERS IF NECESSARY #####
	#		chunkdf = chunkdf.apply(lambda row: MiscFxns.ChunkFlip(row, mdb),1)
	#		##### UPDATE REFERENCE MARKER DATABASE #####
	#		mdb = mdb.apply(lambda row: MiscFxns.MdbUpdate(row,chunkdf),1)
	#		##### APPEND ANY NEW MARKERS TO REFERENCE MARKER DATABASE #####
	#		mdb = MiscFxns.MdbAppend(mdb, chunkdf)
	#	##### FILL IN ANY MISSING ALT ALLELES FROM PREVIOUS RESULTS WITH ALT ALLELES FROM UPDATED MARKER DATABASE #####
	#	for t in cfg['model_order'][0:cfg['model_order'].index(k)]:
	#		if 'a2' in cfg['models'][t]['results'] and cfg['models'][t]['results'][cfg['models'][t]['results']['a2'] == 'NA'].shape[0] > 0:
	#			cfg['models'][t]['results']['a2'][cfg['models'][t]['results']['a2'] == 'NA'] = cfg['models'][t]['results'][cfg['models'][t]['results']['a2'] == 'NA'].apply(lambda row: MiscFxns.UpdateAltAllele(row, mdb), 1)		

	#cpdef calc_meta(meta_tags, metas, meta_incl):
	#	meta_dtypes = []
	#	meta_dtypes.extend([('chr','|S100'),('start','|S100'),('end','|S100'),('id','|S100')])
	#	for meta in meta_tags:
	#		meta_dtypes.extend([(meta + '.err','>f8'),(meta + '.nmiss','>f8'),(meta + '.nsnps','>f8'),(meta + '.cmaf','>f8'),(meta + '.p','>f8'),(meta + '.pmin','>f8'),(meta + '.rho','>f8')])
	#	self.meta_results = np.full((1,1), fill_value=np.nan, dtypes=meta_dtypes)
	#	self.meta_results['chr'][0] = self.biodata.chr
	#	self.meta_results['start'][0] = self.biodata.start
	#	self.meta_results['end'][0] = self.biodata.end
	#	self.meta_results['id'][0] = self.biodata.id














		#self.out=pd.DataFrame({'chr': [self.biodata.chr], 'start': [self.biodata.start], 'end': [self.biodata.end], 'id': [self.biodata.id], 'err': [np.nan], 'nmiss': [np.nan], 'nsnps': [np.nan], 'cmaf': [np.nan], 'p': [np.nan], 'pmin': [np.nan], 'rho': [np.nan]})
		#passed = list(np.where(self.marker_stats['filter'] == 0)[0])
		#if len(passed) > 1:
		#	pheno_df = pd.DataFrame(self.pheno,dtype='object')
		#	biodata_df = pd.DataFrame(self.biodata.marker_data[:,[0] + passed],dtype='object')
		#	biodata_df.columns = [self.iid] + list(self.biodata.marker_info['marker_unique'][np.where(self.marker_stats['filter'] == 0)[0]])
		#	ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(pheno_df.merge(biodata_df, on=self.iid, how='left'), strings_as_factors=False)
		#	#self.results = np.full((self.biodata.marker_info.shape[0],1), fill_value=np.nan, dtype=[('err','>f8'),('nmiss','>f8'),('nsnps','>f8'),('cmaf','>f8'),('p','>f8'),('pmin','>f8'),('rho','>f8')])
		#	#self.results = np.empty((1,1),dtype=[('chr','uint8'),('start','uint32'),('end','uint32'),('id','|S100'),('err','>f8'),('nmiss','>f8'),('nsnps','>f8'),('cmaf','>f8'),('p','>f8'),('pmin','>f8'),('rho','>f8')])
		#	ro.globalenv['markers'] = ro.StrVector(self.biodata.marker_info['marker_unique'][np.where(self.marker_stats['filter'] == 0)[0]])
		#	ro.globalenv['model_cols'] = ro.StrVector(self.model_cols)
		#	ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': self.biodata.marker_info['marker_unique'][np.where(self.marker_stats['filter'] == 0)[0]], 'gene': 'NA'}), strings_as_factors=False)
		#	ro.globalenv['z'] = ro.r('data.matrix(model_df[,names(model_df) %in% markers,drop=FALSE])')
		#	ro.globalenv['pheno'] = ro.r('model_df[model_cols]')
		#	cmd = 'prepScores2(Z=z,formula=' + self.formula + ',SNPInfo=snp_info,data=pheno,family="binomial")'
		#	try:
		#		ro.globalenv['ps'] = ro.r(cmd)
		#	except RRuntimeError as rerr:
		#		self.out['err'][0] = 2
		#		raise Error(rerr.message)
		#	cmd = 'skatOMeta(ps,SNPInfo=snp_info)'
		#	try:
		#		ro.globalenv['result'] = ro.r(cmd)
		#	except RRuntimeError as rerr:
		#		self.out['err'][0] = 3
		#		raise Error(rerr.message)
		#	else:
		#		ro.r('result$err<-0')
		#		ro.r('result$err[! is.finite(result$p) | result$p == 0]<-1')
		#		ro.r('result$nmiss[! is.na(result$err) & result$err == 1]<-NA')
		#		ro.r('result$nsnps[! is.na(result$err) & result$err == 1]<-NA')
		#		ro.r('result$p[! is.na(result$err) & result$err == 1]<-NA')
		#		ro.r('result$pmin[! is.na(result$err) & result$err == 1]<-NA')
		#		ro.r('result$rho[! is.na(result$err) & result$err == 1]<-NA')
		#		ro.r('result$cmaf[! is.na(result$err) & result$err == 1]<-NA')
		#		self.out['err'][0] = np.array(ro.r('result$err'))[:,None]
		#		self.out['nmiss'][0] = np.array(ro.r('result$nmiss'))[:,None]
		#		self.out['nsnps'][0] = np.array(ro.r('result$nsnps'))[:,None]
		#		self.out['cmaf'][0] = np.array(ro.r('result$cmaf'))[:,None]
		#		self.out['pmin'][0] = np.array(ro.r('result$pmin'))[:,None]
		#		self.out['p'][0] = np.array(ro.r('result$p'))[:,None]
		#		self.out['rho'][0] = np.array(ro.r('result$rho'))[:,None]

"""
class BglmFrame(ModelFrame):
	def __init__(self, formula, type):
		super(BglmFrame, self).__init__(formula, type)
		self.header = ['chr','pos','a1','a2','marker','callrate','freq','freq.case','freq.ctrl','mac','rsq','hwe','filter','sample']
		# parse bglm, gglm, bgee, ggee, fskato, gskato, bskato, fskat, gskat, bskat, fburden, gburden, or bburden formulas into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		fields = {}
		for x in [a for a in list(set(re_split('~|\+|\*|:|1|0|-1|factor|\(|\)',self.formula))) if a != '']:
			mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				fields[x] = {'class': 'factor', 'type': mtype}
			else:
				fields[x] = {'class': 'numeric', 'type': mtype}
		self.dict = fields
		for x in self.focus:
			self.header = self.header + [x + '.effect',x + '.stderr',x + '.or',x + '.z',x + '.p']
		self.header = self.header + ['n','status']

class GglmFrame(ModelFrame):
	def __init__(self, formula, type):
		super(GglmFrame, self).__init__(formula, type)
		self.header = ['chr','pos','a1','a2','marker','callrate','freq','freq.case','freq.ctrl','mac','rsq','hwe','filter','sample']
		# parse bglm, gglm, bgee, ggee, fskato, gskato, bskato, fskat, gskat, bskat, fburden, gburden, or bburden formulas into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		fields = {}
		for x in [a for a in list(set(re_split('~|\+|\*|:|1|0|-1|factor|\(|\)',self.formula))) if a != '']:
			mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				fields[x] = {'class': 'factor', 'type': mtype}
			else:
				fields[x] = {'class': 'numeric', 'type': mtype}
		self.dict = fields
		for x in self.focus:
			self.header = self.header + [x + '.effect',x + '.stderr',x + '.or',x + '.z',x + '.p']
		self.header = self.header + ['n','status']

class BgeeFrame(ModelFrame):
	def __init__(self, formula, type):
		super(BgeeFrame, self).__init__(formula, type)
		self.header = ['chr','pos','a1','a2','marker','callrate','freq','freq.case','freq.ctrl','mac','rsq','hwe','filter','sample']
		# parse bglm, gglm, bgee, ggee, fskato, gskato, bskato, fskat, gskat, bskat, fburden, gburden, or bburden formulas into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		fields = {}
		for x in [a for a in list(set(re_split('~|\+|\*|:|1|0|-1|factor|\(|\)',self.formula))) if a != '']:
			mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				fields[x] = {'class': 'factor', 'type': mtype}
			else:
				fields[x] = {'class': 'numeric', 'type': mtype}
		self.dict = fields
		for x in self.focus:
			self.header = self.header + [x + '.effect',x + '.stderr',x + '.or',x + '.z',x + '.p']
		self.header = self.header + ['n','status']

		from rpy2.robjects.packages import importr
		print "loading R package geepack"
		importr('geepack')

class GgeeFrame(ModelFrame):
	def __init__(self, formula, type):
		super(GgeeFrame, self).__init__(formula, type)
		self.header = ['chr','pos','a1','a2','marker','callrate','freq','freq.case','freq.ctrl','mac','rsq','hwe','filter','sample']
		# parse bglm, gglm, bgee, ggee, fskato, gskato, bskato, fskat, gskat, bskat, fburden, gburden, or bburden formulas into dictionary
		#
		# currently valid symbols for right hand side terms in these models
		# SYMBOL			EXAMPLE				ACTION
		# +					+x					estimates: the global effect of x
		# :					x:y					estimates: the global effect of the interaction between x and y
		# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
		# factor			factor(x)			specify x as a categorical variable (factor)
		fields = {}
		for x in [a for a in list(set(re_split('~|\+|\*|:|1|0|-1|factor|\(|\)',self.formula))) if a != '']:
			mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
			if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
				fields[x] = {'class': 'factor', 'type': mtype}
			else:
				fields[x] = {'class': 'numeric', 'type': mtype}
		self.dict = fields
		for x in self.focus:
			self.header = self.header + [x + '.effect',x + '.stderr',x + '.or',x + '.z',x + '.p']
		self.header = self.header + ['n','status']

		from rpy2.robjects.packages import importr
		print "loading R package geepack"
		importr('geepack')
"""



"""
elif self.type in ['bssmeta']:
			# parse bglm, gglm, bgee, ggee, fskato, gskato, bskato, fskat, gskat, bskat, fburden, gburden, or bburden formulas into dictionary
			#
			# currently valid symbols for right hand side terms in these models
			# SYMBOL			EXAMPLE				ACTION
			# +					+x					estimates: the global effect of x
			# :					x:y					estimates: the global effect of the interaction between x and y
			# *					x*y					estimates: the global effect of x, y, and the interaction between x and y
			# factor			factor(x)			specify x as a categorical variable (factor)
			fields = {}
			for x in [a for a in list(set(re_split('~|\+|\*|:|1|0|-1|factor|\(|\)',self.formula))) if a != '']:
				mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
				if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
					fields[x] = {'class': 'factor', 'type': mtype}
				else:
					fields[x] = {'class': 'numeric', 'type': mtype}
			self.dict = fields
			self.header = self.header + ['effect','stderr','or','z','p','nmiss','ntotal','status']
		elif self.type in ['blme','glme']:
			# parse blme or glme formulas into dictionary
			#
			# currently valid symbols for right hand side terms in these models, where x and y are covariates, o is some know offset, and g, g1, and g2 are grouping factors
			# EXAMPLE									ACTION
			# +x										estimates: the global effect of x
			# x:y										estimates: the global effect of the interaction between x and y
			# x*y										estimates: the global effect of x, y, and the interaction between x and y
			# (1|g), 1+(1|g)							estimates: a random intercept with fixed mean
			# 0+offset(o)+(1|g), -1+offset(o)+(1|g)		estimates: a random intercept with a priori means
			# (1|g1)+(1|g2), 1+(1|g1)+(1|g2)			estimates: an intercept varying among g1 and g2
			# (1|g1/g2), (1|g1)+(1|g1:g2)				estimates: an intercept varying among g1 and g2 within g1
			# x+(x|g), 1+x+(1+x|g) 						estimates: a correlated random intercept and slope
			# x+(x||g), 1+x+(1|g)+(0+x|g)				estimates: an uncorrelated random intercept and slope
			# factor			factor(x)				specify x as a categorical variable (factor)
			fields = {}
			for x in [a for a in list(set(re_split('\(|\)|~|\+|1|0|-1|\||\*|:|/|factor',self.formula))) if a != '']:
				mtype = "dependent" if x == re_split('~',self.formula)[0] else "independent"
				if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
					fields[x] = {'class': 'factor', 'type': mtype}
				else:
					fields[x] = {'class': 'numeric', 'type': mtype}
			self.dict = fields
		elif self.type == 'coxph':
			# parse coxph formulas into dictionary
			#
			# currently valid symbols for right hand side terms in these models, where x and y are covariates, t1 is the follow up time variable (for right censored data) or the starting
			# time of the interval (for interval censored data), t2 is the ending time of the interval (for interval censored data), e is the status indicator (eg. case/ctrl) variable
			# 
			# EXAMPLE									ACTION
			# +x										estimates: the global effect of x
			# x:y										estimates: the global effect of the interaction between x and y
			# x*y										estimates: the global effect of x, y, and the interaction between x and y
			# factor(x)									specify x as a categorical variable (factor)
			# cluster(factor(x))						specify x as a cluster variable (returns robust sandwich variance estimators)
			# frailty(factor(x))						specify x as a frailty variable (expects x to be a factor, adds a simple random effect for x)
			#
			# currently valid symbols for left hand side terms in these models
			# Surv(t1,e)								a survival object indicating a time (t1) and status (e) variable for right censored data
			# Surv(t1,t2,e)								a survival object indicating a start time (t1), an end time (t2), and a status (e) variable for interval censored data
			fields = {}
			for x in [a for a in list(set(re_split('Surv|\(|\)|,|~|\+|1|0|-1|\*|:|factor|cluster|frailty',self.formula))) if a != '']:
				surv = re_split('Surv|\(|\)|,',re_split('~',self.formula)[0])
				if x in surv:
					if len(surv) == 2:
						if x == surv[0]:
							mtype = "time1"
						else:
							mtype = "event"
					else:
						if x == surv[0]:
							mtype = "time1"
						elif x == surv[1]:
							mtype = "time2"
						else:
							mtype = "event"
				else:
					mtype = "independent"
				if self.formula[self.formula.find(x)-7:self.formula.find(x)] == 'factor(':
					fields[x] = {'class': 'factor', 'type': mtype}
				else:
					fields[x] = {'class': 'numeric', 'type': mtype}
			self.dict = fields
		else:
			raise Error("model type " + self.type + " unavailable")
			return
"""

