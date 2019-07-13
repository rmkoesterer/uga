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

from itertools import islice, chain
from re import sub as re_sub
import pysam
import pandas as pd
import time
import numpy as np
cimport numpy as np
import numpy.lib.recfunctions as recfxns
cimport cython
import Variant
cimport Variant
from cython.view cimport array
from collections import OrderedDict
import Process
import logging
import resource
import os

module_logger = logging.getLogger("Geno")
ILLEGAL_CHARS = '!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~'
MAX_CELL_LEN = 1000

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef double[:] vcf2dose(np.ndarray genos, hom1, het, hom2, np.int gt):
	cdef double[:] dose = np.ndarray(len(genos))
	for i in xrange(len(genos)):
		geno = genos[i].split(':')[gt]
		dose[i] = float('nan')
		if geno in hom1:
			dose[i] = 2.0
		elif geno in het:
			dose[i] = 1.0
		elif geno in hom2:
			dose[i] = 0.0
	return dose

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef double[:] extract_vcf_dose(np.ndarray genos, np.int field):
	cdef double[:] dose = np.ndarray(len(genos))
	for i in xrange(len(genos)):
		try:
			dose[i] = 2.0 - float(genos[i].split(':')[field])
		except:
			dose[i] = float('nan')
	return dose

cdef class Variants(object):

	def __cinit__(self, filename, sample_filename = None, chr = None, pos = None, id = None, a1 = None, a2 = None):
		logger = logging.getLogger("Geno.Variants.__cinit__")
		logger.debug("initialize variants")
		self.filename = filename
		self.sample_filename = sample_filename
		self.snvgroup_map = None

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef align(self, Variant.Ref ref):
		logger = logging.getLogger("Geno.Variants.align")
		logger.debug("align")
		cdef unsigned int i = 0
		for row in self.info:
			i += 1
			if ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a1'] + row['a2']:
				self.info['a1'][i-1] = ref.db[row['___uid___']]['a1']
				self.info['a2'][i-1] = ref.db[row['___uid___']]['a2']
				self.info['id'][i-1] = ref.db[row['___uid___']]['id']
				self.info['id_unique'][i-1] = ref.db[row['___uid___']]['id_unique']
			elif (((ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a2'][0]) + Variant.complement(row['a1'][0]) or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + row['a1'] or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a2'][0]) + 'NA' or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + 'NA') and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "AT" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "TA" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "GC" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "CG") or 
					((ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "AT" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "TA" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "GC" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "CG") and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + row['a1'])):
				self.info['a1'][i-1] = ref.db[row['___uid___']]['a1']
				self.info['a2'][i-1] = ref.db[row['___uid___']]['a2']
				self.info['id'][i-1] = ref.db[row['___uid___']]['id']
				self.info['id_unique'][i-1] = ref.db[row['___uid___']]['id_unique']
				self.data[:,i] = 2.0 - self.data[:,i].astype(float)

	def load_snvgroup_map(self, snvgroup_map):
		logger = logging.getLogger("Geno.Variants.load_snvgroup_map")
		logger.debug("load_snvgroup_map")
		print "loading snvgroup map " + os.path.basename(snvgroup_map)
		try:
			self.snvgroup_map=pd.read_table(snvgroup_map,names=['chr','pos','id','group_id'], compression='gzip' if snvgroup_map.split('.')[-1] == 'gz' else None)
		except:
			raise Process.Error("failed to load snvgroup map " + os.path.basename(snvgroup_map))

cdef class Vcf(Variants):
	def __cinit__(self, filename, sample_filename = None):
		logger = logging.getLogger("Geno.Vcf.__cinit__")
		logger.debug("initialize vcf")
		super(Vcf, self).__init__(filename)

		print "loading vcf file " + os.path.basename(filename)
		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asVCF())
		except:
			raise Process.Error("failed to load vcf file " + os.path.basename(filename))
		else:
			self.samples = np.array([a for a in self.handle.header][-1].split('\t')[9:])

	def get_region(self, region, group_id = None):
		logger = logging.getLogger("Geno.Vcf.get_region")
		logger.debug("get_region " + region)
		self.group_id = group_id
		if ':' in region:
			self.chr = int(region.split(':')[0])
			self.start = int(region.split(':')[1].split('-')[0])
			self.end = int(region.split(':')[1].split('-')[1])
		else:
			self.chr = int(region.split(':')[0])
			self.start = 0
			self.end = 1000000000
		try:
			self.region_iter = self.handle.fetch(region=region, parser=pysam.asTuple())
		except:
			raise

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_chunk(self, int buffer):
		logger = logging.getLogger("Geno.Vcf.get_chunk")
		logger.debug("get_chunk")
		self.snv_chunk = np.empty((buffer,11 + len(self.samples)), dtype='object')
		cdef unsigned int i = 0
		slice = islice(self.region_iter, buffer)
		try:
			first = slice.next()
		except StopIteration:
			raise
		else:
			slice = chain([first], slice)
			for r in slice:
				record = np.array(r, dtype='|S' + str(MAX_CELL_LEN))
				fields = record[8].split(':')
				if 'DS' in fields:
					field = fields.index('DS')
					self.snv_chunk[i,:9] = record[:9]
					self.snv_chunk[i,11:] = extract_vcf_dose(record[9:], field)
					self.snv_chunk[i,5] = self.group_id
					self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:1000])
					self.snv_chunk[i,10] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
					i += 1
				elif 'GT' in fields:
					gt = fields.index('GT')
					alts = record[4].split(',')
					if len(alts) > 1:
						self.snv_chunk = np.append(self.snv_chunk, np.empty((2*len(alts),11 + len(self.samples)), dtype='object'),axis=0)
						for alt in xrange(len(alts)):
							oth = [alts.index(a) for a in alts if a != alts[alt]]
							het = het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']] + [str(alts.index(k)+1) + sep + str(alt+1) for k in alts if k != alts[alt] for sep in ['/','|']] + [str(alt+1) + sep + str(alts.index(k)+1) for k in alts if k != alts[alt] for sep in ['/','|']]))
							hom1 = list(set(['0/0'] + [str(a+1) + '/' + str(a+1) for a in oth] + ['0|0'] + [str(a+1) + '|' + str(a+1) for a in oth]))
							hom2 = list(set([str(alt+1) + '/' + str(alt+1)] + [str(alt+1) + '|' + str(alt+1)]))
							self.snv_chunk[i,:9] = record[:9]
							self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
							self.snv_chunk[i,3] = '&'.join([record[3]] + [alts[a] for a in oth])
							self.snv_chunk[i,4] = alts[alt]
							self.snv_chunk[i,5] = self.group_id
							self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:1000])
							self.snv_chunk[i,10] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
							het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
							hom1 = ['0/0','0|0']
							hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
							self.snv_chunk[i+1,:9] = record[:9]
							self.snv_chunk[i+1,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
							self.snv_chunk[i+1,3] = record[3]
							self.snv_chunk[i+1,4] = alts[alt]
							self.snv_chunk[i+1,5] = self.group_id
							self.snv_chunk[i+1,9] = 'chr' + self.snv_chunk[i+1,0] + 'bp' + self.snv_chunk[i+1,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,4][0:1000])
							self.snv_chunk[i+1,10] = Variant.get_universal_variant_id(self.snv_chunk[i+1,0],self.snv_chunk[i+1,1],self.snv_chunk[i+1,3],self.snv_chunk[i+1,4],'><')
							i += 2
					else:
						het = ['1/0','0/1','1|0','0|1']
						hom1 = ['0/0','0|0']
						hom2 = ['1/1','1|1']
						self.snv_chunk[i,:9] = record[:9]
						self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
						self.snv_chunk[i,5] = self.group_id
						self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:1000]) + '.' + str(i)
						self.snv_chunk[i,10] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
						i += 1
				else:
					raise Process.Error("failed to load vcf file, GT or DS field required")
			self.snv_chunk = self.snv_chunk[np.where((self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end))]

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Geno.Vcf.get_snvs")
		logger.debug("get_snvs")
		try:
			self.get_chunk(buffer)
		except:
			raise
		self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','id','a1','a2','group_id','id_unique','___uid___']),np.array(['uint8','uint32','|S60','|S1000','|S1000','|S1000','|S1000','|S1000'])))
		self.data = np.column_stack((self.samples, self.snv_chunk[:,11:].transpose()))
		var_ids, var_ids_idx, var_ids_cnt = np.unique(['.'.join(x.split('.')[0:len(x.split('.'))]) for x in self.info['id_unique']], return_inverse=True, return_counts=True)
		self.duplicated = var_ids[var_ids_cnt > 1]
		np.place(self.data, self.data == 'NA',np.nan)

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvgroup(self, int buffer, group_id):
		logger = logging.getLogger("Geno.Vcf.get_snvgroup")
		logger.debug("get_snvgroup " + group_id)
		cdef unsigned int i = 0
		self.snvgroup_chunk = np.empty((0,11 + len(self.samples)), dtype='object')
		if self.snvgroup_map is not None:
			snvgroup_snvs = list(self.snvgroup_map['id'][self.snvgroup_map['group_id'] == group_id])
		while True:
			try:
				self.get_chunk(buffer)
			except:
				break
			i = i + 1
			if self.snvgroup_map is not None:
				self.snv_chunk = self.snv_chunk[np.where(np.in1d(self.snv_chunk[:,2],snvgroup_snvs))]
			self.snvgroup_chunk = np.vstack((self.snvgroup_chunk,self.snv_chunk))
		if i == 0:
			raise
		else:
			self.info = np.array([tuple(row) for row in self.snvgroup_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','id','a1','a2','group_id','id_unique','___uid___']),np.array(['uint8','uint32','|S60','|S1000','|S1000','|S1000','|S1000','|S1000'])))
			self.data = np.column_stack((self.samples, self.snvgroup_chunk[:,11:].transpose()))
			np.place(self.data, self.data == 'NA',np.nan)

cdef class Dos(Variants):
	def __cinit__(self, filename, sample_filename):
		logger = logging.getLogger("Geno.Dos.__cinit__")
		logger.debug("initialize dos")
		super(Dos, self).__init__(filename, sample_filename)

		print "loading dos file " + os.path.basename(filename)
		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asTuple())
		except:
			raise Process.Error("failed to load dos file " + os.path.basename(filename))
		else:
			print "loading dos sample file " + os.path.basename(sample_filename)
			try:
				self.samples=np.genfromtxt(fname=self.sample_filename, dtype='object')
			except:
				raise Process.Error("failed to load dos sample file " + os.path.basename(sample_filename))

	def get_region(self, region, group_id = None):
		logger = logging.getLogger("Geno.Dos.get_region")
		logger.debug("get_region " + region)
		self.group_id = group_id
		if ':' in region:
			self.chr = int(region.split(':')[0])
			self.start = int(region.split(':')[1].split('-')[0])
			self.end = int(region.split(':')[1].split('-')[1])
		else:
			self.chr = int(region.split(':')[0])
			self.start = 0
			self.end = 1000000000
		try:
			self.region_iter = self.handle.fetch(region=region, parser=pysam.asTuple())
		except:
			raise

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_chunk(self, int buffer):
		logger = logging.getLogger("Geno.Dos.get_chunk")
		logger.debug("get_chunk")
		self.snv_chunk = np.empty((buffer,8 + len(self.samples)), dtype='object')
		cdef unsigned int i = 0
		slice = islice(self.region_iter, buffer)
		try:
			first = slice.next()
		except StopIteration:
			raise
		else:
			slice = chain([first], slice)
			for r in slice:
				record = np.array(r, dtype='|S' + str(MAX_CELL_LEN))
				self.snv_chunk[i,:5] = record[:5]
				self.snv_chunk[i,8:] = record[5:]
				self.snv_chunk[i,5] = self.group_id
				self.snv_chunk[i,6] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:1000]) + '.' + str(i)
				self.snv_chunk[i,7] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
				i += 1
			self.snv_chunk = self.snv_chunk[np.where((self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end))]

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Geno.Dos.get_snvs")
		logger.debug("get_snvs")
		try:
			self.get_chunk(buffer)
		except:
			raise
		self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,6,7]]], dtype=zip(np.array(['chr','pos','id','a1','a2','group_id','id_unique','___uid___']),np.array(['uint8','uint32','|S60','|S1000','|S1000','|S1000','|S1000','|S1000'])))
		self.data = np.column_stack((self.samples, self.snv_chunk[:,8:].transpose()))
		var_ids, var_ids_idx, var_ids_cnt = np.unique(['.'.join(x.split('.')[0:len(x.split('.'))]) for x in self.info['id_unique']], return_inverse=True, return_counts=True)
		self.duplicated = var_ids[var_ids_cnt > 1]
		np.place(self.data, self.data == 'NA',np.nan)

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvgroup(self, int buffer, group_id):
		logger = logging.getLogger("Geno.Dos.get_snvgroup")
		logger.debug("get_snvgroup " + group_id)
		cdef unsigned int i = 0
		self.snvgroup_chunk = np.empty((0,8 + len(self.samples)), dtype='object')
		if self.snvgroup_map is not None:
			snvgroup_snvs = list(self.snvgroup_map['id'][self.snvgroup_map['group_id'] == group_id])
		while True:
			try:
				self.get_chunk(buffer)
			except:
				break
			i = i + 1
			if self.snvgroup_map is not None:
				self.snv_chunk = self.snv_chunk[np.where(np.in1d(self.snv_chunk[:,2],snvgroup_snvs))]
			self.snvgroup_chunk = np.vstack((self.snvgroup_chunk,self.snv_chunk))
		if i == 0:
			raise
		else:
			self.info = np.array([tuple(row) for row in self.snvgroup_chunk[:,[0,1,2,3,4,5,6,7]]], dtype=zip(np.array(['chr','pos','id','a1','a2','group_id','id_unique','___uid___']),np.array(['uint8','uint32','|S60','|S1000','|S1000','|S1000','|S1000','|S1000'])))
			self.data = np.column_stack((self.samples, self.snvgroup_chunk[:,8:].transpose()))
			np.place(self.data, self.data == 'NA',np.nan)

cdef class Results(Variants):
	def __cinit__(self, filename, chr, pos, id, a1, a2):
		logger = logging.getLogger("Geno.Results.__cinit__")
		logger.debug("initialize results")
		super(Results, self).__init__(filename, chr, pos, id, a1, a2)

		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asTuple())
		except:
			raise Process.Error("failed to load results file")
		else:
			self.header = [x for x in self.handle.header]
			#self.method = self.header[2].split()[2]
			colsdict = OrderedDict({v:k for k, v in {'chr': chr, 'pos': pos, 'id': id, 'a1': a1, 'a2': a2}.iteritems()})
			self.cols = [x for x in map(colsdict.get, self.header[-1].split(), self.header[-1].split())]
			self.cols = ['chr','pos','id','a1','a2'] + [x for x in self.cols if x not in ['chr','pos','id','a1','a2']]
			self.dtypes = zip([x for x in self.cols if x in ['chr','pos','id','a1','a2']],['uint8','uint32','|S60','|S1000','|S1000']) + zip([x for x in self.cols if x not in ['chr','pos','id','a1','a2']],['|S1000' if x in ['dir','test'] else 'f8' for x in self.cols if x not in ['chr','pos','id','a1','a2']]) + [('id_unique','|S1000'),('___uid___','|S1000')]

	def get_region(self, region):
		logger = logging.getLogger("Geno.Results.get_region")
		logger.debug("get_region " + region)
		if ':' in region:
			self.chr = int(region.split(':')[0])
			self.start = int(region.split(':')[1].split('-')[0])
			self.end = int(region.split(':')[1].split('-')[1])
		else:
			self.chr = int(region.split(':')[0])
			self.start = 0
			self.end = 1000000000
		try:
			self.region_iter = self.handle.fetch(region=region, parser=pysam.asTuple())
		except:
			raise

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_chunk(self, int buffer):
		logger = logging.getLogger("Geno.Results.get_chunk")
		logger.debug("get_chunk")
		self.snv_chunk = np.empty((buffer,len(self.dtypes)), dtype='object')
		cdef unsigned int i = 0
		slice = islice(self.region_iter, buffer)
		try:
			first = slice.next()
		except StopIteration:
			raise
		else:
			slice = chain([first], slice)
			for r in slice:
				record = np.array(r, dtype='object')
				self.snv_chunk[i,:len(self.cols)] = record[:len(self.cols)]
				self.snv_chunk[i,len(self.cols)] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:1000]) + '.'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:1000])
				self.snv_chunk[i,len(self.cols)+1] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
				i += 1
			self.snv_chunk = self.snv_chunk[np.where((self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end))]

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Geno.Results.get_snvs")
		logger.debug("get_snvs")
		cdef unsigned int i = 0
		self.snv_results = np.empty((0,len(self.dtypes)), dtype='object')
		while True:
			try:
				self.get_chunk(buffer)
			except:
				break
			i = i + 1
			self.snv_results = np.vstack((self.snv_results,self.snv_chunk))
		if i == 0:
			raise
		else:
			self.snv_results[self.snv_results == "NA"]=np.nan
			self.snv_results = np.array([tuple(row) for row in self.snv_results], dtype=self.dtypes)
			var_ids, var_ids_idx, var_ids_cnt = np.unique(self.snv_results['id_unique'], return_inverse=True, return_counts=True)
			self.duplicated = var_ids[var_ids_cnt > 1]
			np.place(self.snv_results, self.snv_results == 'NA',np.nan)

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef align_results(self, Variant.Ref ref):
		logger = logging.getLogger("Geno.Results.align_results")
		logger.debug("align_results")
		cdef unsigned int i = 0
		for row in self.snv_results:
			i += 1
			if (ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a1'] + row['a2'] or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a1'] + 'NA' or
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == 'NA' + row['a2'] or
					(ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "AT" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "TA" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "GC" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "CG" and 
					(ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a1'][0]) + Variant.complement(row['a2'][0]) or	
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a1'][0]) + 'NA' or
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == 'NA' + Variant.complement(row['a2'][0])))):
				self.snv_results['a1'][i-1] = ref.db[row['___uid___']]['a1']
				self.snv_results['a2'][i-1] = ref.db[row['___uid___']]['a2']
				self.snv_results['id'][i-1] = ref.db[row['___uid___']]['id']
				self.snv_results['id_unique'][i-1] = ref.db[row['___uid___']]['id_unique']
			elif (ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + row['a1'] or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + 'NA' or
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == 'NA' + row['a1'] or
					((ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a2'][0]) + Variant.complement(row['a1'][0]) or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == Variant.complement(row['a2'][0]) + 'NA' or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == 'NA' + Variant.complement(row['a1'][0])) and
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "AT" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "TA" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "GC" and 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] != "CG") or
					(ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == row['a2'] + row['a1'] and
					(ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "AT" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "TA" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "GC" or 
					ref.db[row['___uid___']]['a1'] + ref.db[row['___uid___']]['a2'] == "CG"))):
				self.snv_results['a1'][i-1] = ref.db[row['___uid___']]['a1']
				self.snv_results['a2'][i-1] = ref.db[row['___uid___']]['a2']
				self.snv_results['id'][i-1] = ref.db[row['___uid___']]['id']
				self.snv_results['id_unique'][i-1] = ref.db[row['___uid___']]['id_unique']
				if 'freq' in self.snv_results.dtype.names:
					self.snv_results['freq'][i-1] = 1.0 - self.snv_results['freq'][i-1]
				if 'freq.case' in self.snv_results.dtype.names:
					self.snv_results['freq.case'][i-1] = 1.0 - self.snv_results['freq.case'][i-1]
				if 'freq.ctrl' in self.snv_results.dtype.names:
					self.snv_results['freq.ctrl'][i-1] = 1.0 - self.snv_results['freq.ctrl'][i-1]
				if 'effect' in self.snv_results.dtype.names:
					self.snv_results['effect'][i-1] = -1.0 * self.snv_results['effect'][i-1]
				if 'or' in self.snv_results.dtype.names:
					self.snv_results['or'][i-1] = 1.0 / self.snv_results['or'][i-1]
				if 'z' in self.snv_results.dtype.names:
					self.snv_results['z'][i-1] = -1.0 * self.snv_results['z'][i-1]
				if 't' in self.snv_results.dtype.names:
					self.snv_results['t'][i-1] = -1.0 * self.snv_results['t'][i-1]
				if 'dir' in self.snv_results.dtype.names:
					self.snv_results['dir'][i-1] = self.snv_results['dir'][i-1].replace('+','!').replace('-','+').replace('!','-')

	def tag_results(self, tag):
		logger = logging.getLogger("Geno.Results.get_tagged_results")
		logger.debug("get_tagged_results")
		self.snv_results_tagged = pd.to_numeric(pd.DataFrame(self.snv_results),errors='coerce')
		self.snv_results_tagged.columns = [tag + '.' + x[0] if x[0] not in ['chr','pos','id','a1','a2','id_unique','___uid___'] else x[0] for x in self.dtypes]
