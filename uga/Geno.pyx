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
cimport cython
import Variant
cimport Variant
from cython.view cimport array
from collections import OrderedDict
import Process
import logging

module_logger = logging.getLogger("Geno")
ILLEGAL_CHARS = '!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~'

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
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

cdef class Variants(object):

	def __cinit__(self, filename, sample_filename = None):
		logger = logging.getLogger("Geno.Variants.__cinit__")
		logger.debug("initialize variants")
		self.filename = filename
		self.sample_filename = sample_filename
		self.snvgroup_map = None

	cpdef align(self, Variant.Ref ref):
		logger = logging.getLogger("Geno.Variants.align")
		logger.debug("align")
		cdef unsigned int i = 0
		for row in self.info:
			i += 1
			if ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == row['a1'] + row['a2']:
				self.info['a1'][i-1] = ref.db[row['uid']]['a1']
				self.info['a2'][i-1] = ref.db[row['uid']]['a2']
				self.info['variant'][i-1] = ref.db[row['uid']]['variant']
				self.info['variant_unique'][i-1] = ref.db[row['uid']]['variant_unique']
			if ((ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == Variant.complement(row['a2'][0]) + Variant.complement(row['a1'][0]) or 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == row['a2'] + row['a1'] or 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == Variant.complement(row['a2'][0]) + 'NA' or 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == row['a2'] + 'NA') and 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] != "AT" and 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] != "TA" and 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] != "GC" and 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] != "CG"):
				self.info['a1'][i-1] = ref.db[row['uid']]['a1']
				self.info['a2'][i-1] = ref.db[row['uid']]['a2']
				self.info['variant'][i-1] = ref.db[row['uid']]['variant']
				self.info['variant_unique'][i-1] = ref.db[row['uid']]['variant_unique']
				self.data[:,i] = 2.0 - self.data[:,i].astype(float)

	def load_snvgroup_map(self, snvgroup_map):
		logger = logging.getLogger("Geno.Variants.load_snvgroup_map")
		logger.debug("load_snvgroup_map")
		print "loading snvgroup map"
		try:
			self.snvgroup_map=pd.read_table(snvgroup_map,names=['chr','pos','variant','id'], compression='gzip' if snvgroup_map.split('.')[-1] == 'gz' else None)
		except:
			raise Process.Error("failed to load snvgroup map")

cdef class Vcf(Variants):
	def __cinit__(self, filename, sample_filename = None):
		logger = logging.getLogger("Geno.Vcf.__cinit__")
		logger.debug("initialize vcf")
		super(Vcf, self).__init__(filename)

		print "loading vcf file"
		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asVCF())
		except:
			raise Process.Error("failed to load vcf file")
		else:
			self.samples = np.array([a for a in self.handle.header][-1].split('\t')[9:])

	def get_region(self, region, id = None):
		logger = logging.getLogger("Geno.Vcf.get_region")
		logger.debug("get_region " + region)
		self.id = id
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
				record = np.array(r)
				gt = record[8].split(':').index('GT')
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
						self.snv_chunk[i,5] = self.id
						self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:20]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:20])
						self.snv_chunk[i,10] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
						het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
						hom1 = ['0/0','0|0']
						hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
						self.snv_chunk[i+1,:9] = record[:9]
						self.snv_chunk[i+1,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
						self.snv_chunk[i+1,3] = record[3]
						self.snv_chunk[i+1,4] = alts[alt]
						self.snv_chunk[i+1,5] = self.id
						self.snv_chunk[i+1,9] = 'chr' + self.snv_chunk[i+1,0] + 'bp' + self.snv_chunk[i+1,1] + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,2][0:60]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,3][0:20]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i+1,4][0:20])
						self.snv_chunk[i+1,10] = Variant.get_universal_variant_id(self.snv_chunk[i+1,0],self.snv_chunk[i+1,1],self.snv_chunk[i+1,3],self.snv_chunk[i+1,4],'><')
						i += 2
				else:
					het = ['1/0','0/1','1|0','0|1']
					hom1 = ['0/0','0|0']
					hom2 = ['1/1','1|1']
					self.snv_chunk[i,:9] = record[:9]
					self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
					self.snv_chunk[i,5] = self.id
					self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:20]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:20])
					self.snv_chunk[i,10] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
					i += 1
			self.snv_chunk = self.snv_chunk[(self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end)]

	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Geno.Vcf.get_snvs")
		logger.debug("get_snvs")
		try:
			self.get_chunk(buffer)
		except:
			raise
		self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
		self.data = np.column_stack((self.samples, self.snv_chunk[:,11:].transpose()))
		np.place(self.data, self.data == 'NA',np.nan)

	cpdef get_snvgroup(self, int buffer, id):
		logger = logging.getLogger("Geno.Vcf.get_snvgroup")
		logger.debug("get_snvgroup " + id)
		cdef unsigned int i = 0
		self.snvgroup_chunk = np.empty((0,11 + len(self.samples)), dtype='object')
		if self.snvgroup_map is not None:
			snvgroup_snvs = list(self.snvgroup_map['variant'][self.snvgroup_map['id'] == id])
		while True:
			try:
				self.get_chunk(buffer)
			except:
				break
			i = i + 1
			self.snv_chunk = self.snv_chunk[np.where(np.in1d(self.snv_chunk[:,2],snvgroup_snvs))]
			if self.snvgroup_map is not None:
				self.snvgroup_chunk = np.vstack((self.snvgroup_chunk,self.snv_chunk))
		if i == 0:
			raise
		else:
			self.info = np.array([tuple(row) for row in self.snvgroup_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
			self.data = np.column_stack((self.samples, self.snvgroup_chunk[:,11:].transpose()))
			np.place(self.data, self.data == 'NA',np.nan)

cdef class Dos(Variants):
	def __cinit__(self, filename, sample_filename):
		logger = logging.getLogger("Geno.Dos.__cinit__")
		logger.debug("initialize dos")
		super(Dos, self).__init__(filename, sample_filename)

		print "loading dos file"
		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asTuple())
		except:
			raise Process.Error("failed to load dos file")
		else:
			print "loading dos sample file"
			try:
				self.samples=np.genfromtxt(fname=self.sample_filename, dtype='object')
			except:
				raise Process.Error("failed to load dos sample file")

	def get_region(self, region, id = None):
		logger = logging.getLogger("Geno.Dos.get_region")
		logger.debug("get_region " + region)
		self.id = id
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
				record = np.array(r)
				self.snv_chunk[i,:5] = record[:5]
				self.snv_chunk[i,8:] = record[5:]
				self.snv_chunk[i,5] = self.id
				self.snv_chunk[i,6] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,3][0:20]) + '_'  + re_sub(ILLEGAL_CHARS,'_',self.snv_chunk[i,4][0:20])
				self.snv_chunk[i,7] = Variant.get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
				i += 1
			self.snv_chunk = self.snv_chunk[(self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end)]

	cpdef get_snvs(self, int buffer):
		logger = logging.getLogger("Geno.Dos.get_snvs")
		logger.debug("get_snvs")
		try:
			self.get_chunk(buffer)
		except:
			raise
		self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,6,7]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
		self.data = np.column_stack((self.samples, self.snv_chunk[:,8:].transpose()))
		np.place(self.data, self.data == 'NA',np.nan)

	cpdef get_snvgroup(self, int buffer, id):
		logger = logging.getLogger("Geno.Dos.get_snvgroup")
		logger.debug("get_snvgroup " + id)
		cdef unsigned int i = 0
		self.snvgroup_chunk = np.empty((0,8 + len(self.samples)), dtype='object')
		if self.snvgroup_map is not None:
			snvgroup_snvs = list(self.snvgroup_map['variant'][self.snvgroup_map['id'] == id])
		while True:
			try:
				self.get_chunk(buffer)
			except:
				break
			i = i + 1
			self.snv_chunk = self.snv_chunk[np.where(np.in1d(self.snv_chunk[:,2],snvgroup_snvs))]
			if self.snvgroup_map is not None:
				self.snvgroup_chunk = np.vstack((self.snvgroup_chunk,self.snv_chunk))
		if i == 0:
			raise
		else:
			self.info = np.array([tuple(row) for row in self.snvgroup_chunk[:,[0,1,2,3,4,5,6,7]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
			self.data = np.column_stack((self.samples, self.snvgroup_chunk[:,8:].transpose()))
			np.place(self.data, self.data == 'NA',np.nan)
