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

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.cdivision(True)


from itertools import islice, chain
from re import sub as re_sub
import pysam
import pandas as pd
import time
import numpy as np
cimport numpy as np
cimport cython
from cython.view cimport array
from collections import OrderedDict
from Process import Error

@cython.boundscheck(False)
@cython.wraparound(False)
def complement(allele):
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
def get_universal_variant_id(chr_py,pos_py,a1_py,a2_py,delim_py):
	cdef str chr = chr_py
	cdef str pos = pos_py
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
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
	cdef public bytes filename, region, id
	cdef public np.ndarray samples
	cdef public object handle, region_iter, snvgroup_map
	cdef public unsigned int chr, start, end
	cdef public np.ndarray genos, data, info, snv_chunk, snvgroup_chunk
	def __cinit__(self, filename):
		self.filename = filename
		self.snvgroup_map = None

	cpdef align(self, VariantRef ref):
		cdef unsigned int i = 0
		for row in self.info:
			i += 1
			if ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == row['a1'] + row['a2']:
				self.info['a1'][i-1] = ref.db[row['uid']]['a1']
				self.info['a2'][i-1] = ref.db[row['uid']]['a2']
				self.info['variant'][i-1] = ref.db[row['uid']]['variant']
				self.info['variant_unique'][i-1] = ref.db[row['uid']]['variant_unique']
			if ((ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == complement(row['a2'][0]) + complement(row['a1'][0]) or 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == row['a2'] + row['a1'] or 
					ref.db[row['uid']]['a1'] + ref.db[row['uid']]['a2'] == complement(row['a2'][0]) + 'NA' or 
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
		print "loading snvgroup map"
		try:
			self.snvgroup_map=pd.read_table(snvgroup_map,names=['chr','pos','variant','id'])
		except:
			raise Error("failed to load snvgroup map")

cdef class Vcf(Variants):
	def __cinit__(self, filename):
		super(Vcf, self).__init__(filename)

		print "loading vcf file"
		try:
			self.handle=pysam.TabixFile(filename=filename,parser=pysam.asVCF())
		except:
			raise Error("failed to load vcf file")
		else:
			self.samples = np.array([a for a in self.handle.header][-1].split('\t')[9:])

	def get_region(self, region, id = None):
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
						self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
						self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
						het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
						hom1 = ['0/0','0|0']
						hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
						self.snv_chunk[i+1,:9] = record[:9]
						self.snv_chunk[i+1,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
						self.snv_chunk[i+1,3] = record[3]
						self.snv_chunk[i+1,4] = alts[alt]
						self.snv_chunk[i+1,5] = self.id
						self.snv_chunk[i+1,9] = 'chr' + self.snv_chunk[i+1,0] + 'bp' + self.snv_chunk[i+1,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i+1,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i+1,3][0:20]) + '_'  + self.snv_chunk[i+1,4][0:20]
						self.snv_chunk[i+1,10] = get_universal_variant_id(self.snv_chunk[i+1,0],self.snv_chunk[i+1,1],self.snv_chunk[i+1,3],self.snv_chunk[i+1,4],'><')
						i += 2
				else:
					het = ['1/0','0/1','1|0','0|1']
					hom1 = ['0/0','0|0']
					hom2 = ['1/1','1|1']
					self.snv_chunk[i,:9] = record[:9]
					self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
					self.snv_chunk[i,5] = self.id
					self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
					self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
					i += 1
			self.snv_chunk = self.snv_chunk[(self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end)]

	cpdef get_snvs(self, int buffer):
		try:
			self.get_chunk(buffer)
		except:
			raise
		self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
		self.data = np.column_stack((self.samples, self.snv_chunk[:,11:].transpose()))
		np.place(self.data, self.data == 'NA',np.nan)

	cpdef get_snvgroup(self, int buffer, id):
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


		#while True:
		#	i = 0
		#	c = np.empty((buffer,11 + len(self.samples)), dtype='object')
		#	slice = islice(self.region_iter, buffer)
		#	try:
		#		first = slice.next()
		#		print first
		#	except StopIteration:
		#		print "HERE"
		#		raise
		#	else:
		#		slice = chain([first], slice)
		#		for r in slice:
		#			record = np.array(r)
		#			if record[2] in gene_snvs:
		#				gt = record[8].split(':').index('GT')
		#				alts = record[4].split(',')
		#				if len(alts) > 1:
		#					c = np.append(c, np.empty((2*len(alts),10 + len(self.samples)), dtype='object'),axis=0)
		#					for alt in xrange(len(alts)):
		#						oth = [alts.index(a) for a in alts if a != alts[alt]]
		#						het = het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']] + [str(alts.index(k)+1) + sep + str(alt+1) for k in alts if k != alts[alt] for sep in ['/','|']] + [str(alt+1) + sep + str(alts.index(k)+1) for k in alts if k != alts[alt] for sep in ['/','|']]))
		#						hom1 = list(set(['0/0'] + [str(a+1) + '/' + str(a+1) for a in oth] + ['0|0'] + [str(a+1) + '|' + str(a+1) for a in oth]))
		#						hom2 = list(set([str(alt+1) + '/' + str(alt+1)] + [str(alt+1) + '|' + str(alt+1)]))
		#						self.snv_chunk[i,:9] = record[:9]
		#						self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
		#						self.snv_chunk[i,3] = '&'.join([record[3]] + [alts[a] for a in oth])
		#						self.snv_chunk[i,4] = alts[alt]
		#						self.snv_chunk[i,5] = self.id
		#						self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
		#						self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
		#						het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
		#						hom1 = ['0/0','0|0']
		#						hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
		#						self.snv_chunk[i+1,:9] = record[:9]
		#						self.snv_chunk[i+1,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
		#						self.snv_chunk[i+1,3] = record[3]
		#						self.snv_chunk[i+1,4] = alts[alt]
		#						self.snv_chunk[i+1,5] = self.id
		#						self.snv_chunk[i+1,9] = 'chr' + self.snv_chunk[i+1,0] + 'bp' + self.snv_chunk[i+1,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i+1,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i+1,3][0:20]) + '_'  + self.snv_chunk[i+1,4][0:20]
		#						self.snv_chunk[i+1,10] = get_universal_variant_id(self.snv_chunk[i+1,0],self.snv_chunk[i+1,1],self.snv_chunk[i+1,3],self.snv_chunk[i+1,4],'><')
		#						i += 2
		#				else:
		#					het = ['1/0','0/1','1|0','0|1']
		#					hom1 = ['0/0','0|0']
		#					hom2 = ['1/1','1|1']
		#					self.snv_chunk[i,:9] = record[:9]
		#					self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
		#					self.snv_chunk[i,5] = self.id
		#					self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
		#					self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
		#					i += 1
		#		c = self.snv_chunk[(self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end)]
		#	c_compiled = np.vstack((c_compiled,c))
		#self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
		#self.data = np.column_stack((self.samples, self.snv_chunk[:,11:].transpose()))
		#np.place(self.data, self.data == 'NA',np.nan)
		

	#cpdef get_gene(self, gene):
	#	c = np.empty((100000,11 + len(self.samples)), dtype='object')
	#	gene_snvs = list(self.gene_map['variant'][self.gene_map['gene'] == gene])
	#	cdef unsigned int i = 0
	#	for r in self.region_iter:
	#		record = np.array(r)
	#		if record[2] in gene_snvs:
	#			gt = record[8].split(':').index('GT')
	#			alts = record[4].split(',')
	#			if len(alts) > 1:
	#				c = np.append(c, np.empty((2*len(alts),10 + len(self.samples)), dtype='object'),axis=0)
	#				for alt in xrange(len(alts)):
	#					oth = [alts.index(a) for a in alts if a != alts[alt]]
	#					het = het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']] + [str(alts.index(k)+1) + sep + str(alt+1) for k in alts if k != alts[alt] for sep in ['/','|']] + [str(alt+1) + sep + str(alts.index(k)+1) for k in alts if k != alts[alt] for sep in ['/','|']]))
	#					hom1 = list(set(['0/0'] + [str(a+1) + '/' + str(a+1) for a in oth] + ['0|0'] + [str(a+1) + '|' + str(a+1) for a in oth]))
	#					hom2 = list(set([str(alt+1) + '/' + str(alt+1)] + [str(alt+1) + '|' + str(alt+1)]))
	#					self.snv_chunk[i,:9] = record[:9]
	#					self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
	#					self.snv_chunk[i,3] = '&'.join([record[3]] + [alts[a] for a in oth])
	#					self.snv_chunk[i,4] = alts[alt]
	#					self.snv_chunk[i,5] = self.id
	#					self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
	#					self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
	#					het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
	#					hom1 = ['0/0','0|0']
	#					hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
	#					self.snv_chunk[i+1,:9] = record[:9]
	#					self.snv_chunk[i+1,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
	#					self.snv_chunk[i+1,3] = record[3]
	#					self.snv_chunk[i+1,4] = alts[alt]
	#					self.snv_chunk[i+1,5] = self.id
	#					self.snv_chunk[i+1,9] = 'chr' + self.snv_chunk[i+1,0] + 'bp' + self.snv_chunk[i+1,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i+1,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i+1,3][0:20]) + '_'  + self.snv_chunk[i+1,4][0:20]
	#					self.snv_chunk[i+1,10] = get_universal_variant_id(self.snv_chunk[i+1,0],self.snv_chunk[i+1,1],self.snv_chunk[i+1,3],self.snv_chunk[i+1,4],'><')
	#					i += 2
	#			else:
	#				het = ['1/0','0/1','1|0','0|1']
	#				hom1 = ['0/0','0|0']
	#				hom2 = ['1/1','1|1']
	#				self.snv_chunk[i,:9] = record[:9]
	#				self.snv_chunk[i,11:] = vcf2dose(record[9:], hom1, het, hom2, gt)
	#				self.snv_chunk[i,5] = self.id
	#				self.snv_chunk[i,9] = 'chr' + self.snv_chunk[i,0] + 'bp' + self.snv_chunk[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',self.snv_chunk[i,2][0:60]) + '_'  + re_sub('&','_',self.snv_chunk[i,3][0:20]) + '_'  + self.snv_chunk[i,4][0:20]
	#				self.snv_chunk[i,10] = get_universal_variant_id(self.snv_chunk[i,0],self.snv_chunk[i,1],self.snv_chunk[i,3],self.snv_chunk[i,4],'><')
	#				i += 1
	#	c = self.snv_chunk[(self.snv_chunk[:i,1].astype(int) >= self.start) & (self.snv_chunk[:i,1].astype(int) <= self.end)]
	#	self.info = np.array([tuple(row) for row in self.snv_chunk[:,[0,1,2,3,4,5,9,10]]], dtype=zip(np.array(['chr','pos','variant','a1','a2','id','variant_unique','uid']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100','|S1000'])))
	#	self.data = np.column_stack((self.samples, self.snv_chunk[:,11:].transpose()))
	#	np.place(self.data, self.data == 'NA',np.nan)

cdef class VariantRef(object):
	cdef public object db
	def __cinit__(self, Variants v):
		self.db = {}
		for row in v.info:
			self.db[row['uid']] = {}
			self.db[row['uid']]['variant'] = row['variant']
			self.db[row['uid']]['variant_unique'] = row['variant_unique']
			self.db[row['uid']]['a1'] = row['a1']
			self.db[row['uid']]['a2'] = row['a2']

	cpdef update(self, Variants v):
		for row in v.info:
			if not row['uid'] in self.db:
				self.db[row['uid']] = {}
				self.db[row['uid']]['variant'] = row['variant']
				self.db[row['uid']]['variant_unique'] = row['variant_unique']
				self.db[row['uid']]['a1'] = row['a1']
				self.db[row['uid']]['a2'] = row['a2']

"""
def ExtractPlink(data_handle, reg, snv_list = None):
	if ':' not in reg:
		try:
			records = (x for x in data_handle if x[0].chromosome == int(reg) and not ',' in x[0].allele2)
		except:
			records = None
	else:
		start = reg.split(':')[1].split('-')[0]
		end = reg.split(':')[1].split('-')[1]
		try:
			records = (x for x in data_handle if x[0].chromosome == int(reg) and x[0].bp_position >= int(start) and x[0].bp_position <= int(end) and not ',' in x[0].allele2)
		except:
			records = None
	if snv_list is not None and records is not None:
		records = (x for x in records if snv_list[(snv_list['chr'] == x[0].chromosome) & (snv_list['pos'] == x[0].bp_position) & (snv_list['a1'] == x[0].allele1) & (snv_list['a2'] == x[0].allele2)].shape[0] > 0)
	return records

def ConvertPlink(chunk, sample_ids, iid_col):
	info=OrderedDict()
	info['chr'] = OrderedDict()
	info['pos'] = OrderedDict()
	info['variant'] = OrderedDict()
	info['a1'] = OrderedDict()
	info['a2'] = OrderedDict()
	data=OrderedDict()
	for locus, row in chunk:
		if locus.allele1 == '0':
			locus.allele1 = locus.allele2
		variant_unique = 'chr' + str(locus.chromosome) + 'bp' + str(locus.bp_position) + '.'  + str(locus.name) + '.' + str(locus.allele2) + '.' + str(locus.allele1)
		info['chr'][variant_unique] = str(locus.chromosome)
		info['variant'][variant_unique] = str(locus.name)
		info['pos'][variant_unique] = str(locus.bp_position)
		info['a1'][variant_unique] = str(locus.allele2)
		info['a2'][variant_unique] = str(locus.allele1)
		for sample, geno in zip(sample_ids, row):
			if not variant_unique in data.keys():
				data[variant_unique] = OrderedDict({sample.iid: geno})
			else:
				data[variant_unique][sample.iid] = geno if geno != 3 else 'NA'
	info = pd.DataFrame(info)
	data = pd.DataFrame(data)
	data = data.convert_objects(convert_numeric=True)
	data[iid_col] = data.index
	chunkdf = info.join(data.transpose())
	chunkdf.index = IndexChunk(chunkdf)
	return chunkdf

def ExtractDos1(data_handle, reg, snv_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[2]) >= int(start) and int(x[2]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if snv_list is not None and records is not None:
		records = (x for x in records if snv_list[(snv_list['chr'] == int(x[0])) & (snv_list['pos'] == int(x[2])) & (snv_list['a1'] == x[3]) & (snv_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDos1(chunk, sample_ids, iid_col):
	chunkdf = pd.DataFrame(chunk)
	chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
	chunkdf.columns = ['chr','pos','variant','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf

def ExtractDos2(data_handle, reg, snv_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[1]) >= int(start) and int(x[1]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if snv_list is not None and records is not None:
		records = (x for x in records if snv_list[(snv_list['chr'] == int(x[0])) & (snv_list['pos'] == int(x[1])) & (snv_list['a1'] == x[3]) & (snv_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDos2(chunk, sample_ids, iid_col):
	chunkdf = pd.DataFrame(chunk)
	chunkdf.columns = ['chr','pos','variant','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf

def ExtractOxford(data_handle, reg, snv_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[2]) >= int(start) and int(x[2]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if snv_list is not None and records is not None:
		records = (x for x in records if snv_list[(snv_list['chr'] == int(x[0])) & (snv_list['pos'] == int(x[2])) & (snv_list['a1'] == x[3]) & (snv_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDosage(row):
	newrow = row[:5]
	a=zip(row[5::3],row[6::3],row[7::3])
	newrow = newrow + [2*float(t[0]) + 1*float(t[1]) if float(t[0]) > 0 or float(t[1]) > 0 or float(t[2]) > 0 else float('nan') for t in a]
	return newrow

def ConvertOxford(chunk, sample_ids, iid_col):
	chunk = [ConvertDosage(row) for row in chunk]
	chunkdf = pd.DataFrame(chunk)
	chunkdf.columns = ['chr','pos','variant','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf
"""