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
import pysam
import pandas as pd
import numpy as np
cimport numpy as np
cimport cython
from cython.view cimport array
from collections import OrderedDict
from Process import Error

#def IndexChunk(chunkdf):
#	return chunkdf['chr'].astype(str) + '><' + chunkdf['pos'].astype(str) + '><'  + chunkdf['a1'].astype(str).str[0:20].replace('-','_') + '><'  + chunkdf['a2'].astype(str).str[0:20].replace('-','_')

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] Vcf2Dose(np.ndarray genos, hom1, het, hom2, np.int gt):
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

cdef class Biodata(object):
	cdef public bytes filename, variant_list_file, region, id
	cdef public np.ndarray samples
	cdef public object handle, variant_list_handle, region_iter
	cdef public unsigned int chr, start, end
	cdef public np.ndarray variant_list, chunk, marker_data, marker_info
	def __cinit__(self, filename, variant_list_file = None):
		self.filename = filename
		self.variant_list_file = variant_list_file

		if self.variant_list_file is not None:
			print "loading variant list"
			try:
				self.variant_list_handle=pysam.TabixFile(filename=self.variant_list_file,parser=pysam.asTuple())
			except:
				raise Error("failed to load variant list iterator")

	#def get_variant_list_region(self, region):
	#	start = int(region.split(':')[1].split('-')[0])
	#	end = int(region.split(':')[1].split('-')[1])
	#	try:
	#		v = self.variant_list_handle.fetch(region=region, parser=pysam.asTuple())
	#	except:
	#		pass
	#	else:
	#		self.variant_list = np.array(['chr' + x[0] + 'bp' + x[1] + '_'  + x[2].replace('.','NA') + '_'  + x[3][0:20] + '_'  + x[4][0:20] for x in v if int(x[1]) >= start and int(x[1]) <= end])

cdef class Vcf(Biodata):
	def __cinit__(self, filename, variant_list_file = None):
		super(Vcf, self).__init__(filename, variant_list_file)

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
			raise Exception
		if self.variant_list_file is not None:
			try:
				v = self.variant_list_handle.fetch(region=region, parser=pysam.asTuple())
			except:
				pass
			else:
				self.variant_list = np.array(['chr' + x[0] + 'bp' + x[1] + '_'  + x[2][0:60].replace('.','NA') + '_'  + x[3][0:20] + '_'  + x[4][0:20] for x in v if int(x[1]) >= self.start and int(x[1]) <= self.end])

	cpdef get_chunk(self, int buffer):
		c = np.empty((buffer,10 + len(self.samples)), dtype='object')
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
				for alt in xrange(len(alts)):
					oth = [alts.index(a) for a in alts if a != alts[alt]]
					het1_unphased = [str(alt+1) + '/0'] + [str(alt+1) + '/' + str(a+1) for a in oth]
					het1_phased = [str(alt+1) + '|0'] + [str(alt+1) + '|' + str(a+1) for a in oth]
					het2_unphased = ['0/' + str(alt+1)] + [str(a+1) + '/' + str(alt+1) for a in oth]
					het2_phased = ['0|' + str(alt+1)] + [str(a+1) + '|' + str(alt+1) for a in oth]
					het = list(set(het1_unphased + het1_phased + het2_unphased + het2_phased))
					hom1_unphased = ['0/0'] + [str(a+1) + '/' + str(b+1) for a in oth for b in oth]
					hom1_phased = ['0|0'] + [str(a+1) + '|' + str(b+1) for a in oth for b in oth]
					hom1 = list(set(hom1_unphased + hom1_phased))
					hom2_unphased = [str(alt+1) + '/' + str(alt+1)]
					hom2_phased = [str(alt+1) + '|' + str(alt+1)]
					hom2 = list(set(hom2_unphased + hom2_phased))
					c[i,:9] = record[:9]
					c[i,10:] = Vcf2Dose(record[9:], hom1, het, hom2, gt)
					c[i,4] = ','.join([alts[alt]] + [alts[a] for a in oth])
					c[i,5] = self.id if self.id is not None else 'NA'
					c[i,9] = 'chr' + c[i,0] + 'bp' + c[i,1] + '_'  + c[i,2][0:60].replace('.','NA') + '_'  + c[i,3][0:20] + '_'  + c[i,4][0:20]
				i += 1
			self.chunk = c[(c[:i,1].astype(int) >= self.start) & (c[:i,1].astype(int) <= self.end)]
			if self.variant_list is not None:
				self.chunk = self.chunk[np.logical_or.reduce([self.chunk[:,9] == x for x in self.variant_list])]
			np.place(self.chunk,self.chunk == '.', 'NA')
			self.marker_info = np.array([tuple(row) for row in self.chunk[:,[0,1,2,3,4,5,9]]], dtype=zip(np.array(['chr','pos','marker','a1','a2','id','marker_unique']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100'])))
			self.marker_data = np.column_stack((self.samples, self.chunk[:,10:].transpose()))
			np.place(self.marker_data, self.marker_data == 'NA',np.nan)

"""
def ChunkVcf(region_iter, sample_ids, buffer, start, end, var_list = None):
	chunk = []
	for record in islice(region_iter, buffer):
		gt = record[8].split(':').index('GT')
		alts = record[4].split(',')
		for alt in xrange(len(alts)):
			oth = [alts.index(a) for a in alts if a != alts[alt]]
			het1_unphased = [str(alt+1) + '/0'] + [str(alt+1) + '/' + str(a+1) for a in oth]
			het1_phased = [str(alt+1) + '|0'] + [str(alt+1) + '|' + str(a+1) for a in oth]
			het2_unphased = ['0/' + str(alt+1)] + [str(a+1) + '/' + str(alt+1) for a in oth]
			het2_phased = ['0|' + str(alt+1)] + [str(a+1) + '|' + str(alt+1) for a in oth]
			het = list(set(het1_unphased + het1_phased + het2_unphased + het2_phased))
			hom1_unphased = ['0/0'] + [str(a+1) + '/' + str(b+1) for a in oth for b in oth]
			hom1_phased = ['0|0'] + [str(a+1) + '|' + str(b+1) for a in oth for b in oth]
			hom1 = list(set(hom1_unphased + hom1_phased))
			hom2_unphased = [str(alt+1) + '/' + str(alt+1)]
			hom2_phased = [str(alt+1) + '|' + str(alt+1)]
			hom2 = list(set(hom2_unphased + hom2_phased))
			newrecord = np.array(record)
			newrecord[9:] = Vcf2Dose(newrecord[9:], hom1, het, hom2, gt)
			newrecord[4] = ','.join([alts[alt]] + [alts[a] for a in oth])
			chunk.append(newrecord)
	out = pd.DataFrame(chunk, columns=['chr','pos','marker','a1','a2','qual','filter','info','format'] + sample_ids)
	if start is not None and end is not None:
		out = out[(out['pos'].astype(np.int) >= int(start)) & (out['pos'].astype(np.int) <= int(end))]
	if var_list:
		out = out[(out['chr'].astype(np.int) == var_list['chr'].astype(np.int)) & (out['pos'].astype(np.int) == var_list['pos'].astype(np.int)) & (out['a1'] == var_list['a1']) & (out['a2'] == var_list['a2'])]
	return out
"""

"""
def LoadVcf(file):
	print "loading vcf file iterator"
	try:
		v=pysam.TabixFile(filename=file,parser=pysam.asVCF())
	except:
		raise Error("failed to load vcf file iterator")
		return
	else:
		sample_ids = [a for a in v.header][-1].split('\t')[9:]
		return v, sample_ids

def ExtractVcf(data_handle, reg, variant_list = None):
	try:
		records = data_handle.fetch(region=reg, parser=pysam.asTuple())
	except:
		records = None
	return records

cdef double[:] Geno2Dose(np.ndarray genos, hom1, het, hom2, np.int gt):
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

def ExtractPlink(data_handle, reg, variant_list = None):
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
	if variant_list is not None and records is not None:
		records = (x for x in records if variant_list[(variant_list['chr'] == x[0].chromosome) & (variant_list['pos'] == x[0].bp_position) & (variant_list['a1'] == x[0].allele1) & (variant_list['a2'] == x[0].allele2)].shape[0] > 0)
	return records

def ConvertPlink(chunk, sample_ids, iid_col):
	marker_info=OrderedDict()
	marker_info['chr'] = OrderedDict()
	marker_info['pos'] = OrderedDict()
	marker_info['marker'] = OrderedDict()
	marker_info['a1'] = OrderedDict()
	marker_info['a2'] = OrderedDict()
	marker_data=OrderedDict()
	for locus, row in chunk:
		if locus.allele1 == '0':
			locus.allele1 = locus.allele2
		marker_unique = 'chr' + str(locus.chromosome) + 'bp' + str(locus.bp_position) + '.'  + str(locus.name) + '.' + str(locus.allele2) + '.' + str(locus.allele1)
		marker_info['chr'][marker_unique] = str(locus.chromosome)
		marker_info['marker'][marker_unique] = str(locus.name)
		marker_info['pos'][marker_unique] = str(locus.bp_position)
		marker_info['a1'][marker_unique] = str(locus.allele2)
		marker_info['a2'][marker_unique] = str(locus.allele1)
		for sample, geno in zip(sample_ids, row):
			if not marker_unique in marker_data.keys():
				marker_data[marker_unique] = OrderedDict({sample.iid: geno})
			else:
				marker_data[marker_unique][sample.iid] = geno if geno != 3 else 'NA'
	marker_info = pd.DataFrame(marker_info)
	marker_data = pd.DataFrame(marker_data)
	marker_data = marker_data.convert_objects(convert_numeric=True)
	marker_data[iid_col] = marker_data.index
	chunkdf = marker_info.join(marker_data.transpose())
	chunkdf.index = IndexChunk(chunkdf)
	return chunkdf

def ExtractDos1(data_handle, reg, variant_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[2]) >= int(start) and int(x[2]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if variant_list is not None and records is not None:
		records = (x for x in records if variant_list[(variant_list['chr'] == int(x[0])) & (variant_list['pos'] == int(x[2])) & (variant_list['a1'] == x[3]) & (variant_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDos1(chunk, sample_ids, iid_col):
	chunkdf = pd.DataFrame(chunk)
	chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf

def ExtractDos2(data_handle, reg, variant_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[1]) >= int(start) and int(x[1]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if variant_list is not None and records is not None:
		records = (x for x in records if variant_list[(variant_list['chr'] == int(x[0])) & (variant_list['pos'] == int(x[1])) & (variant_list['a1'] == x[3]) & (variant_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDos2(chunk, sample_ids, iid_col):
	chunkdf = pd.DataFrame(chunk)
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf

def ExtractOxford(data_handle, reg, variant_list = None):
	try:
		records = data_handle.querys(reg)
		if ':' in reg:
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			records = (x for x in records if int(x[2]) >= int(start) and int(x[2]) <= int(end) and not ',' in str(x[4]))
	except:
		records = None
	if variant_list is not None and records is not None:
		records = (x for x in records if variant_list[(variant_list['chr'] == int(x[0])) & (variant_list['pos'] == int(x[2])) & (variant_list['a1'] == x[3]) & (variant_list['a2'] == x[4])].shape[0] > 0)
	return records

def ConvertDosage(row):
	newrow = row[:5]
	a=zip(row[5::3],row[6::3],row[7::3])
	newrow = newrow + [2*float(t[0]) + 1*float(t[1]) if float(t[0]) > 0 or float(t[1]) > 0 or float(t[2]) > 0 else float('nan') for t in a]
	return newrow

def ConvertOxford(chunk, sample_ids, iid_col):
	chunk = [ConvertDosage(row) for row in chunk]
	chunkdf = pd.DataFrame(chunk)
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf
"""