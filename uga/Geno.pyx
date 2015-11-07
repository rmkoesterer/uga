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
	cdef public bytes filename, region, id
	cdef public np.ndarray samples
	cdef public object handle, region_iter
	cdef public unsigned int chr, start, end
	cdef public np.ndarray chunk, genos, marker_data, marker_info
	def __cinit__(self, filename):
		self.filename = filename

cdef class Vcf(Biodata):
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
			raise Exception

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
				if len(alts) > 1:
					c = np.append(c, np.empty((2*len(alts),10 + len(self.samples)), dtype='object'),axis=0)
					for alt in xrange(len(alts)):
						oth = [alts.index(a) for a in alts if a != alts[alt]]
						het = het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']] + [str(alts.index(k)+1) + sep + str(alt+1) for k in alts if k != alts[alt] for sep in ['/','|']] + [str(alt+1) + sep + str(alts.index(k)+1) for k in alts if k != alts[alt] for sep in ['/','|']]))
						hom1 = list(set(['0/0'] + [str(a+1) + '/' + str(a+1) for a in oth] + ['0|0'] + [str(a+1) + '|' + str(a+1) for a in oth]))
						hom2 = list(set([str(alt+1) + '/' + str(alt+1)] + [str(alt+1) + '|' + str(alt+1)]))
						c[i,:9] = record[:9]
						c[i,10:] = Vcf2Dose(record[9:], hom1, het, hom2, gt)
						c[i,3] = '&'.join([record[3]] + [alts[a] for a in oth])
						c[i,4] = alts[alt]
						c[i,5] = self.id
						c[i,9] = 'chr' + c[i,0] + 'bp' + c[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',c[i,2][0:60]) + '_'  + re_sub('&','_',c[i,3][0:20]) + '_'  + c[i,4][0:20]
						het = list(set(['0' + sep + str(alt+1) for sep in ['/','|']] + [str(alt+1) + sep + '0' for sep in ['/','|']]))
						hom1 = ['0/0','0|0']
						hom2 = [str(alt+1) + '/' + str(alt+1), str(alt+1) + '|' + str(alt+1)]
						c[i+1,:9] = record[:9]
						c[i+1,10:] = Vcf2Dose(record[9:], hom1, het, hom2, gt)
						c[i+1,3] = record[3]
						c[i+1,4] = alts[alt]
						c[i+1,5] = self.id
						c[i+1,9] = 'chr' + c[i+1,0] + 'bp' + c[i+1,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',c[i+1,2][0:60]) + '_'  + re_sub('&','_',c[i+1,3][0:20]) + '_'  + c[i+1,4][0:20]
						i += 2
				else:
					het = ['1/0','0/1','1|0','0|1']
					hom1 = ['0/0','0|0']
					hom2 = ['1/1','1|1']
					c[i,:9] = record[:9]
					c[i,10:] = Vcf2Dose(record[9:], hom1, het, hom2, gt)
					c[i,5] = self.id
					c[i,9] = 'chr' + c[i,0] + 'bp' + c[i,1] + '_'  + re_sub('!|@|#|\$|%|\^|&|\*|\(|\)|-|_|=|\+|\||{|}|\[|\]|;|:|\'|,|<|\.|>|/|\?|~','_',c[i,2][0:60]) + '_'  + re_sub('&','_',c[i,3][0:20]) + '_'  + c[i,4][0:20]
					i += 1
			c = c[(c[:i,1].astype(int) >= self.start) & (c[:i,1].astype(int) <= self.end)]
			self.marker_info = np.array([tuple(row) for row in c[:,[0,1,2,3,4,5,9]]], dtype=zip(np.array(['chr','pos','marker','a1','a2','id','marker_unique']),np.array(['uint8','uint32','|S60','|S20','|S20','|S100','|S100'])))
			self.marker_data = np.column_stack((self.samples, c[:,10:].transpose()))
			np.place(self.marker_data, self.marker_data == 'NA',np.nan)

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
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
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
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
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
	chunkdf.columns = ['chr','pos','marker','a1','a2'] + sample_ids
	chunkdf.index=IndexChunk(chunkdf)
	return chunkdf
"""