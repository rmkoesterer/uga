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

import os
import gzip
from progressbar import ProgressBar, Counter, Timer
from multiprocessing import Process, Manager, cpu_count
import signal
from Bio import bgzf
from re import split as re_split
import pandas as pd
import numpy as np
from glob import glob
from collections import OrderedDict
from plinkio import plinkfile
import tabix
import pysam
from Process import Error
import Fxns

def LoadPlink(file):
	print "loading plink file iterator"
	try:
		p = plinkfile.open(file)
	except:
		raise Error("failed to load plink file iterator")
		return
	else:
		z = zip(p.get_loci(),p)
		sample_ids = [x.iid for x in p.samples]
		return z, sample_ids

def LoadDos1(data, samples):
	print "loading dos1 file iterator"
	try:
		t = tabix.open(data)
	except:
		raise Error("failed to load dos1 file iterator")
		return
	else:
		print "loading dos1 sample file"
		sample_ids = []
		open_func = gzip.open if samples[-3:] == '.gz' else open
		try:
			with open_func(samples) as sf:
				lines = (line.rstrip() for line in sf)
				lines = (line for line in lines if line)
		except:
			raise Error("failed to load dos1 sample file " + samples)
			return
		else:
			for line in lines:
				sample_ids.append(line)
			return t, sample_ids

def LoadDos2(data, samples):
	print "loading dos2 file iterator"
	try:
		t = tabix.open(data)
	except:
		raise Error("failed to load dos2 file iterator")
		return
	else:
		print "loading dos2 sample file"
		sample_ids = []
		open_func = gzip.open if samples[-3:] == '.gz' else open
		try:
			with open_func(samples) as sf:
				lines = (line.rstrip() for line in sf)
				lines = (line for line in lines if line)
		except:
			raise Error("failed to load dos2 sample file " + samples)
			return
		else:
			for line in lines:
				sample_ids.append(line)
			return t, sample_ids

def LoadOxford(data, samples):
	print "loading oxford file iterator"
	try:
		t = tabix.open(data)
	except:
		raise Error("failed to load oxford file iterator")
		return
	else:
		print "loading oxford sample file"
		sample_ids = []
		open_func = gzip.open if samples[-3:] == '.gz' else open
		try:
			with open_func(samples) as sf:
				lines = (line.rstrip() for line in sf)
				lines = (line for line in lines if line)
		except:
			raise Error("failed to load oxford sample file " + samples)
			return
		else:
			for line in lines:
				sample_ids.append(line)
			return t, sample_ids

def LoadVariantList(file):
	print "loading variant list"
	try:
		t = tabix.open(file)
	except:
		raise Error("failed to load variant list iterator")
		return
	else:
		return t

def LoadResults(f):
	print "loading results file iterator"
	try:
		t = tabix.open(f)
	except:
		raise Error("failed to load bgzip format results file " + f)
		return
	else:
		return t

def LoadPheno(pheno,model_fields,fid,iid,matid=None,patid=None,family=None,sex=None,case=None,ctrl=None,sep='tab'):
	print "extracting model fields from pheno file and reducing to complete observations ..."
	sep = Fxns.GetDelimiter(sep)
	vars_df = pd.read_table(pheno,sep=sep,dtype='str')
	vars_df[vars_df == "."] = None
	for x in model_fields.keys():
		if x in list(vars_df.columns):
			print "   %s variable %s found" % (model_fields[x]['type'], x)
		elif x in ['marker','marker1','marker2','marker.interact']:
			print "   %s variable %s skipped" % (model_fields[x]['type'], x)
		else:
			raise Error("column " + x + " not found in phenotype file " + pheno)
			return
	if sex:
		if sex in list(vars_df.columns):
			print "   sex column %s found" % sex
		else:
			raise Error("column " + sex + " not found in phenotype file " + pheno)
			return
	vars_df = vars_df[list(set([a for a in [fid,iid,matid,patid,sex] if a] + list([a for a in model_fields.keys() if a != 'marker'])))]
	vars_df[fid] = vars_df[fid].astype(str)
	vars_df[iid] = vars_df[iid].astype(str)
	if matid:
		vars_df[matid] = vars_df[matid].astype(str)
	if patid:
		vars_df[patid] = vars_df[patid].astype(str)
	for x in model_fields.keys():
		if model_fields[x]['type'] == 'dependent' and family == 'binomial':
			vars_df = vars_df[vars_df[x].isin([str(case),str(ctrl)])]
			vars_df[x] = vars_df[x].map({str(ctrl): '0', str(case): '1'})
	vars_df.dropna(inplace = True)
	if len(vars_df.index) == 0:
		raise Error("no non-missing data left for analysis in phenotype file " + pheno)
		return
	return vars_df

def PrepareChrDirs(regions, directory):
	try:
		os.mkdir(directory.replace('chr[CHR]/',''))
	except OSError:
		pass
	for chr in set([a.split(":")[0] for a in regions]):
		try:
			os.mkdir(directory.replace("[CHR]",chr))
		except OSError:
			continue

def PrepareListDirs(n, directory):
	try:
		os.mkdir(directory.replace('list[LIST]/',''))
	except OSError:
		pass
	for a in range(int(np.ceil(n/100.0))):
		try:
			os.mkdir(directory.replace("[LIST]",str(int(np.floor(a/100.0) + 100*a)) + '-' + str(int(np.floor(a/100.0) + 100*a+ + 99))))
		except OSError:
			continue
			
def GenerateSubFiles(region_df, f, dist_mode, n):
	out_files=OrderedDict()
	for i in range(n):
		if dist_mode in ['region','split-list']:
			of = f.replace('[CHR]',str(region_df['chr'][i])) + '.chr' + str(region_df['chr'][i]) + 'bp' + str(region_df['start'][i]) + '-' + str(region_df['end'][i])
			out_files[str(region_df['chr'][i]) + ':' + str(region_df['start'][i]) + '-' + str(region_df['end'][i])] = of
		elif dist_mode == 'split-list-n':
			of = f.replace('[LIST]',str(int(np.floor(i/100.0) + 99*np.floor(i/100.0))) + '-' + str(int(np.floor(i/100.0) + 99*np.floor(i/100.0) + 99))) + '.list' + str(i)
			out_files[i] = of
		elif dist_mode == 'chr':
			of = f + '.chr' + str(region_df['chr'][i])
			out_files[str(region_df['chr'][i])] = of
		else:
			print Error("invalid dist_mode: " + dist_mode)
			return
	return out_files

def FindSubFiles(f):
	out_files=OrderedDict()
	list_files = glob("list*/" + f + ".list*.gz")
	chr_files = glob("chr*/" + f + ".chr*bp*.gz")
	if len(list_files) > 0 and len(chr_files) > 0:
		print Error("found both list and chr directories containing valid data files")
		return
	elif len(list_files) > 0:
		for t in sorted([(int(k.split('list')[2].split('.')[0]),k.replace('.gz','')) for k in list_files],key=lambda l: l[0]):
			out_files[t[0]] = t[1]
	else:
		for t in sorted([(int(k.split('chr')[2].split('bp')[0]),int(k.split('chr')[2].split('bp')[1].split('-')[0]),int(k.split('chr')[2].split('bp')[1].split('-')[1].split('.')[0]),k.replace('.gz','')) for k in chr_files],key=lambda l: (l[0],l[1],l[2])):
			out_files[str(t[0]) + ':' + str(t[1]) + '-' + str(t[2])] = t[3]
	return out_files

class Regions(object):
	def __init__(self, filename = None, region = None, id=None):
		self.filename = filename
		self.region = region
		self.id = id
		if self.filename is not None and self.region is not None:
			print "region list option provided, region option will be ignored"
		elif self.filename is None and self.region is None:
			self.df = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': [1 for i in range(26)],'end': [1000000000 for i in range(26)],'region': [str(i+1) + ':1-1000000000' for i in range(26)],'id': [None for i in range(26)]})
		elif self.filename is not None:
			print "loading region list"
			try:
				with open(self.filename) as f:
					self.df = pd.read_table(f, header=None)
			except:
				raise Error("unable to load region list " + self.filename)
			else:
				if self.df.shape[1] == 2:
					self.df.columns=['region','id']
				elif self.df.shape[1] == 1:
					self.df.columns=['region']
				else:
					raise Error("region list contains too many columns")
				if not 'id' in self.df.columns:
					self.df['id'] = None
				self.df['chr'] = self.df['region'].apply(lambda row: int(row.split(':')[0]),1)
				self.df['start'] = self.df['region'].apply(lambda row: int(row.split(':')[1].split('-')[0]),1)
				self.df['end'] = self.df['region'].apply(lambda row: int(row.split(':')[1].split('-')[1]),1)
				self.df['chr'] = self.df['chr'].astype(int)
				self.df['start'] = self.df['start'].astype(int)
				self.df['end'] = self.df['end'].astype(int)
				self.df = self.df.sort_index(by=['chr','start'],ascending=[True,True])
		else:
			if len(self.region.split(':')) > 1:
				self.df = pd.DataFrame({'chr': [re_split(':|-',self.region)[0]],'start': [re_split(':|-',self.region)[1]],'end': [re_split(':|-',self.region)[2]],'region': [self.region],'id': [self.id]})
			else:
				self.df = pd.DataFrame({'chr': [self.region],'start': [1],'end': [1000000000],'region': [self.region + ':1-1000000000'],'id': [self.id]})
			self.id = self.id if not self.id is None else 'NA'
