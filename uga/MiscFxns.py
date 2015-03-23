import numpy as np
import pandas as pd
import re
from SystemFxns import Error
from multi_key_dict import multi_key_dict
from itertools import islice,takewhile,ifilter
from plinkio import plinkfile
import pysam
import vcf as VCF
import tabix
import gzip
	
def LoadPlink(data):
	plink_handle = plinkfile.open(data)
	plink_locus_it = plink_handle.get_loci()
	plink_sample_it = plink_handle.get_samples()
	sample_ids = [x.iid for x in plink_handle.samples]
	return plink_handle, plink_locus_it, plink_sample_it, sample_ids

def countPlinkRegion(iter, chr, start, end):
	k = 0
	for record in ifilter(lambda c: c.chromosome == chr and c.bp_position >= start and c.bp_position <= end, iter):
		k += 1
	return k

def LoadVcf(data):
	tb = tabix.open(data)
	v=VCF.Reader(filename=data)
	record=v.next()
	sample_ids = [a.sample for a in record.samples]
	return tb, sample_ids

def LoadDos(data, samples):
	tb = tabix.open(data)
	sample_ids = []
	open_func = gzip.open if samples[-3:] == '.gz' else open
	with open_func(samples) as sf:
		lines = (line.rstrip() for line in sf)
		lines = (line for line in lines if line)
		for line in lines:
			sample_ids.append(line)
	return tb, sample_ids

def countNonPlinkRegion(iter, type, start, end):
	k = 0
	if type in ['vcf','dos2']:
		for record in ifilter(lambda c: int(c[1]) >= start and int(c[1]) <= end, iter):
			k += 1
	else:
		for record in ifilter(lambda c: int(c[2]) >= start and int(c[2]) <= end, iter):
			k += 1
	return k

def Complement(x):
	if x != "NA":
		letters = list(x)
	else:
		letters = "NA"
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
	return ''.join(comp)

def ListCompatibleMarkers(chr,pos,a1,a2,delim):
	markers = []
	markers.append(chr + delim + pos + delim + a1 + delim + a2)
	markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(a2))
	markers.append(chr + delim + pos + delim + Complement(a2) + delim + Complement(a1))
	markers.append(chr + delim + pos + delim + a2 + delim + a1)
	for alt in a2.split(','):
		markers.append(chr + delim + pos + delim + a1 + delim + alt)
		markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(alt))
		markers.append(chr + delim + pos + delim + Complement(alt) + delim + Complement(a1))
		markers.append(chr + delim + pos + delim + alt + delim + a1)
	return list(set(markers))

def ConvertDosage(row):
	newrow = row[:5]
	a=zip(row[5::3],row[6::3],row[7::3])
	newrow = newrow + [2*float(t[0]) + 1*float(t[1]) if float(t[0]) > 0 or float(t[1]) > 0 or float(t[2]) > 0 else float('nan') for t in a]
	return newrow

def CalcCallrate(x):
	xlen = len(x)
	x = x.dropna().astype(float)
	if len(x) > 0:
		return float('%.5g' % (len(x)/float(xlen)))
	else:
		return 0.0

def CalcFreq(marker, chr, male_idx = None, female_idx = None):
	marker = marker.dropna().astype(float)
	if int(chr) == 23:
		if not male_idx is None and not female_idx is None:
			n = len(marker.loc[male_idx]) + 2 * len(marker.loc[female_idx])
			count = (marker.loc[male_idx].sum()/2) + marker.loc[female_idx].sum()
			return float('%.5g' % (count / n)) if len(marker.loc[male_idx]) > 0 and len(marker.loc[female_idx]) > 0 else float('NaN')
		else:
			float('NaN')
	else:
		n = 2 * len(marker)
		count = marker.sum()
		return float('%.5g' % (count / n)) if len(marker) > 0 else float('NaN')

def CalcRsq(x):
	x = x.dropna().astype(float)
	rsq = x.var() / (2 * (x.mean() / 2) * (1 - (x.mean() / 2)))
	rsq = 1 / rsq if rsq > 1 else rsq
	return float('%.5g' % (rsq))

def CalcHWE(marker, chr = None, female_idx = None):
	marker = marker.dropna().astype(float)
	if list(set(marker)) == [0,1,2]:
		if not chr is None and int(chr) == 23:
			if not female_idx is None:
				obs_hets = len(marker.loc[female_idx][marker.loc[female_idx].isin([1])])
				obs_hom1 = len(marker.loc[female_idx][marker.loc[female_idx].isin([2])])
				obs_hom2 = len(marker.loc[female_idx][marker.loc[female_idx].isin([0])])
			else:
				obs_hets = 0
				obs_hom1 = 0
				obs_hom2 = 0
		else:
			obs_hets = len(marker[marker.isin([1])])
			obs_hom1 = len(marker[marker.isin([2])])
			obs_hom2 = len(marker[marker.isin([0])])	
		if not (obs_hets == 0 and obs_hom1 == 0 and obs_hom2 == 0):
			obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
			obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

			rare_copies = 2 * obs_homr + obs_hets
			genotypes   = obs_hets + obs_homc + obs_homr

			het_probs = [0.0] * (rare_copies + 1)

			#start at midpoint
			mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)

			#check to ensure that midpoint and rare alleles have same parity
			if (rare_copies & 1) ^ (mid & 1):
				mid += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			het_probs[mid] = 1.0
			sum = float(het_probs[mid])

			for curr_hets in xrange(mid, 1, -2):
				het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
				sum += het_probs[curr_hets - 2]
				# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
				curr_homr += 1
				curr_homc += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			for curr_hets in xrange(mid, rare_copies - 1, 2):
				het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
				sum += het_probs[curr_hets + 2]
				#add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
				curr_homr -= 1
				curr_homc -= 1

			for i in xrange(0, rare_copies + 1):
				het_probs[i] /= sum

			#alternate p-value calculation for p_hi/p_lo
			p_hi = float(het_probs[obs_hets])
			for i in xrange(obs_hets, rare_copies+1):
				p_hi += het_probs[i]

			p_lo = float(het_probs[obs_hets])
			for i in xrange(obs_hets-1, -1, -1):
				p_lo += het_probs[i]

			p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

			p_hwe = 0.0
			#  p-value calculation for p_hwe
			for i in xrange(0, rare_copies + 1):
				if het_probs[i] > het_probs[obs_hets]:
					continue
				p_hwe += het_probs[i]

			p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

			return float('%.2e' % (p_hwe))
		else:
			return float('NaN')
	else:
		return float('NaN')

def GetRowCalls(row):
	newrow = row
	i = newrow[8].split(':').index('GT')
	newrow[9:] = newrow[9:].apply(lambda x: x.split(':')[i])
	newrow[9:] = newrow[9:].apply(lambda x: 'NaN' if not x in ['0/0','0/1','1/1','1/0'] else x.replace('0/0','2').replace('0/1','1').replace('1/1','0').replace('1/0','1'))
	return newrow

def ExtractModelVars(pheno,model,fid,iid,fxn=None,sex=None,case=None,ctrl=None,pheno_sep='\t'):
	model_vars_dict = {}
	dependent = re.split('\\~',model)[0]
	independent = re.split('\\~',model)[1]
	vars_df = pd.read_table(pheno,sep=pheno_sep,dtype='str')
	for x in [a for a in list(set(re.split('Surv\(|,|\)|~|\+|cluster\(|\(1\||\*|factor\(',model))) if a != '']:
		mtype = ''
		if dependent.find(x) != -1:
			pre = dependent.split(x)[0] if dependent.split(x)[0] != '' else '.'
			post = dependent.split(x)[1] if dependent.split(x)[1] != '' else '.'
			if not pre[-1].isalpha() and not post[0].isalpha():
				mtype = 'dependent'
		if independent.find(x) != -1:
			pre = independent.split(x)[0] if independent.split(x)[0] != '' else '.'
			post = independent.split(x)[1] if independent.split(x)[1] != '' else '.'
			if not pre[-1].isalpha() and not post[0].isalpha():
				if mtype == '':
					mtype = 'independent'
				else:
					print Error("a column in the phenotype file is defined in the model as both an independent and dependent variable")
					return
		if mtype != '':
			if model[model.find(x)-7:model.find(x)] == 'factor(':
				model_vars_dict[x] = {'class': 'factor', 'type': mtype}
			elif model[model.find(x)-15:model.find(x)] == 'ordered(factor(':
				model_vars_dict[x] = {'class': 'orderedfactor', 'type': mtype}
			elif model[model.find(x)-3:model.find(x)] == '(1|':
				model_vars_dict[x] = {'class': 'random', 'type': mtype}
			elif model[model.find(x)-8:model.find(x)] == 'cluster(':
				model_vars_dict[x] = {'class': 'cluster', 'type': mtype}
			else:
				model_vars_dict[x] = {'class': 'numeric', 'type': mtype}

	for x in model_vars_dict.keys():
		if x in list(vars_df.columns):
			print "   %s variable %s found" % (model_vars_dict[x]['type'], x)
		elif x in ['marker','marker1','marker2','marker.interact']:
			print "   %s variable %s skipped" % (model_vars_dict[x]['type'], x)
		else:
			print "   %s variable %s not found" % (model_vars_dict[x]['type'], x)
			print Error("model variable missing from phenotype file")
			return
	
	if sex:
		if sex in list(vars_df.columns):
			print "   sex column %s found" % sex
		else:
			print "   sex column %s not found" % sex
			print Error("sex column missing from phenotype file")
			return

	vars_df = vars_df[list(set([a for a in [fid,iid,sex] if a] + list([a for a in model_vars_dict.keys() if a != 'marker'])))]
	
	vars_df[fid] = vars_df[fid].astype(str)
	vars_df[iid] = vars_df[iid].astype(str)
	
	##### EXTRACT CASE/CTRL IF BINOMIAL fxn #####
	for x in model_vars_dict.keys():
		if model_vars_dict[x]['type'] == 'dependent' and fxn == 'binomial':
			vars_df = vars_df[vars_df[x].isin([str(case),str(ctrl)])]
			vars_df[x] = vars_df[x].map({str(ctrl): '0', str(case): '1'})
	vars_df.dropna(inplace = True)
	if len(vars_df.index) == 0:
		print Error("no data left for analysis")
		return
	return vars_df, model_vars_dict

def GetDelimiter(delimiter):
	if delimiter == 'tab':
		delimiter = '\t'
	elif delimiter == 'space':
		delimiter = ' '
	else:
		delimiter = ','
	return delimiter

def GetFocus(method,model,vars_df):
	focus = ['(Intercept)'] if method.split('_')[0] in ['gee','geeboss','glm','lme'] else []
	for x in re.sub("\+|\~|\-",",",model.split('~')[1]).split(','):
		if not x[0] == '(' and x.find('cluster(') == -1:
			if x.find('factor(') != -1:
				for v in vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64) - min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64)):
					if v != min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64) - min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64))):
						focus.append(x + str(v))
			elif x.find('*') != -1:
				focus.append('*'.join(sorted(x.split('*'))))
			else:
				focus.append(x)
	return focus

class MarkerRefDb(object):

	def __init__(self):
		self.ref = multi_key_dict()
		self.ref_alleles = multi_key_dict()

	def Update(self, row):
		newrow=row
		chr=newrow['chr']
		pos=newrow['pos']
		marker=newrow['marker']
		a1=newrow['a1']
		a2=newrow['a2']
		record_markers = ListCompatibleMarkers(chr,pos,a1,a2,'><')
		if len([r for r in record_markers if self.ref.has_key(r)]) > 0:
			temp_ref = self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			temp_alleles = self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			temp_alleles[1] = ','.join(list(set(','.join([b.split('><')[3] for b in record_markers if b.split('><')[2] == temp_alleles[0]]).split(','))))
			del self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			del self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			if not marker in temp_ref:
				temp_ref.append(marker)
				temp_ref = list(set(temp_ref))
			self.ref[tuple(record_markers)] = temp_ref
			self.ref_alleles[tuple(record_markers)] = temp_alleles
		else:
			self.ref[record_markers] = [marker]
			self.ref_alleles[record_markers] = [a1,a2]
		refmarker = filter(lambda x:'rs' in x, self.ref[chr + '><' + pos + '><' + a1 + '><' + a2])
		if len(refmarker) > 0:
			refmarker = refmarker[0]
		else:
			refmarker = self.ref[chr + '><' + pos + '><' + a1 + '><' + a2][0]
		keymarkers = self.ref[chr + '><' + pos + '><' + a1 + '><' + a2]
		refa1 = self.ref_alleles[chr + '><' + pos + '><' + a1 + '><' + a2][0]
		refa2 = self.ref_alleles[chr + '><' + pos + '><' + a1 + '><' + a2][1]
		if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
			newrow[5:] = 2-newrow[5:].str.replace('NA','nan').astype(float)
		newrow['marker']=refmarker	
		newrow['a1']=refa1
		newrow['a2']=refa2
		return newrow

class ChunkRefDb(object):

	def __init__(self):
		self.ref = multi_key_dict()
		self.ref_alleles = multi_key_dict()

	def ListKeys(self):
		return [b for t in self.ref.keys() for b in t]

	def Update(self,row):
		newrow=row
		chr=newrow['chr']
		pos=newrow['pos']
		marker=newrow['marker']
		a1=newrow['a1']
		a2=newrow['a2']
		record_markers = ListCompatibleMarkers(chr,pos,a1,a2,'><')
		if len([r for r in record_markers if self.ref.has_key(r)]) > 0:
			temp_ref = self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			temp_alleles = self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			temp_alleles[1] = ','.join(list(set(','.join([b.split('><')[3] for b in record_markers if b.split('><')[2] == temp_alleles[0]]).split(','))))
			del self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			del self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			if not marker in temp_ref:
				temp_ref.append(marker)
				temp_ref = list(set(temp_ref))
			self.ref[tuple(record_markers)] = temp_ref
			self.ref_alleles[tuple(record_markers)] = temp_alleles
		else:
			self.ref[record_markers] = [marker]
			self.ref_alleles[record_markers] = [a1,a2]
