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

from Model import *
	
def LoadPlink(data):
	plink_handle = plinkfile.open(data)
	plink_locus_it = plink_handle.get_loci()
	plink_sample_it = plink_handle.get_samples()
	sample_ids = [x.iid for x in plink_handle.samples]
	return plink_handle, plink_locus_it, plink_sample_it, sample_ids

def countPlinkRegion(iter, chr, start, end, ml = None):
	k = 0
	if ml is not None:
		filt = ifilter(lambda c: ml[(ml['chr'] == c.chromosome) & (ml['pos'] == c.bp_position) & (ml['a1'] == c.allele1) & (ml['a2'] == c.allele2)].shape[0] > 0, iter)
	else:
		filt = ifilter(lambda c: c.chromosome == chr and c.bp_position >= start and c.bp_position <= end, iter)
	for record in filt:
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

def LoadResults(data):
	tb = tabix.open(data)
	return tb

def countNonPlinkRegion(iter, type, start, end, ml = None):
	k = 0
	if type in ['vcf','dos2']:
		if ml is not None:
			filt = ifilter(lambda c: ml[(ml['chr'] == int(c[0])) & (ml['pos'] == int(c[1])) & (ml['a1'] == c[3]) & (ml['a2'] == c[4])].shape[0] > 0, iter)
		else:
			filt = ifilter(lambda c: int(c[1]) >= start and int(c[1]) <= end, iter)
		for record in filt:
			k += 1
	else:
		if ml is not None:
			filt = ifilter(lambda c: ml[(ml['chr'] == int(c[0])) & (ml['pos'] == int(c[2])) & (ml['a1'] == c[2]) & (ml['a2'] == c[4])].shape[0] > 0, iter)
		else:
			ifilter(lambda c: int(c[2]) >= start and int(c[2]) <= end, iter)
		for record in filt:
			k += 1
	return k

def Complement(x):
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

def ListCompatibleMarkers(chr,pos,a1,a2,delim):
	markers = []
	markers.append(chr + delim + pos + delim + a1 + delim + a2)
	markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(a2))
	markers.append(chr + delim + pos + delim + Complement(a2) + delim + Complement(a1))
	markers.append(chr + delim + pos + delim + a2 + delim + a1)
	markers.append(chr + delim + pos + delim + a1 + delim + 'NA')
	markers.append(chr + delim + pos + delim + Complement(a1) + delim + 'NA')
	markers.append(chr + delim + pos + delim + Complement(a2) + delim + 'NA')
	markers.append(chr + delim + pos + delim + a2 + delim + 'NA')
	markers.append(chr + delim + pos + delim + 'NA' + delim + a2)
	markers.append(chr + delim + pos + delim + 'NA' + delim + Complement(a2))
	markers.append(chr + delim + pos + delim + 'NA' + delim + Complement(a1))
	markers.append(chr + delim + pos + delim + 'NA' + delim + a1)
	for alt in a2.split(','):
		markers.append(chr + delim + pos + delim + a1 + delim + alt)
		markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(alt))
		markers.append(chr + delim + pos + delim + Complement(alt) + delim + Complement(a1))
		markers.append(chr + delim + pos + delim + alt + delim + a1)
		markers.append(chr + delim + pos + delim + a1 + delim + 'NA')
		markers.append(chr + delim + pos + delim + Complement(a1) + delim + 'NA')
		markers.append(chr + delim + pos + delim + Complement(alt) + delim + 'NA')
		markers.append(chr + delim + pos + delim + alt + delim + 'NA')
		markers.append(chr + delim + pos + delim + 'NA' + delim + alt)
		markers.append(chr + delim + pos + delim + 'NA' + delim + Complement(alt))
		markers.append(chr + delim + pos + delim + 'NA' + delim + Complement(a1))
		markers.append(chr + delim + pos + delim + 'NA' + delim + a1)
	return list(set([m for m in markers if not (m.split('><')[2] == 'NA' and m.split('><')[3] == 'NA')]))

def ListCompatibleMarkersMeta(chr,pos,a1,a2,delim):
	markers = []
	markers.append(chr + delim + pos + delim + a1 + delim + a2)
	markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(a2))
	markers.append(chr + delim + pos + delim + Complement(a2) + delim + Complement(a1))
	markers.append(chr + delim + pos + delim + a2 + delim + a1)
	return "_".join(sorted(markers))

def FlipEffect(refa1, refa2, a1, a2, effect):
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * effect
	else:
		return effect

def FlipFreq(refa1, refa2, a1, a2, freq):
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 - freq
	else:
		return freq

def FlipOR(refa1, refa2, a1, a2, o_r):
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 / o_r
	else:
		return o_r

def FlipZ(refa1, refa2, a1, a2, z):
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * z
	else:
		return z

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
	if list(set(marker)).issubset(set([0,1,2])):
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

def ExtractModelVars(pheno,model,fid,iid,family=None,sex=None,case=None,ctrl=None,sep='\t'):
	model_vars_dict = {}
	dependent = re.split('\\~',model)[0]
	independent = re.split('\\~',model)[1]
	vars_df = pd.read_table(pheno,sep=sep,dtype='str')
	vars_df[vars_df == "."] = None
	for x in [a for a in list(set(re.split('Surv\(|,|\)|~|\+|cluster\(|\(1\||\*|factor\(',dependent))) if a != '']:
		mtype = 'dependent'
		model_vars_dict[x] = {'class': 'numeric', 'type': mtype}
	for x in [a for a in list(set(re.split('Surv\(|,|\)|~|\+|cluster\(|\(1\||\*|factor\(',independent))) if a != '']:
		mtype = 'independent'
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
			print SystemFxns.Error("model variable missing from phenotype file")
			return
	
	if sex:
		if sex in list(vars_df.columns):
			print "   sex column %s found" % sex
		else:
			print "   sex column %s not found" % sex
			print SystemFxns.Error("sex column missing from phenotype file")
			return

	vars_df = vars_df[list(set([a for a in [fid,iid,sex] if a] + list([a for a in model_vars_dict.keys() if a != 'marker'])))]
	
	vars_df[fid] = vars_df[fid].astype(str)
	vars_df[iid] = vars_df[iid].astype(str)
	
	##### EXTRACT CASE/CTRL IF BINOMIAL family #####
	for x in model_vars_dict.keys():
		if model_vars_dict[x]['type'] == 'dependent' and family == 'binomial':
			vars_df = vars_df[vars_df[x].isin([str(case),str(ctrl)])]
			vars_df[x] = vars_df[x].map({str(ctrl): '0', str(case): '1'})
	vars_df.dropna(inplace = True)
	if len(vars_df.index) == 0:
		print SystemFxns.Error("no data left for analysis")
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

def GetFocus(model_fxn,model,vars_df,model_vars_dict):
	focus = ['(Intercept)'] if model_fxn.split('_')[0] in ['gee','geeboss','glm','lme'] else []
	for x in re.sub("\+|\~|\-",",",model.split('~')[1]).split(','):
		if not x[0] == '(' and x.find('cluster(') == -1:
			if x.find('factor(') != -1:
				x = re.sub('factor\(|\)','',x)
			if x.find('*') != -1:
				interact = []
				for a in x.split('*'):
					if model_vars_dict[a]['class'] == 'factor':
						for v in vars_df[a].unique().astype(np.int64) - min(vars_df[a].unique().astype(np.int64)):
							if v != min(vars_df[a].unique().astype(np.int64) - min(vars_df[a].unique().astype(np.int64))):
								interact.append('factor(' + a + ')' + str(v))
					else:
						interact.append(a)
				focus.append(interact[0] + '*' + interact[1])
			elif model_vars_dict[x]['class'] == 'factor':
				for v in vars_df[x].unique().astype(np.int64) - min(vars_df[x].unique().astype(np.int64)):
					if v != min(vars_df[x].unique().astype(np.int64) - min(vars_df[x].unique().astype(np.int64))):
						focus.append('factor(' + x + ')' + str(v))
			elif x.find('*') != -1:
				focus.append('*'.join(sorted(x.split('*'))))
			else:
				focus.append(x)
	return focus

def GenerateFilterCode(marker_info, no_mono=True, miss = None, freq = None, max_freq = None, rsq = None, hwe = None):
	filter = 0
	if (not miss is None and not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) < float(miss)) or (not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) == 0) or (math.isnan(marker_info['callrate'])):
		filter += 1000
	if not math.isnan(marker_info['freq']): 
		if no_mono and (float(marker_info['freq']) == 0 or float(marker_info['freq']) == 1):
			filter += 100
		else:
			if ((	not freq is None
				and 
					(		float(marker_info['freq']) < float(freq)
						 or float(marker_info['freq']) > 1-float(freq)
					)
				and 
					(		float(marker_info['freq.unrel']) < float(freq)
						 or float(marker_info['freq.unrel']) > 1-float(freq)
					)
			   ) 
			  or
			   (	not max_freq is None
				and 
					(	   
						(		float(marker_info['freq']) >= float(max_freq)
							and float(marker_info['freq']) <= 1-float(max_freq)
						)
					)
				and
					(
						(		float(marker_info['freq.unrel']) >= float(max_freq)
							and float(marker_info['freq.unrel']) <= 1-float(max_freq)
						)
						or float(marker_info['freq.unrel']) == 0
						or float(marker_info['freq.unrel']) == 1
					)
			   )):
				filter += 100
	if not rsq is None and not math.isnan(marker_info['rsq']) and (float(marker_info['rsq']) < float(rsq) and float(marker_info['rsq.unrel']) < float(rsq)):
		filter += 10
	if not hwe is None and not math.isnan(marker_info['hwe']) and (float(marker_info['hwe']) < float(hwe) and float(marker_info['hwe.unrel']) < float(hwe)):
		filter += 1
	return filter

def ChunkFlip(row, mdb):
	mdb_match = mdb[mdb.analogs.map(lambda x: set(x).issubset(set(row['analogs'])) or set(row['analogs']).issubset(set(x)))]
	if mdb_match.shape[0] > 0:
		a1=row['a1']
		a2=row['a2']
		refa1=mdb_match['a1'][0]
		refa2=mdb_match['a2'][0]
		row['a1'] = refa1
		row['a2'] = refa2 if refa2 != 'NA' else a2
		if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
			if refa1 + refa2 == Complement(a2) + 'NA':
				row['a2'] = Complement(a1)
			if refa1 + refa2 == a2 + 'NA':	
				row['a2'] = a1
			row[5:len(row)-1] = 2-row[5:len(row)-1].str.replace('NA','nan').astype(float)
	return row

def UpdateAltAllele(row, mdb):
	mdb_match = mdb[mdb.analogs.map(lambda x: '><'.join([row['chr'],row['pos'],row['a1'],row['a2']]) in x)]
	if mdb_match.shape[0] > 0:
		row['a2'] = mdb_match['a2'][0]
	return row['a2']

def MdbUpdate(row, chunkdf):
	chunkdf_match = chunkdf[chunkdf.analogs.map(lambda x: set(x).issubset(set(row['analogs'])) or set(row['analogs']).issubset(set(x)))]
	if chunkdf_match.shape[0] > 0:
		if row['a2'] == 'NA':
			row['a2'] = chunkdf_match['a2'][0]
		row['analogs'] = sorted(list(set(row['analogs'] + chunkdf_match['analogs'][0])))
	return row

def MdbAppend(mdb, chunkdf):
	mdb['check'] = mdb['analogs'].apply(lambda x: '_'.join(x),1)
	chunkdf['check'] = chunkdf['analogs'].apply(lambda x: '_'.join(x),1)
	chunkdf_nomatch = chunkdf[~chunkdf['check'].isin(mdb['check'])]
	if chunkdf_nomatch.shape[0] > 0:
		mdb = mdb.append(pd.DataFrame({'chr': list(chunkdf_nomatch['chr']),'pos': list(chunkdf_nomatch['pos']),'marker': list(chunkdf_nomatch['marker']),'a1': chunkdf_nomatch['a1'],'a2': chunkdf_nomatch['a2'], 
													'analogs': list(chunkdf_nomatch['analogs'])}))
	return mdb
