import os
import sys
import getopt
import subprocess
import gzip
import tabix
import math
import pandas as pd
from rpy2.robjects.packages import importr
from lib.Messages import *
from lib.MarkerCalc import *
from lib.Stats import *
from lib.Coordinates import *
import re
from itertools import islice
from Bio import bgzf

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def analyze(out = None, 
			data = None, 
			samples = None, 
			pheno = None, 
			model = None, 
			fid = None, 
			iid = None, 
			method = None, 
			focus = None, 
			sig = 5, 
			region_list = None, 
			region = None, 
			sex = None, 
			male = None, 
			female = None, 
			buffer = 100, 
			miss = None, 
			freq = None, 
			rsq = None, 
			hwe = None, 
			case = 1, 
			ctrl = 0, 
			nofail = False):

	assert out, usage(Error("an output file must be specified"))
	assert data, usage(Error("a data file must be specified"))
	assert samples, usage(Error("a sample file must be specified"))
	assert pheno, usage(Error("a phenotype file must be specified"))
	assert model, usage(Error("a model must be specified"))
	assert fid, usage(Error("a family ID column must be specified"))
	assert iid, usage(Error("an individual ID column must be specified"))
	assert method, usage(Error("an analysis method must be specified"))
	
	fxn = method.split('_')[1] if method in ["gee_gaussian","gee_binomial","glm_gaussian","glm_binomial","lme_gaussian","lme_binomial"] else ''
	
	if method in ["gee_gaussian","gee_binomial"]: geepack = importr('geepack')
	if method in ["lme_gaussian","lme_binomial"]: lme4 = importr('lme4')
	if method == "coxph": coxph = importr('survival')
	
	if focus: focus = focus.split(',')
	
	##### READ model VARIABLES FROM FILE #####
	print "   ... extracting model variables and removing missing/invalid samples"
	model_vars_dict = {}
	dependent = re.split('\\~',model)[0]
	independent = re.split('\\~',model)[1]
	vars_df = pd.read_table(pheno)
	for x in list(vars_df.columns) + ['marker','marker1','marker2','marker.interact']:
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
					usage(Error("a column in the phenotype file is defined in the model as both an independent and dependent variable"))			
		if mtype != '':
			if x[x.find('FID')-18:x.find('FID')] == 'as.numeric(factor(':
				model_vars_dict[x] = {'class': 'numericfactor', 'type': mtype}
			elif x[x.find('A1_SEX')-7:x.find('A1_SEX')] == 'factor(':
				model_vars_dict[x] = {'class': 'factor', 'type': mtype}
			elif x[x.find('FID')-15:x.find('FID')] == 'ordered(factor(':
				model_vars_dict[x] = {'class': 'orderedfactor', 'type': mtype}
			else:
				model_vars_dict[x] = {'class': 'numeric', 'type': mtype}
	print "   ... analysis model: %s" % model

	for x in model_vars_dict.keys():
		if x in list(vars_df.columns):
			print "          %s variable %s found" % (model_vars_dict[x]['type'], x)
		elif x in ['marker','marker1','marker2','marker.interact']:
			print "          %s variable %s skipped" % (model_vars_dict[x]['type'], x)
		else:
			print "          %s variable %s not found" % (model_vars_dict[x]['type'], x)
			usage(Error("model variable missing from phenotype file"))
	
	if sex:
		if sex in list(vars_df.columns):
			print "          sex column %s found" % sex
		else:
			print "          sex column %s not found" % sex
			usage(Error("sex column missing from phenotype file"))

	vars_df = vars_df[list(set([fid,iid,sex] + list([a for a in model_vars_dict.keys() if a != 'marker'])))]
	
	nmale = len(vars_df[sex][vars_df[sex] == male]) if sex and male and female else 'NA'
	nfemale = len(vars_df[sex][vars_df[sex] == female]) if sex and male and female else 'NA'
	
	##### EXTRACT CASE/CTRL IF BINOMIAL fxn #####
	for x in model_vars_dict.keys():
		if model_vars_dict[x]['type'] == 'dependent' and fxn == 'binomial':
			vars_df = vars_df[vars_df[x].isin([int(case),int(ctrl)])]
			vars_df[x] = vars_df[x].map({int(ctrl): 0, int(case): 1})
	vars_df.dropna(inplace = True)
	for i in model_vars_dict.keys():
		if model_vars_dict[i]['class'] == 'factor':
			vars_df[i] = pd.Categorical.from_array(vars_df[i]).labels
	if len(vars_df.index) == 0:
		usage(Error("no data left for analysis"))
		
	##### DETERMINE MODEL STATS TO BE EXTRACTED #####
	if focus == '':
		focus = ['(Intercept)'] if method.split('_')[0] in ['gee','glm','lme'] else []
		for x in re.sub("\+|\~|\-",",",model.split('~')[1]).split(','):
			if not x[0] == '(' and not x[-1] == ')':
				if x.find('factor(') != -1:
					focus.append(x + str(max(vars_df[re.sub('factor\(|\)','',x)])))
				elif x.find('*') != -1:
					focus.append(':'.join(sorted(x.split('*'))))
				else:
					focus.append(x)
	
	print "   ... reading sample list"
	sample_ids = []
	open_func = gzip.open if samples[-3:] == '.gz' else open
	with open_func(samples) as sf:
		lines = (line.rstrip() for line in sf)
		lines = (line for line in lines if line)
		for line in lines:
			sample_ids.append(line)
	if len(vars_df[vars_df[iid].isin(sample_ids)]) == 0:
		usage(Error("phenotype file and data file contain no common samples"))
	
	print "   ... generating unrelated sample list"
	unrelated_ids = vars_df.drop_duplicates(subset=[fid])[iid]
	
	samples = len(vars_df.index)
	samples_unique = len(vars_df.drop_duplicates(subset=[iid]).index)
	clusters = len(unrelated_ids)
	dep_var = [key for key in model_vars_dict if model_vars_dict[key]['type'] == 'dependent']
	cases = len(vars_df[dep_var[0]][vars_df[dep_var[0]].isin([1])]) if fxn == 'binomial' else 'NA'
	ctrls = len(vars_df[dep_var[0]][vars_df[dep_var[0]].isin([0])]) if fxn == 'binomial' else 'NA'
	print "   ... data summary"
	print "          " + str(samples) + " total observations"
	print "          " + str(samples_unique) + " unique samples"
	print "          " + str(clusters) + " clusters"
	print "          " + str(cases) + " cases"
	print "          " + str(ctrls) + " controls"
	print "          " + str(nmale) + " male"
	print "          " + str(nfemale) + " female"
	
	##### DETERMINE MARKERS TO BE ANALYZED #####
	if region_list:
		print "   ... reading list of regions from file"
		marker_list = Coordinates(region_list).Load()
	if region:
		marker_list = pd.DataFrame({'chr': [re.split(':|-',region)[0]],'start': [re.split(':|-',region)[1]],'end': [re.split(':|-',region)[2]],'region': [region]})
	marker_list['n'] = 0
	for i in range(len(marker_list.index)):
		tb = tabix.open(data)
		try:
			records = tb.querys(marker_list['region'][i])
		except:
			pass
		else:
			for record in records:
				marker_list['n'][i] = marker_list['n'][i] + 1
	if marker_list['n'].sum() == 0:
		print Error("no markers found")
		return()
	else:
		print "          " + str(marker_list['n'].sum()) + " markers in " + str(len(marker_list.index)) + " regions"

	sig = int(sig) if sig else sig
	buffer = int(buffer) if buffer else buffer

	##### RESET BUFFER TO MATCH NUMBER OF MARKERS IF CALC INDEP TESTS #####
	if method == 'indep_tests':
		buffer = marker_list['n'].max()
	vars_df.sort([fid],inplace = True)
	vars_df[fid] = pd.Categorical.from_array(vars_df[fid]).labels
	
	print "   ... starting marker analysis"
	written = False
	k = 0
	markers_noalt = 0
	bgzfile = bgzf.BgzfWriter(out + '.gz', 'wb')
	tb = tabix.open(data)
	i=0
	for r in range(len(marker_list.index)):
		reg = marker_list['region'][r]
		chr = reg.split(':')[0]
		try:
			records = tb.querys(reg)
		except:
			pass
		else:
			while True:
				i = i + 1
				chunk=list(islice(records, buffer))
				if not chunk:
					break
				k = 0
				chunkdf = pd.DataFrame(chunk)
				marker_info = chunkdf.ix[:,:4]
				marker_info.columns = ['chr','pos','marker','a1','a2']
				marker_info['marker_unique'] = 'chr' + marker_info['chr'].astype(str) + 'bp' + marker_info['pos'].astype(str) + '_'  + marker_info['marker'].astype(str) + '_'  + marker_info['a1'].astype(str) + '_'  + marker_info['a2'].astype(str)
				marker_info.index = marker_info['marker_unique']
				marker_data = chunkdf.ix[:,5:].transpose()
				marker_data.columns = marker_info['marker_unique']
				marker_data['IID'] = sample_ids
				model_df = pd.merge(vars_df, marker_data, on = [iid], how='left').sort([fid])
				marker_info['callrate']=model_df.drop_duplicates(subset=[iid])[marker_info['marker_unique']].apply(lambda col: CalcCallrate(col), 0)
				if sex and male and female:
					male_idx = model_df[model_df[sex].isin([male])].index.values
					female_idx = model_df[model_df[sex].isin([female])].index.values
				else:
					male_idx = None
					female_idx = None
				marker_info['freq']=model_df.drop_duplicates(subset=[iid])[marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				marker_info['freq.unrel']=model_df.drop_duplicates(subset=[iid])[model_df[iid].isin(unrelated_ids)][marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				if fxn == 'binomial':
					marker_info['freq.ctrl']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 0][marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.case']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 1][marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.unrel.ctrl']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 0)][marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.unrel.case']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 1)][marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				marker_info['rsq']=model_df.drop_duplicates(subset=[iid])[marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
				marker_info['rsq.unrel']=model_df.drop_duplicates(subset=[iid])[model_df[iid].isin(unrelated_ids)][marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
				if fxn == 'binomial':
					marker_info['rsq.ctrl']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 0][marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.case']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 1][marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.unrel.ctrl']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 0)][marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.unrel.case']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 1)][marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
				marker_info['hwe']=model_df.drop_duplicates(subset=[iid])[marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				marker_info['hwe.unrel']=model_df.drop_duplicates(subset=[iid])[model_df[iid].isin(unrelated_ids)][marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				if fxn == 'binomial':
					marker_info['hwe.ctrl']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 0][marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.case']=model_df.drop_duplicates(subset=[iid])[model_df[dep_var[0]] == 1][marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.unrel.ctrl']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 0)][marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.unrel.case']=model_df.drop_duplicates(subset=[iid])[(model_df[iid].isin(unrelated_ids)) & (model_df[dep_var[0]] == 1)][marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				marker_info['filter']=marker_info.apply(lambda row: GenerateFilterCode(marker_info=row, miss=miss, freq=freq, rsq=rsq, hwe=hwe), 1)
				marker_info['samples'] = str(samples) + '/' + str(samples_unique) + '/' + str(clusters) + '/' + str(cases) + '/' + str(ctrls)
				markercols = [col for col in model_df.columns if 'chr' in col]
				model_df[markercols] = model_df[markercols].astype(float)
				if method.split('_')[0] == 'gee':
					results = marker_info.apply(lambda row: CalcGEE(marker_info=row, model_df=model_df, model_vars_dict=model_vars_dict, model=model, iid=iid, fid=fid, method=method, fxn=fxn, focus=focus, dep_var=dep_var), 1)
				elif method.split('_')[0] == 'glm':
					results = marker_info.apply(lambda row: CalcGLM(marker_info=row, model_df=model_df, model_vars_dict=model_vars_dict, model=model, iid=iid, fid=fid, method=method, fxn=fxn, focus=focus, dep_var=dep_var), 1)
				elif method.split('_')[0] == 'lme':
					results = marker_info.apply(lambda row: CalcLME(marker_info=row, model_df=model_df, model_vars_dict=model_vars_dict, model=model, iid=iid, fid=fid, method=method, fxn=fxn, focus=focus, dep_var=dep_var), 1)
				elif method == 'coxph':
					results = marker_info.apply(lambda row: CalcCoxPH(marker_info=row, model_df=model_df, model_vars_dict=model_vars_dict, model=model, iid=iid, fid=fid, method=method, fxn=fxn, focus=focus, dep_var=dep_var), 1)
				elif method == 'indep_tests':
					n_indep, tot_tests = CalcIndepTests(marker_info=marker_info, model_df=model_df)
					results = pd.DataFrame({'chr': [marker_list['chr'][r]], 'start': [marker_list['start'][r]], 'end': [marker_list['end'][r]], 'reg_id': [marker_list['reg_id'][r]], 'n_total': [tot_tests], 'n_indep': [n_indep]})[['chr','start','end','reg_id','n_total','n_indep']]
				results.fillna('NA', inplace=True)
				if nofail:
					results = results[results['status'] > 0]
				if 'reg_id' in marker_list.columns:
					results['reg_id'] = marker_list['reg_id'][r]
				if not written:
					bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results.columns.values.tolist()]) + '\n')
					bgzfile.flush()
					written = True
				results.to_csv(bgzfile, header=False, sep='\t', index=False)
		pr = '   ... processed ' + str(min(i*buffer,marker_list['n'].sum())) + ' of ' + str(marker_list['n'].sum()) + ' markers' if method != 'indep_tests' else '   ... processed ' + str(marker_list['n'][r]) + ' markers from region ' + str(r+1) + ' of ' + str(len(marker_list.index))
		print pr
	bgzfile.close()
	print "   ... mapping results file"
	cmd = 'tabix -b 2 -e 2 ' + out + '.gz'
	p = subprocess.Popen(cmd, shell=True)
	p.wait()	
	print "   ... process complete"
