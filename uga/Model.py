import os
import sys
import getopt
import subprocess
import gzip
import tabix
import math
import numpy as np
import pandas as pd
from rpy2.robjects.packages import importr
import pandas.rpy.common as py2r
import re
from itertools import islice,takewhile
from Bio import bgzf
import psutil
from MarkerCalc import Complement,ConvertDosage,CalcCallrate,CalcFreq,CalcRsq,CalcHWE,CallToDos
from Messages import Error
from Stats import GenerateFilterCode,CalcGEE,CalcGLM,CalcLME,CalcCoxPH,CalcEffTests,CalcFamSkatO,CalcFamSkat
from Coordinates import Coordinates
from ModelFxns import ExtractModelVars,GetFocus,GetDelimiter
from plinkio import plinkfile
import collections
import pysam
import vcf as VCF
from Iterators import LoadPlink,LoadVcf,LoadDos
		
#from memory_profiler import profile, memory_usage

pd.options.mode.chained_assignment = None

#@profile
def Model(out = None, 
			cfg = None, 
			data = [None], 
			format = ['oxford'], 
			samples = [None], 
			pheno = [None], 
			model = [None], 
			fid = [None], 
			iid = [None], 
			method = [None], 
			focus = [None], 
			sig = 5, 
			region_list = None, 
			region = None, 
			region_id = None, 
			pedigree = [None], 
			sex = [None], 
			male = [None], 
			female = [None], 
			buffer = 100, 
			miss = None, 
			freq = None, 
			rsq = None, 
			hwe = None, 
			case = [1], 
			ctrl = [0], 
			mem = 3, 
			corstr = ['exchangeable'], 
			delimiter = ['tab'], 
			nofail = False):

	print "model options ..."
	if not cfg is None:
		for arg in ['cfg','region','region_list','region_id','mem']:
			if not str(locals()[arg]) in ['None','False']:
				print "   {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])
	else:
		for arg in locals().keys():
			if not str(locals()[arg]) in ['None','False']:
				print "   {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	##### populate configuration #####
	if cfg is None:
		cfg={'out': out, 'buffer': int(buffer), 'hwe': hwe, 'data_order': ['NA'], 'freq': freq, 'miss': freq, 'rsq': freq, 'mem': mem, 'sig': int(sig), 'nofail': nofail,
				'data_info': {'NA': {'data': data[0], 'format': format[0], 'samples': samples[0], 'pheno': pheno[0], 'model': model[0], 'fid': fid[0], 'iid': iid[0],
					'method': method[0], 'focus': focus[0], 'pedigree': pedigree[0], 'sex': sex[0], 'male': male[0], 'female': female[0], 'case': case[0], 'ctrl': ctrl[0], 
						'corstr': corstr[0], 'delimiter': delimiter[0]}}}

	##### READ model VARIABLES FROM FILE #####
	vars_df_dict = {}
	model_vars_dict_dict = {}
	for k in cfg['data_order']:
		if len(cfg['data_info'].keys()) > 1:
			print "extracting model variables for model " + k + " and removing missing/invalid samples ..."
		else:
			print "extracting model variables and removing missing/invalid samples ..."
		cfg['data_info'][k]['delimiter'] = GetDelimiter(cfg['data_info'][k]['delimiter'])
		cfg['data_info'][k]['fxn'] = cfg['data_info'][k]['method'].split('_')[1] if cfg['data_info'][k]['method'] in ["gee_gaussian","gee_binomial","glm_gaussian","glm_binomial","lme_gaussian","lme_binomial"] else None
		cfg['data_info'][k]['vars_df'], cfg['data_info'][k]['model_vars_dict'] = ExtractModelVars(cfg['data_info'][k]['pheno'],cfg['data_info'][k]['model'],cfg['data_info'][k]['fid'],cfg['data_info'][k]['iid'],fxn=cfg['data_info'][k]['fxn'],sex=cfg['data_info'][k]['sex'],delimiter=cfg['data_info'][k]['delimiter'])

	##### DETERMINE MODEL STATS TO BE EXTRACTED #####
	for k in cfg['data_order']:
		if 'focus' in cfg['data_info'][k].keys() and not cfg['data_info'][k]['focus'] is None:
			cfg['data_info'][k]['focus'] = cfg['data_info'][k]['focus'].split(',')
		else:
			cfg['data_info'][k]['focus'] = GetFocus(method=cfg['data_info'][k]['method'],model=cfg['data_info'][k]['model'],vars_df=cfg['data_info'][k]['vars_df'])

	##### LOAD ITERATORS AND SAMPLE LISTS #####
	for k in cfg['data_order']:
		if cfg['data_info'][k]['format'] == 'plink':
			print "reading Plink binary files for model " + k if len(cfg['data_info'].keys()) > 1 else "reading Plink binary files"
			cfg['data_info'][k]['data_it'], cfg['data_info'][k]['sample_it'], cfg['data_info'][k]['sample_ids'] = LoadPlink(cfg['data_info'][k]['data'])
		elif cfg['data_info'][k]['format'] == 'vcf':
			print "reading vcf file for model " + k if len(cfg['data_info'].keys()) > 1 else "reading vcf file"
			cfg['data_info'][k]['data_it'], cfg['data_info'][k]['sample_ids'] = LoadVcf(cfg['data_info'][k]['data'])
		else:
			print "reading data and sample files for model " + k if len(cfg['data_info'].keys()) > 1 else "reading data and sample files"
			cfg['data_info'][k]['data_it'], cfg['data_info'][k]['sample_ids'] = LoadDos(cfg['data_info'][k]['data'], cfg['data_info'][k]['samples'])
		if len(cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].isin(cfg['data_info'][k]['sample_ids'])]) == 0:
			print Error("phenotype file and data file contain no common samples")
			return

	##### REDUCE PHENO DATA TO GENOTYPE IDs #####
	for k in cfg['data_order']:
		cfg['data_info'][k]['vars_df'] = cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].isin(cfg['data_info'][k]['sample_ids'])]

	##### CREATE SUMMARY FOR DATA TO BE ANALYZED #####
	for k in cfg['data_order']:
		cfg['data_info'][k]['samples'] = len(cfg['data_info'][k]['vars_df'].index)
		cfg['data_info'][k]['vars_df_nodup'] = cfg['data_info'][k]['vars_df'].drop_duplicates(subset=[cfg['data_info'][k]['iid']])
		cfg['data_info'][k]['samples_unique'] = len(cfg['data_info'][k]['vars_df_nodup'].index)
		cfg['data_info'][k]['clusters'] = len(cfg['data_info'][k]['vars_df_nodup'].drop_duplicates(subset=[cfg['data_info'][k]['fid']]).index)
		cfg['data_info'][k]['dep_var'] = [key for key in cfg['data_info'][k]['model_vars_dict'] if cfg['data_info'][k]['model_vars_dict'][key]['type'] == 'dependent']
		cfg['data_info'][k]['cases'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['dep_var'][0]][cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['dep_var'][0]].isin(['1'])]) if cfg['data_info'][k]['fxn'] == 'binomial' else 'NA'
		cfg['data_info'][k]['ctrls'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['dep_var'][0]][cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['dep_var'][0]].isin(['0'])]) if cfg['data_info'][k]['fxn'] == 'binomial' else 'NA'
		cfg['data_info'][k]['nmale'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['sex']].isin([str(cfg['data_info'][k]['male'])])].index.values) if not cfg['data_info'][k]['sex'] is None and not cfg['data_info'][k]['male'] is None and not cfg['data_info'][k]['female'] is None else 'NA'
		cfg['data_info'][k]['nfemale'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['sex']].isin([str(cfg['data_info'][k]['female'])])].index.values) if not cfg['data_info'][k]['sex'] is None and not cfg['data_info'][k]['male'] is None and not cfg['data_info'][k]['female'] is None else 'NA'
		if len(cfg['data_info'].keys()) > 1:
			print "data summary for model " + k + " ..."
		else:
			print "data summary ..."
		print "   " + str(cfg['data_info'][k]['samples']) + " total observations"
		print "   " + str(cfg['data_info'][k]['samples_unique']) + " unique samples"
		print "   " + str(cfg['data_info'][k]['clusters']) + " clusters"
		print "   " + str(cfg['data_info'][k]['cases']) + " case"
		print "   " + str(cfg['data_info'][k]['ctrls']) + " control"
		print "   " + str(cfg['data_info'][k]['nmale']) + " male"
		print "   " + str(cfg['data_info'][k]['nfemale']) + " female"

	##### READ PEDIGREE FROM FILE FOR FAMILY BASED SKAT TEST #####
	for k in cfg['data_order']:
		if 'pedigree' in cfg['data_info'][k].keys() and not cfg['data_info'][k]['pedigree'] is None:
			print "extracting pedigree from file for model " + k if len(cfg['data_info'].keys()) > 1 else "extracting pedigree from file"
			cfg['data_info'][k]['ped_df'] = pd.read_table(cfg['data_info'][k]['pedigree'],sep='\t',dtype='str',header=None, names=['FID','IID','PAT','MAT'])
			cfg['data_info'][k]['ped_df'] = cfg['data_info'][k]['ped_df'][cfg['data_info'][k]['ped_df']['IID'].isin(list(cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].values))]

	##### GENERATE REGION LIST #####
	if region_list:
		print "reading list of regions from file"
		marker_list = Coordinates(region_list).Load()
	elif region:
		if len(region.split(':')) > 1:
			marker_list = pd.DataFrame({'chr': [re.split(':|-',region)[0]],'start': [re.split(':|-',region)[1]],'end': [re.split(':|-',region)[2]],'region': [region]})
		else:
			marker_list = pd.DataFrame({'chr': [region],'start': ['NA'],'end': ['NA'],'region': [region]})
		marker_list['reg_id'] = region_id
	else:
		marker_list = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)]})
	
	##### DETERMINE REGION COUNTS #####
	for k in cfg['data_order']:
		marker_list[k] = 0
		if cfg['data_info'][k]['format'] == 'plink':
			for i in range(len(marker_list.index)):
				if marker_list['start'][i] != 'NA':
					marker_list[k][i] = len([c for c in cfg['data_info'][k]['data_it'] if c.chromosome == int(marker_list['chr'][i]) and c.bp_position >= int(marker_list['start'][i]) and c.bp_position <= int(marker_list['end'][i])])
				else:
					marker_list[k][i] = len(cfg['data_info'][k]['data_it'])
		else:
			tb = tabix.open(cfg['data_info'][k]['data'])
			for i in range(len(marker_list.index)):
				try:
					records = tb.querys(marker_list['region'][i])
					records = (x for x in records if int(x[1]) >= int(marker_list['start'][i]) and int(x[1]) <= int(marker_list['end'][i]))
				except:
					pass
				else:
					for record in records:
						if marker_list['start'][i] != 'NA' and int(marker_list['end'][i]) - int(marker_list['start'][i]) <= 10000000:
							marker_list[k][i] = marker_list[k][i] + 1
						else:
							marker_list[k][i] = 'n'
							break
	#marker_list = marker_list[(marker_list[k] > 0].reset_index(drop=True)
	#if marker_list['n'].sum() == 0:
	#	print Error("no markers found")
	#	return()
	#else:
	#	print "   " + str(len(marker_list.index)) + " non-empty regions"

	##### START ANALYSIS #####
	print "modelling data ..."
	for k in cfg['data_order']:
		written = False
		cfg['data_info'][k]['out'] = cfg['out'] + '.' + k + '.gz' if len(cfg['data_info'].keys()) > 1 else cfg['out'] + '.gz'
		bgzfile = bgzf.BgzfWriter(cfg['data_info'][k]['out'], 'wb')
		for r in range(len(marker_list.index)):
			reg = marker_list['region'][r]
			chr = reg.split(':')[0]
			i = 0
			if cfg['data_info'][k]['format'] == 'plink':
				if reg == chr:
					try:
						records = (x for x in cfg['data_info'][k]['data_it'] if x.chromosome == int(chr))
					except:
						break
				else:
					try:
						records = (x for x in cfg['data_info'][k]['data_it'] if x.chromosome == int(chr) and x.bp_position >= int(marker_list['start'][r]) and x.bp_position <= int(marker_list['end'][r]))
					except:
						break
			else:
				try:
					records = cfg['data_info'][k]['data_it'].querys(reg)
					records = (x for x in records if int(x[1]) >= int(marker_list['start'][r]) and int(x[1]) <= int(marker_list['end'][r]))
				except:
					break
			while True:
				i = i + 1
				chunk=list(islice(records, cfg['buffer']))
				if not chunk:
					break
				if cfg['data_info'][k]['format'] == 'plink':
					marker_info=collections.OrderedDict()
					marker_info['chr'] = collections.OrderedDict()
					marker_info['pos'] = collections.OrderedDict()
					marker_info['marker'] = collections.OrderedDict()
					marker_info['a1'] = collections.OrderedDict()
					marker_info['a2'] = collections.OrderedDict()
					marker_data=collections.OrderedDict()
					for locus, row in zip(chunk, cfg['data_info'][k]['data_it']):
						if locus.allele1 == '0':
							locus.allele1 = locus.allele2
						marker_unique = 'chr' + str(locus.chromosome) + 'bp' + str(locus.bp_position) + '.'  + locus.name + '.' + locus.allele2 + '.' + locus.allele1
						marker_info['chr'][marker_unique] = locus.chromosome
						marker_info['marker'][marker_unique] = locus.name
						marker_info['pos'][marker_unique] = locus.bp_position
						marker_info['a1'][marker_unique] = locus.allele2
						marker_info['a2'][marker_unique] = locus.allele1
						for sample, geno in zip(cfg['data_info'][k]['sample_it'], row):
							if not marker_unique in marker_data.keys():
								marker_data[marker_unique] = collections.OrderedDict({sample.iid: geno})
							else:
								marker_data[marker_unique][sample.iid] = geno if geno != 3 else 'NA'
					marker_info = pd.DataFrame(marker_info)
					marker_info['marker_unique'] = marker_info.index
					marker_data = pd.DataFrame(marker_data)
					marker_data = marker_data.convert_objects(convert_numeric=True)
					marker_data[cfg['data_info'][k]['iid']] = marker_data.index
				else:
					chunkdf = pd.DataFrame(chunk)
					marker_info = chunkdf.ix[:,:4]
					marker_info.replace('.','NA',inplace=True)
					marker_info.columns = ['chr','pos','marker','a1','a2'] if cfg['data_info'][k]['format'] in ['dos2','vcf'] else ['chr','marker','pos','a1','a2']
					marker_info = marker_info[['chr','pos','marker','a1','a2']]
					marker_info['marker_unique'] = 'chr' + marker_info['chr'].astype(str) + 'bp' + marker_info['pos'].astype(str) + '.'  + marker_info['marker'].astype(str) + '.'  + marker_info['a1'].astype(str) + '.'  + marker_info['a2'].astype(str)
					marker_info.index = marker_info['marker_unique']
					if cfg['data_info'][k]['format'] == 'vcf':
						marker_data = chunkdf.ix[:,9:].transpose()
						marker_data=marker_data.applymap(CallToDos)
					elif cfg['data_info'][k]['format'] == 'oxford':
						marker_data = chunkdf.ix[:,5:].transpose()
						marker_data = marker_data.apply(lambda col: pd.Series(np.array(ConvertDosage(list(col))).astype(np.float64)),0)
					else:
						marker_data = chunkdf.ix[:,5:].transpose()
					marker_data = marker_data.convert_objects(convert_numeric=True)
					marker_data.columns = marker_info['marker_unique']
					marker_data[cfg['data_info'][k]['iid']] = cfg['data_info'][k]['sample_ids']
				model_df = pd.merge(cfg['data_info'][k]['vars_df'], marker_data, on = [cfg['data_info'][k]['iid']], how='left').sort([cfg['data_info'][k]['fid']])
				model_df_nodup=model_df.drop_duplicates(subset=[cfg['data_info'][k]['iid']]).reset_index(drop=True)
				marker_info['callrate']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: CalcCallrate(col), 0)
				model_df_nodup_unrel=model_df_nodup.drop_duplicates(subset=[cfg['data_info'][k]['fid']]).reset_index(drop=True)
				if not cfg['data_info'][k]['sex'] is None and not cfg['data_info'][k]['male'] is None and not cfg['data_info'][k]['female']:
					male_idx = model_df_nodup[model_df_nodup[cfg['data_info'][k]['sex']].isin([cfg['data_info'][k]['male']])].index.values
					female_idx = model_df_nodup[model_df_nodup[cfg['data_info'][k]['sex']].isin([cfg['data_info'][k]['female']])].index.values
				else:
					male_idx = None
					female_idx = None
				marker_info['freq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				marker_info['freq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['freq.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
					marker_info['freq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcFreq(marker=col, chr = chr, male_idx = male_idx, female_idx = female_idx))
				marker_info['rsq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
				marker_info['rsq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: CalcRsq(col), 0)
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['rsq.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcRsq(col), 0)
					marker_info['rsq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcRsq(col), 0)
				marker_info['hwe']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				marker_info['hwe.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['hwe.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
					marker_info['hwe.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: CalcHWE(marker=col, chr=chr, female_idx=female_idx), 0)
				marker_info['filter']=marker_info.apply(lambda row: GenerateFilterCode(marker_info=row, miss=cfg['miss'], freq=cfg['freq'], rsq=cfg['rsq'], hwe=cfg['hwe']), 1)
				marker_info['samples'] = str(cfg['data_info'][k]['samples']) + '/' + str(cfg['data_info'][k]['samples_unique']) + '/' + str(cfg['data_info'][k]['clusters']) + '/' + str(cfg['data_info'][k]['cases']) + '/' + str(cfg['data_info'][k]['ctrls'])
				markercols = [col for col in model_df.columns if 'chr' in col]
				model_df[markercols] = model_df[markercols].astype(float)
				for x in cfg['data_info'][k]['model_vars_dict'].keys():
					if cfg['data_info'][k]['model_vars_dict'][x]['class'] == 'factor':
						model_df[x] = pd.Categorical.from_array(model_df[x]).codes.astype(np.int64)
				for x in [a for a in cfg['data_info'][k]['model_vars_dict'].keys() if a != 'marker']:
					if cfg['data_info'][k]['model_vars_dict'][x]['class'] not in ['factor','random','cluster']:
						model_df[x] = model_df[x].astype(float)
				if cfg['data_info'][k]['method'].split('_')[0] in ['gee','glm','lme','coxph']:
					if cfg['data_info'][k]['method'].split('_')[0] == 'gee':
						model_df[cfg['data_info'][k]['fid']] = pd.Categorical.from_array(model_df[cfg['data_info'][k]['fid']]).codes.astype(np.int64)
						model_df.sort([cfg['data_info'][k]['fid']],inplace = True)
						results = marker_info.apply(lambda row: CalcGEE(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], dep_var=cfg['data_info'][k]['dep_var'], corstr=cfg['data_info'][k]['corstr']), 1)
					elif cfg['data_info'][k]['method'].split('_')[0] == 'glm':
						results = marker_info.apply(lambda row: CalcGLM(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], dep_var=cfg['data_info'][k]['dep_var']), 1)
					elif cfg['data_info'][k]['method'].split('_')[0] == 'lme':
						results = marker_info.apply(lambda row: CalcLME(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], dep_var=cfg['data_info'][k]['dep_var']), 1)
					elif cfg['data_info'][k]['method'] == 'coxph':
						results = marker_info.apply(lambda row: CalcCoxPH(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], dep_var=cfg['data_info'][k]['dep_var']), 1)
					results.fillna('NA', inplace=True)
					if cfg['nofail']:
						results = results[results['status'] > 0]
					if 'reg_id' in marker_list.columns and not marker_list['reg_id'][r] is None:
						results['reg_id'] = marker_list['reg_id'][r]
					if 'marker_unique' in results.columns.values:
						results.drop('marker_unique',axis=1,inplace=True)
					if not written:
						bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results.columns.values.tolist()]) + '\n')
						bgzfile.flush()
						written = True
					results.to_csv(bgzfile, header=False, sep='\t', index=False)
				else:
					if cfg['data_info'][k]['method'] in ['efftests','famskat_o']:
						if i == 1:
							reg_model_df = model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1)
							reg_marker_info = marker_info[marker_info['filter'] == 0]
						else:
							reg_model_df = pd.merge(reg_model_df,model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1),how='outer',copy=False)
							reg_marker_info = reg_marker_info.append(marker_info[marker_info['filter'] == 0],ignore_index=True)
				cur_markers = str(min(i*cfg['buffer'],marker_list[k][r])) if marker_list['start'][r] != 'NA' else str(i*cfg['buffer'])
				tot_markers = str(marker_list[k][r]) if marker_list['start'][r] != 'NA' else '> 0'
				status = '   processed ' + cur_markers + ' of ' + tot_markers + ' markers from region ' + str(r+1) + ' of ' + str(len(marker_list.index)) if len(cfg['data_info'].keys()) == 1 else '   processed ' + cur_markers + ' of ' + tot_markers + ' markers from region ' + str(r+1) + ' of ' + str(len(marker_list.index)) + " for model " + k
				print status
			if method == 'efftests':
				tot_tests = reg_marker_info.shape[0]
				if tot_tests > 0:
					if (tot_tests * tot_tests) / 100000000.0 <= mem:
						n_eff = CalcEffTests(model_df=reg_model_df[reg_marker_info['marker_unique']], mem=mem)
						status = 0
					else:
						n_eff = float('nan')
						status = 2
				else:
					n_eff, tot_tests, status = (0, 0, 1)
				results = pd.DataFrame({'chr': [marker_list['chr'][r]], 'start': [marker_list['start'][r]], 'end': [marker_list['end'][r]], 'reg_id': [marker_list['reg_id'][r]], 'n_total': [tot_tests], 'n_eff': [n_eff], 'status': [status]})[['chr','start','end','reg_id','n_total','n_eff','status']]
				results.fillna('NA', inplace=True)
				if nofail:
					results = results[results['status'] > 0]
				if 'reg_id' in marker_list.columns:
					results['reg_id'] = marker_list['reg_id'][r]
				if 'marker_unique' in results.columns.values:
					results.drop('marker_unique',axis=1,inplace=True)
				if not written:
					bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results.columns.values.tolist()]) + '\n')
					bgzfile.flush()
					written = True
				results.to_csv(bgzfile, header=False, sep='\t', index=False)
				status = '   processed effective tests calculation for region ' + str(r+1) + ' of ' + str(len(marker_list.index)) if len(cfg['data_info'].keys()) == 1 else '   processed effective tests calculation for region ' + str(r+1) + ' of ' + str(len(marker_list.index)) + " for model " + k
				print status
			elif cfg['data_info'][k]['method'] == 'famskat_o':
				nmarkers = reg_marker_info.shape[0]
				if nmarkers > 0:
					reg_model_df.dropna(inplace=True)
					snp_info = pd.DataFrame({'Name': reg_marker_info['marker_unique'], 'gene': marker_list['reg_id'][r]})
					z = reg_model_df[list(marker_info['marker_unique'])]
					pheno = reg_model_df[list(set(cfg['data_info'][k]['model_vars_dict'].keys() + [cfg['data_info'][k]['iid'],cfg['data_info'][k]['fid']]))]
					results_pre = CalcFamSkatO(snp_info=snp_info, z=z, model=cfg['data_info'][k]['model'], pheno=cfg['data_info'][k]['pheno'], pedigree=cfg['data_info'][k]['ped_df'])
					results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
					results = results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]
					if not written:
						bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results.columns.values.tolist()]) + '\n')
						bgzfile.flush()
						written = True
					results.to_csv(bgzfile, header=False, sep='\t', index=False)
		bgzfile.close()
		print "mapping results file" if len(cfg['data_info'].keys()) == 1 else "mapping results file for model " + k
		cmd = 'tabix -b 2 -e 2 ' + cfg['data_info'][k]['out']
		p = subprocess.Popen(cmd, shell=True)
		p.wait()
	print "process complete"
