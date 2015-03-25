import os
import sys
import getopt
import subprocess
import gzip
import tabix
import math
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import pandas.rpy.common as py2r
import re
from itertools import islice,takewhile,ifilter
from Bio import bgzf
import psutil
from SystemFxns import Error
from MiscFxns import *
from StatsFxns import *
from FileFxns import Coordinates
from plinkio import plinkfile
import collections
import pysam
import vcf as VCF
		
#from memory_profiler import profile, memory_usage

pd.options.mode.chained_assignment = None
ro.r('options(warn=1)')

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
			corstr = ['exchangeable'], 
			pheno_sep = ['tab'], 
			nofail = False):

	print "model options ..."
	if not cfg is None:
		for k in cfg['data_order']:
			if 'pheno_sep' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['pheno_sep'] = GetDelimiter(cfg['data_info'][k]['pheno_sep'])
			else:
				cfg['data_info'][k]['pheno_sep'] = '\t'
			if '_gaussian' in cfg['data_info'][k]['method'] or '_binomial' in cfg['data_info'][k]['method']:
				cfg['data_info'][k]['fxn'] = cfg['data_info'][k]['method'].split('_')[-1]
			else:
				cfg['data_info'][k]['fxn'] = None
			if not 'sex' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['sex'] = None
			if not 'case' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['case'] = None
			if not 'ctrl' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['ctrl'] = None
			if not 'corstr' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['corstr'] = 'exchangeable'
			if not 'pheno_sep' in cfg['data_info'][k].keys():
				cfg['data_info'][k]['pheno_sep'] = '\t'
		for arg in ['cfg','out','sig','buffer','miss','freq','rsq','hwe','nofail','region_list','region','region_id']:
			if not str(locals()[arg]) in ['None','False']:
				print "   {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])
	else:
		for arg in locals().keys():
			if not str(locals()[arg]) in ['None','False']:
				print "   {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	##### populate configuration #####
	if cfg is None:
		cfg={'out': out, 'buffer': int(buffer), 'hwe': hwe, 'data_order': ['NA'], 'freq': freq, 'miss': freq, 'rsq': freq, 'sig': int(sig), 'nofail': nofail, 
				'region': region, 'region_list': region_list, 'region_id': region_id,
				'data_info': {'NA': {'data': data[0], 'format': format[0], 'samples': samples[0], 'pheno': pheno[0], 'model': model[0], 'fid': fid[0], 'iid': iid[0],
					'method': method[0], 'focus': focus[0], 'pedigree': pedigree[0], 'sex': sex[0], 'male': male[0], 'female': female[0], 'case': case[0], 'ctrl': ctrl[0], 
						'corstr': corstr[0], 'pheno_sep': pheno_sep[0]}}}
		if 'pheno_sep' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['pheno_sep'] = GetDelimiter(cfg['data_info']['NA']['pheno_sep'])
		else:
			cfg['data_info']['NA']['pheno_sep'] = '\t'
		if '_gaussian' in cfg['data_info']['NA']['method'] or '_binomial' in cfg['data_info']['NA']['method']:
			cfg['data_info']['NA']['fxn'] = cfg['data_info']['NA']['method'].split('_')[-1]
		else:
			cfg['data_info']['NA']['fxn'] = None
		if not 'sex' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['sex'] = None
		if not 'case' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['case'] = None
		if not 'ctrl' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['ctrl'] = None
		if not 'corstr' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['corstr'] = 'exchangeable'
		if not 'pheno_sep' in cfg['data_info']['NA'].keys():
			cfg['data_info']['NA']['pheno_sep'] = '\t'

	##### DEFINE MODEL TYPES #####
	marker_tests = ['gee_gaussian','gee_binomial','glm_gaussian','glm_binomial','lme_gaussian','lme_binomial','coxph']
	seqmeta_tests = ['famskat_o','skat_o_gaussian','skat_o_binomial','famskat','skat_gaussian','skat_binomial','famburden','burden_gaussian','burden_binomial']
	efftests=['efftests']

	##### READ model VARIABLES FROM FILE #####
	vars_df_dict = {}
	model_vars_dict_dict = {}
	for k in cfg['data_order']:
		if len(cfg['data_info'].keys()) > 1:
			print "extracting model variables for model " + k + " and removing missing/invalid samples ..."
		else:
			print "extracting model variables and removing missing/invalid samples ..."
		cfg['data_info'][k]['vars_df'], cfg['data_info'][k]['model_vars_dict'] = ExtractModelVars(pheno=cfg['data_info'][k]['pheno'], model=cfg['data_info'][k]['model'],
																										fid=cfg['data_info'][k]['fid'], iid=cfg['data_info'][k]['iid'],
																										fxn=cfg['data_info'][k]['fxn'], sex=cfg['data_info'][k]['sex'],
																										case=cfg['data_info'][k]['case'], ctrl=cfg['data_info'][k]['ctrl'],
																										pheno_sep=cfg['data_info'][k]['pheno_sep'])

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
			cfg['data_info'][k]['plink_handle'],cfg['data_info'][k]['plink_locus_it'],cfg['data_info'][k]['plink_sample_it'],cfg['data_info'][k]['sample_ids'] = LoadPlink(cfg['data_info'][k]['data'])
			cfg['data_info'][k]['plink_zip'] = zip(cfg['data_info'][k]['plink_locus_it'],cfg['data_info'][k]['plink_handle'])
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
		cfg['data_info'][k]['families'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['vars_df_nodup'].duplicated(subset=[cfg['data_info'][k]['fid']])].index)
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
		print "   " + str(cfg['data_info'][k]['families']) + " clusters of size > 1"
		print "   " + str(cfg['data_info'][k]['cases']) + " case"
		print "   " + str(cfg['data_info'][k]['ctrls']) + " control"
		print "   " + str(cfg['data_info'][k]['nmale']) + " male"
		print "   " + str(cfg['data_info'][k]['nfemale']) + " female"

	##### READ PEDIGREE FROM FILE FOR FAMILY BASED SKAT TEST #####
	kinship2 = None
	for k in cfg['data_order']:
		if 'pedigree' in cfg['data_info'][k].keys():
			if kinship2 is None:
				kinship2 = importr('kinship2')
			if cfg['data_info'][k]['pedigree'] is not None:
				print "loading pedigree for model " + k if len(cfg['data_info'].keys()) > 1 else "loading pedigree"
				cfg['data_info'][k]['ped_df'] = pd.read_table(cfg['data_info'][k]['pedigree'],sep='\t',dtype='str',usecols=['FID','IID','PAT','MAT'])
				cfg['data_info'][k]['ped_df'] = cfg['data_info'][k]['ped_df'][cfg['data_info'][k]['ped_df']['IID'].isin(list(cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].values))]
				rpedigree = py2r.convert_to_r_dataframe(cfg['data_info'][k]['ped_df'], strings_as_factors=False)
				cfg['data_info'][k]['kins'] = kinship2.makekinship(rpedigree.rx2('FID'),rpedigree.rx2('IID'),rpedigree.rx2('PAT'),rpedigree.rx2('MAT'))

	##### GENERATE REGION LIST #####
	if not region_list is None:
		print "loading region list"
		marker_list = Coordinates(region_list).Load()
	elif not region is None:
		if len(region.split(':')) > 1:
			marker_list = pd.DataFrame({'chr': [re.split(':|-',region)[0]],'start': [re.split(':|-',region)[1]],'end': [re.split(':|-',region)[2]],'region': [region]})
		else:
			marker_list = pd.DataFrame({'chr': [region],'start': ['NA'],'end': ['NA'],'region': [region]})
		marker_list['reg_id'] = region_id if not region_id is None else 'NA'
	else:
		marker_list = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)],'reg_id': ['NA' for i in range(26)]})

	##### DETERMINE REGION COUNTS #####
	for k in cfg['data_order']:
		marker_list[k] = 0
		if cfg['data_info'][k]['format'] == 'plink':
			for i in range(len(marker_list.index)):
				if marker_list['start'][i] != 'NA':
					marker_list[k][i] = countPlinkRegion(cfg['data_info'][k]['plink_handle'].loci, int(marker_list['chr'][i]), int(marker_list['start'][i]), int(marker_list['end'][i]))
				else:
					marker_list[k][i] = 'n'
		else:
			for i in range(len(marker_list.index)):
				if marker_list['start'][i] != 'NA':
					try:
						records = cfg['data_info'][k]['data_it'].querys(marker_list['region'][i])
						marker_list[k][i] = countNonPlinkRegion(records, cfg['data_info'][k]['format'], int(marker_list['start'][i]), int(marker_list['end'][i]))
					except:
						marker_list[k][i] = 'n'
				else:
					marker_list[k][i] = 'n'
	marker_list = marker_list[(marker_list[cfg['data_order']].T != 0).any()].reset_index(drop=True)
	if len(marker_list.index) == 0:
		print Error("no markers found")
		return()
	else:
		print "   " + str(len(marker_list.index)) + " non-empty regions"

	##### DETERMINE ANALYSIS TYPE AND SETUP FILE HANDLES #####
	for k in cfg['data_order']:
		if cfg['data_info'][k]['method'] in marker_tests:
			cfg['data_info'][k]['method_type'] = 'marker'
		else:
			cfg['data_info'][k]['method_type'] = 'gene'
		cfg['data_info'][k]['written'] = False
	if len(list(set([cfg['data_info'][k]['method_type'] for k in cfg['data_order']]))) > 1:
		print Error("mixing gene based and marker based methods not available")
		return()
	bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')
	if 'meta' in cfg.keys():
		cfg['meta_results']={}
		cfg['meta_written']={}
		for meta in cfg['meta']:
			meta_tag = meta.split(':')[0]
			cfg['meta_written'][meta_tag]=False

	##### START ANALYSIS #####
	print "modelling data ..."

	##### LOOP OVER REGIONS #####
	for r in range(len(marker_list.index)):
		mdb = MarkerRefDb()
		reg = marker_list['region'][r]
		chr = reg.split(':')[0]

		##### LOOP OVER DATASETS #####
		for k in cfg['data_order']:
			i = 0
			cfg['reg_model_df'] = None
			cfg['reg_marker_info'] = None
			if cfg['data_info'][k]['format'] == 'plink':
				if reg == chr:
					try:
						records = (x for x in cfg['data_info'][k]['plink_zip'] if x[0].chromosome == int(chr))
					except:
						break
				else:
					try:
						records = (x for x in cfg['data_info'][k]['plink_zip'] if x[0].chromosome == int(chr) and x[0].bp_position >= int(marker_list['start'][r]) and x[0].bp_position <= int(marker_list['end'][r]))
					except:
						break
			else:
				pos_ind = 1 if cfg['data_info'][k]['format'] in ['dos2','vcf'] else 2
				try:
					records = cfg['data_info'][k]['data_it'].querys(reg)
					if reg != chr:
						records = (x for x in records if int(x[pos_ind]) >= int(marker_list['start'][r]) and int(x[pos_ind]) <= int(marker_list['end'][r]))
				except:
					break

			##### CONTINUE SLICING UNTIL NO RECORDS LEFT #####
			chunk_db = ChunkRefDb()
			while True:
				i = i + 1
				chunk=list(islice(records, cfg['buffer']))
				if not chunk:
					break
				##### READ IN CHUNK AND CREATE STANDARD DATA FRAME #####
				if cfg['data_info'][k]['format'] == 'plink':
					marker_info=collections.OrderedDict()
					marker_info['chr'] = collections.OrderedDict()
					marker_info['pos'] = collections.OrderedDict()
					marker_info['marker'] = collections.OrderedDict()
					marker_info['a1'] = collections.OrderedDict()
					marker_info['a2'] = collections.OrderedDict()
					marker_data=collections.OrderedDict()
					for locus, row in chunk:
						if locus.allele1 == '0':
							locus.allele1 = locus.allele2
						marker_unique = 'chr' + str(locus.chromosome) + 'bp' + str(locus.bp_position) + '.'  + str(locus.name) + '.' + str(locus.allele2) + '.' + str(locus.allele1)
						marker_info['chr'][marker_unique] = str(locus.chromosome)
						marker_info['marker'][marker_unique] = str(locus.name)
						marker_info['pos'][marker_unique] = str(locus.bp_position)
						marker_info['a1'][marker_unique] = str(locus.allele2)
						marker_info['a2'][marker_unique] = str(locus.allele1)
						for sample, geno in zip(cfg['data_info'][k]['plink_sample_it'], row):
							if not marker_unique in marker_data.keys():
								marker_data[marker_unique] = collections.OrderedDict({sample.iid: geno})
							else:
								marker_data[marker_unique][sample.iid] = geno if geno != 3 else 'NA'
							#if sample.iid in ['005QOL','0075FS','347']:
							#	print sample.iid, geno
					marker_info = pd.DataFrame(marker_info)
					marker_data = pd.DataFrame(marker_data)
					marker_data = marker_data.convert_objects(convert_numeric=True)
					marker_data[cfg['data_info'][k]['iid']] = marker_data.index
					chunkdf = marker_info.join(marker_data.transpose())
				else:
					if cfg['data_info'][k]['format'] == 'vcf':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf.apply(GetRowCalls,1)
						chunkdf.drop(chunkdf.columns[5:9],axis=1,inplace=True)
					elif cfg['data_info'][k]['format'] == 'oxford':
						chunk = [ConvertDosage(row) for row in chunk]
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['data_info'][k]['format'] == 'dos1':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['data_info'][k]['format'] == 'dos2':
						chunkdf = pd.DataFrame(chunk)
					chunkdf.columns = ['chr','pos','marker','a1','a2'] + cfg['data_info'][k]['sample_ids']
				chunkdf = chunkdf.drop_duplicates(subset=['chr','pos','a1','a2'])
				chunkdf.index=chunkdf['chr'].astype(str) + '><' + chunkdf['pos'].astype(str) + '><'  + chunkdf['a1'].astype(str) + '><'  + chunkdf['a2'].astype(str)
				chunkdf = chunkdf[~chunkdf.index.isin(chunk_db.ListKeys())]
				chunkdf.apply(chunk_db.Update,1)
				chunkdf = chunkdf.apply(mdb.Update,1)

				##### EXTRACT MARKER INFO AND MARKER DATA #####
				marker_info = chunkdf.ix[:,:5]
				marker_info.replace('.','NA',inplace=True)
				marker_info['marker_unique'] = 'chr' + marker_info['chr'].astype(str) + 'bp' + marker_info['pos'].astype(str) + '.'  + marker_info['marker'].astype(str).str.replace('-','_') + '.'  + marker_info['a1'].astype(str) + '.'  + marker_info['a2'].astype(str)
				marker_info.index = marker_info['marker_unique']
				marker_data = chunkdf.ix[:,5:].transpose()
				marker_data = marker_data.convert_objects(convert_numeric=True)
				marker_data.columns = marker_info['marker_unique']
				marker_data[cfg['data_info'][k]['iid']] = marker_data.index

				##### MERGE PHENO DATAFRAME AND MARKER DATA #####
				model_df = pd.merge(cfg['data_info'][k]['vars_df'], marker_data, on = [cfg['data_info'][k]['iid']], how='left').sort([cfg['data_info'][k]['fid']])

				##### EXTRACT UNIQUE SAMPLES AND CALCULATE DESCRIPTIVE STATS #####
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
				if cfg['data_info'][k]['method'] in seqmeta_tests:
					marker_info['filter']=marker_info.apply(lambda row: GenerateFilterCode(marker_info=row, no_mono=True, miss=cfg['miss'], freq=cfg['freq'], rsq=cfg['rsq'], hwe=cfg['hwe']), 1)
				else:
					marker_info['filter']=marker_info.apply(lambda row: GenerateFilterCode(marker_info=row, miss=cfg['miss'], freq=cfg['freq'], rsq=cfg['rsq'], hwe=cfg['hwe']), 1)
				marker_info['samples'] = str(cfg['data_info'][k]['samples']) + '/' + str(cfg['data_info'][k]['samples_unique']) + '/' + str(cfg['data_info'][k]['clusters']) + '/' + str(cfg['data_info'][k]['cases']) + '/' + str(cfg['data_info'][k]['ctrls'])

				##### CONVERT ALL COLUMNS TO APPROPRIATE FORMAT FOR ANALYSIS #####
				markercols = [col for col in model_df.columns if 'chr' in col]
				model_df[markercols] = model_df[markercols].astype(float)
				for x in cfg['data_info'][k]['model_vars_dict'].keys():
					if cfg['data_info'][k]['model_vars_dict'][x]['class'] == 'factor':
						model_df[x] = pd.Categorical.from_array(model_df[x]).codes.astype(np.int64)
				for x in [a for a in cfg['data_info'][k]['model_vars_dict'].keys() if a != 'marker']:
					if cfg['data_info'][k]['model_vars_dict'][x]['class'] not in ['factor','random','cluster']:
						model_df[x] = model_df[x].astype(float)

				##### MARKER ANALYSIS #####
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
					if cfg['data_info'][k]['written'] == False:
						results.columns = [k + '.' + a if not a in ['chr','pos','a1','a2'] and k != 'NA' else a for a in results.columns]
						cfg['data_info'][k]['results'] = results.copy()
						cfg['data_info'][k]['written'] = True
					else:
						results.columns = [k + '.' + a if not a in ['chr','pos','a1','a2'] and k != 'NA' else a for a in results.columns]
						cfg['data_info'][k]['results'] = cfg['data_info'][k]['results'].append(results).reset_index(drop=True)

				##### PREPARE DATA FOR GENE BASED ANALYSIS #####
				else:
					if cfg['data_info'][k]['method_type'] == 'gene':
						if i == 1:
							cfg['reg_model_df'] = model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1)
							cfg['reg_marker_info'] = marker_info[marker_info['filter'] == 0]
						else:
							cfg['reg_model_df'] = pd.merge(cfg['reg_model_df'],model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1),how='outer',copy=False)
							cfg['reg_marker_info'] = cfg['reg_marker_info'].append(marker_info[marker_info['filter'] == 0],ignore_index=True)

				##### UPDATE LOOP STATUS #####
				cur_markers = str(min(i*cfg['buffer'],marker_list[k][r])) if marker_list['start'][r] != 'NA' else str(i*cfg['buffer'])
				tot_markers = str(marker_list[k][r]) if marker_list['start'][r] != 'NA' else '> 0'
				status_reg = marker_list['region'][r] + ': ' + marker_list['reg_id'][r] if 'reg_id' in marker_list.columns and marker_list['reg_id'][r] != 'NA' else marker_list['region'][r]
				status = '   processed ' + cur_markers + '/' + tot_markers + ' markers from region ' + str(r+1) + '/' + str(len(marker_list.index)) + ' (' + status_reg + ')' if len(cfg['data_info'].keys()) == 1 else '   processed ' + cur_markers + '/' + tot_markers + ' markers from region ' + str(r+1) + '/' + str(len(marker_list.index)) + ' (' + status_reg + ')' + " for model " + str(cfg['data_order'].index(k) + 1) + '/' + str(len(cfg['data_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE EFFECTIVE TESTS #####
			if cfg['data_info'][k]['method'] == 'efftests':
				tot_tests = cfg['reg_marker_info'].shape[0]
				if tot_tests > 0:
					n_eff = CalcEffTests(model_df=cfg['reg_model_df'][cfg['reg_marker_info']['marker_unique']])
					status = 0
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
				if cfg['data_info'][k]['written'] == False:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = results.copy()
					cfg['data_info'][k]['written'] = True
				else:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = cfg['data_info'][k]['results'].append(results).reset_index(drop=True)
				status = '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(marker_list.index)) if len(cfg['data_info'].keys()) == 1 else '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(marker_list.index)) + " for model " + str(cfg['data_order'].index(k) + 1) + '/' + str(len(cfg['data_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE PREPSCORES AND INDIVIDUAL MODELS FOR GENE BASED TESTS #####
			elif cfg['data_info'][k]['method'] in seqmeta_tests:
				if not cfg['reg_marker_info'] is None and cfg['reg_marker_info'].shape[0] > 1:
					cfg['data_info'][k]['snp_info'] = pd.DataFrame({'Name': cfg['reg_marker_info']['marker_unique'], 'gene': marker_list['reg_id'][r]})
					z = cfg['reg_model_df'][list(cfg['reg_marker_info']['marker_unique'][cfg['reg_marker_info']['filter'] == 0])]
					pheno = cfg['reg_model_df'][list(set(cfg['data_info'][k]['model_vars_dict'].keys() + [cfg['data_info'][k]['iid'],cfg['data_info'][k]['fid']]))]
					cfg['data_info'][k]['snp_info'].to_csv("test.snpinfo." + k, index=True, header=True, sep="\t")
					z.to_csv("test.z." + k, index=True, header=True, sep="\t")
					pheno.to_csv("test.pheno." + k, index=True, header=True, sep="\t")
					
					##### PREPSCORES #####
					if cfg['data_info'][k]['method'] in ['famskat_o','famskat','famburden']:
						ro.globalenv['ps' + k] = PrepScoresFam(snp_info=cfg['data_info'][k]['snp_info'], z=z, model=cfg['data_info'][k]['model'], pheno=pheno, kinship=cfg['data_info'][k]['kins'])
					else:
						ro.globalenv['ps' + k] = PrepScores(snp_info=cfg['data_info'][k]['snp_info'], z=z, model=cfg['data_info'][k]['model'], pheno=pheno, family=cfg['data_info'][k]['fxn'])

					##### SKAT-O #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results_pre = SkatOMeta('skatOMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'])
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results_pre = SkatMeta('skatMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'])
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'Qmeta': [results_pre['pmin'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]]})
						results = results[['chr','start','end','reg_id','p','Qmeta','cmaf','nmiss','nsnps']]

					##### BURDEN #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results_pre = BurdenMeta('burdenMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'])
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','reg_id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SNP_INFO DF #####
					cfg['data_info'][k]['snp_info'] = pd.DataFrame({'Name': [], 'gene': []})

					##### EMPTY SKAT-O DF #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results = SkatOMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])

					##### EMPTY SKAT DF #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results = SkatMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])

					##### EMPTY BURDEN DF #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results = BurdenMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])
					
				##### APPEND TO RESULTS DF #####
				if cfg['data_info'][k]['written'] == False:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = results.copy()
					cfg['data_info'][k]['written'] = True
				else:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = cfg['data_info'][k]['results'].append(results).reset_index(drop=True)

		##### PERFORM GENE BASED TEST META ANALYSIS #####
		if cfg['data_info'][cfg['data_order'][0]]['method'] in seqmeta_tests and 'meta' in cfg.keys():
			for meta in cfg['meta']:
				meta_tag = meta.split(':')[0]
				meta_incl = meta.split(':')[1].split('+')
				seqmeta_cmd = ''
				snp_info_meta = None
				for k in meta_incl:
					if 'ps' + k in ro.globalenv and cfg['data_info'][k]['results'][k + '.p'][0] != 'NA':
						if snp_info_meta is None:
							snp_info_meta = cfg['data_info'][k]['snp_info']
						else:
							snp_info_meta = snp_info_meta.merge(cfg['data_info'][k]['snp_info'], how='outer')
						if seqmeta_cmd == '':
							seqmeta_cmd = 'ps' + k
						else:
							seqmeta_cmd = seqmeta_cmd + ', ps' + k
				if snp_info_meta is not None and snp_info_meta.shape[0] > 1:

					##### SKAT-O #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results_pre = SkatOMeta('skatOMeta(' + seqmeta_cmd + ', SNPInfo=rsnp_info)', snp_info_meta)
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
													'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
													'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results_pre = SkatMeta('skatMeta(' + seqmeta_cmd + ', SNPInfo=rsnp_info)', snp_info_meta)
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
													'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
													'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### BURDEN #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results_pre = BurdenMeta('burdenMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'])
						results = pd.DataFrame({'chr': [marker_list['chr'][r]],'start': [marker_list['start'][r]],'end': [marker_list['end'][r]],'reg_id': [marker_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','reg_id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SKAT-O DF #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results = SkatOMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])

					##### EMPTY SKAT DF #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results = SkatMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])

					##### EMPTY BURDEN DF #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results = BurdenMetaEmpty(marker_list['chr'][r],marker_list['start'][r],marker_list['end'][r],marker_list['reg_id'][r])

				##### APPEND TO META DF #####
				if cfg['meta_written'][meta_tag] == False:
					results.columns = [meta_tag + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['meta_results'][meta_tag] = results.copy()
					cfg['meta_written'][meta_tag] = True
				else:
					results.columns = [meta_tag + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['meta_results'][meta_tag] = cfg['meta_results'][meta_tag].merge(results, how='outer', copy=False)

	##### COMPILE ALL RESULTS #####
	header = ['chr','start','end','reg_id'] if list(set([cfg['data_info'][k]['method_type'] for k in cfg['data_order']]))[0] == 'gene' else ['chr','pos','a1','a2']
	if 'meta' in cfg.keys():
		for meta in cfg['meta']:
			meta_tag = meta.split(':')[0]
			if meta == cfg['meta'][0]:
				results = cfg['meta_results'][meta_tag]
			else:
				results = results.merge(cfg['meta_results'][meta_tag], how='outer', copy=False)
			header = header + [a for a in cfg['meta_results'][meta_tag].columns.values.tolist() if not a in header]
	for k in cfg['data_order']:
		if 'reg_id' in cfg['data_info'][k]['results'].keys() and len(list(cfg['data_info'][k]['results']['reg_id'].unique())) and list(cfg['data_info'][k]['results']['reg_id'].unique())[0] == 'NA':
			cfg['data_info'][k]['results'].drop('reg_id',axis=1,inplace=True)
		if k + '.reg_id' in cfg['data_info'][k]['results'].keys() and len(list(cfg['data_info'][k]['results'][k + '.reg_id'].unique())) and list(cfg['data_info'][k]['results'][k + '.reg_id'].unique())[0] == 'NA':
			cfg['data_info'][k]['results'].drop(k + '.reg_id',axis=1,inplace=True)
		if cfg['data_info'][k]['fxn'] is not None and cfg['data_info'][k]['fxn'] == 'gaussian':
			cfg['data_info'][k]['results'].drop([x for x in cfg['data_info'][k]['results'].columns if '.or' in x], axis=1,inplace=True)
		if k == cfg['data_order'][0] and 'meta' not in cfg.keys():
			results = cfg['data_info'][k]['results']
		else:
			results = results.merge(cfg['data_info'][k]['results'], how='outer', copy=False)
		header = header + [a for a in cfg['data_info'][k]['results'].columns.values.tolist() if not a in header]

	##### SORT RESULTS AND CONVERT CHR, POS, AND START COLUMNS TO INT #####
	if list(set([cfg['data_info'][k]['method_type'] for k in cfg['data_order']]))[0] == 'marker':
		results[['chr','pos']] = results[['chr','pos']].astype(int)
		results.sort(columns=['chr','pos'],inplace=True)
	else:
		results[['chr','start','end']] = results[['chr','start','end']].astype(int)
		results.sort(columns=['chr','start'],inplace=True)

	##### FILL IN NA's, ORDER HEADER, AND WRITE TO FILE #####
	results.fillna('NA', inplace=True)	
	results = results[header]
	bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results.columns.values.tolist()]) + '\n')
	bgzfile.flush()
	results.to_csv(bgzfile, header=False, sep='\t', index=False)
	bgzfile.close()

	##### MAP OUTPUT FILES FOR TABIX #####
	if list(set([cfg['data_info'][k]['method_type'] for k in cfg['data_order']]))[0] == 'marker':
		cmd = ['tabix','-b','2','-e','2',cfg['out'] + '.gz']
	else:
		cmd = ['tabix','-b','2','-e','3',cfg['out'] + '.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessError:
		print Error("file mapping failed")
	else:
		print "process complete"
