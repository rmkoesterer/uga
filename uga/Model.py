#!/usr/bin/env python
from __main__ import *
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
import pandas.rpy.common as py2r
from itertools import islice,takewhile,ifilter
from Bio import bgzf
import collections
from scipy.stats import norm,t
import vcf as VCF
from plinkio import plinkfile
import MiscFxns
import MiscFxnsCy
import StatsFxns

pandas2ri.activate()
ro.r('options(warn=1)')

def Model(cfg):
	Parse.PrintModelOptions(cfg)
	for k in cfg['models']:
		cfg['models'][k]['sep'] = MiscFxnsCy.GetDelimiterCy(cfg['models'][k]['sep'])

	##### DEFINE MODEL TYPES #####
	marker_tests = ['bgee','ggee','bglm','gglm','blme','glme','cph']
	geepack_tests = ['bgee','ggee']
	lme_tests = ['blme','glme']
	lme_gaussian_tests = ['glme']
	survival_tests = ['cph']
	seqmeta_tests = ['fskato','gskato','bskato','fskat','gskat','bskat','fburden','gburden','bburden']
	seqmeta_famtests = ['fskato','fskat','fburden']
	efftests=['efftests']
	binomial_tests = ['bgee','bglm','blme','bskato','bskat','bburden']

	##### STOP IF MODEL TYPES IN META ANALYSIS ARE INCOMPATIBLE #####
	for x in cfg['meta']:
		if not len(set([cfg['models'][k]['model_fxn'] for k in x.split(':')[1].split('+')])) == 1 or not set([cfg['models'][k]['model_fxn'] for k in x.split(':')[1].split('+')]) & set(seqmeta_tests):
			print SystemFxns.Error("only matching seqMeta models may be meta-analyzed")
			return()

	##### LOAD NECESSARY R PACKAGES #####
	for m in list(set([cfg['models'][k]['model_fxn'] for k in cfg['model_order']])):
		if m in geepack_tests:
			importr('geepack')
		elif m in lme_tests:
			importr('lmerTest')
			if m in lme_gaussian_tests:
				importr('pbkrtest')
		elif m in survival_tests:
			importr('survival')
		elif m in seqmeta_tests:
			importr('seqMeta')
			if m in seqmeta_famtests:
				importr('kinship2')

	##### READ model VARIABLES FROM FILE #####
	for k in cfg['model_order']:
		if len(cfg['models'].keys()) > 1:
			print "extracting model variables for model " + k + " and removing missing/invalid samples ..."
		else:
			print "extracting model variables and removing missing/invalid samples ..."
		cfg['models'][k]['vars_df'], cfg['models'][k]['model_vars_dict'] = MiscFxns.ExtractModelVars(pheno=cfg['models'][k]['pheno'], model=cfg['models'][k]['model'],
																										fid=cfg['models'][k]['fid'], iid=cfg['models'][k]['iid'],
																										family=cfg['models'][k]['family'], 
																										sex=cfg['models'][k]['sex'], case=cfg['models'][k]['case'], 
																										ctrl=cfg['models'][k]['ctrl'], sep=cfg['models'][k]['sep'])

	##### DETERMINE MODEL STATS TO BE EXTRACTED #####
	for k in cfg['model_order']:
		if 'focus' in cfg['models'][k].keys() and not cfg['models'][k]['focus'] is None:
			cfg['models'][k]['focus'] = cfg['models'][k]['focus'].split(',')
		else:
			cfg['models'][k]['focus'] = MiscFxns.GetFocus(model_fxn=cfg['models'][k]['model_fxn'],model=cfg['models'][k]['model'],vars_df=cfg['models'][k]['vars_df'],model_vars_dict=cfg['models'][k]['model_vars_dict'])

	##### LOAD ITERATORS AND SAMPLE LISTS #####
	for k in cfg['model_order']:
		if cfg['models'][k]['format'] == 'plink':
			print "loading Plink binary files for model " + k if len(cfg['models'].keys()) > 1 else "loading Plink binary files"
			cfg['models'][k]['plink_handle'],cfg['models'][k]['plink_locus_it'],cfg['models'][k]['plink_sample_it'],cfg['models'][k]['sample_ids'] = MiscFxns.LoadPlink(cfg['models'][k]['data'])
			cfg['models'][k]['plink_zip'] = zip(cfg['models'][k]['plink_locus_it'],cfg['models'][k]['plink_handle'])
		elif cfg['models'][k]['format'] == 'vcf':
			print "loading vcf file for model " + k if len(cfg['models'].keys()) > 1 else "loading vcf file"
			cfg['models'][k]['data_it'], cfg['models'][k]['sample_ids'] = MiscFxns.LoadVcf(cfg['models'][k]['data'])
		else:
			print "loading data and sample files for model " + k if len(cfg['models'].keys()) > 1 else "loading data and sample files"
			cfg['models'][k]['data_it'], cfg['models'][k]['sample_ids'] = MiscFxns.LoadDos(cfg['models'][k]['data'], cfg['models'][k]['sample'])
		if len(cfg['models'][k]['vars_df'][cfg['models'][k]['vars_df'][cfg['models'][k]['iid']].isin(cfg['models'][k]['sample_ids'])]) == 0:
			print SystemFxns.Error("phenotype file and data file contain no common samples")
			return
		cfg['models'][k]['varlist_handle'] = None
		if 'varlist' in cfg['models'][k]:
			if cfg['models'][k]['varlist'] is not None:
				cfg['models'][k]['varlist_handle'] = tabix.open(cfg['models'][k]['varlist'])

	##### REDUCE PHENO DATA TO GENOTYPE IDs #####
	for k in cfg['model_order']:
		cfg['models'][k]['vars_df'] = cfg['models'][k]['vars_df'][cfg['models'][k]['vars_df'][cfg['models'][k]['iid']].isin(cfg['models'][k]['sample_ids'])]

	##### CREATE SUMMARY FOR DATA TO BE ANALYZED #####
	for k in cfg['model_order']:
		cfg['models'][k]['sample'] = len(cfg['models'][k]['vars_df'].index)
		cfg['models'][k]['vars_df_nodup'] = cfg['models'][k]['vars_df'].drop_duplicates(subset=[cfg['models'][k]['iid']])
		cfg['models'][k]['samples_unique'] = len(cfg['models'][k]['vars_df_nodup'].index)
		cfg['models'][k]['clusters'] = len(cfg['models'][k]['vars_df_nodup'].drop_duplicates(subset=[cfg['models'][k]['fid']]).index)
		cfg['models'][k]['families'] = len(cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['vars_df_nodup'].duplicated(subset=[cfg['models'][k]['fid']])][cfg['models'][k]['fid']].unique())
		cfg['models'][k]['dep_var'] = [v for v in cfg['models'][k]['model_vars_dict'] if cfg['models'][k]['model_vars_dict'][v]['type'] == 'dependent']
		cfg['models'][k]['cases'] = len(cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['dep_var'][0]][cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['dep_var'][0]].isin(['1'])]) if cfg['models'][k]['family'] == 'binomial' else 'NA'
		cfg['models'][k]['ctrls'] = len(cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['dep_var'][0]][cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['dep_var'][0]].isin(['0'])]) if cfg['models'][k]['family'] == 'binomial' else 'NA'
		cfg['models'][k]['nmale'] = len(cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['sex']].isin([str(cfg['models'][k]['male'])])].index.values) if not cfg['models'][k]['sex'] is None and not cfg['models'][k]['male'] is None and not cfg['models'][k]['female'] is None else 'NA'
		cfg['models'][k]['nfemale'] = len(cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['vars_df_nodup'][cfg['models'][k]['sex']].isin([str(cfg['models'][k]['female'])])].index.values) if not cfg['models'][k]['sex'] is None and not cfg['models'][k]['male'] is None and not cfg['models'][k]['female'] is None else 'NA'
		if len(cfg['models'].keys()) > 1:
			print "data summary for model " + k + " ..."
		else:
			print "data summary ..."
		print "   " + str(cfg['models'][k]['sample']) + " total observations"
		print "   " + str(cfg['models'][k]['samples_unique']) + " unique samples"
		print "   " + str(cfg['models'][k]['clusters']) + " clusters"
		print "   " + str(cfg['models'][k]['families']) + " clusters of size > 1"
		print "   " + str(cfg['models'][k]['cases']) + " case"
		print "   " + str(cfg['models'][k]['ctrls']) + " control"
		print "   " + str(cfg['models'][k]['nmale']) + " male"
		print "   " + str(cfg['models'][k]['nfemale']) + " female"

	##### READ PEDIGREE FROM FILE FOR FAMILY BASED SKAT TEST #####
	for k in cfg['model_order']:
		if 'pedigree' in cfg['models'][k].keys():
			if cfg['models'][k]['pedigree'] is not None:
				print "loading pedigree for model " + k if len(cfg['models'].keys()) > 1 else "loading pedigree"
				cfg['models'][k]['ped_df'] = pd.read_table(cfg['models'][k]['pedigree'],sep='\t',dtype='str',usecols=['FID','IID','PAT','MAT'])
				cfg['models'][k]['ped_df'] = cfg['models'][k]['ped_df'][cfg['models'][k]['ped_df']['IID'].isin(list(cfg['models'][k]['vars_df'][cfg['models'][k]['iid']].values))]

	##### GENERATE REGION LIST #####
	if not cfg['reglist'] is None:
		print "loading region list"
		reglist = FileFxns.LoadCoordinates(cfg['reglist'])
	elif not cfg['region'] is None:
		if len(cfg['region'].split(':')) > 1:
			reglist = pd.DataFrame({'chr': [re.split(':|-',cfg['region'])[0]],'start': [re.split(':|-',cfg['region'])[1]],'end': [re.split(':|-',cfg['region'])[2]],'region': [cfg['region']]})
		else:
			reglist = pd.DataFrame({'chr': [cfg['region']],'start': ['NA'],'end': ['NA'],'region': [cfg['region']]})
		reglist['id'] = cfg['id'] if not cfg['id'] is None else 'NA'
	else:
		reglist = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)],'id': ['NA' for i in range(26)]})

	##### DETERMINE ANALYSIS TYPE AND SETUP FILE HANDLES #####
	for k in cfg['model_order']:
		if cfg['models'][k]['model_fxn'] in marker_tests:
			cfg['models'][k]['model_fxn_type'] = 'marker'
		else:
			cfg['models'][k]['model_fxn_type'] = 'gene'
		cfg['models'][k]['written'] = False
	if len(list(set([cfg['models'][k]['model_fxn_type'] for k in cfg['model_order']]))) > 1:
		print SystemFxns.Error("mixing gene based and marker based model_fxns not available")
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
	for r in range(len(reglist.index)):
		mdb = MiscFxns.MarkerRefDb()
		reg = reglist['region'][r]
		chr = int(reg.split(':')[0])

		##### LOOP OVER DATASETS #####
		for k in cfg['model_order']:
			i = 0
			cfg['reg_model_df'] = None
			cfg['reg_marker_info'] = None
			varlist = None
			if cfg['models'][k]['varlist_handle'] is not None:
				varlist = cfg['models'][k]['varlist_handle'].querys(reg)
				varlist = pd.DataFrame([x for x in varlist if int(x[1]) >= int(reglist['start'][r]) and int(x[1]) <= int(reglist['end'][r])])
				varlist.columns = ['chr','pos','marker','a1','a2']
				varlist['chr'] = varlist['chr'].astype(np.int64)
				varlist['pos'] = varlist['pos'].astype(np.int64)
			if cfg['models'][k]['format'] == 'plink':
				if reg == str(chr):
					try:
						records = (x for x in cfg['models'][k]['plink_zip'] if x[0].chromosome == chr and not ',' in x[0].allele2)
					except:
						break
				else:
					try:
						records = (x for x in cfg['models'][k]['plink_zip'] if x[0].chromosome == chr and x[0].bp_position >= int(reglist['start'][r]) and x[0].bp_position <= int(reglist['end'][r]) and not ',' in x[0].allele2)
					except:
						break
				if varlist is not None:
					records = (x for x in records if varlist[(varlist['chr'] == x[0].chromosome) & (varlist['pos'] == x[0].bp_position) & (varlist['a1'] == x[0].allele1) & (varlist['a2'] == x[0].allele2)].shape[0] > 0)
			else:
				pos_ind = 1 if cfg['models'][k]['format'] in ['dos2','vcf'] else 2
				try:
					records = cfg['models'][k]['data_it'].querys(reg)
					if reg != str(chr):
						records = (x for x in records if int(x[pos_ind]) >= int(reglist['start'][r]) and int(x[pos_ind]) <= int(reglist['end'][r]) and not ',' in str(x[4]))
				except:
					break
				if varlist is not None:
					records = (x for x in records if varlist[(varlist['chr'] == int(x[0])) & (varlist['pos'] == int(x[pos_ind])) & (varlist['a1'] == x[3]) & (varlist['a2'] == x[4])].shape[0] > 0)

			##### CONTINUE SLICING UNTIL NO RECORDS LEFT #####
			chunk_db = MiscFxns.ChunkRefDb()

			while True:
				i = i + 1
				chunk=list(islice(records, cfg['buffer']))
				if not chunk:
					break

				##### READ IN CHUNK AND CREATE STANDARD DATA FRAME #####
				if cfg['models'][k]['format'] == 'plink':
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
						for sample, geno in zip(cfg['models'][k]['plink_sample_it'], row):
							if not marker_unique in marker_data.keys():
								marker_data[marker_unique] = collections.OrderedDict({sample.iid: geno})
							else:
								marker_data[marker_unique][sample.iid] = geno if geno != 3 else 'NA'
					marker_info = pd.DataFrame(marker_info)
					marker_data = pd.DataFrame(marker_data)
					marker_data = marker_data.convert_objects(convert_numeric=True)
					marker_data[cfg['models'][k]['iid']] = marker_data.index
					chunkdf = marker_info.join(marker_data.transpose())
				else:
					if cfg['models'][k]['format'] == 'vcf':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf.apply(MiscFxnsCy.GetRowCallsCy,axis=1,raw=True)
						chunkdf.drop(chunkdf.columns[5:9],axis=1,inplace=True)
					elif cfg['models'][k]['format'] == 'oxford':
						chunk = [MiscFxns.ConvertDosage(row) for row in chunk]
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['models'][k]['format'] == 'dos1':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['models'][k]['format'] == 'dos2':
						chunkdf = pd.DataFrame(chunk)
					chunkdf.columns = ['chr','pos','marker','a1','a2'] + cfg['models'][k]['sample_ids']
				if len(cfg['model_order']) > 1:
					chunkdf = chunkdf.drop_duplicates(subset=['chr','pos','a1','a2'])
					chunkdf.index=chunkdf['chr'].astype(str) + '><' + chunkdf['pos'].astype(str) + '><'  + chunkdf['a1'].astype(str) + '><'  + chunkdf['a2'].astype(str)
					chunkdf = chunkdf[~chunkdf.index.isin(chunk_db.ListKeys())]
					chunkdf.apply(chunk_db.Update,1)
					chunkdf = chunkdf.apply(mdb.Update,1)
				
				##### EXTRACT MARKER INFO AND MARKER DATA #####
				marker_info = chunkdf.ix[:,:5]
				marker_info.replace('.','NA',inplace=True)
				marker_info['marker_unique'] = 'chr' + marker_info['chr'].astype(str) + 'bp' + marker_info['pos'].astype(str) + '_'  + marker_info['marker'].astype(str).str.replace('-','_').str.replace(':','.').str.replace(';','.') + '_'  + marker_info['a1'].astype(str).str.replace('-','_') + '_'  + marker_info['a2'].astype(str).str.replace('-','_')
				marker_info.index = marker_info['marker_unique']
				marker_data = chunkdf.ix[:,5:].transpose()
				marker_data.replace('NA',np.nan,inplace=True)
				marker_data.columns = marker_info['marker_unique']
				marker_data[cfg['models'][k]['iid']] = marker_data.index

				##### MERGE PHENO DATAFRAME AND MARKER DATA #####
				model_df = pd.merge(cfg['models'][k]['vars_df'], marker_data, on = [cfg['models'][k]['iid']], how='left').sort([cfg['models'][k]['fid']])

				##### EXTRACT UNIQUE SAMPLES AND CALCULATE DESCRIPTIVE STATS #####
				model_df_nodup=model_df.drop_duplicates(subset=[cfg['models'][k]['iid']]).reset_index(drop=True)

				marker_info['callrate']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcCallrateCy(np.array(col).astype(np.float)), raw=True)
				model_df_nodup_unrel=model_df_nodup.drop_duplicates(subset=[cfg['models'][k]['fid']]).reset_index(drop=True)
				if not cfg['models'][k]['sex'] is None and not cfg['models'][k]['male'] is None and not cfg['models'][k]['female']:
					male_idx = model_df_nodup[model_df_nodup[cfg['models'][k]['sex']].isin([cfg['models'][k]['male']])].index.values
					female_idx = model_df_nodup[model_df_nodup[cfg['models'][k]['sex']].isin([cfg['models'][k]['female']])].index.values
				else:
					male_idx = None
					female_idx = None
				if chr == 23:
					marker_info['freq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreq23Cy(marker_male=np.array(col[male_idx]).astype(np.float), marker_female=np.array(col[female_idx]).astype(np.float)), raw=True)
					marker_info['freq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreq23Cy(marker_male=np.array(col[male_idx]).astype(np.float), marker_female=np.array(col[female_idx]).astype(np.float)), raw=True)
				else:
					marker_info['freq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['freq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
				if cfg['models'][k]['family'] == 'binomial':
					marker_info['freq.ctrl']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['freq.case']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['freq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['freq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(np.array(col).astype(np.float)), raw=True)
				marker_info['rsq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
				marker_info['rsq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
				if cfg['models'][k]['family'] == 'binomial':
					marker_info['rsq.ctrl']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['rsq.case']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['rsq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
					marker_info['rsq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(np.array(col).astype(np.float)), raw=True)
				marker_info['hwe']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
				marker_info['hwe.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
				if cfg['models'][k]['family'] == 'binomial':
					marker_info['hwe.ctrl']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
					marker_info['hwe.case']=model_df_nodup[model_df_nodup[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
					marker_info['hwe.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
					marker_info['hwe.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['models'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(np.array(col).astype(np.float)), raw=True)
				if cfg['models'][k]['model_fxn'] in seqmeta_tests:
					marker_info['filter']=marker_info.apply(lambda row: MiscFxnsCy.GenerateFilterCodeCy(row_callrate=row['callrate'], no_mono=False, 
																											row_freq=row['freq'], row_freq_unrel=row['freq.unrel'], 
																											row_rsq=row['rsq'], row_rsq_unrel=row['rsq.unrel'], 
																											row_hwe=row['hwe'], row_hwe_unrel=row['hwe.unrel'], 
																											miss_thresh=cfg['models'][k]['miss'], maf_thresh=cfg['models'][k]['maf'], 
																											maxmaf_thresh=cfg['models'][k]['maxmaf'], rsq_thresh=cfg['models'][k]['rsq'], 
																											hwe_thresh=cfg['models'][k]['hwe']), axis=1, raw=True)
				else:
					marker_info['filter']=marker_info.apply(lambda row: MiscFxnsCy.GenerateFilterCodeCy(row_callrate=row['callrate'], 
																											row_freq=row['freq'], row_freq_unrel=row['freq.unrel'], 
																											row_rsq=row['rsq'], row_rsq_unrel=row['rsq.unrel'], 
																											row_hwe=row['hwe'], row_hwe_unrel=row['hwe.unrel'], 
																											miss_thresh=cfg['models'][k]['miss'], maf_thresh=cfg['models'][k]['maf'], 
																											maxmaf_thresh=cfg['models'][k]['maxmaf'], rsq_thresh=cfg['models'][k]['rsq'], 
																											hwe_thresh=cfg['models'][k]['hwe']), axis=1, raw=True)
				marker_info['sample'] = str(cfg['models'][k]['sample']) + '/' + str(cfg['models'][k]['samples_unique']) + '/' + str(cfg['models'][k]['clusters']) + '/' + str(cfg['models'][k]['cases']) + '/' + str(cfg['models'][k]['ctrls']) + '/' + str(cfg['models'][k]['nmale']) + '/' + str(cfg['models'][k]['nfemale'])

				##### CONVERT ALL COLUMNS TO APPROPRIATE FORMAT FOR ANALYSIS #####
				markercols = [col for col in model_df.columns if 'chr' in col]
				model_df[markercols] = model_df[markercols].astype(float)
				for x in cfg['models'][k]['model_vars_dict'].keys():
					if cfg['models'][k]['model_vars_dict'][x]['class'] == 'factor':
						model_df[x] = pd.Categorical.from_array(model_df[x]).codes.astype(np.int64)
				for x in [a for a in cfg['models'][k]['model_vars_dict'].keys() if a != 'marker']:
					if cfg['models'][k]['model_vars_dict'][x]['class'] not in ['factor','random','cluster']:
						model_df[x] = model_df[x].astype(float)

				##### MARKER ANALYSIS #####
				if cfg['models'][k]['model_fxn'] in ['bgee','ggee','bglm','gglm','blme','glme','cph']:
					marker_info_passed = marker_info[marker_info['filter'] == 0]
					marker_info_filtered = marker_info[marker_info['filter'] > 0]
					marker_info_filtered['n'] = 0
					marker_info_filtered['status'] = -1
					model_df = model_df[[c for c in model_df.columns if not c in marker_info_filtered['marker_unique']]]
					model_df[cfg['models'][k]['fid']] = pd.Categorical.from_array(model_df[cfg['models'][k]['fid']]).codes.astype(np.int64)
					model_df.sort([cfg['models'][k]['fid']],inplace = True)
					ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(model_df, strings_as_factors=False)
					ro.r('model_df[,names(model_df) == "' + cfg['models'][k]['fid'] + '"]<-as.factor(model_df[,names(model_df) == "' + cfg['models'][k]['fid'] + '"])')
					if cfg['models'][k]['model_fxn'] in ['bgee','ggee']:
						for x in cfg['models'][k]['focus']:
							marker_info_passed[x + '.effect'] = float('nan')
							marker_info_passed[x + '.stderr'] = float('nan')
							marker_info_passed[x + '.or'] = float('nan')
							marker_info_passed[x + '.z'] = float('nan')
							marker_info_passed[x + '.p'] = float('nan')
						marker_info_passed['n'] = 0
						marker_info_passed['status'] = 0
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcGEE(marker_info=row, model_df=model_df, model_vars_dict=cfg['models'][k]['model_vars_dict'], 
																					model=cfg['models'][k]['model'], iid=cfg['models'][k]['iid'], fid=cfg['models'][k]['fid'], 
																					model_fxn=cfg['models'][k]['model_fxn'], family=cfg['models'][k]['family'], focus=cfg['models'][k]['focus'], 
																					dep_var=cfg['models'][k]['dep_var'], corstr=cfg['models'][k]['corstr']), 1)
						
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['models'][k]['model_fxn'] in ['bglm','gglm']:
						for x in cfg['models'][k]['focus']:
							marker_info_passed[x + '.effect'] = float('nan')
							marker_info_passed[x + '.stderr'] = float('nan')
							marker_info_passed[x + '.or'] = float('nan')
							marker_info_passed[x + '.z'] = float('nan')
							marker_info_passed[x + '.p'] = float('nan')
						marker_info_passed['n'] = 0
						marker_info_passed['status'] = 0
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcGLM(marker_info=row, model_df=model_df, model_vars_dict=cfg['models'][k]['model_vars_dict'], 
																								model=cfg['models'][k]['model'], iid=cfg['models'][k]['iid'], fid=cfg['models'][k]['fid'], 
																								family=cfg['models'][k]['family'], 
																								focus=cfg['models'][k]['focus'], dep_var=cfg['models'][k]['dep_var']), 1)
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['models'][k]['model_fxn'] == 'blme':
						for x in cfg['models'][k]['focus']:
							marker_info_passed[x + '.effect'] = float('nan')
							marker_info_passed[x + '.stderr'] = float('nan')
							marker_info_passed[x + '.or'] = float('nan')
							marker_info_passed[x + '.z'] = float('nan')
							marker_info_passed[x + '.p'] = float('nan')
							if cfg['models'][k]['lrt'] and x != '(Intercept)':
								marker_info_passed[x + '.anova.chisq'] = float('nan')
								marker_info_passed[x + '.anova.p'] = float('nan')
						marker_info_passed['n'] = 0
						marker_info_passed['status'] = 0
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcLMEBinomial(marker_info=row, model_df=model_df, model_vars_dict=cfg['models'][k]['model_vars_dict'], 
																					model=cfg['models'][k]['model'], iid=cfg['models'][k]['iid'], fid=cfg['models'][k]['fid'], 
																					model_fxn=cfg['models'][k]['model_fxn'], focus=cfg['models'][k]['focus'], 
																					dep_var=cfg['models'][k]['dep_var'], lmer_ctrl=cfg['models'][k]['lmer_ctrl'],
																					lrt=cfg['models'][k]['lrt']), 1)
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['models'][k]['model_fxn'] == 'glme':
						for x in cfg['models'][k]['focus']:
							marker_info_passed[x + '.effect'] = float('nan')
							marker_info_passed[x + '.stderr'] = float('nan')
							marker_info_passed[x + '.or'] = float('nan')
							marker_info_passed[x + '.z'] = float('nan')
							marker_info_passed[x + '.p'] = float('nan')
							marker_info_passed[x + '.satt.df'] = float('nan')
							marker_info_passed[x + '.satt.t'] = float('nan')
							marker_info_passed[x + '.satt.p'] = float('nan')
							marker_info_passed[x + '.kenrog.p'] = float('nan')
							if cfg['models'][k]['lrt'] and x != '(Intercept)':
								marker_info_passed[x + '.anova.chisq'] = float('nan')
								marker_info_passed[x + '.anova.p'] = float('nan')
						marker_info_passed['n'] = 0
						marker_info_passed['status'] = 0
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcLMEGaussian(marker_info=row, model_df=model_df, model_vars_dict=cfg['models'][k]['model_vars_dict'], 
																					model=cfg['models'][k]['model'], iid=cfg['models'][k]['iid'], fid=cfg['models'][k]['fid'], 
																					model_fxn=cfg['models'][k]['model_fxn'], focus=cfg['models'][k]['focus'], 
																					dep_var=cfg['models'][k]['dep_var'], lmer_ctrl=cfg['models'][k]['lmer_ctrl'],
																					reml=cfg['models'][k]['reml'],lrt=cfg['models'][k]['lrt']), 1)
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['models'][k]['model_fxn'] == 'cph':
						for x in cfg['models'][k]['focus']:
							marker_info_passed[x + '.effect'] = float('nan')
							marker_info_passed[x + '.or'] = float('nan')
							marker_info_passed[x + '.ci_lower'] = float('nan')
							marker_info_passed[x + '.ci_upper'] = float('nan')
							marker_info_passed[x + '.stderr'] = float('nan')
							marker_info_passed[x + '.robust_stderr'] = float('nan')
							marker_info_passed[x + '.z'] = float('nan')
							marker_info_passed[x + '.p'] = float('nan')
						marker_info_passed['n'] = 0
						marker_info_passed['status'] = 0
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcCoxPH(marker_info=row, model_df=model_df, model_vars_dict=cfg['models'][k]['model_vars_dict'], 
																						model=cfg['models'][k]['model'], iid=cfg['models'][k]['iid'], fid=cfg['models'][k]['fid'], 
																						model_fxn=cfg['models'][k]['model_fxn'], focus=cfg['models'][k]['focus'], 
																						dep_var=cfg['models'][k]['dep_var'],cph_ctrl=cfg['models'][k]['cph_ctrl']), 1)
						results = pd.merge(results,marker_info_filtered,how='outer')
					if 'id' in reglist.columns and not reglist['id'][r] is None:
						results['id'] = reglist['id'][r]
					if 'marker_unique' in results.columns.values:
						results.drop('marker_unique',axis=1,inplace=True)
					if len(cfg['model_order']) > 1:
						results.columns = [k + '.' + a if not a in ['chr','pos','a1','a2'] else a for a in results.columns]
					if cfg['models'][k]['written'] == False:
						cfg['models'][k]['results'] = results.copy()
						cfg['models'][k]['written'] = True
					else:
						cfg['models'][k]['results'] = cfg['models'][k]['results'].append(results).reset_index(drop=True)

				##### PREPARE DATA FOR GENE BASED ANALYSIS #####
				else:
					if cfg['models'][k]['model_fxn_type'] == 'gene':
						if i == 1:
							cfg['reg_model_df'] = model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1)
							cfg['reg_marker_info'] = marker_info[marker_info['filter'] == 0]
						else:
							cfg['reg_model_df'] = pd.merge(cfg['reg_model_df'],model_df.drop(marker_info['marker_unique'][marker_info['filter'] != 0],axis=1),how='outer',copy=False)
							cfg['reg_marker_info'] = cfg['reg_marker_info'].append(marker_info[marker_info['filter'] == 0],ignore_index=True)

				##### UPDATE LOOP STATUS #####
				cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + len(marker_info.index)))
				status_reg = reglist['region'][r] + ': ' + reglist['id'][r] if 'id' in reglist.columns and reglist['id'][r] != 'NA' else reglist['region'][r]
				status = '   processed ' + cur_markers + ' markers from region ' + str(r+1) + '/' + str(len(reglist.index)) + ' (' + status_reg + ')' if len(cfg['models'].keys()) == 1 else '   processed ' + cur_markers + ' markers from region ' + str(r+1) + '/' + str(len(reglist.index)) + ' (' + status_reg + ')' + " for cohort " + str(cfg['model_order'].index(k) + 1) + '/' + str(len(cfg['model_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE EFFECTIVE TESTS #####
			if cfg['models'][k]['model_fxn'] == 'neff':
				tot_tests = cfg['reg_marker_info'].shape[0]
				if tot_tests > 0:
					n_eff = StatsFxns.CalcEffTests(model_df=cfg['reg_model_df'][cfg['reg_marker_info']['marker_unique']])
					status = 0
				else:
					n_eff, tot_tests, status = (0, 0, 1)
				results = pd.DataFrame({'chr': [reglist['chr'][r]], 'start': [reglist['start'][r]], 'end': [reglist['end'][r]], 'id': [reglist['id'][r]], 'n_total': [tot_tests], 'n_eff': [n_eff], 'status': [status]})[['chr','start','end','id','n_total','n_eff','status']]
				if 'id' in reglist.columns:
					results['id'] = reglist['id'][r]
				if 'marker_unique' in results.columns.values:
					results.drop('marker_unique',axis=1,inplace=True)
				if len(cfg['model_order']) > 1:
					results.columns = [k + '.' + a if not a in ['chr','start','end','id'] else a for a in results.columns]
				if cfg['models'][k]['written'] == False:
					cfg['models'][k]['results'] = results.copy()
					cfg['models'][k]['written'] = True
				else:
					cfg['models'][k]['results'] = cfg['models'][k]['results'].append(results).reset_index(drop=True)
				status = '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(reglist.index)) if len(cfg['models'].keys()) == 1 else '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(reglist.index)) + " for model " + str(cfg['model_order'].index(k) + 1) + '/' + str(len(cfg['model_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE PREPSCORES AND INDIVIDUAL MODELS FOR GENE BASED TESTS #####
			elif cfg['models'][k]['model_fxn'] in seqmeta_tests:
				if not cfg['reg_marker_info'] is None and cfg['reg_marker_info'].shape[0] > 0:
					cfg['models'][k]['snp_info'] = pd.DataFrame({'Name': cfg['reg_marker_info']['marker_unique'], 'gene': reglist['id'][r]})
					ro.globalenv['snp_info'] = py2r.convert_to_r_dataframe(pd.DataFrame({'Name': cfg['reg_marker_info']['marker_unique'], 'gene': reglist['id'][r]}), strings_as_factors=False)
					ro.globalenv['z'] = py2r.convert_to_r_dataframe(cfg['reg_model_df'][list(cfg['reg_marker_info']['marker_unique'][cfg['reg_marker_info']['filter'] == 0])], strings_as_factors=False)
					ro.globalenv['pheno'] = py2r.convert_to_r_dataframe(cfg['reg_model_df'][list(set(cfg['models'][k]['model_vars_dict'].keys() + [cfg['models'][k]['iid'],cfg['models'][k]['fid']]))], strings_as_factors=False)

					##### PREPSCORES #####
					if cfg['models'][k]['model_fxn'] in ['famskat_o','famskat','famburden']:
						ro.globalenv['pedigree'] = py2r.convert_to_r_dataframe(cfg['models'][k]['ped_df'], strings_as_factors=False)
						ro.globalenv['kins'] = ro.r("makekinship(pedigree.rx2('FID'),pedigree.rx2('IID'),pedigree.rx2('PAT'),pedigree.rx2('MAT'))")
						cmd = 'prepScores(Z=z,formula=' + cfg['models'][k]['model'] + ',SNPInfo=snp_info,data=pheno,kins=kins,sparse=TRUE)'
						ro.globalenv['ps' + k] = ro.r(cmd)
					else:
						ro.globalenv['family'] = ro.r('gaussian()') if cfg['models'][k]['family'] is None or cfg['models'][k]['family'] == 'gaussian' else ro.r('binomial()')
						cmd = 'prepScores(Z=z,formula=' + cfg['models'][k]['model'] + ',SNPInfo=snp_info,data=pheno,family=family)'
						ro.globalenv['ps' + k] = ro.r(cmd)

					##### SKAT-O #####
					if cfg['models'][k]['model_fxn'] in ['fskato','gskato','bskato']:
						ro.globalenv['rho'] = ro.r('seq(0,1,' + str(cfg['models'][k]['rho']) + ')')
						results_pre = StatsFxns.SkatOMeta('skatOMeta(ps' + k + ',SNPInfo=snp_info,rho=rho,skat.wts=' + cfg['models'][k]['skat_wts'] + ',burden.wts=' + cfg['models'][k]['burden_wts'] + ',method="' + cfg['models'][k]['skat_method'] + '")')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],
											'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['models'][k]['model_fxn'] in ['fskat','gskat','bskat']:
						results_pre = StatsFxns.SkatMeta('skatMeta(ps' + k + ',SNPInfo=snp_info,wts=' + cfg['models'][k]['skat_wts'] + ',method="' + cfg['models'][k]['skat_method'] + '")')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],
											'p': [results_pre['p'][1]],'Qmeta': [results_pre['Qmeta'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]]})
						results = results[['chr','start','end','id','p','Qmeta','cmaf','nmiss','nsnps']]

					##### BURDEN #####
					elif cfg['models'][k]['model_fxn'] in ['fburden','gburden','bburden']:
						results_pre = StatsFxns.BurdenMeta('burdenMeta(ps' + k + ',SNPInfo=snp_info,wts=' + cfg['models'][k]['burden_wts'] + ')')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SNP_INFO DF #####
					cfg['models'][k]['snp_info'] = pd.DataFrame({'Name': [], 'gene': []})

					##### EMPTY SKAT-O DF #####
					if cfg['models'][k]['model_fxn'] in ['fskato','gskato','bskato']:
						results = StatsFxns.SkatOMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r])

					##### EMPTY SKAT DF #####
					elif cfg['models'][k]['model_fxn'] in ['fskat','gskat','bskat']:
						results = StatsFxns.SkatMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r])

					##### EMPTY BURDEN DF #####
					elif cfg['models'][k]['model_fxn'] in ['fburden','gburden','bburden']:
						results = StatsFxns.BurdenMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r])
					
				##### APPEND TO RESULTS DF #####
				if len(cfg['model_order']) > 1:
					results.columns = [k + '.' + a if not a in ['chr','start','end','id'] else a for a in results.columns]
				if cfg['models'][k]['written'] == False:
					cfg['models'][k]['results'] = results.copy()
					cfg['models'][k]['written'] = True
				else:
					cfg['models'][k]['results'] = cfg['models'][k]['results'].append(results).reset_index(drop=True)

		##### PERFORM GENE BASED TEST META ANALYSIS #####
		if cfg['models'][cfg['model_order'][0]]['model_fxn'] in seqmeta_tests and len(cfg['meta']) > 0:
			for meta in cfg['meta']:
				meta_tag = meta.split(':')[0]
				meta_incl = meta.split(':')[1].split('+')
				seqmeta_cmd = ''
				snp_info_meta = None
				meta_incl_string = ''
				for k in meta_incl:
					if 'ps' + k in ro.globalenv and cfg['models'][k]['results'].loc[r][k + '.p'] != 'NA':
						meta_incl_string = meta_incl_string + '+'
						if snp_info_meta is None:
							snp_info_meta = cfg['models'][k]['snp_info']
						else:
							snp_info_meta = snp_info_meta.merge(cfg['models'][k]['snp_info'], how='outer')
						if seqmeta_cmd == '':
							seqmeta_cmd = 'ps' + k
						else:
							seqmeta_cmd = seqmeta_cmd + ', ps' + k
					else:
						meta_incl_string = meta_incl_string + 'x'
				if snp_info_meta is not None and snp_info_meta.shape[0] > 0:
					ro.globalenv['snp_info_meta'] = py2r.convert_to_r_dataframe(snp_info_meta, strings_as_factors=False)
					##### SKAT-O #####
					if cfg['models'][cfg['model_order'][0]]['model_fxn'] in ['fskato','gskato','bskato']:
						results_pre = StatsFxns.SkatOMeta('skatOMeta(' + seqmeta_cmd + ',SNPInfo=snp_info_meta,rho=rho,skat.wts=' + cfg['models'][cfg['model_order'][0]]['skat_wts'] + ',burden.wts=' + cfg['models'][cfg['model_order'][0]]['burden_wts'] + ',method="' + cfg['models'][cfg['model_order'][0]]['skat_method'] + '")')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],'incl': [meta_incl_string],
											'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','id','incl','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['models'][cfg['model_order'][0]]['model_fxn'] in ['fskat','gskat','bskat']:
						results_pre = StatsFxns.SkatMeta('skatMeta(' + seqmeta_cmd + ', SNPInfo=snp_info_meta,wts=' + cfg['models'][cfg['model_order'][0]]['skat_wts'] + ',method="' + cfg['models'][cfg['model_order'][0]]['skat_method'] + '")')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],'meta_incl': [meta_incl_string],
													'p': [results_pre['p'][1]],'Qmeta': [results_pre['Qmeta'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
													'nsnps': [results_pre['nsnps'][1]]})
						results = results[['chr','start','end','id','p','Qmeta','cmaf','nmiss','nsnps']]

					##### BURDEN #####
					elif cfg['models'][cfg['model_order'][0]]['model_fxn'] in ['fburden','gburden','bburden']:
						results_pre = StatsFxns.BurdenMeta('burdenMeta(' + seqmeta_cmd + ', SNPInfo=snp_info_meta,wts=' + cfg['models'][cfg['model_order'][0]]['burden_wts'] + ')')
						results = pd.DataFrame({'chr': [reglist['chr'][r]],'start': [reglist['start'][r]],'end': [reglist['end'][r]],'id': [reglist['id'][r]],'meta_incl': [meta_incl_string],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SKAT-O DF #####
					if cfg['models'][k]['model_fxn'] in ['fskato','gskato','bskato']:
						results = StatsFxns.SkatOMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r],meta_incl_string)

					##### EMPTY SKAT DF #####
					elif cfg['models'][k]['model_fxn'] in ['fskat','gskat','bskat']:
						results = StatsFxns.SkatMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r],meta_incl_string)

					##### EMPTY BURDEN DF #####
					elif cfg['models'][k]['model_fxn'] in ['fburden','gburden','bburden']:
						results = StatsFxns.BurdenMetaEmpty(reglist['chr'][r],reglist['start'][r],reglist['end'][r],reglist['id'][r],meta_incl_string)

				##### APPEND TO META DF #####
				if len(cfg['model_order']) > 1:
					results.columns = [meta_tag + '.' + a if not a in ['chr','start','end','id'] else a for a in results.columns]
				if cfg['meta_written'][meta_tag] == False:
					cfg['meta_results'][meta_tag] = results.copy()
					cfg['meta_written'][meta_tag] = True
				else:
					cfg['meta_results'][meta_tag] = cfg['meta_results'][meta_tag].merge(results, how='outer', copy=False)

			##### REMOVE PREPSCORES OBJECTS FROM R GLOBAL ENVIRONMENT #####
			for k in cfg['model_order']:
				if 'ps' + k in ro.globalenv:
					ro.r['rm']('ps' + k)

	##### COMPILE ALL RESULTS #####
	header = ['chr','start','end','id'] if list(set([cfg['models'][k]['model_fxn_type'] for k in cfg['model_order']]))[0] == 'gene' else ['chr','pos','a1','a2']
	if len(cfg['meta']) > 0:
		for meta in cfg['meta']:
			meta_tag = meta.split(':')[0]
			if meta == cfg['meta'][0]:
				results_out = cfg['meta_results'][meta_tag]
			else:
				results_out = results_out.merge(cfg['meta_results'][meta_tag], how='outer', copy=False)
			header = header + [a for a in cfg['meta_results'][meta_tag].columns.values.tolist() if not a in header]

	for k in cfg['model_order']:
		if 'results' in cfg['models'][k]:
			if 'id' in cfg['models'][k]['results'].keys() and len(list(cfg['models'][k]['results']['id'].unique())) == 1 and list(cfg['models'][k]['results']['id'].unique())[0] == 'NA':
				cfg['models'][k]['results'].drop('id',axis=1,inplace=True)
			if k + '.id' in cfg['models'][k]['results'].keys() and len(list(cfg['models'][k]['results'][k + '.id'].unique())) == 1 and list(cfg['models'][k]['results'][k + '.id'].unique())[0] == 'NA':
				cfg['models'][k]['results'].drop(k + '.id',axis=1,inplace=True)
			if cfg['models'][k]['family'] is not None and cfg['models'][k]['family'] == 'gaussian':
				cfg['models'][k]['results'].drop([x for x in cfg['models'][k]['results'].columns if '.or' in x], axis=1,inplace=True)
			if (k == cfg['model_order'][0] and len(cfg['meta']) == 0) or not 'results_out' in locals():
				results_out = cfg['models'][k]['results']
			else:
				results_out = results_out.merge(cfg['models'][k]['results'], how='outer', copy=False)
			header = header + [a for a in cfg['models'][k]['results'].columns.values.tolist() if not a in header]

	##### SORT RESULTS AND CONVERT CHR, POS, AND START COLUMNS TO INT #####
	if list(set([cfg['models'][k]['model_fxn_type'] for k in cfg['model_order']]))[0] == 'marker':
		results_out[['chr','pos']] = results_out[['chr','pos']].astype(int)
		results_out.sort(columns=['chr','pos'],inplace=True)
	else:
		results_out[['chr','start','end']] = results_out[['chr','start','end']].astype(int)
		results_out.sort(columns=['chr','start'],inplace=True)
	sci_cols = [x for x in results_out.columns.values.tolist() if x.endswith(('.p','.pmin')) or x in ['.p','.pmin']]
	float_cols = [x for x in results_out.columns.values.tolist() if x.endswith(('.rho','.cmaf','.Qmeta','.beta','.se','.cmafTotal','.cmafUsed','.hwe','.hwe.unrel','.hwe.ctrl','.hwe.case','.hwe.unrel.ctrl','.hwe.unrel.case','.callrate','.freq','.freq.unrel','.freq.ctrl','.freq.case','.freq.unrel.ctrl','.freq.unrel.case','.rsq','.rsq.unrel','.rsq.ctrl','.rsq.case','.rsq.unrel.ctrl','.rsq.unrel.case','.effect','.stderr','.or','.z')) or x in ['rho','cmaf','Qmeta','beta','se','cmafTotal','cmafUsed','hwe','hwe.unrel','hwe.ctrl','hwe.case','hwe.unrel.ctrl','hwe.unrel.case','callrate','freq','freq.unrel','freq.ctrl','freq.case','freq.unrel.ctrl','freq.unrel.case','rsq','rsq.unrel','rsq.ctrl','rsq.case','rsq.unrel.ctrl','rsq.unrel.case','effect','stderr','or','z']]
	int_cols = [x for x in results_out.columns.values.tolist() if x.endswith(('.filter','.n','.status','.nmiss','.nsnps','.errflag','.nsnpsTotal','.nsnpsUsed')) or x in ['filter','n','status','nmiss','nsnps','errflag','nsnpsTotal','nsnpsUsed']]
	if len(sci_cols) > 0:
		for c in sci_cols:
			results_out[c] = results_out[c].map(lambda x: '%.4e' % (x) if not math.isnan(x) else x)
			results_out[c] = results_out[c].astype(object)
	if len(float_cols) > 0:
		for c in float_cols:
			results_out[c] = results_out[c].map(lambda x: '%.5g' % (x) if not math.isnan(x) else x)
			results_out[c] = results_out[c].astype(object)
	if len(int_cols) > 0:
		for c in int_cols:
			results_out[c] = results_out[c].map(lambda x: '%d' % (x) if not math.isnan(x) else x)
			results_out[c] = results_out[c].astype(object)

	##### FILL IN NA's, ORDER HEADER, AND WRITE TO FILE #####
	results_out.fillna('NA',inplace=True)
	results_out = results_out[header]
	bgzfile.write("\t".join(['#' + x if x == 'chr' else x for x in results_out.columns.values.tolist()]) + '\n')
	bgzfile.flush()
	results_out.to_csv(bgzfile, header=False, sep='\t', index=False)
	bgzfile.close()

	##### MAP OUTPUT FILES FOR TABIX #####
	if list(set([cfg['models'][k]['model_fxn_type'] for k in cfg['model_order']]))[0] == 'marker':
		cmd = ['tabix','-b','2','-e','2',cfg['out'] + '.gz']
	else:
		cmd = ['tabix','-b','2','-e','3',cfg['out'] + '.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessSystemFxns.Error:
		print SystemFxns.Error("file mapping failed")
	else:
		print "process complete"
