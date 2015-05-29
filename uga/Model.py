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

def Model(cfg):
	print "model options ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			print "   {0:>{1}}".format(str(k), len(max(cfg.keys(),key=len))) + ": " + str(cfg[k])
	cfg['data_info']['NA']['pheno_sep'] = MiscFxnsCy.GetDelimiterCy(cfg['data_info']['NA']['pheno_sep'])

	##### DEFINE MODEL TYPES #####
	marker_tests = ['gee_gaussian','gee_binomial','glm_gaussian','glm_binomial','lme_gaussian','lme_binomial','coxph','geeboss_gaussian','geeboss_binomial']
	seqmeta_tests = ['famskat_o','skat_o_gaussian','skat_o_binomial','famskat','skat_gaussian','skat_binomial','famburden','burden_gaussian','burden_binomial']
	efftests=['efftests']

	##### LOAD NECESSARY R PACKAGES #####
	methods = []
	for k in cfg['data_order']:
		methods.append(cfg['data_info'][k]['method'])
	for method in list(set(methods)):
		if method in ['gee_gaussian','gee_binomial','geeboss_gaussian','geeboss_binomial']:
			importr('geepack')
			if method in ['geeboss_gaussian','geeboss_binomial']:
				importr('boss')
		elif method in ['lme_gaussian','lme_binomial']:
			importr('lme4')
		elif method == 'coxph':
			importr('survival')
		elif method in ['famskat_o','skat_o_gaussian','skat_o_binomial','famskat','skat_gaussian','skat_binomial','famburden','burden_gaussian','burden_binomial']:
			importr('seqMeta')
			if method in ['famskat_o','famskat','famburden']:
				importr('kinship2')

	##### READ model VARIABLES FROM FILE #####
	vars_df_dict = {}
	model_vars_dict_dict = {}
	for k in cfg['data_order']:
		if len(cfg['data_info'].keys()) > 1:
			print "extracting model variables for model " + k + " and removing missing/invalid samples ..."
		else:
			print "extracting model variables and removing missing/invalid samples ..."
		cfg['data_info'][k]['vars_df'], cfg['data_info'][k]['model_vars_dict'] = MiscFxns.ExtractModelVars(pheno=cfg['data_info'][k]['pheno'], model=cfg['data_info'][k]['model'],
																										fid=cfg['data_info'][k]['fid'], iid=cfg['data_info'][k]['iid'],
																										fxn=cfg['data_info'][k]['fxn'], sex=cfg['data_info'][k]['sex'],
																										case=cfg['data_info'][k]['case'], ctrl=cfg['data_info'][k]['ctrl'],
																										pheno_sep=cfg['data_info'][k]['pheno_sep'])

	##### DETERMINE MODEL STATS TO BE EXTRACTED #####
	for k in cfg['data_order']:
		if 'focus' in cfg['data_info'][k].keys() and not cfg['data_info'][k]['focus'] is None:
			cfg['data_info'][k]['focus'] = cfg['data_info'][k]['focus'].split(',')
		else:
			cfg['data_info'][k]['focus'] = MiscFxns.GetFocus(method=cfg['data_info'][k]['method'],model=cfg['data_info'][k]['model'],vars_df=cfg['data_info'][k]['vars_df'],model_vars_dict=cfg['data_info'][k]['model_vars_dict'])

	##### REMOVE FACTOR() FROM MODEL STRING #####
	for k in cfg['data_order']:
		for v in cfg['data_info'][k]['model_vars_dict'].keys():
			cfg['data_info'][k]['model'] = cfg['data_info'][k]['model'].replace('factor(' + v + ')',v)

	##### LOAD ITERATORS AND SAMPLE LISTS #####
	for k in cfg['data_order']:
		if cfg['data_info'][k]['format'] == 'plink':
			print "reading Plink binary files for model " + k if len(cfg['data_info'].keys()) > 1 else "reading Plink binary files"
			cfg['data_info'][k]['plink_handle'],cfg['data_info'][k]['plink_locus_it'],cfg['data_info'][k]['plink_sample_it'],cfg['data_info'][k]['sample_ids'] = MiscFxns.LoadPlink(cfg['data_info'][k]['data'])
			cfg['data_info'][k]['plink_zip'] = zip(cfg['data_info'][k]['plink_locus_it'],cfg['data_info'][k]['plink_handle'])
		elif cfg['data_info'][k]['format'] == 'vcf':
			print "reading vcf file for model " + k if len(cfg['data_info'].keys()) > 1 else "reading vcf file"
			cfg['data_info'][k]['data_it'], cfg['data_info'][k]['sample_ids'] = MiscFxns.LoadVcf(cfg['data_info'][k]['data'])
		else:
			print "reading data and sample files for model " + k if len(cfg['data_info'].keys()) > 1 else "reading data and sample files"
			cfg['data_info'][k]['data_it'], cfg['data_info'][k]['sample_ids'] = MiscFxns.LoadDos(cfg['data_info'][k]['data'], cfg['data_info'][k]['samples'])
		if len(cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].isin(cfg['data_info'][k]['sample_ids'])]) == 0:
			print Error("phenotype file and data file contain no common samples")
			return
		cfg['data_info'][k]['marker_list_handle'] = None
		if cfg['data_info'][k]['marker_list'] is not None:
			cfg['data_info'][k]['marker_list_handle'] = tabix.open(cfg['data_info'][k]['marker_list'])

	##### REDUCE PHENO DATA TO GENOTYPE IDs #####
	for k in cfg['data_order']:
		cfg['data_info'][k]['vars_df'] = cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].isin(cfg['data_info'][k]['sample_ids'])]

	##### CREATE SUMMARY FOR DATA TO BE ANALYZED #####
	for k in cfg['data_order']:
		cfg['data_info'][k]['samples'] = len(cfg['data_info'][k]['vars_df'].index)
		cfg['data_info'][k]['vars_df_nodup'] = cfg['data_info'][k]['vars_df'].drop_duplicates(subset=[cfg['data_info'][k]['iid']])
		cfg['data_info'][k]['samples_unique'] = len(cfg['data_info'][k]['vars_df_nodup'].index)
		cfg['data_info'][k]['clusters'] = len(cfg['data_info'][k]['vars_df_nodup'].drop_duplicates(subset=[cfg['data_info'][k]['fid']]).index)
		cfg['data_info'][k]['families'] = len(cfg['data_info'][k]['vars_df_nodup'][cfg['data_info'][k]['vars_df_nodup'].duplicated(subset=[cfg['data_info'][k]['fid']])][cfg['data_info'][k]['fid']].unique())
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
			if cfg['data_info'][k]['pedigree'] is not None:
				print "loading pedigree for model " + k if len(cfg['data_info'].keys()) > 1 else "loading pedigree"
				cfg['data_info'][k]['ped_df'] = pd.read_table(cfg['data_info'][k]['pedigree'],sep='\t',dtype='str',usecols=['FID','IID','PAT','MAT'])
				cfg['data_info'][k]['ped_df'] = cfg['data_info'][k]['ped_df'][cfg['data_info'][k]['ped_df']['IID'].isin(list(cfg['data_info'][k]['vars_df'][cfg['data_info'][k]['iid']].values))]
				ro.globalenv['pedigree'] = py2r.convert_to_r_dataframe(cfg['data_info'][k]['ped_df'], strings_as_factors=False)
				cfg['data_info'][k]['kins'] = ro.r("makekinship(pedigree.rx2('FID'),pedigree.rx2('IID'),pedigree.rx2('PAT'),pedigree.rx2('MAT'))")

	##### GENERATE REGION LIST #####
	if not cfg['region_list'] is None:
		print "loading region list"
		region_list = FileFxns.LoadCoordinates(cfg['region_list'])
	elif not cfg['region'] is None:
		if len(cfg['region'].split(':')) > 1:
			region_list = pd.DataFrame({'chr': [re.split(':|-',cfg['region'])[0]],'start': [re.split(':|-',cfg['region'])[1]],'end': [re.split(':|-',cfg['region'])[2]],'region': [cfg['region']]})
		else:
			region_list = pd.DataFrame({'chr': [cfg['region']],'start': ['NA'],'end': ['NA'],'region': [cfg['region']]})
		region_list['reg_id'] = cfg['region_id'] if not cfg['region_id'] is None else 'NA'
	else:
		region_list = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)],'reg_id': ['NA' for i in range(26)]})

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
	for r in range(len(region_list.index)):
		mdb = MiscFxns.MarkerRefDb()
		reg = region_list['region'][r]
		chr = int(reg.split(':')[0])

		##### LOOP OVER DATASETS #####
		for k in cfg['data_order']:
			i = 0
			cfg['reg_model_df'] = None
			cfg['reg_marker_info'] = None
			marker_list = None
			if cfg['data_info'][k]['marker_list_handle'] is not None:
				marker_list = cfg['data_info'][k]['marker_list_handle'].querys(reg)
				marker_list = pd.DataFrame([x for x in marker_list if int(x[1]) >= int(region_list['start'][r]) and int(x[1]) <= int(region_list['end'][r])])
				marker_list.columns = ['chr','pos','marker','a1','a2']
				marker_list['chr'] = marker_list['chr'].astype(np.int64)
				marker_list['pos'] = marker_list['pos'].astype(np.int64)
			if cfg['data_info'][k]['format'] == 'plink':
				if reg == str(chr):
					try:
						records = (x for x in cfg['data_info'][k]['plink_zip'] if x[0].chromosome == chr and not ',' in x[0].allele2)
					except:
						break
				else:
					try:
						records = (x for x in cfg['data_info'][k]['plink_zip'] if x[0].chromosome == chr and x[0].bp_position >= int(region_list['start'][r]) and x[0].bp_position <= int(region_list['end'][r]) and not ',' in x[0].allele2)
					except:
						break
				if marker_list is not None:
					records = (x for x in records if marker_list[(marker_list['chr'] == x[0].chromosome) & (marker_list['pos'] == x[0].bp_position) & (marker_list['a1'] == x[0].allele1) & (marker_list['a2'] == x[0].allele2)].shape[0] > 0)
			else:
				pos_ind = 1 if cfg['data_info'][k]['format'] in ['dos2','vcf'] else 2
				try:
					records = cfg['data_info'][k]['data_it'].querys(reg)
					if reg != str(chr):
						records = (x for x in records if int(x[pos_ind]) >= int(region_list['start'][r]) and int(x[pos_ind]) <= int(region_list['end'][r]) and not ',' in str(x[4]))
				except:
					break
				if marker_list is not None:
					records = (x for x in records if marker_list[(marker_list['chr'] == int(x[0])) & (marker_list['pos'] == int(x[pos_ind])) & (marker_list['a1'] == x[3]) & (marker_list['a2'] == x[4])].shape[0] > 0)

			##### CONTINUE SLICING UNTIL NO RECORDS LEFT #####
			chunk_db = MiscFxns.ChunkRefDb()
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
					marker_info = pd.DataFrame(marker_info)
					marker_data = pd.DataFrame(marker_data)
					marker_data = marker_data.convert_objects(convert_numeric=True)
					marker_data[cfg['data_info'][k]['iid']] = marker_data.index
					chunkdf = marker_info.join(marker_data.transpose())
				else:
					if cfg['data_info'][k]['format'] == 'vcf':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf.apply(MiscFxnsCy.GetRowCallsCy,axis=1,raw=True)
						chunkdf.drop(chunkdf.columns[5:9],axis=1,inplace=True)
					elif cfg['data_info'][k]['format'] == 'oxford':
						chunk = [MiscFxns.ConvertDosage(row) for row in chunk]
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['data_info'][k]['format'] == 'dos1':
						chunkdf = pd.DataFrame(chunk)
						chunkdf = chunkdf[[0,2,1] + range(3,len(chunkdf.columns))]
					elif cfg['data_info'][k]['format'] == 'dos2':
						chunkdf = pd.DataFrame(chunk)
					chunkdf.columns = ['chr','pos','marker','a1','a2'] + cfg['data_info'][k]['sample_ids']
				if len(cfg['data_order']) > 1:
					chunkdf = chunkdf.drop_duplicates(subset=['chr','pos','a1','a2'])
					chunkdf.index=chunkdf['chr'].astype(str) + '><' + chunkdf['pos'].astype(str) + '><'  + chunkdf['a1'].astype(str) + '><'  + chunkdf['a2'].astype(str)
					chunkdf = chunkdf[~chunkdf.index.isin(chunk_db.ListKeys())]
					chunkdf.apply(chunk_db.Update,1)
					chunkdf = chunkdf.apply(mdb.Update,1)
				
				##### EXTRACT MARKER INFO AND MARKER DATA #####
				marker_info = chunkdf.ix[:,:5]
				marker_info.replace('.','NA',inplace=True)
				marker_info['marker_unique'] = 'chr' + marker_info['chr'].astype(str) + 'bp' + marker_info['pos'].astype(str) + '_'  + marker_info['marker'].astype(str).str.replace('-','_').str.replace(':','.') + '_'  + marker_info['a1'].astype(str) + '_'  + marker_info['a2'].astype(str)
				marker_info.index = marker_info['marker_unique']
				marker_data = chunkdf.ix[:,5:].transpose()
				marker_data = marker_data.convert_objects(convert_numeric=True)
				marker_data.columns = marker_info['marker_unique']
				marker_data[cfg['data_info'][k]['iid']] = marker_data.index

				##### MERGE PHENO DATAFRAME AND MARKER DATA #####
				model_df = pd.merge(cfg['data_info'][k]['vars_df'], marker_data, on = [cfg['data_info'][k]['iid']], how='left').sort([cfg['data_info'][k]['fid']])

				##### EXTRACT UNIQUE SAMPLES AND CALCULATE DESCRIPTIVE STATS #####
				model_df_nodup=model_df.drop_duplicates(subset=[cfg['data_info'][k]['iid']]).reset_index(drop=True)
				marker_info['callrate']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcCallrateCy(col.astype(np.float)), raw=True)
				model_df_nodup_unrel=model_df_nodup.drop_duplicates(subset=[cfg['data_info'][k]['fid']]).reset_index(drop=True)
				if not cfg['data_info'][k]['sex'] is None and not cfg['data_info'][k]['male'] is None and not cfg['data_info'][k]['female']:
					male_idx = model_df_nodup[model_df_nodup[cfg['data_info'][k]['sex']].isin([cfg['data_info'][k]['male']])].index.values
					female_idx = model_df_nodup[model_df_nodup[cfg['data_info'][k]['sex']].isin([cfg['data_info'][k]['female']])].index.values
				else:
					male_idx = None
					female_idx = None
				if chr == 23:
					marker_info['freq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreq23Cy(marker_male=col[male_idx].astype(np.float), marker_female=col[female_idx].astype(np.float)), raw=True)
					marker_info['freq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreq23Cy(marker_male=col[male_idx].astype(np.float), marker_female=col[female_idx].astype(np.float)), raw=True)
				else:
					marker_info['freq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
					marker_info['freq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['freq.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
					marker_info['freq.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
					marker_info['freq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
					marker_info['freq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcFreqCy(col.astype(np.float)), raw=True)
				marker_info['rsq']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
				marker_info['rsq.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['rsq.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
					marker_info['rsq.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
					marker_info['rsq.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
					marker_info['rsq.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcRsqCy(col.astype(np.float)), raw=True)
				marker_info['hwe']=model_df_nodup[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
				marker_info['hwe.unrel']=model_df_nodup_unrel[marker_info['marker_unique']].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
				if cfg['data_info'][k]['fxn'] == 'binomial':
					marker_info['hwe.ctrl']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
					marker_info['hwe.case']=model_df_nodup[model_df_nodup[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
					marker_info['hwe.unrel.ctrl']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '0'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
					marker_info['hwe.unrel.case']=model_df_nodup_unrel[model_df_nodup_unrel[cfg['data_info'][k]['dep_var'][0]] == '1'][list(marker_info['marker_unique'])].apply(lambda col: MiscFxnsCy.CalcHWECy(col.astype(np.float)), raw=True)
				if cfg['data_info'][k]['method'] in seqmeta_tests:
					marker_info['filter']=marker_info.apply(lambda row: MiscFxnsCy.GenerateFilterCodeCy(row_callrate=row['callrate'], no_mono=False, 
																											row_freq=row['freq'], row_freq_unrel=row['freq.unrel'], 
																											row_rsq=row['rsq'], row_rsq_unrel=row['rsq.unrel'], 
																											row_hwe=row['hwe'], row_hwe_unrel=row['hwe.unrel'], 
																											miss_thresh=cfg['miss'], freq_thresh=cfg['freq'], 
																											max_freq_thresh=cfg['max_freq'], rsq_thresh=cfg['rsq'], 
																											hwe_thresh=cfg['hwe']), axis=1, raw=True)
				else:
					marker_info['filter']=marker_info.apply(lambda row: MiscFxnsCy.GenerateFilterCodeCy(row_callrate=row['callrate'], 
																											row_freq=row['freq'], row_freq_unrel=row['freq.unrel'], 
																											row_rsq=row['rsq'], row_rsq_unrel=row['rsq.unrel'], 
																											row_hwe=row['hwe'], row_hwe_unrel=row['hwe.unrel'], 
																											miss_thresh=cfg['miss'], freq_thresh=cfg['freq'], 
																											max_freq_thresh=cfg['max_freq'], rsq_thresh=cfg['rsq'], 
																											hwe_thresh=cfg['hwe']), axis=1, raw=True)
				marker_info['samples'] = str(cfg['data_info'][k]['samples']) + '/' + str(cfg['data_info'][k]['samples_unique']) + '/' + str(cfg['data_info'][k]['clusters']) + '/' + str(cfg['data_info'][k]['cases']) + '/' + str(cfg['data_info'][k]['ctrls']) + '/' + str(cfg['data_info'][k]['nmale']) + '/' + str(cfg['data_info'][k]['nfemale'])

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
				if cfg['data_info'][k]['method'].split('_')[0] in ['gee','geeboss','glm','lme','coxph']:
					marker_info_passed = marker_info[marker_info['filter'] == 0]
					marker_info_filtered = marker_info[marker_info['filter'] > 0]
					marker_info_filtered['n'] = '%d' % (0)
					marker_info_filtered['status'] = '%d' % (-1)
					model_df = model_df[[c for c in model_df.columns if not c in marker_info_filtered['marker_unique']]]
					model_df[cfg['data_info'][k]['fid']] = pd.Categorical.from_array(model_df[cfg['data_info'][k]['fid']]).codes.astype(np.int64)
					model_df.sort([cfg['data_info'][k]['fid']],inplace = True)
					ro.globalenv['model_df'] = py2r.convert_to_r_dataframe(model_df, strings_as_factors=False)
					ro.r('model_df[,names(model_df) == "' + cfg['data_info'][k]['fid'] + '"]<-as.factor(model_df[,names(model_df) == "' + cfg['data_info'][k]['fid'] + '"])')
					for x in cfg['data_info'][k]['focus']:
						marker_info_passed[x + '.effect'] = float('nan')
						marker_info_passed[x + '.stderr'] = float('nan')
						marker_info_passed[x + '.or'] = float('nan')
						marker_info_passed[x + '.z'] = float('nan')
						marker_info_passed[x + '.p'] = float('nan')
					marker_info_passed['n'] = '%d' % (0)
					marker_info_passed['status'] = '%d' % (0)
					if cfg['data_info'][k]['method'].split('_')[0] == 'gee':
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcGEE(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], 
																					model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], 
																					method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], 
																					dep_var=cfg['data_info'][k]['dep_var'], corstr=cfg['data_info'][k]['corstr']), 1)
						
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['data_info'][k]['method'].split('_')[0] == 'geeboss':
						model_df[cfg['data_info'][k]['fid']] = pd.Categorical.from_array(model_df[cfg['data_info'][k]['fid']]).codes.astype(np.int64)
						model_df.sort([cfg['data_info'][k]['fid']],inplace = True)
						model_df['id'] = model_df[cfg['data_info'][k]['fid']] if cfg['data_info'][k]['fid'] != 'id' else model_df['id']
						results = marker_info.apply(lambda row: StatsFxns.CalcGEEBoss(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], 
																						model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], 
																						method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], 
																						dep_var=cfg['data_info'][k]['dep_var'], corstr=cfg['data_info'][k]['corstr'], 
																						thresh=cfg['data_info'][k]['geeboss_thresh'], boss=boss), 1)
					elif cfg['data_info'][k]['method'].split('_')[0] == 'glm':
						results = marker_info_passed.apply(lambda row: StatsFxns.CalcGLM(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], 
																								model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], 
																								fxn=cfg['data_info'][k]['fxn'], 
																								focus=cfg['data_info'][k]['focus'], dep_var=cfg['data_info'][k]['dep_var']), 1)
						results = pd.merge(results,marker_info_filtered,how='outer')
					elif cfg['data_info'][k]['method'].split('_')[0] == 'lme':
						results = marker_info.apply(lambda row: StatsFxns.CalcLME(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], 
																					model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], 
																					method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], 
																					dep_var=cfg['data_info'][k]['dep_var'], lme4=lme4), 1)
					elif cfg['data_info'][k]['method'] == 'coxph':
						results = marker_info.apply(lambda row: StatsFxns.CalcCoxPH(marker_info=row, model_df=model_df, model_vars_dict=cfg['data_info'][k]['model_vars_dict'], 
																						model=cfg['data_info'][k]['model'], iid=cfg['data_info'][k]['iid'], fid=cfg['data_info'][k]['fid'], 
																						method=cfg['data_info'][k]['method'], fxn=cfg['data_info'][k]['fxn'], focus=cfg['data_info'][k]['focus'], 
																						dep_var=cfg['data_info'][k]['dep_var'], survival=survival), 1)
					if 'reg_id' in region_list.columns and not region_list['reg_id'][r] is None:
						results['reg_id'] = region_list['reg_id'][r]
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
				cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + len(marker_info.index)))
				status_reg = region_list['region'][r] + ': ' + region_list['reg_id'][r] if 'reg_id' in region_list.columns and region_list['reg_id'][r] != 'NA' else region_list['region'][r]
				status = '   processed ' + cur_markers + ' markers from region ' + str(r+1) + '/' + str(len(region_list.index)) + ' (' + status_reg + ')' if len(cfg['data_info'].keys()) == 1 else '   processed ' + cur_markers + ' markers from region ' + str(r+1) + '/' + str(len(region_list.index)) + ' (' + status_reg + ')' + " for cohort " + str(cfg['data_order'].index(k) + 1) + '/' + str(len(cfg['data_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE EFFECTIVE TESTS #####
			if cfg['data_info'][k]['method'] == 'efftests':
				tot_tests = cfg['reg_marker_info'].shape[0]
				if tot_tests > 0:
					n_eff = StatsFxns.CalcEffTests(model_df=cfg['reg_model_df'][cfg['reg_marker_info']['marker_unique']])
					status = 0
				else:
					n_eff, tot_tests, status = (0, 0, 1)
				results = pd.DataFrame({'chr': [region_list['chr'][r]], 'start': [region_list['start'][r]], 'end': [region_list['end'][r]], 'reg_id': [region_list['reg_id'][r]], 'n_total': [tot_tests], 'n_eff': [n_eff], 'status': [status]})[['chr','start','end','reg_id','n_total','n_eff','status']]
				if 'reg_id' in region_list.columns:
					results['reg_id'] = region_list['reg_id'][r]
				if 'marker_unique' in results.columns.values:
					results.drop('marker_unique',axis=1,inplace=True)
				if cfg['data_info'][k]['written'] == False:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = results.copy()
					cfg['data_info'][k]['written'] = True
				else:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = cfg['data_info'][k]['results'].append(results).reset_index(drop=True)
				status = '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(region_list.index)) if len(cfg['data_info'].keys()) == 1 else '   processed effective tests calculation for region ' + str(r+1) + '/' + str(len(region_list.index)) + " for model " + str(cfg['data_order'].index(k) + 1) + '/' + str(len(cfg['data_order'])) + ' (' + k + ')'
				print status

			##### CALCULATE PREPSCORES AND INDIVIDUAL MODELS FOR GENE BASED TESTS #####
			elif cfg['data_info'][k]['method'] in seqmeta_tests:
				if not cfg['reg_marker_info'] is None and cfg['reg_marker_info'].shape[0] > 0:
					cfg['data_info'][k]['snp_info'] = pd.DataFrame({'Name': cfg['reg_marker_info']['marker_unique'], 'gene': region_list['reg_id'][r]})
					z = cfg['reg_model_df'][list(cfg['reg_marker_info']['marker_unique'][cfg['reg_marker_info']['filter'] == 0])]
					pheno = cfg['reg_model_df'][list(set(cfg['data_info'][k]['model_vars_dict'].keys() + [cfg['data_info'][k]['iid'],cfg['data_info'][k]['fid']]))]

					##### PREPSCORES #####
					if cfg['data_info'][k]['method'] in ['famskat_o','famskat','famburden']:
						ro.globalenv['ps' + k] = StatsFxns.PrepScoresFam(snp_info=cfg['data_info'][k]['snp_info'], z=z, model=cfg['data_info'][k]['model'], pheno=pheno, kinship=cfg['data_info'][k]['kins'], seqmeta=seqmeta)
					else:
						ro.globalenv['ps' + k] = StatsFxns.PrepScores(snp_info=cfg['data_info'][k]['snp_info'], z=z, model=cfg['data_info'][k]['model'], pheno=pheno, family=cfg['data_info'][k]['fxn'], seqmeta=seqmeta)

					##### SKAT-O #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results_pre = StatsFxns.SkatOMeta('skatOMeta(ps' + k + ', SNPInfo=rsnp_info, rho=r_rho)', cfg['data_info'][k]['snp_info'], cfg['data_info'][k]['skat_o_rho'], seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results_pre = StatsFxns.SkatMeta('skatMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'], seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'Qmeta': [results_pre['Qmeta'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
											'nsnps': [results_pre['nsnps'][1]]})
						results = results[['chr','start','end','reg_id','p','Qmeta','cmaf','nmiss','nsnps']]

					##### BURDEN #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results_pre = StatsFxns.BurdenMeta('burdenMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'], seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','reg_id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SNP_INFO DF #####
					cfg['data_info'][k]['snp_info'] = pd.DataFrame({'Name': [], 'gene': []})

					##### EMPTY SKAT-O DF #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results = StatsFxns.SkatOMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r])

					##### EMPTY SKAT DF #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results = StatsFxns.SkatMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r])

					##### EMPTY BURDEN DF #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results = StatsFxns.BurdenMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r])
					
				##### APPEND TO RESULTS DF #####
				if cfg['data_info'][k]['written'] == False:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = results.copy()
					cfg['data_info'][k]['written'] = True
				else:
					results.columns = [k + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['data_info'][k]['results'] = cfg['data_info'][k]['results'].append(results).reset_index(drop=True)

		##### PERFORM GENE BASED TEST META ANALYSIS #####
		if cfg['data_info'][cfg['data_order'][0]]['method'] in seqmeta_tests and len(cfg['meta']) > 0:
			for meta in cfg['meta']:
				meta_tag = meta.split(':')[0]
				meta_incl = meta.split(':')[1].split('+')
				seqmeta_cmd = ''
				snp_info_meta = None
				meta_incl_string = ''
				for k in meta_incl:
					if 'ps' + k in ro.globalenv and cfg['data_info'][k]['results'].loc[r][k + '.p'] != 'NA':
						meta_incl_string = meta_incl_string + '+'
						if snp_info_meta is None:
							snp_info_meta = cfg['data_info'][k]['snp_info']
						else:
							snp_info_meta = snp_info_meta.merge(cfg['data_info'][k]['snp_info'], how='outer')
						if seqmeta_cmd == '':
							seqmeta_cmd = 'ps' + k
						else:
							seqmeta_cmd = seqmeta_cmd + ', ps' + k
					else:
						meta_incl_string = meta_incl_string + 'x'

				if snp_info_meta is not None and snp_info_meta.shape[0] > 1:
					##### SKAT-O #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results_pre = StatsFxns.SkatOMeta('skatOMeta(' + seqmeta_cmd + ', SNPInfo=rsnp_info, rho=r_rho)', snp_info_meta, cfg['data_info'][k]['skat_o_rho'], seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],'incl': [meta_incl_string],
													'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
													'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','incl','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### SKAT #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results_pre = StatsFxns.SkatMeta('skatMeta(' + seqmeta_cmd + ', SNPInfo=rsnp_info)', snp_info_meta, seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],'meta_incl': [meta_incl_string],
													'p': [results_pre['p'][1]],'pmin': [results_pre['pmin'][1]],'rho': [results_pre['rho'][1]],'cmaf': [results_pre['cmaf'][1]],'nmiss': [results_pre['nmiss'][1]],
													'nsnps': [results_pre['nsnps'][1]],'errflag': [results_pre['errflag'][1]]})
						results = results[['chr','start','end','reg_id','incl','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

					##### BURDEN #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results_pre = StatsFxns.BurdenMeta('burdenMeta(ps' + k + ', SNPInfo=rsnp_info)', cfg['data_info'][k]['snp_info'], seqmeta=seqmeta)
						results = pd.DataFrame({'chr': [region_list['chr'][r]],'start': [region_list['start'][r]],'end': [region_list['end'][r]],'reg_id': [region_list['reg_id'][r]],'meta_incl': [meta_incl_string],
											'p': [results_pre['p'][1]],'beta': [results_pre['beta'][1]],'se': [results_pre['se'][1]],'cmafTotal': [results_pre['cmafTotal'][1]],
											'cmafUsed': [results_pre['cmafUsed'][1]],'nsnpsTotal': [results_pre['nsnpsTotal'][1]],'nsnpsUsed': [results_pre['nsnpsUsed'][1]],
											'nmiss': [results_pre['nmiss'][1]]})
						results = results[['chr','start','end','reg_id','incl','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
				else:

					##### EMPTY SKAT-O DF #####
					if cfg['data_info'][k]['method'] in ['famskat_o','skat_o_gaussian','skat_o_binomial']:
						results = StatsFxns.SkatOMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r],meta_incl_string)

					##### EMPTY SKAT DF #####
					elif cfg['data_info'][k]['method'] in ['famskat','skat_gaussian','skat_binomial']:
						results = StatsFxns.SkatMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r],meta_incl_string)

					##### EMPTY BURDEN DF #####
					elif cfg['data_info'][k]['method'] in ['famburden','burden_gaussian','burden_binomial']:
						results = StatsFxns.BurdenMetaEmpty(region_list['chr'][r],region_list['start'][r],region_list['end'][r],region_list['reg_id'][r],meta_incl_string)

				##### APPEND TO META DF #####
				if cfg['meta_written'][meta_tag] == False:
					results.columns = [meta_tag + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['meta_results'][meta_tag] = results.copy()
					cfg['meta_written'][meta_tag] = True
				else:
					results.columns = [meta_tag + '.' + a if not a in ['chr','start','end','reg_id'] and k != 'NA' else a for a in results.columns]
					cfg['meta_results'][meta_tag] = cfg['meta_results'][meta_tag].merge(results, how='outer', copy=False)

			##### REMOVE PREPSCORES OBJECTS FROM R GLOBAL ENVIRONMENT #####
			for k in cfg['data_order']:
				if 'ps' + k in ro.globalenv:
					ro.r['rm']('ps' + k)

	##### COMPILE ALL RESULTS #####
	header = ['chr','start','end','reg_id'] if list(set([cfg['data_info'][k]['method_type'] for k in cfg['data_order']]))[0] == 'gene' else ['chr','pos','a1','a2']
	if len(cfg['meta']) > 0:
		for meta in cfg['meta']:
			meta_tag = meta.split(':')[0]
			if meta == cfg['meta'][0]:
				results = cfg['meta_results'][meta_tag]
			else:
				results = results.merge(cfg['meta_results'][meta_tag], how='outer', copy=False)
			header = header + [a for a in cfg['meta_results'][meta_tag].columns.values.tolist() if not a in header]

	for k in cfg['data_order']:
		if 'reg_id' in cfg['data_info'][k]['results'].keys() and len(list(cfg['data_info'][k]['results']['reg_id'].unique())) == 1 and list(cfg['data_info'][k]['results']['reg_id'].unique())[0] == 'NA':
			cfg['data_info'][k]['results'].drop('reg_id',axis=1,inplace=True)
		if k + '.reg_id' in cfg['data_info'][k]['results'].keys() and len(list(cfg['data_info'][k]['results'][k + '.reg_id'].unique())) == 1 and list(cfg['data_info'][k]['results'][k + '.reg_id'].unique())[0] == 'NA':
			cfg['data_info'][k]['results'].drop(k + '.reg_id',axis=1,inplace=True)
		if cfg['data_info'][k]['fxn'] is not None and cfg['data_info'][k]['fxn'] == 'gaussian':
			cfg['data_info'][k]['results'].drop([x for x in cfg['data_info'][k]['results'].columns if '.or' in x], axis=1,inplace=True)
		if k == cfg['data_order'][0] and len(cfg['meta']) == 0:
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
	results[[x for x in results.columns if x.endswith(('.p','hwe','hwe.unrel','hwe.ctrl','hwe.case','hwe.unrel.ctrl','hwe.unrel.case','callrate','freq','freq.unrel','freq.ctrl','freq.case','freq.unrel.ctrl','freq.unrel.case','rsq','rsq.unrel','rsq.ctrl','rsq.case','rsq.unrel.ctrl','rsq.unrel.case','effect','stderr','or','z'))]] = results[[x for x in results.columns if x.endswith(('.p','hwe','hwe.unrel','hwe.ctrl','hwe.case','hwe.unrel.ctrl','hwe.unrel.case','callrate','freq','freq.unrel','freq.ctrl','freq.case','freq.unrel.ctrl','freq.unrel.case','rsq','rsq.unrel','rsq.ctrl','rsq.case','rsq.unrel.ctrl','rsq.unrel.case','effect','stderr','or','z'))]].astype(float)
	results[[x for x in results.columns if x.endswith(('filter','n','status'))]] = results[[x for x in results.columns if x.endswith(('filter','n','status'))]].astype(int)
	for c in [x for x in results.columns if x.endswith(('.p','hwe','hwe.unrel','hwe.ctrl','hwe.case','hwe.unrel.ctrl','hwe.unrel.case'))]:
		results[c] = results[c].map(lambda x: '%.4e' % (x) if not math.isnan(x) else x)
		results[c] = results[c].astype(object)
	for c in [x for x in results.columns if x.endswith(('callrate','freq','freq.unrel','freq.ctrl','freq.case','freq.unrel.ctrl','freq.unrel.case','rsq','rsq.unrel','rsq.ctrl','rsq.case','rsq.unrel.ctrl','rsq.unrel.case','effect','stderr','or','z'))]:
		results[c] = results[c].map(lambda x: '%.5g' % (x) if not math.isnan(x) else x)
		results[c] = results[c].astype(object)

	##### FILL IN NA's, ORDER HEADER, AND WRITE TO FILE #####
	results.fillna('NA',inplace=True)
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
