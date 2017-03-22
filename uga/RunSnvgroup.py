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

import pandas as pd
import numpy as np
import Model
import Parse
import Variant
import pysam
import Fxns
from Bio import bgzf
import Process
import multiprocessing as mp
import sys
import os
import resource
import logging
import pickle
import glob

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunSnvgroup")

def process_regions(regions_df, cfg, cpu, log):
	regions_df = regions_df[regions_df['cpu'] == cpu].reset_index(drop=True)

	if log:
		try:
			log_file = open(cfg['out'] + '.cpu' + str(cpu) + '.log','w')
		except:
			print Process.Error("unable to initialize log file " + cfg['out'] + '.cpu' + str(cpu) + '.log').out
			return 1
		stdout_orig = sys.stdout
		sys.stdout = log_file

	models_obj = {}
	variants_found = False
	model_written = {}
	meta_written = {}
	results_final_models = {}
	results_final_models_headers = {}
	results_final_meta = {}
	meta_objs = {}
	variants_files = {}
	variant_ref = Variant.Ref()
	model_loaded = {}
	for n in cfg['model_order']:
		model_written[n] = False
		results_final_models[n] = pd.DataFrame({})
		variants_files[n] = glob.glob(cfg['models'][n]['file'].replace('[CHR]','*'))
		model_loaded[n] = False
	for meta in cfg['meta_order']:
		meta_written[meta] = False
		results_final_meta[meta] = pd.DataFrame({})
		meta_objs[meta] = getattr(Model,cfg['models'][cfg['meta'][meta].split('+')[0]]['fxn'].capitalize() + 'Meta')(tag = meta, meta = cfg['meta'][meta])
	last_chr = None
	for k in xrange(len(regions_df.index)):
		meta_incl = []
		region_written = False
		results_region = pd.DataFrame({})
		print ''
		print 'loading region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['group_id'][k] + ": " + regions_df['region'][k] + ') ...'
		for n in cfg['model_order']:
			if not model_loaded[n] or (last_chr != regions_df['chr'][k] and len(variants_files[n]) > 1):
				if not model_loaded[n]:
					print "\nloading model for " + n if n != '___no_tag___' else "\nloading model"
				else:
					print "\nupdating model for " + n if n != '___no_tag___' else "\nupdating model"
				try:
					models_obj[n] = getattr(Model,cfg['models'][n]['fxn'].capitalize())(fxn=cfg['models'][n]['fxn'], 
																						snvgroup_map=cfg['snvgroup_map'], 
																						dep_var=cfg['models'][n]['dep_var'], 
																						covars=cfg['models'][n]['covars'], 
																						format=cfg['models'][n]['format'], 
																						skat_wts=cfg['models'][n]['skat_wts'], 
																						burden_wts=cfg['models'][n]['burden_wts'], 
																						skat_method=cfg['models'][n]['skat_method'], 
																						cmac=cfg['models'][n]['cmac'], 
																						mafrange=cfg['models'][n]['mafrange'], 
																						timeout=cfg['timeout'], 
																						all_founders=cfg['models'][n]['all_founders'], 
																						case_code=cfg['models'][n]['case_code'], 
																						ctrl_code=cfg['models'][n]['ctrl_code'], 
																						pheno=cfg['models'][n]['pheno'], 
																						variants_file=cfg['models'][n]['file'].replace('[CHR]',str(regions_df['chr'][k])), # variants_file=cfg['models'][n]['file']
																						samples_file=cfg['models'][n]['sample'], 
																						drop_file=cfg['models'][n]['drop'], 
																						keep_file=cfg['models'][n]['keep'], 
																						type=cfg['models'][n]['fxn'], 
																						fid=cfg['models'][n]['fid'], 
																						iid=cfg['models'][n]['iid'], 
																						matid=cfg['models'][n]['matid'], 
																						patid=cfg['models'][n]['patid'], 
																						sex=cfg['models'][n]['sex'], 
																						male=cfg['models'][n]['male'], 
																						female=cfg['models'][n]['female'], 
																						sep=cfg['models'][n]['sep'])
				except Process.Error as err:
					print err.out
					return 1
				model_loaded[n] = True

			try:
				models_obj[n].get_region(regions_df['region'][k], group_id=regions_df['group_id'][k])
			except:
				pass

			try:
				models_obj[n].get_snvgroup(cfg['buffer'], regions_df['group_id'][k])
			except:
				if not variants_found:
					print '   (' + n + ') processed 0 variants'
				pass
			variants_found = True

			if models_obj[n].variants.duplicated is not None:
				print '   WARNING! The following duplicated variant identifiers were generated'
				print '\n'.join(['      ' + d for d in models_obj[n].variants.duplicated])

			if len(cfg['meta_order']) > 0:
				if n == cfg['model_order'][0]:
					variant_ref.load(models_obj[n].variants.info)
				else:
					variant_ref.update(models_obj[n].variants.info)
					models_obj[n].variants.align(variant_ref)

			try:
				models_obj[n].filter(miss_thresh=cfg['models'][n]['miss'], maf_thresh=cfg['models'][n]['maf'], maxmaf_thresh=cfg['models'][n]['maxmaf'], 
								mac_thresh=cfg['models'][n]['mac'], rsq_thresh=cfg['models'][n]['rsq'], hwe_thresh=cfg['models'][n]['hwe'], 
								hwe_maf_thresh=cfg['models'][n]['hwe_maf'], allow_mono=cfg['models'][n]['allow_mono'])
			except:
				pass

			try:
				logger.debug("calc_model")
				models_obj[n].calc_model()
			except Process.Error as err:
				print err.out
				pass

			if len(cfg['meta_order']) > 0:
				if models_obj[n].results['err'][0] == 0:
					meta_incl.append(n)

			if not model_written[n]:
				results_final_models[n] = models_obj[n].out
				results_final_models_headers[n] = models_obj[n].out.columns.values
				model_written[n] = True
			else:
				results_final_models[n] = results_final_models[n].append(models_obj[n].out, ignore_index=True)

			if len(cfg['meta_order']) > 0:
				models_obj[n].tag_results(n)

			if not region_written:
				results_region = models_obj[n].out
				region_written = True
			else:
				results_region = results_region.merge(models_obj[n].out, how='outer')

			analyzed = len(models_obj[n].variant_stats['filter'][models_obj[n].variant_stats['filter'] == 0])

			status = '   (' + n + ') processed ' + str(models_obj[n].variants.info.shape[0]) + ' variants, ' + str(analyzed) + ' passed filters'
			print status
			sys.stdout.flush()

		if len(cfg['meta_order']) > 0:
			for meta in cfg['meta_order']:
				meta_objs[meta].calc_meta(regions_df['chr'][k], regions_df['start'][k], regions_df['end'][k], regions_df['group_id'][k], models_obj[cfg['meta'][meta].split('+')[0]],meta_incl)
				print '   processed meta analysis ' + meta + ' (' + "+".join([x for x in cfg['meta'][meta].split('+') if x in meta_incl]) + ')'
				if not meta_written[meta]:
					results_final_meta[meta] = meta_objs[meta].out.copy()
					meta_written[meta] = True
				else:
					results_final_meta[meta] = results_final_meta[meta].merge(meta_objs[meta].out, how='outer')
		last_chr = regions_df['chr'][k]

	for n in cfg['model_order']:
		pkl = open('/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(cpu) + '.' + n + '.pkl', "wb")
		pickle.dump([results_final_models[n].sort_values(by=['chr','start']),models_obj[n].metadata,results_final_models_headers[n],models_obj[n].tbx_start,models_obj[n].tbx_end],pkl,protocol=2)
		pkl.close()

	if len(cfg['meta_order']) > 0:
		for meta in cfg['meta_order']:
			results_final_meta[meta] = results_final_meta[meta].sort_values(by=['chr','start'])
			results_final_meta[meta]['chr'] = results_final_meta[meta]['chr'].astype(np.int64)
			results_final_meta[meta]['start'] = results_final_meta[meta]['start'].astype(np.int64)
			results_final_meta[meta]['end'] = results_final_meta[meta]['end'].astype(np.int64)
			pkl = open('/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(cpu) + '.' + meta + '.pkl', "wb")
			pickle.dump([results_final_meta[meta],meta_objs[meta].metadata,np.array(results_final_meta[meta].columns.values),meta_objs[meta].tbx_start,meta_objs[meta].tbx_end],pkl,protocol=2)
			pkl.close()

	if log:
		sys.stdout = stdout_orig
		log_file.close()

	if variants_found:
		return 0
	else:
		return -1

def RunSnvgroup(args):
	cfg = Parse.generate_snvgroup_cfg(args)
	Parse.print_snvgroup_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	regions_df = pd.read_table(cfg['region_file'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
	regions_df = regions_df[regions_df['job'] == int(cfg['job'])].reset_index(drop=True)
	return_values = {}
	models_out = {}
	bgzfiles = {}
	print ''
	for m in cfg['model_order']:
		print "initializing out file for model " + m
		models_out[m] = cfg['out'] if m == '___no_tag___' else cfg['out'] + '.' + m
		try:
			bgzfiles[m] = bgzf.BgzfWriter(models_out[m] + '.gz', 'wb')
		except:
			print Process.Error("failed to initialize bgzip format out file " + models_out[m] + '.gz').out
			return 1
	if len(cfg['meta_order']) > 0:
		for m in cfg['meta_order']:
			print "initializing out file for meta " + m
			models_out[m] = cfg['out'] + '.' + m
			try:
				bgzfiles[m] = bgzf.BgzfWriter(models_out[m] + '.gz', 'wb')
			except:
				print Process.Error("failed to initialize bgzip format out file " + models_out[m] + '.gz').out
				return 1

	if cfg['cpus'] > 1:
		pool = mp.Pool(cfg['cpus']-1)
		for i in xrange(1,cfg['cpus']):
			return_values[i] = pool.apply_async(process_regions, args=(regions_df,cfg,i,True,))
			print "submitting job on cpu " + str(i) + " of " + str(cfg['cpus'])
		pool.close()
		print "executing job for cpu " + str(cfg['cpus']) + " of " + str(cfg['cpus']) + " via main process"
		main_return = process_regions(regions_df,cfg,cfg['cpus'],True)
		pool.join()

		if 1 in [return_values[i].get() for i in return_values] or main_return == 1:
			print Process.Error("error detected, see log files").out
			return 1

	else:
		main_return = process_regions(regions_df,cfg,1,True)
		if main_return == 1:
			print Process.Error("error detected, see log files").out
			return 1

	for i in xrange(1,cfg['cpus']+1):
		try:
			logfile = open(cfg['out'] + '.cpu' + str(i) + '.log', 'r')
		except:
			print Process.Error("failed to initialize log file " + cfg['out'] + '.cpu' + str(i) + '.log').out
			return 1
		print logfile.read()
		logfile.close()
		os.remove(cfg['out'] + '.cpu' + str(i) + '.log')

	for m in cfg['model_order']:
		written = False
		for i in xrange(1,cfg['cpus']+1):
			out_model_cpu = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.' + m + '.pkl'
			pkl = open(out_model_cpu,"rb")
			results_final,metadata,results_header,tbx_start,tbx_end = pickle.load(pkl)
			if not written:
				bgzfiles[m].write(metadata)
				bgzfiles[m].write("\t".join(results_header) + '\n')
				written = True
			if results_final.shape[0] > 0:
				results_final.replace({'None': 'NA'}).to_csv(bgzfiles[m], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = results_header, append=True)
			pkl.close()
			os.remove(out_model_cpu)

		bgzfiles[m].close()

		print "indexing out file for model " + m if m != '___no_tag___' else "indexing out file"
		try:
			pysam.tabix_index(models_out[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Process.Error('failed to generate index for file ' + models_out[m] + '.gz').out
			return 1

	if len(cfg['meta_order']) > 0:
		for m in cfg['meta_order']:
			written = False
			for i in xrange(1,cfg['cpus']+1):
				out_model_meta = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.' + m + '.pkl'
				pkl = open(out_model_meta,"rb")
				results_final_meta,metadata,results_header,tbx_start,tbx_end = pickle.load(pkl)
				if not written:
					bgzfiles[m].write(metadata)
					bgzfiles[m].write('\t'.join(results_header) + '\n')
					written = True
				if results_final_meta.shape[0] > 0:
					results_final_meta.replace({'None': 'NA'}).to_csv(bgzfiles[m], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = results_header, append=True)
				pkl.close()
				os.remove(out_model_meta)

			bgzfiles[m].close()

			print "indexing out file for meta " + m
			try:
				pysam.tabix_index(models_out[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
			except:
				print Process.Error('failed to generate index for file ' + models_out[m] + '.gz').out
				return 1

	print "process complete"
	return 0
