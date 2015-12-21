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
import numpy.lib.recfunctions as recfxns
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

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunSnv")

def process_regions(regions_df, cfg, cpu, log):
	regions_df = regions_df[regions_df['cpu'] == cpu].reset_index(drop=True)

	if log:
		try:
			log_file = open(cfg['out'] + '.cpu' + str(cpu) + '.log','w')
		except:
			print Process.Error("unable to initialize log file " + cfg['out'] + '.cpu' + str(cpu) + '.log').out
			return 1
		else:
			stdout_orig = sys.stdout
			stderr_orig = sys.stderr
			sys.stdout = log_file
			sys.stderr = log_file

	models_obj = {}
	for n in cfg['model_order']:
		print "\nloading model " + n if n != '___no_tag___' else "\nloading model"
		try:
			models_obj[n] = getattr(Model,cfg['models'][n]['fxn'].capitalize())(fxn=cfg['models'][n]['fxn'], 
																				format=cfg['models'][n]['format'], 
																				corstr=cfg['models'][n]['corstr'], 
																				pheno=cfg['models'][n]['pheno'], 
																				covars=cfg['models'][n]['covars'], 
																				covars_categorical=cfg['models'][n]['covars_categorical'], 
																				interact=cfg['models'][n]['interact'], 
																				reverse=cfg['models'][n]['reverse'], 
																				all_founders=cfg['models'][n]['all_founders'], 
																				case_code=cfg['models'][n]['case_code'], 
																				ctrl_code=cfg['models'][n]['ctrl_code'], 
																				ped=cfg['models'][n]['ped'], 
																				variants_file=cfg['models'][n]['file'], 
																				samples_file=cfg['models'][n]['sample'], 
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

	variant_ref = Variant.Ref()
	for n in cfg['model_order']:
		written = False
		models_obj[n].out_all = pd.DataFrame({})
		print ''
		print 'loading regions for model ' + n + ' ...'
		for k in xrange(len(regions_df.index)):
			variants_found = False
			i = 0
			try:
				models_obj[n].get_region(regions_df['region'][k])
			except Process.Error as err:
				print err.out
				pass
			else:
				while True:
					i = i + 1
					try:
						models_obj[n].get_snvs(cfg['buffer'])
					except:
						if not variants_found:
							print '   processed 0 variants in region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['region'][k] + ')'
						break
					variants_found = True

					if len(cfg['meta_order']) > 0:
						if n == cfg['model_order'][0]:
							variant_ref.load(models_obj[n].variants.info)
						else:
							variant_ref.update(models_obj[n].variants.info)
							models_obj[n].variants.align(variant_ref)

					try:
						models_obj[n].filter(miss_thresh=cfg['models'][n]['miss'], maf_thresh=cfg['models'][n]['maf'], maxmaf_thresh=cfg['models'][n]['maxmaf'], 
										mac_thresh=cfg['models'][n]['mac'], rsq_thresh=cfg['models'][n]['rsq'], hwe_thresh=cfg['models'][n]['hwe'], 
										hwe_maf_thresh=cfg['models'][n]['hwe_maf'])
					except:
						break
					try:
						models_obj[n].calc_model()
					except Process.Error as err:
						print err.out
						break
					if not written:
						models_obj[n].out_all = models_obj[n].out
						written = True
					else:
						models_obj[n].out_all = models_obj[n].out_all.append(models_obj[n].out, ignore_index=True)
					analyzed = len(models_obj[n].variant_stats['filter'][models_obj[n].variant_stats['filter'] == 0])
					cur_variants = min(i*cfg['buffer'],(i-1)*cfg['buffer'] + models_obj[n].variants.info.shape[0])
					status = '   processed ' + str(cur_variants) + ' variants in region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['region'][k] + '), ' + str(analyzed) + ' passed filters'
					print status
					sys.stdout.flush()

		pkl = open('/'.join(cfg['out'].split('/')[0:-1]) + '/' + (cfg['out'] + '.cpu' + str(cpu) + '.' + n).split('/')[-1] + '.pkl', "wb")
		pickle.dump([models_obj[n].out_all,models_obj[n].metadata,models_obj[n].results_header,models_obj[n].tbx_start,models_obj[n].tbx_end],pkl,protocol=2)
		pkl.close()

	results_all = None
	results_all_written = False
	results_all_header = {0: ['chr','pos','id','a1','a2']}
	print ''
	if len(cfg['meta_order']) > 0:
		print "preparing data for meta analysis ..."
		for n in cfg['model_order']:
			models_obj[n].out_all = models_obj[n].out_all[[x for x in models_obj[n].out_all.columns if x not in ['group_id','id_unique','uid']]]
			models_obj[n].out_all.columns = np.array([n + '.' + x if x not in ['chr','pos','id','a1','a2'] else x for x in models_obj[n].out_all.columns])
			models_obj[n].out_all[n + '.n'] = (models_obj[n].out_all[n + '.callrate'].astype(float) * models_obj[n].nunique).round()
			if n == cfg['model_order'][0]:
				results_all = models_obj[n].out_all
			else:
				results_all = results_all.merge(models_obj[n].out_all,how='outer',on=['chr','pos','id','a1','a2'])
		results_all_header[2] = list(results_all.columns[5:])
		results_all_header[1] = []
		for meta in cfg['meta_order']:
			new_cols, results_all = Model.SnvMeta(results_all, meta, cfg['meta'][meta].split('+'), meta_type = cfg['meta_type'])
			results_all_header[1].extend(new_cols)
			print "   processed meta analysis (" + meta + ")"
		results_all = results_all.sort_values(by=['chr','pos'])
		results_all['chr'] = results_all['chr'].astype(np.int64)
		results_all['pos'] = results_all['pos'].astype(np.int64)
		pkl = open('/'.join(cfg['out'].split('/')[0:-1]) + '/' + (cfg['out'] + '.cpu' + str(cpu)).split('/')[-1] + '.meta.pkl', "wb")
		pickle.dump([results_all,np.array(results_all_header[0] + results_all_header[1] + results_all_header[2]),models_obj[cfg['model_order'][0]].tbx_start,models_obj[cfg['model_order'][0]].tbx_end],pkl,protocol=2)
		pkl.close()

	if log:
		sys.stdout = stdout_orig
		sys.stderr = stderr_orig
		log_file.close()

	if variants_found:
		return 0
	else:
		return -1

def RunSnv(args):
	cfg = Parse.generate_snv_cfg(args)
	Parse.print_snv_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	regions_df = pd.read_table(cfg['region_file'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
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
		print "initializing out file for meta results"
		models_out['___meta___'] = cfg['out'] + '.meta'
		try:
			bgzfiles['___meta___'] = bgzf.BgzfWriter(models_out['___meta___'] + '.gz', 'wb')
		except:
			print Process.Error("failed to initialize bgzip format out file " + models_out['___meta___'] + '.gz').out
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
			regions_cpu_df = regions_df[regions_df['cpu'] == i].reset_index(drop=True)
			out_model_range = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + (cfg['out'] + '.cpu' + str(i) + '.' + m).split('/')[-1] + '.pkl'
			pkl = open(out_model_range,"rb")
			results_final,metadata,results_header,tbx_start,tbx_end = pickle.load(pkl)
			if not written:
				bgzfiles[m].write(metadata)
				bgzfiles[m].write("\t".join(results_header) + '\n')
				written = True
			if results_final.shape[0] > 0:
				results_final.replace({'None': 'NA'}).to_csv(bgzfiles[m], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = results_header, append=True)
			pkl.close()
			os.remove(out_model_range)

		bgzfiles[m].close()
		print "indexing out file for model " + m
		try:
			pysam.tabix_index(models_out[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Process.Error('failed to generate index for file ' + models_out[m] + '.gz').out
			return 1

	if len(cfg['meta_order']) > 0:
		written = False
		for i in xrange(1,cfg['cpus']+1):
			out_model_meta = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.meta.pkl'
			pkl = open(out_model_meta,"rb")
			results_final_meta,results_header,tbx_start,tbx_end = pickle.load(pkl)
			if not written:
				bgzfiles['___meta___'].write('#' + '\t'.join(results_header) + '\n')
				written = True
			if results_final_meta.shape[0] > 0:
				results_final_meta.replace({'None': 'NA'}).to_csv(bgzfiles['___meta___'], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = results_header, append=True)
			pkl.close()
			os.remove(out_model_meta)

		bgzfiles['___meta___'].close()

		print "indexing out file for meta results"
		try:
			pysam.tabix_index(models_out['___meta___'] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Process.Error('failed to generate index for file ' + models_out['___meta___'] + '.gz').out
			return 1

	print "process complete"
	return 0
