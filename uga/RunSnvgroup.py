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
import Geno
import pysam
import Fxns
import time
from Bio import bgzf
from Process import Error
import multiprocessing as mp
import sys
import os
import resource

def process_regions(regions_df, cfg, cpu, log):
	regions_df = regions_df[regions_df['cpu'] == cpu].reset_index(drop=True)

	if log:
		try:
			log_file = open(cfg['out'] + '.cpu' + str(cpu) + '.log','w')
		except:
			print Error("unable to initialize log file " + cfg['out'] + '.cpu' + str(cpu) + '.log')
			return 1
		else:
			stdout_orig = sys.stdout
			sys.stdout = log_file

	models_obj = {}
	for n in cfg['model_order']:
		print "\nloading model " + n if n != '___no_tag___' else "\nloading model"
		try:
			models_obj[n] = getattr(Model,cfg['models'][n]['fxn'].capitalize())(fxn=cfg['models'][n]['fxn'],snvgroup_map=cfg['snvgroup_map'],formula=cfg['models'][n]['formula'],format=cfg['models'][n]['format'],
									skat_wts=cfg['models'][n]['skat_wts'],burden_wts=cfg['models'][n]['burden_wts'],skat_method=cfg['models'][n]['skat_method'],
									burden_mafrange=cfg['models'][n]['burden_mafrange'],
									all_founders=cfg['models'][n]['all_founders'],case_code=cfg['models'][n]['case_code'],ctrl_code=cfg['models'][n]['ctrl_code'],
									pheno_file=cfg['models'][n]['pheno'],variants_file=cfg['models'][n]['file'],type=cfg['models'][n]['fxn'],fid=cfg['models'][n]['fid'],
									iid=cfg['models'][n]['iid'],matid=cfg['models'][n]['matid'],patid=cfg['models'][n]['patid'],sex=cfg['models'][n]['sex'],
									male=cfg['models'][n]['male'],female=cfg['models'][n]['female'],pheno_sep=Fxns.GetDelimiter(cfg['models'][n]['sep']))
		except Error as err:
			print err.out
			return 1

	variants_found = False
	model_written = {}
	final_written = False
	results_final_models = {}
	results_final_models_headers = {}
	for n in cfg['model_order']:
		model_written[n] = False
		results_final_models[n] = pd.DataFrame({})
	results_final_meta = pd.DataFrame({})
	for k in xrange(len(regions_df.index)):
		meta_incl = []
		variants_db = {}
		region_written = False
		results_region = pd.DataFrame({})
		print ''
		print 'loading region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['id'][k] + ": " + regions_df['region'][k] + ') ...'
		for n in cfg['model_order']:
			try:
				models_obj[n].get_region(regions_df['region'][k], id=regions_df['id'][k])
			except:
				pass
			try:
				models_obj[n].get_snvgroup(cfg['buffer'], regions_df['id'][k])
			except:
				pass
			variants_found = True

			if n == cfg['model_order'][0]:
				ref = Geno.VariantRef(models_obj[n].variants)
			else:
				ref.update(models_obj[n].variants)
				models_obj[n].variants.align(ref)

			try:
				models_obj[n].filter(miss_thresh=cfg['models'][n]['miss'], maf_thresh=cfg['models'][n]['maf'], maxmaf_thresh=cfg['models'][n]['maxmaf'], 
								mac_thresh=cfg['models'][n]['mac'], rsq_thresh=cfg['models'][n]['rsq'], hwe_thresh=cfg['models'][n]['hwe'], 
								hwe_maf_thresh=cfg['models'][n]['hwe_maf'])
			except:
				pass

			try:
				models_obj[n].calc_model()
			except Error as err:
				print err.out
				pass

			if models_obj[n].results['err'][0] == 0:
				meta_incl.append(n)

			if not model_written[n]:
				results_final_models[n] = models_obj[n].out
				results_final_models_headers[n] = models_obj[n].out.columns.values
				model_written[n] = True
			else:
				results_final_models[n] = results_final_models[n].append(models_obj[n].out, ignore_index=True)

			if len(cfg['model_order']) > 1:
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

		if len(cfg['meta_order']) > 1:
			h1 = ['chr','start','end','id']
			h2 = [x for x in results_region if x not in ['chr','start','end','id']]
			for meta in cfg['meta_order']:
				meta_result = getattr(Model,cfg['models'][n]['fxn'].capitalize() + 'Meta')(models_obj[cfg['model_order'][0]], meta, cfg['meta'][meta], meta_incl)
				h1 = h1 + [x for x in meta_result.columns.values]
				results_region = pd.concat([results_region,meta_result], axis=1)
			h = h1 + h2
			if not final_written:
				results_final_meta = results_region[h]
				final_written = True
			else:
				results_final_meta = results_final_meta.merge(results_region[h], how='outer')

	for n in cfg['model_order']:
		store = pd.HDFStore('/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(cpu) + '.' + n + '.h5')
		store.put('df',results_final_models[n].sort(columns=['chr','start']))
		store.get_storer('df').attrs.metadata = models_obj[n].metadata
		store.get_storer('df').attrs.results_header = results_final_models_headers[n]
		store.get_storer('df').attrs.tbx_start = models_obj[n].tbx_start
		store.get_storer('df').attrs.tbx_end = models_obj[n].tbx_end
		store.close()

	if len(cfg['meta_order']) > 1:
		results_final_meta = results_final_meta.sort(columns=['chr','start'])
		results_final_meta['chr'] = results_final_meta['chr'].astype(np.int64)
		results_final_meta['start'] = results_final_meta['start'].astype(np.int64)
		results_final_meta['end'] = results_final_meta['end'].astype(np.int64)
		store = pd.HDFStore('/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(cpu) + '.meta.h5')
		store.put('df',results_final_meta)
		store.get_storer('df').attrs.results_header = np.array(results_final_meta.columns.values)
		store.get_storer('df').attrs.tbx_start = models_obj[cfg['model_order'][0]].tbx_start
		store.get_storer('df').attrs.tbx_end = models_obj[cfg['model_order'][0]].tbx_end
		store.close()

	if log:
		sys.stdout = stdout_orig
		log_file.close()

	if variants_found:
		return 0
	else:
		return -1

def RunSnvgroup(args):
	cfg = Parse.GenerateSnvgroupCfg(args)
	Parse.PrintSnvgroupOptions(cfg)

	regions_df = pd.read_table(cfg['region_file'])
	return_values = {}
	models_out = {}
	bgzfiles = {}
	for m in cfg['model_order']:
		print "initializing out file for model " + m + " results"
		models_out[m] = cfg['out'] if m == '___no_tag___' else cfg['out'] + '.' + m
		try:
			bgzfiles[m] = bgzf.BgzfWriter(models_out[m] + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + models_out[m] + '.gz').out
			return 1
	if len(cfg['meta_order']) > 0:
		print "initializing out file for meta results"
		models_out['___meta___'] = cfg['out'] + '.meta'
		try:
			bgzfiles['___meta___'] = bgzf.BgzfWriter(models_out['___meta___'] + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + models_out['___meta___'] + '.gz').out
			return 1

	if cfg['cpus'] > 1:
		pool = mp.Pool(cfg['cpus'])
		for i in xrange(1,cfg['cpus']+1):
			return_values[i] = pool.apply_async(process_regions, args=(regions_df,cfg,i,True,))
			print "submitting job on cpu " + str(i) + " of " + str(cfg['cpus'])
		pool.close()
		pool.join()

		if 1 in [return_values[i].get() for i in return_values]:
			print Error("error detected, see log files").out
			return 1

	else:
		return_values[1] = process_regions(regions_df,cfg,1,True)
		if return_values[1] == -1:
			print Error("error detected, see log files").out
			return 1

	for i in xrange(1,cfg['cpus']+1):
		try:
			logfile = open(cfg['out'] + '.cpu' + str(i) + '.log', 'r')
		except:
			print Error("failed to initialize log file " + cfg['out'] + '.cpu' + str(i) + '.log').out
			return 1
		print logfile.read()
		logfile.close()
		os.remove(cfg['out'] + '.cpu' + str(i) + '.log')

	for m in cfg['model_order']:
		written = False
		for i in xrange(1,cfg['cpus']+1):
			out_model_cpu = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.' + m + '.h5'
			store = pd.HDFStore(out_model_cpu)
			if not written:
				bgzfiles[m].write(store.get_storer('df').attrs.metadata)
				bgzfiles[m].write('\t'.join(store.get_storer('df').attrs.results_header) + '\n')
				tbx_start = store.get_storer('df').attrs.tbx_start
				tbx_end = store.get_storer('df').attrs.tbx_end
				written = True
			if store['df'].shape[0] > 0:
				store['df'].replace({'None': 'NA'}).to_csv(bgzfiles[m], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.results_header, append=True)
			store.close()
			os.remove(out_model_cpu)

		bgzfiles[m].close()

		print "indexing out file for model " + m
		try:
			pysam.tabix_index(models_out[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + models_out[m] + '.gz').out
			return 1

	if len(cfg['meta_order']) > 1:
		written = False
		for i in xrange(1,cfg['cpus']+1):
			out_model_meta = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.meta.h5'
			store = pd.HDFStore(out_model_meta)
			if not written:
				bgzfiles['___meta___'].write('#' + '\t'.join(store.get_storer('df').attrs.results_header) + '\n')
				tbx_start = store.get_storer('df').attrs.tbx_start
				tbx_end = store.get_storer('df').attrs.tbx_end
				written = True
			if store['df'].shape[0] > 0:
				store['df'].replace({'None': 'NA'}).to_csv(bgzfiles['___meta___'], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.results_header, append=True)
			store.close()
			os.remove(out_model_meta)

		bgzfiles['___meta___'].close()

		print "indexing out file for meta results"
		try:
			pysam.tabix_index(models_out['___meta___'] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + models_out['___meta___'] + '.gz').out
			return 1

	print "process complete"
	return 0
