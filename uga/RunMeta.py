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
import Geno
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
logger = logging.getLogger("RunMeta")

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

	results_obj = {}
	for f in cfg['file_order']:
		print "\nloading results file " + f
		try:
			results_obj[f] = Geno.Results(filename=cfg['files'][f])
		except Process.Error as err:
			print err.out
			return 1

	variants_found = False
	meta_written = {}
	results_final_meta = {}
	meta_objs = {}
	variant_ref = Variant.Ref()
	for meta in cfg['meta_order']:
		meta_written[meta] = False
		results_final_meta[meta] = pd.DataFrame({})
		meta_objs[meta] = Model.SnvMeta(tag = meta, meta = cfg['meta'][meta], type = cfg['meta_type'][meta])
	for k in xrange(len(regions_df.index)):
		region_written = False
		print ''
		print 'loading region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['region'][k] + ') ...'
		for f in cfg['file_order']:
			try:
				results_obj[f].get_region(regions_df['region'][k])
			except:
				pass

			try:
				results_obj[f].get_snvs(cfg['buffer'])
			except:
				pass
			variants_found = True

			if f == cfg['file_order'][0]:
				variant_ref.load(results_obj[f].snv_results)
			else:
				variant_ref.update(results_obj[f].snv_results)
				results_obj[f].align_results(variant_ref)

			results_obj[f].tag_results(f)

			if not region_written:
				results_region = results_obj[f].snv_results_tagged.copy()
				region_written = True
			else:
				results_region_cols = [x for x in results_region.columns.values] + [x for x in results_obj[f].snv_results_tagged.columns.values if x not in results_region.columns.values]
				if results_region.empty and not results_obj[f].snv_results_tagged.empty:
					results_region=pd.concat([results_obj[f].snv_results_tagged[['chr','pos','id','a1','a2','id_unique','___uid___']].iloc[[0]],pd.DataFrame(dict(zip([x for x in results_region.columns.values if x not in ['chr','pos','id','a1','a2','id_unique','___uid___']],[np.nan for x in results_region.columns.values if x not in ['chr','pos','id','a1','a2','id_unique','___uid___']])),index=[0])],axis=1)
				results_region = results_region.merge(results_obj[f].snv_results_tagged, how='outer')			
				results_region = results_region[results_region_cols]

			status = '   (' + f + ') processed ' + str(results_obj[f].snv_results.shape[0]) + ' variants'
			print status
			sys.stdout.flush()

		for meta in cfg['meta_order']:
			meta_objs[meta].calc_meta(results_region)
			print '   processed meta analysis ' + meta + ' (' + cfg['meta'][meta] + ')'
			if not meta_written[meta]:
				results_final_meta[meta] = meta_objs[meta].out.copy()
				meta_written[meta] = True
			else:
				results_final_meta[meta] = results_final_meta[meta].merge(meta_objs[meta].out, how='outer')

	for meta in cfg['meta_order']:
		results_final_meta[meta] = results_final_meta[meta].sort_values(by=['chr','pos'])
		results_final_meta[meta]['chr'] = results_final_meta[meta]['chr'].astype(np.int64)
		results_final_meta[meta]['pos'] = results_final_meta[meta]['pos'].astype(np.int64)
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

def RunMeta(args):
	cfg = Parse.generate_meta_cfg(args)
	Parse.print_meta_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	regions_df = pd.read_table(cfg['region_file'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
	regions_df = regions_df[regions_df['job'] == int(cfg['job'])].reset_index(drop=True)
	return_values = {}
	meta_out = {}
	bgzfiles = {}
	print ''
	for m in cfg['meta_order']:
		print "initializing out file for meta " + m
		meta_out[m] = cfg['out'] + '.' + m
		try:
			bgzfiles[m] = bgzf.BgzfWriter(meta_out[m] + '.gz', 'wb')
		except:
			print Process.Error("failed to initialize bgzip format out file " + meta_out[m] + '.gz').out
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
			pysam.tabix_index(meta_out[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Process.Error('failed to generate index for file ' + meta_out[m] + '.gz').out
			return 1

	print "process complete"
	return 0
