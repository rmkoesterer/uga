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
import pysam
import Fxns
import time
from Bio import bgzf
from Process import Error
import multiprocessing as mp
import sys
import os
import resource

def process_regions(regions_df, cfg, models, cpu, log):
	regions_df = regions_df[regions_df['cpu'] == cpu].reset_index(drop=True)

	for n in cfg['model_order']:
		model_out = models[n] + '.cpu' + str(cpu)
		if log:
			try:
				log_file = open(model_out + '.log','w')
			except:
				print Error("unable to initialize log file " + out + '.log')
				return 1
			else:
				stdout_orig = sys.stdout
				sys.stdout = log_file

		print "\nloading model " + n if n != '___no_tag___' else "\nloading model"
		try:
			m = getattr(Model,cfg['models'][n]['fxn'].capitalize())(formula=cfg['models'][n]['formula'],format=cfg['models'][n]['format'], 
									all_founders=cfg['models'][n]['all_founders'],case_code=cfg['models'][n]['case_code'],ctrl_code=cfg['models'][n]['ctrl_code'],
									pheno_file=cfg['models'][n]['pheno'],biodata_file=cfg['models'][n]['file'],type=cfg['models'][n]['fxn'],fid=cfg['models'][n]['fid'],
									iid=cfg['models'][n]['iid'],matid=cfg['models'][n]['matid'],patid=cfg['models'][n]['patid'],sex=cfg['models'][n]['sex'],
									pheno_sep=Fxns.GetDelimiter(cfg['models'][n]['sep']))
		except Error as err:
			print err.msg
			return 1

		variants_found = False
		for k in xrange(len(regions_df.index)):
			i = 0
			written = False
			results_final = pd.DataFrame({})
			print 'loading region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['gene'][k] + ":" + regions_df['region'][k] + ')',

			try:
				m.get_region(regions_df['region'][k])
			except:
				print " <-- chromosome not found"
				pass
			else:
				print '...'
				while True:
					i = i + 1
					try:
						m.get_chunk(cfg['buffer'])
					except:
						break
					variants_found = True

					try:
						m.filter(miss_thresh=cfg['models'][n]['miss'], maf_thresh=cfg['models'][n]['maf'], maxmaf_thresh=cfg['models'][n]['maxmaf'], 
										mac_thresh=cfg['models'][n]['mac'], rsq_thresh=cfg['models'][n]['rsq'], hwe_thresh=cfg['models'][n]['hwe'], 
										hwe_maf_thresh=cfg['models'][n]['hwe_maf'])
					except:
						break

					try:
						m.calc_model()
					except Error as err:
						print err.msg
						break

					if not written:
						results_final = m.out
						written = True
					else:
						results_final = results_final.append(m.out, ignore_index=True)
					analyzed = len(m.marker_stats['filter'][m.marker_stats['filter'] == 0])

					cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + m.biodata.marker_info.shape[0]))
					status = '   processed ' + cur_markers + ' variants, ' + str(analyzed) + ' passed filters'
					print status

			store = pd.HDFStore('/'.join(model_out.split('/')[0:-1]) + '/chr' + str(regions_df['chr'][k]) + '/' + model_out.split('/')[-1] + '.chr' + str(regions_df['chr'][k]) + 'bp' + str(regions_df['start'][k]) + '-' + str(regions_df['end'][k]) + '.h5')
			store.put('df',results_final)
			store.get_storer('df').attrs.results_header_metadata = m.results_header_metadata
			store.get_storer('df').attrs.results_header = m.results_header
			store.get_storer('df').attrs.tbx_start = m.tbx_start
			store.get_storer('df').attrs.tbx_end = m.tbx_end
			store.close()
		if log:
			sys.stdout = stdout_orig
			log_file.close()

	if variants_found:
		return 0
	else:
		return -1

def RunGene(args):
	cfg = Parse.GenerateGeneCfg(args)
	Parse.PrintGeneOptions(cfg)

	regions_df = pd.read_table(cfg['region_file'])
	print regions_df
	return_values = {}
	models = {}
	bgzfiles = {}
	for m in cfg['model_order']:
		print "initializing out file for model " + m
		models[m] = cfg['out'] if m == '___no_tag___' else cfg['out'] + '.' + m
		try:
			bgzfiles[m] = bgzf.BgzfWriter(models[m] + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + models[m] + '.gz').msg
			return 1

	if cfg['cpus'] > 1:
		pool = mp.Pool(cfg['cpus'])
		for i in xrange(1,cfg['cpus']+1):
			return_values[i] = pool.apply_async(process_regions, args=(regions_df,cfg,models,i,True,))
			print "submitting job on cpu " + str(i) + " of " + str(cfg['cpus'])
		pool.close()
		pool.join()

		if 1 in [return_values[i].get() for i in return_values]:
			print Error("error detected, see log files").msg
			return 1

	else:
		return_values[1] = process_regions(regions_df,cfg,models,1,True)
		if return_values[1] == -1:
			print Error("error detected, see log files").msg
			return 1

	for m in cfg['model_order']:
		for i in xrange(1,cfg['cpus']+1):
			try:
				logfile = open(models[m] + '.cpu' + str(i) + '.log', 'r')
			except:
				print Error("failed to initialize log file " + models[m] + '.cpu' + str(i) + '.log').msg
				return 1
			print logfile.read()
			logfile.close()
			os.remove(models[m] + '.cpu' + str(i) + '.log')

		written = False
		for i in xrange(1,cfg['cpus']+1):
			regions_cpu_df = regions_df[regions_df['cpu'] == i].reset_index(drop=True)
			for j in xrange(len(regions_cpu_df.index)):
				out_model_range = '/'.join(models[m].split('/')[0:-1]) + '/chr' + str(regions_cpu_df['chr'][j]) + '/' + models[m].split('/')[-1] + '.cpu' + str(i) + '.chr' + str(regions_cpu_df['chr'][j]) + 'bp' + str(regions_cpu_df['start'][j]) + '-' + str(regions_cpu_df['end'][j]) + '.h5'
				store = pd.HDFStore(out_model_range)
				if not written:
					bgzfiles[m].write(store.get_storer('df').attrs.results_header_metadata)
					bgzfiles[m].write('\t'.join(store.get_storer('df').attrs.results_header) + '\n')
					tbx_start = store.get_storer('df').attrs.tbx_start
					tbx_end = store.get_storer('df').attrs.tbx_end
					written = True
				if store['df'].shape[0] > 0:
					store['df'].replace({'None': 'NA'}).to_csv(bgzfiles[m], index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.results_header, append=True)
				store.close()
				os.remove(out_model_range)

		bgzfiles[m].close()
		print "indexing out file for model " + m
		try:
			pysam.tabix_index(models[m] + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + models[m] + '.gz').msg
			return 1

	print "process complete"
	return 0
