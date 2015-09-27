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
import Model
import IO
import Parse
import pysam
import Fxns
import time
from Bio import bgzf
from Process import Error
import multiprocessing as mp
import sys
import os

def process_regions(regions, cfg, out, r):
	try:
		log_file = open(out + '.log','w')
	except:
		print Error("unable to initialize log file " + out + '.log')
		return 1
	else:
		stdout_orig = sys.stdout
		sys.stdout = log_file

	for k in cfg['model_order']:
		out_model = out if k == '___no_tag___' else out + '.' + k

		##### GENERATE MODEL OBJECT #####
		print "\nloading model " + k if k != '___no_tag___' else "\nloading model"
		try:
			model = getattr(Model,cfg['models'][k]['fxn'].capitalize())(formula=cfg['models'][k]['formula'],format=cfg['models'][k]['format'],variant_list_file=cfg['models'][k]['variant_list'],
									all_founders=cfg['models'][k]['all_founders'],case_code=cfg['models'][k]['case_code'],ctrl_code=cfg['models'][k]['ctrl_code'],
									pheno_file=cfg['models'][k]['pheno'],biodata_file=cfg['models'][k]['file'],type=cfg['models'][k]['fxn'],fid=cfg['models'][k]['fid'],
									iid=cfg['models'][k]['iid'],matid=cfg['models'][k]['matid'],patid=cfg['models'][k]['patid'],sex=cfg['models'][k]['sex'],
									pheno_sep=Fxns.GetDelimiter(cfg['models'][k]['sep']))
		except Error as err:
			print err.msg
			return 1

		##### INITIALIZE OUTPUT FILE HANDLE #####
		store = pd.HDFStore(out_model + '.h5')
	
		##### RUN ANALYSIS #####
		print "running models ..."
		analyzed = 0
		variants_found = False
		i = 0
		print 'loading region ' + str(r+1) + '/' + str(len(regions.df.index)) + ' (' + regions.df['region'][r] + ')' if regions.df['id'] is not None else 'loading region ' + str(r+1) + '/' + str(len(regions.df.index)) + ' (' + regions.df['id'][r] + ': ' + regions.df['region'][r] + ')',
		try:
			model.get_region(regions.df['region'][r], regions.df['id'][r])
		except:
			print " <-- no variants found"
			pass
		else:
			print '...'
			variants_found = True
			while True:
				i = i + 1
				try:
					model.get_chunk(cfg['buffer'])
				except:
					break

				try:
					model.filter(miss_thresh=cfg['models'][k]['miss'], maf_thresh=cfg['models'][k]['maf'], maxmaf_thresh=cfg['models'][k]['maxmaf'], 
									mac_thresh=cfg['models'][k]['mac'], rsq_thresh=cfg['models'][k]['rsq'], hwe_thresh=cfg['models'][k]['hwe'], 
									hwe_maf_thresh=cfg['models'][k]['hwe_maf'])
				except:
					break

				try:
					model.calc_model()
				except Error as err:
					print err.msg
					break

				store.put('df',model.out)
				store.get_storer('df').attrs.metadata = model.results_header_metadata
				store.get_storer('df').attrs.header = model.results_header
				store.get_storer('df').attrs.tbx_start = model.tbx_start
				store.get_storer('df').attrs.tbx_end = model.tbx_end
				analyzed += len(model.marker_stats['filter'][model.marker_stats['filter'] == 0])

				cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + model.biodata.marker_info.shape[0]))
				status = '   processed ' + cur_markers + ' variants, ' + str(analyzed) + ' passed filters'
				print status
	
		store.close()
	sys.stdout = stdout_orig
	log_file.close()
	if variants_found:
		return 0
	else:
		return -1

def RunModels(args):

	##### PARSE ARGS #####
	cfg = Parse.GenerateModelCfg(args)

	##### PRINT OPTIONS #####
	Parse.PrintModelOptions(cfg)

	##### LOAD REGION LIST #####
	try:
		regions = IO.Regions(filename=cfg['region_list'],region=cfg['region'],id=cfg['id'])
	except Error as err:
		print err.msg
		return 1

	cpus = mp.cpu_count() if mp.cpu_count() < int(cfg['cpus']) else int(cfg['cpus'])

	pool = mp.Pool(cpus)
	return_values = {}
	for r in xrange(regions.df.shape[0]):
		return_values[r] = pool.apply_async(process_regions, args=(regions,cfg,cfg['out'] + '.chr' + regions.df['region'][r].replace(':','bp'),r,))
		print "submitted cpu " + str(r) + " of " + str(cpus)
	pool.close()
	pool.join()
	

	##### figure out how to check return values #####
	print [return_values[i].get() for i in return_values]
	if 1 in [return_values[i].get() for i in return_values]:
		print Error("error detected").msg
		return 1

	##### COLLATE LOG FILES #####
	print "collating logs"
	for r in xrange(0,regions.df.shape[0]):
		region_log = cfg['out'] + '.chr' + regions.df['region'][r].replace(':','bp')
		try:
			logfile = open(region_log + '.log', 'r')
		except:
			print Error("failed to initialize log file " + region_log + '.log').msg
			return 1
		print logfile.read()
		os.remove(region_log + '.log')

	##### COLLATE RESULTS #####
	print "collating results"
	for k in cfg['model_order']:
		i = 0
		print "initializing out file for model " + k
		out = cfg['out'] if k == '___no_tag___' else cfg['out'] + '.' + k
		try:
			bgzfile = bgzf.BgzfWriter(out + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + out + '.gz').msg
			return 1
		out_model = cfg['out'] + '.chr' + regions.df['region'][0].replace(':','bp') if k == '___no_tag___' else cfg['out'] + '.chr' + regions.df['region'][0].replace(':','bp') + '.' + k
		store = pd.HDFStore(out_model + '.h5')
		bgzfile.write(store.get_storer('df').attrs.metadata)
		store['df'].to_csv(bgzfile, index=False, sep='\t', header=True, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.header)
		tbx_start = store.get_storer('df').attrs.tbx_start
		tbx_end = store.get_storer('df').attrs.tbx_end
		store.close()
		os.remove(out_model + '.h5')
		for r in xrange(1,regions.df.shape[0]):
			i += 1
			out_model = cfg['out'] + '.chr' + regions.df['region'][r].replace(':','bp') if k == '___no_tag___' else cfg['out'] + '.chr' + regions.df['region'][r].replace(':','bp') + '.' + k
			if return_values[r].get() != -1:
				store = pd.HDFStore(out_model + '.h5')
				store['df'].to_csv(bgzfile, index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.header)
				store.close()
			os.remove(out_model + '.h5')
		bgzfile.close()
		logfile.close()
		print "indexing out file for model " + k
		try:
			pysam.tabix_index(out + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + out + '.gz').msg
			return 1

	print "process complete"
	return 0
