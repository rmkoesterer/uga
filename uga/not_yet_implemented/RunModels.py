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
import resource

def process_regions(regions, cfg, k, out, r, log):
	regions.df = regions.df[regions.df['cpu'] == r].reset_index(drop=True)

	if log:
		try:
			log_file = open(out + '.log','w')
		except:
			print Error("unable to initialize log file " + out + '.log')
			return 1
		else:
			stdout_orig = sys.stdout
			sys.stdout = log_file

	print "\nloading model " + k if k != '___no_tag___' else "\nloading model"
	try:
		m = getattr(Model,cfg['models'][k]['fxn'].capitalize())(formula=cfg['models'][k]['formula'],format=cfg['models'][k]['format'],variant_list_file=cfg['models'][k]['variant_list'],
								all_founders=cfg['models'][k]['all_founders'],case_code=cfg['models'][k]['case_code'],ctrl_code=cfg['models'][k]['ctrl_code'],
								pheno_file=cfg['models'][k]['pheno'],biodata_file=cfg['models'][k]['file'],type=cfg['models'][k]['fxn'],fid=cfg['models'][k]['fid'],
								iid=cfg['models'][k]['iid'],matid=cfg['models'][k]['matid'],patid=cfg['models'][k]['patid'],sex=cfg['models'][k]['sex'],
								pheno_sep=Fxns.GetDelimiter(cfg['models'][k]['sep']))
	except Error as err:
		print err.msg
		return 1


	##### RUN ANALYSIS #####
	analyzed = 0
	variants_found = False
	for r in xrange(len(regions.df.index)):
		i = 0
		written = False
		results_final = pd.DataFrame({})
		print 'loading region ' + str(r+1) + '/' + str(len(regions.df.index)) + ' (' + regions.df['region'][r] + ')' if regions.df['id'] is not None else 'loading region ' + str(r+1) + '/' + str(len(regions.df.index)) + ' (' + regions.df['id'][r] + ': ' + regions.df['region'][r] + ')',

		try:
			m.get_region(regions.df['region'][r], regions.df['id'][r])
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
					m.filter(miss_thresh=cfg['models'][k]['miss'], maf_thresh=cfg['models'][k]['maf'], maxmaf_thresh=cfg['models'][k]['maxmaf'], 
									mac_thresh=cfg['models'][k]['mac'], rsq_thresh=cfg['models'][k]['rsq'], hwe_thresh=cfg['models'][k]['hwe'], 
									hwe_maf_thresh=cfg['models'][k]['hwe_maf'])
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
				analyzed += len(m.marker_stats['filter'][m.marker_stats['filter'] == 0])

				cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + m.biodata.marker_info.shape[0]))
				status = '   processed ' + cur_markers + ' variants, ' + str(analyzed) + ' passed filters'
				print status

		store = pd.HDFStore(out + '.chr' + str(regions.df['chr'][r]) + 'bp' + str(regions.df['start'][r]) + '-' + str(regions.df['end'][r]) + '.h5')
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

	##### CALCULATE REGION SPLIT FACTOR FOR CPUS #####
	cpus = mp.cpu_count() if mp.cpu_count() < int(cfg['cpus']) else int(cfg['cpus'])

	##### SPLIT REGIONS #####
	regions.df['cpu'] = 1
	if cpus > 1:
		if cpus >= regions.df.shape[0]:
			regions.df['cpu'] = np.array(range(regions.df.shape[0]))
		else:
			n = int(np.ceil(regions.df.shape[0] / float(cpus)))
			n_remain = int(regions.df.shape[0] - (n-1) * cpus)
			regions.df['cpu'] = np.append(np.repeat(range(cpus)[:n_remain],n),np.repeat(range(cpus)[n_remain:],n-1)) + 1

	##### GENERATE MODEL OBJECTS AND INITIALIZE OUTPUT FILES #####
	return_values = {}
	for k in cfg['model_order']:
		print "initializing out file for model " + k
		out_model = cfg['out'] if k == '___no_tag___' else cfg['out'] + '.' + k
		try:
			bgzfile = bgzf.BgzfWriter(out_model + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + out_model + '.gz').msg
			return 1

		##### SPLIT REGIONS #####
		if cpus > 1:
			pool = mp.Pool(cpus)
			for r in xrange(cpus):
				return_values[r+1] = pool.apply_async(process_regions, args=(regions,cfg,k,out_model + '.cpu' + str(r+1),r+1,True,))
				print "submitted cpu " + str(r+1) + " of " + str(cpus)
			pool.close()
			pool.join()

			##### CHECK RETURN VALUES #####
			print [return_values[i].get() for i in return_values]
			if 1 in [return_values[i].get() for i in return_values]:
				print Error("error detected").msg
				return 1

			##### COLLATE LOG FILES #####
			print "collating logs"
			for r in xrange(cpus):
				try:
					logfile = open(out_model + '.cpu' + str(r+1) + '.log', 'r')
				except:
					print Error("failed to initialize log file " + out_model + '.cpu' + str(r+1) + '.log').msg
					return 1
				print logfile.read()
				logfile.close()
				os.remove(out_model + '.cpu' + str(r+1) + '.log')

			##### COLLATE RESULTS #####
			print "collating results"
			i = 0
			write_header = True
			written = False
			for r in xrange(cpus):
				regions_df = regions.df[regions.df['cpu'] == r+1].reset_index(drop=True)
				i += 1
				for s in xrange(len(regions_df.index)):
					out_model_range = out_model + '.cpu' + str(r+1) + '.chr' + str(regions_df['chr'][s]) + 'bp' + str(regions_df['start'][s]) + '-' + str(regions_df['end'][s]) + '.h5'
					if return_values[r+1].get() != -1:
						store = pd.HDFStore(out_model_range)
						if not written:
							bgzfile.write(store.get_storer('df').attrs.results_header_metadata)
							bgzfile.write('\t'.join(store.get_storer('df').attrs.results_header))
							tbx_start = store.get_storer('df').attrs.tbx_start
							tbx_end = store.get_storer('df').attrs.tbx_end
							written = True
						if store['df'].shape[0] > 0:
							store['df'].to_csv(bgzfile, index=False, sep='\t', header=write_header, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.results_header, append=True)
							write_header = False
						store.close()
					os.remove(out_model_range)

		else:
			return_values[1] = process_regions(regions,cfg,k,out_model,1,False)
			if return_values[1] != -1:
				for s in xrange(len(regions.df.index)):
					out_model_range = out_model + '.chr' + str(regions.df['chr'][s]) + 'bp' + str(regions.df['start'][s]) + '-' + str(regions.df['end'][s]) + '.h5'
					store = pd.HDFStore(out_model_range)
					bgzfile.write(store.get_storer('df').attrs.results_header_metadata)
					bgzfile.write('\t'.join(store.get_storer('df').attrs.results_header))
					tbx_start = store.get_storer('df').attrs.tbx_start
					tbx_end = store.get_storer('df').attrs.tbx_end
					if store['df'].shape[0] > 0:
						store['df'].to_csv(bgzfile, index=False, sep='\t', header=True, na_rep='NA', float_format='%.5g', columns = store.get_storer('df').attrs.results_header, append=True)
					store.close()
					os.remove(out_model_range)

		bgzfile.close()
		print "indexing out file for model " + k
		try:
			pysam.tabix_index(out_model + '.gz',seq_col=0,start_col=tbx_start,end_col=tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + out_model + '.gz').msg
			return 1

	print "process complete"
	return 0
