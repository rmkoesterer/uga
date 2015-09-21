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

import Model
import IO
import Parse
import pysam
import Fxns
from Bio import bgzf
from Process import Error

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

	##### GENERATE MODEL OBJECT/S #####
	for k in cfg['model_order']:
		out = cfg['out'] if k == '___no_tag___' else cfg['out'] + '.' + k
		print "\nloading model " + k
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
		print "initializing out file"
		try:
			bgzfile = bgzf.BgzfWriter(out + '.gz', 'wb')
		except:
			print Error("failed to initialize bgzip format out file " + out + '.gz').msg
			return 1
	
		##### RUN ANALYSIS #####
		print "running models ..."
		written = False
		write_header=True
		analyzed = 0
		bgzfile.write(model.results_header_metadata)
		for r in range(len(regions.df.index)):
			i = 0
			print "loading region " + str(r+1) + '/' + str(len(regions.df.index)) + ' (' + regions.df['region'][r] + ')',
			try:
				model.get_region(regions.df['region'][r], regions.df['id'][r])
			except:
				print " <-- no variants found"
				pass
			else:
				print '...'
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

					model.out.to_csv(bgzfile, index=False, sep='\t', header=write_header, na_rep='NA', float_format='%.5g', columns = model.results_header)
					write_header = False
					analyzed += len(model.marker_stats['filter'][model.marker_stats['filter'] == 0])
	
					cur_markers = str(min(i*cfg['buffer'],(i-1)*cfg['buffer'] + model.biodata.marker_info.shape[0]))
					status = '   processed ' + cur_markers + ' variants, ' + str(analyzed) + ' passed filters'
					print status
	
		bgzfile.close()
		print "indexing out file"
		try:
			pysam.tabix_index(out + '.gz',seq_col=0,start_col=model.tbx_start,end_col=model.tbx_end,force=True)
		except:
			print Error('failed to generate index for file ' + out + '.gz').msg
			return 1

	print "process complete"
	return 0
