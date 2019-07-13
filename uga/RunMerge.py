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
import Geno
import Parse
import Variant
import pysam
from Bio import bgzf
import Process
import multiprocessing as mp
import sys
import os
import logging
import pickle

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunMerge")

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
			results_obj[f] = Geno.Results(filename=cfg['files'][f], chr=cfg['columns'][f]['chr'], pos=cfg['columns'][f]['pos'], id=cfg['columns'][f]['id'], a1=cfg['columns'][f]['a1'], a2=cfg['columns'][f]['a2'])
		except Process.Error as err:
			print err.out
			return 1

	variants_found = False
	variant_ref = Variant.Ref()
	results_final = None
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

		if k == 0:
			results_final = results_region.copy()
		else:
			results_final = results_final.merge(results_region, how='outer')

	results_final = results_final[[a for a in results_final.columns if a not in ['id_unique','___uid___']]]
	results_final = results_final.sort_values(by=['chr','pos'])
	results_final['chr'] = results_final['chr'].astype(np.int64)
	results_final['pos'] = results_final['pos'].astype(np.int64)
	pkl = open('/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(cpu) + '.pkl', "wb")
	pickle.dump([results_final,np.array(results_final.columns.values)],pkl,protocol=2)
	pkl.close()

	if log:
		sys.stdout = stdout_orig
		log_file.close()

	if variants_found:
		return 0
	else:
		return -1

def RunMerge(args):
	cfg = Parse.generate_merge_cfg(args)
	Parse.print_merge_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	regions_df = pd.read_table(cfg['region_file'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
	regions_df = regions_df[regions_df['job'] == int(cfg['job'])].reset_index(drop=True)
	return_values = {}
	print ''
	try:
		bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')
	except:
		print Process.Error("failed to initialize bgzip format out file " + cfg['out'] + '.gz').out
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

	written = False
	for i in xrange(1,cfg['cpus']+1):
		out = '/'.join(cfg['out'].split('/')[0:-1]) + '/' + cfg['out'].split('/')[-1] + '.cpu' + str(i) + '.pkl'
		pkl = open(out,"rb")
		results_final,results_header = pickle.load(pkl)
		if not written:
			bgzfile.write('#' + '\t'.join(results_header) + '\n')
			written = True
		if results_final.shape[0] > 0:
			results_final.replace({'None': 'NA', 'nan': 'NA'}).to_csv(bgzfile, index=False, sep='\t', header=False, na_rep='NA', float_format='%.5g', columns = results_header, append=True)
		pkl.close()
		os.remove(out)

	bgzfile.close()
	print "indexing out file"
	try:
		pysam.tabix_index(cfg['out'] + '.gz',seq_col=0,start_col=1,end_col=1,force=True)
	except:
		print Process.Error('failed to generate index for file ' + cfg['out'] + '.gz').out
		return 1

	if cfg['snpeff']:
		from ConfigParser import SafeConfigParser
		from pkg_resources import resource_filename
		import subprocess
		import xlsxwriter
		import time
		ini = SafeConfigParser()
		ini.read(resource_filename('uga', 'settings.ini'))

		results_final = pd.read_table(cfg['out'] + '.gz')
		outdf = results_final[['#chr','pos','id','a1','a2']]
		outdf = outdf.rename(columns={'#chr':'#CHROM','pos':'POS','id':'ID','a1':'REF','a2':'ALT'})
		outdf['QUAL'] = None
		outdf['FILTER'] = None
		outdf['INFO'] = None
		outdf.to_csv(cfg['out'] + '.annot1',header=True, index=False, sep='\t')

		time.sleep(1)
		try:
			cmd = 'java -jar ' + ini.get('main','snpeff') + ' -s ' + cfg['out'] + '.annot.summary.html -v -canon GRCh37.75 ' + cfg['out'] + '.annot1 > ' + cfg['out'] + '.annot2'
			print cmd
			p = subprocess.Popen(cmd,shell=True)
			p.wait()
		except KeyboardInterrupt:
			kill_all(p.pid)
			print "canonical annotation process terminated by user"
			sys.exit(1)

		return
		time.sleep(1)
		try:
			cmd = 'java -jar ' + ini.get('main','snpsift') + ' extractFields -s "," -e "NA" ' + cfg['out'] + '.annot2 CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" | sed "s/ANN\[\*\]/ANN/g" > ' + cfg['out'] + '.annot'
			print cmd
			p = subprocess.Popen(cmd,shell=True)
			p.wait()
		except KeyboardInterrupt:
			kill_all(p.pid)
			print "SnpSift annotation process terminated by user"
			sys.exit(1)
		os.remove(cfg['out'] + '.annot1')
		os.remove(cfg['out'] + '.annot2')

		results_final = results_final.rename(columns={'#chr':'#CHROM','pos':'POS','id':'ID','a1':'REF','a2':'ALT'})
		annot = pd.read_table(cfg['out'] + '.annot')
		out = results_final.merge(annot,how='outer')

		out.fillna('NA',inplace=True)

		wkbk = xlsxwriter.Workbook(cfg['out'] + '.annot.xlsx')
		wksht = wkbk.add_worksheet()

		header_format = wkbk.add_format({'bold': True,
											 'align': 'center',
											 'valign': 'vcenter'})
		string_format = wkbk.add_format({'align': 'center', 'valign': 'center'})
		float_format = wkbk.add_format({'align': 'center', 'valign': 'center'})
		float_format.set_num_format('0.000')
		integer_format = wkbk.add_format({'align': 'center', 'valign': 'center'})
		integer_format.set_num_format('0')
		sci_format = wkbk.add_format({'align': 'center', 'valign': 'center'})
		sci_format.set_num_format('0.00E+00')
		i = 0
		for field in out.columns:
			wksht.write(0,i,field,header_format)
			i += 1

		i = 0
		for row in range(out.shape[0]):
			j = 0
			for field in out.columns:
				if field in ['#CHROM','POS'] or field.endswith('.filtered') or field.endswith('.n'):
					wksht.write(row+1,j,out[field][i], integer_format)
				elif field.endswith(('.p','hwe','hwe.unrel')):
					wksht.write(row+1,j,out[field][i], sci_format)
				elif field.endswith(('.effect','.stderr','.or','.z','freq','freq.unrel','rsq','rsq.unrel','callrate')):
					wksht.write(row+1,j,out[field][i], float_format)
				else:
					wksht.write(row+1,j,out[field][i], string_format)
				j += 1
			i += 1
		wksht.freeze_panes(1, 0)
		wkbk.close()

		os.remove(cfg['out'] + '.annot')

	print "process complete"
	return 0
