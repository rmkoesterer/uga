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
import Parse
import pysam
from Bio import bgzf
import scipy.stats as scipy
import math
import Process
import logging
import subprocess
import os
import sys
import xlsxwriter

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunAnnot")

def RunAnnot(args):
	cfg = Parse.generate_annot_cfg(args)
	Parse.print_annot_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	results = pd.read_table(cfg['file'])

	outdf = results[['#chr','pos','id','a1','a2']]
	outdf.rename(columns={'#chr':'#CHROM','pos':'POS','a1':'REF','a2':'ALT','id':'ID'},inplace=True)
	outdf['QUAL'] = None
	outdf['FILTER'] = None
	outdf['INFO'] = None
	outdf.sort(["#CHROM","POS"],inplace=True)
	outdf.to_csv(f + '.annot1',header=True, index=False, sep='\t')

	try:
		cmd = 'java -jar ' + home_dir + '/snpEff.jar -v -canon GRCh37.75 ' + f + '.annot1 > ' + f + '.annot2'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "canonical annotation process terminated by user"
		sys.exit(1)

	try:
		cmd = 'java -jar ' + home_dir + '/SnpSift.jar dbnsfp -db ' + home_dir + '/data/dbNSFP_current -v ' + f + '.annot2 | sed "s/\+//g" > ' + f + '.annot3'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "dbNSFP annotation process terminated by user"
		sys.exit(1)

	try:
		cmd = 'java -jar ' + home_dir + '/SnpSift.jar extractFields -s "," -e "NA" ' + f + '.annot3 CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" "dbNSFP_GERP_RS" "dbNSFP_1000Gp1_ASN_AF" "dbNSFP_Uniprot_acc" "dbNSFP_MutationTaster_pred" "dbNSFP_ESP6500_AA_AF" "dbNSFP_SIFT_pred" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_1000Gp1_EUR_AF" "dbNSFP_1000Gp1_AFR_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_LRT_pred" "dbNSFP_1000Gp1_AMR_AF" "dbNSFP_GERP_NR" "dbNSFP_1000Gp1_AF" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Interpro_domain" | sed "s/ANN\[\*\]/ANN/g" > ' + f + '.annot'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "dbNSFP annotation process terminated by user"
		sys.exit(1)

	os.remove(f + '.annot1')
	os.remove(f + '.annot2')
	os.remove(f + '.annot3')

	if gc is not None:
		for col in gc.keys():
			print 'applying genomic inflation correction ' + str(gc[col]) + ' to column ' + col
			if col.replace('.p','.stderr') in results.columns:
				results[col.replace('.p','.stderr')] = results[col.replace('.p','.stderr')] * math.sqrt(float(gc[col]))
			if col.replace('.p','.z') in results.columns:
				results[col.replace('.p','.z')] = results[col.replace('.p','.z')] / math.sqrt(float(gc[col]))
			if col.replace('.p','.effect') in results.columns and col.replace('.p','.stderr') in results.columns:
				results[col] = 2 * scipy.norm.cdf(-1 * np.abs(results[col.replace('.p','.effect')]) / results[col.replace('.p','.stderr')])
			else:
				results[col] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*results[col]) / math.sqrt(float(gc[col]))))

	results.rename(columns={'#chr':'#CHROM','pos':'POS','a1':'REF','a2':'ALT','id':'ID'},inplace=True)
	annot = pd.read_table(f + '.annot')
	out = results.merge(annot,how='outer')

	out.fillna('NA',inplace=True)

	wkbk = xlsxwriter.Workbook(f + '.annot.xlsx')
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

	os.remove(f + '.annot')

	print "process complete"
	return 0
