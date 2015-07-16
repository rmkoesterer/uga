#!/usr/bin/env python
from __main__ import *
import xlsxwriter

def Annot(cfg):
	Parse.PrintExploreOptions(cfg)

	f = cfg['file'].replace('.gz','') if cfg['file'].endswith('.gz') else cfg['file']

	print 'importing results'
	if cfg['file'][-3:len(cfg['file'])] == '.gz':
		results = pd.read_table(cfg['file'],compression='gzip')
	else:
		results = pd.read_table(cfg['file'])
	outdf = results[['#chr','pos','marker','a1','a2']]
	outdf.rename(columns={'#chr':'#CHROM','pos':'POS','ID':'marker','a1':'REF','a2':'ALT'},inplace=True)
	outdf['QUAL'] = None
	outdf['FILTER'] = None
	outdf['INFO'] = None
	outdf.sort(["#CHROM","POS"],inplace=True)
	outdf.to_csv(f + '.annot1',header=True, index=False, sep='\t')

	print 'generating canonical annotation'
	try:
		cmd = 'java -jar ' + cfg['snpeff'] + ' -canon GRCh37.75 ' + f + '.annot1  -stats ' + f + '.annot.summary.html -nolog > ' + f + '.annot2'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "canonical annotation process terminated by user"
		sys.exit(1)

	print 'generating dbNSFP annotation'
	try:
		cmd = 'java -jar ' + cfg['snpsift'] + ' dbnsfp -db ' + cfg['dbnsfp'] + ' ' + f + '.annot2 -nolog | sed "s/\+//g" > ' + f + '.annot3'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "dbNSFP annotation process terminated by user"
		sys.exit(1)

	print 'parsing annotation'
	try:
		cmd = 'java -jar ' + cfg['snpsift'] + ' extractFields -s "," -e "NA" ' + f + '.annot3 CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" "dbNSFP_GERP_RS" "dbNSFP_1000Gp1_ASN_AF" "dbNSFP_Uniprot_acc" "dbNSFP_MutationTaster_pred" "dbNSFP_ESP6500_AA_AF" "dbNSFP_SIFT_pred" "dbNSFP_phastCons100way_vertebrate" "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_1000Gp1_EUR_AF" "dbNSFP_1000Gp1_AFR_AF" "dbNSFP_ESP6500_EA_AF" "dbNSFP_LRT_pred" "dbNSFP_1000Gp1_AMR_AF" "dbNSFP_GERP_NR" "dbNSFP_1000Gp1_AF" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_Interpro_domain" | sed "s/ANN\[\*\]/ANN/g" > ' + f + '.annot'
		p = subprocess.Popen(cmd,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "SnpSift extraction terminated by user"
		sys.exit(1)

	os.remove(f + '.annot1')
	os.remove(f + '.annot2')
	os.remove(f + '.annot3')

	results.rename(columns={'#chr':'#CHROM','pos':'POS','marker':'ID','a1':'REF','a2':'ALT'},inplace=True)
	annot = pd.read_table(f + '.annot')
	out = results.merge(annot,how='outer')
	out.rename(columns={'#CHROM':'#chr','POS':'pos','ID':'marker','REF':'a1','ALT':'a2'},inplace=True)

	out.fillna('NA',inplace=True)

	print 'generating xlsx file'
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
			if field in ['#chr','pos'] or field.endswith('.filtered') or field.endswith('.n'):
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