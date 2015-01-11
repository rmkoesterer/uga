import os
import sys
import subprocess
import tabix
import math
import gzip
import time
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
import re
from MarkerCalc import Complement
from Coordinates import Coordinates
from Messages import Error
from multi_key_dict import multi_key_dict
from Bio import bgzf
import scipy.stats as scipy

def Meta(cfg=None, 
			region=None, 
			region_list=None, 
			method=None, 
			mem=3):

	print "   ... arguments"
	for arg in locals().keys():
		if not locals()[arg] in [None, False]:
			print "      {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	assert cfg, Error("no configuration file specified")
	assert method, Error("no meta-analysis method specified")

	##### get headers from process files #####
	for tag in cfg['data_info'].keys():
		try:
			p = subprocess.Popen(['tabix','-h',cfg['data_info'][tag]['process_file'],'0'], stdout=subprocess.PIPE)
		except:
			usage(Error("process_file " + cfg['data_info'][tag]['process_file'] + " has incorrect format"))
		cfg['data_info'][tag]['header'] = p.communicate()[0]
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].replace("#","")
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].strip()
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].split()		

	##### initialize out file #####
	bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')
	
	output = {}
	out_cols = []
	for tag in cfg['file_order']:
		for x in ['freq','effect','stderr','or','z','p']:
			if x in cfg['data_info'][tag]:
				out_cols.append(tag + '.' + x)
	out_cols.append(tag + '.filtered')
	
	##### determine markers to be analyzed #####
	print "   ... generating region list"
	if region_list:
		marker_list = Coordinates(region_list).Load()
	if region:
		marker_list = pd.DataFrame({'chr': [re.split(':|-',region)[0]],'start': [re.split(':|-',region)[1]],'end': [re.split(':|-',region)[2]],'region': [region]})
	marker_list['n'] = 0
	for i in range(len(marker_list.index)):
		for key in cfg['data_info'].keys():
			tb = tabix.open(cfg['data_info'][key]['process_file'])
			try:
				records = tb.querys(marker_list['region'][i])
			except:
				pass
			else:
				marker_list['n'][i] = marker_list['n'][i] + 1
				break
	if marker_list['n'].sum() == 0:
		print Error("no markers found")
		return()

	cfg['sig'] = int(cfg['sig']) if 'sig' in cfg.keys() else 5

	for r in range(len(marker_list.index)):
		reg = marker_list['region'][r]
		print "   ... building reference database for " + str(r + 1) + " of " + str(len(marker_list.index)) + " regions"
		delim = "><"
		ref = multi_key_dict()
		ref_alleles = multi_key_dict()
		records_count = {}
		for tag in cfg['file_order']:
			chr = cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['chr'])
			pos = cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['pos'])
			marker = cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['marker'])
			a1 = cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['a1'])
			a2 = cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['a2'])
			n = None if not 'n' in cfg['data_info'][tag] else int(cfg['data_info'][tag]['n'])
			effect = None if not 'effect' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['effect'])
			stderr = None if not 'stderr' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['stderr'])
			freq = None if not 'freq' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['freq'])
			rsq = None if not 'rsq' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['rsq'])
			o_r = None if not 'or' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['or'])
			z = None if not 'z' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['z'])
			p = None if not 'p' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['p'])
			tb = tabix.open(cfg['data_info'][tag]['process_file'])
			records_count[tag] = 0
			try:
				records = tb.querys(reg)
			except:
				pass
			else:
				for record in records:
					filtered = 0
					records_count[tag] = records_count[tag] + 1
					if ref.has_key(record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]):
						if not record[marker] in ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]]:
							ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]].append(record[marker])
							ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]] = list(set(ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]]))
					else:
						ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2], record[chr] + delim + record[pos] + delim + Complement(record[a1]) + delim + Complement(record[a2]),record[chr] + delim + record[pos] + delim + Complement(record[a2]) + delim + Complement(record[a1]),record[chr] + delim + record[pos] + delim + record[a2] + delim + record[a1]] = [record[marker]]
						ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2], record[chr] + delim + record[pos] + delim + Complement(record[a1]) + delim + Complement(record[a2]),record[chr] + delim + record[pos] + delim + Complement(record[a2]) + delim + Complement(record[a1]),record[chr] + delim + record[pos] + delim + record[a2] + delim + record[a1]] = [record[a1],record[a2]]
					refmarker = filter(lambda x:'rs' in x, ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]])
					if len(refmarker) > 0:
						refmarker = refmarker[0]
					else:
						refmarker = ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][0]
					refa1 = ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][0]
					refa2 = ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][1]
					if (refa1 + refa2 == Complement(record[a2]) + Complement(record[a1]) or refa1 + refa2 == record[a2] + record[a1]) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
						if effect:
							record[effect] = str(-1 * float(record[effect])) if record[effect] != "NA" else float('nan')
						if freq:
							record[freq] = str(1-float(record[freq])) if record[freq] != "NA" else float('nan')
						if o_r:
							record[o_r] = str(1/float(record[o_r])) if record[o_r] != "NA" else float('nan')
						if z:
							record[z] = str(-1 * float(record[z])) if record[z] != "NA" else float('nan')
					for f in cfg['data_info'][tag]['filters']:
						col = cfg['data_info'][tag]['header'].index(f.split()[0])
						if f.split()[0] != "status":
							f = f.replace(f.split()[0],str(record[col]))
						else:
							f = f.replace(f.split()[0],'"' + str(record[col]) + '"')
						if eval('"' + str(record[col]) + '" == \'nan\' or "' + str(record[col]) + '" == \'NA\' or not ' + f):
							filtered = 1
					out_vals = {tag + '.filtered': filtered}
					if n:
						out_vals[tag + '.n'] = n
					else:
						out_vals[tag + '.n'] = float('nan')
					if freq:
						out_vals[tag + '.freq'] = float(record[freq]) if record[freq] != "NA" else float('nan')
					if rsq:
						out_vals[tag + '.rsq'] = float(record[rsq]) if record[rsq] != "NA" else float('nan')
					if effect:
						out_vals[tag + '.effect'] = float(record[effect]) if record[effect] != "NA" else float('nan')
					if stderr:
						out_vals[tag + '.stderr'] = float(record[stderr]) if record[stderr] != "NA" else float('nan')
					if o_r:
						out_vals[tag + '.or'] = float(record[o_r]) if record[o_r] != "NA" else float('nan')
					if z:
						out_vals[tag + '.z'] = float(record[z]) if record[z] != "NA" else float('nan')
					if p:
						out_vals[tag + '.p'] = float(record[p]) if record[p] != "NA" else float('nan')
					if not record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2 in output.keys():
						output[record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2] = dict(out_vals)
					else:
						output[record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2].update(dict(out_vals))
				print "          merged cohort " + tag + " : " + str(len(output.keys())) + " markers"
	output_df = pd.DataFrame(output).transpose()
	header=['chr','pos','a1','a2','marker']
	for tag in cfg['file_order']:
		n = None if not 'n' in cfg['data_info'][tag] else int(cfg['data_info'][tag]['n'])
		effect = None if not 'effect' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['effect'])
		stderr = None if not 'stderr' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['stderr'])
		freq = None if not 'freq' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['freq'])
		rsq = None if not 'rsq' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['rsq'])
		o_r = None if not 'or' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['or'])
		z = None if not 'z' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['z'])
		p = None if not 'p' in cfg['data_info'][tag] else cfg['data_info'][tag]['header'].index(cfg['data_info'][tag]['p'])
		for suffix in ['freq','effect','stderr','or','z','p','n','filtered']:
			if not tag + '.' + suffix in output_df.columns.values:
				output_df[tag + '.' + suffix] = float('nan')
			header.append(tag + '.' + suffix)

	output_df['ref'] = output_df.index
	output_df['chr'] = output_df['ref'].apply(lambda x: x.split(delim)[0]).astype(int)
	output_df['pos'] = output_df['ref'].apply(lambda x: x.split(delim)[1]).astype(int)
	output_df['marker'] = output_df['ref'].apply(lambda x: x.split(delim)[2])
	output_df['a1'] = output_df['ref'].apply(lambda x: x.split(delim)[3])
	output_df['a2'] = output_df['ref'].apply(lambda x: x.split(delim)[4])
	
	for tag in cfg['file_order']:
		for suffix in ['freq','effect','stderr','or','z','p','n','filtered']:
			if tag + '.' + suffix in output_df.columns.values:
				output_df[tag + '.' + suffix] = output_df[tag + '.' + suffix].convert_objects(convert_numeric=True)

	##### apply genomic control #####
	if method in ['sample_size','stderr','efftest']:
		for tag in cfg['file_order']:
			if 'gc' in cfg['data_info'][tag].keys():
				print "   ... applying genomic control correction for cohort " + tag
				if 'effect' in cfg['data_info'][tag].keys():
					output_df[tag + '.effect'] = output_df[tag + '.effect'] / math.sqrt(float(cfg['data_info'][tag]['gc']))
				if 'stderr' in cfg['data_info'][tag].keys():
					output_df[tag + '.stderr'] = output_df[tag + '.stderr'] / math.sqrt(float(cfg['data_info'][tag]['gc']))
				if 'or' in cfg['data_info'][tag].keys():
					output_df[tag + '.or'] = np.exp(np.log(output_df[tag + '.or']) / math.sqrt(float(cfg['data_info'][tag]['gc'])))
				if 'z' in cfg['data_info'][tag].keys():
					output_df[tag + '.z'] = output_df[tag + '.z'] / math.sqrt(float(cfg['data_info'][tag]['gc']))
				if 'p' in cfg['data_info'][tag].keys():
					output_df[tag + '.p'] = 2 * scipy.norm.cdf(-1 * np.abs(output_df[tag + '.effect']) / output_df[tag + '.stderr'])
	output_df['chr'] = output_df['chr'].astype(int)
	output_df['pos'] = output_df['pos'].astype(int)
	
	##### apply efftest correction #####
	if method == 'efftest':
		for tag in cfg['file_order']:
			print "   ... applying efftest correction for " + tag
			with gzip.open(cfg['data_info'][tag]['efftest_file']) as f:
				efftests = pd.read_table(f, header=False)
			output_df[tag + '.n_eff'] = 0
			output_df[tag + '.p_eff'] = 0
			header.extend([tag + '.n_eff', tag + '.p_eff'])
			for r in range(len(marker_list.index)):
				if marker_list['reg_id'][r] in list(efftests['reg_id']):
					output_df[tag + '.n_eff'][(output_df['chr'] == marker_list['chr'][r]) & (output_df['pos'] >= marker_list['start'][r]) & (output_df['chr'] <= marker_list['end'][r])] = efftests['n_eff'][efftests['reg_id'] == marker_list['reg_id'][r]].values[0]
			output_df[tag + '.p_eff'] = output_df[tag + '.p'] * output_df[tag + '.n_eff']
			output_df[tag + '.p_eff'] = output_df.apply(lambda x: 0.9999999999 if x[tag + '.p_eff'] >= 1 else x[tag + '.p_eff'],axis=1)

	##### meta analysis #####
	for meta in cfg['meta_order']:
		print "   ... running meta analysis for cohorts " + meta
		if method in ['efftest','sample_size']:
			p_ext = '.p' if method == 'sample_size' else '.p_eff'
			output_df[meta + '.dir'] = ''
			for tag in cfg['meta_info'][meta]:
				output_df[tag + '.filtered'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + p_ext]) or x[tag + p_ext] > 1 or x[tag + p_ext] <= 0 else x[tag + '.filtered'],axis=1)
				filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filtered')][0]
				N_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith(p_ext)][0]
				Eff_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				output_df[meta + '.' + tag + '.dir'] = output_df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				output_df[meta + '.' + tag + '.zi'] = output_df.apply(lambda x: (-1 * scipy.norm.ppf(1 - (x[P_idx]/2)) if x[Eff_idx] < 0 else scipy.norm.ppf(1 - (x[P_idx]/2))) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.n'] = output_df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.wi'] = output_df.apply(lambda x: math.sqrt(x[meta + '.' + tag + '.n']) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.ziwi'] = output_df.apply(lambda x: x[meta + '.' + tag + '.zi'] * x[meta + '.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.dir'] = output_df[meta + '.dir'] + output_df[meta + '.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.wi')]
			ZiWi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.ziwi')]
			output_df[meta + '.n'] = output_df.apply(lambda x: x[N_idx_all].sum(),axis=1)
			output_df[meta + '.z'] = output_df.apply(lambda x: x[ZiWi_idx_all].sum()/math.sqrt(x[N_idx_all].sum()) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
		else:
			output_df[meta + '.dir'] = ''
			for tag in cfg['meta_info'][meta]:
				output_df[tag + '.filtered'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + '.p']) or x[tag + '.p'] > 1 or x[tag + '.p'] <= 0 else x[tag + '.filtered'],axis=1)
				filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filtered')][0]
				N_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.p')][0]
				Eff_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				StdErr_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.stderr')][0]
				output_df[meta + '.' + tag + '.dir'] = output_df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				output_df[meta + '.' + tag + '.n'] = output_df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.wi'] = output_df.apply(lambda x: 1/(x[StdErr_idx]**2) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.biwi'] = output_df.apply(lambda x: x[Eff_idx] * x[meta + '.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.dir'] = output_df[meta + '.dir'] + output_df[meta + '.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.wi')]
			BiWi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.biwi')]
			output_df[meta + '.n'] = output_df.apply(lambda x: x[N_idx_all].sum(),axis=1)
			output_df[meta + '.stderr'] = output_df.apply(lambda x: math.sqrt(1/(x[Wi_idx_all].sum())) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			output_df[meta + '.effect'] = output_df.apply(lambda x: (x[BiWi_idx_all].sum())/(x[Wi_idx_all].sum()) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			output_df[meta + '.z'] = output_df.apply(lambda x: x[meta + '.effect']/x[meta + '.stderr'] if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			header.append(meta + '.effect')
			header.append(meta + '.stderr')
		output_df[meta + '.p'] = output_df.apply(lambda x: 2 * scipy.norm.cdf(-1 * abs(float(x[meta + '.z']))) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
		output_df[meta + '.dir'] = output_df.apply(lambda x: x[meta + '.dir'] if not math.isnan(x[meta + '.p']) else float('nan'), axis=1)
		output_df[meta + '.effect'] = output_df[meta + '.effect'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.stderr'] = output_df[meta + '.stderr'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.z'] = output_df[meta + '.z'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.p'] = output_df[meta + '.p'].map(lambda x: '%.2e' % x if not math.isnan(x) else x)
		header.append(meta + '.z')
		header.append(meta + '.p')
		header.append(meta + '.dir')
		header.append(meta + '.n')
		
	##### format output to significant digits #####
	for tag in cfg['file_order']:
		if 'effect' in cfg['data_info'][tag].keys():
			output_df[tag + '.effect'] = output_df[tag + '.effect'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'stderr' in cfg['data_info'][tag].keys():
			output_df[tag + '.stderr'] = output_df[tag + '.stderr'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'or' in cfg['data_info'][tag].keys():
			output_df[tag + '.or'] = output_df[tag + '.or'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'z' in cfg['data_info'][tag].keys():
			output_df[tag + '.z'] = output_df[tag + '.z'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'p' in cfg['data_info'][tag].keys():
			output_df[tag + '.p'] = output_df[tag + '.p'].map(lambda x: '%.2e' % x if not math.isnan(x) else x)
		if 'n' in cfg['data_info'][tag].keys():
			output_df[tag + '.n'] = output_df[tag + '.n'].map(lambda x: '%.0f' % x if not math.isnan(x) else x)
		output_df[tag + '.filtered'] = output_df[tag + '.filtered'].map(lambda x: '%.0f' % x if not math.isnan(x) else x)
	if meta + '.n' in list(output_df.columns.values):
		output_df[meta + '.n'] = output_df[meta + '.n'].map(lambda x: '%.0f' % x if not math.isnan(x) else x)

	##### extract efftest corrected top p values #####
	if method == 'efftest':
		out_efftest_df = marker_list.copy()
		header = ['chr','start','end','reg_id']
		for tag in cfg['file_order']:
			out_efftest_df[tag + '.min_snp'] = float('nan')
			out_efftest_df[tag + '.min_pos'] = float('nan')
			out_efftest_df[tag + '.min_p'] = float('nan')
			for r in range(len(out_efftest_df.index)):
				min_idx = output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)][tag + '.p_eff'].argmin() if len(output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)].index) > 0 else float('nan')
				out_efftest_df[tag + '.min_snp'][r] = output_df['marker'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[tag + '.min_pos'][r] = output_df['pos'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[tag + '.min_p'][r] = output_df[tag + '.p'][min_idx] if min_idx in output_df.index else float('nan')
			header.append(tag + '.min_snp')
			header.append(tag + '.min_pos')
			header.append(tag + '.min_p')
			out_efftest_df[tag + '.min_pos'] = out_efftest_df[tag + '.min_pos'].map(lambda x: '%.0f' % x if not math.isnan(x) else x)
			out_efftest_df[tag + '.min_p'] = out_efftest_df[tag + '.min_p'].map(lambda x: '%.2e' % x if not math.isnan(x) else x)
		for meta in cfg['meta_order']:
			out_efftest_df[meta + '.min_snp'] = float('nan')
			out_efftest_df[meta + '.min_pos'] = float('nan')
			out_efftest_df[meta + '.min_dir'] = float('nan')
			out_efftest_df[meta + '.min_p'] = float('nan')
			for r in range(len(out_efftest_df.index)):
				min_idx = output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)][meta + '.p'].argmin() if len(output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)].index) > 0 else float('nan')
				out_efftest_df[meta + '.min_snp'][r] = output_df['marker'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_pos'][r] = output_df['pos'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_dir'][r] = output_df[meta + '.dir'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_p'][r] = output_df[meta + '.p'][min_idx] if min_idx in output_df.index else float('nan')
			header.append(meta + '.min_snp')
			header.append(meta + '.min_pos')
			header.append(meta + '.min_dir')
			header.append(meta + '.min_p')
			out_efftest_df[meta + '.min_pos'] = out_efftest_df[meta + '.min_pos'].map(lambda x: '%.0f' % x if not math.isnan(x) else x)
			out_efftest_df[meta + '.min_p'] = out_efftest_df[meta + '.min_p'].map(lambda x: '%.2e' % x if not math.isnan(x) else x)
		out_efftest_df = out_efftest_df[header]
		out_efftest_df.fillna('NA',inplace=True)
		out_efftest_df.rename(columns=lambda x: x.replace('chr','#chr'),inplace=True)
		out_efftest_df.to_csv(bgzfile, header=True, index=False, sep="\t")
	else:
		output_df = output_df[header]	
		output_df.fillna('NA',inplace=True)
		output_df.sort(['chr','pos'],inplace=True)
		output_df.rename(columns=lambda x: x.replace('chr','#chr'),inplace=True)
		output_df.to_csv(bgzfile, header=True, index=False, sep="\t")
	bgzfile.close()

	print "   ... mapping results file"
	cmd = 'tabix -b 2 -e 3 ' + cfg['out'] + '.gz'
	p = subprocess.Popen(cmd, shell=True)
	time.sleep(1)
	p.wait()
	
	print '   ... process complete'
