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

import os
import argparse
from __version__ import version
import textwrap
from datetime import date
from Args import *

def Version():
	print ''
	print 'uga v' + version + ' (c) 2015 Ryan Koesterer   GNU General Public License v3'
	print ''
	print textwrap.fill("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.", initial_indent='   ', subsequent_indent='   ')
	print ''
	print textwrap.fill("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.", initial_indent='   ', subsequent_indent='   ')
	print ''
	print textwrap.fill("You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>", initial_indent='   ', subsequent_indent='   ')
	print ''

def GetParser():
	main_parser = argparse.ArgumentParser()
	main_parser.add_argument('--version', 
						action='version', 
						version=Version(), 
						help='display version information and exit')
	subparsers = main_parser.add_subparsers(title='modules', dest='which')
	global_parser = argparse.ArgumentParser(add_help=False)

	set_parser = SetArgs(subparsers.add_parser('set', help='user definable settings', parents=[global_parser]))
	snv_parser = SnvArgs(subparsers.add_parser('snv', help='run single nucleotide variant association models', parents=[global_parser]))
	gene_parser = GeneArgs(subparsers.add_parser('gene', help='run gene / locus based association models', parents=[global_parser]))

	meta_parser = MetaArgs(subparsers.add_parser('meta', help='meta-analysis', parents=[global_parser]))
	map_parser = MapArgs(subparsers.add_parser('map', help='map non-empty regions in genotype/imputed data files', parents=[global_parser]))
	compile_parser = CompileArgs(subparsers.add_parser('compile', help='verify and compile results files', parents=[global_parser]))
	eval_parser = EvalArgs(subparsers.add_parser('eval', help='evaluate results: filter, plot, list top results, etc.', parents=[global_parser]))
	gc_parser = GcArgs(subparsers.add_parser('gc', help='apply genomic control to 1 or more p-value columns', parents=[global_parser]))
	annot_parser = AnnotArgs(subparsers.add_parser('annot', help='annotate variant results using snpEff and SnpSift', parents=[global_parser]))

	return main_parser
	
def GetArgs(parser):
	args=parser.parse_args()
	if 'ordered_args' in args:
		for k in args.ordered_args:
			vars(args).update({k[0]: k[1]})
	else:
		if args.which in ['snv','gene','meta','map','compile','eval','gc','annot']:
			parser.error("missing argument: no options selected")
	print ''
	print 'active module: ' + args.which
	return args

def GenerateSnvCfg(args):
	config = {'out': None, 'buffer': 1000, 'region': None, 'region_file': None, 'id': None, 'cpus': 1, 'mb': 1, 'qsub': None, 'split': False, 'split_n': None, 'replace': False, 
					'models': {}, 'model_order': [], 'meta': {}, 'meta_order': []}
	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'buffer' and arg[1] is not None:
			config['buffer'] = arg[1]
		if arg[0] == 'region':
			config['region'] = arg[1]
		if arg[0] == 'region_file':
			config['region_file'] = arg[1]
		if arg[0] == 'id':
			config['id'] = arg[1]
		if arg[0] == 'meta':
			config['meta'][arg[1].split(':')[0]] = arg[1].split(':')[1]
		if arg[0] == 'cpus' and arg[1] is not None:
			config['cpus'] = arg[1]
		if arg[0] == 'mb' and arg[1] is not None:
			config['mb'] = arg[1]
		if arg[0] == 'qsub':
			config['qsub'] = arg[1]
		if arg[0] == 'split' and arg[1] is True:
			config['split'] = arg[1]
		if arg[0] == 'split_n':
			config['split_n'] = arg[1]

	tags_idx = [args.index((x,y)) for x, y in args if x == 'tag'] + [len(args)]
	global_args = [x for x in args[:tags_idx[0]] if x[0] not in config]
	config_default = {'formula': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
							'male': 1, 'female': 2, 'miss': 0.0, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
							'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'snv_list': None}
	if len(tags_idx) > 1:
		for i in xrange(len(tags_idx[:-1])):
			config['models'][args[tags_idx[i]][1]] = config_default
			for arg in global_args:
				if arg[0] in ['bssmeta']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['bssmeta']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in global_args:
				if arg[0] in ['bssmeta']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['formula'] = arg[1]
				elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['bssmeta']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['formula'] = arg[1]
				elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			config['model_order'].append(args[tags_idx[i]][1])
	else:
		config['models']['___no_tag___'] = config_default
		for arg in global_args:
			if arg[0] in ['bssmeta']:
				config['models']['___no_tag___']['case_code'] = 1
				config['models']['___no_tag___']['ctrl_code'] = 0
		for arg in global_args:
			if arg[0] in ['bssmeta']:
				config['models']['___no_tag___']['fxn'] = arg[0]
				config['models']['___no_tag___']['formula'] = arg[1]
			elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
				config['models']['___no_tag___']['format'] = arg[0]
				config['models']['___no_tag___']['file'] = arg[1]
			else:
				config['models']['___no_tag___'][arg[0]] = arg[1]
		config['model_order'].append('___no_tag___')
	for arg in args:
		if arg[0] == 'meta':
			config['meta_order'].append(arg[1].split(':')[0])
			config['meta'][arg[1].split(':')[0]] = arg[1].split(':')[1]
	return config

def PrintSnvOptions(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if not k in ['models','model_order','meta','meta_order']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['models']:
		print '   model ' + str(m) + ' ...' if len(cfg['models']) > 1 else '   model ...'
		for n in cfg['models'][m]:
			if cfg['models'][m][n] is not None and cfg['models'][m][n] is not False:
				if cfg['models'][m][n] is True:
					print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len))) + " " + str(cfg['models'][m][n])
	if len(cfg['meta_order']) > 0:
		print '   meta analysis ...'
		for m in cfg['meta_order']:
			print "      {0:>{1}}".format(str('--meta ' + m), len(max(['--' + k for k in cfg['meta_order']],key=len))) + ":" + str(cfg['meta'][m])

def GenerateGeneCfg(args):
	config = {'out': None, 'buffer': 1000, 'id': None, 'cpus': 1, 'qsub': None, 'split': False, 'split_n': None, 'replace': False, 'gene_map': None, 
					'models': {}, 'model_order': [], 'meta': {}, 'meta_order': []}
	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'buffer' and arg[1] is not None:
			config['buffer'] = arg[1]
		if arg[0] == 'meta':
			config['meta'][arg[1].split(':')[0]] = arg[1].split(':')[1]
		if arg[0] == 'cpus' and arg[1] is not None:
			config['cpus'] = arg[1]
		if arg[0] == 'qsub':
			config['qsub'] = arg[1]
		if arg[0] == 'split' and arg[1] is True:
			config['split'] = arg[1]
		if arg[0] == 'split_n':
			config['split_n'] = arg[1]
		if arg[0] == 'gene_map':
			config['gene_map'] = arg[1]

	tags_idx = [args.index((x,y)) for x, y in args if x == 'tag'] + [len(args)]
	global_args = [x for x in args[:tags_idx[0]] if x[0] not in config]
	config_default = {'formula': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
							'male': 1, 'female': 2, 'miss': 0.0, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
							'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None}
	if len(tags_idx) > 1:
		for i in xrange(len(tags_idx[:-1])):
			config['models'][args[tags_idx[i]][1]] = config_default
			for arg in global_args:
				if arg[0] in ['bskato']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['bskato']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in global_args:
				if arg[0] in ['bskato']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['formula'] = arg[1]
				elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['bskato']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['formula'] = arg[1]
				elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			config['model_order'].append(args[tags_idx[i]][1])
	else:
		config['models']['___no_tag___'] = config_default
		for arg in global_args:
			if arg[0] in ['bskato']:
				config['models']['___no_tag___']['case_code'] = 1
				config['models']['___no_tag___']['ctrl_code'] = 0
		for arg in global_args:
			if arg[0] in ['bskato']:
				config['models']['___no_tag___']['fxn'] = arg[0]
				config['models']['___no_tag___']['formula'] = arg[1]
			elif arg[0] in ['plink','vcf','dos1','dos2','oxford']:
				config['models']['___no_tag___']['format'] = arg[0]
				config['models']['___no_tag___']['file'] = arg[1]
			else:
				config['models']['___no_tag___'][arg[0]] = arg[1]
		config['model_order'].append('___no_tag___')
	for arg in args:
		if arg[0] == 'meta':
			config['meta_order'].append(arg[1].split(':')[0])
			config['meta'][arg[1].split(':')[0]] = arg[1].split(':')[1]
	return config

def PrintGeneOptions(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if not k in ['models','model_order','meta','meta_order']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['models']:
		print '   model ' + str(m) + ' ...' if len(cfg['models']) > 1 else '   model ...'
		for n in cfg['models'][m]:
			if cfg['models'][m][n] is not None and cfg['models'][m][n] is not False:
				if cfg['models'][m][n] is True:
					print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len))) + " " + str(cfg['models'][m][n])
	if len(cfg['meta_order']) > 0:
		print '   meta analysis ...'
		for m in cfg['meta_order']:
			print "      {0:>{1}}".format(str('--meta ' + m), len(max(['--' + k for k in cfg['meta_order']],key=len))) + ":" + str(cfg['meta'][m])







def GenerateMetaCfg(args):
	config = {'out': None, 'method': None, 'buffer': 100, 'region_file': None, 'region': None, 'id': None, 'write_header': False, 'write_eof': False,
					'data_info': {}, 'meta_info': {}, 'meta_order': [], 'file_order': []}

	##### add top level variables to config
	top_args = [a for a in args if a[0] in ['out','method','buffer','region_file','region','id']]
	for arg in top_args:
		config[arg[0]] = arg[1]

	# list all possible model level arguments
	model_vars = ['gc','marker_col','freq_col','rsq_col','hwe_col','effect_col','stderr_col','or_col','z_col','p_col','n_col','n','tag','file','maf','rsq','hwe']
	if not 'tag' in [a[0] for a in args]:
		args = [('tag', 'A')] + args
	pre_tag_idx = [a[0] for a in args if a[0] in model_vars].index('tag')

	# extract global model arguments
	stat_args_global = dict([a for a in args if a[0] in model_vars][:pre_tag_idx])
	for k in stat_args_global.keys():
		stat_args_global[k] = stat_args_global[k]

	# extract local model arguments
	stat_args_local = [a for a in args if a[0] in model_vars][pre_tag_idx:]
	local_tag_idx = [a for a, b in enumerate(stat_args_local) if b[0] == 'tag']
	if len(local_tag_idx) == 0:
		local_tag_idx = [0]
	local_args = [stat_args_local[i:j] for i, j in zip([0]+local_tag_idx[1:], local_tag_idx[1:]+[None])]
	# extract meta arguments
	meta_args = [a for a in args if a[0] == 'meta']
	for k in meta_args:
		config['meta_info'][k[1].split(':')[0]] = k[1].split(':')[1].split('+')
		config['meta_order'].append(k[1].split(':')[0])

	# loop over model argument sets and add to config['models']
	for l in local_args:
		tag = l[0][1]
		config['file_order'].append(tag)
		config['data_info'][tag] = {'tag': tag}
		config['data_info'][tag].update(stat_args_global)
		for i in xrange(1,len(l)):
			config['data_info'][tag][l[i][0]] = l[i][1]
	for tag in config['data_info']:
		for k in ['marker_col','freq_col','rsq_col','hwe_col','effect_col','stderr_col','or_col','z_col','p_col','n_col']:
			if not k in config['data_info'][tag]:
				config['data_info'][tag][k] = None
	return config

def PrintMetaOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if not k in ['file_order','meta_order','meta_info','data_info']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['data_info']:
		print '   dataset ' + str(m) + ' ...'
		for n in cfg['data_info'][m]:
			if not n in ['model_fxn','model','family','format','data','write_eof','write_header']:
				if cfg['data_info'][m][n] is not None and cfg['data_info'][m][n] is not False:
					if cfg['data_info'][m][n] is True:
						print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['data_info'][m].keys()],key=len)))
					else:
						print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['data_info'][m].keys()],key=len))) + " " + str(cfg['data_info'][m][n])
	if len(cfg['meta_order']) > 0:
		print '   meta analysis ...'
		for m in cfg['meta_info']:
			print "      {0:>{1}}".format(str('--meta ' + m), len(max(['--' + k for k in cfg['meta_info']],key=len))) + ":" + str('+'.join(cfg['meta_info'][m]))
	print ''

def GenerateEvalCfg(args, ini):
	config = {'file': None, 'qq': False, 'qq_strat': False, 'qq_n': False, 'mht': False, 'color': False, 'plot_gc': False,
				'set_gc': None, 'ext': 'tiff', 'sig': 0.000000054, 'stat': 'marker', 'top': None, 'region': None, 'region_id': None, 
				'region_file': None, 'pmax': 1e-4, 'tag': None, 'unrel': False, 'lz_source': None, 'lz_build': None, 'lz_pop': None, 
				'callrate': None, 'rsq': None, 'maf': None, 'hwe': None, 'effect': None, 'stderr': None, 'odds_ratio': None, 'df': None}
	for arg in args:
		config[arg[0]] = arg[1]
	config['locuszoom'] = ini.get('main','locuszoom')
	return config

def PrintEvalOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	print ''

def GenerateGcCfg(args, ini):
	config = {'file': None, 'gc': {}}
	for arg in args:
		if arg[0] == 'gc':
			config['gc'][arg[1][0]] = arg[1][1]
		else:
			config[arg[0]] = arg[1]
	return config

def PrintGcOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	print ''

def GenerateAnnotCfg(args, ini):
	config = {'file': None, 'build': 'GRCh37.75', 'replace': False, 'qsub': None, 'snpeff': None, 'snpsift': None, 'dbnsfp': None}
	for arg in args:
		config[arg[0]] = arg[1]
	config['snpeff'] = ini.get('main','snpeff')
	config['snpsift'] = ini.get('main','snpsift')
	config['dbnsfp'] = ini.get('main','dbnsfp')
	return config

def PrintAnnotOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	print ''
