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
import __version__
import textwrap
from datetime import date
import Args

def version():
	print ''
	print 'uga v' + __version__.version + ' (c) 2015 Ryan Koesterer   GNU General Public License v3'
	print ''
	print textwrap.fill("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.", initial_indent='   ', subsequent_indent='   ')
	print ''
	print textwrap.fill("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.", initial_indent='   ', subsequent_indent='   ')
	print ''
	print textwrap.fill("You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>", initial_indent='   ', subsequent_indent='   ')
	print ''

def get_parser():
	main_parser = argparse.ArgumentParser()
	main_parser.add_argument('--version', 
						action='version', 
						version=version(), 
						help='display version information and exit')
	subparsers = main_parser.add_subparsers(title='modules', dest='which')
	global_parser = argparse.ArgumentParser(add_help=False)

	set_parser = Args.set_args(subparsers.add_parser('set', help='user definable settings', parents=[global_parser]))
	snv_parser = Args.snv_args(subparsers.add_parser('snv', help='run single nucleotide variant association models', parents=[global_parser]))
	snvgroup_parser = Args.snvgroup_args(subparsers.add_parser('snvgroup', help='run variant group (ie. gene) based association models', parents=[global_parser]))
	meta_parser = Args.meta_args(subparsers.add_parser('meta', help='run meta analysis', parents=[global_parser]))
	compile_parser = Args.compile_args(subparsers.add_parser('compile', help='verify and compile results files', parents=[global_parser]))
	resubmit_parser = Args.resubmit_args(subparsers.add_parser('resubmit', help='resubmit failed results', parents=[global_parser]))
	snvplot_parser = Args.snvplot_args(subparsers.add_parser('snvplot', help='generate qq and manhattan plots for single variant results', parents=[global_parser]))
	snvgroupplot_parser = Args.snvgroupplot_args(subparsers.add_parser('snvgroupplot', help='generate qq and manhattan plots for grouped variant results', parents=[global_parser]))

	return main_parser

def get_args(parser):
	args=parser.parse_args()
	if 'ordered_args' in args:
		for k in args.ordered_args:
			vars(args).update({k[0]: k[1]})
	else:
		if args.which in ['snv','snvgroup','meta','map','compile','eval','gc','annot']:
			parser.error("missing argument: no options selected")
	print ''
	print 'active module: ' + args.which
	return args

def generate_snv_cfg(args):
	config = {'out': None, 'buffer': 100, 'region': None, 'region_file': None, 'cpus': 1, 'mb': 1, 'qsub': None, 'split': False, 'split_n': None, 'replace': False, 
					'debug': False, 'models': {}, 'model_order': [], 'meta': {}, 'meta_order': [], 'meta_type': {}}
	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'buffer' and arg[1] is not None:
			config['buffer'] = arg[1]
		if arg[0] == 'region':
			config['region'] = arg[1]
		if arg[0] == 'region_file':
			config['region_file'] = arg[1]
		if arg[0] == 'meta_stderr':
			config['meta'][arg[1][0]] = arg[1][1]
			config['meta_type'][arg[1][0]] = 'stderr'
			config['meta_order'].append(arg[1][0])
		if arg[0] == 'meta_sample_size':
			config['meta'][arg[1][0]] = arg[1][1]
			config['meta_type'][arg[1][0]] = 'sample_size'
			config['meta_order'].append(arg[1][0])
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
		if arg[0] == 'replace':
			config['replace'] = arg[1]
		if arg[0] == 'debug':
			config['debug'] = arg[1]

	args = [x for x in args if x[0] not in config and not x[0] in ['meta_sample_size','meta_stderr']]

	tags_idx = [args.index((x,y)) for x, y in args if x == 'tag'] + [len(args)]
	global_args = [x for x in args[:tags_idx[0]] if x[0] not in config]
	config_default = {'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
							'male': 1, 'female': 2, 'miss': 0.0, 'maf': 0.0, 'maxmaf': 1.0, 'mac': 0.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
							'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'corstr': None, 
							'pheno': None, 'covars': None, 'covars_categorical': None, 'interact': None, 'reverse': False}
	if len(tags_idx) > 1:
		for i in xrange(len(tags_idx[:-1])):
			config['models'][args[tags_idx[i]][1]] = config_default.copy()
			for arg in global_args:
				if arg[0] in ['score','lm','glm','gee']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['score','lm','glm','gee']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in global_args:
				if arg[0] in ['score','lm','glm','gee']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
				elif arg[0] in ['vcf','dos','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['score','lm','glm','gee']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
				elif arg[0] in ['vcf','dos','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			config['model_order'].append(args[tags_idx[i]][1])
	else:
		config['models']['___no_tag___'] = config_default
		for arg in global_args:
			if arg[0] in ['score','lm','glm','gee']:
				config['models']['___no_tag___']['case_code'] = 1
				config['models']['___no_tag___']['ctrl_code'] = 0
		for arg in global_args:
			if arg[0] in ['score','lm','glm','gee']:
				config['models']['___no_tag___']['fxn'] = arg[0]
			elif arg[0] in ['vcf','dos','oxford']:
				config['models']['___no_tag___']['format'] = arg[0]
				config['models']['___no_tag___']['file'] = arg[1]
			else:
				config['models']['___no_tag___'][arg[0]] = arg[1]
		config['model_order'].append('___no_tag___')
	return config

def print_snv_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if not k in ['models','model_order','meta','meta_order','meta_type']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['model_order']:
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
			if cfg['meta_type'][m] == 'stderr':
				print "      {0:>{1}}".format(str('--meta-stderr'), len(max(['--' + k for k in cfg['meta'].keys()],key=len))) + " " + m + ' ' + str(cfg['meta'][m])
			else:
				print "      {0:>{1}}".format(str('--meta-sample-size'), len(max(['--' + k for k in cfg['meta'].keys()],key=len))) + " " + m + ' ' + str(cfg['meta'][m])

def generate_snvgroup_cfg(args):
	config = {'out': None, 'buffer': 100, 'region': None, 'region_file': None, 'cpus': 1, 'qsub': None, 'split': False, 'split_n': None, 'replace': False, 'snvgroup_map': None, 
					'debug': False, 'timeout': 3600, 'models': {}, 'model_order': [], 'meta': {}, 'meta_order': []}
	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'buffer' and arg[1] is not None:
			config['buffer'] = arg[1]
		if arg[0] == 'meta':
			config['meta'][arg[1][0]] = arg[1][1]
			config['meta_order'].append(arg[1][0])
		if arg[0] == 'cpus' and arg[1] is not None:
			config['cpus'] = arg[1]
		if arg[0] == 'qsub':
			config['qsub'] = arg[1]
		if arg[0] == 'split' and arg[1] is True:
			config['split'] = arg[1]
		if arg[0] == 'split_n':
			config['split_n'] = arg[1]
		if arg[0] == 'snvgroup_map':
			config['snvgroup_map'] = arg[1]
		if arg[0] == 'region':
			config['region'] = arg[1]
		if arg[0] == 'region_file':
			config['region_file'] = arg[1]
		if arg[0] == 'replace':
			config['replace'] = arg[1]
		if arg[0] == 'debug':
			config['debug'] = arg[1]
		if arg[0] == 'timeout':
			config['timeout'] = arg[1]

	args = [x for x in args if x[0] not in config]

	tags_idx = [args.index((x,y)) for x, y in args if x == 'tag'] + [len(args)]
	global_args = [x for x in args[:tags_idx[0]] if x[0] not in config]
	config_default = {'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
							'male': 1, 'female': 2, 'miss': 0.0, 'maf': 0.0, 'maxmaf': 1.0, 'mac': 0.0, 'snvgroup_mac': 0.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
							'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'skat_wts': None, 'burden_wts': None, 'skat_method': None,
							'pheno': None, 'covars': None, 'covars_categorical': None, 'mafrange': None, 'skato_rho': None}
	if len(tags_idx) > 1:
		for i in xrange(len(tags_idx[:-1])):
			config['models'][args[tags_idx[i]][1]] = config_default.copy()
			for arg in global_args:
				if arg[0] in ['skat','skato','burden']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['skat','skato','burden']:
					config['models'][args[tags_idx[i]][1]]['case_code'] = 1
					config['models'][args[tags_idx[i]][1]]['ctrl_code'] = 0
			for arg in global_args:
				if arg[0] in ['skat','skato','burden']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
				elif arg[0] in ['vcf','dos','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			for arg in args[tags_idx[i]+1:tags_idx[i+1]]:
				if arg[0] in ['skat','skato','burden']:
					config['models'][args[tags_idx[i]][1]]['fxn'] = arg[0]
				elif arg[0] in ['vcf','dos','oxford']:
					config['models'][args[tags_idx[i]][1]]['format'] = arg[0]
					config['models'][args[tags_idx[i]][1]]['file'] = arg[1]
				else:
					config['models'][args[tags_idx[i]][1]][arg[0]] = arg[1]
			config['model_order'].append(args[tags_idx[i]][1])
	else:
		config['models']['___no_tag___'] = config_default
		for arg in global_args:
			if arg[0] in ['skat','skato','burden']:
				config['models']['___no_tag___']['case_code'] = 1
				config['models']['___no_tag___']['ctrl_code'] = 0
		for arg in global_args:
			if arg[0] in ['skat','skato','burden']:
				config['models']['___no_tag___']['fxn'] = arg[0]
			elif arg[0] in ['vcf','dos','oxford']:
				config['models']['___no_tag___']['format'] = arg[0]
				config['models']['___no_tag___']['file'] = arg[1]
			else:
				config['models']['___no_tag___'][arg[0]] = arg[1]
		config['model_order'].append('___no_tag___')
	return config

def print_snvgroup_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if not k in ['models','model_order','meta','meta_order']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['model_order']:
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
			print "      {0:>{1}}".format(str('--meta'), len(max(['--' + k for k in cfg['meta'].keys()],key=len))) + " " + m + ' ' + str(cfg['meta'][m])

def generate_meta_cfg(args):
	config = {'out': None, 'region': None, 'region_file': None, 'buffer': 100, 'cpus': 1, 'mb': 1, 'qsub': None, 'split': False, 'split_n': None, 'replace': False, 
					'debug': False, 'files': {}, 'file_order': [], 'meta': {}, 'meta_order': [], 'meta_type': {}}

	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'region':
			config['region'] = arg[1]
		if arg[0] == 'region_file':
			config['region_file'] = arg[1]
		if arg[0] == 'meta_stderr':
			config['meta'][arg[1][0]] = arg[1][1]
			config['meta_type'][arg[1][0]] = 'stderr'
			config['meta_order'].append(arg[1][0])
		if arg[0] == 'meta_sample_size':
			config['meta'][arg[1][0]] = arg[1][1]
			config['meta_type'][arg[1][0]] = 'sample_size'
			config['meta_order'].append(arg[1][0])
		if arg[0] == 'buffer':
			config['buffer'] = arg[1]
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
		if arg[0] == 'replace':
			config['replace'] = arg[1]
		if arg[0] == 'debug':
			config['debug'] = arg[1]
		if arg[0] == 'file':
			config['files'][arg[1][0]] = arg[1][1]
			config['file_order'].append(arg[1][0])
	return config

def print_meta_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if not k in ['files','file_order','meta','meta_order','meta_type']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for f in cfg['file_order']:
		print '   file ' + str(f) + ' ...' if len(cfg['files']) > 1 else '   file ...'
		print "      {0:>{1}}".format('--file ' + f, len('--file')) + " " + cfg['files'][f]
	if len(cfg['meta_order']) > 0:
		print '   meta analysis ...'
		for m in cfg['meta_order']:
			if cfg['meta_type'][m] == 'stderr':
				print "      {0:>{1}}".format(str('--meta-stderr'), len(max(['--' + k for k in cfg['meta'].keys()],key=len))) + " " + m + ' ' + str(cfg['meta'][m])
			else:
				print "      {0:>{1}}".format(str('--meta-sample-size'), len(max(['--' + k for k in cfg['meta'].keys()],key=len))) + " " + m + ' ' + str(cfg['meta'][m])

def generate_compile_cfg(args):
	config = {'dir': None, 'replace': False}
	for arg in args:
		if arg[0] == 'dir':
			config['dir'] = arg[1]
		if arg[0] == 'replace':
			config['replace'] = arg[1]
	return config

def print_compile_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])

def generate_snvplot_cfg(args):
	config = {'file': None, 'qsub': None, 'replace': False, 'debug': False, 'out': None, 'ext': 'tiff', 'nogc': False, 'color': False, 'qq': False, 'qq_strat': False, 
				'mht': False, 'crop': 10, 'pcol': None}
	for arg in args:
		if arg[0] == 'file':
			config['file'] = arg[1]
		if arg[0] == 'qsub':
			config['qsub'] = arg[1]
		if arg[0] == 'replace':
			config['replace'] = arg[1]
		if arg[0] == 'debug':
			config['debug'] = arg[1]
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'ext':
			config['ext'] = arg[1]
		if arg[0] == 'nogc':
			config['nogc'] = arg[1]
		if arg[0] == 'color':
			config['color'] = arg[1]
		if arg[0] == 'qq':
			config['qq'] = arg[1]
		if arg[0] == 'qq_strat':
			config['qq_strat'] = arg[1]
		if arg[0] == 'mht':
			config['mht'] = arg[1]
		if arg[0] == 'crop':
			config['crop'] = arg[1]
		if arg[0] == 'pcol':
			config['pcol'] = arg[1]
	return config

def print_snvplot_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])

def generate_snvgroupplot_cfg(args):
	config = {'file': None, 'qsub': None, 'replace': False, 'debug': False, 'out': None, 'ext': 'tiff', 'nogc': False, 'color': False, 'qq': False, 
				'mht': False, 'crop': 10, 'pcol': None, 'cmaf': 0.0}
	for arg in args:
		if arg[0] == 'file':
			config['file'] = arg[1]
		if arg[0] == 'qsub':
			config['qsub'] = arg[1]
		if arg[0] == 'replace':
			config['replace'] = arg[1]
		if arg[0] == 'debug':
			config['debug'] = arg[1]
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'ext':
			config['ext'] = arg[1]
		if arg[0] == 'nogc':
			config['nogc'] = arg[1]
		if arg[0] == 'color':
			config['color'] = arg[1]
		if arg[0] == 'qq':
			config['qq'] = arg[1]
		if arg[0] == 'mht':
			config['mht'] = arg[1]
		if arg[0] == 'crop':
			config['crop'] = arg[1]
		if arg[0] == 'pcol':
			config['pcol'] = arg[1]
		if arg[0] == 'cmaf':
			config['cmaf'] = arg[1]
	return config

def print_snvgroupplot_options(cfg):
	print ''
	print "main options ..."
	for k in cfg:
		if cfg[k] is not None and cfg[k] is not False:
			if cfg[k] is True:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
			else:
				print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
