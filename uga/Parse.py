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

	##### Define Subparsers #####
	set_parser = SetArgs(subparsers.add_parser('set', help='user definable settings', parents=[global_parser]))

	##### Define Stat Parsers #####
	#bglm_parser = BGlmArgs(StatArgs(subparsers.add_parser('bglm', help='binomial generalized linear model', parents=[global_parser])))
	#gglm_parser = GGlmArgs(StatArgs(subparsers.add_parser('gglm', help='gaussian generalized linear model', parents=[global_parser])))
	#bssmeta_parser = BSsMetaArgs(StatArgs(subparsers.add_parser('bssmeta', help='binomial singlesnpMeta model', parents=[global_parser])))
	model_parser = ModelArgs(subparsers.add_parser('model', help='run association models', parents=[global_parser]))

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
		if args.which in ['model','meta','map','compile','eval','gc','annot']:
			parser.error("missing argument: no options selected")
	if args.which in ['model','bglm']:
		if args.region:
			assert not args.split, parser.error("argument -s/--split: not allowed with argument --region")
			assert not args.split_n, parser.error("argument -n/--split-n: not allowed with argument --region")
			assert not args.job, parser.error("argument -j/--job: not allowed with argument --region")
			assert not args.jobs, parser.error("argument --jobs: not allowed with argument --region")
			assert not args.split_chr, parser.error("argument --split-chr: not allowed with argument --region")
		if args.region_list:
			assert os.path.exists(args.region_list), parser.error("argument --region-list: file does not exist")
		if args.jobs:
			assert os.path.exists(args.jobs), parser.error("argument --jobs: file does not exist")
	if args.which in ['bglm'] and len([x for x in vars(args) if x in ['plink','vcf','dos1','dos2','oxford'] and vars(args)[x] is not None]) > 1:
		parser.error("multiple data files detected: only 1 of --plink, --vcf, --dos1, --dos2, and --oxford allowed for module " + args.which)
	if args.which == 'map' and not (args.b or args.kb or args.mb or args.n):
		parser.error("missing argument: --b, --kb, --mb, or --n required in module map")
	if args.which == 'stat' and not (args.fskato is None or args.fskat is None or args.fburden is None) and args.ped is None:
		parser.error("missing argument: --ped required for fskato, fskat, or fburden models using data with family structure")
	print ''
	print 'active module: ' + args.which
	return args
"""
def GenerateStatCfg(args):
	config = {'out': None, 'buffer': 100, 'region_list': None, 'region': None,'id': None, 'write_header': False, 'write_eof': False,
					'models': {}, 'model_order': [], 'meta': []}

	##### add top level variables to config
	top_args = [a for a in args if a[0] in ['out','buffer','region_list','region','id']]
	for arg in top_args:
		config[arg[0]] = arg[1]

	# list all possible model level arguments
	model_vars = ['tag','sample','pheno','variant_list','fid','iid','patid','matid','all_founders','focus','ped','sex', 
					'male','female','miss','maf','maxmaf','mac','rsq','hwe','hwe_maf','rho','case_code','ctrl_code','corstr','sep','lmer_ctrl',
					'reml','lrt','cph_ctrl','skat_wts','burden_wts','skat_method',
					'ggee','bgee','gglm','bglm','glme','blme','cph','neff','bssmeta',
					'fskato','gskato','bskato','fskat','gskat','bskat','fburden','gburden','bburden',
					'oxford','vcf','plink','dos1','dos2']
	if not 'tag' in [a[0] for a in args]:
		args = [('tag', 'NA')] + args
	pre_tag_idx = [a[0] for a in args if a[0] in model_vars].index('tag')

	# extract global model arguments
	stat_args_global = dict([a for a in args if a[0] in model_vars][:pre_tag_idx])
	for k in stat_args_global.keys():
		if k in ['ggee','bgee','gglm','bglm','glme','blme','cph','neff','bssmeta',
						'fskato','gskato','bskato','fskat','gskat','bskat','fburden','gburden','bburden'] and stat_args_global[k] is not None:
			stat_args_global['model'] = stat_args_global[k]
			stat_args_global['model_fxn'] = k
			if k in ['bgee','bglm','blme','bssmeta','bskato','bskat','bburden']:
				stat_args_global['family'] = "binomial"
			elif k in ['ggee','gglm','glme','gskato','gskat','gburden']:
				stat_args_global['family'] = "gaussian"
			else:
				stat_args_global['family'] = None
		elif k in ['oxford','vcf','plink','dos1','dos2'] and stat_args_global[k] is not None:
			stat_args_global['data'] = stat_args_global[k]
			stat_args_global['format'] = k
			if not 'variant_list' in stat_args_global:
				stat_args_global['variant_list'] = None
			if not 'sample' in stat_args_global:
				stat_args_global['sample'] = None
		else:
			stat_args_global[k] = stat_args_global[k]

	# extract local model arguments
	stat_args_local = [a for a in args if a[0] in model_vars][pre_tag_idx:]
	local_tag_idx = [a for a, b in enumerate(stat_args_local) if b[0] == 'tag']
	if len(local_tag_idx) == 0:
		local_tag_idx = [0]
	local_args = [stat_args_local[i:j] for i, j in zip([0]+local_tag_idx[1:], local_tag_idx[1:]+[None])]
	# extract meta arguments
	meta_args = [a for a in args if a[0] == 'meta']

	# loop over model argument sets and add to config['models']
	for l in local_args:
		if l[0][0] != 'tag':
			tag = 'NA'
		else:
			tag = l[0][1]
		config['model_order'].append(tag)
		config['models'][tag] = {'tag': tag}
		config['models'][tag].update(stat_args_global)
		for i in xrange(1,len(l)):
			config['models'][tag][l[i][0]] = l[i][1]
			if l[i][0] in ['ggee','bgee','gglm','bglm','glme','blme','cph','neff','bssmeta',
							'fskato','gskato','bskato','fskat','gskat','bskat','fburden','gburden','bburden'] and l[i][1] is not None:
				for k in [x for x in ['ggee','bgee','gglm','bglm','glme','blme','cph','neff','bssmeta',
							'fskato','gskato','bskato','fskat','gskat','bskat','fburden','gburden','bburden'] if x != l[i][0]]:
					if k in config['models'][tag]:
						del config['models'][tag][k]
				config['models'][tag]['model'] = l[i][1]
				config['models'][tag]['model_fxn'] = l[i][0]
				if l[i][0] in ['bgee','bglm','blme','bssmeta','bskato','bskat','bburden']:
					config['models'][tag]['family'] = "binomial"
				elif l[i][0] in ['ggee','gglm','glme','gskato','gskat','gburden']:
					config['models'][tag]['family'] = "gaussian"
				else:
					config['models'][tag]['family'] = None
			elif l[i][0] in ['oxford','vcf','plink','dos1','dos2'] and l[i][1] is not None:
				for k in [x for x in ['oxford','vcf','plink','dos1','dos2'] if x != l[i][0]]:
					if k in config['models'][tag]:
						del config['models'][tag][k]
				config['models'][tag]['data'] = l[i][1]
				config['models'][tag]['format'] = l[i][0]
		if not 'variant_list' in config['models'][tag]:
			config['models'][tag]['variant_list'] = None
		if not 'sample' in config['models'][tag]:
			config['models'][tag]['sample'] = None

	# loop over meta arguments and add to config['meta']
	for m in meta_args:
		config['meta'].append(m[1])

	# update all defaults and remove any extra settings
	for x in config['models']:
		if config['models'][x]['model_fxn'] in ['ggee','bgee']:
			if 'corstr' not in config['models'][x]:
				config['models'][x]['corstr'] = 'exchangeable'
		if config['models'][x]['model_fxn'] in ['gskato','bskato','fskato']:
			if not 'rho' in config['models'][x]:
				config['models'][x]['rho'] = 1
		else:
			if 'rho' in config['models'][x]:
				del config['models'][x]['rho']
		if config['models'][x]['model_fxn'] in ['gskato','bskato','fskato','fskat','gskat','bskat']:
			if not 'skat_wts' in config['models'][x]:
				config['models'][x]['skat_wts'] = 'function(maf){dbeta(maf,1,25)}'
			if not 'skat_method' in config['models'][x]:
				config['models'][x]['skat_method'] = 'saddlepoint'
		if config['models'][x]['model_fxn'] in ['gskato','bskato','fskato','fburden','gburden','bburden']:
			if not 'burden_wts' in config['models'][x]:
				config['models'][x]['burden_wts'] = 'function(maf){as.numeric(maf<0.01)}'
		if config['models'][x]['model_fxn'] in ['gskato','bskato','fskato','fburden','gburden','bburden','neff'] and config['region'] is not None:
			if config['id'] is None:
				config['id'] = 'NA'
		if config['models'][x]['model_fxn'] in ['glme','blme']:
			if not 'lmer_ctrl' in config['models'][x]:
				config['models'][x]['lmer_ctrl'] = "check.nobs.vs.rankZ='stop',check.nlev.gtreq.5='stop',check.rankX='stop.deficient',check.scaleX='stop',check.conv.grad=.makeCC('stop',tol=1e-3,relTol=NULL),check.conv.singular=.makeCC(action='stop',tol=1e-4),check.conv.hess=.makeCC(action='stop',tol=1e-6)"
			if not 'lrt' in config['models'][x]:
				config['models'][x]['lrt'] = False
		if config['models'][x]['model_fxn'] == 'glme':
			if not 'reml' in config['models'][x]:
				config['models'][x]['reml'] = False
		if config['models'][x]['model_fxn'] == 'cph':
			if not 'cph_ctrl' in config['models'][x]:
				config['models'][x]['cph_ctrl'] = "eps=1e-09,toler.chol=.Machine$double.eps^0.75,iter.max=20,toler.inf=sqrt(1e-09),outer.max=10"
		if not 'miss' in config['models'][x]:
			config['models'][x]['miss'] = None
		if not 'maf' in config['models'][x]:
			config['models'][x]['maf'] = None
		if not 'maxmaf' in config['models'][x]:
			config['models'][x]['maxmaf'] = None
		if not 'mac' in config['models'][x]:
			config['models'][x]['mac'] = None
		if not 'hwe' in config['models'][x]:
			config['models'][x]['hwe'] = None
		if not 'hwe_maf' in config['models'][x]:
			config['models'][x]['hwe_maf'] = None
		if not 'rsq' in config['models'][x]:
			config['models'][x]['rsq'] = None
		if not 'sep' in config['models'][x]:
			config['models'][x]['sep'] = 'tab'
		if not 'case_code' in config['models'][x]:
			config['models'][x]['case_code'] = 1
		if not 'ctrl_code' in config['models'][x]:
			config['models'][x]['ctrl_code'] = 0
		if not 'sex' in config['models'][x]:
			config['models'][x]['sex'] = None
		if not 'male' in config['models'][x]:
			config['models'][x]['male'] = 1
		if not 'female' in config['models'][x]:
			config['models'][x]['female'] = 2
	return config

def PrintStatOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if not k in ['models','model_order','meta']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	for m in cfg['models']:
		print '   model ' + str(m) + ' ...' if len(cfg['models']) > 1 else '   model ...'
		cfg['models'][m][cfg['models'][m]['model_fxn']] = cfg['models'][m]['model']
		for n in cfg['models'][m]:
			if not n in ['model_fxn','model','family','format','data','write_eof','write_header']:
				if cfg['models'][m][n] is not None and cfg['models'][m][n] is not False:
					if cfg['models'][m][n] is True:
						print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len)))
					else:
						print "      {0:>{1}}".format(str('--' + n.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg['models'][m].keys()],key=len))) + " " + str(cfg['models'][m][n])
	if len(cfg['meta']) > 0:
		print '   meta analysis ...'
		metas = [(x.split(':')[0],x.split(':')[1]) for x in cfg['meta']]
		for m in metas:
			print "      {0:>{1}}".format(str('--meta ' + m[0]), len(max(['--' + k[0] for k in metas],key=len))) + ":" + str(m[1])
	print ''

def GenerateBglmCfg(args):
	config = {'out': None, 'pheno': None, 'formula': None, 'variant_list': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sample': None, 
				'sex': None, 'male': 1, 'female': 2, 'buffer': 100, 'miss': 0.5, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None, 'id': None, 'region': None, 'region_list': None, 
				'oxford': None, 'dos1': None, 'dos2': None, 'plink': None, 'vcf': None, 'case_code': 1, 'ctrl_code': 0, 'write_header': True, 'write_eof': True, 'fxn': 'bglm'}
	for arg in args:
		if arg[0] in ['plink','vcf','dos1','dos2','oxford']:
			config['format'] = arg[0]
			config['file'] = arg[1]
		config[arg[0]] = arg[1]
	return config

def PrintBglmOptions(cfg):
	print ''
	print "options in effect ..."
	for k in cfg:
		if not k in ['format','file','write_eof','write_header']:
			if cfg[k] is not None and cfg[k] is not False:
				if cfg[k] is True:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
				else:
					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
	print ''

bssmeta_defaults = {'formula': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
														'male': 1, 'female': 2, 'miss': 0.5, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
														'case_code': 1, 'ctrl_code': 0, 'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'variant_list': None}
"""
def GenerateModelCfg(args):
	#config = {'out': None, 'pheno': None, 'formula': None, 'variant_list': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sample': None, 
	#			'sex': None, 'male': 1, 'female': 2, 'buffer': 100, 'miss': 0.5, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None, 'region': None, 'region_list': None, 
	#			'oxford': None, 'dos1': None, 'dos2': None, 'plink': None, 'vcf': None, 'case_code': 1, 'ctrl_code': 0, 'write_header': True, 'write_eof': True, 'fxn': 'bssmeta',
	config = {'out': None, 'buffer': 1000, 'region': None, 'region_list': None, 'id': None, 'models': {}, 'model_order': [], 'meta': {}, 'meta_order': [], 'replace': False}
	#config = {'out': None, 'formula': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 'male': 1, 'female': 2,
	#			'buffer': 100, 'miss': 0.5, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None, 'region': None, 'region_list': None, 
	#			'case_code': 1, 'ctrl_code': 0, 'write_header': True, 'write_eof': True, 'fxn': 'bssmeta',
	#			'tag': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'variant_list': None, 'metas': []}
	for arg in args:
		if arg[0] == 'out':
			config['out'] = arg[1]
		if arg[0] == 'buffer':
			config['buffer'] = arg[1]
		if arg[0] == 'region':
			config['region'] = arg[1]
		if arg[0] == 'region_list':
			config['region_list'] = arg[1]
		if arg[0] == 'meta':
			config['meta'][arg[1].split(':')[0]] = arg[1].split(':')[1]

	tags_idx = [args.index((x,y)) for x, y in args if x == 'tag'] + [len(args)]
	global_args = [x for x in args[:tags_idx[0]] if x[0] not in config]
	config_default = {'formula': None, 'fid': 'FID', 'iid': 'IID', 'patid': None, 'matid': None, 'all_founders': False, 'sep': 'tab', 'sex': None, 
							'male': 1, 'female': 2, 'miss': 0.5, 'maf': 0.000001, 'maxmaf': 1.0, 'mac': 3.0, 'rsq': 0.0, 'hwe': None, 'hwe_maf': None,
							'fxn': None, 'format': None, 'file': None, 'sample': None, 'pheno': None, 'variant_list': None}
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
		for arg in args:
			if arg[0] in ['bssmeta']:
				config['models']['___no_tag___']['case_code'] = 1
				config['models']['___no_tag___']['ctrl_code'] = 0
		for arg in args:
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

def PrintModelOptions(cfg):
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


#def PrintBssmetaOptions(cfg):
#	print ''
#	print "options in effect ..."
#	for k in cfg:
#		if not k in ['format','file']:
#			if cfg[k] is not None and cfg[k] is not False:
#				if cfg[k] is True:
#					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len)))
#				else:
#					print "      {0:>{1}}".format(str('--' + k.replace('_','-')), len(max(['--' + key.replace('_','-') for key in cfg.keys()],key=len))) + " " + str(cfg[k])
#	print ''

def GenerateMetaCfg(args):
	config = {'out': None, 'method': None, 'buffer': 100, 'region_list': None, 'region': None, 'id': None, 'write_header': False, 'write_eof': False,
					'data_info': {}, 'meta_info': {}, 'meta_order': [], 'file_order': []}

	##### add top level variables to config
	top_args = [a for a in args if a[0] in ['out','method','buffer','region_list','region','id']]
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
				'region_list': None, 'pmax': 1e-4, 'tag': None, 'unrel': False, 'lz_source': None, 'lz_build': None, 'lz_pop': None, 
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
