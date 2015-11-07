#!/usr/bin/env python

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
import numpy as np
import pandas as pd
from collections import OrderedDict
from glob import glob
from ConfigParser import SafeConfigParser
from pkg_resources import resource_filename
import shutil
import Parse
import Process
import Map

def main(args=None):
	args = Parse.GetArgs(Parse.GetParser())
	cfg = getattr(Parse, 'Generate' + args.which.capitalize() + 'Cfg')(args.ordered_args)

	##### read settings file #####
	ini = SafeConfigParser()
	ini.read(resource_filename('uga', 'settings.ini'))

	##### distribute jobs #####
	run_type = 0
	if cfg['cpus'] is not None and cfg['cpus'] > 1:
		run_type = run_type + 1
	if cfg['split']:
		run_type = run_type + 10
	if cfg['split_n']:
		run_type = run_type + 100

	if args.which == 'snv':
		#	generate regions dataframe with M rows, either from --snv-map or by splitting data file or --snv-region according to --mb
		#	run_type = 0:   run as single job
		#	run_type = 1:   --cpus C (distribute M regions over C cpus and run single job, 1 job C cpus)
		#	run_type = 10:  --split (split M regions into single region jobs, M jobs 1 cpu)
		#	run_type = 100: --split-n N (distribute M regions over N jobs, N jobs 1 cpu)
		#	run_type = 11:  --split, --cpus C (split M regions into chunks of size M / C and run M jobs, M jobs C cpus)
		#	run_type = 101: --split-n N, --cpus C (distribute M regions over N jobs and distribute each over C cpus, N jobs C cpus)

		if cfg['region_file']:
			regions_df = pd.read_table(cfg['region_file'],header=None,names=['region','id'])
			regions_df['chr'] = [x.split(':')[0] for x in regions_df['region']]
			regions_df['start'] = [x.split(':')[1].split('-')[0] for x in regions_df['region']]
			regions_df['end'] = [x.split(':')[1].split('-')[1] for x in regions_df['region']]
			regions_df['job'] = 1
			regions_df['cpu'] = 1
		else:
			snv_map = []
			data_files = []
			for m in cfg['models']:
				if cfg['models'][m]['file'] not in data_files:
					snv_map.extend(Map.Map(out=cfg['out'] + '.' + m + '.regions', file=cfg['models'][m]['file'], mb = cfg['mb'], region = cfg['region'], format=cfg['models'][m]['format']))
					data_files.append(cfg['models'][m]['file'])
			snv_map = list(set(snv_map))
			regions_df = pd.DataFrame({'region': snv_map, 'chr': [x.split(':')[0] for x in snv_map], 'start': [int(x.split(':')[1].split('-')[0]) for x in snv_map], 'end': [int(x.split(':')[1].split('-')[1]) for x in snv_map]})
			regions_df.sort(columns=['chr','start'],inplace=True)
			regions_df.reset_index(drop=True,inplace=True)
			regions_df['job'] = 1
			regions_df['cpu'] = 1
			del data_files
			del snv_map
		regions_df = regions_df[['chr','start','end','region','job','cpu']]

	if args.which == 'gene':
		#	generate regions dataframe with M rows from --gene-map
		#	run_type = 0:   run as single job
		#	run_type = 1:   --cpus C (distribute M genes over C cpus and run single job, 1 job C cpus)
		#	run_type = 10:  --split (split M genes into single region jobs, M jobs 1 cpu)
		#	run_type = 100: --split-n N (distribute M genes over N jobs, N jobs 1 cpu)
		#	run_type = 101: --split-n N, --cpus C (distribute M genes over N jobs and distribute each job over C cpus, N jobs C cpus)
		gene_map = pd.read_table(cfg['gene_map'],header=None,names=['chr','pos','marker','gene'])
		regions_df = gene_map
		regions_df['start'] = [min(x) for x in [regions_df['pos'][regions_df['gene'] == y] for y in regions_df['gene']]]
		regions_df['end'] = [max(x) for x in [regions_df['pos'][regions_df['gene'] == y] for y in regions_df['gene']]]
		regions_df['region'] = regions_df.chr.map(str) + ':' + regions_df.start.map(str) + '-' + regions_df.end.map(str)
		regions_df['job'] = 1
		regions_df['cpu'] = 1
		regions_df = regions_df[['chr','start','end','region','gene','job','cpu']]
		regions_df.drop_duplicates(inplace=True)
		regions_df.reset_index(drop=True,inplace=True)

	if run_type == 1:
		n = int(np.ceil(regions_df.shape[0] / float(cfg['cpus'])))
		n_remain = int(regions_df.shape[0] - (n-1) * cfg['cpus'])
		regions_df['cpu'] = np.append(np.repeat(range(cfg['cpus'])[:n_remain],n),np.repeat(range(cfg['cpus'])[n_remain:],n-1)) + 1
	elif run_type == 10:
		regions_df['job'] = regions_df.index.values + 1
	elif run_type == 100:
		n = int(np.ceil(regions_df.shape[0] / float(cfg['split_n'])))
		n_remain = int(regions_df.shape[0] - (n-1) * cfg['split_n'])
		regions_df['job'] = np.append(np.repeat(range(cfg['split_n'])[:n_remain],n),np.repeat(range(cfg['split_n'])[n_remain:],n-1)) + 1
	elif run_type == 11 and args.which != 'gene':
		cfg['split_n'] = int(np.ceil(regions_df.shape[0] / float(cfg['cpus'])))
		n = int(np.ceil(regions_df.shape[0] / float(cfg['split_n'])))
		n_remain = int(regions_df.shape[0] - (n-1) * cfg['split_n'])
		regions_df['job'] = np.append(np.repeat(range(cfg['split_n'])[:n_remain],n),np.repeat(range(cfg['split_n'])[n_remain:],n-1)) + 1
		for i in range(1,int(max(regions_df['job'])) + 1):
			n = int(np.ceil(regions_df[regions_df['job'] == i].shape[0] / float(cfg['cpus'])))
			n_remain = int(regions_df[regions_df['job'] == i].shape[0] - (n-1) * cfg['cpus'])
			regions_df.loc[regions_df['job'] == i,'cpu'] = np.append(np.repeat(range(cfg['cpus'])[:n_remain],n),np.repeat(range(cfg['cpus'])[n_remain:],n-1)) + 1
		cfg['split'] = None
	elif run_type == 101:
		n = int(np.ceil(regions_df.shape[0] / float(cfg['split_n'])))
		n_remain = int(regions_df.shape[0] - (n-1) * cfg['split_n'])
		regions_df['job'] = np.append(np.repeat(range(cfg['split_n'])[:n_remain],n),np.repeat(range(cfg['split_n'])[n_remain:],n-1)) + 1
		for i in range(1,int(max(regions_df['job'])) + 1):
			n = int(np.ceil(regions_df[regions_df['job'] == i].shape[0] / float(cfg['cpus'])))
			n_remain = int(regions_df[regions_df['job'] == i].shape[0] - (n-1) * cfg['cpus'])
			regions_df.loc[regions_df['job'] == i,'cpu'] = np.append(np.repeat(range(cfg['cpus'])[:n_remain],n),np.repeat(range(cfg['cpus'])[n_remain:],n-1)) + 1
	if int(max(regions_df['job'])) + 1 > 100000:
		print Process.PrintError('number of jobs exceeds 100,000, consider using --split-n to reduce the total number of jobs')
		return

	##### generate output directories #####
	directory = os.getcwd()

	if args.which in ['snv','gene','meta']:
		directory = directory + '/' + cfg['out']
		if os.path.exists(directory):
			if args.replace:
				print 'replacing existing results'
				try:
					shutil.rmtree(directory)
				except OSError:
					print Process.PrintError('unable to replace results directory' + directory)
			else:
				print Process.PrintError('results directory ' + directory + ' already exists, use --replace to overwrite existing results')
				return
		try:
			os.mkdir(directory)
		except OSError:
			pass

		if run_type in [10,11,100,101]:
			for j in range(1, int(max(regions_df['job'])) + 1):
				try:
					os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100))
				except OSError:
					pass
				try:
					os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j))
				except OSError:
					pass
				for chr in np.unique(regions_df['chr'][regions_df['job'] == j]):
					try:
						os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/chr' + str(chr))
					except OSError:
						continue
		else:
			for chr in np.unique(regions_df['chr']):
				try:
					os.mkdir(directory + '/chr' + str(chr))
				except OSError:
					continue

	##### locate qsub wrapper #####
	qsub_wrapper = resource_filename('uga', 'Qsub.py')

	if args.which == 'set':
		if 'ordered_args' in args:
			for k in args.ordered_args:
				ini.set('main',k[0],k[1])
			with open(resource_filename('uga', 'settings.ini'), 'w') as f:
				ini.write(f)
		print 'main settings ...'
		for s in ini.sections():
			for k in ini.options(s):
				print '   ' + k + ' = ' + ini.get(s,k)

	elif args.which == 'snv':
		if cfg['qsub']:
			print "submitting jobs\n"
		out = cfg['out']
		for j in range(1, int(max(regions_df['job'])) + 1):
			regions_job_df = regions_df[regions_df['job'] == j].reset_index(drop=True)
			if int(max(regions_df['job'])) > 1:
				cfg['out'] = directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + out + '.job' + str(j)
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			else:
				cfg['out'] = directory + '/' + out
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			args.ordered_args = [('out',cfg['out']),('region_file',cfg['out'] + '.regions'),('cpus',int(max(regions_job_df['cpu'])))] + [x for x in args.ordered_args if x[0] not in ['out','region_file','cpus']]
			cmd = 'RunSnv(' + str(args.ordered_args) + ')'
			if cfg['qsub']:
				Process.Qsub(['qsub'] + cfg['qsub'].split() + ['-o',cfg['out'] + '.' + args.which + '.log',qsub_wrapper],'\"' + cmd + '\"')
			else:
				Process.Interactive(qsub_wrapper, cmd, cfg['out'] + '.' + args.which + '.log')

	elif args.which == 'gene':
		if cfg['qsub']:
			print "submitting jobs\n"
		out = cfg['out']
		for j in range(1, int(max(regions_df['job'])) + 1):
			regions_job_df = regions_df[regions_df['job'] == j].reset_index(drop=True)
			if int(max(regions_df['job'])) > 1:
				cfg['out'] = directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + out + '.job' + str(j)
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			else:
				cfg['out'] = directory + '/' + out
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			args.ordered_args = [('out',cfg['out']),('region_file',cfg['out'] + '.regions'),('cpus',int(max(regions_job_df['cpu'])))] + [x for x in args.ordered_args if x[0] not in ['out','region_file','cpus']]
			cmd = 'RunGene(' + str(args.ordered_args) + ')'
			if cfg['qsub']:
				Process.Qsub(['qsub'] + cfg['qsub'].split() + ['-o',cfg['out'] + '.' + args.which + '.log',qsub_wrapper],'\"' + cmd + '\"')
			else:
				Process.Interactive(qsub_wrapper, cmd, cfg['out'] + '.' + args.which + '.log')

	#elif args.which == 'map':
	#	config = '{\'out\': \'' + cfg['out'] + '\''
	#	for x in ['oxford','dos1','dos2','plink','vcf','region','b','kb','mb','n','chr','shift_mb','shift_kb','shift_b']:
	#		if x in vars(args).keys() and not vars(args)[x] in [False,None]:
	#			if type(vars(args)[x]) is str:
	#				config = config + ',\'' + x + '\': \'' + str(vars(args)[x]) + '\''
	#			else:
	#				config = config + ',\'' + x + '\': ' + str(vars(args)[x])
	#	config = config + '}'
	#	cmd = 'memory_usage((' + args.which.capitalize() + ', (),' + config + '), interval=0.1)'
	#	if args.replace:
	#		for f in [cfg['out'], cfg['out'] + '.map.log']:
	#			try:
	#				os.remove(f)
	#			except OSError:
	#				continue
	#	else:
	#		for f in [cfg['out'], cfg['out'] + '.map.log']:
	#			if os.path.exists(f):
	#				print Process.PrintError("1 or more output files already exists (use --replace flag to replace)")
	#				return
	#	if cfg['qsub']:
	#		Process.Qsub('qsub ' + cfg['qsub'] + ' -o ' + cfg['out'] + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
	#	else:
	#		Process.Interactive(qsub_wrapper, cmd, cfg['out'] + '.' + args.which + '.log')

	#elif args.which == 'compile':
	#	out_files = {}
	#	if args.which == "compile":
	#		out_files = FindSubFiles(f = args.file)
	#	if len(out_files.keys()) > 1:
	#		existing_files = glob(args.file + '.gz*') + glob(args.file + '.log')
	#		if len(existing_files) > 0:
	#			if not args.replace:
	#				print Process.PrintError("1 or more output files or files with similar basename already exists (use --replace flag to replace)")
	#				return
	#			else:
	#				for f in existing_files:
	#					try:
	#						os.remove(f)
	#					except OSError:
	#						continue
	#		complete, complete_reg = Compile.CheckResults(file_dict=out_files, out=args.file + '.verify')
	#		if not complete:
	#			print Process.PrintError("results could not be verified")
	#			return
	#		out_files = OrderedDict([(x, out_files[x]) for x in out_files.keys() if str(x) in complete_reg])
	#		if cfg['split']:
	#			compile_fxn = 'Compile.CompileResultsSplit'
	#		elif cfg['split']_chr:
	#			compile_fxn = 'Compile.CompileResultsSplitChr'
	#		else:
	#			compile_fxn = 'Compile.CompileResults'
	#		if not globals()[compile_fxn](out_files, args.file):
	#			print Process.PrintError("results could not be compiled")
	#			return
	#	else:
	#		print Process.PrintError("no split results to compile")
	#		return

	#elif args.which == 'eval':
	#	check_files = [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.eval.log']
	#	check_files = [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.top_results']
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.tiff'] if 'qq' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.tiff'] if 'qq_strat' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.eps'] if 'qq' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.eps'] if 'qq_strat' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.pdf'] if 'qq' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.pdf'] if 'qq_strat' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.tiff'] if 'mht' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.eps'] if 'mht' in vars(args).keys() else check_files
	#	check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.pdf'] if 'mht' in vars(args).keys() else check_files
	#	existing_files = []
	#	for f in check_files:
	#		if os.path.exists(f):
	#			if not args.replace:
	#				existing_files = existing_files + [f]
	#				print "found file " + str(f)
	#			else:
	#				try:
	#					os.remove(f)
	#				except OSError:
	#					continue
	#	if len(existing_files) > 0:
	#		print Process.PrintError("above files already exist (use --replace flag to replace)")
	#		return
	#	cmd = 'memory_usage((' + args.which.capitalize() + ', (' + str(config) + ',)), interval=0.1)'
	#	if cfg['qsub']:
	#		Process.Qsub('qsub ' + cfg['qsub'] + ' -o ' + config['file'].replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
	#	else:
	#		Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.eval.log')

	#elif args.which == 'gc':
	#	check_files = [args.file.replace('.gz','') + '.gc.log']
	#	check_files = check_files + [args.file.replace('.gz','') + '.gc.gz']
	#	check_files = check_files + [args.file.replace('.gz','') + '.gc.gz.tbi']
	#	existing_files = []
	#	for f in check_files:
	#		if os.path.exists(f):
	#			if not args.replace:
	#				existing_files = existing_files + [f]
	#				print "found file " + str(f)
	#			else:
	#				try:
	#					os.remove(f)
	#				except OSError:
	#					continue
	#	if len(existing_files) > 0:
	#		print Process.PrintError("above files already exist (use --replace flag to replace)")
	#		return
	#	cmd = 'memory_usage((' + args.which.upper() + ', (' + str(config) + ',)), interval=0.1)'
	#	if cfg['qsub']:
	#		Process.Qsub('qsub ' + cfg['qsub'] + ' -o ' + config['file'].replace('.gz','') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
	#	else:
	#		Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.gc.log')

	#elif args.which == 'annot':
	#	check_files = [args.file.replace('.gz','') + '.annot.xlsx']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot.log']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot1']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot2']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot3']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot.summary.genes.txt']
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot.summary.html']
	#	existing_files = []
	#	for f in check_files:
	#		if os.path.exists(f):
	#			if not args.replace:
	#				existing_files = existing_files + [f]
	#				print "found file " + str(f)
	#			else:
	#				try:
	#					os.remove(f)
	#				except OSError:
	#					continue
	#	if len(existing_files) > 0:
	#		print Process.PrintError("above files already exist (use --replace flag to replace)")
	#		return
	#	cmd = 'memory_usage((' + args.which.capitalize() + ', (' + str(config) + ',)), interval=0.1)'
	#	if cfg['qsub']:
	#		Process.Qsub('qsub ' + cfg['qsub'] + ' -o ' + config['file'].replace('.gz','') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
	#	else:
	#		Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.annot.log')

	else:
		print Process.PrintError(args.which + " not a module")

	print ''

if __name__ == "__main__":
	main()
	os._exit(0)

