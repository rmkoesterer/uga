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
import glob
from ConfigParser import SafeConfigParser
from pkg_resources import resource_filename
import signal
import subprocess
import shutil
import Parse
import Process
import Map
import Fxns
import pickle
from Bio import bgzf

def main(args=None):
	rerun = []
	args = Parse.GetArgs(Parse.GetParser())

	if args.which in ['snv','snvgroup','meta','resubmit']:
		if args.which == 'resubmit':
			with open(args.dir + '/' + os.path.basename(args.dir) + '.args.pkl', 'rb') as p:
				args,cfg = pickle.load(p)
			with open(cfg['out'] + '/' + os.path.basename(cfg['out']) + '.rerun.pkl', 'rb') as p:
				rerun = pickle.load(p)
			cfg['replace'] = True
		else:
			cfg = getattr(Parse, 'Generate' + args.which.capitalize() + 'Cfg')(args.ordered_args)

	##### read settings file #####
	ini = SafeConfigParser()
	ini.read(resource_filename('uga', 'settings.ini'))

	##### distribute jobs #####
	if args.which in ['snv','snvgroup','meta']:
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
				regions_df = pd.read_table(cfg['region_file'],header=None,names=['region'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
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

		if args.which == 'snvgroup':
			#	generate regions dataframe with M rows from --snvgroup-map
			#	run_type = 0:   run as single job
			#	run_type = 1:   --cpus C (distribute M snvgroups over C cpus and run single job, 1 job C cpus)
			#	run_type = 10:  --split (split M snvgroups into single region jobs, M jobs 1 cpu)
			#	run_type = 100: --split-n N (distribute M snvgroups over N jobs, N jobs 1 cpu)
			#	run_type = 101: --split-n N, --cpus C (distribute M snvgroups over N jobs and distribute each job over C cpus, N jobs C cpus)

			if cfg['region_file']:
				regions_df = pd.read_table(cfg['region_file'],header=None,names=['region','id'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
				regions_df['chr'] = [x.split(':')[0] for x in regions_df['region']]
				regions_df['start'] = [x.split(':')[1].split('-')[0] for x in regions_df['region']]
				regions_df['end'] = [x.split(':')[1].split('-')[1] for x in regions_df['region']]
				regions_df['job'] = 1
				regions_df['cpu'] = 1
				regions_df = regions_df[['chr','start','end','region','id','job','cpu']].reset_index(drop=True)
			elif cfg['region']:
				snv_map = []
				data_files = []
				for m in cfg['models']:
					if cfg['models'][m]['file'] not in data_files:
						snv_map.extend(Map.Map(out=cfg['out'] + '.' + m + '.regions', file=cfg['models'][m]['file'], mb = 1000, region = cfg['region'], format=cfg['models'][m]['format']))
						data_files.append(cfg['models'][m]['file'])
				snv_map = list(set(snv_map))
				regions_df = pd.DataFrame({'region': snv_map, 'chr': [x.split(':')[0] for x in snv_map], 'start': [int(x.split(':')[1].split('-')[0]) for x in snv_map], 'end': [int(x.split(':')[1].split('-')[1]) for x in snv_map]})
				regions_df.sort(columns=['chr','start'],inplace=True)
				regions_df.reset_index(drop=True,inplace=True)
				regions_df['id'] = cfg['region']
				regions_df['job'] = 1
				regions_df['cpu'] = 1
				del data_files
				del snv_map
				regions_df = regions_df[['chr','start','end','region','id','job','cpu']].reset_index(drop=True)
			else:
				if cfg['snvgroup_map']:
					snvgroup_map = pd.read_table(cfg['snvgroup_map'],header=None,names=['chr','pos','marker','id'], compression='gzip' if cfg['snvgroup_map'].split('.')[-1] == 'gz' else None)
					regions_df = snvgroup_map[['chr','pos','id']]
					regions_df=regions_df.groupby(['chr','id'])
					regions_df = regions_df.agg({'pos': [np.min,np.max]})
					regions_df.columns = ['start','end']
					regions_df['chr'] = regions_df.index.get_level_values('chr')
					regions_df['id'] = regions_df.index.get_level_values('id')
					regions_df.reset_index(drop=True,inplace=True)
					regions_df['region'] = regions_df.chr.map(str) + ':' + regions_df.start.map(str) + '-' + regions_df.end.map(str)
					regions_df['job'] = 1
					regions_df['cpu'] = 1
					regions_df = regions_df[['chr','start','end','region','id','job','cpu']].reset_index(drop=True)
					regions_df.drop_duplicates(inplace=True)
					regions_df.sort(columns=['chr','start'],inplace=True)
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
		elif run_type == 11 and args.which != 'snvgroup':
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

	if args.which in ['snv','snvgroup','meta']:
		if len(rerun) == 0:
			if int(max(regions_df['job'])) > 1 and cfg['qsub'] is not None:
				print 'detected run type ' + str(run_type) + ' ...'
				if 'mb' in cfg:
					print '   ' + str(regions_df.shape[0]) + ' regions of size ' + str(cfg['mb']) + 'mb detected'
				else:
					print '   ' + str(regions_df.shape[0]) + ' regions detected'
				print '   ' + str(int(max(regions_df['job']))) + ' jobs will be submitted'
				print '   <= ' + str(max(np.bincount(regions_df['job']))) + ' regions per job'
				print '   <= '  + str(int(max(regions_df['cpu']))) + ' cpus per job'
				print '   qsub options: ' + cfg['qsub']
				print '   output directory: ' + directory
				print '   replace: ' + str(cfg['replace'])
				input_var = None
				while input_var not in ['y','n','Y','N']:
					input_var = raw_input('\nsubmit jobs (yY/nN)? ')
				if input_var.lower() == 'n':
					print 'canceled by user'
					return

			directory = directory + '/' + os.path.basename(cfg['out'])
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

			with open(directory + '/' + os.path.basename(cfg['out']) + '.args.pkl', 'wb') as p:
				pickle.dump([args, cfg], p)

			if run_type in [10,11,100,101] and regions_df.shape[0] > 1:
				for j in range(1, int(max(regions_df['job'])) + 1):
					try:
						os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100))
					except OSError:
						pass
					try:
						os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j))
					except OSError:
						pass
					if args.which in ['snv','meta']:
						for chr in np.unique(regions_df['chr'][regions_df['job'] == j]):
							try:
								os.mkdir(directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/chr' + str(chr))
							except OSError:
								continue
				with open(directory + '/' + cfg['out'] + '.files', 'w') as jlist:
					for j in range(1, int(max(regions_df['job'])) + 1):
						if 'model_order' in cfg:
							for m in cfg['model_order']:
								jlist.write(str(j) + '\t' + cfg['out'] + '.' + m + '.gz' + '\t' + 'jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + cfg['out'] + '.job' + str(j) + '.' + m + '.gz\n')
						else:
							jlist.write(str(j) + '\t' + cfg['out'] + '.gz' + '\t' + 'jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + cfg['out'] + '.job' + str(j) + '.gz\n')
						if 'meta_order' in cfg:
							if len(cfg['meta_order']) > 0:
								jlist.write(str(j) + '\t' + cfg['out'] + '.meta.gz' + '\t' + 'jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + cfg['out'] + '.job' + str(j) + '.meta.gz\n')
			else:
				if args.which in ['snv','meta']:
					for chr in np.unique(regions_df['chr']):
						try:
							os.mkdir(directory + '/chr' + str(chr))
						except OSError:
							continue
		else:
			directory = directory + '/' + os.path.basename(cfg['out'])
			if int(max(regions_df['job'])) > 1 and cfg['qsub'] is not None:
				print 'detected resubmit ...'
				print '   ' + str(len(rerun)) + ' jobs will be submitted'
				print '   <= ' + str(max(np.bincount(regions_df['job']))) + ' regions per job'
				print '   <= '  + str(int(max(regions_df['cpu']))) + ' cpus per job'
				print '   qsub options: ' + cfg['qsub']
				print '   output directory: ' + directory
				print '   replace: ' + str(cfg['replace'])
				input_var = None
				while input_var not in ['y','n','Y','N']:
					input_var = raw_input('\nresubmit jobs (yY/nN)? ')
				if input_var.lower() == 'n':
					print 'canceled by user'
					return

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

	elif args.which in ['snv','snvgroup','resubmit']:
		if cfg['qsub']:
			print "submitting jobs\n"
		out = cfg['out']
		joblist = range(1, int(max(regions_df['job'])) + 1) if len(rerun) == 0 else rerun
		for j in joblist:
			regions_job_df = regions_df[regions_df['job'] == j].reset_index(drop=True)
			if int(max(regions_df['job'])) > 1:
				cfg['out'] = directory + '/jobs' + str(100 * ((j-1) / 100) + 1) + '-' + str(100 * ((j-1) / 100) + 100) + '/job' + str(j) + '/' + out + '.job' + str(j)
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			else:
				cfg['out'] = directory + '/' + out
				regions_job_df.to_csv(cfg['out'] + '.regions', index = False, header = True, sep='\t', na_rep='None')
			args.ordered_args = [('out',cfg['out']),('region_file',cfg['out'] + '.regions'),('cpus',int(max(regions_job_df['cpu'])))] + [x for x in args.ordered_args if x[0] not in ['out','region_file','cpus']]
			cmd = 'Run' + args.which.capitalize() + '(' + str(args.ordered_args) + ')'
			if cfg['qsub']:
				Process.Qsub(['qsub'] + cfg['qsub'].split() + ['-o',cfg['out'] + '.' + args.which + '.log',qsub_wrapper],'\"' + cmd + '\"')
			else:
				Process.Interactive(qsub_wrapper, cmd, cfg['out'] + '.' + args.which + '.log')

	elif args.which == 'compile':
		files = pd.read_table(args.dir + '/' + os.path.basename(args.dir) + '.files', names=['job','out','file'])
		complete, rerun = Fxns.VerifyResults(args.dir,files)
		if len(rerun) > 0:
			print Process.PrintError('detected ' + str(len(rerun)) + ' failed jobs\n       writing failed job numbers to file ' + args.dir + '/' + os.path.basename(args.dir) + '.rerun\n       execute  script with the --resubmit option to resubmit only the failed jobs')
			with open(args.dir + '/' + os.path.basename(args.dir) + '.rerun.pkl', 'wb') as p:
				pickle.dump(rerun, p)
		else:
			complete = Fxns.CompileResults(args.dir,files)
			if complete:
				input_var = None
				while input_var not in ['y','n','Y','N']:
					input_var = raw_input('delete job data directories (yY/nN)? ')
				if input_var.lower() == 'n':
					print 'canceled by user'
				else:
					print 'deleting job data directories'
					for d in glob.glob(args.dir + '/jobs*-*'):
						try:
							shutil.rmtree(d)
						except OSError:
							print Process.PrintError('unable to delete job data directory ' + d)
			else:
				print Process.PrintError('file compilation incomplete')

		#jobs['complete'] = 0
		#for f in jobs['file']:
		#	if os.path.exists(f):
		#		jobs.loc[jobs['file'] == f,'complete'] = 1
		#incomplete = np.unique(jobs.loc[jobs['complete'] == 0,'job'])
		#if len(incomplete) > 0:
		#	with open(directory + '/' + args.dir + '.rerun', 'w') as rerun:
		#		for x in incomplete:
		#			rerun.write(str(x) + '\n')
		#else:
		#	bgzfile = {}
		#	for o in np.unique(jobs['out']):
		#		bgzfile[o] = bgzf.BgzfWriter(o, 'wb')
		#	for j in np.unique(jobs['job']):
		#		if j == 1:
		#			for f in jobs.loc[jobs['job'] == j,'file']:
		#				p1 = subprocess.Popen(['zcat',f], stdout=subprocess.PIPE, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
		#				p2 = subprocess.Popen(['head','-1'], stdin=p1.stdout, stdout=subprocess.PIPE)
		#				bgzfile[jobs.loc[jobs['file'] == f,'out']].write(p2.communicate()[0])
		#				p1.wait()
		#				p2.wait()
		#		else:
		#			for f in jobs.loc[jobs['job'] == j,'file']:
		#				p1 = subprocess.Popen(['zcat',f], stdout=subprocess.PIPE, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
		#				p2 = subprocess.Popen(['sed','1d'], stdin=p1.stdout, stdout=subprocess.PIPE)
		#				bgzfile[jobs.loc[jobs['file'] == f,'out']].write(p2.communicate()[0])
		#				p1.wait()
		#				p2.wait()
		#	for o in np.unique(jobs['out']):
		#		bgzfile[o].close()
		#		
		#print jobs; return
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
	#	check_files = check_files + [args.file.replace('.gz','') + '.annot.summary.features.txt']
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

