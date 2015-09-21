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
from collections import OrderedDict
from glob import glob
from ConfigParser import SafeConfigParser
from pkg_resources import resource_filename
import Parse
import Process
import Stdout
import IO
import Compile

def main(args=None):
	args = Parse.GetArgs(Parse.GetParser())

	##### read settings file #####
	ini = SafeConfigParser()
	ini.read(resource_filename('uga', 'settings.ini'))

	##### read cfg file into dictionary #####
	#if args.which == 'meta':
	#	print "reading configuration from file"
	#	from Parse import GenerateMetaCfg
	#	config = GenerateMetaCfg(args.ordered_args)
	#	args.cfg = 1
	#elif args.which == 'stat':
	#	print "preparing stat configuration"
	#	from parse import GenerateStatCfg
	#	config = GenerateStatCfg(args.ordered_args)
	#	args.cfg = 1
	#elif args.which == 'eval':
	#	print "preparing eval configuration"
	#	from Parse import GenerateEvalCfg
	#	config = GenerateEvalCfg(args.ordered_args, ini)
	#	args.cfg = 1
	#elif args.which == 'gc':
	#	print "preparing gc configuration"
	#	from Parse import GenerateGcCfg
	#	config = GenerateGcCfg(args.ordered_args, ini)
	#	args.cfg = 1
	#elif args.which == 'annot':
	#	print "preparing annot configuration"
	#	from Parse import GenerateAnnotCfg
	#	config = GenerateAnnotCfg(args.ordered_args, ini)
	#	args.cfg = 1
	#elif args.which == 'bglm':
	#	print "preparing bglm configuration"
	#	from Parse import GenerateBglmCfg
	#	config = GenerateBglmCfg(args.ordered_args)
	#	args.cfg = 1
	#elif args.which == 'bssmeta':
	#	print "preparing bssmeta configuration"
	#	from Parse import GenerateBssmetaCfg
	#	#config = GenerateBssmetaCfg(args.ordered_args)
	#	#args.cfg = 1
	#print config; return
	##### define region list #####
	if args.which in ['model','meta']:
		n = 1
		dist_mode = 'full'
		regions = IO.Regions(filename=args.region_list,region=args.region,id=None)
		if args.region_list:
			if args.split or args.split_n:
				if not args.split_n or args.split_n > len(regions.df.index):
					n = len(regions.df.index)
					dist_mode = 'split-list'
				else:
					n = args.split_n
					dist_mode = 'split-list-n'
			else:
				dist_mode = 'list'
			print " " + str(len(regions.df.index)) + " regions found"
		elif args.region:
			if len(args.region.split(':')) > 1:
				dist_mode = 'region'
			else:
				dist_mode = 'chr'
			n = 1
		else:
			n = 1

		##### get job list from file #####
		if 'jobs' in vars(args).keys() and args.jobs is not None:
			jobs = []
			with open(args.jobs) as f:
				lines = (line.rstrip() for line in f)
				lines = (line for line in lines if line)
				for line in lines:
					if line.find(':') != -1:
						jobs.append(regions.df['region'][regions.df['region'] == line].index[0])
					else:
						jobs.append(int(line))
			print "" + str(len(jobs)) + " jobs read from job list file"

		##### define output directory and update out file name #####
		directory = os.getcwd() + '/'
		if n > 1:
			if dist_mode == 'split-list':
				directory = directory + 'chr[CHR]/'
			elif dist_mode == 'split-list-n':
				directory = directory + 'list[LIST]/'
		if 'out' in vars(args).keys() and args.which in ['model','meta']:
			args.out = directory + args.out

		##### generate out file names for split jobs or chr/region specific jobs #####
		out_files = {}
		if dist_mode in ['chr','region','split-list','split-list-n']:
			out_files = IO.GenerateSubFiles(region_df = regions.df, f = args.out, dist_mode = dist_mode, n = n)

	##### get user home directory #####
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

	elif args.which == 'map':
		config = '{\'out\': \'' + args.out + '\''
		for x in ['oxford','dos1','dos2','plink','vcf','region','b','kb','mb','n','chr','shift_mb','shift_kb','shift_b']:
			if x in vars(args).keys() and not vars(args)[x] in [False,None]:
				if type(vars(args)[x]) is str:
					config = config + ',\'' + x + '\': \'' + str(vars(args)[x]) + '\''
				else:
					config = config + ',\'' + x + '\': ' + str(vars(args)[x])
		config = config + '}'
		cmd = 'memory_usage((' + args.which.capitalize() + ', (),' + config + '), interval=0.1)'
		if args.replace:
			for f in [args.out, args.out + '.map.log']:
				try:
					os.remove(f)
				except OSError:
					continue
		else:
			for f in [args.out, args.out + '.map.log']:
				if os.path.exists(f):
					print Stdout.PrintError("1 or more output files already exists (use --replace flag to replace)")
					return
		if args.qsub:
			Process.Qsub('qsub ' + args.qsub + ' -o ' + args.out + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			Process.Interactive(qsub_wrapper, cmd, args.out + '.' + args.which + '.log')

	elif args.which == 'compile':
		out_files = {}
		if args.which == "compile":
			out_files = FindSubFiles(f = args.file)
		if len(out_files.keys()) > 1:
			existing_files = glob(args.file + '.gz*') + glob(args.file + '.log')
			if len(existing_files) > 0:
				if not args.replace:
					print Stdout.PrintError("1 or more output files or files with similar basename already exists (use --replace flag to replace)")
					return
				else:
					for f in existing_files:
						try:
							os.remove(f)
						except OSError:
							continue
			complete, complete_reg = Compile.CheckResults(file_dict=out_files, out=args.file + '.verify')
			if not complete:
				print Stdout.PrintError("results could not be verified")
				return
			out_files = OrderedDict([(x, out_files[x]) for x in out_files.keys() if str(x) in complete_reg])
			if args.split:
				compile_fxn = 'Compile.CompileResultsSplit'
			elif args.split_chr:
				compile_fxn = 'Compile.CompileResultsSplitChr'
			else:
				compile_fxn = 'Compile.CompileResults'
			if not globals()[compile_fxn](out_files, args.file):
				print Stdout.PrintError("results could not be compiled")
				return
		else:
			print Stdout.PrintError("no split results to compile")
			return

	elif args.which == 'eval':
		check_files = [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.eval.log']
		check_files = [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.top_results']
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.tiff'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.tiff'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.eps'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.eps'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq.pdf'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.qq_strat.pdf'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.tiff'] if 'mht' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.eps'] if 'mht' in vars(args).keys() else check_files
		check_files = check_files + [args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.mht.pdf'] if 'mht' in vars(args).keys() else check_files
		existing_files = []
		for f in check_files:
			if os.path.exists(f):
				if not args.replace:
					existing_files = existing_files + [f]
					print "found file " + str(f)
				else:
					try:
						os.remove(f)
					except OSError:
						continue
		if len(existing_files) > 0:
			print Stdout.PrintError("above files already exist (use --replace flag to replace)")
			return
		cmd = 'memory_usage((' + args.which.capitalize() + ', (' + str(config) + ',)), interval=0.1)'
		if args.qsub:
			Process.Qsub('qsub ' + args.qsub + ' -o ' + config['file'].replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.' + args.stat.replace('*','_inter_') + '.eval.log')

	elif args.which == 'gc':
		check_files = [args.file.replace('.gz','') + '.gc.log']
		check_files = check_files + [args.file.replace('.gz','') + '.gc.gz']
		check_files = check_files + [args.file.replace('.gz','') + '.gc.gz.tbi']
		existing_files = []
		for f in check_files:
			if os.path.exists(f):
				if not args.replace:
					existing_files = existing_files + [f]
					print "found file " + str(f)
				else:
					try:
						os.remove(f)
					except OSError:
						continue
		if len(existing_files) > 0:
			print Stdout.PrintError("above files already exist (use --replace flag to replace)")
			return
		cmd = 'memory_usage((' + args.which.upper() + ', (' + str(config) + ',)), interval=0.1)'
		if args.qsub:
			Process.Qsub('qsub ' + args.qsub + ' -o ' + config['file'].replace('.gz','') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.gc.log')

	elif args.which == 'annot':
		check_files = [args.file.replace('.gz','') + '.annot.xlsx']
		check_files = check_files + [args.file.replace('.gz','') + '.annot.log']
		check_files = check_files + [args.file.replace('.gz','') + '.annot1']
		check_files = check_files + [args.file.replace('.gz','') + '.annot2']
		check_files = check_files + [args.file.replace('.gz','') + '.annot3']
		check_files = check_files + [args.file.replace('.gz','') + '.annot.summary.genes.txt']
		check_files = check_files + [args.file.replace('.gz','') + '.annot.summary.html']
		existing_files = []
		for f in check_files:
			if os.path.exists(f):
				if not args.replace:
					existing_files = existing_files + [f]
					print "found file " + str(f)
				else:
					try:
						os.remove(f)
					except OSError:
						continue
		if len(existing_files) > 0:
			print Stdout.PrintError("above files already exist (use --replace flag to replace)")
			return
		cmd = 'memory_usage((' + args.which.capitalize() + ', (' + str(config) + ',)), interval=0.1)'
		if args.qsub:
			Process.Qsub('qsub ' + args.qsub + ' -o ' + config['file'].replace('.gz','') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			Process.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.annot.log')

	elif args.which in ['model','meta']:
		print "preparing output directories"
		if dist_mode == 'split-list' and n > 1:
			PrepareChrDirs(regions.df['region'], directory)
		elif dist_mode == 'split-list-n' and n > 1:
			PrepareListDirs(n, directory)
		if args.qsub:
			print "submitting jobs\n"
		joblist = []
		if not args.job is None:
			joblist.append(args.job)
		elif not args.jobs is None:
			joblist.extend(jobs)
		else:
			joblist.extend(range(n))
		for i in joblist:
			if dist_mode in ['split-list', 'region']:
				args.out = out_files['%s:%s-%s' % (str(regions.df['chr'][i]), str(regions.df['start'][i]), str(regions.df['end'][i]))]
				args.ordered_args = [x if x[0] != 'out' else ('out',args.out) for x in args.ordered_args]
				if n > 1:
					args.region = '%s:%s-%s' % (str(regions.df['chr'][i]), str(regions.df['start'][i]), str(regions.df['end'][i]))
					args.region_list = None
					args.ordered_args = [x if x[0] != 'region' else ('region',args.region) for x in args.ordered_args]
					args.ordered_args = [x if x[0] != 'region_list' else ('region_list',args.region_list) for x in args.ordered_args]
			elif dist_mode == 'split-list-n':
				args.out = out_files[i]
				args.ordered_args = [x if x[0] != 'out' else ('out',args.out) for x in args.ordered_args]
				rlist = args.out + '.region_list'
				regions.df.loc[np.array_split(np.array(regions.df.index), n)[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'id'])
				args.region_list = rlist
				args.ordered_args = [x if x[0] != 'region_list' else ('region_list',args.region_list) for x in args.ordered_args]
			elif dist_mode == 'chr':
				args.out = out_files['%s' % (str(regions.df['chr'][i]))]
				args.ordered_args = [x if x[0] != 'out' else ('out',args.out) for x in args.ordered_args]
			if args.replace:
				for f in [args.out, args.out + '.' + args.which + '.log',args.out + '.gz', args.out + '.gz.tbi']:
					try:
						os.remove(f)
					except OSError:
						continue
			else:
				for f in [args.out,args.out + '.log',args.out + '.gz', args.out + '.gz.tbi']:
					if not os.path.exists(f):
						print Stdout.PrintError("1 or more output files already exists (use --replace flag to replace)")
						return
			cmd = 'RunModels(' + str(args.ordered_args) + ')'
			if args.qsub:
				Process.Qsub(['qsub'] + args.qsub.split() + ['-o',args.out + '.' + args.which + '.log',qsub_wrapper],'\"' + cmd + '\"')
			else:
				Process.Interactive(qsub_wrapper, cmd, args.out + '.' + args.which + '.log')
	else:
		print Stdout.PrintError(args.which + " not a module")

	print ''

if __name__ == "__main__":
	main()
	os._exit(0)

