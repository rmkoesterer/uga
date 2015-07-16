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
import gzip
import sys
import string
import math
import subprocess
import numpy as np
import pandas as pd
import collections
import re
import glob
import tabix
from Bio import bgzf
import Parse
import multi_key_dict
import FileFxns
import SystemFxns
from ConfigParser import SafeConfigParser
from pkg_resources import resource_filename
pd.options.mode.chained_assignment = None

def main(args=None):
	parser=Parse.Parser()
	args=Parse.Parse(parser)

	##### read settings file #####
	ini = SafeConfigParser()
	ini.read(resource_filename('uga', 'settings.ini'))

	##### read cfg file into dictionary #####
	if args.which == 'meta':
		print "reading configuration from file"
		config = Parse.GenerateMetaCfg(args.ordered_args)
		args.cfg = 1
	elif args.which == 'model':
		print "preparing model configuration"
		config = Parse.GenerateModelCfg(args.ordered_args)
		args.cfg = 1
	elif args.which == 'explore':
		print "preparing explore configuration"
		config = Parse.GenerateExploreCfg(args.ordered_args, ini)
		args.cfg = 1
	elif args.which == 'gc':
		print "preparing gc configuration"
		config = Parse.GenerateGcCfg(args.ordered_args, ini)
		args.cfg = 1
	elif args.which == 'annot':
		print "preparing annot configuration"
		config = Parse.GenerateAnnotCfg(args.ordered_args, ini)
		args.cfg = 1

	##### define region list #####
	if args.which in ['model','meta']:
		n = 1
		dist_mode = 'full'
		if args.reglist:
			print "generating list of genomic regions ...", 
			region_df = FileFxns.LoadCoordinates(args.reglist)
			if args.split or args.split_n:
				if not args.split_n or args.split_n > len(region_df.index):
					n = len(region_df.index)
					dist_mode = 'split-list'
				else:
					n = args.split_n
					dist_mode = 'split-list-n'
			else:
				dist_mode = 'list'
			print " " + str(len(region_df.index)) + " regions found"
		elif args.region:
			if len(args.region.split(':')) > 1:
				region_df = pd.DataFrame({'chr': [re.split(':|-', args.region)[0]], 'start': [re.split(':|-', args.region)[1]], 'end': [re.split(':|-', args.region)[2]], 'region': [args.region]})
				dist_mode = 'region'
			else:
				region_df = pd.DataFrame({'chr': [args.region],'start': ['NA'],'end': ['NA'],'region': [args.region]})
				dist_mode = 'chr'
			n = 1
		else:
			region_df = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)]})
			n = 1

		##### get job list from file #####
		if 'jobs' in vars(args).keys() and args.jobs is not None:
			jobs = []
			with open(args.jobs) as f:
				lines = (line.rstrip() for line in f)
				lines = (line for line in lines if line)
				for line in lines:
					if line.find(':') != -1:
						jobs.append(region_df['region'][region_df['region'] == line].index[0])
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
			out_files = FileFxns.GenerateSubFiles(region_df = region_df, f = args.out, dist_mode = dist_mode, n = n)

	##### get user home directory #####
	#home_dir = os.path.expanduser("~")
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
		cmd = args.which.capitalize() + '(out=\'' + args.out + '\''
		for x in ['oxford','dos1','dos2','plink','vcf','region','b','kb','mb','n','chr','shift_mb','shift_kb','shift_b']:
			if x in vars(args).keys() and not vars(args)[x] in [False,None]:
				if type(vars(args)[x]) is str:
					cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
				else:
					cmd = cmd + ',' + x + '=' + str(vars(args)[x])
		cmd = cmd + ')'
		if args.replace:
			for f in [args.out, args.out + '.map.log']:
				try:
					os.remove(f)
				except OSError:
					continue
		else:
			for f in [args.out, args.out + '.map.log']:
				if os.path.exists(f):
					print SystemFxns.Error("1 or more output files already exists (use --replace flag to replace)")
					return
		if args.qsub:
			SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + args.out + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			SystemFxns.Interactive(qsub_wrapper, cmd, args.out + '.' + args.which + '.log')

	elif args.which == 'compile':
		out_files = {}
		if args.which == "compile":
			out_files = FileFxns.FindSubFiles(f = args.data)
		if len(out_files.keys()) > 1:
			existing_files = glob.glob(args.out + '.gz*') + glob.glob(args.out + '.log')
			if len(existing_files) > 0:
				if not args.replace:
					print SystemFxns.Error("1 or more output files or files with similar basename already exists (use --replace flag to replace)")
					return
				else:
					for f in existing_files:
						try:
							os.remove(f)
						except OSError:
							continue
			complete, complete_reg = FileFxns.CheckResults(file_dict=out_files, out=args.out + '.verify')
			if not complete:
				print SystemFxns.Error("results could not be verified")
				return
			out_files = collections.OrderedDict([(x, out_files[x]) for x in out_files.keys() if str(x) in complete_reg])
			if not FileFxns.CompileResults(out_files, args.out):
				print SystemFxns.Error("results could not be compiled")
				return
		else:
			print SystemFxns.Error("no split results to compile")
			return

	elif args.which == 'explore':
		check_files = [args.out + '.explore.log']
		check_files = [args.out + '.top_results']
		check_files = check_files + [args.out + '.qq.tiff'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.qq_strat.tiff'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.qq.eps'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.qq_strat.eps'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.qq.pdf'] if 'qq' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.qq_strat.pdf'] if 'qq_strat' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.mht.tiff'] if 'mht' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.mht.eps'] if 'mht' in vars(args).keys() else check_files
		check_files = check_files + [args.out + '.mht.pdf'] if 'mht' in vars(args).keys() else check_files
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
			print SystemFxns.Error("above files already exist (use --replace flag to replace)")
			return
		cmd = args.which.capitalize() + '(cfg=' + str(config) + ')'
		if args.qsub:
			SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + config['out'] + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			SystemFxns.Interactive(qsub_wrapper, cmd, args.out + '.explore.log')

	elif args.which == 'gc':
		check_files = [args.out + '.gc.log']
		check_files = check_files + [args.out + '.gc.gz']
		check_files = check_files + [args.out + '.gc.gz.tbi']
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
			print SystemFxns.Error("above files already exist (use --replace flag to replace)")
			return
		cmd = args.which.upper() + '(cfg=' + str(config) + ')'
		if args.qsub:
			SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + config['out'] + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			SystemFxns.Interactive(qsub_wrapper, cmd, args.out + '.gc.log')
		#cmd = 'GC(data="' + args.data + '",out="' + args.out + '"'
		#if args.gc:
		#	cmd = cmd + ',gc=' + str(dict(args.gc)) + ')'
		#SystemFxns.Interactive(qsub_wrapper, cmd, args.out + '.gc.log')

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
			print SystemFxns.Error("above files already exist (use --replace flag to replace)")
			return
		cmd = args.which.capitalize() + '(cfg=' + str(config) + ')'
		if args.qsub:
			SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + config['file'].replace('.gz','') + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
		else:
			SystemFxns.Interactive(qsub_wrapper, cmd, args.file.replace('.gz','') + '.annot.log')

	elif args.which in ['model','meta']:
		print "preparing output directories"
		if dist_mode == 'split-list' and n > 1:
			FileFxns.PrepareChrDirs(region_df['region'], directory)
		elif dist_mode == 'split-list-n' and n > 1:
			FileFxns.PrepareListDirs(n, directory)
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
			if i == joblist[0] and (args.jobs is None or i == range(n)[0]):
				config['write_header'] = True
			else:
				config['write_header'] = False
			if i == joblist[len(joblist)-1] and (args.jobs is None or i == range(n)[len(range(n))-1]):
				config['write_eof'] = True
			else:
				config['write_eof'] = False
			if dist_mode in ['split-list', 'region']:
				config['out'] = out_files['%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))]
				if n > 1:
					config['region'] = '%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))
					config['reglist'] = None
			elif dist_mode == 'split-list-n':
				config['out'] = out_files[i]
				rlist = config['out'] + '.reglist'
				region_df.loc[np.array_split(np.array(region_df.index), n)[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'id'])
				config['reglist'] = rlist
			elif dist_mode == 'chr':
				config['out'] = out_files['%s' % (str(region_df['chr'][i]))]
			if args.replace:
				for f in [config['out'], config['out'] + '.' + args.which + '.log',config['out'] + '.gz', config['out'] + '.gz.tbi']:
					try:
						os.remove(f)
					except OSError:
						continue
			else:
				for f in [config['out'],config['out'] + '.log',config['out'] + '.gz', config['out'] + '.gz.tbi']:
					if os.path.exists(f):
						print SystemFxns.Error("1 or more output files already exists (use --replace flag to replace)")
						return
			cmd = args.which.capitalize() + '(cfg=' + str(config) + ')'
			if args.qsub:
				SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + config['out'] + '.' + args.which + '.log ' + qsub_wrapper + ' \"' + cmd + '\"')
			else:
				SystemFxns.Interactive(qsub_wrapper, cmd, config['out'] + '.' + args.which + '.log')
	else:
		print SystemFxns.Error(args.which + " not a module")

	print ''

if __name__ == "__main__":
	main()
	os._exit(0)

