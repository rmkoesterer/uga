#!/usr/bin/env python
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
from Parse import Parse,Parser
import multi_key_dict
import FileFxns
import SystemFxns
pd.options.mode.chained_assignment = None

def main(args=None):
	parser=Parser()
	args=Parse(parser)

	##### read cfg file into dictionary #####
	if args.which == 'meta':
		print "reading configuration from file"
		config = FileFxns.LoadMetaCfg(args.cfg, args.which, args.vars)
		args.out = config['out']
	elif  args.which == 'model' and not args.cfg is None:
		print "reading configuration from file"
		config = FileFxns.LoadModelCfg(args.cfg, args.which)
		args.out = config['out']
	elif args.which == 'model':
		print "preparing configuration"
		config = FileFxns.GenerateSingleModelCfg(vars(args))
		args.cfg = 1

	##### define region list #####
	if args.which in ['model','meta','compile']:
		n = 1
		dist_mode = 'full'
		if args.region_list:
			print "generating list of genomic regions ...", 
			region_df = FileFxns.LoadCoordinates(args.region_list)
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
		if 'job_list' in vars(args).keys() and args.job_list is not None:
			jobs = []
			with open(args.job_list) as f:
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
		if 'data' in vars(args).keys() and args.which == 'compile':
			args.data = directory + args.data
		if 'out' in vars(args).keys() and args.which in ['model','meta']:
			args.out = directory + args.out

		##### generate out file names for split jobs or chr/region specific jobs #####
		out_files = {}
		if dist_mode in ['chr','region','split-list','split-list-n']:
			if args.which == "compile":
				out_files = FileFxns.GenerateSubFiles(region_df = region_df, f = args.data, dist_mode = dist_mode, n = n)
			else:
				out_files = FileFxns.GenerateSubFiles(region_df = region_df, f = args.out, dist_mode = dist_mode, n = n)

	##### get user home directory #####
	home_dir = os.path.expanduser("~")

	if args.which == 'map':
		cmd = args.which.capitalize() + '(out=\'' + args.out + '\''
		for x in ['oxford','dos1','dos2','plink','vcf','b','kb','mb','n','chr','shift_mb','shift_kb','shift_b']:
			if x in vars(args).keys() and not vars(args)[x] in [False,None]:
				if type(vars(args)[x]) is str:
					cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
				else:
					cmd = cmd + ',' + x + '=' + str(vars(args)[x])
		cmd = cmd + ')'
		if args.overwrite:
			for f in [args.out, args.out + '.map.log']:
				try:
					os.remove(f)
				except OSError:
					continue
		else:
			for f in [args.out, args.out + '.map.log']:
				if os.path.exists(f):
					print SystemFxns.Error("1 or more output files already exists (use --overwrite flag to replace)")
					return
		SystemFxns.Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.' + args.which + '.log')

	elif args.which == 'compile':
		if len(out_files.keys()) > 1:
			existing_files = glob.glob(args.out + '*')
			if len(existing_files) > 0:
				if not args.overwrite:
					print SystemFxns.Error("1 or more output files or files with similar basename already exists (use --overwrite flag to replace)")
					return
				else:
					for f in existing_files:
						try:
							os.remove(f)
						except OSError:
							continue
			complete, complete_reg = FileFxns.CheckResults(file_dict=out_files, out=args.out + '.verify', overwrite=args.overwrite)
			if not complete:
				print SystemFxns.Error("results could not be verified")
				return
			out_files = collections.OrderedDict([(x, out_files[x]) for x in out_files.keys() if x in complete_reg])
			if not FileFxns.CompileResults(out_files, args.out):
				print SystemFxns.Error("results could not be compiled")
				return
		else:
			print SystemFxns.Error("no split results to compile")
			return

	elif args.which == 'explore':
		check_files = [args.out + '.explore.log']
		check_files = [args.out + '.top_results']
		check_files = check_files + [args.out + '.qq.tiff'] if 'qq' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'tiff' else check_files
		check_files = check_files + [args.out + '.qq_strat.tiff'] if 'qq_strat' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'tiff' else check_files
		check_files = check_files + [args.out + '.qq.eps'] if 'qq' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'eps' else check_files
		check_files = check_files + [args.out + '.qq_strat.eps'] if 'qq_strat' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'eps' else check_files
		check_files = check_files + [args.out + '.qq.pdf'] if 'qq' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'pdf' else check_files
		check_files = check_files + [args.out + '.qq_strat.pdf'] if 'qq_strat' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'pdf' else check_files
		check_files = check_files + [args.out + '.mht.tiff'] if 'mht' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'tiff' else check_files
		check_files = check_files + [args.out + '.mht.eps'] if 'mht' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'eps' else check_files
		check_files = check_files + [args.out + '.mht.pdf'] if 'mht' in vars(args).keys() and 'ext' in vars(args).keys() and args.ext == 'pdf' else check_files
		#check_files = check_files + glob.glob(args.out + '.rgnl.*') if 'region_list' in vars(args).keys() or 'region' in vars(args).keys() or 'regional_n' in vars(args).keys() else check_files
		existing_files = []
		for f in check_files:
			if os.path.exists(f):
				if not args.overwrite:
					existing_files = existing_files + [f]
					print "found file " + str(f)
				else:
					try:
						os.remove(f)
					except OSError:
						continue
		if len(existing_files) > 0:
			print SystemFxns.Error("above files already exist (use --overwrite flag to replace)")
			return
		cmd = 'Explore(data="' + args.data + '",out="' + args.out + '"'
		for x in ['qq','qq_n','qq_strat','mht','color','ext','sig','gc','set_gc','lz_source','lz_build','lz_pop','regional_n','region_list','region_id','region','stat','top_p','tag','unrel','f_dist_dfn','f_dist_dfd','callrate_thresh','rsq_thresh','freq_thresh','hwe_thresh','effect_thresh','stderr_thresh','or_thresh','df_thresh']:
			if x in vars(args).keys() and not vars(args)[x] in [False,None]:
				if type(vars(args)[x]) is str:
					cmd = cmd + ',' + x + '="' + str(vars(args)[x]) + '"'
				else:
					cmd = cmd + ',' + x + '=' + str(vars(args)[x])
		cmd = cmd + ')'
		SystemFxns.Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.explore.log')

	elif args.which == 'gc':
		check_files = [args.out + '.gc.log']
		check_files = check_files + [args.out + '.gc.gz']
		check_files = check_files + [args.out + '.gc.gz.tbi']
		existing_files = []
		for f in check_files:
			if os.path.exists(f):
				if not args.overwrite:
					existing_files = existing_files + [f]
					print "found file " + str(f)
				else:
					try:
						os.remove(f)
					except OSError:
						continue
		if len(existing_files) > 0:
			print SystemFxns.Error("above files already exist (use --overwrite flag to replace)")
			return
		cmd = 'GC(data="' + args.data + '",out="' + args.out + '"'
		if args.gc:
			cmd = cmd + ',gc=' + str(dict(args.gc)) + ')'
		SystemFxns.Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.gc.log')

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
		elif not args.job_list is None:
			joblist.extend(jobs)
		else:
			joblist.extend(range(n))
		for i in joblist:
			if dist_mode in ['split-list', 'region']:
				config['out'] = out_files['%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))]
				if n > 1:
					config['region'] = '%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))
					config['region_list'] = None
			elif dist_mode == 'split-list-n':
				config['out'] = out_files[i]
				rlist = config['out'] + '.regions'
				region_df.loc[np.array_split(np.array(region_df.index), n)[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'reg_id'])
				config['region_list'] = rlist
			elif dist_mode == 'chr':
				config['out'] = out_files['%s' % (str(region_df['chr'][i]))]
			if args.overwrite:
				for f in [config['out'], config['out'] + '.' + args.which + '.log',config['out'] + '.gz', config['out'] + '.gz.tbi']:
					try:
						os.remove(f)
					except OSError:
						continue
			else:
				for f in [config['out'],config['out'] + '.log',config['out'] + '.gz', config['out'] + '.gz.tbi']:
					if os.path.exists(f):
						print SystemFxns.Error("1 or more output files already exists (use --overwrite flag to replace)")
						return
			if args.which == 'model':
					cmd = args.which.capitalize() + '(cfg=' + str(config) + ')'
			elif args.which == 'meta':
				cmd = args.which.capitalize() + '(cfg=' + str(config)
				for x in ['region_list', 'region']:
					if x in vars(args).keys() and not vars(args)[x] in [False,None]:
						if type(vars(args)[x]) is str:
							cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
						else:
							cmd = cmd + ',' + x + '=' + str(vars(args)[x])
				cmd = cmd + ')'
			if args.qsub:
				SystemFxns.Qsub('qsub ' + args.qsub + ' -o ' + config['out'] + '.' + args.which + '.log ' + home_dir + '/.uga_wrapper.py \"' + cmd + '\"')
			else:
				SystemFxns.Interactive(home_dir + '/.uga_wrapper.py', cmd, config['out'] + '.' + args.which + '.log')
	else:
		print SystemFxns.Error(args.which + " not a module")

	print ''

if __name__ == "__main__":
	main()
	os._exit(0)

