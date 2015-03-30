#!/usr/bin/env python
import os
import gzip
import sys
import string
import math
import time
import subprocess
import os.path
import numpy as np
import pandas as pd
import re
import glob
from FileFxns import *
from SystemFxns import *
from Parse import Parse,Parser

def main(args=None):
	parser=Parser()
	args=Parse(parser)

	##### read cfg file into dictionary #####
	if args.which == 'meta':
		print "reading configuration from file"
		config = Cfg(args.cfg, args.which, args.vars).Load()
		args.out = config['out']
	elif  args.which == 'model' and not args.cfg is None:
		print "reading configuration from file"
		config = Cfg(args.cfg, args.which).Load()
		args.out = config['out']

	##### define region list #####
	if args.which in ['model','meta','compile']:
		n = 1
		dist_mode = 'full'
		if args.region_list:
			print "generating list of genomic regions ...", 
			region_df = Coordinates(args.region_list).Load()
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
				out_files = GenerateSubFiles(region_df = region_df, f = args.data, dist_mode = dist_mode, n = n)
			else:
				out_files = GenerateSubFiles(region_df = region_df, f = args.out, dist_mode = dist_mode, n = n)

	##### get user home directory #####
	home_dir = os.path.expanduser("~")

	if args.which == 'map':
		if args.split_chr:
			for i in range(26):
				cmd = args.which.capitalize() + '(out=\'' + args.out + '.chr' + str(i+1) + '\',chr=' + str(i+1)
				for x in ['oxford','dos1','dos2','plink','vcf','b','kb','mb','n']:
					if x in vars(args).keys() and not vars(args)[x] in [False,None]:
						if type(vars(args)[x]) is str:
							cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
						else:
							cmd = cmd + ',' + x + '=' + str(vars(args)[x])
				cmd = cmd + ')'
				if args.overwrite:
					for f in [args.out + '.chr' + str(i+1), args.out + '.chr' + str(i+1) + '.map.log']:
						try:
							os.remove(f)
						except OSError:
							continue
				else:
					for f in [args.out + '.chr' + str(i+1), args.out + '.chr' + str(i+1) + '.map.log']:
						if os.path.exists(f):
							print Error("1 or more output files already exists (use --overwrite flag to replace)")
							return
				Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.' + args.which + '.log')
		else:
			cmd = args.which.capitalize() + '(out=\'' + args.out + '\''
			for x in ['oxford','dos1','dos2','plink','vcf','b','kb','mb','n','chr']:
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
						print Error("1 or more output files already exists (use --overwrite flag to replace)")
						return
			Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.' + args.which + '.log')

	elif args.which == 'compile':
		if len(out_files.keys()) > 1:
			existing_files = glob.glob(args.out + '*')
			if len(existing_files) > 0:
				if not args.overwrite:
					print Error("1 or more output files or files with similar basename already exists (use --overwrite flag to replace)")
					return
				else:
					for f in existing_files:
						try:
							os.remove(f)
						except OSError:
							continue		
			if not CheckResults(file_dict=out_files, out=args.out + '.verify', overwrite=args.overwrite):
				print Error("results could not be verified")
				return
			if not CompileResults(out_files, args.out, args.overwrite):
				print Error("results could not be compiled")
				return
		else:
			print Error("no split results to compile")
			return

	elif args.which == 'explore':
		existing_files = glob.glob(args.out + '*')
		if len(existing_files) > 0:
			if not args.overwrite:
				print Error("1 or more output files already exists (use --overwrite flag to replace)")
				return
			else:	
				for f in existing_files:
					try:
						os.remove(f)
					except OSError:
						continue
		cmd = 'Explore(data="' + args.data + '",out="' + args.out + '"'
		for x in ['qq','manhattan','color','ext','sig','gc','top_p','regional_n','stats_prefix','tag','unrel','f_dist_dfn','f_dist_dfd','callrate_thresh','rsq_thresh','freq_thresh','hwe_thresh','effect_thresh','stderr_thresh','or_thresh','df_thresh']:
			if x in vars(args).keys() and not vars(args)[x] in [False,None]:
				if type(vars(args)[x]) is str:
					cmd = cmd + ',' + x + '="' + str(vars(args)[x]) + '"'
				else:
					cmd = cmd + ',' + x + '=' + str(vars(args)[x])
		cmd = cmd + ')'
		Interactive(home_dir + '/.uga_wrapper.py', cmd, args.out + '.' + args.which + '.log')

	elif args.which in ['model','meta']:
		print "preparing output directories"
		if dist_mode == 'split-list' and n > 1:
			PrepareChrDirs(region_df['region'], directory)
		elif dist_mode == 'split-list-n' and n > 1:
			PrepareListDirs(n, directory)
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
				out = out_files['%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))]
				if n > 1:
					vars(args)['region'] = '%s:%s-%s' % (str(region_df['chr'][i]), str(region_df['start'][i]), str(region_df['end'][i]))
					vars(args)['region_list'] = None
			elif dist_mode == 'split-list-n':
				out = out_files[i]
				rlist = out + '.regions'
				region_df.loc[np.array_split(np.array(region_df.index), n)[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'reg_id'])
				vars(args)['region_list'] = rlist
			elif dist_mode == 'chr':
				out = out_files['%s' % (str(region_df['chr'][i]))]
			else:
				out = args.out
			if args.overwrite:
				for f in [out, out + '.' + args.which + '.log',out + '.gz', out + '.gz.tbi']:
					try:
						os.remove(f)
					except OSError:
						continue
			else:
				for f in [out, out + '.log',out + '.gz', out + '.gz.tbi']:
					if os.path.exists(f):
						print Error("1 or more output files already exists (use --overwrite flag to replace)")
						return
			if args.which == 'model':
				if not args.cfg is None:
					config['out'] = out
					cmd = args.which.capitalize() + '(cfg=' + str(config)
					for x in ['region_list', 'region', 'region_id']:
						if x in vars(args).keys() and not vars(args)[x] in [False,None]:
							if type(vars(args)[x]) is str:
								cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
							else:
								cmd = cmd + ',' + x + '=' + str(vars(args)[x])
					cmd = cmd + ')'
				else:
					cmd = args.which.capitalize() + '(out=\'' + out + '\''
					for x in ['oxford','dos1','dos2','plink','vcf','samples','pheno','fid','iid','focus','sig','region_list','gee_gaussian','gee_binomial',
								'glm_gaussian','glm_binomial','lme_gaussian','lme_binomial','coxph','efftests','famskat_o','skat_o_gaussian','skat_o_binomial',
								'famskat','skat_gaussian','skat_binomial','famburden','burden_gaussian','burden_binomial',
								'region','region_id','sex','male','female','buffer','corstr','miss','freq','rsq','hwe','case','ctrl','nofail',
								'pedigree','pheno_sep']:
						if x in vars(args).keys() and not str(vars(args)[x]) in ['False','None']:
							if x in ['oxford','dos1','dos2','plink','vcf']:
								cmd = cmd + ',data=[\'' + str(vars(args)[x]) + '\'],format=[\'' + x + '\']'
							elif x in ['gee_gaussian','gee_binomial','glm_gaussian','glm_binomial','lme_gaussian','lme_binomial','coxph','efftests',
											'famskat_o','skat_o_gaussian','skat_o_binomial','famskat','skat_gaussian','skat_binomial','famburden','burden_gaussian','burden_binomial']:
								cmd = cmd + ',model=[\'' + str(vars(args)[x]) + '\'],method=[\'' + x + '\']'
							elif type(vars(args)[x]) is str:
								cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
							else:
								cmd = cmd + ',' + x + '=' + str(vars(args)[x])
					cmd = cmd + ')'
			elif args.which == 'meta':
				config['out'] = out
				cmd = args.which.capitalize() + '(cfg=' + str(config)
				for x in ['region_list', 'region']:
					if x in vars(args).keys() and not vars(args)[x] in [False,None]:
						if type(vars(args)[x]) is str:
							cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
						else:
							cmd = cmd + ',' + x + '=' + str(vars(args)[x])
				cmd = cmd + ')'
			if args.qsub:
				Qsub('qsub ' + args.qsub + ' -o ' + out + '.' + args.which + '.log ' + home_dir + '/.uga_wrapper.py \"' + cmd + '\"')
			else:
				Interactive(home_dir + '/.uga_wrapper.py', cmd, out + '.' + args.which + '.log')

	else:
		print Error(args.which + " not a module")

	print ''

if __name__ == "__main__":
	main()
