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
from Messages import *
from File import *
from Cfg import *
from Coordinates import *
from Process import *
from Parse import *
from __init__ import __version__

def main(args=None):
	parser = Args(__version__)
	parser.AddAll()
	args=parser.Initiate()

	##### read cfg file into dictionary #####
	if args.which == 'meta':
		print "   ... reading configuration from file"
		config = Cfg(args.cfg, args.which, args.vars).Load()
		args.out = config['out']

	##### define region list #####
	n = 1
	dist_mode = 'full'
	if args.region_list:
		print "   ... reading list of regions from file ...", 
		regionlist = Coordinates(args.region_list).Load()
		print " " + str(len(regionlist.index)) + " regions found"
		if args.split or args.split_n:
			if not args.split_n:
				n = len(regionlist.index)
				dist_mode = 'split-list'
			else:
				n = args.split_n
				dist_mode = 'split-list-n'
		else:
			dist_mode = 'list'
	if args.region:
		regionlist = pd.DataFrame({'chr': [re.split(':|-', args.region)[0]], 'start': [re.split(':|-', args.region)[1]], 'end': [re.split(':|-', args.region)[2]], 'region': [args.region]})
		n = 1
		dist_mode = 'region'
	
	##### get job list from file #####
	if args.job_list:
		jobs = []
		with open(args.job_list) as f:
			lines = (line.rstrip() for line in f)
			lines = (line for line in lines if line)
			for line in lines:
				jobs.append(int(line))
		print "   ... " + str(len(jobs)) + " jobs read from job list file"

	if dist_mode == 'split-list':
		regionlist_idx = regionlist.index
	if dist_mode == 'split-list-n':
		regionlist_idx = np.array_split(np.array(regionlist.index), n)

	##### define output directory #####
	directory = os.path.dirname(args.out) if not args.directory else args.directory
	directory = directory + '/' if args.directory else directory
	if n > 100:
		if dist_mode == 'split-list':
			directory = directory + 'chr[CHR]/'
		elif dist_mode == 'split-list-n':
			directory = directory + 'list[LIST]/'
	args.out = directory + args.out
	
	out_files = {}
	if dist_mode in ['region','split-list','split-list-n']:
		out_files = GenerateSubFiles(regionlist = regionlist, f = args.out, dist_mode = dist_mode, n = n)

	##### define script library path #####
	script_path = '/'.join(sys.argv[0].split('/')[:-1])
	
	if args.which == 'verify':
		complete_string = '   ... process complete' if not args.complete_string else args.complete_string
		if dist_mode in ['split-list', 'split-list-n']:
			CheckResults(out_files, args.verify_out + '.verify', args.cpus, complete_string, args.overwrite)
		else:
			print Error("--check option available only if both --list and --split or --split-n are used")
			parser.print_help()
	elif args.which == 'compile':
		if dist_mode in ['split-list', 'split-list-n'] and len(out_files.keys()) > 1:
			CompileResults(out_files, args.compile_out, args.overwrite)
		else:
			print Error("single file results, nothing to compile")
			parser.print_help()
	else:
		if not os.path.exists(args.directory):
			try:
				os.mkdir(args.directory)
			except OSError:
				print Error("unable to create output directory")
				parser.print_help()
		print "   ... preparing output directories"
		if dist_mode == 'split-list' and n > 100:
			PrepareChrDirs(regionlist['region'], directory)
		elif dist_mode == 'split-list-n' and n > 100:
			PrepareListDirs(n, directory)
		print "   ... submitting analysis jobs\n" if args.qsub else "   ... starting analysis\n"
		joblist = []
		if args.job:
			joblist.append(args.job)
		elif args.job_list:
			joblist.extend(jobs)
		else:
			joblist.extend(range(n))
		name = '.'.join(os.path.basename(args.out).split('.')[:-1]) if not args.name else args.name
		for i in joblist:
			if dist_mode in ['split-list', 'region']:
				out = out_files['%s:%s-%s' % (str(regionlist['chr'][i]), str(regionlist['start'][i]), str(regionlist['end'][i]))]
				vars(args)['region'] = '%s:%s-%s' % (str(regionlist['chr'][i]), str(regionlist['start'][i]), str(regionlist['end'][i]))
				vars(args)['region_list'] = None
			elif dist_mode == 'split-list-n':
				out = out_files[i]
				rlist = out + '.regions'
				regionlist.loc[regionlist_idx[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'reg_id'])
				vars(args)['region_list'] = rlist
			else:
				out = args.out
			if args.overwrite:
				RemoveExistingFiles(out, args.which)
			else:
				CheckExistingFiles(out, args.which)
			if args.which == 'analyze':
				cmd = args.which.capitalize() + '(out="' + out + '"'
				for x in ['data', 'samples', 'pheno', 'model', 'fid', 'iid', 'method', 'focus', 'sig', 'region_list', 'region', 'sex', 'male', 'female', 'buffer', 'miss', 'freq', 'rsq', 'hwe', 'case', 'ctrl', 'nofail']:
					if x in vars(args).keys() and not vars(args)[x] in [False,None]:
						if type(vars(args)[x]) is str:
							cmd = cmd + ',' + x + '="' + str(vars(args)[x]) + '"'
						else:
							cmd = cmd + ',' + x + '=' + str(vars(args)[x])
				cmd = cmd + ',mem=' + str(args.mem) + ')'
			elif args.which == 'meta':
				config['out'] = out
				cmd = args.which.capitalize() + '(cfg=' + str(config)
				for x in ['region_list', 'region', 'method']:
					if x in vars(args).keys() and not vars(args)[x] in [False,None]:
						if type(vars(args)[x]) is str:
							cmd = cmd + ',' + x + '=\'' + str(vars(args)[x]) + '\''
						else:
							cmd = cmd + ',' + x + '=' + str(vars(args)[x])
				cmd = cmd + ',mem=' + str(args.mem) + ')'
			if args.qsub:
				Qsub('qsub -P ' + args.qsub + ' -l mem_free=' + str(args.mem) + 'g -N ' + name + ' -o ' + out + '.log ' + script_path + '/../../bin/submit.py --qsub ' + args.qsub + ' --cmd \"' + cmd + '\"')
			else:
				Interactive(script_path + '/../../bin/submit.py', cmd, out + '.log')
		"""
				if dist_mode in ['list', 'split-list-n']:
					cmd_arg = cmd_arg + " --list " + rlist
				if dist_mode == 'region':
					cmd_arg = cmd_arg + " --region " + regionlist['region'][i]
				cmd="qsub -P " + qsub + " -l mem_free=1g -N " + cfg.split("/")[-1].replace('.', '_') + " -o " + config_temp['out'] + ".log" + " " + scripts + "/Submit.sh \"" + cmd_arg + "\""
				Qsub(cmd)
			else:
				cmd = scripts + "/Submit.sh"
				cmd_arg = scripts + "/Analysis.py " + ' '.join("--%s %s" % (key, val) if key != 'model' else "--%s \"%s\"" % (key, val) for (key, val) in config_temp.items())
				if dist_mode in ['list', 'split-list-n']:
					cmd_arg = cmd_arg + " --list " + rlist
				if dist_mode == 'region':
					cmd_arg = cmd_arg + " --region " + regionlist['region'][i]
				Interactive(cmd, cmd_arg, config_temp['out'] + '.log')
	elif module == 'me':
		print "   ... preparing output directories"
		if dist_mode in ['split-full', 'split-list']:
			PrepareDirs(regionlist['region'], directory)
		print "   ... submitting meta analysis jobs\n" if qsub != '' else "   ... starting meta analysis\n"
		for i in range(n):
			config_temp = config.copy()
			if dist_mode in ['split-list', 'region']:
				config_temp['out'] = config_temp['out'].replace('[CHR]', regionlist['chr'][i]) + '.chr' + regionlist['chr'][i] + 'bp' + regionlist['start'][i] + '-' + regionlist['end'][i]
			if overwrite:
				RemoveExistingFiles(config_temp['out'], module)
			else:
				CheckExistingFiles(config_temp['out'], module)
			cmd_arg = scripts + "/Meta.py --cfg " + cfg
			if dist_mode == 'list':
				cmd_arg = cmd_arg + " --list " + list
			if dist_mode == 'region':
				cmd_arg = cmd_arg + " --region " + regionlist['region'][i]
			cmd_arg = cmd_arg + " --directory " + directory if directory != '' else cmd_arg
			cmd_arg = cmd_arg + " --vars " + vars if vars != '' else cmd_arg
			if qsub != "":
				cmd="qsub -P " + qsub + " -l mem_free=1g -N " + cfg.split("/")[-1].replace('.', '_') + " -o " + config_temp['out'] + ".log" + " " + scripts + "/Submit.sh \"" + cmd_arg + "\""
				Qsub(cmd)
			else:
				cmd = scripts + "/Submit.sh"
				Interactive(cmd, cmd_arg, config_temp['out'] + '.log')
	"""
	print ''

if __name__ == "__main__":
	main()
