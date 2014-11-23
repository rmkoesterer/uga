#!/usr/bin/env python

import os
import gzip
import sys
import argparse
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
import submit
from __init__ import __version__


def main(args=None):
	parser = argparse.ArgumentParser('parent', add_help=False)
	parser.add_argument('--version', 
						action='version', 
						version='%(prog)s v0.6 pre-release')
	parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	parser.add_argument('-q', '--qsub', 
						action='store', 
						help='a group ID under which to submit jobs to the queue')
	parser.add_argument('-v', '--vars', 
						action='append', 
						help='a declaration of the form A=B, C=D, E=F, ... to replace [A] with B, [C] with D, [E] with F, ... in any line of the cfg file')
	parser.add_argument('-d', '--directory', 
						action='store', 
						default=os.getcwd(), 
						help='an output directory path')
	parser.add_argument('--cpus', 
						action='store', 
						type=int, 
						default=1, 
						help='number of cpus (limited module availability)')
	parser.add_argument('--verify', 
						action='store_true', 
						help='verify results for any split-job process')
	parser.add_argument('--compile', 
						action='store', 
						help='a filename for compiled output')

	parser_split_group1 = parser.add_mutually_exclusive_group()
	parser_split_group1.add_argument('-r', '--region', 
						action='store', 
						help='a region specified in tabix format (ie. 1:10583-1010582).')
	parser_split_group1.add_argument('--region-list', 
						action='store', 
						help='a filename for a list of tabix format regions')
	parser_split_group2 = parser.add_mutually_exclusive_group()
	parser_split_group2.add_argument('-s', '--split', 
						action='store_true', 
						help='split region list entirely into separate jobs (requires --region-list)')
	parser_split_group2.add_argument('-n', '--split-n', 
						action='store', 
						type=int, 
						help='split region list into n separate jobs (requires --region-list)')
	parser_split_group2.add_argument('--chr', 
						action='store', 
						type=int, 
						help='chromosome number (1-26)')
	parser_split_group2.add_argument('--split-chr', 
						action='store_true', 
						help='split by chromosome (will generate up to 26 separate jobs depending on chromosome coverage)')
	parser_split_group3 = parser.add_mutually_exclusive_group()
	parser_split_group3.add_argument('-j', '--job', 
						action='store', 
						type=int, 
						help='run a particular job number (requires --split-n)')
	parser_split_group3.add_argument('--job-list', 
						action='store', 
						help='a filename for a list of job numbers (requires --split-n)')

	top_parser = argparse.ArgumentParser(parents=[parser])
	subparsers = top_parser.add_subparsers(title='modules', dest='which')

	analyze_parser = subparsers.add_parser('analyze', help='marker and gene-based analysis', parents=[parser])
	analyze_required = analyze_parser.add_argument_group('required arguments')
	analyze_required.add_argument('--data', 
						action='store', 
						required=True, 
						help='a genomic data file')
	analyze_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='an output file name (basename only: do not include path)')
	analyze_required.add_argument('--samples', 
						action='store', 
						required=True, 
						help='a sample file (single column list of IDs in order of data)')
	analyze_required.add_argument('--pheno', 
						action='store', 
						required=True, 
						help='a tab delimited phenotype file')
	analyze_required.add_argument('--model', 
						action='store', 
						required=True, 
						help='a comma separated list of models in the format "phenotype~age+factor(sex)+pc1+pc2+pc3+marker"')
	analyze_required.add_argument('--fid', 
						action='store', 
						required=True, 
						help='the column name with family ID')
	analyze_required.add_argument('--iid', 
						action='store', 
						required=True, 
						help='the column name with sample ID')
	analyze_required.add_argument('--method', 
						action='store', 
						required=True, 
						choices=['gee_gaussian', 'gee_binomial', 'glm_gaussian', 'glm_binomial', 'lme_gaussian', 'lme_binomial', 'coxph'], 
						help='the analysis method')
	analyze_parser.add_argument('--focus', 
						action='store', 
						help='a comma separated list of variables for which stats will be output (default: marker)')
	analyze_parser.add_argument('--sig', 
						action='store', 
						type=int, 
						default=5, 
						help='significant digits to include in output (default: 5)')
	analyze_parser.add_argument('--sex', 
						action='store', 
						help='name of the column containing male/female status (requires MALE_CODE and female)')
	analyze_parser.add_argument('--male', 
						action='store', 
						type=int, 
						default=1, 
						help='the code for a male in sex (requires --sex and --female)')
	analyze_parser.add_argument('--female', 
						action='store', 
						type=int, 
						default=2, 
						help='the code for a male in sex (requires --sex and --female)')
	analyze_parser.add_argument('--buffer', 
						action='store', 
						type=int, 
						default=100, 
						help='a value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	analyze_parser.add_argument('--miss', 
						action='store', 
						type=float, 
						help='a threshold value for missingness')
	analyze_parser.add_argument('--freq', 
						action='store', 
						type=float, 
						help='a threshold value for allele frequency')
	analyze_parser.add_argument('--rsq', 
						action='store', 
						type=float, 
						help='a threshold value for r-squared (imputation quality)')
	analyze_parser.add_argument('--hwe', 
						action='store', 
						type=float, 
						help='a threshold value for Hardy Weinberg p-value')
	analyze_parser.add_argument('--case', 
						action='store', 
						type=int, 
						default=1, 
						help='the code for a case in the dependent variable column (requires ctrl; binomial fxn family only; default: 1)')
	analyze_parser.add_argument('--ctrl', 
						action='store', 
						type=int, 
						default=0, 
						help='the code for a control in the dependent variable column (requires case; binomial fxn family only; default: 0)')
	analyze_parser.add_argument('--nofail', 
						action='store_true', 
						help='exclude filtered/failed analyses from results')

	meta_parser = subparsers.add_parser('meta', help='meta-analysis', parents=[parser])
	meta_required = meta_parser.add_argument_group('required arguments')	
	meta_required.add_argument('-c', '--cfg', 
						action='store', 
						required=True, 
						help='a configuration file name')

	annot_parser = subparsers.add_parser('annot', help='annotation', parents=[parser])
	annot_required = annot_parser.add_argument_group('required arguments')

	plot_parser = subparsers.add_parser('plot', help='plot generation', parents=[parser])
	plot_required = plot_parser.add_argument_group('required arguments')

	map_parser = subparsers.add_parser('map', help='map non-empty regions in a file', parents=[parser])
	map_required = map_parser.add_argument_group('required arguments')

	map_impute_parser = subparsers.add_parser('map-impute', help='map non-empty imputation regions overlapping with file', parents=[parser])
	map_impute_required = map_impute_parser.add_argument_group('required arguments')

	efftests_parser = subparsers.add_parser('efftests', help='determine effective number of tests (Li and Ji method)', parents=[parser])
	efftests_required = efftests_parser.add_argument_group('required arguments')

	args=top_parser.parse_args()

	if args.region:
		assert not args.split, parser.error("argument -s/--split: not allowed with argument --region")
		assert not args.split_n, parser.error("argument -n/--split-n: not allowed with argument --region")
		assert not args.job, parser.error("argument -j/--job: not allowed with argument --region")
		assert not args.job_list, parser.error("argument --job-list: not allowed with argument --region")
		assert not args.chr, parser.error("argument --chr: not allowed with argument --region")
		assert not args.split_chr, parser.error("argument --split-chr: not allowed with argument --region")
	if args.region_list:
		assert os.path.exists(args.region_list), parser.error("argument --region-list: file does not exist")
		if args.split:
			assert not args.job, parser.error("argument --job: not allowed with argument -s/--split")
			assert not args.job_list, parser.error("argument --job-list: not allowed with argument -s/--split")
		if args.split_n:
			if args.job_list:
				assert os.path.exists(args.job_list), parser.error("argument --job-list: file does not exist")

	print "   " + " ".join(sys.argv)
	opts=vars(args)
	maxopt = max([len(a) for a in opts.keys() if not a == 'which'])
	print "      {0:>{1}}".format('module', maxopt) + ": " + str(args.which)
	for arg in [a for a in opts.keys() if not a == 'which']:
		if opts[arg]:
			print "      {0:>{1}}".format(str(arg), maxopt) + ": " + str(opts[arg])

	##### define script library path #####
	scripts = os.path.dirname(os.path.abspath(__file__))

	##### read cfg file into dictionary #####
	if args.which == 'me':
		print "   ... reading configuration from file"
		config = Cfg(opts['cfg'], opts['which'], opts['vars']).Load()

	cmd = ''
	for x in ['out', 'data', 'samples', 'pheno', 'model', 'fid', 'iid', 'method', 'focus', 'sig', 'region_list', 'region', 'sex', 'male', 'female', 'buffer', 'miss', 'freq', 'rsq', 'hwe', 'case', 'ctrl', 'nofail']:
		if opts[x]:
			cmd = cmd + ' --' + x + ' ' + str(opts[x])
	
	Interactive('submit.py','\"' + cmd + '\"', opts['out'] + '.log')
	"""
	##### define region list #####
	n = 1
	if region_list != '':
		print "   ... reading list of regions from file ...", 
		regionlist = Coordinates(region_list).Load()
		print " " + str(len(regionlist.index)) + " regions found"
		if split:
			if split_n == '':
				n = len(regionlist.index)
				dist_mode = 'split-list'
			else:
				n = split_n
				dist_mode = 'split-list-n'
		else:
			dist_mode = 'list'
	if region != '':
		regionlist = pd.DataFrame({'chr': [re.split(':|-', region)[0]], 'start': [re.split(':|-', region)[1]], 'end': [re.split(':|-', region)[2]], 'region': [region]})
		n = 1
		dist_mode = 'region'
	
		
	##### get job list from file #####
	if job_list != '':
		jobs = []
		with open(job_list) as f:
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
	directory = os.path.dirname(config['out']) if directory == '' else directory
	directory = directory + '/' if directory != '' else directory
	if n > 100:
		if dist_mode == 'split-list':
			directory = directory + 'chr[CHR]/'
		elif dist_mode == 'split-list-n':
			directory = directory + 'list[LIST]/'
	config['out'] = directory + config['out']
	
	out_files = GenerateSubFiles(regionlist = regionlist, f = config['out'], dist_mode = dist_mode, n = split_n)

	if check:
		if module in ['ma', 'ind']:
			complete_string = '   ... process complete' if complete_string == '' else complete_string
		elif module == 'me':
			complete_string = '   ... meta analysis complete' if complete_string == '' else complete_string
		if dist_mode in ['split-list', 'split-list-n']:
			CheckResults(out_files, cfg + '.check', cpus, complete_string, overwrite)
		else:
			usage(Error("--check option available only if both --list and --split or --split-n are used"))
	elif compile != '':
		if dist_mode in ['split-list', 'split-list-n']:
			CompileResults(out_files,  compile, cpus, overwrite)
		else:
			usage(Error("single file results, nothing to compile"))
	else:
		if module in ['ma', 'ind']:
			print "   ... preparing output directories"
			if dist_mode == 'split-list' and n > 100:
				PrepareChrDirs(regionlist['region'], directory)
			elif dist_mode == 'split-list-n' and n > 100:
				PrepareListDirs(n, directory)
			print "   ... submitting analysis jobs\n" if qsub != '' else "   ... starting analysis\n"
			joblist = []
			if job != '':
				joblist.append(job)
			elif job_list != '':
				joblist.extend(jobs)
			else:
				joblist.extend(range(n))
			for i in joblist:
				config_temp = config.copy()
				if dist_mode in ['split-list', 'region']:
					config_temp['out'] = out_files['%d:%d-%d' % (regionlist['chr'][i], regionlist['start'][i], regionlist['end'][i])]
				if dist_mode == 'split-list-n':
					#config_temp['out'] = config_temp['out'].replace('[LIST]', str(int(np.floor(i/100.0) + 99*np.floor(i/100.0))) + '-' + str(int(np.floor(i/100.0) + 99*np.floor(i/100.0) + 99))) + '.list' + str(i)
					config_temp['out'] = out_files[i]
					rlist = config_temp['out'] + '.regions'
					regionlist.loc[regionlist_idx[i]].to_csv(rlist, header=False, index=False, sep='\t', columns=['region', 'reg_id'])
				if overwrite:
					RemoveExistingFiles(config_temp['out'], module)
				else:
					CheckExistingFiles(config_temp['out'], module)
				if qsub != "":
					cmd_arg = scripts + "/Analysis.py " + ' '.join("--%s %s" % (key, val) if key != 'model' else "--%s '%s'" % (key, val) for (key, val) in config_temp.items())
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
	print 'CAFGAP v' + __version__
	main()
