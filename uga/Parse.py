import argparse
import os
import sys
from __init__ import __version__

def Parser():
	parser = argparse.ArgumentParser(add_help=False)
	top_parser = argparse.ArgumentParser(parents=[parser])
	subparsers = top_parser.add_subparsers(title='modules', dest='which')

	parser.add_argument('--version', 
						action='version', 
						version='Universal Genome Analyst: %(prog)s v' + __version__, 
						help='display version information and exit')
	parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	parser.add_argument('-q', '--qsub', 
						action='store', 
						help='a group ID under which to submit jobs to the queue')
	parser.add_argument('--name', 
						action='store', 
						help='a job name (only used with --qsub; if not set, --out basename will be used')
	parser.add_argument('-d', '--directory', 
						action='store', 
						default=os.getcwd(), 
						help='an output directory path')
	parser.add_argument('--cpu', 
						action='store', 
						type=int, 
						default=1, 
						help='number of cpus (limited module availability; default: 1)')
	parser.add_argument('--mem', 
						action='store', 
						type=int, 
						default=3, 
						help='amount of ram memory to request for queued job (in gigabytes)')
	parser.add_argument('--region-id', 
						action='store', 
						help='a region id (for use with --region option)')
	parser.add_argument('--pedigree', 
						action='store', 
						help='a pedigree filename (for use with --method famskat_o)')

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

	model_parser = subparsers.add_parser('model', help='marker and gene-based models', parents=[parser])
	model_required = model_parser.add_argument_group('required arguments')
	model_required.add_argument('--data', 
						action='store', 
						required=True, 
						help='a genomic data file')
	model_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='an output file name (basename only: do not include path)')
	model_required.add_argument('--pheno', 
						action='store', 
						required=True, 
						help='a tab delimited phenotype file')
	model_required.add_argument('--fid', 
						action='store', 
						required=True, 
						help='the column name with family ID')
	model_required.add_argument('--iid', 
						action='store', 
						required=True, 
						help='the column name with sample ID')
	model_required.add_argument('--model', 
						action='store', 
						required=True, 
						help='a model string in the format "phenotype~age+factor(sex)+pc1+pc2+pc3+marker"')
	model_required.add_argument('--method', 
						action='store', 
						required=True, 
						choices = ["gee_gaussian","gee_binomial","glm_gaussian","glm_binomial","lme_gaussian","lme_binomial","coxph","efftests","skat_o","famskat_o"], 
						help='an analysis method')
	model_parser.add_argument('--samples', 
						action='store', 
						help='a sample file (single column list of IDs in same order as --data file, only for format != plink)')
	model_parser.add_argument('--focus', 
						action='store', 
						help='a comma separated list of variables for which stats will be output (default: marker)')
	model_parser.add_argument('--sig', 
						action='store', 
						type=int, 
						default=5, 
						help='significant digits to include in output (default: 5)')
	model_parser.add_argument('--sex', 
						action='store', 
						help='name of the column containing male/female status (requires --male and --female)')
	model_parser.add_argument('--male', 
						action='store', 
						type=int, 
						default=1, 
						help='the code for a male(default: 1; requires --sex and --female)')
	model_parser.add_argument('--female', 
						action='store', 
						type=int, 
						default=2, 
						help='the code for a female(default: 2; requires --sex and --male)')
	model_parser.add_argument('--buffer', 
						action='store', 
						type=int, 
						default=100, 
						help='a value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	model_parser.add_argument('--miss', 
						action='store', 
						type=float, 
						help='a threshold value for missingness')
	model_parser.add_argument('--freq', 
						action='store', 
						type=float, 
						help='a threshold value for allele frequency')
	model_parser.add_argument('--rsq', 
						action='store', 
						type=float, 
						help='a threshold value for r-squared (imputation quality)')
	model_parser.add_argument('--hwe', 
						action='store', 
						type=float, 
						help='a threshold value for Hardy Weinberg p-value')
	model_parser.add_argument('--case', 
						action='store', 
						type=int, 
						default=1, 
						help='the code for a case in the dependent variable column (requires --ctrl; binomial fxn family only; default: 1)')
	model_parser.add_argument('--ctrl', 
						action='store', 
						type=int, 
						default=0, 
						help='the code for a control in the dependent variable column (requires --case; binomial fxn family only; default: 0)')
	model_parser.add_argument('--format', 
						action='store', 
						default='oxford', 
						choices=['oxford', 'dos1','dos2','plink'], 
						help='the format of the data file, either oxford (3 genotype probabilities), dos1 (single allele dosage with oxford order), dos2 (single allele dosage with [chr pos markername a1 a2] format), or Plink binary (default: oxford)')
	model_parser.add_argument('--nofail', 
						action='store_true', 
						help='exclude filtered/failed analyses from results')

	meta_parser = subparsers.add_parser('meta', help='meta-analysis', parents=[parser])
	meta_required = meta_parser.add_argument_group('required arguments')	
	meta_required.add_argument('-c', '--cfg', 
						action='store', 
						required=True, 
						help='a configuration file name')
	meta_parser.add_argument('-v', '--vars', 
						action='append', 
						help='a declaration of the form A=B, C=D, E=F, ... to replace [A] with B, [C] with D, [E] with F, ... in any line of the cfg file')	
	meta_parser.add_argument('--method', 
						action='store', 
						default='sample_size', 
						choices=['sample_size', 'stderr', 'efftest'], 
						help='the meta-analysis method (default: sample_size)')

	annot_parser = subparsers.add_parser('annot', help='annotation (inactive)', parents=[parser])
	annot_required = annot_parser.add_argument_group('required arguments')

	map_parser = subparsers.add_parser('map', help='map non-empty regions in a file', parents=[parser])
	map_required = map_parser.add_argument_group('required arguments')
	map_required.add_argument('--file', 
						action='store', 
						required=True, 
						help='a file name')
	map_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='an output file name (basename only: do not include path)')
	map_split_group1 = map_parser.add_mutually_exclusive_group()
	map_split_group1.add_argument('--mb', 
						action='store', 
						help='region size in megabases')
	map_split_group1.add_argument('--kb', 
						action='store', 
						help='region size in kilobases')
	map_split_group1.add_argument('--b', 
						action='store', 
						help='region size in kilobases')

	map_impute_parser = subparsers.add_parser('map-impute', help='map non-empty regions in a file overlapping imputation reference file (inactive)', parents=[parser])
	map_impute_required = map_impute_parser.add_argument_group('required arguments')

	summary_parser = subparsers.add_parser('summary', help='verify, compile and/or plot results files', parents=[parser])
	summary_required = summary_parser.add_argument_group('required arguments')
	summary_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='results filename (basename only: do not include path)')
	summary_parser.add_argument('--out-rename', 
						action='store', 
						help='a name for compiled summary files (basename only: do not include path)')
	summary_parser.add_argument('--complete-string', 
						action='store', 
						help='a string indicating completeness in the result log file')
	summary_parser.add_argument('--verify', 
						action='store_true', 
						help='verify results')
	summary_parser.add_argument('--compile', 
						action='store_true', 
						help='compile results')
	summary_parser.add_argument('--qq', 
						action='store_true', 
						help='print qq plot')
	summary_parser.add_argument('--manhattan', 
						action='store_true', 
						help='print manhattan plot')
	summary_parser.add_argument('--chr', 
						action='store', 
						help='a column name for chromosome (default: chr)')
	summary_parser.add_argument('--pos', 
						action='store', 
						help='a column name for position (default: pos)')
	summary_parser.add_argument('--p', 
						action='store', 
						help='a column name for p-value (default: p)')
	summary_parser.add_argument('--rsq', 
						action='store', 
						help='a column name for imputation quality (default: rsq)')
	summary_parser.add_argument('--freq', 
						action='store', 
						help='a column name for allele frequency (default: freq)')
	summary_parser.add_argument('--hwe', 
						action='store', 
						help='a column name for Hardy Weinberg p-value (default: hwe)')
	summary_parser.add_argument('--meta-dir', 
						action='store', 
						help='a column name meta analysis direction (default: meta.dir)')
	summary_parser.add_argument('--rsq-thresh', 
						action='store', 
						type=float, 
						help='a threshold for imputation quality (requires --rsq)')
	summary_parser.add_argument('--freq-thresh', 
						action='store', 
						type=float, 
						help='a threshold for allele frequency (requires --freq)')
	summary_parser.add_argument('--hwe-thresh', 
						action='store', 
						type=float, 
						help='a threshold for Hardy Weinberg p-value (requires --hwe)')
	summary_parser.add_argument('--df-thresh', 
						action='store', 
						type=int, 
						help='a threshold for meta analysis degrees of freedom (requires --meta-dir)')
	summary_parser.add_argument('--sig', 
						action='store', 
						type=float, 
						help='a significance level for filtered results')
	summary_parser.add_argument('--calc-sig', 
						action='store_true', 
						help='calculate significance level from number of filtered markers (0.05 / n)')
	summary_parser.add_argument('--ext', 
						action='store', 
						default='tiff', 
						choices=['tiff','eps','pdf'], 
						help='a file type extension')

	
	return top_parser

def Parse(top_parser):
	args=top_parser.parse_args()
	if args.region:
		assert not args.split, top_parser.error("argument -s/--split: not allowed with argument --region")
		assert not args.split_n, top_parser.error("argument -n/--split-n: not allowed with argument --region")
		assert not args.job, top_parser.error("argument -j/--job: not allowed with argument --region")
		assert not args.job_list, top_parser.error("argument --job-list: not allowed with argument --region")
		assert not args.split_chr, top_parser.error("argument --split-chr: not allowed with argument --region")
	if args.region_list:
		assert os.path.exists(args.region_list), top_parser.error("argument --region-list: file does not exist")
		if args.split:
			assert not args.job, top_parser.error("argument --job: not allowed with argument -s/--split")
			assert not args.job_list, top_parser.error("argument --job-list: not allowed with argument -s/--split")
		if args.split_n:
			if args.job_list:
				assert os.path.exists(args.job_list), top_parser.error("argument --job-list: file does not exist")
	if args.which == 'meta' and not args.region and not args.region_list:
		top_parser.error("missing argument: --region or --region-list required in module meta")
	if args.which == 'meta' and args.method == 'efftest' and not args.region_list:
		top_parser.error("missing argument: --region-list required in module meta with --method efftest")
	print ''
	print 'Universal Genome Analyst v' + __version__
	print ''
	print '   ' + ' '.join(sys.argv)
	print ''
	return args
