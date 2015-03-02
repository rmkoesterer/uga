import argparse
import os
import sys
from __init__ import __version__

def Parser():
	parser = argparse.ArgumentParser(add_help=False)
	top_parser = argparse.ArgumentParser(parents=[parser])
	subparsers = top_parser.add_subparsers(title='modules', dest='which')

	top_parser.add_argument('--version', 
						action='version', 
						version='Universal Genome Analyst: %(prog)s v' + __version__, 
						help='display version information and exit')

	model_parser = subparsers.add_parser('model', help='marker and locus-based statistical modeling', parents=[parser])
	model_required = model_parser.add_argument_group('required arguments')
	model_required.add_argument('--out', 
						action='store', 
						help='output file name (basename only: do not include path)')
	model_required.add_argument('--pheno', 
						nargs=1, 
						action='store', 
						help='phenotype file (see documentation for required formatting)')
	model_required.add_argument('--fid', 
						nargs=1, 
						action='store', 
						help='column name with family ID')
	model_required.add_argument('--iid', 
						action='store', 
						nargs=1, 
						help='column name with sample ID (The IDs in this column must match the --samples file)')
	model_required.add_argument('--model', 
						action='store', 
						nargs=1, 
						help='model string in the format "phenotype~age+factor(sex)+pc1+pc2+pc3+marker" (see documentation)')
	model_required.add_argument('--method', 
						action='store', 
						nargs=1, 
						choices = ["gee_gaussian","gee_binomial","glm_gaussian","glm_binomial","lme_gaussian","lme_binomial","coxph","efftests","skat_o","famskat_o"], 
						help='analysis method identifier (see documentation  for descriptions)')
	model_parser.add_argument('--cfg', 
						action='store', 
						help='configuration file name (see documentation)')
	model_parser.add_argument('--delimiter', 
						action='store', 
						nargs=1, 
						choices=['tab','space','comma'], 
						help='phenotype file delimiter (default: tab)')
	model_parser.add_argument('--samples', 
						action='store', 
						nargs=1, 
						help='sample file (single column list of IDs in same order as --data file, only required if format != plink)')
	model_parser.add_argument('--focus', 
						action='store', 
						nargs=1, 
						help='comma separated list of variables for which stats will be reported (default: report all stats)')
	model_parser.add_argument('--sig', 
						action='store', 
						type=int, 
						default=5, 
						help='significant digits reported for float type stats (default: SIG=5)')
	model_parser.add_argument('--sex', 
						action='store', 
						nargs=1, 
						help='name of the column containing male/female status (requires --male and --female)')
	model_parser.add_argument('--male', 
						action='store', 
						nargs=1, 
						type=int, 
						default=1, 
						help='code for male (default: MALE=1; requires --sex and --female)')
	model_parser.add_argument('--female', 
						nargs=1, 
						action='store', 
						type=int, 
						default=2, 
						help='code for female (default: FEMALE=2; requires --sex and --male)')
	model_parser.add_argument('--buffer', 
						action='store', 
						type=int, 
						default=100, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: BUFFER=100)')
	model_parser.add_argument('--miss', 
						action='store', 
						type=float, 
						help='threshold value for missingness (ie. MISS=0.95 allows for up to 5%% missingness')
	model_parser.add_argument('--freq', 
						action='store', 
						type=float, 
						help='threshold value for allele frequency (ie. FREQ=0.03 filters out markers with MAF < 0.03')
	model_parser.add_argument('--rsq', 
						action='store', 
						type=float, 
						help='threshold value for imputation quality (ie. RSQ=0.8 filters out markers with r-squared < 0.8)')
	model_parser.add_argument('--hwe', 
						action='store', 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. HWE=1e-6 filters out markers with Hardy Weinberg p-value < 1e-6')
	model_parser.add_argument('--case', 
						nargs=1, 
						action='store', 
						type=int, 
						default=[1], 
						help='code for case in the dependent variable column (requires --ctrl; binomial fxn family only; default: CASE=1)')
	model_parser.add_argument('--ctrl', 
						nargs=1, 
						action='store', 
						type=int, 
						default=[0], 
						help='code for control in the dependent variable column (requires --case; binomial fxn family only; default: CTRL=0)')
	model_parser.add_argument('--corstr', 
						nargs=1, 
						action='store', 
						choices=['exchangeable','independence','ar1','unstructured'], 
						default='exchangeable', 
						help='correlation structure for gee analyses (default: exchangeable)')
	model_parser.add_argument('--nofail', 
						action='store_true', 
						help='exclude filtered/failed analyses from results (if not set, full results are reported with filtered marker stats set to NA')
	model_parser.add_argument('--pedigree', 
						nargs=1, 
						action='store', 
						help='pedigree file')
	model_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	model_parser.add_argument('-q', '--qsub', 
						action='store', 
						help='group ID under which to submit jobs to the queue')
	model_parser.add_argument('--name', 
						action='store', 
						help='job name (only used with --qsub; if not set, --out basename will be used')
	model_parser.add_argument('-d', '--directory', 
						action='store', 
						default=os.getcwd(), 
						help='output directory path (default: current working directory)')
	model_parser.add_argument('--mem', 
						action='store', 
						type=int, 
						default=3, 
						help='amount of ram memory to request for queued job in GB (default: MEM=3')
	model_parser.add_argument('--region-id', 
						action='store', 
						help='add region id to results (for use with --region option)')
	model_parser_split_group1 = model_parser.add_mutually_exclusive_group()
	model_parser_split_group1.add_argument('-r', '--region', 
						action='store', 
						help='region specified in Tabix format (ie. 1:10583-1010582).')
	model_parser_split_group1.add_argument('--region-list', 
						action='store', 
						help='filename for a list of tabix format regions')
	model_parser_split_group2 = model_parser.add_mutually_exclusive_group()
	model_parser_split_group2.add_argument('-s', '--split', 
						action='store_true', 
						help='split region list into a single job for each line in file (requires --region-list)')
	model_parser_split_group2.add_argument('-n', '--split-n', 
						action='store', 
						type=int, 
						help='split region list into SPLIT_N jobs (requires --region-list)')
	model_parser_split_group2.add_argument('--split-chr', 
						action='store_true', 
						help='split into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	model_parser_split_group3 = model_parser.add_mutually_exclusive_group()
	model_parser_split_group3.add_argument('-j', '--job', 
						action='store', 
						type=int, 
						help='run a particular job number (requires --split-n)')
	model_parser_split_group3.add_argument('--job-list', 
						action='store', 
						help='filename for a list of job numbers (requires --split-n)')
	model_parser_split_group4 = model_parser.add_mutually_exclusive_group()
	model_parser_split_group4.add_argument('--oxford', 
						action='store', 
						help='oxford format genotype data file')
	model_parser_split_group4.add_argument('--dos1', 
						action='store', 
						help='dos1 format genotype data file (see documentation)')
	model_parser_split_group4.add_argument('--dos2', 
						action='store', 
						help='dos2 format genotype data file (see documentation)')
	model_parser_split_group4.add_argument('--plink', 
						action='store', 
						help='plink binary format genotype data file (without extension)')
	model_parser_split_group4.add_argument('--vcf', 
						action='store', 
						help='vcf 4.0/4.1 format genotype data file')

	meta_parser = subparsers.add_parser('meta', help='meta-analysis', parents=[parser])
	meta_required = meta_parser.add_argument_group('required arguments')	
	meta_required.add_argument('-c', '--cfg', 
						action='store', 
						required=True, 
						help='configuration file name (see documentation)')
	meta_parser.add_argument('-v', '--vars', 
						action='append', 
						help='declaration of the form A=B, C=D, E=F, ... to replace [A] with B, [C] with D, [E] with F, ... in any line of the cfg file')	
	meta_parser.add_argument('--method', 
						action='store', 
						default='sample_size', 
						choices=['sample_size', 'stderr', 'efftest'], 
						help='meta-analysis method (default: sample_size)')
	meta_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	meta_parser.add_argument('-q', '--qsub', 
						action='store', 
						help='group ID under which to submit jobs to the queue')
	meta_parser.add_argument('--name', 
						action='store', 
						help='job name (only used with --qsub; if not set, --out basename will be used')
	meta_parser.add_argument('-d', '--directory', 
						action='store', 
						default=os.getcwd(), 
						help='output directory path (default: current working directory)')
	meta_parser.add_argument('--mem', 
						action='store', 
						type=int, 
						default=3, 
						help='amount of ram memory to request for queued job in GB (default: MEM=3')
	meta_parser.add_argument('--region-id', 
						action='store', 
						help='add region id to results (for use with --region option)')
	meta_parser_split_group1 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group1.add_argument('-r', '--region', 
						action='store', 
						help='region specified in Tabix format (ie. 1:10583-1010582).')
	meta_parser_split_group1.add_argument('--region-list', 
						action='store', 
						help='filename for a list of tabix format regions')
	meta_parser_split_group2 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group2.add_argument('-s', '--split', 
						action='store_true', 
						help='split region list into 1 job for each line in file (requires --region-list)')
	meta_parser_split_group2.add_argument('-n', '--split-n', 
						action='store', 
						type=int, 
						help='split region list into SPLIT_N jobs (requires --region-list)')
	meta_parser_split_group2.add_argument('--split-chr', 
						action='store_true', 
						help='split jobs into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	meta_parser_split_group3 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group3.add_argument('-j', '--job', 
						action='store', 
						type=int, 
						help='run a particular job number (requires --split-n)')
	meta_parser_split_group3.add_argument('--job-list', 
						action='store', 
						help='filename for a list of job numbers (requires --split-n)')

	map_parser = subparsers.add_parser('map', help='map non-empty regions in genotype data files', parents=[parser])
	map_required = map_parser.add_argument_group('required arguments')
	map_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='output file name')
	map_parser.add_argument('-c','--chr', 
						action='store', 
						type=int,  
						help='chromosome number from 1-26')
	map_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing out file')
	map_parser.add_argument('-q', '--qsub', 
						action='store', 
						help='group ID under which to submit jobs to the queue')
	map_parser.add_argument('--name', 
						action='store', 
						help='job name (only used with --qsub; if not set, --out basename will be used')
	map_parser.add_argument('--split-chr', 
						action='store_true', 
						help='split jobs into 26 chromosomes')
	map_split_group1 = map_parser.add_mutually_exclusive_group()
	map_split_group1.add_argument('--mb', 
						action='store', 
						help='region size (megabase)')
	map_split_group1.add_argument('--kb', 
						action='store', 
						help='region size (kilobase)')
	map_split_group1.add_argument('--b', 
						action='store', 
						help='region size (base)')
	map_split_group1.add_argument('--n', 
						action='store', 
						help='number of markers to be included in each region')
	map_split_group2 = map_parser.add_mutually_exclusive_group()
	map_split_group2.add_argument('--oxford', 
						action='store', 
						help='oxford format genotype data file')
	map_split_group2.add_argument('--dos1', 
						action='store', 
						help='dos1 format genotype data file (see documentation)')
	map_split_group2.add_argument('--dos2', 
						action='store', 
						help='dos2 format genotype data file (see documentation)')
	map_split_group2.add_argument('--plink', 
						action='store', 
						help='plink binary format genotype data file (without extension)')					
	map_split_group2.add_argument('--vcf', 
						action='store', 
						help='vcf 4.0/4.1 format genotype data file')					

	summary_parser = subparsers.add_parser('summary', help='verify, compile, filter and/or plot results files', parents=[parser])
	summary_required = summary_parser.add_argument_group('required arguments')
	summary_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='filename of existing results (basename only: do not include path)')
	summary_parser.add_argument('--out-rename', 
						action='store', 
						help='name for compiled summary files (basename only: do not include path)')
	summary_parser.add_argument('--complete-string', 
						action='store', 
						help='string indicating completeness in the result log file (default: COMPLETE_STRING=\'   ... process complete\')')
	summary_parser.add_argument('--verify', 
						action='store_true', 
						help='verify results before compiling')
	summary_parser.add_argument('--compile', 
						action='store_true', 
						help='compile results')
	summary_parser.add_argument('--qq', 
						action='store_true', 
						help='print qq plot')
	summary_parser.add_argument('--manhattan', 
						action='store_true', 
						help='print manhattan plot')
	summary_parser.add_argument('--gc', 
						action='store_true', 
						help='print plots with genomic inflation corrected p-values')
	summary_parser.add_argument('--chr', 
						action='store', 
						help='column name for chromosome (default: CHR=chr)')
	summary_parser.add_argument('--pos', 
						action='store', 
						help='column name for position (default: POS=pos)')
	summary_parser.add_argument('--p', 
						action='store', 
						help='column name for p-value (default: P=p)')
	summary_parser.add_argument('--rsq', 
						action='store', 
						help='column name for imputation quality (default: RSQ=rsq)')
	summary_parser.add_argument('--freq', 
						action='store', 
						help='column name for allele frequency (default: FREQ=freq)')
	summary_parser.add_argument('--hwe', 
						action='store', 
						help='column name for Hardy Weinberg p-value (default: HWE=hwe)')
	summary_parser.add_argument('--meta-dir', 
						action='store', 
						help='column name meta analysis direction (default: META_DIR=meta.dir)')
	summary_parser.add_argument('--rsq-thresh', 
						action='store', 
						type=float, 
						help='threshold for imputation quality (ie. RSQ_THRESH=0.8 filters out markers with r-squared < 0.8; requires --rsq)')
	summary_parser.add_argument('--freq-thresh', 
						action='store', 
						type=float, 
						help='threshold for allele frequency (ie. FREQ_THRESH=0.03 filters out markers with MAF < 0.03; requires --freq)')
	summary_parser.add_argument('--hwe-thresh', 
						action='store', 
						type=float, 
						help='threshold for Hardy Weinberg p-value (ie. HWE_THRESH=1e-6 filters out markers with Hardy Weinberg p-value < 1e-6; requires --hwe)')
	summary_parser.add_argument('--df-thresh', 
						action='store', 
						type=int, 
						help='threshold for meta analysis degrees of freedom (ie. DF_THRESH=4 filters out markers less than 5 datasets included in the meta analysis; requires --meta-dir)')
	summary_parser.add_argument('--sig', 
						action='store', 
						type=float, 
						help='significant digits reported for float type stats (default: SIG=5)')
	summary_parser.add_argument('--calc-sig', 
						action='store_true', 
						help='calculate significance level from number of filtered markers reported (0.05 / n)')
	summary_parser.add_argument('--ext', 
						action='store', 
						default='tiff', 
						choices=['tiff','eps','pdf'], 
						help='file type extension for plot files')
	summary_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')

	return top_parser
	
def Parse(top_parser):
	args=top_parser.parse_args()
	if args.which in ['model','meta']:
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
	if args.which == 'meta' and args.method == 'efftest' and not args.region_list:
		top_parser.error("missing argument: --region-list required in module meta with --method efftest")
	if args.which == 'map' and not (args.b or args.kb or args.mb or args.n):
		top_parser.error("missing argument: --b, --kb, --mb, or --n required in module map")
	if args.which == 'model' and args.cfg is None and (args.out is None or args.pheno is None or args.fid is None or args.iid is None or args.model is None or args.method is None):
		top_parser.error("missing argument: --out, --pheno, --fid, --iid, --model and --method required in module model without --cfg")
	print ''
	print 'Universal Genome Analyst v' + __version__
	print ''
	print 'active module: ' + args.which
	return args
