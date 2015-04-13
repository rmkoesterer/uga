import argparse
import os
import sys
from __version__ import version

def Parser():
	parser = argparse.ArgumentParser(add_help=False)
	top_parser = argparse.ArgumentParser(parents=[parser])
	subparsers = top_parser.add_subparsers(title='modules', dest='which')

	top_parser.add_argument('--version', 
						action='version', 
						version='Universal Genome Analyst: %(prog)s v' + version, 
						help='display version information and exit')

	model_parser = subparsers.add_parser('model', help='marker and locus-based statistical modeling', parents=[parser])
	model_parser.add_argument('--out', 
						action='store', 
						help='base output file name (basename only: do not include path)')
	model_parser.add_argument('--pheno', 
						nargs=1, 
						action='store', 
						help='phenotype file (see documentation for required formatting)')
	model_parser.add_argument('--fid', 
						nargs=1, 
						action='store', 
						help='column name with family ID')
	model_parser.add_argument('--iid', 
						action='store', 
						nargs=1, 
						help='column name with sample ID (The IDs in this column must match the --samples file)')
	model_parser.add_argument('--cfg', 
						action='store', 
						help='configuration file name (see documentation)')
	model_parser.add_argument('--pheno-sep', 
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
						help='threshold value for missingness (ie. MISS=0.95 allows for up to 5%% missingness)')
	model_parser.add_argument('--freq', 
						action='store', 
						type=float, 
						help='threshold value for allele frequency (ie. FREQ=0.03 filters out markers with MAF < 0.03)')
	model_parser.add_argument('--rsq', 
						action='store', 
						type=float, 
						help='threshold value for imputation quality (ie. RSQ=0.8 filters out markers with r-squared < 0.8)')
	model_parser.add_argument('--hwe', 
						action='store', 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. HWE=1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
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
						default=['exchangeable'], 
						help='correlation structure for gee analyses (default: exchangeable)')
	model_parser.add_argument('--nofail', 
						action='store_true', 
						help='exclude filtered/failed analyses from results (if not set, full results are reported with filtered marker stats set to NA)')
	model_parser.add_argument('--geeboss-thresh', 
						nargs=1, 
						action='store', 
						type=float, 
						default=1e-7, 
						help='p-value threshold for boss.fit thresh option (see CRAN R boss package documentation)')
	model_parser.add_argument('--pedigree', 
						nargs=1, 
						action='store', 
						help='pedigree file')
	model_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	model_parser.add_argument('-q', '--qsub', 
						action='store', 
						help='string indicating all qsub options to be added to the qsub command')
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
						help='run a particular job (use with --region-list and --split/--split-n)')
	model_parser_split_group3.add_argument('--job-list', 
						action='store', 
						help='filename for a list of jobs (use with --region-list and --split/--split-n)')
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
						help='vcf 4.1/4.2 format genotype data file')
	model_parser_split_group5 = model_parser.add_mutually_exclusive_group()
	model_parser_split_group5.add_argument('--gee-gaussian', 
						action='store', 
						help='model string for gee gaussian analysis')
	model_parser_split_group5.add_argument('--gee-binomial', 
						action='store', 
						help='model string for gee binomial analysis')
	model_parser_split_group5.add_argument('--geeboss-gaussian', 
						action='store', 
						help='model string for gee boss.fit gaussian analysis')
	model_parser_split_group5.add_argument('--geeboss-binomial', 
						action='store', 
						help='model string for gee boss.fit binomial analysis')
	model_parser_split_group5.add_argument('--glm-gaussian', 
						action='store', 
						help='model string for glm gaussian analysis')
	model_parser_split_group5.add_argument('--glm-binomial', 
						action='store', 
						help='model string for glm binomial analysis')
	model_parser_split_group5.add_argument('--lme-gaussian', 
						action='store', 
						help='model string for lme gaussian analysis')
	model_parser_split_group5.add_argument('--lme-binomial', 
						action='store', 
						help='model string for lme binomial analysis')
	model_parser_split_group5.add_argument('--coxph', 
						action='store', 
						help='model string for coxph analysis')
	model_parser_split_group5.add_argument('--efftests', 
						action='store', 
						help='model string for efftests analysis')
	model_parser_split_group5.add_argument('--famskat-o', 
						action='store', 
						help='model string for family based skat-o analysis')
	model_parser_split_group5.add_argument('--skat-o-gaussian', 
						action='store', 
						help='model string for skat-o gaussian analysis')
	model_parser_split_group5.add_argument('--skat-o-binomial', 
						action='store', 
						help='model string for skat-o binomial analysis')
	model_parser_split_group5.add_argument('--famskat', 
						action='store', 
						help='model string for family based skat analysis')
	model_parser_split_group5.add_argument('--skat-gaussian', 
						action='store', 
						help='model string for skat gaussian analysis')
	model_parser_split_group5.add_argument('--skat-binomial', 
						action='store', 
						help='model string for skat binomial analysis')
	model_parser_split_group5.add_argument('--famburden', 
						action='store', 
						help='model string for family based burden analysis')
	model_parser_split_group5.add_argument('--burden-gaussian', 
						action='store', 
						help='model string for burden  gaussian analysis')
	model_parser_split_group5.add_argument('--burden-binomial', 
						action='store', 
						help='model string for burden binomial analysis')

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
						help='string indicating all qsub options to be added to the qsub command')
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
						help='run a particular job (use with --region-list and --split/--split-n)')
	meta_parser_split_group3.add_argument('--job-list', 
						action='store', 
						help='filename for a list of jobs (use with --region-list and --split/--split-n)')

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
						help='vcf 4.1/4.2 format genotype data file')					

	compile_parser = subparsers.add_parser('compile', help='verify and compile plot results files', parents=[parser])
	compile_required = compile_parser.add_argument_group('required arguments')
	compile_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='filename for compiled results (basename only: do not include path)')
	compile_required.add_argument('--data', 
						action='store', 
						required=True, 
						help='base filename of existing results (basename only: do not include path)')
	compile_required.add_argument('--region-list', 
						action='store', 
						required=True, 
						help='filename for a list of tabix format regions')
	compile_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	compile_parser_split_group = compile_parser.add_mutually_exclusive_group()
	compile_parser_split_group.add_argument('-s', '--split', 
						action='store_true', 
						help='split region list into 1 job for each line in file (requires --region-list)')
	compile_parser_split_group.add_argument('-n', '--split-n', 
						action='store', 
						type=int, 
						help='split region list into SPLIT_N jobs (requires --region-list)')
	compile_parser_split_group.add_argument('--split-chr', 
						action='store_true', 
						help='split jobs into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	

	explore_parser = subparsers.add_parser('explore', help='filter and/or plot results files', parents=[parser])
	explore_required = explore_parser.add_argument_group('required arguments')
	explore_required.add_argument('--data', 
						action='store', 
						required=True, 
						help='filename for compiled results')
	explore_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='filename for compiled results (basename only: do not include path)')
	explore_parser.add_argument('--qq', 
						action='store_true', 
						help='print qq plot')
	explore_parser.add_argument('--manhattan', 
						action='store_true', 
						help='print manhattan plot')
	explore_parser.add_argument('--color', 
						action='store_true', 
						help='plot in color')
	explore_parser.add_argument('--gc', 
						action='store_true', 
						help='print manhattan plots with genomic inflation corrected p-values')
	explore_parser.add_argument('--stats-prefix', 
						action='store', 
						help='string indicating prefix for stats to be summarized, not including tag (default: STATS_PREFIX=\'marker\')')
	explore_parser.add_argument('--top-p', 
						action='store', 
						type=float, 
						help='threshold for p-values to export to file (ie. TOP_P=1e-4 prints all markers with p value < 1e-4 to a top results file')
	explore_parser.add_argument('--regional-n', 
						action='store', 
						type=int, 
						help='print regional plots for top REGIONAL_N markers')
	explore_parser.add_argument('--tag', 
						action='store', 
						help='string indicating tag for stats to be summarized, if tag exists')
	explore_parser.add_argument('--unrel', 
						action='store_true', 
						help='filter based on unrel columns')	
	explore_parser.add_argument('--meta-dir', 
						action='store', 
						help='column name meta analysis direction (default: META_DIR=meta.dir)')
	explore_parser.add_argument('--rsq-thresh', 
						action='store', 
						type=float, 
						help='threshold for imputation quality (ie. RSQ_THRESH=0.8 filters out markers with r-squared < 0.8)')
	explore_parser.add_argument('--freq-thresh', 
						action='store', 
						type=float, 
						help='threshold for allele frequency (ie. FREQ_THRESH=0.03 filters out markers with MAF < 0.03)')
	explore_parser.add_argument('--hwe-thresh', 
						action='store', 
						type=float, 
						help='threshold for Hardy Weinberg p-value (ie. HWE_THRESH=1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	explore_parser.add_argument('--callrate-thresh', 
						action='store', 
						type=float, 
						help='threshold for callrate (ie. CALLRATE_THRESH=0.95 filters out markers with callrate < 0.95)')
	explore_parser.add_argument('--effect-thresh', 
						action='store', 
						type=float, 
						help='threshold for effect estimate (ie. EFFECT_THRESH=1.7 filters out markers with effect estimate > 1.7 and < -1.7)')
	explore_parser.add_argument('--stderr-thresh', 
						action='store', 
						type=float, 
						help='threshold for standard error (ie. STDERR_THRESH=5 filters out markers with standard error > 5)')
	explore_parser.add_argument('--or-thresh', 
						action='store', 
						type=float, 
						help='threshold for odds ratio (ie. OR_THRESH=1.3 filters out markers with odds ratio > 1.25 and < 1/1.25 = 0.8)')
	explore_parser.add_argument('--df-thresh', 
						action='store', 
						type=int, 
						help='threshold for meta analysis degrees of freedom (ie. DF_THRESH=4 filters out markers less than 5 datasets included in the meta analysis; requires --meta-dir)')
	explore_parser.add_argument('--sig', 
						action='store', 
						type=float, 
						help='significant digits reported for float type stats (default: SIG=5)')
	explore_parser.add_argument('--ext', 
						action='store', 
						default='tiff', 
						choices=['tiff','eps','pdf'], 
						help='file type extension for plot files')
	explore_parser.add_argument('-o', '--overwrite', 
						action='store_true', 
						help='overwrite existing output files')
	explore_parser.add_argument('--f-dist-dfn', 
						action='store', 
						type=int, 
						help='f-distribution p-value dfn (see python scipy.stats.f.cdf documentation)')
	explore_parser.add_argument('--f-dist-dfd', 
						action='store', 
						type=int, 
						help='f-distribution p-value dfd (see python scipy.stats.f.cdf documentation)')
	explore_parser_split_group = explore_parser.add_mutually_exclusive_group()
	explore_parser_split_group.add_argument('-r', '--region', 
						action='store', 
						help='region specified in Tabix format (ie. 1:10583-1010582).')
	explore_parser_split_group.add_argument('--region-list', 
						action='store', 
						help='filename for a list of tabix format regions')

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
			#if args.split:
			#	assert not args.job, top_parser.error("argument --job: not allowed with argument -s/--split")
			#	assert not args.job_list, top_parser.error("argument --job-list: not allowed with argument -s/--split")
		if args.job_list:
			assert os.path.exists(args.job_list), top_parser.error("argument --job-list: file does not exist")
	if args.which == 'meta' and args.method == 'efftest' and not args.region_list:
		top_parser.error("missing argument: --region-list required in module meta with --method efftest")
	if args.which == 'map' and not (args.b or args.kb or args.mb or args.n):
		top_parser.error("missing argument: --b, --kb, --mb, or --n required in module map")
	if args.which == 'model' and args.cfg is None and (args.out is None or args.pheno is None or args.fid is None or args.iid is None or 
															(args.gee_gaussian is None and args.gee_binomial is None and
															args.glm_gaussian is None and args.glm_binomial is None and
															args.lme_gaussian is None and args.lme_binomial is None and
															args.coxph is None and args.efftests is None and
															args.famskat_o is None and args.skat_o_gaussian is None and args.skat_o_binomial is None and
															args.famskat is None and args.skat_gaussian is None and args.skat_binomial is None and
															args.famburden is None and args.burden_gaussian is None and args.burden_binomial is None and
															args.geeboss_gaussian is None and args.geeboss_binomial is None)):
		top_parser.error("missing argument: --out, --pheno, --fid, --iid, and a model string (ie. --gee-gaussian, etc.) required in module model without --cfg")
	if args.which == 'model' and not (args.famskat_o is None or args.famskat is None or args.famburden is None) and args.pedigree is None:
		top_parser.error("missing argument: --pedigree required for gene based modelling of family data")
	if args.which == 'compile' and not (args.split or args.split_n or args.split_chr):
		top_parser.error("missing argument: --split, --split-n, or --split-chr required for compile module")
	print ''
	print 'Universal Genome Analyst v' + version
	print ''
	print 'active module: ' + args.which
	return args
