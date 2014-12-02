import argparse
import os

def uga_parse():	
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
	parser.add_argument('--name', 
						action='store', 
						help='a job name (only used with --qsub; if not set, --out basename will be used')
	parser.add_argument('-d', '--directory', 
						action='store', 
						default=os.getcwd(), 
						help='an output directory path')
	parser.add_argument('--cpus', 
						action='store', 
						type=int, 
						default=1, 
						help='number of cpus (limited module availability)')
	parser.add_argument('--mem', 
						action='store', 
						type=int, 
						default=3, 
						help='amount of ram memory to request for queued job (in gigabytes)')

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
						choices=['gee_gaussian', 'gee_binomial', 'glm_gaussian', 'glm_binomial', 'lme_gaussian', 'lme_binomial', 'coxph', 'efftests'], 
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
	meta_parser.add_argument('-v', '--vars', 
						action='append', 
						help='a declaration of the form A=B, C=D, E=F, ... to replace [A] with B, [C] with D, [E] with F, ... in any line of the cfg file')	
	meta_parser.add_argument('--method', 
						action='store', 
						default='sample_size', 
						choices=['sample_size', 'stderr', 'efftest'], 
						help='the meta-analysis method')

	annot_parser = subparsers.add_parser('annot', help='annotation', parents=[parser])
	annot_required = annot_parser.add_argument_group('required arguments')

	plot_parser = subparsers.add_parser('plot', help='plot generation', parents=[parser])
	plot_required = plot_parser.add_argument_group('required arguments')

	map_parser = subparsers.add_parser('map', help='map non-empty regions in a file', parents=[parser])
	map_required = map_parser.add_argument_group('required arguments')

	map_impute_parser = subparsers.add_parser('map-impute', help='map non-empty imputation regions overlapping with file', parents=[parser])
	map_impute_required = map_impute_parser.add_argument_group('required arguments')

	verify_parser = subparsers.add_parser('verify', help='verify output files', parents=[parser])
	verify_required = verify_parser.add_argument_group('required arguments')
	verify_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='an output file name (basename only: do not include path)')
	verify_required.add_argument('--verify-out', 
						action='store', 
						required=True, 
						help='an verification output file name (basename only: including path)')
	verify_required.add_argument('--complete-string', 
						action='store', 
						help='a string indicating completeness in the log file (used only with verify module)')

	compile_parser = subparsers.add_parser('compile', help='compile output files', parents=[parser])
	compile_required = compile_parser.add_argument_group('required arguments')
	compile_required.add_argument('--out', 
						action='store', 
						required=True, 
						help='an output file name (basename only: do not include path)')
	compile_required.add_argument('--compile-out', 
						action='store', 
						required=True, 
						help='a compiled output file name (basename only: including path)')
	return top_parser