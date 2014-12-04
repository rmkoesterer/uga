import argparse
import os
import sys

class Args(object):

	def __init__(self, version):
		self.parser = argparse.ArgumentParser('parent', add_help=False)
		self.parser.add_argument('--version', 
						action='version', 
						version='Universal Genome Analyst: %(prog)s v' + version, 
						help='display version information and exit')
		self.parser.add_argument('-o', '--overwrite', 
							action='store_true', 
							help='overwrite existing output files')
		self.parser.add_argument('-q', '--qsub', 
							action='store', 
							help='a group ID under which to submit jobs to the queue')
		self.parser.add_argument('--name', 
							action='store', 
							help='a job name (only used with --qsub; if not set, --out basename will be used')
		self.parser.add_argument('-d', '--directory', 
							action='store', 
							default=os.getcwd(), 
							help='an output directory path')
		self.parser.add_argument('--cpus', 
							action='store', 
							type=int, 
							default=1, 
							help='number of cpus (limited module availability)')
		self.parser.add_argument('--mem', 
							action='store', 
							type=int, 
							default=3, 
							help='amount of ram memory to request for queued job (in gigabytes)')

		parser_split_group1 = self.parser.add_mutually_exclusive_group()
		parser_split_group1.add_argument('-r', '--region', 
							action='store', 
							help='a region specified in tabix format (ie. 1:10583-1010582).')
		parser_split_group1.add_argument('--region-list', 
							action='store', 
							help='a filename for a list of tabix format regions')
		parser_split_group2 = self.parser.add_mutually_exclusive_group()
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
		parser_split_group3 = self.parser.add_mutually_exclusive_group()
		parser_split_group3.add_argument('-j', '--job', 
							action='store', 
							type=int, 
							help='run a particular job number (requires --split-n)')
		parser_split_group3.add_argument('--job-list', 
							action='store', 
							help='a filename for a list of job numbers (requires --split-n)')
		self.top_parser = argparse.ArgumentParser(parents=[self.parser])
		self.subparsers = self.top_parser.add_subparsers(title='modules', dest='which')

	def AddModel(self):
		model_parser = self.subparsers.add_parser('model', help='marker and gene-based models', parents=[self.parser])
		model_required = model_parser.add_argument_group('required arguments')
		model_required.add_argument('--data', 
							action='store', 
							required=True, 
							help='a genomic data file')
		model_required.add_argument('--out', 
							action='store', 
							required=True, 
							help='an output file name (basename only: do not include path)')
		model_required.add_argument('--samples', 
							action='store', 
							required=True, 
							help='a sample file (single column list of IDs in order of data)')
		model_required.add_argument('--pheno', 
							action='store', 
							required=True, 
							help='a tab delimited phenotype file')
		model_required.add_argument('--model', 
							action='store', 
							required=True, 
							help='a comma separated list of models in the format "phenotype~age+factor(sex)+pc1+pc2+pc3+marker"')
		model_required.add_argument('--fid', 
							action='store', 
							required=True, 
							help='the column name with family ID')
		model_required.add_argument('--iid', 
							action='store', 
							required=True, 
							help='the column name with sample ID')
		model_required.add_argument('--method', 
							action='store', 
							required=True, 
							choices=['gee_gaussian', 'gee_binomial', 'glm_gaussian', 'glm_binomial', 'lme_gaussian', 'lme_binomial', 'coxph', 'efftests'], 
							help='the analysis method')
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
							help='name of the column containing male/female status (requires MALE_CODE and female)')
		model_parser.add_argument('--male', 
							action='store', 
							type=int, 
							default=1, 
							help='the code for a male in sex (requires --sex and --female)')
		model_parser.add_argument('--female', 
							action='store', 
							type=int, 
							default=2, 
							help='the code for a male in sex (requires --sex and --female)')
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
							help='the code for a case in the dependent variable column (requires ctrl; binomial fxn family only; default: 1)')
		model_parser.add_argument('--ctrl', 
							action='store', 
							type=int, 
							default=0, 
							help='the code for a control in the dependent variable column (requires case; binomial fxn family only; default: 0)')
		model_parser.add_argument('--nofail', 
							action='store_true', 
							help='exclude filtered/failed analyses from results')

	def AddMeta(self):
		meta_parser = self.subparsers.add_parser('meta', help='meta-analysis', parents=[self.parser])
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

	def AddAnnot(self):
		annot_parser = self.subparsers.add_parser('annot', help='annotation', parents=[self.parser])
		annot_required = annot_parser.add_argument_group('required arguments')

	def AddPlot(self):
		plot_parser = self.subparsers.add_parser('plot', help='plot generation', parents=[self.parser])
		plot_required = plot_parser.add_argument_group('required arguments')

	def AddMap(self):
		map_parser = self.subparsers.add_parser('map', help='map non-empty regions in a file', parents=[self.parser])
		map_required = map_parser.add_argument_group('required arguments')

	def AddMapImpute(self):
		map_impute_parser = self.subparsers.add_parser('map-impute', help='map non-empty imputation regions overlapping with file', parents=[self.parser])
		map_impute_required = map_impute_parser.add_argument_group('required arguments')

	def AddVerify(self):
		verify_parser = self.subparsers.add_parser('verify', help='verify output files', parents=[self.parser])
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

	def AddCompile(self):
		compile_parser = self.subparsers.add_parser('compile', help='compile output files', parents=[self.parser])
		compile_required = compile_parser.add_argument_group('required arguments')
		compile_required.add_argument('--out', 
							action='store', 
							required=True, 
							help='an output file name (basename only: do not include path)')
		compile_required.add_argument('--compile-out', 
							action='store', 
							required=True, 
							help='a compiled output file name (basename only: including path)')

	def AddAll(self):
		self.AddModel()
		self.AddMeta()
		self.AddAnnot()
		self.AddPlot()
		self.AddMap()
		self.AddMapImpute()
		self.AddVerify()
		self.AddCompile()

	def Parse(self):
		return self.top_parser.parse_args()

	def Initiate(self):
		args=self.Parse()
		if args.region:
			assert not args.split, self.top_parser.error("argument -s/--split: not allowed with argument --region")
			assert not args.split_n, self.top_parser.error("argument -n/--split-n: not allowed with argument --region")
			assert not args.job, self.top_parser.error("argument -j/--job: not allowed with argument --region")
			assert not args.job_list, self.top_parser.error("argument --job-list: not allowed with argument --region")
			assert not args.chr, self.top_parser.error("argument --chr: not allowed with argument --region")
			assert not args.split_chr, self.top_parser.error("argument --split-chr: not allowed with argument --region")
		if args.region_list:
			assert os.path.exists(args.region_list), self.top_parser.error("argument --region-list: file does not exist")
			if args.split:
				assert not args.job, self.top_parser.error("argument --job: not allowed with argument -s/--split")
				assert not args.job_list, self.top_parser.error("argument --job-list: not allowed with argument -s/--split")
			if args.split_n:
				if args.job_list:
					assert os.path.exists(args.job_list), self.top_parser.error("argument --job-list: file does not exist")
		print "   " + " ".join(sys.argv)
		return args
	