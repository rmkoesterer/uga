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

import argparse

class AddString(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		if not 'ordered_args' in namespace:
			setattr(namespace, 'ordered_args', [])
		previous = namespace.ordered_args
		previous.append((self.dest, values))
		setattr(namespace, 'ordered_args', previous)

class AddTrue(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		if not 'ordered_args' in namespace:
			setattr(namespace, 'ordered_args', [])
		previous = namespace.ordered_args
		previous.append((self.dest, True))
		setattr(namespace, 'ordered_args', previous)

class AddFalse(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		if not 'ordered_args' in namespace:
			setattr(namespace, 'ordered_args', [])
		previous = namespace.ordered_args
		previous.append((self.dest, False))
		setattr(namespace, 'ordered_args', previous)

def set_args(set_parser):
	set_parser.add_argument('--snpeff', 
						action=AddString, 
						help='set full path to snpEff executable')
	set_parser.add_argument('--snpsift', 
						action=AddString, 
						help='set full path to SnpSift executable')
	set_parser.add_argument('--dbnsfp', 
						action=AddString, 
						help='set full path to dbNSFP database')
	set_parser.add_argument('--locuszoom', 
						action=AddString, 
						help='set full path to locuszoom executable')
	return set_parser

def snv_args(snv_parser):
	snv_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	snv_parser.add_argument('--ped', 
						action=AddString, 
						help='phenotype file')
	snv_parser.add_argument('--cpus', 
						action=AddString, 
						type=int, 
						help='number of cpus')
	snv_parser.add_argument('--fid', 
						action=AddString, 
						help='column name with family ID')
	snv_parser.add_argument('--iid', 
						action=AddString, 
						help='column name with sample ID (The IDs in this column must match the --samples file)')
	snv_parser.add_argument('--matid', 
						action=AddString, 
						help='column name with maternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	snv_parser.add_argument('--patid', 
						action=AddString, 
						help='column name with paternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	snv_parser.add_argument('--all-founders', 
						nargs=0, 
						action=AddTrue, 
						help='use all samples to calculate variant statistics regardless of --matid and --patid column information')
	snv_parser.add_argument('--sep', 
						action=AddString, 
						choices=['tab','space','comma'], 
						help='phenotype file delimiter (default: tab)')
	snv_parser.add_argument('--pheno', 
						action=AddString, 
						help='phenotype')
	snv_parser.add_argument('--covars', 
						action=AddString, 
						help='comma separated list of numeric covariates')
	snv_parser.add_argument('--covars-categorical', 
						action=AddString, 
						help='comma separated list of categorical covariates')
	snv_parser.add_argument('--interact', 
						action=AddString, 
						help='a variable to include in a snv interaction term (if used, the p value for the interaction will be reported: --interact AGE -> SNV*AGE)')
	snv_parser.add_argument('--reverse', 
						nargs=0, 
						action=AddTrue, 
						help='reverse the model so that snv is used as the dependent variable and return a p-value for the phenotype')
	snv_parser.add_argument('--sample', 
						action=AddString, 
						help='sample file (not required for vcf format files)')
	snv_parser.add_argument('--sex', 
						action=AddString, 
						help='name of the column containing male/female status (requires --male and --female)')
	snv_parser.add_argument('--male', 
						action=AddString, 
						type=int, 
						help='code for male (default: 1; requires --sex and --female)')
	snv_parser.add_argument('--female', 
						action=AddString, 
						type=int, 
						help='code for female (default: 2; requires --sex and --male)')
	snv_parser.add_argument('--buffer', 
						action=AddString, 
						type=int, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	snv_parser.add_argument('--miss', 
						action=AddString, 
						type=float, 
						help='threshold value for missingness (ie. 0.95 allows for up to 5%% missingness)')
	snv_parser.add_argument('--maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency (ie. 0.03 filters out markers with maf < 0.03)')
	snv_parser.add_argument('--maxmaf', 
						action=AddString, 
						type=float, 
						help='threshold value for maximum minor allele frequency (ie. 0.01 filters out markers with maf >= 0.01)')
	snv_parser.add_argument('--mac', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele count (ie. 3 filters out markers with mac < 3)')
	snv_parser.add_argument('--rsq', 
						action=AddString, 
						type=float, 
						help='threshold value for imputation quality (ie. 0.8 filters out markers with rsq < 0.8)')
	snv_parser.add_argument('--hwe', 
						action=AddString, 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. 1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	snv_parser.add_argument('--hwe-maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency for which Hardy Weinberg p-value threshold is evaluated (ie. 0.01 only filters markers with maf >= 0.01 that do not pass hwe threshold)')
	snv_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	snv_parser.add_argument('--mb', 
						action=AddString, 
						help='region size in megabases to use for split analyses (default: 1)')
	snv_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	snv_parser.add_argument('--tag', 
						action=AddString, 
						help='tag for individual model')
	snv_parser_split_group1 = snv_parser.add_mutually_exclusive_group()
	snv_parser_split_group1.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format of any size up to an entire chromosome (ie. 1:1-1000000, 21).')
	snv_parser_split_group1.add_argument('--region-file', 
						action=AddString, 
						help='filename for a list of tabix format regions of any size up to an entire chromosome (ie. 1:1-1000000, 21)')
	snv_parser_split_group2 = snv_parser.add_mutually_exclusive_group()
	snv_parser_split_group2.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='split --region-file into an individual job for each line in file (requires --region-file)')
	snv_parser_split_group2.add_argument('--split-n', 
						action=AddString, 
						type=int, 
						help='split --region-file into n individual jobs each with a subset of regions in the file (requires --region-file)')
	snv_parser_split_group2.add_argument('--split-chr', 
						nargs=0, 
						action=AddTrue,  
						help='split data into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	snv_parser_split_group3 = snv_parser.add_mutually_exclusive_group()
	snv_parser_split_group3.add_argument('--job', 
						action=AddString, 
						type=int, 
						help='run a particular job (use with --region-file and --split with value a tabix format region or --split-n with value a number from 1..n)')
	snv_parser_split_group3.add_argument('--jobs', 
						action=AddString, 
						help='filename for a list of jobs to run (use with --region-file and --split with a column of tabix format regions or --split-n with a column of numbers from 1..n)')
	snv_parser.add_argument('--vcf', 
						action=AddString, 
						help='vcf 4.1/4.2 format genotype data file')
	snv_parser.add_argument('--dos', 
						action=AddString, 
						help='single allele dosage format genotype data file')
	snv_parser.add_argument('--case-code', 
						action=AddString, 
						type=int, 
						help='code for case in the dependent variable column (requires --ctrl-code; binomial fxn family only; default: 1)')
	snv_parser.add_argument('--ctrl-code', 
						action=AddString, 
						type=int, 
						help='code for control in the dependent variable column (requires --case-code; binomial fxn family only; default: 0)')
	snv_parser.add_argument('--score', 
						nargs=0, 
						action=AddTrue,
						help='score test (singlesnpMeta)')
	snv_parser.add_argument('--gee', 
						nargs=0, 
						action=AddTrue,
						help='gee test (R geepack geeglm function)')
	snv_parser.add_argument('--corstr', 
						action=AddString, 
						help='correlation structure for gee test (default: exchangeable)')
	snv_parser.add_argument('--glm', 
						nargs=0, 
						action=AddTrue,
						help='glm test')
	snv_parser.add_argument('--lm', 
						nargs=0, 
						action=AddTrue,
						help='lm test')
	snv_parser.add_argument('--meta-sample-size', 
						nargs=2, 
						action=AddString, 
						help='a sample size weighted meta analysis (ie. --meta-sample-size meta m1+m2: meta-analysis with tag "meta", including tags "m1" and "m2"; requires sample size, p-value, and effect size (for direction of effect))')
	snv_parser.add_argument('--meta-stderr', 
						nargs=2, 
						action=AddString, 
						help='a standard error weighted meta analysis (ie. --meta-stderr meta m1+m2: meta-analysis with tag "meta", including tags "m1" and "m2"; requires p-value, effect size, standard error, and consistent units for effect size estimates and standard errors (ie. same trait and model))')
	snv_parser.add_argument('--debug', 
						nargs=0, 
						action=AddTrue, 
						help='enable debug mode (prints debug info to log file)')
	return snv_parser

def snvgroup_args(snvgroup_parser):
	snvgroup_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	snvgroup_parser.add_argument('--ped', 
						action=AddString, 
						help='phenotype file')
	snvgroup_parser.add_argument('--cpus', 
						action=AddString, 
						type=int, 
						help='number of cpus')
	snvgroup_parser.add_argument('--fid', 
						action=AddString, 
						help='column name with family ID')
	snvgroup_parser.add_argument('--iid', 
						action=AddString, 
						help='column name with sample ID (The IDs in this column must match the --samples file)')
	snvgroup_parser.add_argument('--matid', 
						action=AddString, 
						help='column name with maternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	snvgroup_parser.add_argument('--patid', 
						action=AddString, 
						help='column name with paternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	snvgroup_parser.add_argument('--all-founders', 
						nargs=0, 
						action=AddTrue, 
						help='use all samples to calculate variant statistics regardless of --matid and --patid column information')
	snvgroup_parser.add_argument('--sep', 
						action=AddString, 
						choices=['tab','space','comma'], 
						help='phenotype file delimiter (default: tab)')
	snvgroup_parser.add_argument('--sample', 
						action=AddString, 
						help='sample file (not required for vcf format files)')
	snvgroup_parser.add_argument('--pheno', 
						action=AddString, 
						help='phenotype')
	snvgroup_parser.add_argument('--covars', 
						action=AddString, 
						help='comma separated list of numeric covariates')
	snvgroup_parser.add_argument('--covars-categorical', 
						action=AddString, 
						help='comma separated list of categorical covariates')
	snvgroup_parser.add_argument('--sex', 
						action=AddString, 
						help='name of the column containing male/female status (requires --male and --female)')
	snvgroup_parser.add_argument('--male', 
						action=AddString, 
						type=int, 
						help='code for male (default: 1; requires --sex and --female)')
	snvgroup_parser.add_argument('--female', 
						action=AddString, 
						type=int, 
						help='code for female (default: 2; requires --sex and --male)')
	snvgroup_parser.add_argument('--buffer', 
						action=AddString, 
						type=int, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	snvgroup_parser.add_argument('--miss', 
						action=AddString, 
						type=float, 
						help='threshold value for missingness (ie. 0.95 allows for up to 5%% missingness)')
	snvgroup_parser.add_argument('--maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency (ie. 0.03 filters out markers with maf < 0.03)')
	snvgroup_parser.add_argument('--maxmaf', 
						action=AddString, 
						type=float, 
						help='threshold value for maximum minor allele frequency (ie. 0.01 filters out markers with maf >= 0.01)')
	snvgroup_parser.add_argument('--mac', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele count (ie. 3 filters out markers with mac < 3)')
	snvgroup_parser.add_argument('--snvgroup-mac', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele count in snv group (ie. 3 filters out snv groups with mac < 3)')
	snvgroup_parser.add_argument('--rsq', 
						action=AddString, 
						type=float, 
						help='threshold value for imputation quality (ie. 0.8 filters out markers with rsq < 0.8)')
	snvgroup_parser.add_argument('--hwe', 
						action=AddString, 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. 1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	snvgroup_parser.add_argument('--hwe-maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency for which Hardy Weinberg p-value threshold is evaluated (ie. 0.01 only filters markers with maf >= 0.01 that do not pass hwe threshold)')
	snvgroup_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	snvgroup_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	snvgroup_parser.add_argument('--tag', 
						action=AddString, 
						help='tag for individual model')
	snvgroup_parser.add_argument('--snvgroup-map', 
						action=AddString, 
						help='filename for list of tabix formatted single variant regions (ie. 1:100-100) and snvgroup name (or any categorical naming scheme) assignments')
	snvgroup_parser_split_group1 = snvgroup_parser.add_mutually_exclusive_group()
	snvgroup_parser_split_group1.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format of any size up to an entire chromosome (ie. 1:1-1000000, 21).')
	snvgroup_parser_split_group1.add_argument('--region-file', 
						action=AddString, 
						help='filename for a list of tabix format regions of any size up to an entire chromosome (ie. 1:1-1000000, 21)')
	snvgroup_parser_split_group2 = snvgroup_parser.add_mutually_exclusive_group()
	snvgroup_parser_split_group2.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='split --region-file into an individual job for each line in file (requires --region-file)')
	snvgroup_parser_split_group2.add_argument('--split-n', 
						action=AddString, 
						type=int, 
						help='split --region-file into n individual jobs each with a subset of regions in the file (requires --region-file)')
	snvgroup_parser_split_group2.add_argument('--split-chr', 
						nargs=0, 
						action=AddTrue,  
						help='split data into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	snvgroup_parser_split_group3 = snvgroup_parser.add_mutually_exclusive_group()
	snvgroup_parser_split_group3.add_argument('--job', 
						action=AddString, 
						type=int, 
						help='run a particular job (use with --region-file and --split with value a tabix format region or --split-n with value a number from 1..n)')
	snvgroup_parser_split_group3.add_argument('--jobs', 
						action=AddString, 
						help='filename for a list of jobs to run (use with --region-file and --split with a column of tabix format regions or --split-n with a column of numbers from 1..n)')
	snvgroup_parser.add_argument('--vcf', 
						action=AddString, 
						help='vcf 4.1/4.2 format genotype data file')
	snvgroup_parser.add_argument('--dos', 
						action=AddString, 
						help='single allele dosage format genotype data file')
	snvgroup_parser.add_argument('--case-code', 
						action=AddString, 
						type=int, 
						help='code for case in the dependent variable column (requires --ctrl-code; binomial fxn family only; default: 1)')
	snvgroup_parser.add_argument('--ctrl-code', 
						action=AddString, 
						type=int, 
						help='code for control in the dependent variable column (requires --case-code; binomial fxn family only; default: 0)')
	snvgroup_parser.add_argument('--skat', 
						nargs=0, 
						action=AddTrue,
						help='skat test (skatMeta)')
	snvgroup_parser.add_argument('--skato', 
						nargs=0, 
						action=AddTrue,
						help='skato test (skatOMeta)')
	snvgroup_parser.add_argument('--skat-wts', 
						action=AddString, 
						help='skat weights (default: beta weights: function(maf){dbeta(maf,1,25)})')
	snvgroup_parser.add_argument('--burden-wts', 
						action=AddString, 
						help='burden weights (default: T1 weights for skato test: function(maf){maf < 0.01}; constant weight for burden test: 1)')
	snvgroup_parser.add_argument('--skat-method', 
						action=AddString, 
						help='skat method for p-value calculation (default: saddlepoint)')
	snvgroup_parser.add_argument('--skato-rho', 
						action=AddString, 
						help='skato rho parameter (default: seq(0,1,0.1)))')
	snvgroup_parser.add_argument('--timeout', 
						action=AddString, 
						type=int, 
						help='timeout in seconds for model test function (default: 3600)')
	snvgroup_parser.add_argument('--burden', 
						nargs=0, 
						action=AddTrue,
						help='burden test (burdenMeta)')
	snvgroup_parser.add_argument('--mafrange', 
						action=AddString, 
						help='maf range for snv group tests which is different from --maxmaf and --maf since it is to be calculated before and during meta analysis (default: all snvs c(0,0.5))')
	snvgroup_parser.add_argument('--meta', 
						nargs=2, 
						action=AddString, 
						help='a meta analysis (ie. --meta meta m1+m2: meta-analysis with tag "meta", including tags "m1" and "m2")')
	snvgroup_parser.add_argument('--debug', 
						nargs=0, 
						action=AddTrue, 
						help='enable debug mode (prints debug info to log file)')
	return snvgroup_parser

def meta_args(meta_parser):
	meta_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	meta_parser.add_argument('--file', 
						nargs=2,
						action=AddString, 
						help='results file tag and name (ie. --file X file.gz: add file "file.gz" with tag "X"')
	meta_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	meta_parser.add_argument('--cpus', 
						action=AddString, 
						type=int, 
						help='number of cpus')
	meta_parser.add_argument('--mb', 
						action=AddString, 
						help='region size in megabases to use for split analyses (default: 1)')
	meta_parser.add_argument('--buffer', 
						action=AddString, 
						type=int, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	meta_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	meta_parser_split_group1 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group1.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format of any size up to an entire chromosome (ie. 1:1-1000000, 21).')
	meta_parser_split_group1.add_argument('--region-file', 
						action=AddString, 
						help='filename for a list of tabix format regions of any size up to an entire chromosome (ie. 1:1-1000000, 21)')
	meta_parser_split_group2 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group2.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='split --region-file into an individual job for each line in file (requires --region-file)')
	meta_parser_split_group2.add_argument('--split-n', 
						action=AddString, 
						type=int, 
						help='split --region-file into n individual jobs each with a subset of regions in the file (requires --region-file)')
	meta_parser_split_group2.add_argument('--split-chr', 
						nargs=0, 
						action=AddTrue,  
						help='split data into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	meta_parser_split_group3 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group3.add_argument('--job', 
						action=AddString, 
						type=int, 
						help='run a particular job (use with --region-file and --split with value a tabix format region or --split-n with value a number from 1..n)')
	meta_parser_split_group3.add_argument('--jobs', 
						action=AddString, 
						help='filename for a list of jobs to run (use with --region-file and --split with a column of tabix format regions or --split-n with a column of numbers from 1..n)')
	meta_parser.add_argument('--meta-sample-size', 
						nargs=2, 
						action=AddString, 
						help='a sample size weighted meta analysis (ie. --meta-sample-size meta m1+m2: meta-analysis with tag "meta", including tags "m1" and "m2"; requires sample size, p-value, and effect size (for direction of effect))')
	meta_parser.add_argument('--meta-stderr', 
						nargs=2, 
						action=AddString, 
						help='a standard error weighted meta analysis (ie. --meta-stderr meta m1+m2: meta-analysis with tag "meta", including tags "m1" and "m2"; requires p-value, effect size, standard error, and consistent units for effect size estimates and standard errors (ie. same trait and model))')
	meta_parser.add_argument('--debug', 
						nargs=0, 
						action=AddTrue, 
						help='enable debug mode (prints debug info to log file)')

def compile_args(compile_parser):
	compile_required = compile_parser.add_argument_group('required arguments')
	compile_required.add_argument('--dir', 
						action=AddString, 
						required=True, 
						help='base directory name of existing results')
	compile_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	return compile_parser

def resubmit_args(resubmit_parser):
	resubmit_required = resubmit_parser.add_argument_group('required arguments')
	resubmit_required.add_argument('--dir', 
						action=AddString, 
						required=True, 
						help='base directory name of existing results')
	resubmit_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	return resubmit_parser

def snvplot_args(snvplot_parser):
	snvplot_required = snvplot_parser.add_argument_group('required arguments')
	snvplot_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='filename of existing results')
	snvplot_required.add_argument('--pcol', 
						action=AddString, 
						help='a comma separated list of p value column names or a single p value column name (default: p)')
	snvplot_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	snvplot_parser.add_argument('--ext', 
						action=AddString, 
						choices=['tiff','eps','pdf'], 
						help='file type extension for plot files (default: tiff)')
	snvplot_parser.add_argument('--freqcol', 
						action=AddString, 
						help='effect allele frequency column (default: freq)')
	snvplot_parser.add_argument('--crop', 
						action=AddString, 
						type=float, 
						help='crop extreme values at this -log10(p) (default: 10)')
	snvplot_parser.add_argument('--qq', 
						nargs=0, 
						action=AddTrue, 
						help='enable qq plot')
	snvplot_parser.add_argument('--qq-strat', 
						nargs=0, 
						action=AddTrue, 
						help='enable frequency stratified qq plot')
	snvplot_parser.add_argument('--mht', 
						nargs=0, 
						action=AddTrue, 
						help='enable manhattan plot')
	snvplot_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	snvplot_parser.add_argument('--gc', 
						nargs=0, 
						action=AddTrue, 
						help='apply genomic control adjustment for manhattan plots')
	snvplot_parser.add_argument('--color', 
						nargs=0, 
						action=AddTrue, 
						help='plot unique color for each chromosome in manhattan plot (default: two-tone blue color manhattan plot)')
	snvplot_parser.add_argument('--debug', 
						nargs=0, 
						action=AddTrue, 
						help='enable debug mode (prints debug info to log file)')
	return snvplot_parser

def snvgroupplot_args(snvgroupplot_parser):
	snvgroupplot_required = snvgroupplot_parser.add_argument_group('required arguments')
	snvgroupplot_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='filename of existing results')
	snvgroupplot_required.add_argument('--pcol', 
						action=AddString, 
						help='a comma separated list of p value column names or a single p value column name (default: p)')
	snvgroupplot_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	snvgroupplot_parser.add_argument('--ext', 
						action=AddString, 
						choices=['tiff','eps','pdf'], 
						help='file type extension for plot files (default: tiff)')
	snvgroupplot_parser.add_argument('--cmaccol', 
						action=AddString, 
						help='group cumulative minor allele count column (default: cmac)')
	snvgroupplot_parser.add_argument('--crop', 
						action=AddString, 
						type=float, 
						help='crop extreme values at this -log10(p) (default: 10)')
	snvgroupplot_parser.add_argument('--qq', 
						nargs=0, 
						action=AddTrue, 
						help='enable qq plot')
	snvgroupplot_parser.add_argument('--qq-strat', 
						nargs=0, 
						action=AddTrue, 
						help='enable cumulative minor allele count stratified qq plot')
	snvgroupplot_parser.add_argument('--mht', 
						nargs=0, 
						action=AddTrue, 
						help='enable manhattan plot')
	snvgroupplot_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	snvgroupplot_parser.add_argument('--gc', 
						nargs=0, 
						action=AddTrue, 
						help='apply genomic control adjustment for manhattan plots')
	snvgroupplot_parser.add_argument('--color', 
						nargs=0, 
						action=AddTrue, 
						help='plot unique color for each chromosome in manhattan plot (default: two-tone blue color manhattan plot)')
	snvgroupplot_parser.add_argument('--debug', 
						nargs=0, 
						action=AddTrue, 
						help='enable debug mode (prints debug info to log file)')
	return snvgroupplot_parser
