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

def SetArgs(set_parser):
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

def SnvArgs(snv_parser):
	snv_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	snv_parser.add_argument('--pheno', 
						action=AddString, 
						help='phenotype file')
	snv_parser.add_argument('--cpus', 
						action=AddString, 
						type=int, 
						help='number of cpus')
	#snv_parser.add_argument('--snv-list', 
	#					action=AddString, 
	#					help='variant list file')
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
	snv_parser.add_argument('--sample', 
						action=AddString, 
						help='sample file (not required for Plink format files)')
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
	#snv_parser.add_argument('--id', 
	#					action=AddString, 
	#					help='add id column to indicate region / locus')
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
	snv_parser.add_argument('--oxford', 
						action=AddString, 
						help='oxford format genotype data file')
	snv_parser.add_argument('--dos1', 
						action=AddString, 
						help='dos1 format genotype data file')
	snv_parser.add_argument('--dos2', 
						action=AddString, 
						help='dos2 format genotype data file')
	snv_parser.add_argument('--plink', 
						action=AddString, 
						help='Plink binary format genotype data file (without extension)')
	snv_parser.add_argument('--vcf', 
						action=AddString, 
						help='vcf 4.1/4.2 format genotype data file')
	snv_parser.add_argument('--case-code', 
						action=AddString, 
						type=int, 
						help='code for case in the dependent variable column (requires --ctrl-code; binomial fxn family only; default: 1)')
	snv_parser.add_argument('--ctrl-code', 
						action=AddString, 
						type=int, 
						help='code for control in the dependent variable column (requires --case-code; binomial fxn family only; default: 0)')
	snv_parser.add_argument('--bssmeta', 
						action=AddString, 
						help='model string for bssmeta (binomial singlesnpMeta)')
	snv_parser.add_argument('--meta', 
						action=AddString, 
						help='a meta analysis string')
	return snv_parser

def GeneArgs(gene_parser):
	gene_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	gene_parser.add_argument('--pheno', 
						action=AddString, 
						help='phenotype file')
	gene_parser.add_argument('--cpus', 
						action=AddString, 
						type=int, 
						help='number of cpus')
	gene_parser.add_argument('--fid', 
						action=AddString, 
						help='column name with family ID')
	gene_parser.add_argument('--iid', 
						action=AddString, 
						help='column name with sample ID (The IDs in this column must match the --samples file)')
	gene_parser.add_argument('--matid', 
						action=AddString, 
						help='column name with maternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	gene_parser.add_argument('--patid', 
						action=AddString, 
						help='column name with paternal ID (used to determine founders for marker stats and for analyses requiring a pedigree)')
	gene_parser.add_argument('--all-founders', 
						nargs=0, 
						action=AddTrue, 
						help='use all samples to calculate variant statistics regardless of --matid and --patid column information')
	gene_parser.add_argument('--sep', 
						action=AddString, 
						choices=['tab','space','comma'], 
						help='phenotype file delimiter (default: tab)')
	gene_parser.add_argument('--sample', 
						action=AddString, 
						help='sample file (not required for Plink format files)')
	gene_parser.add_argument('--sex', 
						action=AddString, 
						help='name of the column containing male/female status (requires --male and --female)')
	gene_parser.add_argument('--male', 
						action=AddString, 
						type=int, 
						help='code for male (default: 1; requires --sex and --female)')
	gene_parser.add_argument('--female', 
						action=AddString, 
						type=int, 
						help='code for female (default: 2; requires --sex and --male)')
	gene_parser.add_argument('--buffer', 
						action=AddString, 
						type=int, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	gene_parser.add_argument('--miss', 
						action=AddString, 
						type=float, 
						help='threshold value for missingness (ie. 0.95 allows for up to 5%% missingness)')
	gene_parser.add_argument('--maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency (ie. 0.03 filters out markers with maf < 0.03)')
	gene_parser.add_argument('--maxmaf', 
						action=AddString, 
						type=float, 
						help='threshold value for maximum minor allele frequency (ie. 0.01 filters out markers with maf >= 0.01)')
	gene_parser.add_argument('--mac', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele count (ie. 3 filters out markers with mac < 3)')
	gene_parser.add_argument('--rsq', 
						action=AddString, 
						type=float, 
						help='threshold value for imputation quality (ie. 0.8 filters out markers with rsq < 0.8)')
	gene_parser.add_argument('--hwe', 
						action=AddString, 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. 1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	gene_parser.add_argument('--hwe-maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency for which Hardy Weinberg p-value threshold is evaluated (ie. 0.01 only filters markers with maf >= 0.01 that do not pass hwe threshold)')
	gene_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	gene_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	gene_parser.add_argument('--tag', 
						action=AddString, 
						help='tag for individual model')
	gene_parser_split_group1 = gene_parser.add_mutually_exclusive_group()
	gene_parser_split_group1.add_argument('--gene-map', 
						action=AddString, 
						help='filename for list of tabix formatted single variant regions (ie. 1:100-100) and gene name (or any categorical naming scheme) assignments')
	gene_parser_split_group2 = gene_parser.add_mutually_exclusive_group()
	gene_parser_split_group2.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='split --region-file into an individual job for each line in file (requires --region-file)')
	gene_parser_split_group2.add_argument('--split-n', 
						action=AddString, 
						type=int, 
						help='split --region-file into n individual jobs each with a subset of regions in the file (requires --region-file)')
	gene_parser_split_group2.add_argument('--split-chr', 
						nargs=0, 
						action=AddTrue,  
						help='split data into chromosomes (will generate up to 26 separate jobs depending on chromosome coverage)')
	gene_parser_split_group3 = gene_parser.add_mutually_exclusive_group()
	gene_parser_split_group3.add_argument('--job', 
						action=AddString, 
						type=int, 
						help='run a particular job (use with --region-file and --split with value a tabix format region or --split-n with value a number from 1..n)')
	gene_parser_split_group3.add_argument('--jobs', 
						action=AddString, 
						help='filename for a list of jobs to run (use with --region-file and --split with a column of tabix format regions or --split-n with a column of numbers from 1..n)')
	gene_parser.add_argument('--oxford', 
						action=AddString, 
						help='oxford format genotype data file')
	gene_parser.add_argument('--dos1', 
						action=AddString, 
						help='dos1 format genotype data file')
	gene_parser.add_argument('--dos2', 
						action=AddString, 
						help='dos2 format genotype data file')
	gene_parser.add_argument('--plink', 
						action=AddString, 
						help='Plink binary format genotype data file (without extension)')
	gene_parser.add_argument('--vcf', 
						action=AddString, 
						help='vcf 4.1/4.2 format genotype data file')
	gene_parser.add_argument('--case-code', 
						action=AddString, 
						type=int, 
						help='code for case in the dependent variable column (requires --ctrl-code; binomial fxn family only; default: 1)')
	gene_parser.add_argument('--ctrl-code', 
						action=AddString, 
						type=int, 
						help='code for control in the dependent variable column (requires --case-code; binomial fxn family only; default: 0)')
	gene_parser.add_argument('--bskato', 
						action=AddString, 
						help='model string for bskato (binomial skatOMeta)')
	gene_parser.add_argument('--meta', 
						action=AddString, 
						help='a meta analysis string')
	return gene_parser

def MetaArgs(meta_parser):
	meta_required = meta_parser.add_argument_group('required arguments')
	meta_parser.add_argument('--out', 
						action=AddString, 
						help='output file basename (do not include path or extension)')
	meta_parser.add_argument('--marker-col', 
						action=AddString, 
						help='variant(marker) name column')
	meta_parser.add_argument('--freq-col', 
						action=AddString, 
						help='effect allele frequency column')
	meta_parser.add_argument('--rsq-col', 
						action=AddString, 
						help='imputation quality column')
	meta_parser.add_argument('--hwe-col', 
						action=AddString, 
						help='hardy-weinberg p-value column')
	meta_parser.add_argument('--effect-col', 
						action=AddString, 
						help='effect column')
	meta_parser.add_argument('--stderr-col', 
						action=AddString, 
						help='standard error column')
	meta_parser.add_argument('--or-col', 
						action=AddString, 
						help='odds ratio column')
	meta_parser.add_argument('--z-col', 
						action=AddString, 
						help='z statistic column')
	meta_parser.add_argument('--p-col', 
						action=AddString, 
						help='p-value column')
	meta_parser.add_argument('--n-col', 
						action=AddString, 
						help='sample size column')
	meta_parser.add_argument('--gc', 
						action=AddString, 
						type=float, 
						help='set genomic inflation value instead of calculating it')
	meta_parser.add_argument('--n', 
						action=AddString, 
						help='set sample size')
	meta_parser.add_argument('--tag', 
						action=AddString, 
						help='tag for data file')
	meta_parser.add_argument('--file', 
						action=AddString, 
						help='results file name')
	meta_parser.add_argument('--meta', 
						action=AddString, 
						help='a meta analysis string')
	meta_parser.add_argument('--maf', 
						action=AddString, 
						type=float, 
						help='threshold value for minimum minor allele frequency (ie. 0.03 filters out markers with maf < 0.03)')
	meta_parser.add_argument('--rsq', 
						action=AddString, 
						type=float, 
						help='threshold value for imputation quality (ie. 0.8 filters out markers with rsq < 0.8)')
	meta_parser.add_argument('--hwe', 
						action=AddString, 
						type=float, 
						help='threshold value for Hardy Weinberg p-value (ie. 1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	meta_parser.add_argument('--method', 
						action=AddString, 
						choices=['sample_size', 'stderr'], 
						help='meta-analysis method (default: sample_size)')
	meta_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace existing output files')
	meta_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	meta_parser.add_argument('--buffer', 
						action=AddString, 
						type=int, 
						help='value for number of markers calculated at a time (WARNING: this argument will affect RAM memory usage; default: 100)')
	meta_parser.add_argument('--id', 
						action=AddString, 
						help='add region id to results (for use with --region option)')
	meta_parser_split_group1 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group1.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format (ie. 1:10583-1010582).')
	meta_parser_split_group1.add_argument('--region-file', 
						action=AddString, 
						help='filename for a list of tabix format regions')
	meta_parser_split_group2 = meta_parser.add_mutually_exclusive_group()
	meta_parser_split_group2.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='split region_file into an individual job for each line in file (requires --region-file)')
	meta_parser_split_group2.add_argument('--split-n', 
						action=AddString, 
						type=int, 
						help='split region_file into n individual jobs each with a subset of regions in the file (requires --region-file)')
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
	return meta_parser

def MapArgs(map_parser):
	map_required = map_parser.add_argument_group('required arguments')
	map_required.add_argument('--out', 
						action=AddString, 
						required=True, 
						help='output file name')
	map_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing out files')
	map_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	map_split_group1 = map_parser.add_mutually_exclusive_group()
	map_split_group1.add_argument('--mb', 
						action=AddString, 
						help='region size (megabase)')
	map_split_group1.add_argument('--kb', 
						action=AddString, 
						help='region size (kilobase)')
	map_split_group1.add_argument('--b', 
						action=AddString, 
						help='region size (base)')
	map_split_group1.add_argument('--n', 
						action=AddString, 
						help='number of markers to be included in each region')
	map_split_group2 = map_parser.add_mutually_exclusive_group()
	map_split_group2.add_argument('--oxford', 
						action=AddString, 
						help='oxford format genotype data file')
	map_split_group2.add_argument('--dos1', 
						action=AddString, 
						help='dos1 format genotype data file')
	map_split_group2.add_argument('--dos2', 
						action=AddString, 
						help='dos2 format genotype data file')
	map_split_group2.add_argument('--plink', 
						action=AddString, 
						help='Plink binary format genotype data file (without extension)')					
	map_split_group2.add_argument('--vcf', 
						action=AddString, 
						help='vcf 4.1/4.2 format genotype data file')
	map_split_group3 = map_parser.add_mutually_exclusive_group()
	map_split_group3.add_argument('--shift-mb', 
						action=AddString, 
						help='shift size (megabase)')
	map_split_group3.add_argument('--shift-kb', 
						action=AddString, 
						help='shift size (kilobase)')
	map_split_group3.add_argument('--shift-b', 
						action=AddString, 
						help='shift size (base)')
	map_split_group4 = map_parser.add_mutually_exclusive_group()
	map_split_group4.add_argument('--chr', 
						action=AddString, 
						type=int,  
						help='chromosome number from 1-26')
	map_split_group4.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format (ie. 1:10583-1010582).')
	return map_parser

def CompileArgs(compile_parser):
	compile_required = compile_parser.add_argument_group('required arguments')
	compile_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='base filename of existing results (basename only: do not include path or extension or list / region portion of filename; ex. set to X if out file names are of the form chr2/X.chr2bp1-2.gz or list0-99/X.list1.gz)')
	compile_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	compile_parser.add_argument('--split', 
						nargs=0, 
						action=AddTrue, 
						help='for use if all files to compile contain a header and an end of file marker')
	compile_parser.add_argument('--split-chr', 
						nargs=0, 
						action=AddTrue,  
						help='for use if each chromosome was run independently (chromosome directories contain a file with a header and a file with an end of file marker')
	return compile_parser

def EvalArgs(eval_parser):
	eval_required = eval_parser.add_argument_group('required arguments')
	eval_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='filename for results')
	eval_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	eval_parser.add_argument('--qq', 
						nargs=0, 
						action=AddTrue, 
						help='print qq plot')
	eval_parser.add_argument('--qq-strat', 
						nargs=0, 
						action=AddTrue, 
						help='print frequency stratified qq plot')
	eval_parser.add_argument('--qq-n', 
						nargs=0, 
						action=AddTrue, 
						help='print number of markers on qq plot')
	eval_parser.add_argument('--mht', 
						nargs=0, 
						action=AddTrue, 
						help='print manhattan plot')
	eval_parser.add_argument('--color', 
						nargs=0, 
						action=AddTrue, 
						help='plot in color')
	eval_parser.add_argument('--plot-gc', 
						nargs=0, 
						action=AddTrue, 
						help='print manhattan plots with genomic inflation corrected p-values')
	eval_parser.add_argument('--set-gc', 
						action=AddString, 
						type=float, 
						help='set genomic inflation value instead of calculating it')
	eval_parser.add_argument('--pmax', 
						action=AddString, 
						type=float, 
						help='set maximum p-value for top results file (default = 1e-4; will print at least top 100')
	eval_parser.add_argument('--stat', 
						action=AddString, 
						help='string indicating prefix for statistics to be summarized, not including tag (default: STAT=\'marker\')')
	eval_parser.add_argument('--top', 
						action=AddString, 
						type=int, 
						help='an integer; print regional plots for top n regions')
	eval_parser.add_argument('--tag', 
						action=AddString, 
						help='string indicating tag for stats to be summarized, if tag exists (example: TAG=aa and STAT=marker -> aa.marker.p)')
	eval_parser.add_argument('--unrel', 
						nargs=0, 
						action=AddTrue, 
						help='filter based on unrel columns')	
	eval_parser.add_argument('--rsq', 
						action=AddString, 
						type=float, 
						help='threshold for imputation quality (ie. 0.8 filters out markers with r-squared < 0.8)')
	eval_parser.add_argument('--maf', 
						action=AddString, 
						type=float, 
						help='threshold for allele frequency (ie. 0.03 filters out markers with MAF < 0.03)')
	eval_parser.add_argument('--hwe', 
						action=AddString, 
						type=float, 
						help='threshold for Hardy Weinberg p-value (ie. 1e-6 filters out markers with Hardy Weinberg p-value < 1e-6)')
	eval_parser.add_argument('--callrate', 
						action=AddString, 
						type=float, 
						help='threshold for callrate (ie. 0.95 filters out markers with callrate < 0.95)')
	eval_parser.add_argument('--effect', 
						action=AddString, 
						type=float, 
						help='threshold for effect estimate (ie. 1.7 filters out markers with effect estimate > 1.7 and < -1.7)')
	eval_parser.add_argument('--stderr', 
						action=AddString, 
						type=float, 
						help='threshold for standard error (ie. 5 filters out markers with standard error > 5)')
	eval_parser.add_argument('--oddsratio', 
						action=AddString, 
						type=float, 
						help='threshold for odds ratio (ie. 1.3 filters out markers with odds ratio > 1.25 and < 1/1.25 = 0.8)')
	eval_parser.add_argument('--df', 
						action=AddString, 
						type=int, 
						help='threshold for meta analysis degrees of freedom (ie. 4 filters out markers less than 5 datasets included in the meta analysis; requires --meta-dir)')
	eval_parser.add_argument('--sig', 
						action=AddString, 
						type=float, 
						help='line of significance p value (default: 5.4e-8)')
	eval_parser.add_argument('--ext', 
						action=AddString, 
						choices=['tiff','eps','pdf'], 
						help='file type extension for plot files (default: tiff)')
	eval_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	eval_parser_split_group = eval_parser.add_mutually_exclusive_group()
	eval_parser_split_group.add_argument('--region', 
						action=AddString, 
						help='genomic region specified in Tabix format (ie. 1:10583-1010582).')
	eval_parser_split_group.add_argument('--region-file', 
						action=AddString, 
						help='filename for a list of tabix format regions')
	eval_parser.add_argument('--lz-source', 
						action=AddString, 
						help='locuszoom source option')
	eval_parser.add_argument('--lz-build', 
						action=AddString, 
						help='locuszoom build option')
	eval_parser.add_argument('--lz-pop', 
						action=AddString, 
						help='locuszoom pop option')
	return eval_parser

def GcArgs(gc_parser):
	gc_required = gc_parser.add_argument_group('required arguments')
	gc_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='filename for results')
	gc_required.add_argument('--gc', 
						nargs=2, 
						required=True, 
						action=AddString, 
						help='apply genomic control to a 1 or more p-value columns (ex. --gc meta.p 1.0123 --gc meta.aa.p 1.002123)')
	gc_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	gc_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	return gc_parser

def AnnotArgs(annot_parser):
	annot_required = annot_parser.add_argument_group('required arguments')
	annot_required.add_argument('--file', 
						action=AddString, 
						required=True, 
						help='file name for single variant results')
	annot_parser.add_argument('--build', 
						action=AddString, 
						help='genomic build (default: GRCh37.75)')
	annot_parser.add_argument('--replace', 
						nargs=0, 
						action=AddTrue, 
						help='replace any existing output files')
	annot_parser.add_argument('--qsub', 
						action=AddString, 
						help='string indicating all qsub options to be added to the qsub command (triggers submission of all jobs to the cluster)')
	return annot_parser
