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

import pandas as pd
import numpy as np
import Parse
import pysam
from Bio import bgzf
import scipy.stats as scipy
import math
import Process
import logging

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunFilter")

def RunFilter(args):
	cfg = Parse.generate_filter_cfg(args)
	Parse.print_filter_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	print ''
	print "loading file header"
	try:
		handle=pysam.TabixFile(filename=cfg['file'],parser=pysam.asTuple())
	except:
		print Process.Error("unable to load file header").out
		return 1
	header = [x for x in handle.header]
	cols = header[-1].split()

	found = True
	if not cfg['pcol'] in cols:
		print Process.Error("p-value column, --pcol, not found").out
		found = False
	if not cfg['bpcol'] in cols:
		print Process.Error("genomic position column, --bpcol, not found").out
		found = False
	if cfg['miss'] is not None and not cfg['misscol'] in cols:
		print Process.Error("callrate column, --misscol, not found; required for --miss option").out
		found = False
	if cfg['maf'] is not None and not cfg['freqcol'] in cols:
		print Process.Error("allele frequency column, --freqcol, not found; required for --maf option").out
		found = False
	if cfg['mac'] is not None and not cfg['maccol'] in cols:
		print Process.Error("minor allele count column, --maccol, not found; required for --mac option").out
		found = False
	if cfg['cmac'] is not None and not cfg['cmaccol'] in cols:
		print Process.Error("cumulative minor allele count column, --cmaccol, not found; required for --cmac option").out
		found = False
	if cfg['rsq'] is not None and not cfg['rsqcol'] in cols:
		print Process.Error("imputation quality (rsq) column, --rsqcol, not found; required for --rsq option").out
		found = False
	if cfg['hwe'] is not None and not cfg['hwecol'] in cols:
		print Process.Error("Hardy Weinberg p-value column, --hwecol, not found; required for --hwe option").out
		found = False
	if cfg['hwe_maf'] is not None and (not cfg['hwecol'] in cols or not cfg['freqcol'] in cols):
		print Process.Error("either Hardy Weinberg p-value or allele frequency column, --hwecol or --freqcol, not found; both required for --hwe-maf option").out
		found = False
	if not found:
		return 1

	print "reading data from file"
	skip_rows = len(header)-1
	cols = header[-1].split()
	try:
		r = pd.read_table(cfg['file'],sep='\t',skiprows=skip_rows,compression='gzip')
	except:
		print Process.Error("unable to read data from file " + cfg['file']).out
		return 1
	r = r.loc[~ np.isnan(r[cfg['pcol']])]
	print str(r.shape[0]) + " results found with valid p-values"

	nsnps = r.shape[0]
	if cfg['miss'] is not None:
		r = r.loc[r[cfg['misscol']] >= cfg['miss']]
	if cfg['maf'] is not None:
		r = r.loc[(r[cfg['freqcol']] >= cfg['maf']) & (r[cfg['freqcol']] <= 1-cfg['maf'])]
	if cfg['mac'] is not None:
		r = r.loc[r[cfg['maccol']] >= cfg['mac']]
	if cfg['cmac'] is not None:
		r = r.loc[r[cfg['cmaccol']] >= cfg['cmac']]
	if cfg['rsq'] is not None:
		r = r.loc[(~ np.isnan(r[cfg['rsqcol']])) & (r[cfg['rsqcol']] >= cfg['rsq'])]
	if cfg['hwe'] is not None:
		if cfg['hwe_maf'] is not None:
			r = r.loc[(~ np.isnan(r[cfg['hwecol']])) & (~ np.isnan(r[cfg['freqcol']])) & (~ (r[cfg['freqcol']] >= cfg['hwe_maf']) & (r[cfg['hwecol']] < cfg['hwe']))]
		else:
			r = r.loc[(~ np.isnan(r[cfg['hwecol']])) & (r[cfg['hwecol']] >= cfg['hwe'])]
	print str(r.shape[0]) + " results remain after filtering, " + str(nsnps - r.shape[0]) + " removed"

	if cfg['gc']:
		l = np.median(scipy.chi2.ppf([1-x for x in r.loc[~ np.isnan(r[cfg['pcol']]),cfg['pcol']].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
		print "genomic inflation = " + str(l)

		if cfg['stderrcol'] in r.columns:
			print "adjusting stderr"
			r[cfg['stderrcol']] = r[cfg['stderrcol']] * math.sqrt(l)
		if cfg['waldcol'] in r.columns:
			print "adjusting wald statistic"
			r[cfg['waldcol']] = r[cfg['waldcol']] / math.sqrt(l)
			print "calculating corrected p-value from wald statistic"
			r[cfg['pcol']] = scipy.chisqprob(r[cfg['waldcol']],1)
		elif cfg['zcol'] in r.columns:
			print "adjusting z statistic"
			r[cfg['zcol']] = r[cfg['zcol']] / math.sqrt(l)
			print "calculating corrected p-value from z statistic"
			r[cfg['pcol']] = 2 * scipy.norm.cdf(-1 * np.abs(r[cfg['zcol']]))
		elif cfg['effectcol'] in r.columns and cfg['stderrcol'] in r.columns:
			print "calculating corrected p-value from effect and stderr using a calculated z statistic"
			r[cfg['pcol']] = 2 * scipy.norm.cdf(-1 * np.abs(r[cfg['effectcol']]) / r[cfg['stderrcol']])
		else:
			print "calculating corrected p-value from existing p-value using an estimated z statistic"
			r[cfg['pcol']] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*r[cfg['pcol']]) / math.sqrt(l)))

	print "writing filtered results to file"
	try:
		bgzfile = bgzf.BgzfWriter(cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz'), 'wb')
	except:
		print Process.Error("unable to initialize out file " + cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz')).out
		return 1
	bgzfile.write('\n'.join([x for x in handle.header]) + '\n')
	r[cols].to_csv(bgzfile,header=False,index=False,sep="\t",na_rep='NA', float_format='%.5g')
	bgzfile.close()
	handle.close()

	print "indexing out file"
	try:
		pysam.tabix_index(cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz'),seq_col=0,start_col=r.columns.get_loc(cfg['bpcol']),end_col=r.columns.get_loc(cfg['bpcol']),force=True)
	except:
		print Process.Error('failed to generate index for file ' + cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz')).out
		return 1

	print "process complete"
	return 0
