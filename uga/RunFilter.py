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
	if cfg['miss'] is not None and not 'callrate' in cols:
		print Process.Error("callrate column not found, required for --miss option").out
		found = False
	if cfg['maf'] is not None and not 'freq' in cols:
		print Process.Error("freq column not found, required for --maf option").out
		found = False
	if cfg['mac'] is not None and not 'mac' in cols:
		print Process.Error("mac column not found, required for --mac option").out
		found = False
	if cfg['cmac'] is not None and not 'cmac' in cols:
		print Process.Error("cmac column not found, required for --cmac option").out
		found = False
	if cfg['rsq'] is not None and not 'rsq' in cols:
		print Process.Error("rsq column not found, required for --rsq option").out
		found = False
	if cfg['hwe'] is not None and not 'hwe' in cols:
		print Process.Error("hwe column not found, required for --hwe option").out
		found = False
	if cfg['hwe_maf'] is not None and (not 'hwe' in cols or not 'freq' in cols):
		print Process.Error("either hwe or freq column not found, both required for --hwe option").out
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
	r = r.loc[~ np.isnan(r['p'])]
	print str(r.shape[0]) + " results found with valid p-values"

	nsnps = r.shape[0]
	if cfg['miss'] is not None:
		r = r.loc[r['callrate'] >= cfg['miss']]
	if cfg['maf'] is not None:
		r = r.loc[(r['freq'] >= cfg['maf']) & (r['freq'] <= 1-cfg['maf'])]
	if cfg['mac'] is not None:
		r = r.loc[r['mac'] >= cfg['mac']]
	if cfg['cmac'] is not None:
		r = r.loc[r['cmac'] >= cfg['cmac']]
	if cfg['rsq'] is not None:
		r = r.loc[(~ np.isnan(r['rsq'])) & (r['rsq'] >= cfg['rsq'])]
	if cfg['hwe'] is not None:
		if cfg['hwe_maf'] is not None:
			r = r.loc[(~ np.isnan(r['hwe'])) & (~ np.isnan(r['freq'])) & (~ (r['freq'] >= cfg['hwe_maf']) & (r['hwe'] < cfg['hwe']))]
		else:
			r = r.loc[(~ np.isnan(r['hwe'])) & (r['hwe'] >= cfg['hwe'])]
	print str(r.shape[0]) + " results remain after filtering, " + str(nsnps - r.shape[0]) + " removed"

	if cfg['gc']:
		l = np.median(scipy.chi2.ppf([1-x for x in r.loc[~ np.isnan(r['p']),'p'].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
		print "genomic inflation = " + str(l)

		print "adjusting stderr"
		if 'stderr' in r.columns:
			r['stderr'] = r['stderr'] * math.sqrt(l)

		if 'wald' in r.columns:
			print "adjusting wald statistic"
			r['wald'] = r['wald'] / math.sqrt(l)
			print "calculating corrected p-value from wald statistic"
			r['p'] = scipy.chisqprob(r['wald'],1)
		elif 'z' in r.columns:
			print "adjusting z statistic"
			r['z'] = r['z'] / math.sqrt(l)
			print "calculating corrected p-value from z statistic"
			r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(r['z']))
		elif 'effect' in r.columns and 'stderr' in r.columns:
			print "calculating corrected p-value from effect and stderr using a calculated z statistic"
			r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(r['effect']) / r['stderr'])
		else:
			print "calculating corrected p-value from existing p-value using an estimated z statistic"
			r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*r['p']) / math.sqrt(l)))

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
		pysam.tabix_index(cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz'),seq_col=0,start_col=r.columns.get_loc('pos'),end_col=r.columns.get_loc('pos'),force=True)
	except:
		print Process.Error('failed to generate index for file ' + cfg['file'].replace('.gz','.' + cfg['tag'] + '.gz')).out
		return 1

	print "process complete"
	return 0
