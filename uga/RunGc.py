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
logger = logging.getLogger("RunGc")

def RunGc(args):
	cfg = Parse.generate_gc_cfg(args)
	Parse.print_gc_options(cfg)

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

	print "reading data from file"
	skip_rows = len(header)-1
	cols = header[-1].split()
	try:
		r = pd.read_table(cfg['file'],sep='\t',skiprows=skip_rows,compression='gzip')
	except:
		print Process.Error("unable to read data from file " + cfg['file']).out
		return 1
	print str(r.shape[0]) + " variants found"

	l = np.median(scipy.chi2.ppf([1-x for x in r.loc[~ np.isnan(r['p']),'p'].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
	print "   genomic inflation = " + str(l)

	print "   adjusting stderr"
	if 'stderr' in r.columns:
		r['stderr'] = r['stderr'] * math.sqrt(l)

	if 'wald' in r.columns:
		print "   adjusting wald statistic"
		r['wald'] = r['wald'] / math.sqrt(l)
		print "   calculating corrected p-value from wald statistic"
		r['p'] = scipy.chisqprob(r['wald'],1)
	elif 'z' in r.columns:
		print "   adjusting z statistic"
		r['z'] = r['z'] / math.sqrt(l)
		print "   calculating corrected p-value from z statistic"
		r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(r['z']))
	elif 'effect' in r.columns and 'stderr' in r.columns:
		print "   calculating corrected p-value from effect and stderr using a calculated z statistic"
		r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(r['effect']) / r['stderr'])
	else:
		print "   calculating corrected p-value from existing p-value using an estimated z statistic"
		r['p'] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*r['p']) / math.sqrt(l)))

	print "writing inflation corrected results to file"
	try:
		bgzfile = bgzf.BgzfWriter(cfg['file'].replace('.gz','.gc.gz'), 'wb')
	except:
		print Process.Error("unable to initialize out file " + cfg['file'].replace('.gz','.gc.gz')).out
		return 1
	bgzfile.write('\n'.join([x for x in handle.header]) + '\n')
	r[cols].to_csv(bgzfile,header=False,index=False,sep="\t",na_rep='NA', float_format='%.5g')
	bgzfile.close()
	handle.close()

	print "indexing out file"
	try:
		pysam.tabix_index(cfg['file'].replace('.gz','.gc.gz'),seq_col=0,start_col=r.columns.get_loc('pos'),end_col=r.columns.get_loc('pos'),force=True)
	except:
		print Process.Error('failed to generate index for file ' + cfg['file'].replace('.gz','.gc.gz')).out
		return 1

	print "process complete"
	return 0
