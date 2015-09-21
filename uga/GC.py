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

from __main__ import *
import scipy.stats as scipy
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

def GC(cfg):
	Parse.PrintGcOptions(cfg)

	##### read data from file #####
	print "loading results from file"
	reader = pd.read_table(cfg['file'], sep='\t', chunksize=1000000,compression='gzip',dtype=object)
	bgzfile = bgzf.BgzfWriter(cfg['file'].replace('.gz','') + '.gc.gz', 'wb')
	i = 0
	for results in reader:
		i = i+1
		lines = len(results)
		if len(results) > 0:
			h = True if i == 1 else False
			for col in cfg['gc'].keys():
				results[[x for x in results.columns if col.replace('.p','.effect') in x or col.replace('.p','.stderr') in x or col.replace('.p','.z') in x or col in x]] = results[[x for x in results.columns if col.replace('.p','.effect') in x or col.replace('.p','.stderr') in x or col.replace('.p','.z') in x or col in x]].astype(float)
				if col.replace('.p','.stderr') in results.columns:
					results[col.replace('.p','.stderr')] = results[col.replace('.p','.stderr')] * math.sqrt(float(cfg['gc'][col]))
				if col.replace('.p','.z') in results.columns:
					results[col.replace('.p','.z')] = results[col.replace('.p','.z')] / math.sqrt(float(cfg['gc'][col]))
				if col.replace('.p','.effect') in results.columns and col.replace('.p','.stderr') in results.columns:
					results[col] = 2 * scipy.norm.cdf(-1 * np.abs(results[col.replace('.p','.effect')]) / results[col.replace('.p','.stderr')])
				else:
					results[col] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*results[col]) / math.sqrt(float(cfg['gc'][col]))))
				if col.replace('.p','.stderr') in results.columns:
					results[col.replace('.p','.stderr')] = results[col.replace('.p','.stderr')].map(lambda x: '%.5g' % x if not math.isnan(x) else 'NA')
				if col.replace('.p','.z') in results.columns:
					results[col.replace('.p','.z')] = results[col.replace('.p','.z')].map(lambda x: '%.5g' % x if not math.isnan(x) else 'NA')
				results[col] = results[col].map(lambda x: '%.2e' % x if not math.isnan(x) else 'NA')
			results.fillna('NA',inplace=True)
			results.to_csv(bgzfile, header=h, index=False, sep="\t")
			del results
		print "   processed " + str((i-1)*1000000 + lines) + " lines"

	bgzfile.close()

	print "mapping results file"
	cmd = ['tabix','-b','2','-e','2',cfg['file'].replace('.gz','') + '.gc.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessError:
		print SystemFxns.Error("file mapping failed")
	else:
		print 'process complete'
