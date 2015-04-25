import pandas as pd
import numpy as np
import subprocess
import psutil
import resource
import math
from SystemFxns import Error
from FileFxns import Coordinates
import scipy.stats as scipy
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import re
import os
from Bio import bgzf

pd.options.mode.chained_assignment = None

#from memory_profiler import profile, memory_usage
#@profile
def GC(data, 
			out, 
			gc, 
			):

	print "gc options ..."
	for arg in locals().keys():
		if not locals()[arg] in [None, False]:
			print "   {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	##### read data from file #####
	print "loading results from file"
	reader = pd.read_table(data, sep='\t', chunksize=1000000,compression='gzip',dtype=object)
	bgzfile = bgzf.BgzfWriter(out + '.gc.gz', 'wb')
	i = 0
	for results in reader:
		i = i+1
		lines = len(results)
		if len(results) > 0:
			h = True if i == 1 else False
			for col in gc.keys():
				results[[x for x in results.columns if col.replace('.p','.effect') in x or col.replace('.p','.stderr') in x or col.replace('.p','.z') in x or col in x]] = results[[x for x in results.columns if col.replace('.p','.effect') in x or col.replace('.p','.stderr') in x or col.replace('.p','.z') in x or col in x]].astype(float)
				if col.replace('.p','.stderr') in results.columns:
					results[col.replace('.p','.stderr')] = results[col.replace('.p','.stderr')] * math.sqrt(float(gc[col]))
				if col.replace('.p','.z') in results.columns:
					results[col.replace('.p','.z')] = results[col.replace('.p','.z')] / math.sqrt(float(gc[col]))
				if col.replace('.p','.effect') in results.columns and col.replace('.p','.stderr') in results.columns:
					results[col] = 2 * scipy.norm.cdf(-1 * np.abs(results[col.replace('.p','.effect')]) / results[col.replace('.p','.stderr')])
				else:
					results[col] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*results[col]) / math.sqrt(float(gc[col]))))
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
	cmd = ['tabix','-b','2','-e','2',out + '.gc.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessError:
		print Error("file mapping failed")
	else:
		print 'process complete'
