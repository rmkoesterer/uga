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

from progressbar import ProgressBar, Counter, Timer
import subprocess
import signal
import pandas as pd
import numpy as np
import glob
import os
from Bio import bgzf
import pysam
import gzip

def get_delimiter(d):
	if d == 'tab':
		d = '\t'
	elif d == 'space':
		d = ' '
	else:
		d = ','
	return d

def verify_results(directory, files):
	print "verifying results"
	reg_complete = []
	reg_rerun = []
	pbar = ProgressBar(maxval=files.shape[0], widgets = ['   processed ', Counter(), ' of ' + str(files.shape[0]) + ' files (', Timer(), ')'])
	pbar.start()
	for i, row in files.iterrows():
		f = directory + '/' + '/'.join(row['file'].split('/')[1:])
		j = row['job']
		logfile = glob.glob('/'.join(f.split('/')[0:len(f.split('/'))-1]) + "/*.log")
		if os.path.exists(f):
			if len(logfile) > 0:
				p = subprocess.Popen(['grep','-cw','process complete',logfile[0]], stdout=subprocess.PIPE)	
				complete = p.communicate()[0]
				complete = int(complete.strip())
				if complete == 0:
					reg_rerun.append(j)
				else:
					reg_complete.append(j)
			else:
				reg_rerun.append(j)
		else:
			reg_rerun.append(j)
		pbar.update(i)
	pbar.finish()
	return list(set(reg_complete)), list(set(reg_rerun))

def compile_results(directory, files):
	out = np.unique(files['out'])
	bgzfile = {}
	for o in out:
		print "compiling results to file " + o
		files_o = files[files['out'] == o].reset_index(drop=True)
		pbar = ProgressBar(maxval=files_o.shape[0], widgets = ['   processed ', Counter(), ' of ' + str(files_o.shape[0]) + ' files (', Timer(), ')'])
		pbar.start()
		bgzfile[o] = bgzf.BgzfWriter(directory + '/' + o, 'wb')
		for j, row in files_o.iterrows():
			f = directory + '/' + '/'.join(row['file'].split('/')[1:])
			sed = ['awk','{print $0}'] if j+1 == 1 else ['grep','-v','^#']
			p1 = subprocess.Popen(['zcat',f], stdout=subprocess.PIPE, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
			p2 = subprocess.Popen(sed, stdin=p1.stdout, stdout=subprocess.PIPE)
			bgzfile[o].write(p2.communicate()[0])
			p1.wait()
			p2.wait()
			pbar.update(j)
		pbar.finish()
		bgzfile[o].close()

	print "compiling log files"
	summary = file(directory + '/' + os.path.basename(directory) + '.summary', 'w')
	logs = file(directory + '/' + os.path.basename(directory) + '.logs', 'w')
	files_o = files[files['out'] == out[0]].reset_index(drop=True)
	pbar = ProgressBar(maxval=files_o.shape[0], widgets = ['   processed ', Counter(), ' of ' + str(files_o.shape[0]) + ' results (', Timer(), ')'])
	pbar.start()
	for j, row in files_o.iterrows():
		f = directory + '/' + '/'.join(row['file'].split('/')[1:])
		lf = glob.glob('/'.join(f.split('/')[0:len(f.split('/'))-1]) + "/*.log")
		if j+1 == 1:
			p1 = subprocess.Popen(['cat',lf[0]], stdout=subprocess.PIPE)
			p2 = subprocess.Popen(['awk','{print \"      \"$0}'], stdin=p1.stdout, stdout=subprocess.PIPE)
			summary.write('Sample log file from ' + lf[0] + '\n\n')
			summary.write(p2.communicate()[0])
			p1.wait()
			p2.wait()
			summary.write('\nElapsed time and max memory used for each job in list\n\n')
		p = subprocess.Popen(['grep','time elapsed:',lf[0]], stdout=subprocess.PIPE)
		elap = p.communicate()[0].strip()
		p.wait()
		p = subprocess.Popen(['grep','max memory used by main process:',lf[0]], stdout=subprocess.PIPE)
		maxmem = p.communicate()[0].strip()
		p.wait()
		p = subprocess.Popen(['grep','max memory used by any subprocess:',lf[0]], stdout=subprocess.PIPE)
		maxmemsub = p.communicate()[0].strip()
		p.wait()
		summary.write('job ' + str(j+1) + ' - ' + elap + ' - ' + maxmem + ' - ' + maxmemsub + '\n')
		p = subprocess.Popen(['cat',lf[0]], stdout=subprocess.PIPE)
		logs.write(p.communicate()[0] + '\n\n')
		p.wait()
		pbar.update(j)
	pbar.finish()
	summary.close()
	logs.close()
	for o in out:
		print "mapping compiled file " + o
		files_o = files[files['out'] == o].reset_index(drop=True)
		h=pysam.TabixFile(filename=directory + '/' + '/'.join(files.iloc[0]['file'].split('/')[1:]),parser=pysam.asTuple())
		header = [x for x in h.header]
		cols = header[-1].split()
		source = header[0]
		if '## source: uga' in source or "#chr" in source:
			if cols[1] == 'pos':
				b = 1
				e = 1
			else:
				b = 1
				e = 2
			pysam.tabix_index(directory + '/' + o,seq_col=0,start_col=b,end_col=e,force=True)
		elif '##fileformat=VCF' in source or "#CHROM\tBEGIN\tEND\tMARKER_ID" in source or "#CHROM\tBEG\tEND\tMARKER_ID" in source: # if labeled as VCF or EPACTS results
			pysam.tabix_index(directory + '/' + o,preset='vcf',force=True)
		else:
			print "compiled file source not recognized"
			return False
	print "file compilation complete"
	return True
