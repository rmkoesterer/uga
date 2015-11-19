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
		f = directory + '/' + row['file']
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
			f = directory + '/' + row['file']
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
	logfile = file(directory + '/' + directory + '.log', 'w')
	files_o = files[files['out'] == out[0]].reset_index(drop=True)
	pbar = ProgressBar(maxval=files_o.shape[0], widgets = ['   processed ', Counter(), ' of ' + str(files_o.shape[0]) + ' results (', Timer(), ')'])
	pbar.start()
	for j, row in files_o.iterrows():
		f = directory + '/' + row['file']
		lf = glob.glob('/'.join(f.split('/')[0:len(f.split('/'))-1]) + "/*.log")
		if j+1 == 1:
			p3 = subprocess.Popen(['cat',lf[0]], stdout=subprocess.PIPE)
			p4 = subprocess.Popen(['awk','{print \"      \"$0}'], stdin=p3.stdout, stdout=subprocess.PIPE)
			logfile.write('Sample log file from ' + lf[0] + '\n\n')
			logfile.write(p4.communicate()[0])
			p3.wait()
			p4.wait()
			logfile.write('\nElapsed time and max memory used for each job in list\n\n')
		p5 = subprocess.Popen(['grep','time elapsed:',lf[0]], stdout=subprocess.PIPE)
		elap = p5.communicate()[0].strip()
		p5.wait()
		p6 = subprocess.Popen(['grep','max memory used by main process:',lf[0]], stdout=subprocess.PIPE)
		maxmem = p6.communicate()[0].strip()
		p6.wait()
		p7 = subprocess.Popen(['grep','max memory used by any subprocess:',lf[0]], stdout=subprocess.PIPE)
		maxmemsub = p7.communicate()[0].strip()
		p7.wait()
		logfile.write('job ' + str(j+1) + ' - time elapsed ' + elap + ' max mem main ' + maxmem + ' max mem sub ' + maxmemsub + '\n')
		pbar.update(j)
	pbar.finish()
	logfile.close()
	for o in out:
		print "mapping compiled file " + o
		files_o = files[files['out'] == o].reset_index(drop=True)
		p7 = subprocess.Popen(['tabix','-H',directory + '/' + files.iloc[0]['file']], stdout=subprocess.PIPE)
		header_full = p7.communicate()[0].splitlines()
		p7.wait()
		header = header_full[-1].split()
		if header[1] == 'pos':
			b = 1
			e = 1
		else:
			b = 1
			e = 2
		pysam.tabix_index(directory + '/' + o,seq_col=0,start_col=b,end_col=e,force=True)
	print "file compilation complete"
	return True
