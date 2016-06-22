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
from Bio import bgzf
import Process
import subprocess
import multiprocessing as mp
import sys
import os
import resource
import logging
import pickle
import glob
import pysam

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunSnvgroup")

def process_regions(regions_df, cfg, cpu, log):
	regions_df = regions_df[regions_df['cpu'] == cpu].reset_index(drop=True)

	if log:
		try:
			log_file = open(cfg['out'] + '.cpu' + str(cpu) + '.log','w')
		except:
			print Process.Error("unable to initialize log file " + cfg['out'] + '.cpu' + str(cpu) + '.log').out
			return 1
		stdout_orig = sys.stdout
		sys.stdout = log_file

	tool_error = False
	print ''
	print 'sourcing bash script ' + cfg['source']
	with open(cfg['source'],'r') as s:
		base_cmd = s.read()
	for k in xrange(len(regions_df.index)):
		print ''
		print 'loading region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['region'][k] + ') ...'
		cmd = base_cmd.replace('UGA_FILE',cfg['file']).replace('UGA_OUT',cfg['out'] + '.cpu' + str(cpu) + '.chr' + regions_df['region'][k].replace(':','bp')).replace('UGA_REGION_BP',regions_df['region'][k].replace(':','bp')).replace('UGA_REGION',regions_df['region'][k]).replace('UGA_OUT',cfg['out'])
		with open(cfg['out'] + '.source','w') as s:
			s.write(cmd)
		try:
			p = subprocess.Popen(['sh','-c', cmd],stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
			for line in iter(p.stdout.readline, ''):
				sys.stdout.write(line)
			p.wait()
		except KeyboardInterrupt:
			kill_all(p.pid)
			tool_error = True
			pass

		status = 'processed region ' + str(k+1) + '/' + str(len(regions_df.index)) + ' (' + regions_df['region'][k] + ') ...'
		print status
		sys.stdout.flush()

	if log:
		sys.stdout = stdout_orig
		log_file.close()

	if tool_error:
		return -1
	else:
		return 0

def RunTools(args):
	cfg = Parse.generate_tools_cfg(args)
	Parse.print_tools_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	regions_df = pd.read_table(cfg['region_file'], compression='gzip' if cfg['region_file'].split('.')[-1] == 'gz' else None)
	regions_df = regions_df[regions_df['job'] == int(cfg['job'])].reset_index(drop=True)
	return_values = {}
	print ''
	print "initializing out file"
	try:
		bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')
	except:
		print Process.Error("failed to initialize bgzip format out file " + cfg['out'] + '.gz').out
		return 1

	if cfg['cpus'] > 1:
		pool = mp.Pool(cfg['cpus']-1)
		for i in xrange(1,cfg['cpus']):
			return_values[i] = pool.apply_async(process_regions, args=(regions_df,cfg,i,True,))
			print "submitting job on cpu " + str(i) + " of " + str(cfg['cpus'])
		pool.close()
		print "executing job for cpu " + str(cfg['cpus']) + " of " + str(cfg['cpus']) + " via main process"
		main_return = process_regions(regions_df,cfg,cfg['cpus'],True)
		pool.join()

		if 1 in [return_values[i].get() for i in return_values] or main_return == 1:
			print Process.Error("error detected, see log files").out
			return 1

	else:
		main_return = process_regions(regions_df,cfg,1,True)
		if main_return == 1:
			print Process.Error("error detected, see log files").out
			return 1

	for i in xrange(1,cfg['cpus']+1):
		try:
			logfile = open(cfg['out'] + '.cpu' + str(i) + '.log', 'r')
		except:
			print Process.Error("failed to initialize log file " + cfg['out'] + '.cpu' + str(i) + '.log').out
			return 1
		print logfile.read()
		logfile.close()
		os.remove(cfg['out'] + '.cpu' + str(i) + '.log')

	written = False
	for i in xrange(1,cfg['cpus']+1):
		cpu_regions_df = regions_df[regions_df['cpu'] == i].reset_index()
		for j in xrange(0,len(cpu_regions_df.index)):
			f_temp=glob.glob(cfg['out'] + '.cpu' + str(i) + '.chr' + cpu_regions_df['region'][j].replace(':','bp') + '*.gz')[0]
			try:
				h=pysam.TabixFile(filename=f_temp,parser=pysam.asVCF())
			except:
				print Process.Error("failed to load vcf file " + f_temp)
				return 1
			if not written:
				for row in h.header:
					bgzfile.write(str(row) + '\n')
				written = True
			h_iter = h.fetch(region=str(cpu_regions_df['chr'][j]))
			for row in h_iter:
				bgzfile.write(str(row) + '\n')
			for f in glob.glob(cfg['out'] + '.cpu' + str(i) + '.chr' + cpu_regions_df['region'][j].replace(':','bp') + '.*'):
				os.remove(f)

	bgzfile.close()

	print "indexing out file"
	try:
		pysam.tabix_index(cfg['out'] + '.gz',preset="vcf",force=True)
	except:
		print Process.Error('failed to generate index').out
		return 1

	print "process complete"
	return 0
