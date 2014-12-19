from progressbar import ProgressBar, Counter, Timer
import sys
import os
import subprocess
from collections import OrderedDict
import numpy as np
from multiprocessing import Process, Manager, cpu_count
import math
from Bio import bgzf
import tabix
from Messages import Error

def RemoveExistingFiles(file, module):
	for f in [file, file + '.log', file + '.gz', file + '.gz.tbi']:
		try:
			os.remove(f)
		except OSError:
			continue
				
def CheckExistingFiles(file, module):
	for f in [file, file + '.log', file + '.gz', file + '.gz.tbi']:
		if os.path.exists(f):
			print Error("1 or more output files already exists (use --overwrite flag to replace)")
			sys.exit()
		else:
			pass

def PrepareChrDirs(regions, directory):
	try:
		os.mkdir(directory.replace('chr[CHR]/',''))
	except OSError:
		pass
	for chr in set([a.split(":")[0] for a in regions]):
		try:
			os.mkdir(directory.replace("[CHR]",chr))
		except OSError:
			continue

def PrepareListDirs(n, directory):
	try:
		os.mkdir(directory.replace('list[LIST]/',''))
	except OSError:
		pass
	for a in range(int(np.ceil(n/100.0))):
		try:
			os.mkdir(directory.replace("[LIST]",str(int(np.floor(a/100.0) + 100*a)) + '-' + str(int(np.floor(a/100.0) + 100*a+ + 99))))
		except OSError:
			continue
			
def GenerateSubFiles(region_df, f, dist_mode, n):
	out_files=OrderedDict()
	for i in range(n):
		if dist_mode in ['region','split-list']:
			of = f.replace('[CHR]',str(region_df['chr'][i])) + '.chr' + str(region_df['chr'][i]) + 'bp' + str(region_df['start'][i]) + '-' + str(region_df['end'][i])
			out_files[str(region_df['chr'][i]) + ':' + str(region_df['start'][i]) + '-' + str(region_df['end'][i])] = of
		elif dist_mode == 'split-list-n':
			of = f.replace('[LIST]',str(int(np.floor(i/100.0) + 99*np.floor(i/100.0))) + '-' + str(int(np.floor(i/100.0) + 99*np.floor(i/100.0) + 99))) + '.list' + str(i)
			out_files[i] = of
		elif dist_mode == 'chr':
			of = f + '.chr' + str(region_df['chr'][i])
			out_files[str(region_df['chr'][i])] = of
		else:
			print Error("   ... invalid dist_mode: " + dist_mode)
			return
	return out_files

def CheckResults(file_dict, out, cpus, complete_string, overwrite):
	if overwrite:
		print "   ... removing previous check results files"
	if os.path.exists(out + '.files.incomplete'):
		if overwrite:
			os.remove(out + '.files.incomplete')
		else :
			print "   ... file " + out + ".files.incomplete already exists (use --overwrite flag to replace the existing file)"
			return
	if os.path.exists(out + '.files.missing'):
		if overwrite:
			os.remove(out + '.files.missing')
		else :
			print "   ... file " + out + ".files.missing already exists (use --overwrite flag to replace the existing file)"
			return
	if os.path.exists(out + '.regions.rerun'):
		if overwrite:
			os.remove(out + '.regions.rerun')
		else :
			print "   ... file " + out + ".regions.rerun already exists (use --overwrite flag to replace the existing file)"
			return
	cpus = cpu_count() if cpu_count() < cpus else cpus	
	n = len(file_dict.keys())
	if n > 1:
		manager = Manager()
		jobs = []
		reg_complete_dict = manager.dict()
		reg_rerun_dict = manager.dict()
		reg_incompletefile_dict = manager.dict()
		reg_missingfile_dict = manager.dict()
		def CheckReg(i, joblist, reg_complete_dict,reg_rerun_dict,reg_incompletefile_dict,reg_missingfile_dict):
			reg_complete = []
			reg_rerun = []
			reg_incompletefile = []
			reg_missingfile = []
			for j in joblist:
				logfile = file_dict[j] + ".log"
				resfile = file_dict[j] + ".gz"
				if os.path.exists(resfile):
					if os.path.exists(logfile):
						p = subprocess.Popen(['grep','-cw',complete_string,logfile], stdout=subprocess.PIPE)	
						complete = p.communicate()[0]
						complete = int(complete.strip())
						if complete == 0:
							reg_rerun.append(str(j))
							reg_incompletefile.append(resfile)
						else:
							reg_complete.append(str(j))
					else:
						reg_missingfile.append(resfile)
						reg_rerun.append(str(j))
				else:
					reg_missingfile.append(resfile)
					reg_rerun.append(str(j))
			reg_complete_dict[i] = reg_complete
			reg_rerun_dict[i] = reg_rerun
			reg_incompletefile_dict[i] = reg_incompletefile
			reg_missingfile_dict[i] = reg_missingfile
			print "   ... finished processing on cpu " + str(i+1)
		joblist = [file_dict.keys()[i:i+int(math.ceil(float(len(file_dict.keys()))/cpus))] for i in range(0,len(file_dict.keys()),int(math.ceil(float(len(file_dict.keys()))/cpus)))]
		for i in range(len(joblist)):
			print "   ... submitting " + str(len(joblist[i])) + " result files to verify on cpu " + str(i+1)
			reg_complete_dict[i] = []
			reg_rerun_dict[i] = []
			reg_incompletefile_dict[i] = []
			reg_missingfile_dict[i] = []
			p = Process(target=CheckReg, args=(i, joblist[i],reg_complete_dict,reg_rerun_dict,reg_incompletefile_dict,reg_missingfile_dict))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()
		reg_complete = []
		reg_rerun = []
		reg_incompletefile = []
		reg_missingfile = []
		for i in range(len(joblist)):
			reg_complete.extend(reg_complete_dict[i])
			reg_rerun.extend(reg_rerun_dict[i])
			reg_incompletefile.extend(reg_incompletefile_dict[i])
			reg_missingfile.extend(reg_missingfile_dict[i])
		print "   ... verification results"
		print "          " + str(len(reg_complete)) + " of " + str(n) + " complete"
		if len(reg_incompletefile) > 0:
			print "          " + str(len(reg_incompletefile)) + " of " + str(n) + " incomplete"
			f = open(out + '.files.incomplete', 'w')
			f.write('\n'.join(reg_incompletefile))
			f.close()
			print "             written to file " + out + ".files.incomplete"
		if len(reg_missingfile) > 0:
			print "          " + str(len(reg_missingfile)) + " of " + str(n) + " files missing"
			f = open(out + '.files.missing', 'w')
			f.write('\n'.join(reg_missingfile))
			f.close()
			print "             written to file " + out + ".files.missing"
		if len(reg_rerun) > 0:
			print "          " + str(len(reg_rerun)) + " of " + str(n) + " total regions to rerun"
			f = open(out + '.regions.rerun', 'w')
			f.write('\n'.join(reg_rerun))
			f.close()
			print "             written to file " + out + ".regions.rerun"
	else:
		resfile = file_basename + '.chr' + file_dict.values()[0].split(":")[0] + 'bp' + file_dict.values()[0].split(":")[1].split("-")[0] + '-' + file_dict.values()[0].split(":")[1].split("-")[1]
		logfile = resfile + ".log"
		resfile = resfile + ".gz"
		if os.path.exists(resfile):
			if os.path.exists(logfile):
				p = subprocess.Popen(['grep','-cw',complete_string,logfile], stdout=subprocess.PIPE)	
				complete = p.communicate()[0]
				complete = int(complete.strip())
				if complete == 0:
					print "          analysis complete"
				else:
					print "          analysis incomplete"
			else:
				print "          analysis incomplete ... no log file found"
		else:
			print "          analysis incomplete ... no out file found"

def CompileResults(out_files, out, overwrite):
	print "   ... removing previous compiled results files"
	if os.path.exists(out + '.gz'):
		if overwrite:
			os.remove(out + '.gz')
		else :
			print "   ... file " + out + ".gz already exists (use --overwrite flag to replace the existing file)"
			return
	if os.path.exists(out + '.gz.tbi'):
		if overwrite:
			os.remove(out + '.gz.tbi')
		else:
			print "   ... file " + out + ".gz.tbi already exists (use --overwrite flag to replace the existing file)"
			return
	print "   ... compiling results"
	bgzfile = bgzf.BgzfWriter(out + '.gz', 'wb')
	logfile = file(out + '.log', 'w')
	pbar = ProgressBar(maxval=len(out_files.keys()), widgets = ['   ... processed ', Counter(), ' of ' + str(len(out_files.keys())) + ' regions (', Timer(), ')'])
	i=0
	pbar.start()
	for reg in out_files.keys():
		i = i + 1
		sed = ['awk','{print $0}'] if i == 1 else ['sed','1d']
		resfile = out_files[reg]
		resfile = resfile + ".gz"
		p1 = subprocess.Popen(['zcat',resfile], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(sed, stdin=p1.stdout, stdout=subprocess.PIPE)
		bgzfile.write(p2.communicate()[0])
		p1.wait()
		p2.wait()
		if i == 1:
			p3 = subprocess.Popen(['cat',out_files[reg] + '.log'], stdout=subprocess.PIPE)
			p4 = subprocess.Popen(['awk','{print \"      \"$0}'], stdin=p3.stdout, stdout=subprocess.PIPE)
			logfile.write('Sample log file from ' + out_files[reg] + '.log\n\n')
			logfile.write(p4.communicate()[0])
			p3.wait()
			p4.wait()
			logfile.write('\nElapsed time and max memory used for each job in list\n\n')
		p5 = subprocess.Popen(['grep','time elapsed:',out_files[reg] + '.log'], stdout=subprocess.PIPE)
		elap = p5.communicate()[0].strip()
		p5.wait()
		p6 = subprocess.Popen(['grep','memory used:',out_files[reg] + '.log'], stdout=subprocess.PIPE)
		maxmem = p6.communicate()[0].strip()
		p6.wait()
		logfile.write('      job ' + str(reg) + '   -   ' + elap + '   ' + maxmem + '\n')
		pbar.update(i)
	pbar.finish()
	#p1.wait()
	#p2.wait()
	bgzfile.close()
	logfile.close()
	print "   ... mapping compiled file"
	cmd = 'tabix -b 2 -e 2 ' + out + '.gz'
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print "   ... file compilation complete\n"
