from __main__ import *
from progressbar import ProgressBar, Counter, Timer
from multiprocessing import Process, Manager, cpu_count
import signal

def LoadCoordinates(filename):
	with open(filename) as f:
		regions = pd.read_table(f, header=None)
	if regions.shape[1] == 2:
		regions.columns=['region','id']
	elif regions.shape[1] == 1:
		regions.columns=['region']
	else:
		print SystemFxns.Error("too many columns in region list")
		return
	if not 'id' in regions.columns:
		regions['id'] = 'NA'
	regions['chr'] = regions['region'].apply(lambda row: int(row.split(':')[0]),1)
	regions['start'] = regions['region'].apply(lambda row: int(row.split(':')[1].split('-')[0]),1)
	regions['end'] = regions['region'].apply(lambda row: int(row.split(':')[1].split('-')[1]),1)
	regions['chr'] = regions['chr'].astype(int)
	regions['start'] = regions['start'].astype(int)
	regions['end'] = regions['end'].astype(int)
	regions = regions.sort_index(by=['chr','start'],ascending=[True,True])
	return regions

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
	out_files=collections.OrderedDict()
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
			print SystemFxns.Error("invalid dist_mode: " + dist_mode)
			return
	return out_files

def CheckResults(file_dict, out, cpus=1):
	print "verifying results"
	cpus = cpu_count() if cpu_count() < cpus else cpus	
	n = len(file_dict.keys())
	if n > 1:
		manager = Manager()
		jobs = []
		reg_complete_dict = manager.dict()
		reg_rerun_dict = manager.dict()
		reg_incompletefile_dict = manager.dict()
		reg_missingfile_dict = manager.dict()
		reg_nomarkers_dict = manager.dict()
		joblist = [file_dict.keys()[i:i+int(math.ceil(float(len(file_dict.keys()))/cpus))] for i in range(0,len(file_dict.keys()),int(math.ceil(float(len(file_dict.keys()))/cpus)))]
		if cpus == 1:
			pbar = ProgressBar(maxval=len(joblist[i]), widgets = ['   processed ', Counter(), ' of ' + str(len(joblist[i])) + ' regions (', Timer(), ')'])
		def CheckReg(i, joblist, reg_complete_dict,reg_rerun_dict,reg_incompletefile_dict,reg_missingfile_dict,reg_nomarkers_dict):
			reg_complete = []
			reg_rerun = []
			reg_incompletefile = []
			reg_missingfile = []
			reg_nomarkers = []
			k = 0
			for j in joblist:
				k = k + 1
				logfile = glob.glob(file_dict[j] + ".*.log")
				resfile = file_dict[j] + ".gz"
				if os.path.exists(resfile):
					if os.path.exists(logfile[0]):
						p = subprocess.Popen(['grep','-cw','process complete',logfile[0]], stdout=subprocess.PIPE)	
						complete = p.communicate()[0]
						complete = int(complete.strip())
						if complete == 0:
							p = subprocess.Popen(['grep','-cw','no markers found',logfile[0]], stdout=subprocess.PIPE)	
							no_makers = p.communicate()[0]
							no_makers = int(no_makers.strip())
							if no_makers == 0:
								reg_rerun.append(str(j))
								reg_incompletefile.append(resfile)
							else:
								reg_nomarkers.append(str(j))
						else:
							reg_complete.append(str(j))
					else:
						reg_missingfile.append(resfile)
						reg_rerun.append(str(j))
				else:
					p = subprocess.Popen(['grep','-cw','no markers found',logfile[0]], stdout=subprocess.PIPE)	
					no_makers = p.communicate()[0]
					no_makers = int(no_makers.strip())
					if no_makers == 0:
						reg_rerun.append(str(j))
						reg_missingfile.append(resfile)
					else:
						reg_nomarkers.append(str(j))
				if cpus == 1:
					pbar.update(k)
			reg_complete_dict[i] = reg_complete
			reg_rerun_dict[i] = reg_rerun
			reg_incompletefile_dict[i] = reg_incompletefile
			reg_missingfile_dict[i] = reg_missingfile
			reg_nomarkers_dict[i] = reg_nomarkers
			if cpus > 1:
				print "   finished processing on cpu " + str(i+1)
		for i in range(len(joblist)):
			print "submitting " + str(len(joblist[i])) + " result files for verification on cpu " + str(i+1)
			reg_complete_dict[i] = []
			reg_rerun_dict[i] = []
			reg_incompletefile_dict[i] = []
			reg_missingfile_dict[i] = []
			if cpus == 1:
				pbar.start()
			p = Process(target=CheckReg, args=(i, joblist[i],reg_complete_dict,reg_rerun_dict,reg_incompletefile_dict,reg_missingfile_dict,reg_nomarkers_dict))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()
		reg_complete = []
		reg_rerun = []
		reg_incompletefile = []
		reg_missingfile = []
		reg_nomarkers = []
		for i in range(len(joblist)):
			reg_complete.extend(reg_complete_dict[i])
			reg_rerun.extend(reg_rerun_dict[i])
			reg_incompletefile.extend(reg_incompletefile_dict[i])
			reg_missingfile.extend(reg_missingfile_dict[i])
			reg_nomarkers.extend(reg_nomarkers_dict[i])
		if cpus == 1:
			pbar.finish()
		print "verification status ... "
		print "   " + str(len(reg_complete) + len(reg_nomarkers)) + " of " + str(n) + " complete"
		if len(reg_incompletefile) > 0:
			print "   " + str(len(reg_incompletefile)) + " of " + str(n) + " incomplete"
			f = open(out + '.files.incomplete', 'w')
			f.write('\n'.join(reg_incompletefile))
			f.close()
			print "      written to file " + out + ".files.incomplete"
		if len(reg_missingfile) > 0:
			print "   " + str(len(reg_missingfile)) + " of " + str(n) + " files missing"
			f = open(out + '.files.missing', 'w')
			f.write('\n'.join(reg_missingfile))
			f.close()
			print "      written to file " + out + ".files.missing"
		if len(reg_nomarkers) > 0:
			print "   " + str(len(reg_nomarkers)) + " of " + str(n) + " contained no markers"
			f = open(out + '.no_markers', 'w')
			f.write('\n'.join(reg_nomarkers))
			f.close()
			print "      written to file " + out + ".no_markers"
		if len(reg_rerun) > 0:
			print "   " + str(len(reg_rerun)) + " of " + str(n) + " total regions to rerun"
			f = open(out + '.regions.rerun', 'w')
			f.write('\n'.join(reg_rerun))
			f.close()
			print "      written to file " + out + ".regions.rerun"
		if len(reg_complete) + len(reg_nomarkers) != n:
			return False, reg_complete
		else:
			return True, reg_complete
	else:
		reg = file_dict.values()[0].split(":")[0] + ':' + file_dict.values()[0].split(":")[1].split("-")[0] + '-' + file_dict.values()[0].split(":")[1].split("-")[1]
		resfile = file_basename + '.chr' + file_dict.values()[0].split(":")[0] + 'bp' + file_dict.values()[0].split(":")[1].split("-")[0] + '-' + file_dict.values()[0].split(":")[1].split("-")[1]
		logfile = resfile + ".log"
		resfile = resfile + ".gz"
		if os.path.exists(resfile):
			if os.path.exists(logfile):
				p = subprocess.Popen(['grep','-cw','process complete',logfile], stdout=subprocess.PIPE)	
				complete = p.communicate()[0]
				complete = int(complete.strip())
				if complete == 0:
					p = subprocess.Popen(['grep','-cw','no markers found',logfile[0]], stdout=subprocess.PIPE)	
					no_makers = p.communicate()[0]
					no_makers = int(no_makers.strip())
					if no_makers != 0:
						print "analysis complete"
						return True, reg
					else:
						print "analysis incomplete"
						return False, reg
				else:
					print "analysis incomplete"
					return False, reg
			else:
				print "analysis incomplete ... no log file found"
				return False, reg
		else:
			print "analysis incomplete ... no out file found"
			return False, reg

def CompileResults(out_files, out):
	print "compiling results"
	bgzfile = bgzf.BgzfWriter(out + '.gz', 'wb')
	logfile = file(out + '.log', 'w')
	pbar = ProgressBar(maxval=len(out_files.keys()), widgets = ['   processed ', Counter(), ' of ' + str(len(out_files.keys())) + ' regions (', Timer(), ')'])
	i=0
	pbar.start()
	for reg in out_files.keys():
		i = i + 1
		sed = ['awk','{print $0}'] if i == 1 else ['sed','1d']
		resfile = out_files[reg]
		resfile = resfile + ".gz"
		p1 = subprocess.Popen(['zcat',resfile], stdout=subprocess.PIPE, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
		p2 = subprocess.Popen(sed, stdin=p1.stdout, stdout=subprocess.PIPE)
		bgzfile.write(p2.communicate()[0])
		p1.wait()
		p2.wait()
		outlog = glob.glob(out_files[reg] + ".*.log")
		if i == 1:
			p3 = subprocess.Popen(['cat',outlog[0]], stdout=subprocess.PIPE)
			p4 = subprocess.Popen(['awk','{print \"      \"$0}'], stdin=p3.stdout, stdout=subprocess.PIPE)
			logfile.write('Sample log file from ' + outlog[0] + '\n\n')
			logfile.write(p4.communicate()[0])
			p3.wait()
			p4.wait()
			logfile.write('\nElapsed time and max memory used for each job in list\n\n')
		p5 = subprocess.Popen(['grep','time elapsed:',outlog[0]], stdout=subprocess.PIPE)
		elap = p5.communicate()[0].strip()
		p5.wait()
		p6 = subprocess.Popen(['grep','memory used:',outlog[0]], stdout=subprocess.PIPE)
		maxmem = p6.communicate()[0].strip()
		p6.wait()
		logfile.write('      job ' + str(reg) + '   -   ' + elap + '   ' + maxmem + '\n')
		pbar.update(i)
	pbar.finish()
	bgzfile.close()
	logfile.close()
	print "mapping compiled file"
	p7 = subprocess.Popen(['zcat',out + '.gz'], stdout=subprocess.PIPE, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
	p8 = subprocess.Popen(['head','-1'], stdin=p7.stdout, stdout=subprocess.PIPE)
	header = p8.communicate()[0]
	header = header.split()

	if 'pos' in header:
		b = header.index('pos') + 1
		e = b
	elif 'chr.pos' in header:
		b = header.index('chr.pos') + 1
		e = b
	else:
		b = header.index('start') + 1
		e = header.index('end') + 1
	cmd = 'tabix -b ' + str(b) + ' -e ' + str(e) + ' ' + out + '.gz'
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print "file compilation complete"
	return True
