import subprocess
import tabix
import math
from Messages import Error
from multiprocessing import Process, Manager, cpu_count

def Map(file, out, cpu, mb = None, kb = None, b = None):

	##### list possible regions #####
	s = int(b) if b else ''
	s = int(kb) * 1000 if kb else s
	s = int(mb) * 1000000 if mb else s
	starts = range(1,300000000,s)
	regions = []
	for i in range(1,27):
		for rp in starts:
			regions.append(str(i) + ":" + str(rp) + "-" + str(rp+s-1)) 
	
	##### start multiprocessing jobs #####
	cpus = cpu_count() if cpu_count() < int(cpu) else int(cpu)
	manager = Manager()
	return_dict = manager.dict()
	jobs = []
	
	##### define function to check chunk of regions for markers #####
	def CheckRegion(i, joblist, return_dict):
		for x in joblist:
			found = False
			tb = tabix.open(file)
			try:
				records = tb.querys(x)
			except:
				pass
			else:
				for record in records:
					found = True
					break
				return_dict[x]=found
		print "   ...processed " + str(len(joblist)) + " regions on cpu " + str(i+1)

	##### submit parallel processing jobs on available cores #####
	joblist = [regions[i:i+int(math.ceil(float(len(regions))/cpus))] for i in range(0,len(regions),int(math.ceil(float(len(regions))/cpus)))]

	for i in range(len(joblist)):
		print "   ...submitting jobs on cpu " + str(i+1)
		p = Process(target=CheckRegion, args=(i, joblist[i], return_dict))
		jobs.append(p)
		p.start()
	for job in jobs:
		job.join()
	
	r_write = []
	for r in return_dict.keys():
		if(return_dict[r]):
			r_write.append(r)
	r_write.sort(key=lambda a: regions.index(a))
	
	with open(out, "w") as fout:
		for i in r_write:
			fout.write(i + "\n")
	
	if b:
		print "   " + str(len(r_write)) + " non-empty regions of size " + str(s) + "bp written to file " + out
	elif kb:
		print "   " + str(len(r_write)) + " non-empty regions of size " + str(s/1000) + "kb written to file " + out
	elif mb:
		print "   " + str(len(r_write)) + " non-empty regions of size " + str(s/1000000) + "mb written to file " + out
