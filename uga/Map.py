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
from multiprocessing import Process, Manager, cpu_count
from itertools import islice,groupby
from operator import attrgetter
import pysam
import glob
import os
from progressbar import ProgressBar, Counter, Timer

def map(file, 
		region = None, 
		mb = '1', 
		shift_mb = None):

	s = int(float(mb) * 1000000) if mb else s
	shift = int(shift_mb) * 1000000 if shift_mb else 0

	regions_out = []

	files_glob = glob.glob(file.replace('[CHR]','*'))
	for f in files_glob:
		print "mapping file " + os.path.basename(f)
		v = pysam.TabixFile(filename=f, parser=pysam.asTuple())
		regions = []
		start = 1
		end = 1000000000
		if not region is None:
			if len(region.split(':')) > 1:
				start = region.split(':')[1].split('-')[0]
				end = region.split(':')[1].split('-')[1]
			if start != end:
				starts = range(int(start),int(end),s-shift)
			else:
				starts = [int(start)]
			for rp in starts:
				regions.append(region.split(':')[0] + ":" + str(rp) + "-" + str(min(rp+s-1,int(end))))
		else:
			starts = range(start,end,s-shift)
			for chr in v.contigs:
				for rp in starts:
					regions.append(chr + ":" + str(rp) + "-" + str(rp+s-1))
		prev_chr=''
		for reg in regions:
			chr = reg.split(':')[0]
			if chr != prev_chr:
				print 'mapping chromosome ' + chr
			prev_chr = chr
			start = reg.split(':')[1].split('-')[0]
			end = reg.split(':')[1].split('-')[1]
			if not s is None:
				try:
					records = v.fetch(region=reg, parser=pysam.asTuple())
				except:
					pass
				else:
					for record in records:
						if int(record[1]) >= int(start) and int(record[1]) <= int(end):
							regions_out.append(reg)
							break
			else:
				try:
					records = v.fetch(region=reg, parser=pysam.asTuple())
				except:
					pass
				else:
					last=0
					j=0
					while True:
						j += 1
						if j == 1:
							for record in islice(records, 0,n):
								pass
							regions_out.append(chr + ':' + start + '-' + record[1])
							last = int(record[1])
						else:
							chunk=tuple(islice(records, n - 1,n))
							if not chunk:
								break
							regions_out.append(chr + ':' + str(last+1) + '-' + chunk[0][1])
							last = int(chunk[0][1])
					if last+1 < end:
						records = v.fetch(region=chr + ':' + str(last+1) + '-' + end, parser=pysam.asTuple())
						for record in records:
							regions_out.append(chr + ':' + start + '-' + record[1])
							break

	return regions_out
