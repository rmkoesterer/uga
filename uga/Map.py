from __main__ import *
from multiprocessing import Process, Manager, cpu_count
from plinkio import plinkfile
from itertools import islice,groupby
from operator import attrgetter
import vcf as VCF

def Map(out, 
		oxford = None, 
		dos1 = None, 
		dos2 = None, 
		plink = None, 
		vcf = None, 
		region = None, 
		chr = None, 
		mb = None, 
		kb = None, 
		b = None, 
		n = None, 
		shift_mb = None,
		shift_kb = None, 
		shift_b = None):

	assert not oxford is None or not dos1 is None or not dos2 is None or not plink is None or not vcf is None, SystemFxns.Error("a genotype data file must be specified")
	assert not b is None or not kb is None or not mb is None or not n is None, SystemFxns.Error("a region size or number of markers must be specified")
	
	s = int(b) if b else None
	s = int(kb) * 1000 if kb else s
	s = int(mb) * 1000000 if mb else s
	shift = int(shift_b) if shift_b else 0
	shift = int(shift_kb) * 1000 if shift_kb else shift
	shift = int(shift_mb) * 1000000 if shift_mb else shift
	if shift == None:
		shift = 0

	if not n is None:
		n = int(n)

	start = 1
	end = 300000000
	if not chr is None:
		if not int(chr) in range(1,27):
			print SystemFxns.Error("chromosome " + str(chr) + " is an invalid choice")
			return
		chrs = [int(chr)]
	elif not region is None:
		chrs = [int(region.split(':')[0])]
		start = int(region.split(':')[1].split('-')[0])
		end = int(region.split(':')[1].split('-')[1])
	else:
		chrs = range(1,27)

	def chunk(seq, size):
		return [seq[pos:pos + size] for pos in xrange(0, len(seq), size)]
	
	##### map bim file #####
	if not plink is None:
		total = 0
		print "mapping Plink format bim file"
		bim=pd.read_table(plink + '.bim',header=None,names=['chr','marker','x','pos','a1','a2'])
		with open(out, "w") as fout:
			for i in chrs:
				if i in bim['chr'].values:
					bim_chr=bim[bim['chr'] == i]
					bim_chr.reset_index(inplace=True,drop=True)
					if not s is None:
						chr_write=pd.DataFrame({'chr': [i for x in range(start,min(bim_chr.iloc[-1]['pos']+1,end),s-shift)], 'start': range(start,min(bim_chr.iloc[-1]['pos']+1,end),s-shift), 'end': range(s+start,min(bim_chr.iloc[-1]['pos']+s+1,end+s),s-shift)})
						chr_write['reg']=chr_write.apply(lambda a: str(a['chr']) + ':' + str(a['start']) + '-' + str(a['end']-1),axis=1)
						for j in chr_write.index:
							if bim[(bim['chr'] == i) & (bim['pos'] >= chr_write.iloc[j]['start']) & (bim['pos'] <= chr_write.iloc[j]['end'])].shape[0] > 0:
								total += 1
								fout.write(chr_write.iloc[j]['reg'] + "\n")
					else:
						bim_chr = bim_chr[(bim_chr['pos'] >= start) & (bim_chr['pos'] <= end)]
						chr_write=pd.DataFrame({'chr': [i for x in range(bim_chr.iloc[::n,3].shape[0])], 'start': np.array(bim_chr.iloc[::n,3]), 'end': np.append(np.array(bim_chr.iloc[::n,3])[1:],np.array(bim_chr['pos'])[-1])})
						chr_write['reg']=chr_write.apply(lambda a: str(a['chr']) + ':' + str(a['start']) + '-' + str(a['end']-1),axis=1)
						chr_write['reg'].to_csv(fout,header=False,index=False,mode='a')
						total += len(chr_write['reg'])
					print "   processed chromosome " + str(i)
	elif not vcf is None:
		total = 0
		with open(out, "w") as fout:
			for i in chrs:
				if not s is None:
					v=VCF.Reader(filename=vcf)
					starts = range(start,end,s-shift)
					regions = []
					for rp in starts:
						regions.append(str(rp) + "-" + str(rp+s-1))
					for reg in regions:
						try:
							records = v.fetch(str(i),int(reg.split('-')[0]),int(reg.split('-')[1]))
						except:
							pass
						else:
							for record in records:
								if record.POS >= int(reg.split('-')[0]) and record.POS <= int(reg.split('-')[1]):
									total += 1
									fout.write(str(i) + ':' + reg + '\n')
									break
					print "   processed chromosome " + str(i)
				else:
					tb = tabix.open(vcf)
					try:
						if not region is None:
							records = tb.querys(region)
						else:
							records = tb.querys(str(i))
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
								total += 1
								fout.write(str(i) + ':' + str(start) + '-' + str(record[1]) + '\n')
								last = int(record[1])
							else:
								chunk=tuple(islice(records, n - 1,n))
								if not chunk:
									break
								total += 1
								fout.write(str(i) + ':' + str(last+1) + '-' + str(chunk[0][1]) + '\n')
								last = int(chunk[0][1])
						if last+1 < end:
							records = tb.querys(str(i) + ':' + str(last+1) + '-' + str(end))
							for record in records:
								total += 1
								fout.write(str(i) + ':' + str(last+1) + '-' + str(end) + '\n')
								break
					print "   processed chromosome " + str(i)
	else:
		if not oxford is None:
			file = oxford
		elif not dos1 is None:
			file = dos1
		else:
			file = dos2
		total = 0
		tb = tabix.open(file)
		with open(out, "w") as fout:
			for i in chrs:
				if not s is None:
					starts = range(start,end,s-shift)
					regions = []
					for rp in starts:
						regions.append(str(i) + ":" + str(rp) + "-" + str(rp+s-1)) 
					for reg in regions:
						try:
							records = tb.querys(reg)
						except:
							pass
						else:
							for record in records:
								total += 1
								fout.write(reg + '\n')
								break
				else:
					try:
						if not region is None:
							records = tb.querys(region)
						else:
							records = tb.querys(str(i))
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
								total += 1
								fout.write(str(i) + ':' + str(start) + '-' + str(record[1]) + '\n')
								last = int(record[2]) if not oxford is None or not dos1 is None else int(record[1])
							else:
								chunk=tuple(islice(records, n - 1,n))
								if not chunk:
									break
								total += 1
								end_temp = int(chunk[0][2]) if not oxford is None or not dos1 is None else int(chunk[0][1])
								fout.write(str(i) + ':' + str(last+1) + '-' + str(end_temp) + '\n')
								last = end_temp
						if last+1 < end:
							records = tb.querys(str(i) + ':' + str(last+1) + '-' + str(end))
							for record in records:
								total += 1
								fout.write(str(i) + ':' + str(last+1) + '-' + str(end) + '\n')
								break
				print "   processed chromosome " + str(i)
	if b:
		print str(total) + " non-empty regions of size " + str(s) + "bp written to file " + out
	elif kb:
		print str(total) + " non-empty regions of size " + str(s/1000) + "kb written to file " + out
	elif mb:
		print str(total) + " non-empty regions of size " + str(s/1000000) + "mb written to file " + out
	elif n:
		print str(total) + " non-empty regions of length " + str(n) + " written to file " + out
