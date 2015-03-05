from plinkio import plinkfile
import pysam
import vcf as VCF
import tabix
import gzip
	
def LoadPlink(data):
	plink_bed_it = plinkfile.open(data)
	plink_sample_it = plink_bed_it.get_samples()
	plink_locus_it = plink_bed_it.get_loci()
	sample_ids = [x.iid for x in plink_sample_it]
	return plink_bed_it, plink_locus_it, plink_sample_it, sample_ids

def LoadVcf(data):
	tb = tabix.open(data)
	v=VCF.Reader(filename=data)
	record=v.next()
	sample_ids = [a.sample for a in record.samples]
	return tb, sample_ids

def LoadDos(data, samples):
	tb = tabix.open(data)
	sample_ids = []
	open_func = gzip.open if samples[-3:] == '.gz' else open
	with open_func(samples) as sf:
		lines = (line.rstrip() for line in sf)
		lines = (line for line in lines if line)
		for line in lines:
			sample_ids.append(line)
	return tb, sample_ids
