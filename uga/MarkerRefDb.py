from MarkerCalc import ListCompatibleMarkers,Complement
from multi_key_dict import multi_key_dict
import numpy as np
import pandas as pd

class MarkerRefDb(object):

	def __init__(self):
		self.ref = multi_key_dict()
		self.ref_alleles = multi_key_dict()

	def Update(self, row):
		newrow=row
		chr=newrow['chr']
		pos=newrow['pos']
		marker=newrow['marker']
		a1=newrow['a1']
		a2=newrow['a2']
		record_markers = ListCompatibleMarkers(chr,pos,a1,a2,'><')
		if len([r for r in record_markers if self.ref.has_key(r)]) > 0:
			temp_ref = self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			temp_alleles = self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			temp_alleles[1] = ','.join(list(set(','.join([b.split('><')[3] for b in record_markers if b.split('><')[2] == temp_alleles[0]]).split(','))))
			del self.ref[[r for r in record_markers if self.ref.has_key(r)][1]]
			del self.ref_alleles[[r for r in record_markers if self.ref_alleles.has_key(r)][1]]
			if not marker in temp_ref:
				temp_ref.append(marker)
				temp_ref = list(set(temp_ref))
			self.ref[tuple(record_markers)] = temp_ref
			self.ref_alleles[tuple(record_markers)] = temp_alleles
		else:
			self.ref[record_markers] = [marker]
			self.ref_alleles[record_markers] = [a1,a2]
		refmarker = filter(lambda x:'rs' in x, self.ref[chr + '><' + pos + '><' + a1 + '><' + a2])
		if len(refmarker) > 0:
			refmarker = refmarker[0]
		else:
			refmarker = self.ref[chr + '><' + pos + '><' + a1 + '><' + a2][0]
		keymarkers = self.ref[chr + '><' + pos + '><' + a1 + '><' + a2]
		refa1 = self.ref_alleles[chr + '><' + pos + '><' + a1 + '><' + a2][0]
		refa2 = self.ref_alleles[chr + '><' + pos + '><' + a1 + '><' + a2][1]
		if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
			newrow[5:] = 1-newrow[5:]
		newrow['marker']=refmarker	
		newrow['a1']=refa1
		newrow['a2']=refa2
		return newrow
					