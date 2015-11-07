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
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def Complement(allele):
	cdef str x = allele
	if x != "NA":
		letters = list(x)
		comp = []
		for l in letters:
			if l == 'T': 
				c = 'A'
			elif l == 'A':
				c = 'T'
			elif l == 'G': 
				c = 'C'
			elif l == 'C':
				c = 'G'
			elif l == '0':
				c = '0'
			elif l == ',':
				c = ','
			elif l == 'NA':
				c = 'NA'
			elif l == '-':
				c = '-'
			elif l == 'I':
				c = 'D'
			elif l == 'D':
				c = 'I'
			elif l in ['1','2','3','4','5','6','7','8','9','0']:
				c = l
			else:
				c = 'X'
			comp.append(c)
	else:
		comp = ['NA']
	return ''.join(comp)

@cython.boundscheck(False)
@cython.wraparound(False)
def ListCompatibleMarkers(chr_py,pos_py,a1_py,a2_py,delim_py):
	cdef str chr = chr_py
	cdef str pos = pos_py
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef str delim = delim_py
	analogs = [chr + delim + pos + delim + a1 + delim + a2, 
					chr + delim + pos + delim + Complement(a1) + delim + Complement(a2)]
	if a2 != 'NA':	
		analogs = analogs + [chr + delim + pos + delim + Complement(a2) + delim + Complement(a1),
							chr + delim + pos + delim + a2 + delim + a1, 
							chr + delim + pos + delim + a1 + delim + 'NA', 
							chr + delim + pos + delim + Complement(a1) + delim + 'NA', 
							chr + delim + pos + delim + Complement(a2) + delim + 'NA', 
							chr + delim + pos + delim + a2 + delim + 'NA']
	return sorted(list(set(analogs)))

@cython.boundscheck(False)
@cython.wraparound(False)
def ListCompatibleMarkersMeta(chr_py,pos_py,a1_py,a2_py,delim_py):
	cdef str chr = chr_py
	cdef str pos = pos_py
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef str delim = delim_py
	analogs = [chr + delim + pos + delim + a1 + delim + a2, 
					chr + delim + pos + delim + Complement(a1) + delim + Complement(a2)]
	if a2 != 'NA':	
		analogs = analogs + [chr + delim + pos + delim + Complement(a2) + delim + Complement(a1),
							chr + delim + pos + delim + a2 + delim + a1, 
							chr + delim + pos + delim + a1 + delim + 'NA', 
							chr + delim + pos + delim + Complement(a1) + delim + 'NA', 
							chr + delim + pos + delim + Complement(a2) + delim + 'NA', 
							chr + delim + pos + delim + a2 + delim + 'NA']
	return "_".join(sorted(list(set(analogs))))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipEffect(refa1_py, refa2_py, a1_py, a2_py, effect_py):
	cdef str refa1 = refa1_py[0:20]
	cdef str refa2 = refa2_py[0:20]
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef float effect = effect_py
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * effect
	else:
		return effect

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipFreq(refa1_py, refa2_py, a1_py, a2_py, freq_py):
	cdef str refa1 = refa1_py[0:20]
	cdef str refa2 = refa2_py[0:20]
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef float freq = freq_py
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 - freq
	else:
		return freq

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipOR(refa1_py, refa2_py, a1_py, a2_py, o_r_py):
	cdef str refa1 = refa1_py[0:20]
	cdef str refa2 = refa2_py[0:20]
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef float o_r = o_r_py
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return 1 / o_r
	else:
		return o_r

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def FlipZ(refa1_py, refa2_py, a1_py, a2_py, z_py):
	cdef str refa1 = refa1_py[0:20]
	cdef str refa2 = refa2_py[0:20]
	cdef str a1 = a1_py[0:20]
	cdef str a2 = a2_py[0:20]
	cdef float z = z_py
	if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
		return -1 * z
	else:
		return z

def ChunkFlip(row, vdb):
	vdb_match = vdb[vdb.analogs.map(lambda x: set(x).issubset(set(row['analogs'])) or set(row['analogs']).issubset(set(x)))]
	if vdb_match.shape[0] > 0:
		a1=row['a1']
		a2=row['a2']
		refa1=vdb_match['a1'][0]
		refa2=vdb_match['a2'][0]
		row['a1'] = refa1
		row['a2'] = refa2 if refa2 != 'NA' else a2
		if (refa1 + refa2 == Complement(a2) + Complement(a1) or refa1 + refa2 == a2 + a1 or refa1 + refa2 == Complement(a2) + 'NA' or refa1 + refa2 == a2 + 'NA') and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
			if refa1 + refa2 == Complement(a2) + 'NA':
				row['a2'] = Complement(a1)
			if refa1 + refa2 == a2 + 'NA':	
				row['a2'] = a1
			row[5:len(row)-1] = 2-row[5:len(row)-1].str.replace('NA','nan').astype(float)
	return row

def UpdateAltAllele(row, vdb):
	vdb_match = vdb[vdb.analogs.map(lambda x: '><'.join([row['chr'],row['pos'],row['a1'],row['a2']]) in x)]
	if vdb_match.shape[0] > 0:
		row['a2'] = vdb_match['a2'][0]
	return row['a2']

def VdbUpdate(row, chunkdf):
	chunkdf_match = chunkdf[chunkdf.analogs.map(lambda x: set(x).issubset(set(row['analogs'])) or set(row['analogs']).issubset(set(x)))]
	if chunkdf_match.shape[0] > 0:
		if row['a2'] == 'NA':
			row['a2'] = chunkdf_match['a2'][0]
		row['analogs'] = sorted(list(set(row['analogs'] + chunkdf_match['analogs'][0])))
	return row

def VdbAppend(vdb, chunkdf):
	vdb['check'] = vdb['analogs'].apply(lambda x: '_'.join(x),1)
	chunkdf['check'] = chunkdf['analogs'].apply(lambda x: '_'.join(x),1)
	chunkdf_nomatch = chunkdf[~chunkdf['check'].isin(vdb['check'])]
	if chunkdf_nomatch.shape[0] > 0:
		vdb = vdb.append(pd.DataFrame({'chr': list(chunkdf_nomatch['chr']),'pos': list(chunkdf_nomatch['pos']),'marker': list(chunkdf_nomatch['marker']),'a1': chunkdf_nomatch['a1'],'a2': chunkdf_nomatch['a2'], 
													'analogs': list(chunkdf_nomatch['analogs'])}))
	return vdb
