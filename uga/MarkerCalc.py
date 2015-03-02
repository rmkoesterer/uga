def Complement(x):
	if x != "NA":
		letters = list(x)
	else:
		letters = "NA"
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
	return ''.join(comp)

def ListCompatibleMarkers(chr,pos,a1,a2,delim):
	markers = []
	markers.append(chr + delim + pos + delim + a1 + delim + a2)
	markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(a2))
	markers.append(chr + delim + pos + delim + Complement(a2) + delim + Complement(a1))
	markers.append(chr + delim + pos + delim + a2 + delim + a1)
	for alt in a2.split(','):
		markers.append(chr + delim + pos + delim + a1 + delim + alt)
		markers.append(chr + delim + pos + delim + Complement(a1) + delim + Complement(alt))
		markers.append(chr + delim + pos + delim + Complement(alt) + delim + Complement(a1))
		markers.append(chr + delim + pos + delim + alt + delim + a1)
	return list(set(markers))

def ConvertDosage(x):
	a=zip(x[0::3],x[1::3],x[2::3])
	return [2*t[0] + 1*t[1] if t[0] > 0 or t[1] > 0 or t[2] > 0 else float('nan') for t in a]
	
def CalcCallrate(x):
	xlen = len(x)
	x = x.dropna().astype(float)
	if len(x) > 0:
		return float('%.5g' % (len(x)/float(xlen)))
	else:
		return 0.0

def CalcFreq(marker, chr, male_idx = None, female_idx = None):
	marker = marker.dropna().astype(float)
	if int(chr) == 23:
		if not male_idx is None and not female_idx is None:
			n = len(marker.loc[male_idx]) + 2 * len(marker.loc[female_idx])
			count = (marker.loc[male_idx].sum()/2) + marker.loc[female_idx].sum()
			return float('%.5g' % (count / n)) if len(marker.loc[male_idx]) > 0 and len(marker.loc[female_idx]) > 0 else float('NaN')
		else:
			float('NaN')
	else:
		n = 2 * len(marker)
		count = marker.sum()
		return float('%.5g' % (count / n)) if len(marker) > 0 else float('NaN')

def CalcRsq(x):
	x = x.dropna().astype(float)
	rsq = x.var() / (2 * (x.mean() / 2) * (1 - (x.mean() / 2))) if (len(x) > 0 and x.mean()/2 != 1 and x.mean()/2 != 0) else None
	if rsq != None: 
		rsq = 1 / rsq if rsq > 1 else rsq
	return float('%.5g' % (rsq)) if rsq != None else float('NaN')

def CalcHWE(marker, chr = None, female_idx = None):
	marker = marker.dropna().astype(float)
	if list(set(marker)) == [0,1,2]:
		if not chr is None and int(chr) == 23:
			if not female_idx is None:
				obs_hets = len(marker.loc[female_idx][marker.loc[female_idx].isin([1])])
				obs_hom1 = len(marker.loc[female_idx][marker.loc[female_idx].isin([2])])
				obs_hom2 = len(marker.loc[female_idx][marker.loc[female_idx].isin([0])])
			else:
				obs_hets = 0
				obs_hom1 = 0
				obs_hom2 = 0
		else:
			obs_hets = len(marker[marker.isin([1])])
			obs_hom1 = len(marker[marker.isin([2])])
			obs_hom2 = len(marker[marker.isin([0])])	
		if not (obs_hets == 0 and obs_hom1 == 0 and obs_hom2 == 0):
			obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
			obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

			rare_copies = 2 * obs_homr + obs_hets
			genotypes   = obs_hets + obs_homc + obs_homr

			het_probs = [0.0] * (rare_copies + 1)

			#start at midpoint
			mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)

			#check to ensure that midpoint and rare alleles have same parity
			if (rare_copies & 1) ^ (mid & 1):
				mid += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			het_probs[mid] = 1.0
			sum = float(het_probs[mid])

			for curr_hets in xrange(mid, 1, -2):
				het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
				sum += het_probs[curr_hets - 2]
				# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
				curr_homr += 1
				curr_homc += 1

			curr_hets = mid
			curr_homr = (rare_copies - mid) / 2
			curr_homc = genotypes - curr_hets - curr_homr

			for curr_hets in xrange(mid, rare_copies - 1, 2):
				het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
				sum += het_probs[curr_hets + 2]
				#add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
				curr_homr -= 1
				curr_homc -= 1

			for i in xrange(0, rare_copies + 1):
				het_probs[i] /= sum

			#alternate p-value calculation for p_hi/p_lo
			p_hi = float(het_probs[obs_hets])
			for i in xrange(obs_hets, rare_copies+1):
				p_hi += het_probs[i]

			p_lo = float(het_probs[obs_hets])
			for i in xrange(obs_hets-1, -1, -1):
				p_lo += het_probs[i]

			p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

			p_hwe = 0.0
			#  p-value calculation for p_hwe
			for i in xrange(0, rare_copies + 1):
				if het_probs[i] > het_probs[obs_hets]:
					continue
				p_hwe += het_probs[i]

			p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

			return float('%.2e' % (p_hwe))
		else:
			return float('NaN')
	else:
		return float('NaN')

def CallToDos(call):
	c = call.split(':')[0]
	c = 'NaN' if not c in ['0/0','0/1','1/1','1/0'] else c
	return float(c.replace('0/0','2').replace('0/1','1').replace('1/1','0').replace('1/0','1'))