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
	

def countPlinkRegion(iter, chr, start, end, ml = None):
	k = 0
	if ml is not None:
		filt = ifilter(lambda c: ml[(ml['chr'] == c.chromosome) & (ml['pos'] == c.bp_position) & (ml['a1'] == c.allele1) & (ml['a2'] == c.allele2)].shape[0] > 0, iter)
	else:
		filt = ifilter(lambda c: c.chromosome == chr and c.bp_position >= start and c.bp_position <= end, iter)
	for record in filt:
		k += 1
	return k

def countNonPlinkRegion(iter, type, start, end, ml = None):
	k = 0
	if type in ['vcf','dos2']:
		if ml is not None:
			filt = ifilter(lambda c: ml[(ml['chr'] == int(c[0])) & (ml['pos'] == int(c[1])) & (ml['a1'] == c[3]) & (ml['a2'] == c[4])].shape[0] > 0, iter)
		else:
			filt = ifilter(lambda c: int(c[1]) >= start and int(c[1]) <= end, iter)
		for record in filt:
			k += 1
	else:
		if ml is not None:
			filt = ifilter(lambda c: ml[(ml['chr'] == int(c[0])) & (ml['pos'] == int(c[2])) & (ml['a1'] == c[2]) & (ml['a2'] == c[4])].shape[0] > 0, iter)
		else:
			ifilter(lambda c: int(c[2]) >= start and int(c[2]) <= end, iter)
		for record in filt:
			k += 1
	return k

def ConvertDosage(row):
	newrow = row[:5]
	a=zip(row[5::3],row[6::3],row[7::3])
	newrow = newrow + [2*float(t[0]) + 1*float(t[1]) if float(t[0]) > 0 or float(t[1]) > 0 or float(t[2]) > 0 else float('nan') for t in a]
	return newrow

def GetFocus(model_fxn,model,vars_df,model_vars_dict):
	focus = ['(Intercept)'] if model_fxn.split('_')[0] in ['gee','geeboss','glm','lme'] else []
	for x in re.sub("\+|\~|\-",",",model.split('~')[1]).split(','):
		if not x[0] == '(' and x.find('cluster(') == -1:
			if x.find('factor(') != -1:
				x = re.sub('factor\(|\)','',x)
			if x.find('*') != -1:
				interact = []
				for a in x.split('*'):
					if model_vars_dict[a]['class'] == 'factor':
						for v in vars_df[a].unique().astype(np.int64) - min(vars_df[a].unique().astype(np.int64)):
							if v != min(vars_df[a].unique().astype(np.int64) - min(vars_df[a].unique().astype(np.int64))):
								interact.append('factor(' + a + ')' + str(v))
					else:
						interact.append(a)
				focus.append(interact[0] + '*' + interact[1])
			elif model_vars_dict[x]['class'] == 'factor':
				for v in vars_df[x].unique().astype(np.int64) - min(vars_df[x].unique().astype(np.int64)):
					if v != min(vars_df[x].unique().astype(np.int64) - min(vars_df[x].unique().astype(np.int64))):
						focus.append('factor(' + x + ')' + str(v))
			elif x.find('*') != -1:
				focus.append('*'.join(sorted(x.split('*'))))
			else:
				focus.append(x)
	return focus

def GenerateFilterCode(marker_info, no_mono=True, miss = None, freq = None, max_freq = None, rsq = None, hwe = None):
	filter = 0
	if (not miss is None and not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) < float(miss)) or (not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) == 0) or (math.isnan(marker_info['callrate'])):
		filter += 1000
	if not math.isnan(marker_info['freq']): 
		if no_mono and (float(marker_info['freq']) == 0 or float(marker_info['freq']) == 1):
			filter += 100
		else:
			if ((	not freq is None
				and 
					(		float(marker_info['freq']) < float(freq)
						 or float(marker_info['freq']) > 1-float(freq)
					)
				and 
					(		float(marker_info['freq.unrel']) < float(freq)
						 or float(marker_info['freq.unrel']) > 1-float(freq)
					)
			   ) 
			  or
			   (	not max_freq is None
				and 
					(	   
						(		float(marker_info['freq']) >= float(max_freq)
							and float(marker_info['freq']) <= 1-float(max_freq)
						)
					)
				and
					(
						(		float(marker_info['freq.unrel']) >= float(max_freq)
							and float(marker_info['freq.unrel']) <= 1-float(max_freq)
						)
						or float(marker_info['freq.unrel']) == 0
						or float(marker_info['freq.unrel']) == 1
					)
			   )):
				filter += 100
	if not rsq is None and not math.isnan(marker_info['rsq']) and (float(marker_info['rsq']) < float(rsq) and float(marker_info['rsq.unrel']) < float(rsq)):
		filter += 10
	if not hwe is None and not math.isnan(marker_info['hwe']) and (float(marker_info['hwe']) < float(hwe) and float(marker_info['hwe.unrel']) < float(hwe)):
		filter += 1
	return filter
