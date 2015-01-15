import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import statsmodels.api as sm
import pandas.rpy.common as py2r
import math

def GenerateFilterCode(marker_info, miss = None, freq = None, rsq = None, hwe = None):
	filter = 0
	if not miss is None and not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) < miss:
		filter += 100
	if not freq is None and not math.isnan(marker_info['freq']) and ((float(marker_info['freq']) < freq or float(marker_info['freq']) > 1-freq or float(marker_info['freq']) == 0 or float(marker_info['freq']) == 1) and (float(marker_info['freq.unrel']) < freq or float(marker_info['freq.unrel']) > 1-freq or float(marker_info['freq.unrel']) == 0 or float(marker_info['freq.unrel']) == 1)):
		filter += 10
	if not rsq is None and not math.isnan(marker_info['rsq']) and (float(marker_info['rsq']) < rsq and float(marker_info['rsq.unrel']) < rsq):
		filter += 1
	if not hwe is None and not math.isnan(marker_info['hwe']) and (float(marker_info['hwe']) < hwe and float(marker_info['hwe.unrel']) < hwe):
		filter += 1
	return filter

"""def CalcGEEPython(marker_info, model_df, model_vars_dict, MODEL, IID_COL, FID_COL, METHOD, FXN_FAMILY, FOCUS_VARS, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	#if FXN_FAMILY == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(rmodel_df.rx2(marker_info['marker_unique'])).labels) < 3 and 0 in rxtab(robjects.r('~' + marker_info['marker_unique'] + '+' + dep_var), data=rmodel_df))):
	if FXN_FAMILY == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df[marker_info['marker_unique']]).labels) < 3 and 0 in pd.crosstab(model_df[marker_info['marker_unique']],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] != 0:
			model = sm.GEE.from_formula(MODEL, groups=FID_COL, data=model_df, family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Exchangeable())
			result = model.fit()
	print marker_info['marker_unique']
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info"""

def CalcGEE(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).labels) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))])
			model_out=rtry(rsummary(geepack.geeglm(ro.r(model),data=rmodel_df,id=rmodel_df.rx2(fid),family=fxn,corstr='exchangeable')),silent=ro.r('TRUE'))
			if 'try-error' in rclass(model_out):
				model_out=rtry(rsummary(geepack.geeglm(ro.r(model),data=rmodel_df,id=rmodel_df.rx2(fid),family=fxn,corstr='independence')),silent=ro.r('TRUE'))
				if 'try-error' in rclass(model_out):
					status = -4
				else:
					status = 2 if model_out.rx2('error')[0] != 1 else -2
					valid = True if model_out.rx2('error')[0] != 1 else False
			else:
				status = 1 if model_out.rx2('error')[0] != 1 else -1
				valid = True if model_out.rx2('error')[0] != 1 else False
		else:
			status = -5
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':')
			marker_info[x + '.beta'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std.err'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate']))
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err'])
			marker_info[x + '.p'] = '%.2e' % (2 * scipy.norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err'])))
	else:
		for x in focus:
			marker_info[x + '.beta'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info.drop('marker_unique',axis=1,inplace=True)

def CalcGLM(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).labels) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))])
			model_out=rtry(rsummary(rglm(ro.r(model),data=rmodel_df,family=fxn)),silent=ro.r('TRUE'))
			if 'try-error' in rclass(model_out):
				status = -4
			else:
				status = 1
				valid = True
		else:
			status = -5
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':')
			marker_info[x + '.beta'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate']))
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error'])
			marker_info[x + '.p'] = '%.2e' % (2 * scipy.norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error'])))
	else:
		for x in focus:
			marker_info[x + '.beta'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info.drop('marker_unique',axis=1,inplace=True)

def CalcLME(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).labels) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))])
			model_out=rtry(rsummary(rlme4.glmer(ro.r(model),data=rmodel_df,REML=ro.r('FALSE'),family=fxn)),silent=ro.r('TRUE'))
			if 'try-error' in rclass(model_out):
				status = -4
			else:
				status = 1
				valid = True
		else:
			status = -5
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':')
			marker_info[x + '.beta'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate']))
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error'])
			marker_info[x + '.p'] = '%.2e' % (2 * scipy.norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error'])))
	else:
		for x in focus:
			marker_info[x + '.beta'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info.drop('marker_unique',axis=1,inplace=True)

def CalcCoxPH(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	if marker_info['filter'] == 0:
		rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))])
		model_out=rtry(rsummary(survival.coxph(ro.r(model),data=rmodel_df,control=ro.r('coxph.control(iter.max = 100)'))),silent=ro.r('TRUE'))
		if 'try-error' in rclass(model_out):
			status = -4
		else:
			status = 1
			valid = True
	else:
		status = -5
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		conf_int = py2r.convert_robj(model_out.rx('conf.int'))['conf.int']
		for x in focus:
			xt = x.replace('*',':')
			marker_info[x + '.n'] = '%d' % (np.array(model_out.rx('n')[0])[0])
			marker_info[x + '.beta'] = '%.5g' % (coef.loc[xt,'coef'])
			marker_info[x + '.or'] = '%.5g' % (coef.loc[xt,'exp(coef)'])
			marker_info[x + '.ci_lower'] = '%.5g' % (conf_int.loc[xt,'lower .95'])
			marker_info[x + '.ci_upper'] = '%.5g' % (conf_int.loc[xt,'upper .95'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'se(coef)'])
			marker_info[x + '.robust_stderr'] = '%.5g' % (coef.loc[xt,'robust se'])
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'z'])
			marker_info[x + '.p'] = '%.2e' % (coef.loc[xt,'Pr(>|z|)'])
	else:
		for x in focus:
			marker_info[x + '.n'] = float('NaN')
			marker_info[x + '.beta'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.ci_lower'] = float('NaN')
			marker_info[x + '.ci_upper'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.robust_stderr'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info.drop('marker_unique',axis=1,inplace=True)

def CalcIndepTests(marker_info, model_df):
	cols = marker_info['marker_unique'][marker_info['filter'] == 0]
	tot_tests = len(cols)
	if len(cols) > 1:
		markers = model_df[cols]
		markers_cor = markers.corr()
		markers_cor_eigvals = np.linalg.eigvalsh(markers_cor)
		markers_cor_eigvalsnew = [x if x > 0 else 0 for x in markers_cor_eigvals]
		n_indep = sum([1 if x else 0 for x in markers_cor_eigvals>1]+(markers_cor_eigvalsnew-np.floor(markers_cor_eigvalsnew)))
	elif len(cols) == 1:
		n_indep = 1
	else:
		n_indep = 0
	return ('%.5g' % n_indep, '%d' % tot_tests)
