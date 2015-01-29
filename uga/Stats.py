#import pickle
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
#import statsmodels.api as sm
import pandas.rpy.common as py2r
import math
from scipy.stats import norm

geepack = importr('geepack')
lme4 = importr('lme4')
survival = importr('survival')
seqmeta = importr('seqMeta')
kinship2 = importr('kinship2')

rtry = ro.r('try')
rsummary = ro.r('summary')
rclass = ro.r('class')
rglm = ro.r('glm')
base=py2r.importr('base')

def GenerateFilterCode(marker_info, miss = None, freq = None, rsq = None, hwe = None):
	filter = 0
	if (not miss is None and not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) < miss) or (not math.isnan(marker_info['callrate']) and float(marker_info['callrate']) == 0) or (math.isnan(marker_info['callrate'])):
		filter += 1000
	if (not freq is None and not math.isnan(marker_info['freq']) and ((float(marker_info['freq']) < freq or float(marker_info['freq']) > 1-freq or float(marker_info['freq']) == 0 or float(marker_info['freq']) == 1) and (float(marker_info['freq.unrel']) < freq or float(marker_info['freq.unrel']) > 1-freq or float(marker_info['freq.unrel']) == 0 or float(marker_info['freq.unrel']) == 1))) or (math.isnan(marker_info['freq'])):
		filter += 100
	if (not rsq is None and not math.isnan(marker_info['rsq']) and (float(marker_info['rsq']) < rsq and float(marker_info['rsq.unrel']) < rsq)) or (math.isnan(marker_info['rsq'])):
		filter += 10
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
	n = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).codes) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(), strings_as_factors=False)
			rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
			n = len(py2r.convert_robj(ro.r.unique(rmodel_df.rx(iid))))
			for x in model_vars_dict.keys():
				if model_vars_dict[x]['class'] == 'factor':
					rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
					rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
					rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
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
			marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std.err'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			marker_info[x + '.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
	else:
		for x in focus:
			marker_info[x + '.effect'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcGLM(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).codes) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(), strings_as_factors=False)
			rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
			n = len(py2r.convert_robj(ro.r.unique(rmodel_df.rx(iid))))
			for x in model_vars_dict.keys():
				if model_vars_dict[x]['class'] == 'factor':
					rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
					rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
					rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
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
			marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			marker_info[x + '.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
	else:
		for x in focus:
			marker_info[x + '.effect'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcLME(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(pd.Categorical.from_array(model_df['marker']).codes) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var]))):
		status = -3
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(),strings_as_factors=False)
			rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
			n = len(py2r.convert_robj(ro.r.unique(rmodel_df.rx(iid))))
			for x in model_vars_dict.keys():
				if model_vars_dict[x]['class'] in ['factor','random']:
					rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
					rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
					rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
			model_out=rtry(rsummary(lme4.glmer(ro.r(model),data=rmodel_df,REML=ro.r('FALSE'),family=fxn)),silent=ro.r('TRUE'))
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
			marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
			marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			marker_info[x + '.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
	else:
		for x in focus:
			marker_info[x + '.effect'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcCoxPH(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	valid = False
	if marker_info['filter'] == 0:
		rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(), strings_as_factors=False)
		for x in model_vars_dict.keys():
			if model_vars_dict[x]['class'] == 'cluster':
				rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
				rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
				rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
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
			marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'coef'])
			marker_info[x + '.or'] = '%.5g' % (coef.loc[xt,'exp(coef)'])
			marker_info[x + '.ci_lower'] = '%.5g' % (conf_int.loc[xt,'lower .95'])
			marker_info[x + '.ci_upper'] = '%.5g' % (conf_int.loc[xt,'upper .95'])
			marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'se(coef)'])
			marker_info[x + '.robust_stderr'] = '%.5g' % (coef.loc[xt,'robust se'])
			marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'z'])
			marker_info[x + '.p'] = '%.2e' % (coef.loc[xt,'Pr(>|z|)'])
		marker_info['n'] = '%d' % (np.array(model_out.rx('n')[0])[0])
	else:
		for x in focus:
			marker_info[x + '.effect'] = float('NaN')
			marker_info[x + '.or'] = float('NaN')
			marker_info[x + '.ci_lower'] = float('NaN')
			marker_info[x + '.ci_upper'] = float('NaN')
			marker_info[x + '.stderr'] = float('NaN')
			marker_info[x + '.robust_stderr'] = float('NaN')
			marker_info[x + '.z'] = float('NaN')
			marker_info[x + '.p'] = float('NaN')
		marker_info['n'] = 0
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcEffTests(model_df, mem):
	if model_df.shape[1] > 1:
		markers_cor = model_df.corr()
		markers_cor_eigvals = np.linalg.eigvalsh(markers_cor)
		markers_cor_eigvalsnew = [x if x > 0 else 0 for x in markers_cor_eigvals]
		n_eff = sum([1 if x else 0 for x in markers_cor_eigvals>1]+(markers_cor_eigvalsnew-np.floor(markers_cor_eigvalsnew)))
	else:
		n_eff = 1
	return '%.5g' % n_eff

def CalcFamSkatO(snp_info, z, model, pheno, pedigree):
	rsnp_info = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	rz = ro.r('as.matrix')(py2r.convert_to_r_dataframe(z, strings_as_factors=False))
	rpheno = py2r.convert_to_r_dataframe(pheno, strings_as_factors=False)
	rpedigree = py2r.convert_to_r_dataframe(pedigree, strings_as_factors=False)
	kins = kinship2.makekinship(rpedigree.rx2('FID'),rpedigree.rx2('IID'),rpedigree.rx2('PAT'),rpedigree.rx2('MAT'))
	ro.globalenv['ps'] = seqmeta.prepScores(Z = rz, formula = ro.r(model), SNPInfo = rsnp_info, data = rpheno, kins = kins, sparse=ro.r('FALSE'))
	result = rtry(seqmeta.skatOMeta(base.as_symbol('ps'), SNPInfo = rsnp_info),silent=ro.r('TRUE'))
	result_df = py2r.convert_robj(result)
	result_df['p'] = '%.2e' % (result_df['p'])
	result_df['pmin'] = '%.2e' % (result_df['pmin'])
	result_df['rho'] = '%.5g' % (result_df['rho'])
	result_df['cmaf'] = '%.5g' % (result_df['cmaf'])
	result_df['nmiss'] = '%d' % (result_df['nmiss'])
	result_df['nsnps'] = '%d' % (result_df['nsnps'])
	result_df['errflag'] = '%d' % (result_df['errflag'])
	return result_df

def CalcFamSkat(snp_info, z, model, pheno, pedigree):
	rsnp_info = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	rz = ro.r('as.matrix')(py2r.convert_to_r_dataframe(z, strings_as_factors=False))
	rpheno = py2r.convert_to_r_dataframe(pheno, strings_as_factors=False)
	rpedigree = py2r.convert_to_r_dataframe(pedigree, strings_as_factors=False)
	kins = kinship2.makekinship(rpedigree.rx2('FID'),rpedigree.rx2('IID'),rpedigree.rx2('PAT'),rpedigree.rx2('MAT'))
	ro.globalenv['ps'] = seqmeta.prepScores(Z = rz, formula = ro.r(model), SNPInfo = rsnp_info, data = rpheno, kins = kins, sparse=ro.r('FALSE'))
	result = rtry(seqmeta.skatMeta(base.as_symbol('ps'), SNPInfo = rsnp_info),silent=ro.r('TRUE'))
	result_df = py2r.convert_robj(result)
	result_df['p'] = '%.2e' % (result_df['p'])
	result_df['Qmeta'] = '%.5g' % (result_df['Qmeta'])
	result_df['cmaf'] = '%.5g' % (result_df['cmaf'])
	result_df['nmiss'] = '%d' % (result_df['nmiss'])
	result_df['nsnps'] = '%d' % (result_df['nsnps'])
	result_df['errflag'] = '%d' % (result_df['errflag'])
	return result_df
