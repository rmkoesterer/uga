from Model import *
rtry = ro.r('try')
rsummary = ro.r('summary')
rclass = ro.r('class')

def CalcGEE(marker_info, model_df, model_vars_dict, model, method, iid, fid, fxn, focus, dep_var, corstr):
	valid = False
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'geeglm(' + model.replace('marker',marker_info['marker_unique']) + ',id=' + fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + fxn + ',corstr="' + corstr + '")'
	try:
		ro.globalenv['model_out']=rtry(rsummary(ro.r(cmd)))
	except RRuntimeError:
		cmd = 'geeglm(' + model.replace('marker',marker_info['marker_unique']) + ',id=' + fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + fxn + ',corstr="independence")'
		try:
			ro.globalenv['model_out']=rtry(rsummary(ro.r(cmd)))
		except RRuntimeError:
			marker_info['status'] = '%d' % (-5)
		else:
			if ro.r('model_out$error') != 1:
				marker_info['status'] = '%d' % (2)
				valid = True
			else:
				marker_info['status'] = '%d' % (-4)
	else:
		if ro.r('model_out$error') != 1:
			marker_info['status'] = '%d' % (1)
			valid = True
		else:
			marker_info['status'] = '%d' % (-3)
	if valid:
		coef = py2r.convert_robj(ro.r('model_out$coefficients'))
		coef.index.values[coef.index.values == marker_info['marker_unique']] = 'marker'
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std.err'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
		marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
	return marker_info
"""
def CalcGEE(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var, geepack, corstr = 'exchangeable'):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if model.find('*') != -1:
		model_df['ugaInter'] = model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[0]]*model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[1]]			
	if (fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(model_df['marker'].unique()) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var])))) or (model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().stack().value_counts()[1] != model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().shape[0]):
		status = -2
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
			model_out=rtry(rsummary(geepack.geeglm(ro.r(model),data=rmodel_df,id=rmodel_df.rx2(fid),family=fxn,corstr=corstr)),silent=ro.r('TRUE'))
			if 'try-error' in rclass(model_out):
				model_out=rtry(rsummary(geepack.geeglm(ro.r(model),data=rmodel_df,id=rmodel_df.rx2(fid),family=fxn,corstr='independence')),silent=ro.r('TRUE'))
				if 'try-error' in rclass(model_out):
					status = -5
				else:
					status = 2 if model_out.rx2('error')[0] != 1 else -4
					valid = True if model_out.rx2('error')[0] != 1 else False
			else:
				status = 1 if model_out.rx2('error')[0] != 1 else -3
				valid = True if model_out.rx2('error')[0] != 1 else False
		else:
			status = -1
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std.err'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			else:
				marker_info[x + '.effect'] = 'NA'
				marker_info[x + '.stderr'] = 'NA'
				marker_info[x + '.or'] = 'NA'
				marker_info[x + '.z'] = 'NA'
				marker_info[x + '.p'] = 'NA'
	else:
		for x in focus:
			marker_info[x + '.effect'] = 'NA'
			marker_info[x + '.stderr'] = 'NA'
			marker_info[x + '.or'] = 'NA'
			marker_info[x + '.z'] = 'NA'
			marker_info[x + '.p'] = 'NA'
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info
"""
def CalcGLM(marker_info, model_df, model_vars_dict, model, iid, fid, fxn, focus, dep_var):
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'glm(' + model.replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),family="' + fxn + '")'
	model_out=ro.r(cmd)
	if model_out.rx2('converged')[0] and not model_out.rx2('boundary')[0] and not 'try-error' in rclass(rsummary(model_out)):
		marker_info['status'] = '%d' % (1)
		coef = py2r.convert_robj(rsummary(model_out).rx('coefficients'))['coefficients']
		coef.index.values[coef.index.values == marker_info['marker_unique']] = 'marker'
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
		marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
	else:
		marker_info['status'] = '%d' % (-3)
	return marker_info

"""
def CalcGLM(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var, rglm):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	valid = False
	if model.find('*') != -1:
		model_df['ugaInter'] = model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[0]]*model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[1]]			
	if (fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(model_df['marker'].unique()) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var])))) or (model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().stack().value_counts()[1] != model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().shape[0]):
		marker_info['status'] = -2
	else:
		rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(), strings_as_factors=False)
		rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
		for x in model_vars_dict.keys():
			if model_vars_dict[x]['class'] == 'factor':
				rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
				rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
				rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
		model_out=rglm(ro.r(model),data=rmodel_df,family=fxn)
		if model_out.rx2('converged')[0] and not model_out.rx2('boundary')[0] and not 'try-error' in rclass(rsummary(model_out)):
			valid = True
			marker_info['status'] = 1
		else:
			marker_info['status'] = -3
	if valid:
		coef = py2r.convert_robj(rsummary(model_out).rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info['n'] = len(py2r.convert_robj(ro.r.unique(rmodel_df.rx(iid))))
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info
"""

def CalcLME(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var, lme4):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if model.find('*') != -1:
		model_df['ugaInter'] = model_df[[x for x in re.split('\+|-',model) if x.find('*') != -1][0].split('*')[0]]*model_df[[x for x in re.split('\+|-',model) if x.find('*') != -1][0].split('*')[1]]			
	if (fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(model_df['marker'].unique()) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var])))) or (model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter']))]].corr().abs().stack().value_counts()[1] > model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter']))]].corr().abs().shape[0]):
		status = -2
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
				status = -3
			else:
				status = 1
				valid = True
		else:
			status = -1
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and fxn == 'binomial' else 'NA'
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else 'NA'
				marker_info[x + '.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else 'NA'
			else:
				marker_info[x + '.effect'] = 'NA'
				marker_info[x + '.stderr'] = 'NA'
				marker_info[x + '.or'] = 'NA'
				marker_info[x + '.z'] = 'NA'
				marker_info[x + '.p'] = 'NA'
	else:
		for x in focus:
			marker_info[x + '.effect'] = 'NA'
			marker_info[x + '.stderr'] = 'NA'
			marker_info[x + '.or'] = 'NA'
			marker_info[x + '.z'] = 'NA'
			marker_info[x + '.p'] = 'NA'
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcCoxPH(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var, survival):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if marker_info['filter'] == 0:
		rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,fid]))].dropna(), strings_as_factors=False)
		rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
		for x in model_vars_dict.keys():
			if model_vars_dict[x]['class'] == 'cluster':
				rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
				rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
				rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
		model_out=rtry(rsummary(survival.coxph(ro.r(model),data=rmodel_df,control=ro.r('coxph.control(iter.max = 100)'))),silent=ro.r('TRUE'))
		if 'try-error' in rclass(model_out):
			status = -2
		else:
			status = 1
			valid = True
	else:
		status = -1
	if valid:
		coef = py2r.convert_robj(model_out.rx('coefficients'))['coefficients']
		conf_int = py2r.convert_robj(model_out.rx('conf.int'))['conf.int']
		for x in focus:
			xt = x.replace('*',':') if x.replace('*',':') in coef.index.values else x.replace('*',':').split(':')[1] + ':' + x.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'coef'])
				marker_info[x + '.or'] = '%.5g' % (coef.loc[xt,'exp(coef)'])
				marker_info[x + '.ci_lower'] = '%.5g' % (conf_int.loc[xt,'lower .95'])
				marker_info[x + '.ci_upper'] = '%.5g' % (conf_int.loc[xt,'upper .95'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'se(coef)'])
				marker_info[x + '.robust_stderr'] = '%.5g' % (coef.loc[xt,'robust se'])
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'z'])
				marker_info[x + '.p'] = '%.2e' % (coef.loc[xt,'Pr(>|z|)'])
			else:
				marker_info[x + '.effect'] = 'NA'
				marker_info[x + '.or'] = 'NA'
				marker_info[x + '.ci_lower'] = 'NA'
				marker_info[x + '.ci_upper'] = 'NA'
				marker_info[x + '.stderr'] = 'NA'
				marker_info[x + '.robust_stderr'] = 'NA'
				marker_info[x + '.z'] = 'NA'
				marker_info[x + '.p'] = 'NA'
		marker_info['n'] = '%d' % (np.array(model_out.rx('n')[0])[0])
	else:
		for x in focus:
			marker_info[x + '.effect'] = 'NA'
			marker_info[x + '.or'] = 'NA'
			marker_info[x + '.ci_lower'] = 'NA'
			marker_info[x + '.ci_upper'] = 'NA'
			marker_info[x + '.stderr'] = 'NA'
			marker_info[x + '.robust_stderr'] = 'NA'
			marker_info[x + '.z'] = 'NA'
			marker_info[x + '.p'] = 'NA'
		marker_info['n'] = 0
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcGEEBoss(marker_info, model_df, model_vars_dict, model, iid, fid, method, fxn, focus, dep_var, boss, corstr = 'exchangeable', thresh = 1e-7):
	model_df.rename(columns={marker_info['marker_unique']: 'marker'}, inplace=True)
	fxn = ro.r.gaussian() if fxn == 'gaussian' else ro.r.binomial()
	interact = ro.r('NULL')
	notes = 'NA'
	status = 0
	n = 0
	valid = False
	if model.find('*') != -1:
		interact_vars = [x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].split('*')
		interact = [x for x in interact_vars if x != 'marker'][0]
		model_df['ugaInter'] = model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[0]]*model_df[[x for x in re.split('\+|-',model) if x.replace('factor(','').replace(')','').find('*') != -1][0].replace('factor(','').replace(')','').split('*')[1]]
	if (fxn == 'binomial' and (marker_info['freq.ctrl'] == 'NA' or marker_info['freq.ctrl'] < 0.001 or marker_info['freq.ctrl'] > 0.999 or marker_info['freq.case'] < 0.001 or marker_info['freq.case'] > 0.999 or (len(model_df['marker'].unique()) < 3 and 0 in pd.crosstab(model_df['marker'],model_df[dep_var])))) or (model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().stack().value_counts()[1] != model_df[[x for x in model_df if x in list(set(model_vars_dict.keys() + ['ugaInter'])) and (x == 'ugaInter' or model_vars_dict[x]['type'] != "dependent")]].corr().abs().shape[0]):
		status = -2
	else:
		if marker_info['filter'] == 0:
			rmodel_df = py2r.convert_to_r_dataframe(model_df[list(set(model_vars_dict.keys() + [iid,'id']))].dropna(), strings_as_factors=False)
			rmodel_df = ro.r.subset(rmodel_df,rmodel_df.rx('marker').ro != "NA")
			n = len(py2r.convert_robj(ro.r.unique(rmodel_df.rx(iid))))
			for x in model_vars_dict.keys():
				if model_vars_dict[x]['class'] == 'factor':
					rmodel_df.colnames=ro.StrVector([x + '_ugaFactored' if a == x else a for a in list(rmodel_df.colnames)])
					rmodel_df=ro.r.cbind(rmodel_df,ugaConvert=ro.r('factor')(rmodel_df.rx2(x + '_ugaFactored')))
					rmodel_df.colnames=ro.StrVector([x if a == 'ugaConvert' else a for a in list(rmodel_df.colnames)])
			for col in rmodel_df:
				col.rclass = None
			try:
				bs=rtry(boss.boss_set(ro.r.formula(model.replace('+marker','').replace('marker+','')),id=rmodel_df.rx2('id'),type="gee",E_name=interact,family=fxn,corstr=corstr,data=rmodel_df),silent=ro.r('TRUE'))
			except RRuntimeError:
				status = -4
			else:
				try:
					model_out=rtry(boss.boss_fit(rmodel_df.rx2('marker'),bs,thresh=thresh,sattdf=ro.r('TRUE')),silent=ro.r('TRUE'))
				except RRuntimeError:
					status = -5
				else:
					status = 1
					valid = True
		else:
			status = -1
	if valid:
		marker_info['marker.effect'] = '%.5g' % (py2r.convert_robj(model_out.rx2('beta.main'))[0]) if py2r.convert_robj(model_out.rx2('beta.main')[0]) else 'NA'
		marker_info['marker.v'] = '%.5g' % (py2r.convert_robj(model_out.rx2('v.main'))[0]) if py2r.convert_robj(model_out.rx2('v.main')[0]) else 'NA'
		if interact != ro.r('NULL'):
			marker_info['inter.effect'] = '%.5g' % (py2r.convert_robj(model_out.rx2('beta.inter'))[0]) if py2r.convert_robj(model_out.rx2('beta.inter')[0]) else 'NA'
			marker_info['inter.v'] = '%.5g' % (py2r.convert_robj(model_out.rx2('v.inter'))[0]) if py2r.convert_robj(model_out.rx2('v.inter')[0]) else 'NA'
			marker_info['inter.cov'] = '%.5g' % (py2r.convert_robj(model_out.rx2('cov.inter'))[0]) if py2r.convert_robj(model_out.rx2('cov.inter')[0]) else 'NA'
		marker_info['satt.df'] = '%.5g' % (py2r.convert_robj(model_out.rx2('df.satt'))[0]) if py2r.convert_robj(model_out.rx2('df.satt')[0]) else 'NA'
		marker_info['marker.stderr'] = '%.5g' % (math.sqrt(py2r.convert_robj(model_out.rx2('v.main'))[0])) if py2r.convert_robj(model_out.rx2('v.main')[0]) else 'NA'
		marker_info['marker.or'] = '%.5g' % (math.exp(py2r.convert_robj(model_out.rx2('beta.main'))[0])) if py2r.convert_robj(model_out.rx2('beta.main')[0]) and not py2r.convert_robj(model_out.rx2('beta.main')[0]) > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.main')[0]) < -709.782712893384 and fxn == 'binomial' else 'NA'
		marker_info['marker.z'] = '%.5g' % (py2r.convert_robj(model_out.rx2('beta.main'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.main'))[0])) if py2r.convert_robj(model_out.rx2('beta.main'))[0] and py2r.convert_robj(model_out.rx2('v.main'))[0] and not py2r.convert_robj(model_out.rx2('beta.main'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.main'))[0] < -709.782712893384 else 'NA'
		marker_info['marker.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(py2r.convert_robj(model_out.rx2('beta.main'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.main'))[0])))) if py2r.convert_robj(model_out.rx2('beta.main'))[0] and py2r.convert_robj(model_out.rx2('v.main'))[0] and not py2r.convert_robj(model_out.rx2('beta.main'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.main'))[0] < -709.782712893384 else 'NA'
		marker_info['marker.sattdf.p'] = '%.2e' % (2 * t.sf(abs(py2r.convert_robj(model_out.rx2('beta.main'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.main'))[0])),py2r.convert_robj(model_out.rx2('df.satt'))[0])) if py2r.convert_robj(model_out.rx2('df.satt')[0]) and py2r.convert_robj(model_out.rx2('beta.main'))[0] and py2r.convert_robj(model_out.rx2('v.main'))[0] and not py2r.convert_robj(model_out.rx2('beta.main'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.main'))[0] < -709.782712893384 else 'NA'
		if interact != ro.r('NULL'):
			marker_info['inter.stderr'] = '%.5g' % (math.sqrt(py2r.convert_robj(model_out.rx2('v.inter'))[0])) if py2r.convert_robj(model_out.rx2('v.inter')[0]) else 'NA'
			marker_info['inter.or'] = '%.5g' % (math.exp(py2r.convert_robj(model_out.rx2('beta.inter'))[0])) if py2r.convert_robj(model_out.rx2('beta.inter')[0]) and not py2r.convert_robj(model_out.rx2('beta.inter')[0]) > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.inter')[0]) < -709.782712893384 and fxn == 'binomial' else 'NA'
			marker_info['inter.z'] = '%.5g' % (py2r.convert_robj(model_out.rx2('beta.inter'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.inter'))[0])) if py2r.convert_robj(model_out.rx2('beta.inter'))[0] and py2r.convert_robj(model_out.rx2('v.inter'))[0] and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] < -709.782712893384 else 'NA'
			marker_info['inter.p'] = '%.2e' % (2 * norm.cdf(-1 * abs(py2r.convert_robj(model_out.rx2('beta.inter'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.inter'))[0])))) if py2r.convert_robj(model_out.rx2('beta.inter'))[0] and py2r.convert_robj(model_out.rx2('v.inter'))[0] and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] < -709.782712893384 else 'NA'
			marker_info['inter.sattdf.p'] = '%.2e' % (2 * t.sf(abs(py2r.convert_robj(model_out.rx2('beta.inter'))[0] / math.sqrt(py2r.convert_robj(model_out.rx2('v.inter'))[0])),py2r.convert_robj(model_out.rx2('df.satt'))[0])) if py2r.convert_robj(model_out.rx2('df.satt')[0]) and py2r.convert_robj(model_out.rx2('beta.inter'))[0] and py2r.convert_robj(model_out.rx2('v.inter'))[0] and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] > 709.782712893384 and not py2r.convert_robj(model_out.rx2('beta.inter'))[0] < -709.782712893384 else 'NA'
	else:
		marker_info['marker.effect'] = 'NA'
		marker_info['marker.v'] = 'NA'
		if interact != ro.r('NULL'):
			marker_info['inter.effect'] = 'NA'
			marker_info['inter.v'] = 'NA'
			marker_info['inter.cov'] = 'NA'
		marker_info['satt.df'] = 'NA'
		marker_info['marker.stderr'] = 'NA'
		marker_info['marker.or'] = 'NA'
		marker_info['marker.z'] = 'NA'
		marker_info['marker.p'] = 'NA'
		marker_info['marker.sattdf.p'] = 'NA'
		if interact != ro.r('NULL'):
			marker_info['inter.stderr'] = 'NA'
			marker_info['inter.or'] = 'NA'
			marker_info['inter.z'] = 'NA'
			marker_info['inter.p'] = 'NA'
			marker_info['inter.sattdf.p'] = 'NA'
	marker_info['n'] = n
	marker_info['status'] = status
	model_df.rename(columns={'marker': marker_info['marker_unique']}, inplace=True)
	return marker_info

def CalcEffTests(model_df):
	if model_df.shape[1] > 1:
		markers_cor = model_df.corr()
		markers_cor_eigvals = np.linalg.eigvalsh(markers_cor)
		markers_cor_eigvalsnew = [x if x > 0 else 0 for x in markers_cor_eigvals]
		n_eff = sum([1 if x else 0 for x in markers_cor_eigvals>1]+(markers_cor_eigvalsnew-np.floor(markers_cor_eigvalsnew)))
	else:
		n_eff = 1
	return '%.5g' % n_eff

def PrepScores(snp_info, z, model, pheno, seqmeta, family=None):
	family = ro.r('gaussian()') if family is None or family == 'gaussian' else ro.r('binomial()')
	rsnp_info = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	rz = ro.r('as.matrix')(py2r.convert_to_r_dataframe(z, strings_as_factors=False))
	rpheno = py2r.convert_to_r_dataframe(pheno, strings_as_factors=False)
	ro.globalenv['ps'] = seqmeta.prepScores(Z = rz, formula = ro.r(model), SNPInfo = rsnp_info, data = rpheno, family=family)
	return ro.globalenv['ps']

def PrepScoresFam(snp_info, z, model, pheno, seqmeta, kinship):
	rsnp_info = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	rz = ro.r('as.matrix')(py2r.convert_to_r_dataframe(z, strings_as_factors=False))
	rpheno = py2r.convert_to_r_dataframe(pheno, strings_as_factors=False)
	ro.globalenv['ps'] = seqmeta.prepScores(Z = rz, formula = ro.r(model), SNPInfo = rsnp_info, data = rpheno, kins = kinship, sparse=ro.r('FALSE'))
	return ro.globalenv['ps']

def SkatOMeta(cmd, snp_info, seqmeta, rho = 1):
	ro.globalenv['rsnp_info'] = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	ro.globalenv['r_rho'] = ro.r('seq(0,1,' + str(rho) + ')')
	try:
		result = rtry(ro.reval(cmd),silent=ro.r('TRUE'))
		f1=open('test.result', 'w+')
		f1.write(str(result))
		f1.close()
	except RRuntimeError:
		result_df = pd.DataFrame({'p': ['NA'], 'pmin': ['NA'], 'rho': ['NA'], 'cmaf': ['NA'], 'nmiss': ['NA'], 'nsnps': ['NA'], 'errflag': ['100']})
		result_df.index = [1]
	else:
		result_df = py2r.convert_robj(result)
		result_df['p'] = '%.2e' % (result_df['p']) if result_df['p'][1] != 0 and result_df['errflag'][1] == 0 else 'NA'
		result_df['pmin'] = '%.2e' % (result_df['pmin']) if result_df['pmin'][1] != 0 and result_df['errflag'][1] == 0 else 'NA'
		result_df['rho'] = '%.5g' % (result_df['rho'])
		result_df['cmaf'] = '%.5g' % (result_df['cmaf'])
		result_df['nmiss'] = '%d' % (result_df['nmiss'])
		result_df['nsnps'] = '%d' % (result_df['nsnps'])
		result_df['errflag'] = '%d' % (result_df['errflag'])
	return result_df

def SkatMeta(cmd, snp_info, seqmeta):
	ro.globalenv['rsnp_info'] = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	try:
		result = rtry(ro.reval(cmd),silent=ro.r('TRUE'))
	except RRuntimeError:
		result_df = pd.DataFrame({'p': ['NA'], 'Qmeta': ['NA'], 'cmaf': ['NA'], 'nmiss': ['NA'], 'nsnps': ['NA']})
		result_df.index = [1]
	else:
		result_df = py2r.convert_robj(result)
		result_df['p'] = '%.2e' % (result_df['p']) if result_df['p'][1] != 0 else 'NA'
		result_df['Qmeta'] = '%.5g' % (result_df['Qmeta'])
		result_df['cmaf'] = '%.5g' % (result_df['cmaf'])
		result_df['nmiss'] = '%d' % (result_df['nmiss'])
		result_df['nsnps'] = '%d' % (result_df['nsnps'])
	return result_df

def BurdenMeta(cmd, snp_info, seqmeta):
	ro.globalenv['rsnp_info'] = py2r.convert_to_r_dataframe(snp_info, strings_as_factors=False)
	try:
		result = rtry(ro.reval(cmd),silent=ro.r('TRUE'))
	except RRuntimeError:
		result_df = pd.DataFrame({'p': ['NA'], 'beta': ['NA'], 'se': ['NA'], 'cmafTotal': ['NA'], 'cmafUsed': ['NA'], 'nsnpsTotal': ['NA'], 'nsnpsUsed': ['NA'], 'nmiss': ['NA']})
		result_df.index = [1]
	else:
		result_df = py2r.convert_robj(result)
		result_df['p'] = '%.2e' % (result_df['p']) if result_df['p'][1] != 0 else 'NA'
		result_df['beta'] = '%.5g' % (result_df['beta'])
		result_df['se'] = '%.5g' % (result_df['se'])
		result_df['cmafTotal'] = '%.5g' % (result_df['cmafTotal'])
		result_df['cmafUsed'] = '%.5g' % (result_df['cmafUsed'])
		result_df['nsnpsTotal'] = '%d' % (result_df['nsnpsTotal'])
		result_df['nsnpsUsed'] = '%d' % (result_df['nsnpsUsed'])
		result_df['nmiss'] = '%d' % (result_df['nmiss'])
	return result_df


def SkatOMetaEmpty(chr, start, end, reg_id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],'incl': [meta_incl_string],
								'p': ['NA'],'pmin': ['NA'],'rho': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA'],'errflag': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],
								'p': ['NA'],'pmin': ['NA'],'rho': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA'],'errflag': ['NA']})
	return results[['chr','start','end','reg_id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

def SkatMetaEmpty(chr, start, end, reg_id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],'incl': [meta_incl_string],
								'p': ['NA'],'Qmeta': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],
								'p': ['NA'],'Qmeta': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA']})
	return results[['chr','start','end','reg_id','p','Qmeta','cmaf','nmiss','nsnps']]

def BurdenMetaEmpty(chr, start, end, reg_id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],'incl': [meta_incl_string],
								'p': ['NA'],'beta': ['NA'],'se': ['NA'],'cmafTotal': ['NA'],'cmafUsed': ['NA'],
								'nsnpsTotal': ['NA'],'nsnpsUsed': ['NA'],'nmiss': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'reg_id': [reg_id],
								'p': ['NA'],'beta': ['NA'],'se': ['NA'],'cmafTotal': ['NA'],'cmafUsed': ['NA'],
								'nsnpsTotal': ['NA'],'nsnpsUsed': ['NA'],'nmiss': ['NA']})
	return results[['chr','start','end','reg_id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
