from Model import *
rtry = ro.r('try')
rsummary = ro.r('summary')
rclass = ro.r('class')

def CalcGEE(marker_info, model_df, model_vars_dict, model, model_fxn, iid, fid, family, focus, dep_var, corstr):
	valid = False
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'geeglm(' + model.replace('marker',marker_info['marker_unique']) + ',id=' + fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + family + ',corstr="' + corstr + '")'
	try:
		ro.globalenv['model_out']=rtry(rsummary(ro.r(cmd)))
	except RRuntimeError:
		print "      " + marker_info['marker_unique'] + ": RRuntimeError for function geeglm with corstr=" + corstr
		cmd = 'geeglm(' + model.replace('marker',marker_info['marker_unique']) + ',id=' + fid + ',data=na.omit(model_df[,names(model_df) %in% cols]),family=' + family + ',corstr="independence")'
		try:
			ro.globalenv['model_out']=rtry(rsummary(ro.r(cmd)))
		except RRuntimeError:
			print "      " + marker_info['marker_unique'] + ": RRuntimeError for function geeglm with corstr=independence"
			marker_info['status'] = '%d' % (-4)
		else:
			if ro.r('model_out$error') != 1:
				marker_info['status'] = '%d' % (2)
				valid = True
			else:
				marker_info['status'] = '%d' % (-3)
	else:
		if ro.r('model_out$error') != 1:
			marker_info['status'] = '%d' % (1)
			valid = True
		else:
			marker_info['status'] = '%d' % (-2)
	if valid:
		coef = py2r.convert_robj(ro.r('model_out$coefficients'))
		for x in focus:
			xt = x.replace('marker',marker_info['marker_unique'])
			xt = xt.replace('*',':') if xt.replace('*',':') in coef.index.values else xt.replace('*',':').split(':')[1] + ':' + xt.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std.err'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and family == 'binomial' else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std.err']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
		marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
	return marker_info

def CalcGLM(marker_info, model_df, model_vars_dict, model, iid, fid, family, focus, dep_var):
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'glm(' + model.replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),family="' + family + '")'
	try:
		model_out=ro.r(cmd)
	except RRuntimeError:
		print "      " + marker_info['marker_unique'] + ": RRuntimeError for function glm"
		marker_info['status'] = '%d' % (-3)
	else:
		if model_out.rx2('converged')[0] and not model_out.rx2('boundary')[0]:
			marker_info['status'] = '%d' % (1)
			coef = py2r.convert_robj(rsummary(model_out).rx('coefficients'))['coefficients']
			coef.index.values[coef.index.values == marker_info['marker_unique']] = 'marker'
			for x in focus:
				xt = x.replace('marker',marker_info['marker_unique'])
				xt = xt.replace('*',':') if xt.replace('*',':') in coef.index.values else xt.replace('*',':').split(':')[1] + ':' + xt.replace('*',':').split(':')[0]
				if xt in coef.index.values:
					marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
					marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
					marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 and family == 'binomial' else float('nan')
					marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
					marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
		else:
			marker_info['status'] = '%d' % (-2)
	return marker_info

def CalcLMEBinomial(marker_info,model_df,model_vars_dict,model,model_fxn,iid,fid,focus,dep_var,lmer_ctrl,lrt):
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'glmer(' + model.replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),family="binomial",control=glmerControl(' + lmer_ctrl + '))'
	try:
		ro.globalenv['model_out']=ro.r(cmd)
	except RRuntimeError:
		print '      ' + marker_info['marker_unique'] + ': RRuntimeError for function glmer with family="binomial"'
		marker_info['status'] = '%d' % (-2)
	else:
		marker_info['status'] = '%d' % (1)
		coef = py2r.convert_robj(ro.r('summary(model_out)$coefficients'))
		coef.index.values[coef.index.values == marker_info['marker_unique']] = 'marker'
		for x in focus:
			xt = x.replace('marker',marker_info['marker_unique'])
			xt = xt.replace('*',':') if xt.replace('*',':') in coef.index.values else xt.replace('*',':').split(':')[1] + ':' + xt.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'z value']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (coef.loc[xt,'Pr(>|z|)']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
			if lrt and x != '(Intercept)':
				xx = x.split(')')[0] + ')' if x.find('factor(') != -1 else x
				cmd2 = 'glmer(' + model.replace(xx + '+','').replace('+' + xx,'').replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),family="binomial",control=glmerControl(' + lmer_ctrl + '))'
				try:
					ro.globalenv['model_out2']=ro.r(cmd2)
				except RRuntimeError:
					print '      ' + marker_info['marker_unique'] + ': RRuntimeError for function glmer with family="binomial" and reduced model (model - ' + xx + ')'
					marker_info['status'] = '%d' % (-3)
				else:
					anova_out = py2r.convert_robj(ro.r('anova(model_out,model_out2)'))
					marker_info[x + '.anova.chisq'] = '%.5g' % (anova_out['Chisq'][1])
					marker_info[x + '.anova.p'] = '%.4e' % (anova_out['Pr(>Chisq)'][1])
		marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
	return marker_info

def CalcLMEGaussian(marker_info,model_df,model_vars_dict,model,model_fxn,iid,fid,focus,dep_var,lmer_ctrl,reml,lrt):
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'lmer(' + model.replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),REML=' + str(reml).upper() + ',control=lmerControl(' + lmer_ctrl + '))'
	try:
		ro.globalenv['model_out']=ro.r(cmd)
	except RRuntimeError:
		print '      ' + marker_info['marker_unique'] + ': RRuntimeError for function lmer with family="gaussian"'
		marker_info['status'] = '%d' % (-2)
	else:
		marker_info['status'] = '%d' % (1)
		coef = py2r.convert_robj(ro.r('summary(model_out)$coefficients'))
		coef.index.values[coef.index.values == marker_info['marker_unique']] = 'marker'
		for x in focus:
			xt = x.replace('marker',marker_info['marker_unique'])
			xt = xt.replace('*',':') if xt.replace('*',':') in coef.index.values else xt.replace('*',':').split(':')[1] + ':' + xt.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'Estimate'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'Std. Error'])
				marker_info[x + '.or'] = '%.5g' % (math.exp(coef.loc[xt,'Estimate'])) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.p'] = '%.4e' % (2 * norm.cdf(-1 * abs(coef.loc[xt,'Estimate'] / coef.loc[xt,'Std. Error']))) if not coef.loc[xt,'Estimate'] > 709.782712893384 and not coef.loc[xt,'Estimate'] < -709.782712893384 else float('nan')
				marker_info[x + '.satt.df'] = '%.5g' % (coef.loc[xt,'df'])
				marker_info[x + '.satt.t'] = '%.5g' % (coef.loc[xt,'t value'])
				marker_info[x + '.satt.p'] = '%.4e' % (coef.loc[xt,'Pr(>|t|)'])
				marker_info[x + '.kenrog.p'] = '%.4e' % (t.sf(coef.loc[xt,'t value'], ro.r('get_ddf_Lb(model_out, fixef(model_out))')))
			if lrt and x != '(Intercept)':
				xx = x.split(')')[0] + ')' if x.find('factor(') != -1 else x
				cmd2 = 'lmer(' + model.replace(xx + '+','').replace('+' + xx,'').replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),REML=' + str(reml).upper() + ',control=lmerControl(' + lmer_ctrl + '))'
				try:
					ro.globalenv['model_out2']=ro.r(cmd2)
				except RRuntimeError:
					print '      ' + marker_info['marker_unique'] + ': RRuntimeError for function glmer with family="binomial" and reduced model (model - ' + xx + ')'
					marker_info['status'] = '%d' % (-3)
				else:
					anova_out = py2r.convert_robj(ro.r('anova(model_out,model_out2)'))
					marker_info[x + '.anova.chisq'] = '%.5g' % (anova_out['Chisq'][1])
					marker_info[x + '.anova.p'] = '%.4e' % (anova_out['Pr(>Chisq)'][1])
		marker_info['n'] = '%d' % (len(model_df[list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))].dropna()[iid].unique()))
	return marker_info

def CalcCoxPH(marker_info, model_df, model_vars_dict, model, iid, fid, model_fxn, focus, dep_var, cph_ctrl):
	ro.globalenv['cols'] = list(set([a for a in model_vars_dict.keys() if a != 'marker'] + [iid,fid] + [marker_info['marker_unique']]))
	cmd = 'coxph(' + model.replace('marker',marker_info['marker_unique']) + ',data=na.omit(model_df[,names(model_df) %in% cols]),control=coxph.control(' + cph_ctrl + '))'
	try:
		model_out=ro.r(cmd)
	except RRuntimeError:
		print "      " + marker_info['marker_unique'] + ": RRuntimeError for function coxph"
		marker_info['status'] = '%d' % (-3)
	else:
		coef = py2r.convert_robj(rsummary(model_out).rx('coefficients'))['coefficients']
		conf_int = py2r.convert_robj(rsummary(model_out).rx('conf.int'))['conf.int']
		for x in focus:
			xt = x.replace('marker',marker_info['marker_unique'])
			xt = xt.replace('*',':') if xt.replace('*',':') in coef.index.values else xt.replace('*',':').split(':')[1] + ':' + xt.replace('*',':').split(':')[0]
			if xt in coef.index.values:
				marker_info[x + '.effect'] = '%.5g' % (coef.loc[xt,'coef'])
				marker_info[x + '.or'] = '%.5g' % (coef.loc[xt,'exp(coef)'])
				marker_info[x + '.ci_lower'] = '%.5g' % (conf_int.loc[xt,'lower .95'])
				marker_info[x + '.ci_upper'] = '%.5g' % (conf_int.loc[xt,'upper .95'])
				marker_info[x + '.stderr'] = '%.5g' % (coef.loc[xt,'se(coef)'])
				marker_info[x + '.robust_stderr'] = '%.5g' % (coef.loc[xt,'robust se'])
				marker_info[x + '.z'] = '%.5g' % (coef.loc[xt,'z'])
				marker_info[x + '.p'] = '%.2e' % (coef.loc[xt,'Pr(>|z|)'])
		marker_info['n'] = '%d' % (np.array(rsummary(model_out).rx('n')[0])[0])
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


def SkatOMetaEmpty(chr, start, end, id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],'incl': [meta_incl_string],
								'p': ['NA'],'pmin': ['NA'],'rho': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA'],'errflag': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],
								'p': ['NA'],'pmin': ['NA'],'rho': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA'],'errflag': ['NA']})
	return results[['chr','start','end','id','p','pmin','rho','cmaf','nmiss','nsnps','errflag']]

def SkatMetaEmpty(chr, start, end, id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],'incl': [meta_incl_string],
								'p': ['NA'],'Qmeta': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],
								'p': ['NA'],'Qmeta': ['NA'],'cmaf': ['NA'],'nmiss': ['NA'],
								'nsnps': ['NA']})
	return results[['chr','start','end','id','p','Qmeta','cmaf','nmiss','nsnps']]

def BurdenMetaEmpty(chr, start, end, id, meta_incl_string = None):
	if meta_incl_string is not None:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],'incl': [meta_incl_string],
								'p': ['NA'],'beta': ['NA'],'se': ['NA'],'cmafTotal': ['NA'],'cmafUsed': ['NA'],
								'nsnpsTotal': ['NA'],'nsnpsUsed': ['NA'],'nmiss': ['NA']})
	else:
		results = pd.DataFrame({'chr': [chr],'start': [start],'end': [end],'id': [id],
								'p': ['NA'],'beta': ['NA'],'se': ['NA'],'cmafTotal': ['NA'],'cmafUsed': ['NA'],
								'nsnpsTotal': ['NA'],'nsnpsUsed': ['NA'],'nmiss': ['NA']})
	return results[['chr','start','end','id','p','beta','se','cmafTotal','cmafUsed','nsnpsTotal','nsnpsUsed','nmiss']]
