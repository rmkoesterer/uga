import numpy as np
import pandas as pd
import re
from Messages import Error

def ExtractModelVars(pheno,model,fid,iid,fxn=None,sex=None,pheno_sep='\t'):
	model_vars_dict = {}
	dependent = re.split('\\~',model)[0]
	independent = re.split('\\~',model)[1]
	vars_df = pd.read_table(pheno,sep=pheno_sep,dtype='str')
	for x in [a for a in list(set(re.split('Surv\(|,|\)|~|\+|cluster\(|\(1\||\*|factor\(',model))) if a != '']:
		mtype = ''
		if dependent.find(x) != -1:
			pre = dependent.split(x)[0] if dependent.split(x)[0] != '' else '.'
			post = dependent.split(x)[1] if dependent.split(x)[1] != '' else '.'
			if not pre[-1].isalpha() and not post[0].isalpha():
				mtype = 'dependent'
		if independent.find(x) != -1:
			pre = independent.split(x)[0] if independent.split(x)[0] != '' else '.'
			post = independent.split(x)[1] if independent.split(x)[1] != '' else '.'
			if not pre[-1].isalpha() and not post[0].isalpha():
				if mtype == '':
					mtype = 'independent'
				else:
					print Error("a column in the phenotype file is defined in the model as both an independent and dependent variable")
					return
		if mtype != '':
			if model[model.find(x)-7:model.find(x)] == 'factor(':
				model_vars_dict[x] = {'class': 'factor', 'type': mtype}
			elif model[model.find(x)-15:model.find(x)] == 'ordered(factor(':
				model_vars_dict[x] = {'class': 'orderedfactor', 'type': mtype}
			elif model[model.find(x)-3:model.find(x)] == '(1|':
				model_vars_dict[x] = {'class': 'random', 'type': mtype}
			elif model[model.find(x)-8:model.find(x)] == 'cluster(':
				model_vars_dict[x] = {'class': 'cluster', 'type': mtype}
			else:
				model_vars_dict[x] = {'class': 'numeric', 'type': mtype}

	for x in model_vars_dict.keys():
		if x in list(vars_df.columns):
			print "   %s variable %s found" % (model_vars_dict[x]['type'], x)
		elif x in ['marker','marker1','marker2','marker.interact']:
			print "   %s variable %s skipped" % (model_vars_dict[x]['type'], x)
		else:
			print "   %s variable %s not found" % (model_vars_dict[x]['type'], x)
			print Error("model variable missing from phenotype file")
			return
	
	if sex:
		if sex in list(vars_df.columns):
			print "   sex column %s found" % sex
		else:
			print "   sex column %s not found" % sex
			print Error("sex column missing from phenotype file")
			return

	vars_df = vars_df[list(set([a for a in [fid,iid,sex] if a] + list([a for a in model_vars_dict.keys() if a != 'marker'])))]
	
	vars_df[fid] = vars_df[fid].astype(str)
	vars_df[iid] = vars_df[iid].astype(str)
	
	##### EXTRACT CASE/CTRL IF BINOMIAL fxn #####
	for x in model_vars_dict.keys():
		if model_vars_dict[x]['type'] == 'dependent' and fxn == 'binomial':
			vars_df = vars_df[vars_df[x].isin([str(case),str(ctrl)])]
			vars_df[x] = vars_df[x].map({str(ctrl): '0', str(case): '1'})
	vars_df.dropna(inplace = True)
	if len(vars_df.index) == 0:
		print Error("no data left for analysis")
		return
	return vars_df, model_vars_dict

def GetDelimiter(delimiter):
	if delimiter == 'tab':
		delimiter = '\t'
	elif delimiter == 'space':
		delimiter = ' '
	else:
		delimiter = ','
	return delimiter

def GetFocus(method,model,vars_df):
	focus = ['(Intercept)'] if method.split('_')[0] in ['gee','glm','lme'] else []
	for x in re.sub("\+|\~|\-",",",model.split('~')[1]).split(','):
		if not x[0] == '(' and x.find('cluster(') == -1:
			if x.find('factor(') != -1:
				for v in vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64) - min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64)):
					if v != min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64) - min(vars_df[re.sub('factor\(|\)','',x)].unique().astype(np.int64))):
						focus.append(x + str(v))
			elif x.find('*') != -1:
				focus.append('*'.join(sorted(x.split('*'))))
			else:
				focus.append(x)
	return focus
