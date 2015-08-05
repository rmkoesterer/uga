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

from __main__ import *
import scipy.stats as scipy
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

def Explore(cfg):
	Parse.PrintExploreOptions(cfg)

	if cfg['qq'] or cfg['mht']:
		import rpy2.robjects.lib.ggplot2 as ggplot2
		rgrid = importr('grid')
		grdevices = importr('grDevices')

	fcols = ['#chr','pos','a1','a2','marker','callrate','freq']
	fcols.append('freq.unrel')
	fcols.append('rsq')
	fcols.append('rsq.unrel')
	fcols.append('hwe')
	fcols.append('hwe.unrel')
	fcols.append('samples')
	fcols.append('n')
	fcols.append(cfg['stat'] + '.effect')
	fcols.append(cfg['stat'] + '.stderr')
	fcols.append(cfg['stat'] + '.or')
	fcols.append(cfg['stat'] + '.z')
	fcols.append(cfg['stat'] + '.p')
	fcols.append(cfg['stat'] + '.dir')
	chr_col = '#chr'
	pos_col = 'pos'
	rsq_col = 'rsq.unrel' if cfg['unrel'] else 'rsq'
	freq_col = 'freq.unrel' if cfg['unrel'] else 'freq'
	hwe_col = 'hwe.unrel' if cfg['unrel'] else 'hwe'
	effect_col = cfg['stat'] + '.effect'
	stderr_col = cfg['stat'] + '.stderr'
	oddsratio_col = cfg['stat'] + '.or'
	z_col = cfg['stat'] + '.z'
	p_col = cfg['stat'] + '.p'
	meta_dir_col = cfg['stat'] + '.dir'
	if cfg['tag'] is not None:
		fcols = [cfg['tag'] + '.' + x if x not in ['#chr','pos','a1','a2'] else x for x in fcols]
		rsq_col = cfg['tag'] + '.' + rsq_col
		freq_col = cfg['tag'] + '.' + freq_col
		hwe_col = cfg['tag'] + '.' + hwe_col
		effect_col = cfg['tag'] + '.' + effect_col
		stderr_col = cfg['tag'] + '.' + stderr_col
		oddsratio_col = cfg['tag'] + '.' + oddsratio_col
		z_col = cfg['tag'] + '.' + z_col
		p_col = cfg['tag'] + '.' + p_col
		meta_dir_col = cfg['tag'] + '.' + meta_dir_col

	##### GENERATE REGION LIST #####
	regions = None
	if not cfg['region_list'] is None:
		print "loading region list"
		variant_list = FileFxns.LoadCoordinates(cfg['region_list'])
		regions = list(variant_list['region'])
	elif not cfg['region'] is None:
		if len(cfg['region'].split(':')) > 1:
			variant_list = pd.DataFrame({'chr': [re.split(':|-',cfg['region'])[0]],'start': [re.split(':|-',cfg['region'])[1]],'end': [re.split(':|-',cfg['region'])[2]],'region': [cfg['region']]})
		else:
			variant_list = pd.DataFrame({'chr': [cfg['region']],'start': ['NA'],'end': ['NA'],'region': [cfg['region']]})
		variant_list['id'] = 'NA'
		regions = list(variant_list['region'])

	##### read data from file #####
	print "loading results from file"
	if regions is None:
		reader = pd.read_table(cfg['file'], sep='\t', chunksize=1000000,compression='gzip',dtype=object)
		i = 0
		for chunk in reader:
			i = i+1
			chunk = chunk[[x for x in chunk.columns if x in fcols]]
			lines = len(chunk)
			chunk = chunk[pd.notnull(chunk[p_col])]
			chunk = chunk[chunk[p_col] != 0]
			if len(chunk) > 0:
				if i == 1:
					pvals = chunk
				else:
					pvals = pvals.append(chunk)
			print "   processed " + str((i-1)*1000000 + lines) + " lines: " + str(len(chunk)) + " markers added: " + str(len(pvals.index)) + " total"
	else:
		try:
			h = subprocess.Popen(['tabix','-h',cfg['file'],'0'], stdout=subprocess.PIPE)
		except:
			usage(SystemFxns.Error("process_file " + cfg['file'] + " has incorrect format"))
		header = h.communicate()[0]
		header = header.replace("#","")
		header = header.strip()
		header = header.split()		
		reader = tabix.open(cfg['file'])
		for r in range(len(variant_list.index)):
			reg = variant_list['region'][r]
			try:
				records = reader.querys(reg)
			except:
				pass
			else:
				chunk = pd.DataFrame([record for record in records])
				chunk.columns = header
				chunk = chunk[chunk[p_col] != 'NA']
				chunk = chunk[pd.notnull(chunk[p_col])]
				chunk = chunk[chunk[p_col] != 0]
				if len(chunk) > 0:
					if r == 0:
						pvals = chunk
					else:
						pvals = pvals.append(chunk)
			print "   processed region " + str(r+1) + "/" + str(len(variant_list.index)) + " (" + str(variant_list['id'][r]) + " " + str(variant_list['region'][r]) + "): " + str(len(chunk)) + " markers added: " + str(len(pvals.index)) + " total"
		pvals[pvals == 'NA'] = float('nan')
	pvals.columns = [a.replace('#','') for a in pvals.columns]
	chr_col = 'chr'
	conv_cols = [chr_col,pos_col] + [x for x in pvals.columns if 'rsq' in x or 'hwe' in x or 'freq' in x or '.effect' in x or '.stderr' in x or '.or' in x or '.z' in x or '.p' in x]
	pvals[conv_cols] = pvals[conv_cols].convert_objects(convert_numeric=True)

	##### filter data #####
	if cfg['rsq']:
		print "filtering data for imputation quality"
		pvals = pvals[(pvals[rsq_col] >= cfg['rsq']) & (pvals[rsq_col] <= 1/cfg['rsq'])]
	if cfg['maf']:
		print "filtering data for frequency"
		pvals = pvals[(pvals[freq_col] >= cfg['maf']) & (pvals[freq_col] <= 1 - cfg['maf'])]
	if cfg['hwe']:
		print "filtering data for Hardy Weinberg p-value"
		pvals = pvals[pvals[hwe_col] > cfg['hwe']]
	if cfg['callrate']:
		print "filtering data for callrate"
		pvals = pvals[pvals[callrate_col] >= cfg['callrate']]
	if cfg['effect']:
		print "filtering data for effect estimate"
		pvals = pvals[(pvals[effect_col] <= cfg['effect']) & (pvals[effect_col] >= -1 * cfg['effect'])]
	if cfg['stderr']:
		print "filtering data for standard error"
		pvals = pvals[pvals[stderr_col] <= cfg['stderr']]
	if cfg['odds_ratio']:
		print "filtering data for odds ratio"
		pvals = pvals[(pvals[oddsratio_col] <= cfg['odds_ratio']) & (pvals[oddsratio_col] >= 1 / cfg['odds_ratio'])]
	def count_df(x):
		return len(re.findall('\\+|-',x))
	if cfg['df']:
		print "filtering data for degrees of freedom"
		pvals = pvals[pvals[meta_dir_col].apply(count_df) >= int(cfg['df']) + 1]
	print str(len(pvals)) + " markers left after filtering"
	
	if regions is None:
		l=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
		print "genomic inflation (all markers) = " + str(l)
		if freq_col in pvals:
			lA='NA'
			lB='NA'
			lC='NA'
			lD='NA'
			lE='NA'
			lE_n=len(pvals[p_col][(pvals[freq_col] < 0.01) | (pvals[freq_col] > 0.99)])
			lD_n=len(pvals[p_col][((pvals[freq_col] >= 0.01) & (pvals[freq_col] < 0.03)) | ((pvals[freq_col] <= 0.99) & (pvals[freq_col] > 0.97))])
			lC_n=len(pvals[p_col][((pvals[freq_col] >= 0.03) & (pvals[freq_col] < 0.05)) | ((pvals[freq_col] <= 0.97) & (pvals[freq_col] > 0.95))])
			lB_n=len(pvals[p_col][((pvals[freq_col] >= 0.05) & (pvals[freq_col] < 0.1)) | ((pvals[freq_col] <= 0.95) & (pvals[freq_col] > 0.9))])
			lA_n=len(pvals[p_col][(pvals[freq_col] >= 0.1) & (pvals[freq_col] <= 0.9)])
			if lE_n > 0:
				lE=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col][(pvals[freq_col] < 0.01) | (pvals[freq_col] > 0.99)].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			if lD_n > 0:
				lD=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col][((pvals[freq_col] >= 0.01) & (pvals[freq_col] < 0.03)) | ((pvals[freq_col] <= 0.99) & (pvals[freq_col] > 0.97))].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			if lC_n > 0:
				lC=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col][((pvals[freq_col] >= 0.03) & (pvals[freq_col] < 0.05)) | ((pvals[freq_col] <= 0.97) & (pvals[freq_col] > 0.95))].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			if lB_n > 0:
				lB=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col][((pvals[freq_col] >= 0.05) & (pvals[freq_col] < 0.1)) | ((pvals[freq_col] <= 0.95) & (pvals[freq_col] > 0.9))].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			if lA_n > 0:
				lA=np.median(scipy.chi2.ppf([1-x for x in pvals[p_col][(pvals[freq_col] >= 0.1) & (pvals[freq_col] <= 0.9)].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			print "genomic inflation (MAF >= 10%, n=" + str(lA_n) + ") = " + str(lA)
			print "genomic inflation (5% <= MAF < 10%, n=" + str(lB_n) + ") = " + str(lB)
			print "genomic inflation (3% <= MAF < 5%, n=" + str(lC_n) + ") = " + str(lC)
			print "genomic inflation (1% <= MAF < 3%, n=" + str(lD_n) + ") = " + str(lD)
			print "genomic inflation (MAF < 1%, n=" + str(lE_n) + ") = " + str(lE)
		else:
			if cfg['qq_strat']:
				print freq_col + " not found in data file, skipping frequency stratified qq plot"

		if cfg['qq'] and freq_col in pvals and cfg['qq_strat']:
			pvals['logp'] = -1 * np.log10(pvals[p_col])
			pvals.sort(columns=['logp'], inplace=True)
			pvals['MAF'] = 'E'
			pvals['MAF'][(pvals[freq_col] >= 0.01) & (pvals[freq_col] <= 0.99)] = 'D'
			pvals['MAF'][(pvals[freq_col] >= 0.03) & (pvals[freq_col] <= 0.97)] = 'C'
			pvals['MAF'][(pvals[freq_col] >= 0.05) & (pvals[freq_col] <= 0.95)] = 'B'
			pvals['MAF'][(pvals[freq_col] >= 0.1) & (pvals[freq_col] <= 0.9)] = 'A'
			a = np.array([])
			b = np.array([])
			c = np.array([])
			if len(pvals[pvals['MAF'] == 'E'].index) > 0:
				aa = -1 * np.log10(ro.r('ppoints(' + str(len(pvals[pvals['MAF'] == 'E'].index)) + ')'))
				aa.sort()
				bb = pvals['logp'][pvals['MAF'] == 'E']
				bb.sort()
				cc = pvals['MAF'][pvals['MAF'] == 'E']
				a = np.append(a,aa)
				b = np.append(b,bb)
				c = np.append(c,cc)
			if len(pvals[pvals['MAF'] == 'D'].index) > 0:
				aa = -1 * np.log10(ro.r('ppoints(' + str(len(pvals[pvals['MAF'] == 'D'].index)) + ')'))
				aa.sort()
				bb = pvals['logp'][pvals['MAF'] == 'D']
				bb.sort()
				cc = pvals['MAF'][pvals['MAF'] == 'D']
				a = np.append(a,aa)
				b = np.append(b,bb)
				c = np.append(c,cc)
			if len(pvals[pvals['MAF'] == 'C'].index) > 0:
				aa = -1 * np.log10(ro.r('ppoints(' + str(len(pvals[pvals['MAF'] == 'C'].index)) + ')'))
				aa.sort()
				bb = pvals['logp'][pvals['MAF'] == 'C']
				bb.sort()
				cc = pvals['MAF'][pvals['MAF'] == 'C']
				a = np.append(a,aa)
				b = np.append(b,bb)
				c = np.append(c,cc)
			if len(pvals[pvals['MAF'] == 'B'].index) > 0:
				aa = -1 * np.log10(ro.r('ppoints(' + str(len(pvals[pvals['MAF'] == 'B'].index)) + ')'))
				aa.sort()
				bb = pvals['logp'][pvals['MAF'] == 'B']
				bb.sort()
				cc = pvals['MAF'][pvals['MAF'] == 'B']
				a = np.append(a,aa)
				b = np.append(b,bb)
				c = np.append(c,cc)
			if len(pvals[pvals['MAF'] == 'A'].index) > 0:
				aa = -1 * np.log10(ro.r('ppoints(' + str(len(pvals[pvals['MAF'] == 'A'].index)) + ')'))
				aa.sort()
				bb = pvals['logp'][pvals['MAF'] == 'A']
				bb.sort()
				cc = pvals['MAF'][pvals['MAF'] == 'A']
				a = np.append(a,aa)
				b = np.append(b,bb)
				c = np.append(c,cc)
			df = ro.DataFrame({'a': ro.FloatVector(a), 'b': ro.FloatVector(b), 'MAF': ro.StrVector(c)})

			print "generating frequency stratified qq plot"
			if cfg['ext'] == 'tiff':
				grdevices.tiff(cfg['file'].replace('.gz','') + '.qq_strat.' + cfg['ext'],width=4,height=4,units="in",bg="white",compression="lzw",res=300)
			elif cfg['ext'] == 'eps':
				grdevices.postscript(cfg['file'].replace('.gz','') + '.qq_strat.' + cfg['ext'],width=4,height=4,bg="white",horizontal=False)
			else:
				grdevices.pdf(cfg['file'].replace('.gz','') + '.qq_strat.' + cfg['ext'],width=4,height=4,bg="white")
			gp = ggplot2.ggplot(df)
			pp = gp + \
					ggplot2.aes_string(x='a',y='b') + \
					ggplot2.geom_point(ggplot2.aes_string(color='MAF'), size=2) + \
					ggplot2.scale_colour_manual(values=ro.r('c("E"="#a8ddb5", "D"="#7bccc4", "C"="#4eb3d3", "B"="#2b8cbe", "A"="#08589e")'), labels=ro.r('c("E"="MAF < 1%","D"="1% <= MAF < 3%","C"="3% <= MAF < 5%","B"="5% <= MAF < 10%","A"="MAF >= 10%")')) + \
					ggplot2.geom_abline(intercept=0, slope=1, alpha=0.5) + \
					ggplot2.scale_x_continuous(ro.r('expression(Expected~~-log[10](italic(p)))')) + \
					ggplot2.scale_y_continuous(ro.r('expression(Observed~~-log[10](italic(p)))')) + \
					ggplot2.theme_bw(base_size = 12)
			pp = pp + ggplot2.theme(**{'axis.title.x': ggplot2.element_text(vjust=-0.5,size=14), 'axis.title.y': ggplot2.element_text(vjust=1,angle=90,size=14), 'legend.title': ggplot2.element_blank(), 'legend.key.height': ro.r.unit(0.1,"in"), 'legend.text': ggplot2.element_text(size=5), 'legend.key': ggplot2.element_blank(), 'legend.justification': ro.r('c(0,1)'), 'legend.position': ro.r('c(0,1)'), 'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.text': ggplot2.element_text(size=12)})
			pp.plot()
			grdevices.dev_off()

		if cfg['qq']:
			a = -1 * np.log10(ro.r('ppoints(' + str(len(pvals.index)) + ')'))
			a.sort()
			
			pvals['logp'] = -1 * np.log10(pvals[p_col]) + 0.0
			pvals.sort(columns=['logp'], inplace=True)
			
			ci_upper = -1 * np.log10(scipy.beta.ppf(0.95, range(1,len(pvals[p_col]) + 1), range(len(pvals[p_col]),0,-1)))
			ci_upper.sort()
			ci_lower = -1 * np.log10(scipy.beta.ppf(0.05, range(1,len(pvals[p_col]) + 1), range(len(pvals[p_col]),0,-1)))
			ci_lower.sort()

			df = ro.DataFrame({'a': ro.FloatVector(a), 'b': ro.FloatVector(pvals['logp']), 'ci_upper': ro.FloatVector(ci_upper), 'ci_lower': ro.FloatVector(ci_lower)})
			dftext_label = 'lambda %~~% ' + str(round(l,3))
			dftext = ro.DataFrame({'x': ro.r('Inf'), 'y': 0.5, 'lab': dftext_label})

			print "generating qq plot"
			if cfg['ext'] == 'tiff':
				grdevices.tiff(cfg['file'].replace('.gz','') + '.qq.' + cfg['ext'],width=4,height=4,units="in",bg="white",compression="lzw",res=300)
			elif cfg['ext'] == 'eps':
				grdevices.postscript(cfg['file'].replace('.gz','') + '.qq.' + cfg['ext'],width=4,height=4,bg="white",horizontal=False)
			else:
				grdevices.pdf(cfg['file'].replace('.gz','') + '.qq.' + cfg['ext'],width=4,height=4,bg="white")
			gp = ggplot2.ggplot(df)
			pp = gp + \
					ggplot2.aes_string(x='a',y='b') + \
					ggplot2.geom_ribbon(ggplot2.aes_string(x='a',ymin='ci_lower',ymax='ci_upper'), data=df, alpha=0.15, fill='black') + \
					ggplot2.geom_point(size=2) + \
					ggplot2.geom_abline(intercept=0, slope=1, alpha=0.5) + \
					ggplot2.scale_x_continuous(ro.r('expression(Expected~~-log[10](italic(p)))')) + \
					ggplot2.scale_y_continuous(ro.r('expression(Observed~~-log[10](italic(p)))')) + \
					ggplot2.theme_bw(base_size = 12) + \
					ggplot2.geom_text(ggplot2.aes_string(x='x', y='y', label='lab'), data = dftext, colour="black", vjust=0, hjust=1, size = 4, parse=ro.r('TRUE'))
			if cfg['qq_n']:
				dftext2_label = '~~~ n == ' + str(len(pvals))
				dftext2 = ro.DataFrame({'x': ro.r('Inf'), 'y': 0, 'lab': dftext2_label})
				pp = pp + ggplot2.geom_text(ggplot2.aes_string(x='x', y='y', label='lab'), data = dftext2, colour="black", vjust=0, hjust=1, size = 4, parse=ro.r('TRUE'))
			pp = pp + ggplot2.theme(**{'axis.title.x': ggplot2.element_text(vjust=-0.5,size=14), 'axis.title.y': ggplot2.element_text(vjust=1,angle=90,size=14), 'legend.position': 'none', 'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.text': ggplot2.element_text(size=12)})
			pp.plot()
			grdevices.dev_off()
		
		if cfg['plot_gc']:
			print "adjusting p-values for genomic inflation"
			pvals[p_col]=2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*pvals[p_col]) / math.sqrt(l)))

		if cfg['mht']:
			print "calculating genomic positions for manhattan plot"
			df = pvals[[chr_col,pos_col,p_col]].reset_index(drop=True)
			df.sort(columns=[chr_col,pos_col], inplace=True)
			ticks = []
			lastbase = 0
			df['gpos'] = 0
			nchr = len(list(np.unique(df[chr_col].values)))
			chrs = np.unique(df[chr_col].values)
			if cfg['color']:
				colours = ["#08306B","#41AB5D","#000000","#F16913","#3F007D","#EF3B2C","#08519C","#238B45","#252525","#D94801","#54278F","#CB181D","#2171B5","#006D2C","#525252","#A63603","#6A51A3","#A50F15","#4292C6","#00441B","#737373","#7F2704","#807DBA","#67000D"]
			else:
				colours = ["#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969"]
			if nchr == 1:
				df['gpos'] = df[pos_col]
				df['colours'] = "#000000"
				if df['gpos'].max() - df['gpos'].min() <= 1000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 100 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 10000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 1000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 100000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 10000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 200000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 20000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 300000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 30000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 400000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 40000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 500000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 50000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 600000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 60000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 700000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 70000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 800000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 80000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 900000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 90000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 1000000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 100000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 10000000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 1000000 == 0]
				elif df['gpos'].max() - df['gpos'].min() <= 100000000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 10000000 == 0]
				elif df['gpos'].max() - df['gpos'].min() > 100000000:
					ticks = [x for x in range(df['gpos'].min(),df['gpos'].max()) if x % 25000000 == 0]
			else:
				df['colours'] = "#000000"
				for i in range(len(chrs)):
					print "   processed chromosome " + str(int(chrs[i]))
					if i == 0:
						df['gpos'][df[chr_col] == chrs[i]] = df[pos_col][df[chr_col] == chrs[i]]
					else:
						lastbase = lastbase + df[pos_col][df[chr_col] == chrs[i-1]].iget(-1)
						df['gpos'][df[chr_col] == chrs[i]] = (df[pos_col][df[chr_col] == chrs[i]]) + lastbase
					ticks.append(df['gpos'][df[chr_col] == chrs[i]].iloc[(int(math.floor(len(df['gpos'][df[chr_col] == chrs[i]]))/2)) + 1])
					df['colours'][df[chr_col] == chrs[i]] = colours[int(chrs[i])]
			df['logp'] = -1 * np.log10(df[p_col])
			maxy=int(max(np.ceil(-1 * np.log10(cfg['sig'])),np.ceil(df['logp'].max())))
			if maxy > 20:
				y_breaks = range(0,maxy,5)
				y_labels = range(0,maxy,5)
			else:
				y_breaks = range(0,maxy)
				y_labels = range(0,maxy)
			rdf = ro.DataFrame({'gpos': ro.FloatVector(df['gpos']), 'logp': ro.FloatVector(df['logp']), 'colours': ro.FactorVector(df['colours'])})
			if cfg['ext'] == 'tiff':
				grdevices.tiff(cfg['file'].replace('.gz','') + '.mht.' + cfg['ext'],width=12,height=4,units="in",bg="white",compression="lzw",res=300)
			elif cfg['ext'] == 'eps':
				grdevices.postscript(cfg['file'].replace('.gz','') + '.mht.' + cfg['ext'],width=12,height=4,bg="white",horizontal=False)
			else:
				grdevices.pdf(cfg['file'].replace('.gz','') + '.mht.' + cfg['ext'],width=12,height=4,bg="white")
			print "generating manhattan plot"
			if nchr == 1:
				gp = ggplot2.ggplot(rdf)
				pp = gp + \
						ggplot2.aes_string(x='gpos',y='logp') + \
						ggplot2.geom_hline(yintercept = -1 * np.log10(cfg['sig']),colour="#B8860B", linetype=5, size = 0.25) + \
						ggplot2.geom_point(size=1.5) + \
						ggplot2.scale_x_continuous(ro.r('expression(Chromosome~~' + str(df[chr_col][0]) + '~~(kb))'),breaks=ro.FloatVector(ticks),labels=ro.Vector(["{:,}".format(x/1000) for x in ticks])) + \
						ggplot2.scale_y_continuous(ro.r('expression(-log[10](italic(p)))'),limits=ro.r('c(0,' + str(maxy) + ')')) + \
						ggplot2.theme_bw(base_size = 8) + \
						ggplot2.theme(**{'axis.title.x': ggplot2.element_text(vjust=-0.5,size=14), 'axis.title.y': ggplot2.element_text(vjust=1,angle=90,size=14), 'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.title': ggplot2.element_text(size=10), 'axis.text': ggplot2.element_text(size=8), 'legend.position': 'none', 'axis.text': ggplot2.element_text(size=12)})
				pp.plot()
			else:
				gp = ggplot2.ggplot(rdf)
				pp = gp + \
						ggplot2.aes_string(x='gpos',y='logp',colour='colours') + \
						ggplot2.geom_hline(yintercept = -1 * np.log10(cfg['sig']),colour="#B8860B", linetype=5, size = 0.25) + \
						ggplot2.geom_point(size=1.5) + \
						ggplot2.scale_colour_manual(values=ro.StrVector(colours)) + \
						ggplot2.scale_x_continuous(ro.r('expression(Chromosome)'),breaks=ro.FloatVector(ticks),labels=ro.FloatVector(chrs)) + \
						ggplot2.scale_y_continuous(ro.r('expression(-log[10](italic(p)))'),limits=ro.r('c(0,' + str(maxy) + ')')) + \
						ggplot2.theme_bw(base_size = 8) + \
						ggplot2.theme(**{'axis.title.x': ggplot2.element_text(vjust=-0.5,size=14), 'axis.title.y': ggplot2.element_text(vjust=1,angle=90,size=14), 'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.title': ggplot2.element_text(size=8), 'axis.text': ggplot2.element_text(size=6), 'legend.position': 'none', 'axis.text': ggplot2.element_text(size=12)})
				pp.plot()
			grdevices.dev_off()

		##### DETERMINE TOP REGIONS #####
		regions = []
		pvals_top = pvals.copy()
		if cfg['top'] is not None:
			print "determining top regions for regional manhattan plots"
			pvals_top.sort(columns=[p_col], inplace=True)
			while len(regions) < cfg['top']:
				region_temp = str(pvals_top['chr'].iloc[0]) + ':' + str(pvals_top['pos'].iloc[0] - 100000) + '-' + str(pvals_top['pos'].iloc[0] + 100000)
				regions.append(region_temp)
				pvals_top = pvals_top[~((pvals_top['chr'] == int(region_temp.split(':')[0])) & (pvals_top['pos'] >= int(region_temp.split(':')[1].split('-')[0])) & (pvals_top['pos'] <= int(region_temp.split(':')[1].split('-')[1])))]

		##### PRINT TOP RESULTS TO FILE #####
		print 'writing top results file'
		pvals.fillna('NA',inplace=True)
		if pvals[pvals[p_col] < cfg['pmax']].shape[0] > 100:
			pvals_out = pvals[pvals[p_col] < cfg['pmax']].sort(columns=[p_col])
		else:
			pvals_out = pvals.sort(columns=[p_col]).head(100)
		pvals_out.rename(columns={'chr':'#chr'},inplace=True)
		pvals_out[p_col] = pvals_out[p_col].map('{:,.4e}'.format)
		if 'logp' in pvals_out:
			pvals_out.drop('logp',axis=1,inplace=True)
		if 'MAF' in pvals_out:
			pvals_out.drop('MAF',axis=1,inplace=True)
		pvals_out.fillna('NA',inplace=True)
		pvals_out.to_csv(cfg['file'].replace('.gz','') + '.top_results', header=True, index=False, sep="\t")

	else:
		if cfg['set_gc'] is not None:
			print "adjusting p-values for genomic inflation"
			pvals[p_col]=2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*pvals[p_col]) / math.sqrt(cfg['set_gc'])))

	if len(regions) > 0:
		print "generating regional plots"
		home_dir = os.path.expanduser("~")
		for reg in regions:
			chr = int(reg.split(':')[0])
			start = int(reg.split(':')[1].split('-')[0])
			end = int(reg.split(':')[1].split('-')[1])
			pvals_region = pvals[(pvals['chr'] == chr) & (pvals['pos'] >= start) & (pvals['pos'] <= end)]
			pvals_region['MarkerName'] = pvals_region['marker'].values
			pvals_region['MarkerName'] = pvals_region[['chr','pos','MarkerName']].apply(lambda row: str(row[0]) + ':' + str(row[1]) if not 'rs' in row[2] else row[2], axis=1)
			pvals_region['P-value'] = pvals_region[p_col].values
			pvals_region.sort(columns=['P-value'], inplace=True)
			pvals_region = pvals_region[['MarkerName','P-value','pos']]
			pvals_region.to_csv(cfg['file'].replace('.gz','') + '.rgnl.chr' + reg.replace(':','bp') + '.plotdata',header=True, index=False, sep='\t')
			cmd = cfg['locuszoom'] + ' --metal ' + cfg['file'].replace('.gz','') + '.rgnl.chr' + reg.replace(':','bp') + '.plotdata --chr ' + str(chr) + ' --start ' + str(start) + ' --end ' + str(end) + ' --plotonly --cache None --prefix ' + cfg['file'].replace('.gz','')
			if cfg['lz_pop'] is not None:
				cmd = cmd + ' --source ' + cfg['lz_source'] + ' --build ' + cfg['lz_build'] + ' --pop ' + cfg['lz_pop'] 
			else:
				cmd = cmd + ' --no-ld'
			try:
				pr = subprocess.Popen(cmd,shell=True)
				pr.wait()
			except KeyboardInterrupt:
				kill_all(pr.pid)
				print Highlight("process terminated by user")
				sys.exit(1)
			print "   regional plot " + reg + " finished"

	print "process complete"
