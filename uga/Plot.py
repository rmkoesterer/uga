import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import math
from Messages import Error
from scipy.stats import chi2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.lib.ggplot2 as ggplot2

rplot = ro.r('plot')
rhist = ro.r('hist')
grdevices = importr('grDevices')
rgrid = importr('grid')

def Plot(data, out, qq = False, manhattan = False, color = True, ext = 'tiff', chr = '#chr', pos = 'pos', p = 'marker.p', rsq = None, freq = None, hwe = None, rsq_thresh = None, freq_thresh = None, hwe_thresh = None):

	print "   ... arguments"
	for arg in locals().keys():
		if not locals()[arg] in [None, False]:
			print "      {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	fcols = [chr,pos,p]
	if rsq:
		fcols.append(rsq)
	if freq:
		fcols.append(freq)
	if hwe:
		fcols.append(hwe)
	
	##### read data from file #####
	pvals = pd.DataFrame(columns=fcols)
	reader = pd.read_table(data, sep='\t', chunksize=1000000,compression='gzip',usecols=fcols)
	i = 0
	for chunk in reader:
		i = i+1
		lines = len(chunk)
		chunk = chunk[pd.notnull(chunk[p])]
		chunk = chunk[chunk[p] != 0]
		if len(chunk) > 0:
			pvals = pvals.append(chunk)
		print "   ... loaded " + str(i*1000000) + " lines: " + str(len(chunk)) + " markers added: " + str(len(pvals.index)) + " total"

	##### filter data #####
	if rsq and rsq_thresh:
		print "   ... filtering data for imputation quality"
		pvals = pvals[(pvals[rsq] >= rsq_thresh) & (pvals[rsq] <= 1/rsq_thresh)]
	if freq and freq_thresh:
		print "   ... filtering data for frequency"
		pvals = pvals[(pvals[freq] >= freq_thresh) & (pvals[freq] <= 1 - freq_thresh)]
	if hwe and hwe_thresh:
		print "   ... filtering data for Hardy Weinberg p-value"
		pvals = pvals[pvals[hwe] > hwe_thresh]
	print "   ... " + str(len(pvals)) + " markers left after filtering"
	
	print "   ... calculating genomic inflation factor"
	l=np.median(chi2.ppf([1-x for x in pvals[p].tolist()], df=1))/0.455
	print "   ... lambda = " + str(l)
	
	a = -1 * np.log10(ro.r('ppoints(' + str(len(pvals.index)) + ')'))
	a.sort()
	b = -1 * np.log10(pvals[p])
	b.sort()
	
	df = ro.DataFrame({'a': ro.FloatVector(a), 'b': ro.FloatVector(b)})
	dftext = ro.DataFrame({'x': a[-1] * 0.10, 'y': b.irow(-1) * 0.80, 'lab': 'lambda==' + str(l)})

	if qq:
		print "   ... generating qq plot"
		if ext == 'tiff':
			grdevices.tiff(out + '.qq.' + ext,width=4,height=4,units="in",bg="white",compression="lzw",res=300)
		elif ext == 'eps':
			grdevices.postscript(out + '.qq.' + ext,width=4,height=4,bg="white",horizontal=False)
		else:
			grdevices.pdf(out + '.qq.' + ext,width=4,height=4,bg="white")
		gp = ggplot2.ggplot(df)
		pp = gp + \
				ggplot2.aes_string(x='a',y='b') + \
				ggplot2.geom_point(size=1.5) + \
				ggplot2.geom_abline(intercept=0, linetype=2)	 + \
				ggplot2.scale_x_continuous(ro.r('expression(Expected~~-log[10](italic(p)))')) + \
				ggplot2.scale_y_continuous(ro.r('expression(Observed~~-log[10](italic(p)))')) + \
				ggplot2.theme_bw(base_size = 8) + \
				ggplot2.geom_text(ggplot2.aes_string(x='x', y='y', label='lab'), data = dftext, colour="black", size = 2, parse=ro.r('TRUE')) + \
				ggplot2.theme(**{'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.title': ggplot2.element_text(size=8), 'axis.text': ggplot2.element_text(size=6)})
		pp.plot()
		grdevices.dev_off()
	
	if manhattan:
		print "   ... ordering by genomic position for manhattan plot"
		df = pvals[[chr,pos,p]].reset_index(drop=True)
		df.sort(columns=[chr,pos], inplace=True)
		ticks = []
		lastbase = 0
		df['gpos'] = 0
		nchr = len(list(np.unique(df[chr].values)))
		chrs = np.unique(df[chr].values)
		if color:
			colours = ["#08306B","#41AB5D","#000000","#F16913","#3F007D","#EF3B2C","#08519C","#238B45","#252525","#D94801","#54278F","#CB181D","#2171B5","#006D2C","#525252","#A63603","#6A51A3","#A50F15","#4292C6","#00441B","#737373","#7F2704","#807DBA","#67000D"]
		else:
			colours = ["#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969","#000000","#696969"]
		if nchr == 1:
			df['gpos'] = df[pos]
			ticks.append(df['gpos'][int((math.floor(len(df['gpos']))/2)) + 1])
		else:
			df['colours'] = "#000000"
			for i in range(len(chrs)):
				print "   ... calculating genomic positions for chromosome " + str(int(chrs[i]))
				if i == 0:
					df['gpos'][df[chr] == chrs[i]] = df[pos][df[chr] == chrs[i]]
				else:
					lastbase = lastbase + df[pos][df[chr] == chrs[i-1]].iget(-1)
					df['gpos'][df[chr] == chrs[i]] = (df[pos][df[chr] == chrs[i]]) + lastbase
				ticks.append(df['gpos'][df[chr] == chrs[i]].iloc[(int(math.floor(len(df['gpos'][df[chr] == chrs[i]]))/2)) + 1])
				df['colours'][df[chr] == chrs[i]] = colours[int(chrs[i])]
		df['logp'] = -1 * np.log10(df[p])
		maxy=int(max(np.ceil(-1 * np.log10(0.000000054)),np.ceil(df['logp'].max())))
		if maxy > 20:
			y_breaks = range(0,maxy,by=5)
			y_labels = range(0,maxy,by=5)
		else:
			y_breaks = range(0,maxy)
			y_labels = range(0,maxy)
		rdf = ro.DataFrame({'gpos': ro.FloatVector(df['gpos']), 'logp': ro.FloatVector(df['logp']), 'colours': ro.FactorVector(df['colours'])})
		if ext == 'tiff':
			grdevices.tiff(out + '.mht.' + ext,width=12,height=4,units="in",bg="white",compression="lzw",res=300)
		elif ext == 'eps':
			grdevices.postscript(out + '.mht.' + ext,width=12,height=4,bg="white",horizontal=False)
		else:
			grdevices.pdf(out + '.mht.' + ext,width=12,height=4,bg="white")
		print "   ... generating manhattan plot"
		if nchr == 1:
			gp = ggplot2.ggplot(rdf)
			pp = gp + \
					ggplot2.aes_string(x='gpos',y='logp') + \
					ggplot2.geom_point(size=1.5) + \
					ggplot2.geom_hline(yintercept = -1 * np.log10(0.000000054),colour="#B8860B") + \
					ggplot2.scale_x_continuous(ro.r('expression(Chromosome)')) + \
					ggplot2.scale_y_continuous(ro.r('expression(-log[10](italic(p)))'),limits=ro.r('c(0,' + str(maxy) + ')')) + \
					ggplot2.theme_bw(base_size = 8) + \
					ggplot2.theme(**{'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.title': ggplot2.element_text(size=8), 'axis.text': ggplot2.element_text(size=6), 'legend.position': 'none'})
			pp.plot()
		else:
			gp = ggplot2.ggplot(rdf)
			pp = gp + \
					ggplot2.aes_string(x='gpos',y='logp',colour='colours') + \
					ggplot2.geom_point(size=1.5) + \
					ggplot2.geom_hline(yintercept = -1 * np.log10(0.000000054),colour="#B8860B") + \
					ggplot2.scale_colour_manual(values=ro.StrVector(colours)) + \
					ggplot2.scale_x_continuous(ro.r('expression(Chromosome)'),breaks=ro.FloatVector(ticks),labels=ro.FloatVector(chrs)) + \
					ggplot2.scale_y_continuous(ro.r('expression(-log[10](italic(p)))'),limits=ro.r('c(0,' + str(maxy) + ')')) + \
					ggplot2.theme_bw(base_size = 8) + \
					ggplot2.theme(**{'panel.background': ggplot2.element_blank(), 'panel.border': ggplot2.element_blank(), 'panel.grid.minor': ggplot2.element_blank(), 'panel.grid.major': ggplot2.element_blank(), 'axis.line': ro.r('element_line(colour="black")'), 'axis.title': ggplot2.element_text(size=8), 'axis.text': ggplot2.element_text(size=6), 'legend.position': 'none'})
			pp.plot()
		grdevices.dev_off()
	
	"   ... process complete"
					