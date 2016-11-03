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
import numpy as np
import scipy.stats as scipy
import Parse
import pysam
import math
import Process
import readline
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import logging
import re
pd.options.mode.chained_assignment = None
pandas2ri.activate()

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunSnvplot")

def RunSnvplot(args):
	cfg = Parse.generate_snvplot_cfg(args)
	Parse.print_snvplot_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	ro.r('suppressMessages(library(ggplot2))')
	ro.r('suppressMessages(library(grid))')
	ro.r('suppressMessages(library(RColorBrewer))')

	handle=pysam.TabixFile(filename=cfg['file'],parser=pysam.asVCF())
	header = [x for x in handle.header]
	skip_rows = len(header)-1
	cols = header[-1].split()
	pcols = cfg['pcol'].split(',')
	cols_extract = [cfg['chrcol'],cfg['bpcol']] + pcols
	if cfg['qq_strat_freq']:
		if cfg['freqcol'] not in cols:
			print Process.Error("frequency column " + cfg['freqcol'] + " not found, unable to proceed with frequency stratified plots").out
			return 1
		else:
			cols_extract = cols_extract + [cfg['freqcol']]
			print "frequency column " + cfg['freqcol'] + " found"
	if cfg['qq_strat_mac']:
		if cfg['maccol'] not in cols:
			print Process.Error("minor allele count column " + cfg['maccol'] + " not found, unable to proceed with minor allele count stratified plots").out
			return 1
		else:
			cols_extract = cols_extract + [cfg['maccol']]
			print "minor allele count column " + cfg['maccol'] + " found"

	print "importing data"
	r = pd.read_table(cfg['file'],sep='\t',skiprows=skip_rows,usecols=cols_extract,compression='gzip')
	print str(r.shape[0]) + " total variants found"

	for pcol in pcols:
		print "plotting p-values for column " + pcol + " ..."
		extract_cols = [cfg['chrcol'],cfg['bpcol'],pcol]
		if cfg['freqcol'] in r:
			extract_cols = extract_cols + [cfg['freqcol']]
		if cfg['maccol'] in r:
			extract_cols = extract_cols + [cfg['maccol']]
		results = r[extract_cols]
		results.dropna(inplace=True)
		results = results[(results[pcol] > 0) & (results[pcol] <= 1)].reset_index(drop=True)
		print "   " + str(results.shape[0]) + " variants with plottable p-values"

		results['logp'] = -1 * np.log10(results[pcol]) + 0.0

		ro.globalenv['results'] = results
		l = np.median(scipy.chi2.ppf([1-x for x in results[pcol].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
		# in R: median(qchisq(results$p, df=1, lower.tail=FALSE))/qchisq(0.5,1)
		print "   genomic inflation (all variants) = " + str(l)

		if cfg['qq']:
			print "   generating standard qq plot"
			print "   minimum p-value: " + str(np.min(results[pcol]))
			a = -1 * np.log10(ro.r('ppoints(' + str(len(results.index)) + ')'))
			a.sort()
			results.sort_values(by=['logp'], inplace=True)
			print "   maximum -1*log10(p-value): " + str(np.max(results['logp']))

			ci_upper = -1 * np.log10(scipy.beta.ppf(0.95, range(1,len(results[pcol]) + 1), range(len(results[pcol]),0,-1)))
			ci_upper.sort()
			ci_lower = -1 * np.log10(scipy.beta.ppf(0.05, range(1,len(results[pcol]) + 1), range(len(results[pcol]),0,-1)))
			ci_lower.sort()
			
			ro.globalenv['df'] = ro.DataFrame({'a': ro.FloatVector(a), 'b': ro.FloatVector(results['logp']), 'ci_lower': ro.FloatVector(ci_lower), 'ci_upper': ro.FloatVector(ci_upper)})
			dftext_label = 'lambda %~~% ' + str(round(l,3))
			ro.globalenv['dftext'] = ro.DataFrame({'x': ro.r('Inf'), 'y': 0.5, 'lab': dftext_label})

			if cfg['ext'] == 'tiff':
				ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq.tiff')
			elif cfg['ext'] == 'png':
				ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq.png')
			elif cfg['ext'] == 'eps':
				ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.eps')
			else:
				ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.pdf')
			ro.r("""
				gp<-ggplot(df)
				pp<-gp + 
					aes_string(x='a',y='b') +
					geom_ribbon(aes_string(x='a',ymin='ci_lower',ymax='ci_upper'), data=df, alpha=0.25, fill='black') + 
					geom_point(size=2) +
					geom_abline(intercept=0, slope=1, alpha=0.5) + 
					scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
					scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
					coord_fixed() +
					theme_bw(base_size = 12) + 
					geom_text(aes_string(x='x', y='y', label='lab'), data = dftext, colour="black", vjust=0, hjust=1, size = 4, parse=TRUE) +
					theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.position = 'none', 
						panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
						panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
				%s
				""" % (ggsave))

			if np.max(results['logp']) > cfg['crop']:
				print "   generating cropped standard qq plot"
				ro.r('df$b[df$b > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
				ro.r('df$shape<-0')
				ro.r('df$shape[df$b == ' + str(cfg['crop']) + ']<-1')
				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq.cropped.tiff')
				elif cfg['ext'] == 'png':
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq.cropped.png')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.cropped.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.cropped.pdf')
				ro.r("""
					gp<-ggplot(df)
					pp<-gp + 
						aes_string(x='a',y='b') +
						geom_ribbon(aes_string(x='a',ymin='ci_lower',ymax='ci_upper'), data=df, alpha=0.25, fill='black') + 
						geom_point(aes(shape=factor(shape)),size=2) +
						geom_abline(intercept=0, slope=1, alpha=0.5) + 
						scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
						scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
						coord_fixed() +
						theme_bw(base_size = 12) + 
						geom_text(aes_string(x='x', y='y', label='lab'), data = dftext, colour="black", vjust=0, hjust=1, size = 4, parse=TRUE) +
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.position = 'none', 
							panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
							panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
					%s
					""" % (ggsave))

		def ppoints(n, a):
			try:
				n = np.float(len(n))
			except TypeError:
				n = np.float(n)
			return (np.arange(n) + 1 - a)/(n + 1 - 2*a)

		if cfg['qq_strat_freq']:
			print "   generating frequency stratified qq plot"

			strat_ticks = np.sort([np.float(x) for x in cfg['freq_ticks'].split(',')])
			results['UGA___QQ_BIN___'] = 0
			for i in xrange(len(strat_ticks)):
				results.loc[(results[cfg['freqcol']] >= strat_ticks[i]) & (results[cfg['freqcol']] <= 1-strat_ticks[i]),'UGA___QQ_BIN___'] = i+1
			bin_values = results['UGA___QQ_BIN___'].value_counts()
			for i in xrange(len(strat_ticks)+1):
				if i not in bin_values.index:
					bin_values[i] = 0
			counts = pd.DataFrame(bin_values)
			counts['lambda'] = np.nan
			results['description'] = 'NA'
			for i in xrange(len(strat_ticks)+1):
				if counts.loc[i,'UGA___QQ_BIN___'] > 0:
					counts.loc[i,'lambda'] = np.median(scipy.chi2.ppf([1-x for x in results[pcol][results['UGA___QQ_BIN___'] == i].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				else:
					counts.loc[i,'lambda'] = np.nan
				if i == 0:
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "(0," + str(strat_ticks[i]) + ") ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAF (0," + str(strat_ticks[i]) + "): n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
				elif i < len(strat_ticks):
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "[" + str(strat_ticks[i-1]) + "," + str(strat_ticks[i]) + ") ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAF [" + str(strat_ticks[i-1]) + "," + str(strat_ticks[i]) + "): n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
				else:
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "[" + str(strat_ticks[i-1]) + ",0.5]  ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAF [" + str(strat_ticks[i-1]) + ",0.5]: n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
			results.sort_values(['UGA___QQ_BIN___','logp'],inplace=True)
			results['expected'] = 0
			for i in counts.index:
				if counts.loc[i,'UGA___QQ_BIN___'] > 0:
					results.loc[results['UGA___QQ_BIN___'] == i,'expected'] = np.sort(-1 * np.log10(ppoints(len(results.loc[results['UGA___QQ_BIN___'] == i,'expected']),0)))
			ro.globalenv['df'] = ro.DataFrame({'expected': ro.FloatVector(results['expected']), 'logp': ro.FloatVector(results['logp']), 'UGA___QQ_BIN___': ro.IntVector(results['UGA___QQ_BIN___']), 'description': ro.StrVector(results['description'])})
			ro.r("df<-df[order(df$UGA___QQ_BIN___),]")
			ro.r("df$description<-ordered(df$description,levels=unique(df$description))")

			if cfg['ext'] == 'tiff':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.tiff')
			elif cfg['ext'] == 'png':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.png')
			elif cfg['ext'] == 'eps':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.eps')
			else:
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.pdf')
			ro.r("""
				gp<-ggplot(df, aes_string(x='expected',y='logp')) +
					geom_point(aes_string(color='description'), size=2) +
					scale_colour_manual(values=colorRampPalette(brewer.pal(9,"Blues"))(length(unique(df$description))+2)[3:(length(unique(df$description))+2)]) +
					geom_abline(intercept=0, slope=1, alpha=0.5) + 
					scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
					scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
					coord_fixed() +
					theme_bw(base_size = 12) + 
					theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
						legend.key.height = unit(0.1,"in"), legend.text = element_text(size=6), legend.key = element_blank(), legend.justification = c(0,1), 
						legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
						panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
				%s
				""" % (ggsave))

			if np.max(results['logp']) > cfg['crop']:
				print "   generating cropped frequency stratified qq plot"
				ro.r('df$logp[df$logp > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
				ro.r('df$shape<-0')
				ro.r('df$shape[df$logp == ' + str(cfg['crop']) + ']<-1')
				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.cropped.tiff')
				elif cfg['ext'] == 'png':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.cropped.png')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.cropped.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_freq.cropped.pdf')
				ro.r("""
					gp<-ggplot(df, aes_string(x='expected',y='logp')) +
						geom_point(aes(shape=factor(shape), color=description), size=2) +
						scale_colour_manual(values=colorRampPalette(brewer.pal(9,"Blues"))(length(unique(df$description))+2)[3:(length(unique(df$description))+2)]) +
						geom_abline(intercept=0, slope=1, alpha=0.5) + 
						scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
						scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
						coord_fixed() +
						theme_bw(base_size = 12) + 
						guides(shape=FALSE) + 
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
							legend.key.height = unit(0.1,"in"), legend.text = element_text(size=6), legend.key = element_blank(), legend.justification = c(0,1), 
							legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
							panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
					%s
					""" % (ggsave))

		if cfg['qq_strat_mac']:
			print "   generating minor allele count stratified qq plot"

			strat_ticks = np.sort([np.float(x) for x in cfg['mac_ticks'].split(',')])
			results['UGA___QQ_BIN___'] = 0
			for i in xrange(len(strat_ticks)):
				results.loc[results[cfg['maccol']] >= strat_ticks[i],'UGA___QQ_BIN___'] = i+1
			bin_values = results['UGA___QQ_BIN___'].value_counts()
			for i in xrange(len(strat_ticks)+1):
				if i not in bin_values.index:
					bin_values[i] = 0
			counts = pd.DataFrame(bin_values)
			counts['lambda'] = 0
			results['description'] = 'NA'
			for i in np.sort(counts.index):
				if counts.loc[i,'UGA___QQ_BIN___'] > 0:
					counts.loc[i,'lambda'] = np.median(scipy.chi2.ppf([1-x for x in results[pcol][results['UGA___QQ_BIN___'] == i].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				else:
					counts.loc[i,'lambda'] = np.nan
				if i == 0:
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "(0," + str(int(strat_ticks[i])) + ") ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAC (0," + str(int(strat_ticks[i])) + "): n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
				elif i < len(strat_ticks):
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "[" + str(int(strat_ticks[i-1])) + "," + str(int(strat_ticks[i])) + ") ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAC [" + str(int(strat_ticks[i-1])) + "," + str(int(strat_ticks[i])) + "): n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
				else:
					results.loc[results['UGA___QQ_BIN___'] == i,'description'] = "[" + str(int(strat_ticks[i-1])) + ",...] ~" + str(round(counts.loc[i,'lambda'],3))
					print "   MAC [" + str(int(strat_ticks[i-1])) + ",...]: n=" + str(np.int(counts.loc[i,'UGA___QQ_BIN___'])) + ", lambda=" + str(counts.loc[i,'lambda'])
			results.sort_values(['UGA___QQ_BIN___','logp'],inplace=True)
			results['expected'] = 0
			for i in counts.index:
				results.loc[results['UGA___QQ_BIN___'] == i,'expected'] = np.sort(-1 * np.log10(ppoints(len(results.loc[results['UGA___QQ_BIN___'] == i,'expected']),0)))

			ro.globalenv['df'] = ro.DataFrame({'expected': ro.FloatVector(results['expected']), 'logp': ro.FloatVector(results['logp']), 'UGA___QQ_BIN___': ro.IntVector(results['UGA___QQ_BIN___']), 'description': ro.StrVector(results['description'])})
			ro.r("df<-df[order(df$UGA___QQ_BIN___),]")
			ro.r("df$description<-ordered(df$description,levels=unique(df$description))")

			if cfg['ext'] == 'tiff':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.tiff')
			elif cfg['ext'] == 'png':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.png')
			elif cfg['ext'] == 'eps':
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.eps')
			else:
				ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.pdf')
			ro.r("""
				gp<-ggplot(df, aes_string(x='expected',y='logp')) +
					geom_point(aes_string(color='description'), size=2) +
					scale_colour_manual(values=colorRampPalette(brewer.pal(9,"Blues"))(length(unique(df$description))+2)[3:(length(unique(df$description))+2)]) +
					geom_abline(intercept=0, slope=1, alpha=0.5) + 
					scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
					scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
					coord_fixed() +
					theme_bw(base_size = 12) + 
					theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
						legend.key.height = unit(0.1,"in"), legend.text = element_text(size=6), legend.key = element_blank(), legend.justification = c(0,1), 
						legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
						panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
				%s
				""" % (ggsave))
        
			if np.max(results['logp']) > cfg['crop']:
				print "   generating cropped frequency stratified qq plot"
				ro.r('df$logp[df$logp > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
				ro.r('df$shape<-0')
				ro.r('df$shape[df$logp == ' + str(cfg['crop']) + ']<-1')
				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.cropped.tiff')
				elif cfg['ext'] == 'png':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.cropped.png')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.cropped.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat_mac.cropped.pdf')
				ro.r("""
					gp<-ggplot(df, aes_string(x='expected',y='logp')) +
						geom_point(aes(shape=factor(shape), color=description), size=2) +
						scale_colour_manual(values=colorRampPalette(brewer.pal(9,"Blues"))(length(unique(df$description))+2)[3:(length(unique(df$description))+2)]) +
						geom_abline(intercept=0, slope=1, alpha=0.5) + 
						scale_x_discrete(expression(Expected~~-log[10](italic(p)))) +
						scale_y_discrete(expression(Observed~~-log[10](italic(p)))) +
						coord_fixed() +
						theme_bw(base_size = 12) + 
						guides(shape=FALSE) + 
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
							legend.key.height = unit(0.1,"in"), legend.text = element_text(size=6), legend.key = element_blank(), legend.justification = c(0,1), 
							legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
							panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
					%s
					""" % (ggsave))

		if cfg['mht']:
			print "   generating standard manhattan plot"
			print "   minimum p-value: " + str(np.min(results[pcol]))
			print "   maximum -1*log10(p-value): " + str(np.max(results['logp']))
			if cfg['gc'] and l > 1:
				print "   adjusting p-values for genomic inflation for p-value column " + pcol
				results[pcol]=2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*results[pcol]) / math.sqrt(l)))
				print "   minimum post-gc adjustment p-value: " + str(np.min(results[pcol]))
				print "   maximum post-gc adjustment -1*log10(p-value): " + str(np.max(results['logp']))
			else:
				print "   skipping genomic inflation correction"

			print "   calculating genomic positions"
			results.sort_values(by=[cfg['chrcol'],cfg['bpcol']], inplace=True)
			ticks = []
			lastbase = 0
			results['gpos'] = 0
			nchr = len(list(np.unique(results[cfg['chrcol']].values)))
			chrs = np.unique(results[cfg['chrcol']].values)
			if cfg['color']:
				colours = ["#08306B","#41AB5D","#000000","#F16913","#3F007D","#EF3B2C","#08519C","#238B45","#252525","#D94801","#54278F","#CB181D","#2171B5","#006D2C","#525252","#A63603","#6A51A3","#A50F15","#4292C6","#00441B","#737373","#7F2704","#807DBA","#67000D"]
			else:
				colours = ["#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3"]
			if nchr == 1:
				results['gpos'] = results[cfg['bpcol']]
				results['colours'] = "#08589e"
				if results['gpos'].max() - results['gpos'].min() <= 1000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 100 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 10000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 1000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 100000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 10000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 200000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 20000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 300000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 30000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 400000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 40000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 500000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 50000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 600000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 60000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 700000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 70000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 800000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 80000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 900000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 90000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 1000000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 100000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 10000000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 1000000 == 0]
				elif results['gpos'].max() - results['gpos'].min() <= 100000000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 10000000 == 0]
				elif results['gpos'].max() - results['gpos'].min() > 100000000:
					ticks = [x for x in range(results['gpos'].min(),results['gpos'].max()) if x % 25000000 == 0]
			else:
				results['colours'] = "#000000"
				for i in range(len(chrs)):
					print "      processed chromosome " + str(int(chrs[i]))
					if i == 0:
						results.loc[results[cfg['chrcol']] == chrs[i],'gpos'] = results.loc[results[cfg['chrcol']] == chrs[i],cfg['bpcol']]
					else:
						lastbase = lastbase + results.loc[results[cfg['chrcol']] == chrs[i-1],cfg['bpcol']].iloc[-1]
						results.loc[results[cfg['chrcol']] == chrs[i],'gpos'] = (results.loc[results[cfg['chrcol']] == chrs[i],cfg['bpcol']]) + lastbase
					if results.loc[results[cfg['chrcol']] == chrs[i]].shape[0] > 1:
						ticks.append(results.loc[results[cfg['chrcol']] == chrs[i],'gpos'].iloc[0] + (results.loc[results[cfg['chrcol']] == chrs[i],'gpos'].iloc[-1] - results.loc[results[cfg['chrcol']] == chrs[i],'gpos'].iloc[0])/2)
					else:
						ticks.append(results.loc[results[cfg['chrcol']] == chrs[i],'gpos'].iloc[0])
					results.loc[results[cfg['chrcol']] == chrs[i],'colours'] = colours[int(chrs[i])]
			results['logp'] = -1 * np.log10(results[pcol])
			if results.shape[0] >= 1000000:
				sig = 5.4e-8
			else:
				sig = 0.05 / results.shape[0]
			print "   significance level set to p-value = " + str(sig) + " (-1*log10(p-value) = " + str(-1 * np.log10(sig)) + ")"
			print "   " + str(len(results[pcol][results[pcol] <= sig])) + " genome wide significant variants"
			chr = results[cfg['chrcol']][0]
			maxy=int(max(np.ceil(-1 * np.log10(sig)),np.ceil(results['logp'].max())))
			if maxy > 20:
				y_breaks = range(0,maxy,5)
				y_labels = range(0,maxy,5)
			else:
				y_breaks = range(0,maxy)
				y_labels = range(0,maxy)
			ro.globalenv['df'] = ro.DataFrame({'gpos': ro.FloatVector(results['gpos']), 'logp': ro.FloatVector(results['logp']), 'colours': ro.FactorVector(results['colours'])})
			ro.globalenv['ticks'] = ro.FloatVector(ticks)
			ro.globalenv['labels'] = ro.Vector(["{:,}".format(x/1000) for x in ticks])
			ro.globalenv['colours'] = ro.StrVector(colours)
			ro.globalenv['chrs'] = ro.FloatVector(chrs)

			print "   generating manhattan plot"
			if cfg['ext'] == 'tiff':
				ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.tiff')
			elif cfg['ext'] == 'png':
				ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.png')
			elif cfg['ext'] == 'eps':
				ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.eps')
			else:
				ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.pdf')
			if nchr == 1:
				ro.r("""
					gp<-ggplot(df, aes_string(x='gpos',y='logp')) +
						geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
						geom_point(size=1.5) + 
						scale_x_continuous(expression(Chromosome~~%d~~(kb))'),breaks=ticks,labels=labels) + \
						scale_y_continuous(expression(-log[10](italic(p))),breaks=seq(0,%d,1),limits=c(0,%d)) + \
						theme_bw(base_size = 8) + \
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
								panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
								panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=10), 
								axis.text = element_text(size=12), legend.position = 'none')
					%s
					""" % (sig, chr, maxy, maxy, ggsave))
			else:
				ro.r("""
					gp = ggplot(df, aes_string(x='gpos',y='logp',colour='colours')) + 
						geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
						geom_point(size=1.5) + 
						scale_colour_manual(values=colours) + 
						scale_x_continuous(expression(Chromosome),breaks=ticks,labels=chrs) + 
						scale_y_continuous(expression(-log[10](italic(p))),breaks=seq(0,%d,1),limits=c(0,%d)) + 
						theme_bw(base_size = 8) + 
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
								panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
								panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=10), 
								axis.text = element_text(size=12), legend.position = 'none')
					%s
					""" % (sig, maxy, maxy, ggsave))

			if maxy > cfg['crop']:
				maxy = cfg['crop']
				ro.r('df$logp[df$logp > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
				ro.r('df$shape<-0')
				ro.r('df$shape[df$logp == ' + str(cfg['crop']) + ']<-1')
				print "   generating cropped manhattan plot"
				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.cropped.tiff')
				elif cfg['ext'] == 'png':
					ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,units="in",bg="white",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.cropped.png')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.cropped.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=gp,width=16,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.cropped.pdf')
				if nchr == 1:
					ro.r("""
						gp<-ggplot(df, aes_string(x='gpos',y='logp')) +
							geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
							geom_point(aes(shape=factor(shape)),size=1.5) + 
							scale_x_continuous(expression(Chromosome~~%d~~(kb))'),breaks=ticks,labels=labels) + 
							scale_y_continuous(expression(-log[10](italic(p))),breaks=seq(0,%d,1),limits=c(0,%d)) + 
							theme_bw(base_size = 8) + 
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
									panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
									panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=10), 
									axis.text = element_text(size=12), legend.position = 'none')
						%s
						""" % (sig, chr, maxy, maxy, ggsave))
				else:
					ro.r("""
						gp = ggplot(df, aes_string(x='gpos',y='logp',colour='colours')) + 
							geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
							geom_point(aes(shape=factor(shape)),size=1.5) + 
							scale_colour_manual(values=colours) + 
							scale_x_continuous(expression(Chromosome),breaks=ticks,labels=chrs) + 
							scale_y_continuous(expression(-log[10](italic(p))),breaks=seq(0,%d,1),limits=c(0,%d)) + 
							theme_bw(base_size = 8) + 
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
									panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
									panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=8), 
									axis.text = element_text(size=12), legend.position = 'none')
						%s
						""" % (sig, maxy, maxy, ggsave))

	print "process complete"
	return 0
