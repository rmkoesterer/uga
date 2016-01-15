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
pandas2ri.activate()

logging.basicConfig(format='%(asctime)s - %(processName)s - %(name)s - %(message)s',level=logging.DEBUG)
logger = logging.getLogger("RunSnvplot")

def RunSnvgroupplot(args):
	cfg = Parse.generate_snvgroupplot_cfg(args)
	Parse.print_snvgroupplot_options(cfg)

	if not cfg['debug']:
		logging.disable(logging.CRITICAL)

	ro.r('suppressMessages(library(ggplot2))')
	ro.r('suppressMessages(library(grid))')

	handle=pysam.TabixFile(filename=cfg['file'],parser=pysam.asVCF())
	header = [x for x in handle.header]
	skip_rows = len(header)-1
	cols = header[-1].split()
	pcols = cfg['pcol'].split(',')
	if cfg['qq_strat']:
		if cfg['cmaccol'] not in cols:
			cols_extract = ['#chr','start','end','id'] + pcols
			print Process.Error("minor allele count column " + cfg['cmaccol'] + " not found, unable to proceed with minor allele count stratified plots").out
			return 1
		else:
			cols_extract = ['#chr','start','end','id',cfg['cmaccol']] + pcols
			print "minor allele count column " + cfg['cmaccol'] + " found"
	else:
		cols_extract = ['#chr','start','end','id'] + pcols
	print "importing data"
	r = pd.read_table(cfg['file'],sep='\t',skiprows=skip_rows,usecols=cols_extract,compression='gzip')
	print str(r.shape[0]) + " total groups found"
	for pcol in pcols:
		print "plotting p-values for column " + pcol + " ..."
		results = r[['#chr','start','end','id',cfg['cmaccol'],pcol]] if cfg['cmaccol'] in cols else r[['#chr','start','end','id',pcol]]
		results.dropna(inplace=True)
		results = results[(results[pcol] > 0) & (results[pcol] <= 1)]
		results.reset_index(drop=True, inplace=True)
		print "   " + str(results.shape[0]) + " groups with plottable p-values"

		if results.shape[0] > 1:
			results['logp'] = -1 * np.log10(results[pcol]) + 0.0
			results['pos'] = results.start + (results.end - results.start) / 2

			ro.globalenv['results'] = results
			l = np.median(scipy.chi2.ppf([1-x for x in results[pcol].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
			# in R: median(qchisq(results$p, df=1, lower.tail=FALSE))/qchisq(0.5,1)
			print "   genomic inflation (all groups) = " + str(l)

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
				dftext_label = 'lambda %~~% ' + str(l)
				ro.globalenv['dftext'] = ro.DataFrame({'x': ro.r('Inf'), 'y': 0.5, 'lab': dftext_label})

				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq.tiff')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.qq.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.pdf')
				ro.r("""
					gp<-ggplot(df)
					pp<-gp + 
						aes_string(x='a',y='b') +
						geom_ribbon(aes_string(x='a',ymin='ci_lower',ymax='ci_upper'), data=df, alpha=0.25, fill='black') + 
						geom_point(size=2) +
						geom_abline(intercept=0, slope=1, alpha=0.5) + 
						scale_x_continuous(expression(Expected~~-log[10](italic(p)))) +
						scale_y_continuous(expression(Observed~~-log[10](italic(p)))) +
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
					elif cfg['ext'] == 'eps':
						ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.qq.cropped.eps')
					else:
						ggsave = 'ggsave(filename="%s",plot=pp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq.cropped.pdf')
					ro.r("""
						gp<-ggplot(df)
						pp<-gp + 
							aes_string(x='a',y='b') +
							geom_ribbon(aes_string(x='a',ymin='ci_lower',ymax='ci_upper'), data=df, alpha=0.25, fill='black') + 
							geom_point(aes(shape=factor(shape)),size=2) +
							geom_abline(intercept=0, slope=1, alpha=0.5) + 
							scale_x_continuous(expression(Expected~~-log[10](italic(p)))) +
							scale_y_continuous(expression(Observed~~-log[10](italic(p)))) +
							theme_bw(base_size = 12) + 
							geom_text(aes_string(x='x', y='y', label='lab'), data = dftext, colour="black", vjust=0, hjust=1, size = 4, parse=TRUE) +
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.position = 'none', 
								panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
								panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
						%s
						""" % (ggsave))

			if cfg['qq_strat']:
				print "   generating frequency stratified qq plot"

				results['CMAC'] = 'E'
				results.loc[results[cfg['cmaccol']] >= 100,'CMAC'] = 'D'
				results.loc[(results[cfg['cmaccol']] >= 50) & (results[cfg['cmaccol']] < 100),'CMAC'] = 'C'
				results.loc[(results[cfg['cmaccol']] >= 20) & (results[cfg['cmaccol']] < 50),'CMAC'] = 'B'
				results.loc[(results[cfg['cmaccol']] >= 10) & (results[cfg['cmaccol']] < 20),'CMAC'] = 'A'
				lA='NA'
				lB='NA'
				lC='NA'
				lD='NA'
				lE='NA'
				lE_n=len(results[pcol][results[cfg['cmaccol']] < 10])
				lD_n=len(results[pcol][(results[cfg['cmaccol']] >= 10) & (results[cfg['cmaccol']] < 20)])
				lC_n=len(results[pcol][(results[cfg['cmaccol']] >= 20) & (results[cfg['cmaccol']] < 50)])
				lB_n=len(results[pcol][(results[cfg['cmaccol']] >= 50) & (results[cfg['cmaccol']] < 100)])
				lA_n=len(results[pcol][results[cfg['cmaccol']] >= 100])
				if lE_n > 0:
					lE=np.median(scipy.chi2.ppf([1-x for x in results[pcol][results[cfg['cmaccol']] < 10].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				if lD_n > 0:
					lD=np.median(scipy.chi2.ppf([1-x for x in results[pcol][(results[cfg['cmaccol']] >= 10) & (results[cfg['cmaccol']] < 20)].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				if lC_n > 0:
					lC=np.median(scipy.chi2.ppf([1-x for x in results[pcol][(results[cfg['cmaccol']] >= 20) & (results[cfg['cmaccol']] < 50)].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				if lB_n > 0:
					lB=np.median(scipy.chi2.ppf([1-x for x in results[pcol][(results[cfg['cmaccol']] >= 50) & (results[cfg['cmaccol']] < 100)].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				if lA_n > 0:
					lA=np.median(scipy.chi2.ppf([1-x for x in results[pcol][results[cfg['cmaccol']] >= 100].tolist()], df=1))/scipy.chi2.ppf(0.5,1)
				print "   genomic inflation (CMAC > 100, n=" + str(lA_n) + ") = " + str(lA)
				print "   genomic inflation (50 <= CMAC < 100, n=" + str(lB_n) + ") = " + str(lB)
				print "   genomic inflation (20 <= CMAC < 50, n=" + str(lC_n) + ") = " + str(lC)
				print "   genomic inflation (10 <= CMAC < 20, n=" + str(lD_n) + ") = " + str(lD)
				print "   genomic inflation (CMAC < 10, n=" + str(lE_n) + ") = " + str(lE)

				a = np.array([])
				b = np.array([])
				c = np.array([])
				results.sort_values(by=['logp'], inplace=True)
				if len(results[results['CMAC'] == 'E'].index) > 0:
					aa = -1 * np.log10(ro.r('ppoints(' + str(len(results[results['CMAC'] == 'E'].index)) + ')'))
					aa.sort()
					bb = results['logp'][results['CMAC'] == 'E']
					#bb.sort()
					cc = results['CMAC'][results['CMAC'] == 'E']
					a = np.append(a,aa)
					b = np.append(b,bb)
					c = np.append(c,cc)
					print "   minimum p-value (CMAC < 10): " + str(np.min(results[pcol][results['CMAC'] == 'E']))
					print "   maximum -1*log10(p-value) (CMAC < 10): " + str(np.max(results['logp'][results['CMAC'] == 'E']))
				if len(results[results['CMAC'] == 'D'].index) > 0:
					aa = -1 * np.log10(ro.r('ppoints(' + str(len(results[results['CMAC'] == 'D'].index)) + ')'))
					aa.sort()
					bb = results['logp'][results['CMAC'] == 'D']
					#bb.sort()
					cc = results['CMAC'][results['CMAC'] == 'D']
					a = np.append(a,aa)
					b = np.append(b,bb)
					c = np.append(c,cc)
					print "   minimum p-value (10 <= CMAC < 20): " + str(np.min(results[pcol][results['CMAC'] == 'D']))
					print "   maximum -1*log10(p-value) (10 <= CMAC < 20): " + str(np.max(results['logp'][results['CMAC'] == 'D']))
				if len(results[results['CMAC'] == 'C'].index) > 0:
					aa = -1 * np.log10(ro.r('ppoints(' + str(len(results[results['CMAC'] == 'C'].index)) + ')'))
					aa.sort()
					bb = results['logp'][results['CMAC'] == 'C']
					#bb.sort()
					cc = results['CMAC'][results['CMAC'] == 'C']
					a = np.append(a,aa)
					b = np.append(b,bb)
					c = np.append(c,cc)
					print "   minimum p-value (20 <= CMAC < 50): " + str(np.min(results[pcol][results['CMAC'] == 'C']))
					print "   maximum -1*log10(p-value) (20 <= CMAC < 50): " + str(np.max(results['logp'][results['CMAC'] == 'C']))
				if len(results[results['CMAC'] == 'B'].index) > 0:
					aa = -1 * np.log10(ro.r('ppoints(' + str(len(results[results['CMAC'] == 'B'].index)) + ')'))
					aa.sort()
					bb = results['logp'][results['CMAC'] == 'B']
					#bb.sort()
					cc = results['CMAC'][results['CMAC'] == 'B']
					a = np.append(a,aa)
					b = np.append(b,bb)
					c = np.append(c,cc)
					print "   minimum p-value (50 <= CMAC < 100): " + str(np.min(results[pcol][results['CMAC'] == 'B']))
					print "   maximum -1*log10(p-value) (50 <= CMAC < 100): " + str(np.max(results['logp'][results['CMAC'] == 'B']))
				if len(results[results['CMAC'] == 'A'].index) > 0:
					aa = -1 * np.log10(ro.r('ppoints(' + str(len(results[results['CMAC'] == 'A'].index)) + ')'))
					aa.sort()
					bb = results['logp'][results['CMAC'] == 'A']
					#bb.sort()
					cc = results['CMAC'][results['CMAC'] == 'A']
					a = np.append(a,aa)
					b = np.append(b,bb)
					c = np.append(c,cc)
					print "   minimum p-value (CMAC >= 100): " + str(np.min(results[pcol][results['CMAC'] == 'A']))
					print "   maximum -1*log10(p-value) (CMAC >= 100): " + str(np.max(results['logp'][results['CMAC'] == 'A']))
			
				ro.globalenv['df'] = ro.DataFrame({'a': ro.FloatVector(a), 'b': ro.FloatVector(b), 'CMAC': ro.StrVector(c)})
			
				if cfg['ext'] == 'tiff':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat.tiff')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.qq_strat.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat.pdf')
				ro.r("""
					gp<-ggplot(df, aes_string(x='a',y='b')) +
						geom_point(aes_string(color='CMAC'), size=2) +
						scale_colour_manual(values=c("E"="#a8ddb5", "D"="#7bccc4", "C"="#4eb3d3", "B"="#2b8cbe", "A"="#08589e"), labels=c("E"="CMAC < 10","D"="10 <= CMAC < 20","C"="20 <= CMAC < 50","B"="50 <= CMAC < 100","A"="CMAC >= 100")) +
						geom_abline(intercept=0, slope=1, alpha=0.5) + 
						scale_x_continuous(expression(Expected~~-log[10](italic(p)))) +
						scale_y_continuous(expression(Observed~~-log[10](italic(p)))) +
						theme_bw(base_size = 12) + 
						theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
							legend.key.height = unit(0.1,"in"), legend.text = element_text(size=5), legend.key = element_blank(), legend.justification = c(0,1), 
							legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
							panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
					%s
					""" % (ggsave))
			
				if np.max(results['logp']) > cfg['crop']:
					print "   generating cropped frequency stratified qq plot"
					ro.r('df$b[df$b > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
					ro.r('df$shape<-0')
					ro.r('df$shape[df$b == ' + str(cfg['crop']) + ']<-1')
					if cfg['ext'] == 'tiff':
						ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.qq_strat.cropped.tiff')
					elif cfg['ext'] == 'eps':
						ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.qq_strat.cropped.eps')
					else:
						ggsave = 'ggsave(filename="%s",plot=gp,width=4,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.qq_strat.cropped.pdf')
					ro.r("""
						gp<-ggplot(df, aes_string(x='a',y='b')) +
							geom_point(aes(shape=factor(shape), color=CMAC), size=2) +
							scale_colour_manual(values=c("E"="#a8ddb5", "D"="#7bccc4", "C"="#4eb3d3", "B"="#2b8cbe", "A"="#08589e"), labels=c("E"="CMAC < 10","D"="10 <= CMAC < 20","C"="20 <= CMAC < 50","B"="50 <= CMAC < 100","A"="CMAC >= 100")) +
							geom_abline(intercept=0, slope=1, alpha=0.5) + 
							scale_x_continuous(expression(Expected~~-log[10](italic(p)))) +
							scale_y_continuous(expression(Observed~~-log[10](italic(p)))) +
							theme_bw(base_size = 12) + 
							guides(shape=FALSE) + 
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), legend.title = element_blank(), 
								legend.key.height = unit(0.1,"in"), legend.text = element_text(size=5), legend.key = element_blank(), legend.justification = c(0,1), 
								legend.position = c(0,1), panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
								panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.text = element_text(size=12))
						%s
						""" % (ggsave))

			if cfg['mht']:
				print "   generating standard manhattan plot"
				print "   minimum p-value: " + str(np.min(results[pcol]))
				print "   maximum -1*log10(p-value): " + str(np.max(results['logp']))
				if cfg['gc']:
					print "   adjusting p-values for genomic inflation for p-value column " + pcol
					results[pcol]=2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*results[pcol]) / math.sqrt(l)))
					print "   minimum post-gc adjustment p-value: " + str(np.min(results[pcol]))
					print "   maximum post-gc adjustment -1*log10(p-value): " + str(np.max(results['logp']))
				else:
					print "   skipping genomic inflation correction"

				print "   calculating genomic positions"
				results.sort_values(by=['#chr','pos'], inplace=True)
				ticks = []
				lastbase = 0
				results['gpos'] = 0
				nchr = len(list(np.unique(results['#chr'].values)))
				chrs = np.unique(results['#chr'].values)
				if cfg['color']:
					colours = ["#08306B","#41AB5D","#000000","#F16913","#3F007D","#EF3B2C","#08519C","#238B45","#252525","#D94801","#54278F","#CB181D","#2171B5","#006D2C","#525252","#A63603","#6A51A3","#A50F15","#4292C6","#00441B","#737373","#7F2704","#807DBA","#67000D"]
				else:
					colours = ["#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3","#08589e","#4eb3d3"]
				if nchr == 1:
					results['gpos'] = results['pos']
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
							results.loc[results['#chr'] == chrs[i],'gpos'] = results.loc[results['#chr'] == chrs[i],'pos']
						else:
							if results.loc[results['#chr'] == chrs[i-1]].shape[0] > 1:
								lastbase = lastbase + results.loc[results['#chr'] == chrs[i-1],'pos'].iloc[-1]
							else:
								lastbase = lastbase + results.loc[results['#chr'] == chrs[i-1],'pos'].iloc[0]
							results.loc[results['#chr'] == chrs[i],'gpos'] = (results.loc[results['#chr'] == chrs[i],'pos']) + lastbase
						if results.loc[results['#chr'] == chrs[i]].shape[0] > 1:
							ticks.append(results.loc[results['#chr'] == chrs[i],'gpos'].iloc[0] + (results.loc[results['#chr'] == chrs[i],'gpos'].iloc[-1] - results.loc[results['#chr'] == chrs[i],'gpos'].iloc[0])/2)
						else:
							ticks.append(results.loc[results['#chr'] == chrs[i],'gpos'].iloc[0])
						results.loc[results['#chr'] == chrs[i],'colours'] = colours[int(chrs[i])]

				results['logp'] = -1 * np.log10(results[pcol])
				if results.shape[0] >= 1000000:
					sig = 5.4e-8
				else:
					sig = 0.05 / results.shape[0]
				print "   significance level set to p-value = " + str(sig) + " (-1*log10(p-value) = " + str(-1 * np.log10(sig)) + ")"
				chr = results['#chr'][0]
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
					ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.tiff')
				elif cfg['ext'] == 'eps':
					ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.mht.eps')
				else:
					ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.pdf')
				if nchr == 1:
					ro.r("""
						gp<-ggplot(df, aes_string(x='gpos',y='logp')) +
							geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
							geom_point(size=1.5) + 
							scale_x_continuous(expression(Chromosome~~%d~~(kb))'),breaks=ticks,labels=labels) + \
							scale_y_continuous(expression(-log[10](italic(p))),limits=c(0,%d)) + \
							theme_bw(base_size = 8) + \
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
									panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
									panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=10), 
									axis.text = element_text(size=8), legend.position = 'none', axis.text = element_text(size=12))
						%s
						""" % (sig, chr, maxy, ggsave))
				else:
					ro.r("""
						gp = ggplot(df, aes_string(x='gpos',y='logp',colour='colours')) + 
							geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
							geom_point(size=1.5) + 
							scale_colour_manual(values=colours) + 
							scale_x_continuous(expression(Chromosome),breaks=ticks,labels=chrs) + 
							scale_y_continuous(expression(-log[10](italic(p))),limits=c(0,%d)) + 
							theme_bw(base_size = 8) + 
							theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
									panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
									panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=8), 
									axis.text = element_text(size=6), legend.position = 'none', axis.text = element_text(size=12))
						%s
						""" % (sig, maxy, ggsave))

				if maxy > cfg['crop']:
					maxy = cfg['crop']
					ro.r('df$logp[df$logp > ' + str(cfg['crop']) + ']<-' + str(cfg['crop']))
					ro.r('df$shape<-0')
					ro.r('df$shape[df$logp == ' + str(cfg['crop']) + ']<-1')
					print "   generating cropped manhattan plot"
					if cfg['ext'] == 'tiff':
						ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,units="in",bg="white",compression="lzw",dpi=300)' % (cfg['out'] + '.' + pcol + '.mht.cropped.tiff')
					elif cfg['ext'] == 'eps':
						ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,bg="white",horizontal=True)' % (cfg['out'] + '.' + pcol + '.mht.cropped.eps')
					else:
						ggsave = 'ggsave(filename="%s",plot=gp,width=12,height=4,bg="white")' % (cfg['out'] + '.' + pcol + '.mht.cropped.pdf')
					if nchr == 1:
						ro.r("""
							gp<-ggplot(df, aes_string(x='gpos',y='logp')) +
								geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
								geom_point(aes(shape=factor(shape)),size=1.5) + 
								scale_x_continuous(expression(Chromosome~~%d~~(kb))'),breaks=ticks,labels=labels) + 
								scale_y_continuous(expression(-log[10](italic(p))),limits=c(0,%d)) + 
								theme_bw(base_size = 8) + 
								theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
										panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
										panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=10), 
										axis.text = element_text(size=8), legend.position = 'none', axis.text = element_text(size=12))
							%s
							""" % (sig, chr, maxy, ggsave))
					else:
						ro.r("""
							gp = ggplot(df, aes_string(x='gpos',y='logp',colour='colours')) + 
								geom_hline(yintercept = -1 * log10(%g),colour="#B8860B", linetype=5, size = 0.25) + 
								geom_point(aes(shape=factor(shape)),size=1.5) + 
								scale_colour_manual(values=colours) + 
								scale_x_continuous(expression(Chromosome),breaks=ticks,labels=chrs) + 
								scale_y_continuous(expression(-log[10](italic(p))),limits=c(0,%d)) + 
								theme_bw(base_size = 8) + 
								theme(axis.title.x = element_text(vjust=-0.5,size=14), axis.title.y = element_text(vjust=1,angle=90,size=14), 
										panel.background = element_blank(), panel.border = element_blank(), panel.grid.minor = element_blank(), 
										panel.grid.major = element_blank(), axis.line = element_line(colour="black"), axis.title = element_text(size=8), 
										axis.text = element_text(size=6), legend.position = 'none', axis.text = element_text(size=12))
							%s
							""" % (sig, maxy, ggsave))

	print "process complete"
	return 0
