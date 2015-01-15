import os
import sys
import subprocess
import tabix
import math
import pandas as pd
from multi_key_dict import multi_key_dict
from Bio import bgzf

pd.options.mode.chained_assignment = None
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
	
def Meta(config=None, 
			region=None, 
			region_list=None):

	print "   ... arguments"
	"""
	for arg in locals().keys():
		if not locals()[arg] in [None, False]:
			print "      {0:>{1}}".format(str(arg), len(max(locals().keys(),key=len))) + ": " + str(locals()[arg])

	assert cfg, Error("no configuration file specified")
	assert os.path.exists(cfg), Error("configuration file not found")
	
	##### get headers from process files #####
	for tag in config['data_info'].keys():
		try:
			p = subprocess.Popen(['tabix','-h',config['data_info'][tag]['process_file'],'0'], stdout=subprocess.PIPE)
		except:
			usage(Error("process_file " + config['data_info'][tag]['process_file'] + " has incorrect format"))
		config['data_info'][tag]['header'] = p.communicate()[0]
		config['data_info'][tag]['header'] = config['data_info'][tag]['header'].replace("#","")
		config['data_info'][tag]['header'] = config['data_info'][tag]['header'].strip()
		config['data_info'][tag]['header'] = config['data_info'][tag]['header'].split()		

	print "   ... starting meta analysis"
	written = False
	bgzfile = bgzf.BgzfWriter(config['out'] + '.gz', 'wb')
	
	output = {}
	out_cols = []
	for tag in config['file_order']:
		for x in ['freq','effect','stderr','or','z','p']:
			if x in config['data_info'][tag]:
				out_cols.append(tag + '.' + x)
	out_cols.append(tag + '.filtered')

	for r in range(len(config['marker_list'])):
		reg = configmarker_list['region'][r]
		print "   ... building reference database for " + str(r + 1) + " of " + str(len(marker_list.index)) + " regions"
		delim = "><"
		ref = multi_key_dict()
		ref_alleles = multi_key_dict()
		records_count = {}
		for tag in config['file_order']:
			chr = config['data_info'][tag]['header'].index(config['data_info'][tag]['chr'])
			pos = config['data_info'][tag]['header'].index(config['data_info'][tag]['pos'])
			marker = config['data_info'][tag]['header'].index(config['data_info'][tag]['marker'])
			a1 = config['data_info'][tag]['header'].index(config['data_info'][tag]['a1'])
			a2 = config['data_info'][tag]['header'].index(config['data_info'][tag]['a2'])
			tb = tabix.open(config['data_info'][tag]['process_file'])
			records_count[tag] = 0
			try:
				records = tb.querys(reg)
			except:
				pass
			else:
				for record in records:
					records_count[tag] = records_count[tag] + 1
					if ref.has_key(record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]):
						if not record[marker] in ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]]:
							ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]].append(record[marker])
							ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]] = list(set(ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]]))
					else:
						ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2], record[chr] + delim + record[pos] + delim + Complement(record[a1]) + delim + Complement(record[a2]),record[chr] + delim + record[pos] + delim + Complement(record[a2]) + delim + Complement(record[a1]),record[chr] + delim + record[pos] + delim + record[a2] + delim + record[a1]] = [record[marker]]
						ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2], record[chr] + delim + record[pos] + delim + Complement(record[a1]) + delim + Complement(record[a2]),record[chr] + delim + record[pos] + delim + Complement(record[a2]) + delim + Complement(record[a1]),record[chr] + delim + record[pos] + delim + record[a2] + delim + record[a1]] = [record[a1],record[a2]]
					chr = config['data_info'][tag]['header'].index(config['data_info'][tag]['chr'])
					pos = config['data_info'][tag]['header'].index(config['data_info'][tag]['pos'])
					marker = config['data_info'][tag]['header'].index(config['data_info'][tag]['marker'])
					a1 = config['data_info'][tag]['header'].index(config['data_info'][tag]['a1'])
					a2 = config['data_info'][tag]['header'].index(config['data_info'][tag]['a2'])
					effect = config['data_info'][tag]['header'].index(config['data_info'][tag]['effect'])
					stderr = config['data_info'][tag]['header'].index(config['data_info'][tag]['stderr'])
					freq = config['data_info'][tag]['header'].index(config['data_info'][tag]['freq'])
					oddsratio = config['data_info'][tag]['header'].index(config['data_info'][tag]['or'])
					z = config['data_info'][tag]['header'].index(config['data_info'][tag]['z'])
					p = config['data_info'][tag]['header'].index(config['data_info'][tag]['p'])
					filtered = 0
					refmarker = filter(lambda x:'rs' in x, ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]])
					if len(refmarker) > 0:
						refmarker = refmarker[0]
					else:
						refmarker = ref[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][0]
					refa1 = ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][0]
					refa2 = ref_alleles[record[chr] + delim + record[pos] + delim + record[a1] + delim + record[a2]][1]
					if (refa1 + refa2 == Complement(record[a2]) + Complement(record[a1]) or refa1 + refa2 == record[a2] + record[a1]) and refa1 + refa2 != "AT" and refa1 + refa2 != "TA" and refa1 + refa2 != "GC" and refa1 + refa2 != "CG":
						record[effect] = str(-1 * float(record[effect])) if record[effect] != "NA" else "NA"
						record[freq] = str(1-float(record[freq])) if record[freq] != "NA" else "NA"
						record[oddsratio] = str(1/float(record[oddsratio])) if record[oddsratio] != "NA" else "NA"
						record[z] = str(-1 * float(record[z])) if record[z] != "NA" else "NA"
					for f in config['data_info'][tag]['filters']:
						col = config['data_info'][tag]['header'].index(f.split()[0])
						if f.split()[0] != "status":
							f = f.replace(f.split()[0],str(record[col]))
						else:
							f = f.replace(f.split()[0],'"' + str(record[col]) + '"')
						exec 'if "' + str(record[col]) + '" == "NA" or not ' + f + ': filtered = 1' in {}
					out_vals = {tag + '.freq': record[freq], tag + '.effect': record[effect], tag + '.stderr': record[stderr], tag + '.or': record[oddsratio], tag + '.z': record[z], tag + '.p': record[p], tag + '.filtered': filtered}
					if not record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2 in output.keys():
						output[record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2] = dict(out_vals)
					else:
						output[record[chr] + delim + record[pos] + delim + refmarker + delim + refa1 + delim + refa2].update(dict(out_vals))
				print "          merged cohort " + tag + " : " + str(len(output.keys())) + " markers"
	output_df = pd.DataFrame(output).transpose()
	header=['chr','pos','a1','a2','marker']
	header.extend(output_df.columns.values)
	
	output_df['ref'] = output_df.index
	output_df['chr'] = output_df['ref'].apply(lambda x: x.split(delim)[0])
	output_df['pos'] = output_df['ref'].apply(lambda x: x.split(delim)[1])
	output_df['marker'] = output_df['ref'].apply(lambda x: x.split(delim)[2])
	output_df['a1'] = output_df['ref'].apply(lambda x: x.split(delim)[3])
	output_df['a2'] = output_df['ref'].apply(lambda x: x.split(delim)[4])
	##### apply genomic control #####
	for tag in config['file_order']:
		if 'gc' in config['data_info'][tag].keys():
			if 'effect' in config['data_info'][tag].keys():
				output_df[tag + '.effect'] = output_df[tag + '.effect'].apply(lambda x: '%.5f' % (float(x) / math.sqrt(float(config['data_info'][tag]['gc']))) if x != 'NA' else 'NA')
			if 'stderr' in config['data_info'][tag].keys():
				output_df[tag + '.stderr'] = output_df[tag + '.stderr'].apply(lambda x: '%.5f' % (float(x) / math.sqrt(float(config['data_info'][tag]['gc']))) if x != 'NA' else 'NA')
			if 'or' in config['data_info'][tag].keys():
				output_df[tag + '.or'] = output_df[tag + '.or'].apply(lambda x: '%.5f' % (math.exp(math.log(float(x)) / math.sqrt(float(config['data_info'][tag]['gc'])))) if x != 'NA' else 'NA')
			if 'z' in config['data_info'][tag].keys():
				output_df[tag + '.z'] = output_df[tag + '.z'].apply(lambda x: '%.5f' % (float(x) / math.sqrt(float(config['data_info'][tag]['gc']))) if x != 'NA' else 'NA')
			if 'p' in config['data_info'][tag].keys():
				output_df[tag + '.p'] = output_df[[tag + '.effect',tag + '.stderr']].apply(lambda x: '%.2e' % (2 * scipy.norm.cdf(-1 * abs(float(x[0]) / float(x[1])))) if x[0] != 'NA' and x[1] != 'NA' else 'NA', axis=1)
	
	#for meta in config['meta_order']
		
	#
				marker.data$zi<-NA
				marker.data$wi<-NA
				marker.data$zi_wi<-NA
				marker.data$bi_wi<-NA
				marker.data$fi_wi<-NA
				marker.data$zi[marker.data$keep == 1] <- ifelse(exp(as.numeric(marker.data$beta[marker.data$keep == 1])) > 1, qnorm(as.numeric(marker.data$p[marker.data$keep == 1])/2,lower.tail=FALSE), -qnorm(as.numeric(marker.data$p[marker.data$keep == 1])/2,lower.tail=FALSE))
				marker.data$wi[marker.data$keep == 1] <- ifelse(! is.na(marker.data$lambda[marker.data$keep == 1]) & as.numeric(marker.data$lambda[marker.data$keep == 1]) > 1, ifelse(marker.data$method[marker.data$keep == 1] == "SampleSize", sqrt(as.numeric(marker.data$samples[marker.data$keep == 1]))/as.numeric(marker.data$lambda[marker.data$keep == 1]), (1/(as.numeric(marker.data$stderr[marker.data$keep == 1]))^2)/as.numeric(marker.data$lambda[marker.data$keep == 1])), ifelse(marker.data$method[marker.data$keep == 1] == "SampleSize", sqrt(as.numeric(marker.data$samples[marker.data$keep == 1])), 1/(as.numeric(marker.data$stderr[marker.data$keep == 1]))^2))
				marker.data$zi_wi[marker.data$keep == 1]<-as.numeric(marker.data$zi[marker.data$keep == 1])*as.numeric(marker.data$wi[marker.data$keep == 1])
				marker.data$bi_wi[marker.data$keep == 1]<-as.numeric(marker.data$beta[marker.data$keep == 1])*as.numeric(marker.data$wi[marker.data$keep == 1])
				marker.data$fi_wi[marker.data$keep == 1]<-as.numeric(marker.data$freq[marker.data$keep == 1])*as.numeric(marker.data$wi[marker.data$keep == 1])
				dist.out$freq<-sum(marker.data$fi_wi[marker.data$keep == 1])/sum(marker.data$wi[marker.data$keep == 1]) 
				dist.out$samples<-sum(marker.data$samples[marker.data$keep == 1])
				dist.out$cases<-ifelse(length(marker.data$cases[marker.data$keep == 1]) == length(marker.data$cases[marker.data$keep == 1 & ! is.na(marker.data$cases)]), sum(marker.data$cases[marker.data$keep == 1]), NA)
				dist.out$ctrls<-ifelse(length(marker.data$ctrls[marker.data$keep == 1]) == length(marker.data$ctrls[marker.data$keep == 1 & ! is.na(marker.data$ctrls)]), sum(marker.data$ctrls[marker.data$keep == 1]), NA)
				dist.out$stderr<-ifelse(method == "SampleSize", NA, sqrt(1/(sum(marker.data$wi[marker.data$keep == 1]))))
				dist.out$beta<-ifelse(method == "SampleSize", NA, sum(marker.data$bi_wi[marker.data$keep == 1])/sum(marker.data$wi[marker.data$keep == 1]))
				dist.out$zscore<-ifelse(method == "SampleSize", formatC(sum(marker.data$zi_wi[marker.data$keep == 1])/sqrt(sum((marker.data$wi[marker.data$keep == 1])^2)),format="g",digits=sig.digits), formatC(dist.out$beta/dist.out$stderr,format="g",digits=sig.digits))
				dist.out$p<-ifelse(method == "SampleSize",formatC(2*pnorm(-abs(sum(marker.data$zi_wi[marker.data$keep == 1])/sqrt(sum((marker.data$wi[marker.data$keep == 1])^2)))),format="e",digits=sig.digits-1), formatC(2*pnorm(-abs(dist.out$beta/dist.out$stderr)),format="e",digits=sig.digits-1))
				dist.out$dir<-paste(marker.data$dir.keep[marker.data$tag %in% sets], collapse="")
				dist.out$hetchisq<-ifelse(method == "SampleSize", NA, sum(marker.data$wi[marker.data$keep == 1]*(marker.data$beta[marker.data$keep == 1])^2)-(((sum(marker.data$beta[marker.data$keep == 1]*marker.data$wi[marker.data$keep == 1]))^2)/sum(marker.data$wi[marker.data$keep == 1])))
				dist.out$hetdf<-ifelse(method == "SampleSize", NA, length(marker.data$keep[marker.data$keep == 1])-1)
				dist.out$hetisq<-ifelse(method == "SampleSize", NA, ifelse(dist.out$hetchisq > dist.out$hetdf, ((dist.out$hetchisq-dist.out$hetdf)/dist.out$hetchisq)*100, 0))
				dist.out$hetp<-ifelse(method == "SampleSize", NA, pchisq(dist.out$hetchisq,df=dist.out$hetdf,lower.tail=FALSE))
	#
	output_df = output_df[header]		
	output_df.fillna('NA',inplace=True)
	#write_header = True if r == 0 else False
	#output_df.to_csv(bgzfile, header=write_header, index=False, sep="\t")
	output_df.to_csv(bgzfile, header=True, index=False, sep="\t")
	bgzfile.close()
	#
		print "   ... writing final results to out file"
		bgzfile = bgzf.BgzfWriter(config['out'] + '.gz', 'wb')
		bgzfile.write("#chr" + "\t" + "chr.pos" + "\t" + "marker" + "\t" + "a1" + "\t" + "a2" + "\t" + "\t".join(out_cols) + "\n")
		for i, row in out_order_df.iterrows():
			k = str(row[0]) + delim + str(row[1]) + delim + delim.join([str(a) for a in row[2:] if not a is None])
			for s in out_cols:
				if not s in output[k].keys():
					if "filtered" in s:
						output[k][s] = 1
					else:
						output[k][s] = "NA"
			bgzfile.write(k.split(delim)[0] + "\t" + k.split(delim)[1] + "\t" + delim.join(k.split(delim)[2:(len(k.split(delim))-2)]) + "\t" + k.split(delim)[len(k.split(delim))-2] + "\t" + k.split(delim)[len(k.split(delim))-1] + "\t" + "\t".join([str(output[k][j]) for j in out_cols]) + "\n")

	bgzfile.close()
	print "   ... mapping results file"
	cmd = 'tabix -b 2 -e 2 ' + config['out'] + '.gz'
	p = subprocess.popen(cmd, shell=True)
	time.sleep(1)
	p.wait()
	"""
	
	print '   ... meta analysis complete'
