from __main__ import *
import MiscFxns
import scipy.stats as scipy
from itertools import islice

def Meta(cfg):
	Parse.PrintMetaOptions(cfg)

	assert cfg, SystemFxns.Error("no configuration file specified")

	##### get headers from process files #####
	for tag in cfg['data_info'].keys():
		try:
			p = subprocess.Popen(['tabix','-h',cfg['data_info'][tag]['file'],'0'], stdout=subprocess.PIPE)
		except:
			print SystemFxns.Error("file " + cfg['data_info'][tag]['file'] + " has incorrect format")
			return
		cfg['data_info'][tag]['header'] = p.communicate()[0]
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].replace("#","")
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].strip()
		cfg['data_info'][tag]['header'] = cfg['data_info'][tag]['header'].split()		

	##### initialize out file #####
	bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')
	
	output = {}
	out_cols = []
	for tag in cfg['file_order']:
		for x in ['freq','effect','stderr','or','z','p']:
			if x in cfg['data_info'][tag]:
				out_cols.append(tag + '.' + x)
	out_cols.append(tag + '.filtered')
	
	##### determine markers to be analyzed #####
	if not cfg['reglist'] is None:
		print "loading region list"
		reglist = FileFxns.LoadCoordinates(cfg['reglist'])
	elif not cfg['region'] is None:
		if len(cfg['region'].split(':')) > 1:
			reglist = pd.DataFrame({'chr': [re.split(':|-',cfg['region'])[0]],'start': [re.split(':|-',cfg['region'])[1]],'end': [re.split(':|-',cfg['region'])[2]],'region': [cfg['region']]})
		else:
			reglist = pd.DataFrame({'chr': [cfg['region']],'start': ['NA'],'end': ['NA'],'region': [cfg['region']]})
		reglist['id'] = cfg['id'] if not cfg['id'] is None else 'NA'
	else:
		reglist = pd.DataFrame({'chr': [str(i+1) for i in range(26)],'start': ['NA' for i in range(26)],'end': ['NA' for i in range(26)],'region': [str(i+1) for i in range(26)],'id': ['NA' for i in range(26)]})
	reglist['n'] = 0
	for i in range(len(reglist.index)):
		for key in cfg['data_info'].keys():
			tb = tabix.open(cfg['data_info'][key]['file'])
			try:
				records = tb.querys(reglist['region'][i])
			except:
				pass
			else:
				reglist['n'][i] = reglist['n'][i] + 1
				break
	if reglist['n'].sum() == 0:
		print SystemFxns.Error("no markers found")
		return()

	##### load iterators #####
	for k in cfg['file_order']:
		print "loading data files for cohort " + k
		cfg['data_info'][k]['file_it'] = MiscFxns.LoadResults(cfg['data_info'][k]['file'])

	##### LOW MEMORY VERSION #####
	for r in range(len(reglist.index)):
		reg = reglist['region'][r]
		print "building reference database for " + str(r + 1) + " of " + str(len(reglist.index)) + " regions"
		i = 0
		for k in cfg['file_order']:
			try:
				records = cfg['data_info'][k]['file_it'].querys(reg)
			except:
				break
			while True:
				i = i + 1
				chunk=list(islice(records, None))
				if not chunk:
					break

				##### READ IN CHUNK AND CREATE STANDARD DATA FRAME #####
				print "   aligning " + k + " in region " + str(r + 1) + " of " + str(len(reglist.index)) + " (" + reg + ")",
				chunkdf = pd.DataFrame(chunk)
				chunkdf.columns = cfg['data_info'][k]['header']
				cols = [x for x in ['chr','pos','a1','a2',cfg['data_info'][k]['marker_col'],cfg['data_info'][k]['effect_col'],cfg['data_info'][k]['stderr_col'],cfg['data_info'][k]['freq_col'],
									cfg['data_info'][k]['rsq_col'],cfg['data_info'][k]['hwe_col'],cfg['data_info'][k]['or_col'],cfg['data_info'][k]['z_col'],
									cfg['data_info'][k]['p_col']] if not x is None]
				chunkdf = chunkdf[cols]
				chunkdf.rename(columns=dict([x for x in [(cfg['data_info'][k]['marker_col'],k + '.marker'),(cfg['data_info'][k]['effect_col'],k + '.effect'),(cfg['data_info'][k]['stderr_col'],k + '.stderr'),
							(cfg['data_info'][k]['freq_col'],k + '.freq'),(cfg['data_info'][k]['rsq_col'],k + '.rsq'),(cfg['data_info'][k]['hwe_col'],k + '.hwe'),
							(cfg['data_info'][k]['or_col'],k + '.or'),(cfg['data_info'][k]['z_col'],k + '.z'),(cfg['data_info'][k]['p_col'],k + '.p')] if not x[0] is None]),inplace=True)
				chunkdf = chunkdf.drop_duplicates(subset=['chr','pos','a1','a2'])
				chunkdf.index = pd.MultiIndex.from_tuples(tuple(chunkdf.apply(lambda row: MiscFxns.ListCompatibleMarkersMeta(row['chr'],row['pos'],row['a1'],row['a2'],'><'),axis=1)),names=['x1','x2','x3','x4'])
				chunkdf.ix[:,5:][chunkdf.ix[:,5:] == 'NA'] = float('nan')
				if 'n' in cfg['data_info'][k]:
					chunkdf[k + '.n'] = int(cfg['data_info'][k]['n'])
				chunkdf[k + '.filter'] = 0
				if 'rsq' in cfg['data_info'][k] and cfg['data_info'][k]['rsq_col'] is not None:
					chunkdf[k + '.filter'][(math.isnan(chunkdf[cfg['data_info'][k]['rsq_col']])) | (chunkdf[cfg['data_info'][k]['rsq_col']] < cfg['data_info'][k]['rsq'])] = 1
				if 'maf' in cfg['data_info'][k] and cfg['data_info'][k]['freq_col'] is not None:
					chunkdf[k + '.filter'][(math.isnan(chunkdf[cfg['data_info'][k]['freq_col']])) | (chunkdf[cfg['data_info'][k]['freq_col']] < cfg['data_info'][k]['maf']) | (chunkdf[cfg['data_info'][k]['freq_col']] > 1 - cfg['data_info'][k]['maf'])] = 1
				if 'hwe' in cfg['data_info'][k] and cfg['data_info'][k]['hwe_col'] is not None:
					chunkdf[k + '.filter'][(math.isnan(chunkdf[cfg['data_info'][k]['hwe_col']])) | (chunkdf[cfg['data_info'][k]['hwe_col']] < cfg['data_info'][k]['hwe'])] = 1
				if i == 1:
					output_df = chunkdf.copy()
				else:
					chunkdf.rename(columns={'chr': k + '.chr', 'pos': k + '.pos', 'a1': k + '.a1', 'a2': k + '.a2'}, inplace=True)
					output_df = output_df.join(chunkdf,how='outer')
					output_df['a1'][output_df['a1'].isnull()] = output_df[k + '.a1'][output_df['a1'].isnull()]
					output_df['a2'][output_df['a2'].isnull()] = output_df[k + '.a2'][output_df['a2'].isnull()]
					output_df['chr'][output_df['chr'].isnull()] = output_df[k + '.chr'][output_df['chr'].isnull()]
					output_df['pos'][output_df['pos'].isnull()] = output_df[k + '.pos'][output_df['pos'].isnull()]
					if k + '.effect' in output_df:
						output_df[k + '.effect'][~output_df[k + '.effect'].isnull()] = output_df[~output_df[k + '.effect'].isnull()].apply(lambda row: MiscFxns.FlipEffect(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.effect']),1)
					if k + '.freq' in output_df:
						output_df[k + '.freq'][~output_df[k + '.freq'].isnull()] = output_df[~output_df[k + '.freq'].isnull()].apply(lambda row: MiscFxns.FlipFreq(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.freq']),1)
					if k + '.or' in output_df:
						output_df[k + '.or'][~output_df[k + '.or'].isnull()] = output_df[~output_df[k + '.or'].isnull()].apply(lambda row: MiscFxns.FlipOR(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.or']),1)
					if k + '.z' in output_df:
						output_df[k + '.z'][~output_df[k + '.z'].isnull()] = output_df[~output_df[k + '.z'].isnull()].apply(lambda row: MiscFxns.FlipZ(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.z']),1)
					output_df.drop(labels=[k + '.chr',k + '.pos',k + '.a1',k + '.a2'],axis=1,inplace=True)
				print ": " + str(output_df.shape[0]) + " markers"

		print output_df.head(); print output_df.shape; return

		##### DO META FOR REGION AND WRITE TO FILE #####
		##### LOOP TO NEXT REGION #####


	header=output_df.columns.values
	output_df['chr'] = output_df['chr'].astype(int)
	output_df['pos'] = output_df['pos'].astype(int)
	for tag in cfg['file_order']:
		for suffix in ['freq','effect','stderr','rsq','hwe','or','z','p','n','filtered']:
			if tag + '.' + suffix in output_df.columns.values:
				output_df[tag + '.' + suffix] = output_df[tag + '.' + suffix].convert_objects(convert_numeric=True)
	print header; print output_df.dtypes; print output_df.shape; sys.exit()

	##### apply genomic control #####
	for tag in cfg['file_order']:
		if 'gc' in cfg['data_info'][tag].keys():
			print "applying genomic control correction for cohort " + tag
			if 'stderr' in cfg['data_info'][tag].keys():
				output_df[tag + '.stderr'] = output_df[tag + '.stderr'] * math.sqrt(float(cfg['data_info'][tag]['gc']))
			if 'z' in cfg['data_info'][tag].keys():
				output_df[tag + '.z'] = output_df[tag + '.z'] / math.sqrt(float(cfg['data_info'][tag]['gc']))
			if 'p' in cfg['data_info'][tag].keys():
				if 'stderr' in cfg['data_info'][tag].keys():
					output_df[tag + '.p'] = 2 * scipy.norm.cdf(-1 * np.abs(output_df[tag + '.effect']) / output_df[tag + '.stderr'])
				else:
					output_df[tag + '.p'] = 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*output_df[tag + '.p']) / math.sqrt(float(cfg['data_info'][tag]['gc']))))

	##### meta analysis #####
	for meta in cfg['meta_order']:
		print "   running meta analysis for cohorts " + meta
		if cfg['method'] in ['efftest','sample_size']:
			p_ext = '.p' if cfg['method'] == 'sample_size' else '.p_eff'
			output_df[meta + '.dir'] = ''
			for tag in cfg['meta_info'][meta]:
				output_df[tag + '.filtered'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + p_ext]) or x[tag + p_ext] > 1 or x[tag + p_ext] <= 0 else x[tag + '.filtered'],axis=1)
				filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filtered')][0]
				N_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith(p_ext)][0]
				Eff_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				output_df[meta + '.' + tag + '.dir'] = output_df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				output_df[meta + '.' + tag + '.zi'] = output_df.apply(lambda x: (-1 * scipy.norm.ppf(1 - (x[P_idx]/2)) if x[Eff_idx] < 0 else scipy.norm.ppf(1 - (x[P_idx]/2))) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.n'] = output_df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.wi'] = output_df.apply(lambda x: math.sqrt(x[meta + '.' + tag + '.n']) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.ziwi'] = output_df.apply(lambda x: x[meta + '.' + tag + '.zi'] * x[meta + '.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.dir'] = output_df[meta + '.dir'] + output_df[meta + '.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.wi')]
			ZiWi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.ziwi')]
			output_df[meta + '.n'] = output_df.apply(lambda x: x[N_idx_all].sum(),axis=1)
			output_df[meta + '.z'] = output_df.apply(lambda x: x[ZiWi_idx_all].sum()/math.sqrt(x[N_idx_all].sum()) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			output_df[meta + '.stderr'] = output_df.apply(lambda x: float('nan'), axis=1)
			output_df[meta + '.effect'] = output_df.apply(lambda x: float('nan'), axis=1)
		else:
			output_df[meta + '.dir'] = ''
			for tag in cfg['meta_info'][meta]:
				output_df[tag + '.filtered'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + '.p']) or x[tag + '.p'] > 1 or x[tag + '.p'] <= 0 else x[tag + '.filtered'],axis=1)
				filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filtered')][0]
				N_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
				P_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.p')][0]
				Eff_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.effect')][0]
				StdErr_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.stderr')][0]
				output_df[meta + '.' + tag + '.dir'] = output_df.apply(lambda x: ('-' if x[Eff_idx] < 0 else '+') if x[filter_idx] == 0 and x[P_idx] <= 1 else 'x',axis=1)
				output_df[meta + '.' + tag + '.n'] = output_df.apply(lambda x: x[tag + '.n'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.wi'] = output_df.apply(lambda x: 1/(x[StdErr_idx]**2) if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.' + tag + '.biwi'] = output_df.apply(lambda x: x[Eff_idx] * x[meta + '.' + tag + '.wi'] if x[filter_idx] == 0 and x[P_idx] <= 1 else float('nan'),axis=1)
				output_df[meta + '.dir'] = output_df[meta + '.dir'] + output_df[meta + '.' + tag + '.dir']
			N_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.n')]
			Wi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.wi')]
			BiWi_idx_all=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(meta) and s.endswith('.biwi')]
			output_df[meta + '.n'] = output_df.apply(lambda x: x[N_idx_all].sum(),axis=1)
			output_df[meta + '.stderr'] = output_df.apply(lambda x: math.sqrt(1/(x[Wi_idx_all].sum())) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			output_df[meta + '.effect'] = output_df.apply(lambda x: (x[BiWi_idx_all].sum())/(x[Wi_idx_all].sum()) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			output_df[meta + '.z'] = output_df.apply(lambda x: x[meta + '.effect']/x[meta + '.stderr'] if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
			header.append(meta + '.effect')
			header.append(meta + '.stderr')
		output_df[meta + '.p'] = output_df.apply(lambda x: 2 * scipy.norm.cdf(-1 * abs(float(x[meta + '.z']))) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
		output_df[meta + '.dir'] = output_df.apply(lambda x: x[meta + '.dir'] if not math.isnan(x[meta + '.p']) else float('nan'), axis=1)
		output_df[meta + '.effect'] = output_df[meta + '.effect'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.stderr'] = output_df[meta + '.stderr'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.z'] = output_df[meta + '.z'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		output_df[meta + '.p'] = output_df[meta + '.p'].map(lambda x: '%.4e' % x if not math.isnan(x) else x)
		header.append(meta + '.z')
		header.append(meta + '.p')
		header.append(meta + '.dir')
		header.append(meta + '.n')
		
	##### format output to significant digits #####
	for tag in cfg['file_order']:
		if 'effect' in cfg['data_info'][tag].keys():
			output_df[tag + '.effect'] = output_df[tag + '.effect'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'stderr' in cfg['data_info'][tag].keys():
			output_df[tag + '.stderr'] = output_df[tag + '.stderr'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'or' in cfg['data_info'][tag].keys():
			output_df[tag + '.or'] = output_df[tag + '.or'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'z' in cfg['data_info'][tag].keys():
			output_df[tag + '.z'] = output_df[tag + '.z'].map(lambda x: '%.5g' % x if not math.isnan(x) else x)
		if 'p' in cfg['data_info'][tag].keys():
			output_df[tag + '.p'] = output_df[tag + '.p'].map(lambda x: '%.4e' % x if not math.isnan(x) else x)
		if 'n' in cfg['data_info'][tag].keys():
			output_df[tag + '.n'] = output_df[tag + '.n'].map(lambda x: '%d' % x if not math.isnan(x) else x)
		output_df[tag + '.filtered'] = output_df[tag + '.filtered'].map(lambda x: '%d' % x if not math.isnan(x) else x)
	for meta in cfg['meta_order']:
		if meta + '.n' in list(output_df.columns.values):
			output_df[meta + '.n'] = output_df[meta + '.n'].map(lambda x: '%d' % x if not math.isnan(x) else x)

	##### extract efftest corrected top p values #####
	if cfg['method'] == 'efftest':
		out_efftest_df = reglist.copy()
		header = ['chr','start','end','id']
		for tag in cfg['file_order']:
			out_efftest_df[tag + '.min_snp'] = float('nan')
			out_efftest_df[tag + '.min_pos'] = float('nan')
			out_efftest_df[tag + '.min_p'] = float('nan')
			for r in range(len(out_efftest_df.index)):
				min_idx = output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)][tag + '.p_eff'].argmin() if len(output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)].index) > 0 else float('nan')
				out_efftest_df[tag + '.min_snp'][r] = output_df['marker'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[tag + '.min_pos'][r] = output_df['pos'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[tag + '.min_p'][r] = output_df[tag + '.p'][min_idx] if min_idx in output_df.index else float('nan')
			header.append(tag + '.min_snp')
			header.append(tag + '.min_pos')
			header.append(tag + '.min_p')
			out_efftest_df[tag + '.min_pos'] = out_efftest_df[tag + '.min_pos'].map(lambda x: '%d' % x if not math.isnan(x) else x)
			out_efftest_df[tag + '.min_p'] = out_efftest_df[tag + '.min_p'].map(lambda x: '%.4e' % x if not math.isnan(x) else x)
		for meta in cfg['meta_order']:
			out_efftest_df[meta + '.min_snp'] = float('nan')
			out_efftest_df[meta + '.min_pos'] = float('nan')
			out_efftest_df[meta + '.min_dir'] = float('nan')
			out_efftest_df[meta + '.min_p'] = float('nan')
			for r in range(len(out_efftest_df.index)):
				min_idx = output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)][meta + '.p'].argmin() if len(output_df[(output_df['chr'] == out_efftest_df['chr'][r]) & (output_df['pos'] >= out_efftest_df['start'][r]) & (output_df['pos'] <= out_efftest_df['end'][r]) & (output_df[tag + '.filtered'] != 1)].index) > 0 else float('nan')
				out_efftest_df[meta + '.min_snp'][r] = output_df['marker'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_pos'][r] = output_df['pos'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_dir'][r] = output_df[meta + '.dir'][min_idx] if min_idx in output_df.index else float('nan')
				out_efftest_df[meta + '.min_p'][r] = output_df[meta + '.p'][min_idx] if min_idx in output_df.index else float('nan')
			header.append(meta + '.min_snp')
			header.append(meta + '.min_pos')
			header.append(meta + '.min_dir')
			header.append(meta + '.min_p')
			out_efftest_df[meta + '.min_pos'] = out_efftest_df[meta + '.min_pos'].map(lambda x: '%d' % x if not math.isnan(x) else x)
			out_efftest_df[meta + '.min_p'] = out_efftest_df[meta + '.min_p'].map(lambda x: '%.4e' % x if not math.isnan(x) else x)
		out_efftest_df = out_efftest_df[header]
		out_efftest_df.fillna('NA',inplace=True)
		out_efftest_df.rename(columns=lambda x: x.replace('chr','#chr'),inplace=True)
		out_efftest_df.to_csv(bgzfile, header=True, index=False, sep="\t")
	else:
		output_df = output_df[header]	
		output_df.fillna('NA',inplace=True)
		output_df.sort(['chr','pos'],inplace=True)
		output_df.rename(columns=lambda x: x.replace('chr','#chr'),inplace=True)
		output_df.to_csv(bgzfile, header=True, index=False, sep="\t")
	bgzfile.close()

	print "mapping results file"
	cmd = ['tabix','-b','2','-e','3',cfg['out'] + '.gz'] if cfg['method'] == 'efftest' else ['tabix','-b','2','-e','2',cfg['out'] + '.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessError:
		print SystemFxns.Error("file mapping failed")
	else:
		print 'process complete'
