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
import MiscFxns
import MiscFxnsCy
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

	##### GENERATE REGION LIST #####
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
		print "loading data files for " + k
		cfg['data_info'][k]['file_it'] = MiscFxns.LoadResults(cfg['data_info'][k]['file'])

	##### INITIALIZE OUT FILE #####
	bgzfile = bgzf.BgzfWriter(cfg['out'] + '.gz', 'wb')

	##### ITERATE OVER REGIONS AND ALIGN VARIANTS #####
	written = False
	header = ['chr','pos','a1','a2','marker']
	for k in cfg['meta_order']:
		header = header + [k + '.z',k + '.p',k + '.dir',k + '.n']
	for k in cfg['file_order']:
		header = header + [x[1] for x in [(cfg['data_info'][k]['marker_col'],k + '.marker'),(cfg['data_info'][k]['freq_col'],k + '.freq'),
									(cfg['data_info'][k]['rsq_col'],k + '.rsq'),(cfg['data_info'][k]['hwe_col'],k + '.hwe'),
									(cfg['data_info'][k]['effect_col'],k + '.effect'),(cfg['data_info'][k]['stderr_col'],k + '.stderr'),
									(cfg['data_info'][k]['or_col'],k + '.or'),(cfg['data_info'][k]['z_col'],k + '.z'),(cfg['data_info'][k]['p_col'],k + '.p')] if not x[0] is None]
		header = header + [k + '.n',k + '.filter']

	for r in range(len(reglist.index)):
		output_df = pd.DataFrame({})
		reg = reglist['region'][r]
		i = 0
		i_all = 0
		for k in cfg['file_order']:
			i_all += 1
			try:
				records = cfg['data_info'][k]['file_it'].querys(reg)
			except:
				break
			while True:
				i += 1
				chunk=list(islice(records, None))
				if not chunk:
					break
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
				chunkdf.index = chunkdf.apply(lambda row: MiscFxnsCy.ListCompatibleMarkersMetaCy(row['chr'],row['pos'],row['a1'],row['a2'],'><'),axis=1)
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
					output_df = chunkdf
				else:
					chunkdf.rename(columns={'chr': k + '.chr', 'pos': k + '.pos', 'a1': k + '.a1', 'a2': k + '.a2'}, inplace=True)
					output_df = output_df.join(chunkdf,how='outer')
					output_df['a1'][output_df['a1'].isnull()] = output_df[k + '.a1'][output_df['a1'].isnull()]
					output_df['a2'][output_df['a2'].isnull()] = output_df[k + '.a2'][output_df['a2'].isnull()]
					output_df['chr'][output_df['chr'].isnull()] = output_df[k + '.chr'][output_df['chr'].isnull()]
					output_df['pos'][output_df['pos'].isnull()] = output_df[k + '.pos'][output_df['pos'].isnull()]
					output_df[[x for x in output_df.columns if not x in ['chr','pos','a1','a2']]] = output_df[[x for x in output_df.columns if not x in ['chr','pos','a1','a2']]].convert_objects(convert_numeric=True, copy=False)
					if k + '.effect' in output_df:
						output_df[k + '.effect'][~output_df[k + '.effect'].isnull()] = output_df[~output_df[k + '.effect'].isnull()].apply(lambda row: MiscFxnsCy.FlipEffectCy(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.effect']),1)
					if k + '.freq' in output_df:
						output_df[k + '.freq'][~output_df[k + '.freq'].isnull()] = output_df[~output_df[k + '.freq'].isnull()].apply(lambda row: MiscFxnsCy.FlipFreqCy(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.freq']),1)
					if k + '.or' in output_df:
						output_df[k + '.or'][~output_df[k + '.or'].isnull()] = output_df[~output_df[k + '.or'].isnull()].apply(lambda row: MiscFxnsCy.FlipORCy(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.or']),1)
					if k + '.z' in output_df:
						output_df[k + '.z'][~output_df[k + '.z'].isnull()] = output_df[~output_df[k + '.z'].isnull()].apply(lambda row: MiscFxnsCy.FlipZCy(row['a1'], row['a2'], row[k + '.a1'], row[k + '.a2'], row[k + '.z']),1)
					output_df.drop(labels=[k + '.chr',k + '.pos',k + '.a1',k + '.a2'],axis=1,inplace=True)
				print "region " + str(r + 1) + "/" + str(len(reglist.index)) + " (" + reg + "): adding variants from cohort " + str(i_all) + "/" + str(len(cfg['file_order'])) + " (" + k + "): " + str(chunkdf.shape[0])

		if output_df.shape[0] > 0:
			##### APPLY GC #####
			for k in cfg['file_order']:
				if 'gc' in cfg['data_info'][k]:
					print "region " + str(r + 1) + "/" + str(len(reglist.index)) + " (" + reg + "): applying genomic control for " + k + ": " + str(cfg['data_info'][k]['gc'])
					if 'stderr_col' in cfg['data_info'][k]:
						output_df[k + '.stderr'][~output_df[k + '.stderr'].isnull()] = output_df[~output_df[k + '.stderr'].isnull()].apply(lambda row: float(row[k + '.stderr']) * math.sqrt(float(cfg['data_info'][k]['gc'])),1)
					if 'z' in cfg['data_info'][k].keys():
						output_df[k + '.z'][~output_df[k + '.z'].isnull()] = output_df[~output_df[k + '.z'].isnull()].apply(lambda row: float(row[k + '.z']) / math.sqrt(float(cfg['data_info'][k]['gc'])),1)
					if 'p' in cfg['data_info'][k].keys():
						if 'stderr' in cfg['data_info'][k].keys():
							output_df[k + '.p'][~((output_df[k + '.effect'].isnull()) | (output_df[k + '.stderr'].isnull()))] = output_df[~((output_df[k + '.effect'].isnull()) | (output_df[k + '.stderr'].isnull()))].apply(lambda row: 2 * scipy.norm.cdf(-1 * np.abs(float(row[k + '.effect'])/float(row[k + '.stderr']))),1)
						else:
							output_df[k + '.p'][~output_df[k + '.p'].isnull()] = output_df[~output_df[k + '.p'].isnull()].apply(lambda row: 2 * scipy.norm.cdf(-1 * np.abs(scipy.norm.ppf(0.5*float(row[k + '.p'])) / math.sqrt(float(cfg['data_info'][k]['gc'])))),1)

			##### META ANALYSIS #####
			for meta in cfg['meta_order']:
				print "region " + str(r + 1) + "/" + str(len(reglist.index)) + " (" + reg + "): running " + meta
				if cfg['method'] == 'sample_size':
					output_df[meta + '.dir'] = ''
					for tag in cfg['meta_info'][meta]:
						output_df[tag + '.filter'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + '.p']) or x[tag + '.p'] > 1 or x[tag + '.p'] <= 0 else x[tag + '.filter'],axis=1)
						filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filter')][0]
						N_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.n')][0]
						P_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.p')][0]
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
						output_df[tag + '.filter'] = output_df.apply(lambda x: 1 if math.isnan(x[tag + '.p']) or x[tag + '.p'] > 1 or x[tag + '.p'] <= 0 else x[tag + '.filter'],axis=1)
						filter_idx=[i for i, s in enumerate(list(output_df.columns.values)) if s.startswith(tag) and s.endswith('.filter')][0]
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
				output_df[meta + '.p'] = output_df.apply(lambda x: 2 * scipy.norm.cdf(-1 * abs(float(x[meta + '.z']))) if len(x[meta + '.dir'].replace('x','')) > 1 else float('nan'), axis=1)
				output_df[meta + '.dir'] = output_df.apply(lambda x: x[meta + '.dir'] if not math.isnan(x[meta + '.p']) else float('nan'), axis=1)

			##### ADD GLOBAL MARKER TO OUTPUT #####
			output_df['marker'] = output_df[cfg['file_order'][0] + '.marker']
			for k in [x for x in cfg['file_order'] if not x == cfg['file_order'][0]]:
				output_df['marker'][(~output_df['marker'].str.startswith('rs')) & (output_df[k + '.marker'].str.startswith('rs'))]=output_df[k + '.marker'][(~output_df['marker'].str.startswith('rs')) & (output_df[k + '.marker'].str.startswith('rs'))]
			for k in [x for x in cfg['file_order'] if not x == cfg['file_order'][0]]:
				output_df['marker'][(output_df['marker'].isnull()) & (~output_df[k + '.marker'].isnull())]=output_df[k + '.marker'][(output_df['marker'].isnull()) & (~output_df[k + '.marker'].isnull())]

			##### SORT RESULTS AND CONVERT CHR, POS, AND START COLUMNS TO INT #####
			output_df[['chr','pos']] = output_df[['chr','pos']].astype(int)
			output_df.sort(columns=['chr','pos'],inplace=True)
			for c in [x for x in output_df.columns if x.endswith(('.p','.pmin'))]:
				output_df[c] = output_df[c].map(lambda x: '%.4e' % (x) if not math.isnan(x) else x)
				output_df[c] = output_df[c].astype(object)
			for c in [x for x in output_df.columns if x.endswith(('.rho','.cmaf','.Qmeta','.beta','.se','.cmafTotal','.cmafUsed','hwe','hwe.unrel','hwe.ctrl','hwe.case','hwe.unrel.ctrl','hwe.unrel.case','callrate','freq','freq.unrel','freq.ctrl','freq.case','freq.unrel.ctrl','freq.unrel.case','rsq','rsq.unrel','rsq.ctrl','rsq.case','rsq.unrel.ctrl','rsq.unrel.case','effect','stderr','or','z'))]:
				output_df[c] = output_df[c].map(lambda x: '%.5g' % (x) if not math.isnan(x) else x)
				output_df[c] = output_df[c].astype(object)
			for c in [x for x in output_df.columns if x.endswith(('filter','n','status'))]:
				output_df[c] = output_df[c].map(lambda x: '%d' % (x) if not math.isnan(x) else x)
				output_df[c] = output_df[c].astype(object)

			##### FILL IN NA's, ORDER HEADER, AND WRITE TO FILE #####
			output_df.fillna('NA',inplace=True)
			output_df.sort(['chr','pos'],inplace=True)
			for col in [x for x in header if not x in output_df.columns]:
				output_df[col] = float('nan')
			output_df = output_df[header]
			if not written:
				output_df.rename(columns=lambda x: x.replace('chr','#chr'),inplace=True)
				output_df.to_csv(bgzfile, header=True, index=False, sep="\t")
				bgzfile.flush()
				written = True
			else:
				output_df.to_csv(bgzfile, header=False, index=False, sep="\t")
				bgzfile.flush()
	bgzfile.close()

	print "mapping results file"
	cmd = ['tabix','-b','2','-e','2',cfg['out'] + '.gz']
	try:
		p = subprocess.check_call(cmd)
	except subprocess.CalledProcessError:
		print SystemFxns.Error("file mapping failed")
	else:
		print 'process complete'
