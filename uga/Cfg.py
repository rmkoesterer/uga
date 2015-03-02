import sys
from Messages import Error

class Cfg(object):

	def __init__(self, filename, module, vars = None):
		self.filename = filename
		self.module = module
		self.vars = dict(item.split('=') for item in vars.split(',')) if vars else {}
		
	def getFilename(self):
		return self.filename
	
	def getModule(self):
		return self.module
		
	def getVars(self):
		return self.vars
		
	def __str__(self):
		return "%s is a configuration file for module %s with variables %s" % (self.filename, self.module, self.vars)
		
	def Load(self):
		if self.module == 'meta':
			config = {'out': None, 'sig': 5, 'method': None, 'data_info': {}, 'meta_info': {}, 'meta_order': [], 'file_order': []}
			config_temp = {'filters': []}
			with open(self.filename) as f:
				lines = (line.rstrip() for line in f)
				lines = (line for line in lines if line)
				i = 0
				for line in lines:
					for k in self.vars.keys():
						line = line.replace('[' + k + ']', self.vars[k])
					key = str(line.split()[0])
					val = " ".join(line.split()[1:])
					if key in ["out","sig","method"]:
						config[key] = val
					elif key == "remove_filters":
						config_temp['filters'] = []
					elif key == "filter":
						config_temp['filters'].append(val)
					elif key == "process_meta":
						config['meta_order'].append(val.split(':')[0])
						config['meta_info'][val.split(':')[0]] = val.split(':')[1].split('+')
					elif key == "process_file":
						i = i + 1
						if not 'tag' in config_temp.keys():
							config_temp['tag']='FILE' + str(i)
						config_temp['process_file'] = val
						config['data_info'][config_temp['tag']] = dict(config_temp)
						config['file_order'].append(config_temp['tag'])
					else:
						config_temp[key] = val
			return config
		elif self.module == 'model':
			config = {'out': None, 'sig': 5, 'buffer': 100, 'miss': None, 'freq': None, 'rsq': None, 'hwe': None, 'mem': 3, 'nofail': False, 'data_info': {}, 'data_order': []}
			config_temp = {}
			with open(self.filename) as f:
				lines = (line.rstrip() for line in f)
				lines = (line for line in lines if line)
				i = 0
				for line in lines:
					for k in self.vars.keys():
						line = line.replace('[' + k + ']', self.vars[k])
					key = str(line.split()[0])
					val = " ".join(line.split()[1:])
					if key in ["out","sig","buffer","miss","freq","rsq","hwe","mem","nofail"]:
						config[key] = val
					elif key == "process_data":
						i = i + 1
						if not 'tag' in config_temp.keys():
							config_temp['tag']='FILE' + str(i)
						config_temp['process_data'] = val
						config['data_info'][config_temp['tag']] = dict(config_temp)
						config['data_order'].append(config_temp['tag'])
					else:
						config_temp[key] = val
			return config
		else:
			print Error('module ' + model + ' cannot be used with cfg file')
