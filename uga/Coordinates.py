import re
import pandas as pd

class Coordinates(object):
	def __init__(self, filename):
		self.filename = filename
		
	def getFilename(self):
		return self.filename
		
	def __str__(self):
		return "%s is a coordinate file" % self.filename
		
	def Load(self):
		with open(self.filename) as f:
			regions = pd.read_table(f, header=None)
		if regions.shape[1] == 2:
			regions.columns=['region','reg_id']
		elif regions.shape[1] == 1:
			regions.columns=['region']
		else:
			print Error("too many columns in region list")
			return
		regions['chr'] = regions['region'].apply(lambda row: int(row.split(':')[0]),1)
		regions['start'] = regions['region'].apply(lambda row: int(row.split(':')[1].split('-')[0]),1)
		regions['end'] = regions['region'].apply(lambda row: int(row.split(':')[1].split('-')[1]),1)
		regions['chr'] = regions['chr'].astype(int)
		regions['start'] = regions['start'].astype(int)
		regions['end'] = regions['end'].astype(int)
		regions = regions.sort_index(by=['chr','start'],ascending=[True,True])
		return regions
