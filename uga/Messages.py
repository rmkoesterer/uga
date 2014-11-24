import sys

def Banner():
	print ''
	with open(sys.path[0] + '/README.rst') as banner:
		lines = (line.rstrip() for line in banner)
		for line in lines:
			if "**Version**" in line:
				vline=line.split(" ")
				print "   * BU_BiomedicalGenetics v" + " ".join(vline[2:])

def Error(m):
	return "\n   *** Error: " + m + "\n"
	
def Bold(m):
	return "\033[1m" + m + "\033[0m"
