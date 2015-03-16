import psutil
import subprocess
import sys

def Error(m):
	return "\n   *** Error: " + m + "\n"

def Highlight(m):
	return "\n   ... " + m + "\n"
	
def Bold(m):
	return "\033[1m" + m + "\033[0m"

def Qsub(command):
	try:
		p = subprocess.Popen(command,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print Highlight("process terminated by user")
		sys.exit(1)

def Banner():
	print ''
	with open(sys.path[0] + '/README.rst') as banner:
		lines = (line.rstrip() for line in banner)
		for line in lines:
			if "**Version**" in line:
				vline=line.split(" ")
				print "   * BU_BiomedicalGenetics v" + " ".join(vline[2:])

def Interactive(submit, cmd, log_file = None):
	try:
		p = subprocess.Popen([submit,cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
		if log_file:
			log = open(log_file, 'w')
		for line in iter(p.stdout.readline, ''):
			sys.stdout.write(line)
			if log_file:	
				log.write(line)
		p.wait()
		if log_file:
			log.close()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print Highlight("process terminated by user")
		sys.exit(1)

def kill_all(pid):
	parent = psutil.Process(pid)
	for child in parent.children(recursive=True):
		child.kill()
	parent.kill()
