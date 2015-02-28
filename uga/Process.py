import psutil
import subprocess
import sys
from Messages import Highlight

def Qsub(command):
	try:
		p = subprocess.Popen(command,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print Highlight("terminated by user")
		sys.exit(1)

def Interactive(submit, cmd, log_file = None):
	try:
		p = subprocess.Popen([submit,'--internal','--cmd', cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
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
		print Highlight("terminated by user")
		sys.exit(1)

def kill_all(pid):
	parent = psutil.Process(pid)
	for child in parent.children(recursive=True):
		child.kill()
	parent.kill()
