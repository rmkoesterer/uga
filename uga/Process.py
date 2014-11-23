import psutil
import subprocess
import sys

def Qsub(command):
	try:
		p = subprocess.Popen(command,shell=True)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "\n   ... process terminated\n"
		sys.exit(1)

def Interactive(command, command_arg, log_file):
	try:
		log = open(log_file, 'w')
		p = subprocess.Popen([command,command_arg], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
		for line in iter(p.stdout.readline, ''):
			sys.stdout.write(line)
			log.write(line)
		p.wait()
		log.close()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "\n   ... process terminated\n"
		sys.exit(1)

def kill_all(pid):
	parent = psutil.Process(pid)
	for child in parent.children(recursive=True):
		child.kill()
	parent.kill()
