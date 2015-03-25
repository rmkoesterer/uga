#!/usr/bin/env python
#$ -cwd
#$ -j y
#$ -V

import os
import pwd
import psutil
import resource
import subprocess
import sys
from time import strftime, localtime, time, gmtime

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def main(argv):

	start_time = (localtime(), time())
	env_vars = os.environ.copy()
	local=False

	user_name=pwd.getpwuid(os.getuid()).pw_name
	group_name=subprocess.check_output(['id','-ng']).rstrip()

	if not 'REQNAME' in env_vars.keys():
		local=True
		env_vars['REQNAME'] = env_vars['HOSTNAME'] + '_' + strftime('%Y_%m_%d_%H_%M_%S', start_time[0]) if 'HOSTNAME' in env_vars.keys() else strftime('%Y_%m_%d_%H_%M_%S', start_time[0])
	if not 'JOB_ID' in env_vars.keys():
		env_vars['JOB_ID'] = '1'
	if not 'SGE_TASK_ID' in env_vars.keys():
		env_vars['SGE_TASK_ID'] = 'None'

	print "start time: " + strftime("%Y-%m-%d %H:%M:%S", start_time[0])
	if 'HOSTNAME' in env_vars.keys():
		print "compute node: " + env_vars['HOSTNAME']
	if 'PWD' in env_vars.keys():
		print "current directory: " + env_vars['PWD']
	print "user name: " + user_name
	print "group id: " + group_name
	if 'JOB_ID' in env_vars.keys():
		print "job id: " + env_vars['JOB_ID']
	if 'REQNAME' in env_vars.keys():
		print "job name: " + env_vars['REQNAME']
	if 'SGE_TASK_ID' in env_vars.keys():
		print "task index number: " + env_vars['SGE_TASK_ID']

	if argv[1].split('(')[0] == "Model":
		from uga.Model import Model
	elif argv[1].split('(')[0] == "Meta":
		from uga.Meta import Meta
	elif argv[1].split('(')[0] == "Map":
		from uga.Map import Map
	elif argv[1].split('(')[0] == "Explore":
		from uga.Explore import Explore

	eval(argv[1])

	end_time = (localtime(), time())
	process = psutil.Process(os.getpid())
	mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000.0
	print 'finish time: ' + strftime("%Y-%m-%d %H:%M:%S", end_time[0])
	print 'time elapsed: ' + strftime('%H:%M:%S', gmtime(end_time[1] - start_time[1]))
	print 'memory used: ' + str('%.2f' % mem) + ' MB'

if __name__ == "__main__":
	main(sys.argv)
