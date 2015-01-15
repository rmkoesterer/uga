#!/usr/bin/env python
#$ -cwd
#$ -j y
#$ -V

import os
import psutil
import subprocess
import argparse
import sys
from uga.Analyze import Analyze
from time import strftime, localtime, mktime
#from memory_profiler import profile, memory_usage


#@profile
def main(args=None):
	parser = argparse.ArgumentParser('submit.py')
	parser.add_argument('--qsub', 
						action='store', 
						help='a group/project id for cluster group identification')
	parser_required = parser.add_argument_group('required arguments')
	parser_required.add_argument('--cmd', 
						action='store',
						required=True, 
						help='python function command')

	args=parser.parse_args()

	start_time = localtime()
	env_vars = os.environ.copy()
	local=False
	env_vars['PROJ_ID'] = 'None' if not args.qsub else args.qsub
	if not 'REQNAME' in env_vars.keys():
		local=True
		env_vars['REQNAME'] = env_vars['HOSTNAME'] + '_' + strftime('%Y_%m_%d_%H_%M_%S', start_time)
	if not 'JOB_ID' in env_vars.keys():
		env_vars['JOB_ID'] = '1'
	if not 'SGE_TASK_ID' in env_vars.keys():
		env_vars['SGE_TASK_ID'] = 'None'

	print "start time: " + strftime("%Y-%m-%d %H:%M:%S", start_time)
	print "compute node: " + env_vars['HOSTNAME']
	print "current directory: " + env_vars['PWD']
	print "current group/project id: " + env_vars['PROJ_ID']
	print "job id: " + env_vars['JOB_ID']
	print "job name: " + env_vars['REQNAME']
	print "task index number: " + env_vars['SGE_TASK_ID']

	##### FUNCTION TO RUN #####
	eval(args.cmd)
	
	end_time = localtime()
	process = psutil.Process(os.getpid())
	mem = process.get_memory_info()[0] / float(2 ** 20)
	print 'finish time: ' + strftime("%Y-%m-%d %H:%M:%S", end_time)
	print 'time elapsed: ' + str(int(((mktime(end_time)-mktime(start_time))/3600)%60)) + ':' + str(int(((mktime(end_time)-mktime(start_time))/60)%60)) + ':' + str((mktime(end_time)-mktime(start_time))%60)
	print 'memory used: ' + str(mem) + 'MB'

if __name__ == "__main__":
	main()
