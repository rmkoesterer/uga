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

import sys
import psutil
import subprocess
import logging
import os

qsub_script_prefix = """#!/bin/bash
#$ -cwd

if [ ! -z "$SGE_TASK_ID" ]
then
taskid=$SGE_TASK_ID
else
taskid="None"
fi

if [ ! -z "$REQNAME" ]
then
reqname=$REQNAME
else
reqname="None"
fi

if [ ! -z "$HOSTNAME" ]
then
hostname=$HOSTNAME
else
hostname="None"
fi

if [ ! -z "$JOB_ID" ]
then
jobid=$JOB_ID
else
jobid="None"
fi

if [ ! -z "$SGE_CLUSTER_NAME" ]
then
clustername=$SGE_CLUSTER_NAME
else
clustername="None"
fi

if [ ! -z "$HOST" ]
then
host=$HOST
else
host="None"
fi

if [ ! -z "$Queue" ]
then
queue=$Queue
else
queue="None"
fi

if [ ! -z "$PWD" ]
then
pwd=$PWD
else
pwd="None"
fi

if [ ! -z "$SGE_STDOUT_PATH" ]
then
stdoutpath=$SGE_STDOUT_PATH
else
stdoutpath="None"
fi

user=$USER
"""

qsub_script_postfix = """
EXIT_CODE=$?
	
exit $EXIT_CODE
"""

class Error(Exception):
	def __init__(self, msg):
		self.out = 'ERROR: ' + msg
		self.msg = msg

def print_error(e):
	return "\nERROR: " + e + "\n"

def write_qsub_script(singularity_cmd, qsub_wrapper, cmd, qsub_script, jobs_run_file, log_file):
	text = [qsub_script_prefix, singularity_cmd + " " + qsub_wrapper + " --user $user --taskid $taskid --reqname $reqname --hostname $hostname --jobid $jobid --clustername $clustername --host $host --queue $queue --pwd $pwd --stdoutpath $stdoutpath --cmd " + '\"' + cmd + '\"' + " --job-list " + jobs_run_file + " --log " + log_file, qsub_script_postfix]
	with open(qsub_script, 'w') as f:
		f.write("\n".join(text).encode('utf-8'))
	os.chmod(qsub_script, 0775)

def qsub(qsub_pre, singularity_cmd, qsub_wrapper, cmd, qsub_script, jobs_run_file, log_file):
	write_qsub_script(singularity_cmd, qsub_wrapper, cmd, qsub_script, jobs_run_file, log_file)
	cmd_list = qsub_pre + [os.path.abspath(qsub_script)]
	try:
		p = subprocess.Popen(cmd_list,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
		for line in iter(p.stdout.readline, ''):
			sys.stdout.write(line)
		p.wait()
	except KeyboardInterrupt:
		kill_all(p.pid)
		print "   ... process terminated by user"
		sys.exit(1)

def interactive(submit, cmd, log_file = None):
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
		print "   ... process terminated by user"
		sys.exit(1)

def kill_all(pid):
	parent = psutil.Process(pid)
	for child in parent.children(recursive=True):
		child.kill()
	parent.kill()
