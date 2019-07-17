#!/usr/bin/env python
#$ -cwd
#$ -j y
#$ -V

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

import os
import pwd
import psutil
import resource
import sys
import argparse
from time import strftime, localtime, time, gmtime

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def main(args=None):

	start_time = (localtime(), time())

	if args.cmd.split('(')[0] in ["RunSnv","RunSnvgroup","RunMeta","RunMerge","RunTools"] and args.taskid != 'None':
		with open(args.job_list) as f:
			joblist = [line.rstrip() for line in f]
		job = joblist[int(args.taskid)-1]
		args.cmd = args.cmd.replace("UGA_JOB_ID",job)
		args.log = args.log.replace("UGA_JOB_ID",job)
		args.cmd = args.cmd.replace("UGA_JOB_RANGE",str((100 * ((int(job)-1) / 100) + 1)) + "-" + str((100 * ((int(job)-1) / 100) + 100)))
		args.log = args.log.replace("UGA_JOB_RANGE",str((100 * ((int(job)-1) / 100) + 1)) + "-" + str((100 * ((int(job)-1) / 100) + 100)))
		try:
			lf = open(args.log,'w')
		except(IOError, OSError):
			return
		sys.stdout = lf
		sys.stderr = lf
    
	if args.reqname != 'None':
		args.reqname = args.hostname + '_' + strftime('%Y_%m_%d_%H_%M_%S', start_time[0]) if args.hostname != 'None' else strftime('%Y_%m_%d_%H_%M_%S', start_time[0])

	from uga.__version__ import version
    
	print "uga v" + version
	print "start time: " + strftime("%Y-%m-%d %H:%M:%S", start_time[0])
	print "sge cluster: " + args.clustername
	print "host: " + args.host
	print "queue: " + args.queue
	print "compute node: " + args.hostname
	print "current directory: " + args.pwd
	print "user name: " + args.user
	print "job id: " + args.jobid
	print "job name: " + args.reqname
	print "task index number: " + args.taskid
	print "sge log file: " + args.stdoutpath
    
	if args.cmd.split('(')[0] == "RunSnv":
		from uga.RunSnv import RunSnv
	if args.cmd.split('(')[0] == "RunSnvgroup":
		from uga.RunSnvgroup import RunSnvgroup
	if args.cmd.split('(')[0] == "RunMeta":
		from uga.RunMeta import RunMeta
	if args.cmd.split('(')[0] == "RunMerge":
		from uga.RunMerge import RunMerge
	if args.cmd.split('(')[0] == "RunSnvplot":
		from uga.RunSnvplot import RunSnvplot
	if args.cmd.split('(')[0] == "RunSnvgroupplot":
		from uga.RunSnvgroupplot import RunSnvgroupplot
	if args.cmd.split('(')[0] == "RunFilter":
		from uga.RunFilter import RunFilter
	if args.cmd.split('(')[0] == "RunTools":
		from uga.RunTools import RunTools
    
	print ""
	print "command entered: " + args.cmd
	exec('r=' + args.cmd)
	if r == 0:
		end_time = (localtime(), time())
		process = psutil.Process(os.getpid())
		mem=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000.0
		mem_children=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1000.0
		print 'finish time: ' + strftime("%Y-%m-%d %H:%M:%S", end_time[0])
		print 'time elapsed: ' + strftime('%H:%M:%S', gmtime(end_time[1] - start_time[1]))
		print 'max memory used by main process: ' + str('%.2f' % mem) + ' MB'
		print 'max memory used by any subprocess: ' + str('%.2f' % mem_children) + ' MB'

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--user', help='USER')
	parser.add_argument('--taskid', help='SGE_TASK_ID')
	parser.add_argument('--reqname', help='REQNAME')
	parser.add_argument('--hostname', help='HOSTNAME')
	parser.add_argument('--jobid', help='JOB_ID')
	parser.add_argument('--clustername', help='SGE_CLUSTER_NAME')
	parser.add_argument('--host', help='HOST')
	parser.add_argument('--queue', help='Queue')
	parser.add_argument('--pwd', help='PWD')
	parser.add_argument('--stdoutpath', help='SGE_STDOUT_PATH')
	requiredArgs = parser.add_argument_group('required arguments')
	requiredArgs.add_argument('--cmd', help='a command', required=True)
	requiredArgs.add_argument('--job-list', help='a job list', required=True)
	requiredArgs.add_argument('--log', help='a log file name', required=True)
	args = parser.parse_args()
	main(args)
