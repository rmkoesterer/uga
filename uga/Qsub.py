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

	from uga.__version__ import version

	print "uga v" + version
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

	if argv[1].split('(')[0] == "RunSnv":
		from uga.RunSnv import RunSnv
	if argv[1].split('(')[0] == "RunSnvgroup":
		from uga.RunSnvgroup import RunSnvgroup
	if argv[1].split('(')[0] == "RunMeta":
		from uga.RunMeta import RunMeta
	if argv[1].split('(')[0] == "RunSnvplot":
		from uga.RunSnvplot import RunSnvplot
	if argv[1].split('(')[0] == "RunSnvgroupplot":
		from uga.RunSnvgroupplot import RunSnvgroupplot
	if argv[1].split('(')[0] == "RunFilter":
		from uga.RunFilter import RunFilter

	exec('r=' + argv[1])
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
	main(sys.argv)
