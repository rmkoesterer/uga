#!/bin/bash
#$ -cwd
#$ -j y
#$ -V

source ~/.bashrc
use UGER
echo "===================================================="
echo "Started on: $(date)"
echo "Node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Current Group/Project ID: $(id -ng)"
echo "Job ID: $JOB_ID"
echo "Job name: $JOB_NAME"
echo "Task index number: $SGE_TASK_ID"
echo "===================================================="
stime=$(date '+%s')

for arg in "$@"; do
	echo -e "...command submitted: ${arg}\n"
	eval ${arg}
done

maxvmem=`qstat -j $JOB_ID | grep maxvmem | awk '{print $8}' | awk -F'=' '{print $2}'`
etime=$(date '+%s')
dtime=$((etime - stime))
dsec=$((dtime % 60))
dmin=$(((dtime / 60) % 60))
dhour=$(((dtime / 3600) % 60))
echo -e "\n===================================================================="
echo "Finished on: $(date)"
echo "Time elapsed: $dhour hours, $dmin minutes, $dsec seconds"
echo "Maximum memory used: $maxvmem"
echo "===================================================================="
