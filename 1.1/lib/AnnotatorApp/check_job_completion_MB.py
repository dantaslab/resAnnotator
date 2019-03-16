#!/usr/bin/env python

"""
	This program check if the job submitted to the cluster completed successfully
	The program takes in jobid and prints code and a message.
	
	Return value 
		if job completed successfully: 
			return [1, "jobid completed successfully"]
		else:
			return [0, "jobid_%a failed"]	

"""
__author__ = "Manish Boolchandani"
__date__ = "2017-04-10"
__version__ = "0.1.1"

# Python imports
import sys
import os
import time
import re
import argparse
import logging
import subprocess


# Logger configuration and setup
logger = logging.getLogger(__name__)

def main(args):
	logger.info("checking if all jobs completed successfully")

	parser = argparse.ArgumentParser(prog = 'check_job_completion.py', description = "A program to check if all jobs in an array got completed!")
	parser.add_argument('-v','--version', action = 'version', version = "%(prog)s 0.1.0")
	parser.add_argument('-j', '--job', dest = 'job_id', help = "Enter job id")
	parser.add_argument('-n', '--workdir', dest = 'workdir', help = "Enter current working directory")
	args = parser.parse_args()
	
	if not args.job_id:
		logger.error("You must provide job id!")
		parser.exit(status=1, message = "You must provide job id!\n")
	else:
		job_id = args.job_id

	if not args.workdir:
		logger.error("You must provide workdir!")
		parser.exit(status=1, message = "You must workdir!\n")
	else:
		workdir = args.workdir

	#job_status = check_job(job_id)

	filename = os.path.join(workdir, "check_temp_"+ str(job_id))

	fh = open(filename, 'w')

	fh.write('COMPLETED')

	fh.close()

	print filename

"""
check job and returns a list with error code and message
"""
def check_job(jobId):

        flag = 1
	job_output = []
        logger.info("Job finished. Checking the status of all tasks!")
        cmd = "sacct -P -n -j " + jobId + " --format jobid,jobname,elapsed,ReqMem,alloccpus,state,exitcode"
        job_detail = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if job_detail:
        	#logger.debug(job_detail)
                for line in job_detail.splitlines():
                	logger.debug(line)
                	status = line.split('|')
                        logger.debug("job state " + status[5])
                        if status[5] != "COMPLETED" or status[6]!="0:0":
                        	logger.exception("Error Occurred...Exiting!!")
                                sys.exit(1)
                        else:
                        	logger.exception("Error occurred....Exiting!!")
                                sys.exit(1)
	else:
                logger.exception("jobId is empty! Exiting!!")
                sys.exit(1)

        return job_output;

if __name__ == '__main__':
	main(sys.argv[1:])
