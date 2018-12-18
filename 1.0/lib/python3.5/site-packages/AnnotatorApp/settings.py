import os
import sys
import re
import time
import json
import logging
import subprocess

from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import SeqIO



LOCAL_DATABASE = os.path.join(os.getcwd(), "localDB")

APP_NAME="ResAnnotator"
SOFTWARE_VERSION = "0.5.0"
SOFTWARE_SUMMARY = 'A program to annotate antibiotic resistance genes/proteins in contigs/protein sequences \
using BLAST(Uniprot, CARD, Resfinder) and/orHMM (Resfams, PFam-A) databases'

# ====================================================================================
# FUNCTIONS
# ====================================================================================

def determine_path():
	try:
		root = __file__
		if os.path.islink(root):
			root = os.path.realpath(root)
		
		return os.path.dirname(os.path.abspath(root))
	except:
		print("I'm sorry, but something is wrong.")
		print("There is no __file__ variable. Please contact the author.")
		sys.exit()


def file_exists_readable(file, raise_IOError=None):
	"""
	Exit with error if file does not exist or is not readable
	Or raise an IOerror if selected
	"""

	if not os.path.isfile(file):
		message="Can not find file "+ file
		logger.critical(message)
		if raise_IOError:
			print("CRITICAL ERROR: " + message)
			raise IOError
		else:
			sys.exit("CRITICAL ERROR: " + message)

	if not os.access(file, os.R_OK):
		message="Not able to read file " + file
		logger.critical(message)
		if raise_IOError:
			print("CRITICAL ERROR: " + message)
			raise IOError
		else:
			sys.exit("CRITICAL ERROR: " + message)


def find_exe_in_path(exe):
	"""
	Check that an executable exists in $PATH
	"""
	paths = os.environ["PATH"].split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path,exe)
		if os.path.exists(fullexe):
			if os.access(fullexe,os.X_OK):
				return fullexe	
	return


def is_fasta(filename):
	"""Determine if the given file is in Fasta format"""
	with open(filename, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		return any(fasta)


def check_headers(filepath):
    """
        Check if FASTA headers contain white spaces.
        Return True, if header has spaces else returns False.
    """
    header_pattern = re.compile('^>.*')
    for record in SeqIO.parse(filepath,"fasta"):
        header_line = str(record.description.rstrip())
        if re.match('.*\s.*', header_line):
            return False
    
    return True

def fix_headers(infile, outfile):
    """Remove white spaces from the headers of a FASTA file."""
    
    count = 0
    cntSeqs = 0
    fasta_seqs = {}
    with open(outfile, 'w') as f_out:
        for seq_record in SeqIO.parse(open(infile, mode='r'), 'fasta'):
            # remove .id from .description record (remove all before first space)
            seq_record.id = re.sub('\s+', '_', seq_record.description.rstrip())
            seq_record.description = seq_record.id
            #print('SequenceID = '  + seq_record.id)
            #print('Description = ' + seq_record.description + '\n')
            # write new fasta file
            count=SeqIO.write(seq_record, f_out, 'fasta')
            if count!=1: 
                print('Error while writing sequence:  ' + seq_record.id)
            else:
                cntSeqs +=1
   
    return cntSeqs 

def is_hmm(filename):
    """Determine if the file is in hmm format"""
    #print("Inside the hmm function {}".format(filename))
    flag = True 
    return flag 


def run_sbatch_script(script_file, split_cnt, workdir, sbatch_params):
	"""
	This function will submit jobs to cluster
	"""
	logger.info("Submitting jobs to cluster")
	job_status = "FAIL"
	sbatch_script = os.path.join(script_path, 'submit_job_array_slurm.sh')
	chkjob_script = os.path.join(script_path, 'check_job_completion_MB.py')
	job_status = []
	jobId = "" 

	if not sbatch_params:
		sbatch_params = '--mem=8000M'

	if not split_cnt:
		split_cnt = 1

	cmd = "{} {} {} | sbatch {} --workdir {}".format(sbatch_script, script_file, str(split_cnt), str(sbatch_params), str(workdir));
	#print(cmd)
	try:
		output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		#print(output)
		matchObj = re.search('^Submitted batch job (\d+)', output.decode())

		if matchObj:
			jobId = matchObj.group(1)
			logger.debug("jobId = " + jobId)
		else:
			print("No match!")
			logger.exception("Error occured: No match!!. Could not get the job Id")
			sys.exit(1)
    
	except subprocess.CalledProcessError as e:
		logger.error("Error submitting script to cluster = {}".format(cmd))
		logger.error("Error Message = {}".format(str(e)))
		sys.exit(1)

	if jobId:
		cmd2 = "sbatch --dependency=afterany:{} {} --job {} --workdir {}".format(str(jobId), chkjob_script, str(jobId), str(workdir))
		#print(cmd2)
		try:
			stat = subprocess.check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
			filename = os.path.join(str(workdir), "check_temp_" + str(jobId))
			logger.debug("filename = {}".format(str(filename)))
			
			while not os.path.exists(filename):
				#logger.debug("Waiting for a job to finish!")
				time.sleep(5)

			if os.path.exists(filename):
				job_status = "SUCCESS"
			else:
				logger.info("job failed!")
				sys.exit(1)

		except subprocess.CalledProcessError as e:
				logger.error("Error submitting script to cluster = {}".format(cmd2))
				logger.error("Error Message = {}".format(str(e)))
				sys.exit(1)
	else:
		logger.exception("Invalid job number!")
		sys.exit(1)
	
	return job_status;

# ====================================================================================
# FILEPATHS
# ====================================================================================

script_path = determine_path()

path = os.path.join(script_path, "_db/")
data_path = os.path.join(script_path, "_data/")

# ====================================================================================
# LOGGING CONFIG
# ====================================================================================
level = logging.DEBUG
logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s: %(message)s','%m/%d/%Y:%I:%M:%S %p')
#logging.basicConfig(
#	format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', 
#	filemode='w', 
#	datefmt='%m/%d/%Y %I:%M:%S %p')
#
#stream_handler = logging.StreamHandler()
#stream_handler.setFormatter(formatter)
#logger.addHandler(stream_handler)


