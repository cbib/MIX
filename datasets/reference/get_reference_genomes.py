#!/usr/bin/env python
# encoding: utf-8

import ftplib
from subprocess import call

import os
import logging

from Queue import *
import threading
n_parallel_connections=4
genomes_to_download=[
# NCBI Directory URL 																	# Local name
["genomes/Bacteria/Aeromonas_hydrophila_ATCC_7966_uid58617"						,"Aeromonas_hydrophila"],
["genomes/Bacteria/Bacillus_cereus_ATCC_10987_uid57673"						,"Bacillus_cereus"],
["genomes/Bacteria/Bacteroides_fragilis_638R_uid84217"						,"Bacteroides_fragilis"],
["genomes/Bacteria/Mycobacterium_abscessus_uid61613"						,"Mycobacterium_abscessus"],
["genomes/Bacteria/Rhodobacter_sphaeroides_2_4_1_uid57653"						,"Rhodobacter_sphaeroides"],
["genomes/Bacteria/Staphylococcus_aureus_USA300_TCH1516_uid58925"						,"Staphylococcus_aureus"],
["genomes/Bacteria/Vibrio_cholerae_O1_biovar_El_Tor_N16961_uid57623"						,"Vibrio_cholerae"],
["genomes/Bacteria/Xanthomonas_axonopodis_citrumelo_F1_uid73179"						,"Xanthomonas_axonopodis"]
]

#genomes_to_download=[{"ncbi_url":x[0],"local_name":x[1],"extension":".fna"} for x in genomes_to_download]
genomes_to_download=[{"ncbi_url":x[0],"local_name":x[1],"extension":".gff"} for x in genomes_to_download]

logger = logging.getLogger('simple_example')
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - THR %(thread)d - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)




class Worker(threading.Thread):
	def __init__(self, function, in_queue, out_queue):
		self.function = function

		self.in_queue, self.out_queue = in_queue, out_queue
		super(Worker, self).__init__()

	def run(self):
		while True:
			try:
				if self.in_queue.empty(): 
					break
				data = self.in_queue.get()
				result = self.function(data)
				self.out_queue.put(result)
				self.in_queue.task_done()
				logger.info("Still %d to do",self.in_queue.qsize())
			except Exception as e:
				logger.critical('something happened!: %s',repr(e))
				self.out_queue.put("result")
				self.in_queue.task_done()
#				print "exception in thread:",repr(e)
				


def process(data, function, num_workers=1):
	in_queue = Queue()
	for item in data:
		in_queue.put(item)
	out_queue = Queue(maxsize=in_queue.qsize())
	workers = [Worker(function, in_queue, out_queue) for i in xrange(num_workers)]
	for worker in workers: 
		worker.setDaemon(True)
		worker.start()

	in_queue.join()

	return out_queue


def download_genome(entry):
	"""Entry is an entry of the genomes_to_download list"""
	logger.info("Will get %s  to %s",entry['ncbi_url'],entry['local_name'])
	ftp = ftplib.FTP("ftp.ncbi.nih.gov")
	ftp.login("anonymous","massyah@gmail.com")
	ftp.cwd(entry['ncbi_url'])
	logger.info("Initiated connection")

	content=[]
	ftp.dir(content.append)
	sequences=[x.split()[-1] for x in content if x.endswith(entry["extension"])]
	logger.info("Listed directory: %d sequences to download",len(sequences))

	content=[]
	ftp.dir(content.append)

	try:
		os.mkdir(entry['local_name'])
	except OSError:
		pass

	try:
		os.mkdir(entry['local_name'])
	except OSError:
		pass

	for seqname in sequences:
		file=open(entry['local_name']+"/"+seqname,"w")
		ftp.retrbinary("RETR "+seqname,file.write)
		file.close()
		logger.info("Downloaded %s",seqname)
	


	return "OK"+entry['local_name']


def concatenate_genome(entry):
	local_dir=entry['local_name']
	logger.info("Concatenating %s",local_dir)
	sequences=[x for x in os.listdir(local_dir) if x.endswith(entry["extension"])]
	if entry["extension"]==".fna":
		target_name=local_dir+"_ref.fa"
	else:
		target_name=local_dir+"_ref"+entry["extension"]

	# call(['cat']+sequences+['> %s'%(target_name)])
	CMD="cd %s && cat %s > %s"%(local_dir," ".join(sequences),target_name)
	# logger.info("Will run CMD:\n%s",CMD )
	os.system(CMD)
	logger.info("Finished concatenatioon %s",local_dir)


process(genomes_to_download,download_genome,n_parallel_connections)
process(genomes_to_download,concatenate_genome,n_parallel_connections)
