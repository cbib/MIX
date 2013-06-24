#! /usr/bin/python
# encoding: utf-8


import sys
import N50 
import collections
from Bio.Seq import Seq
from Bio.SeqUtils import CheckSum
import scipy
from Bio import SeqUtils
from Bio import SeqIO

from Bio.Alphabet import generic_dna

import sys
import os
import argparse
import logging
from scipy.stats.mstats import mquantiles

logger = logging.getLogger('concat')
logger.setLevel(logging.DEBUG)

# while len(logger.handlers()) > 0:
#  	logger.pop()

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(funcName)s - %(message)s',"%Y-%m-%d %H:%M:%S")
# formatter = logging.Formatter('%(asctime)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)



help="""Usage 
%(tool)s  FASTAFILE [-o]

Options
========
-o OUTPUTFILE: Indicate output file (default stdout)

Example usage
==============
%(tool)s genome1.fa -o 

Todo
======
"""%{"tool":sys.argv[1]}

def main(argv=None):
	parser=argparse.ArgumentParser(description="Generate a single contig by stitching together all contigs of FASTAFILE")
#	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")
	parser.add_argument('FASTAFILE',action='append',nargs="+",help='List of fasta files to merge')
	args=parser.parse_args()
	logger.info("Will output to %s",str(args.outfile))
	concatenated_sequence=""
	sequence_id=""
	for afile in args.FASTAFILE:
		afile=afile[0]
		for record in SeqIO.parse(afile, "fasta", generic_dna):
			concatenated_sequence+=str(record.seq)
		this_file_id = afile.replace("/","_")
		this_file_id = this_file_id.replace(" ","_")
		sequence_id+=this_file_id
	logger.info("Concatenated %d bp",len(concatenated_sequence))
	print >>args.outfile,">%s\n"%(sequence_id)
	print >>args.outfile,concatenated_sequence
	print >>args.outfile,"\n"
	args.outfile.close()
if __name__ == "__main__":
	sys.exit(main())
