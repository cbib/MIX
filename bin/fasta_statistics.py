#! /usr/bin/env python
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
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('mix')


help="""Usage 
%(tool)s  FASTAFILE FASTAFILE 

Options
========
-p: Turn on pretty table output, require Pretty Table module (https://pypi.python.org/pypi/PrettyTable/)
-d : Indicate output separator, default to space 
Example usage
==============
%(tool)s genome1.fa genome2.fa

python bin/fasta_statistics.py datasets/*/*/*.fa datasets/*/*/*.fasta

Todo
======
"""%{"tool":sys.argv[1]}

class SequenceStat(object):
	"""docstring for SequenceStat"""
	def __init__(self, file,fastaRecord):
		super(SequenceStat, self).__init__()
		self.file = file
		self.length = len(fastaRecord.seq)
		self.description=fastaRecord.description
		self.gc=SeqUtils.GC(fastaRecord.seq)
		self.crc32=CheckSum.crc32(fastaRecord.seq)
		
def main(argv=None):
	parser=argparse.ArgumentParser(description="Compute various statistics related to the sequences either in the provided fasta files or for the sequences piped in")
#	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('-p',dest="pretty",action="store_true",help="Pretty print using PrettyTable module")
	parser.add_argument('-d',dest="delimiter",help="Colum separator for output, default to whitespace",default=" ")
	parser.add_argument('-t',dest="min_length",help="Minimun length threshold to filter fasta file",default=0,type=int)
	parser.add_argument('-r',dest="reference_length",help="(Not yet implemented)Reference length used to compute corrected Nx values",default=0)
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")
	parser.add_argument('FASTAFILE',action='append',nargs="+",help='List of fasta files to keep. Use "*" to keep them all')
	args=parser.parse_args()
	all_records=[]
	FASTAFILE=args.FASTAFILE[0]
	if args.pretty:
		import prettytable

	for f in FASTAFILE: 
		for record in SeqIO.parse(f, "fasta", generic_dna):
			if len(record.seq)<=args.min_length:
				continue
			all_records.append(SequenceStat(f,record))
	# Display summary statistics per file
	sequences_per_files=collections.defaultdict(list)
	for s in all_records:
		sequences_per_files[s.file].append(s)
	if args.pretty:
		table=prettytable.PrettyTable(["File","#Seqs","Avg GC","Avg Length(kb)", "Quant","min","max",  "Sum Length(kb)","N50(kb)","L50"])
		table.align["File"] = "l" 

		for file,seqs in sequences_per_files.items():
			lengths=[x.length for x in seqs]
			table.add_row([file,len(seqs),round(scipy.average([x.gc for x in seqs]),2),\
				round(scipy.average(lengths)/1000,2),mquantiles(lengths),min(lengths),max(lengths),round(sum(lengths)/1000,2),round(N50.N50(lengths)/1000,2),N50.L50(lengths)])
		print >>args.outfile,table.get_string(sortby="N50(kb)")

	else:
		for file,seqs in sequences_per_files.items():
			lengths=[x.length for x in seqs]

			print >>args.outfile," ".join(map(str,[\
				file,len(seqs),scipy.average([x.gc for x in seqs]),\
				scipy.average(lengths),sum(lengths),N50.N50(lengths),N50.L50(lengths)
				]))
if __name__ == "__main__":
	sys.exit(main())
