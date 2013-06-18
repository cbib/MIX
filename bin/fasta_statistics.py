#! /usr/bin/python
# encoding: utf-8
# Todo : Load, parse and display nucmer alignments in terminal, in a weaving fashion

import sys
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
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")
	parser.add_argument('FASTAFILE',action='append',nargs="+",help='List of fasta files to keep. Use "*" to keep them all')
	args=parser.parse_args()
	all_records=[]
	FASTAFILE=args.FASTAFILE[0]
	if args.pretty:
		import prettytable

	for f in FASTAFILE: 
		for record in SeqIO.parse(f, "fasta", generic_dna):
			all_records.append(SequenceStat(f,record))
	# Display summary statistics per file
	sequences_per_files=collections.defaultdict(list)
	for s in all_records:
		sequences_per_files[s.file].append(s)
	if args.pretty:
		table=prettytable.PrettyTable(["File","#Seqs","Avg GC","Avg Length","Total Length"])
		table.align["File"] = "l" 

		for file,seqs in sequences_per_files.items():
			table.add_row([file,len(seqs),scipy.average([x.gc for x in seqs]),\
				scipy.average([x.length for x in seqs]),sum([x.length for x in seqs])])
		print >>args.outfile,table.get_string(sortby="Total Length")

	else:
		for file,seqs in sequences_per_files.items():
			print >>args.outfile," ".join(map(str,[\
				file,len(seqs),scipy.average([x.gc for x in seqs]),\
				scipy.average([x.length for x in seqs]),sum([x.length for x in seqs])
				]))
if __name__ == "__main__":
	sys.exit(main())
