#! /usr/bin/python
# encoding: utf-8
# Todo : Load, parse and display nucmer alignments in terminal, in a weaving fashion

import sys

import sys
import os
import argparse
import mummerParser 

# Logger setup
import logging
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('mix')


help="""Usage 
%(tool)s [-er] CONTIGID CONTIGID CONTIGID ... NUCMERCOORDFILE



Options
=========
-e: 
-r: 

Example usage
==============
%(tool)s A1 A2 nuc.coords

Todo
======
"""%{"tool":sys.argv[1]}

def main(argv=None):
	parser=argparse.ArgumentParser(description="Pipe out the alignments involving or excluding provided contigs ID")
	parser.add_argument('-e',dest="inverse_match",action="store_true",help="Inverse the criterion, filter out any conting in CONTIGID list (not implemented yet)")
	parser.add_argument('-r',dest="use_re",action="store_true",help="Consider CONTIGID as python regexp (not implemented yet)")
	parser.add_argument('-s',dest="remove_self",action="store_true",help="remove self alignments")

	parser.add_argument('-p',dest="pretty",action="store_true",help="Pretty print")
	parser.add_argument('-v',dest="verbose",action="store_true",help="Verbose output")
	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")
	parser.add_argument('CONTIGS',action='append',nargs="+",help='List of CONTIGS to keep. Use "*" to keep them all')
	args=parser.parse_args()

	contig_ids=args.CONTIGS[0]
	alignments=mummerParser.parse_mummerFile(args.infile)
	kept_alignments=[]
	for a in alignments:
		if "*" in contig_ids and ((args.remove_self and a['TAGQ']!=a['TAGR'])or not args.remove_self):
			kept_alignments.append(a)
		if a['TAGQ'] in contig_ids and a['TAGR'] in contig_ids:
			kept_alignments.append(a)
		elif args.verbose:
			print >>args.outfile,"Skipped",a
	for a in kept_alignments:
		mummerParser.print_alignment(a,args.outfile,args.pretty)

if __name__ == "__main__":
	sys.exit(main())
