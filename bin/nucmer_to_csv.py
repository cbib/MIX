#! /usr/bin/python
# encoding: utf-8
# Todo : Load, parse and display nucmer alignments in terminal, in a weaving fashion

import sys

import sys
import os
import argparse
from operator import itemgetter
import mummerParser 

# Logger setup
import logging
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('mix')


help="""Usage 
%(tool)s  NUCMERCOORDFILE



Options
=========
-e: 
-k: 
Example usage
==============
%(tool)s nuc.coords

Todo
======
"""%{"tool":sys.argv[1]}

def get_assembler_name (contig_name):
	if "_" in contig_name :
		return contig_name.split("_")[0].replace(">", "")
	else :
		return "None"


def is_kept(aln):
	idy_threshold = 75 #98
	cov_threshold = 99 #95
	aln_threshold = 125 
	ctg_threshold= 0
	# skip self
	if aln['TAGQ']==aln['TAGR']:
		return False

	# Skip under threshold 
	threshold_aln_OK = (int(aln["LEN1"]) > aln_threshold and int(aln["LEN2"]) > aln_threshold)
	threshold_ctg_OK = (int(aln["LENR"]) > ctg_threshold and int(aln["LENQ"]) > ctg_threshold)
	if (not threshold_ctg_OK) or (not threshold_aln_OK):
		return False

	A1, A2 = get_assembler_name (aln["TAGR"]), get_assembler_name (aln["TAGQ"])
	if (A1==A2):
		return False

	# if (aln["TAGR"] != aln["TAGQ"]) and ((aln["COVR"] > cov_threshold) or ((aln["LENR"] - aln["LEN1"]) < ctg_threshold)):
	# 	return False
	# elif (aln["TAGR"] != aln["TAGQ"]) and ((aln["COVQ"] > cov_threshold) or ((aln["LENQ"] - aln["LEN2"]) < ctg_threshold)):
	# 	return False
	if (aln["TAGR"] != aln["TAGQ"]) and (aln["IDY"] > idy_threshold) and threshold_aln_OK and threshold_ctg_OK and (A1 != A2):
		return True
	else:
		return False
	mummerParser.print_alignment(aln,sys.stderr,pretty=True)
	assert(False)
	# return False

def is_extremal(a):
	boudary_set_off = 52
	ref_boudary_set_off = boudary_set_off
	query_boudary_set_off = boudary_set_off
	#ref_boudary_set_off = (a["LENR"]/10000)+2
	#query_boudary_set_off = (a["LENQ"]/10000)+2
	ref_forward = a["S1"] < a["E1"]
	query_forward = a["S2"] < a["E2"]
	if ref_forward : 
		ref_aln_5 = a["S1"] < ref_boudary_set_off
		ref_aln_3 = a["E1"] > a["LENR"] - ref_boudary_set_off
	else : 
		ref_aln_5 = a["E1"] < ref_boudary_set_off
		ref_aln_3 = a["S1"] > a["LENR"] - ref_boudary_set_off
	if query_forward : 
		query_aln_5 = a["S2"] < query_boudary_set_off
		query_aln_3 = a["E2"] > a["LENQ"] - query_boudary_set_off
	else : 
		query_aln_5 = a["E2"] < query_boudary_set_off
		query_aln_3 = a["S2"] > a["LENQ"] - query_boudary_set_off
	if ref_forward and query_forward : 
		extends = (ref_aln_5 and query_aln_3) or (query_aln_5 and ref_aln_3)
	else : 
		extends = (ref_aln_5 and query_aln_5) or (query_aln_3 and ref_aln_3)	
	return extends

def main(argv=None):
	parser=argparse.ArgumentParser(description="Pipe out the alignments in CSV format")
	parser.add_argument('-e',dest="extremal_only",action="store_true",help="Only keep extremal alignments")
	parser.add_argument('-k',dest="florence_selection",action="store_true",help="Remove alignments not satisfying Florence criteria")

	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,dest="outfile")
	args=parser.parse_args()

	alignments=mummerParser.parse_mummerFile(args.infile)
	# Get col names 
	a=alignments[0]	
	items = a.items()
	items.sort(key=itemgetter(0))
	print>>args.outfile, "\t".join(map(str,[x[0] for x in items]))

	for a in alignments:
		if args.extremal_only:
			if not is_extremal(a):
				continue
		if args.florence_selection:
			if not is_kept(a):
				continue

		items = a.items()
		items.sort(key=itemgetter(0))
		print>>args.outfile, "\t".join(map(str,[x[1] for x in items]))

if __name__ == "__main__":
	sys.exit(main())
