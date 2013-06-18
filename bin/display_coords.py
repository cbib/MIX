#! /usr/bin/python
# encoding: utf-8
# Todo : Load, parse and display nucmer alignments in terminal, in a weaving fashion

import sys
import os
# Logger setup
import logging
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('mix')


import argparse
import mummerParser 



CODE={
    'ENDC':0,  # RESET COLOR
    'BOLD':1,
    'UNDERLINE':4,
    'BLINK':5,
    'INVERT':7,
    'CONCEALD':8,
    'STRIKE':9,
    'GREY30':90,
    'GREY40':2,
    'GREY65':37,
    'GREY70':97,
    'GREY20_BG':40,
    'GREY33_BG':100,
    'GREY80_BG':47,
    'GREY93_BG':107,
    'DARK_RED':31,
    'RED':91,
    'RED_BG':41,
    'LIGHT_RED_BG':101,
    'DARK_YELLOW':33,
    'YELLOW':93,
    'YELLOW_BG':43,
    'LIGHT_YELLOW_BG':103,
    'DARK_BLUE':34,
    'BLUE':94,
    'BLUE_BG':44,
    'LIGHT_BLUE_BG':104,
    'DARK_MAGENTA':35,
    'PURPLE':95,
    'MAGENTA_BG':45,
    'LIGHT_PURPLE_BG':105,
    'DARK_CYAN':36,
    'AUQA':96,
    'CYAN_BG':46,
    'LIGHT_AUQA_BG':106,
    'DARK_GREEN':32,
    'GREEN':92,
    'GREEN_BG':42,
    'LIGHT_GREEN_BG':102,
    'BLACK':30,
}

def termcode(num):
    return '\033[%sm'%num

def colorstr(astr,color):
    return termcode(CODE[color])+astr+termcode(CODE['ENDC'])

def parse_field(field, nb_arg):
	##
	# @return the value(s) contained in the MUMmer cell
	field=field.split()
	i=0
	while i <len(field):
		if field[i]=='' or "":
			field.pop(i)
		else:
			i+=1
	if nb_arg==2:
		return field[0], field[1]
	if nb_arg==1:
		return field[0]

def parse_line (line):
	## 
	# @brief Treatment of one line of a COORD file. 
	# @param line string representing one line of a COORD file
	# @return all the characteristics of the corresponding alignment
	tab=[x for x in line.split() if x!="|"]
	if len(tab)!=13:
		print "Unexpected number of items in line", line
		print "Got",len(tab),tab
		sys.exit(1)
	return tab 

def parse_mummerFile (file_adr):
	##
	# @brief Read the COORD file 'file_adr' to fill and return a table containing all the alignments between two different contigs.
	# @param file_adr file containing the alignments between two assemblies (MUMmer COORD file)
	# @return a table containing all the alignments between two different contigs
	alignments=[]
	results_file = open(file_adr, "r")
	line=results_file.readline()
	while line and ((len(line.strip())==0) or not line.split()[0].isdigit()):
		line=results_file.readline()

	while line :
		aln={}
		# try:
		aln["S1"],aln["E1"],aln["S2"],aln["E2"],aln["LEN1"],aln["LEN2"],aln["IDY"],aln["LENR"],aln["LENQ"],aln["COVR"],aln["COVQ"],aln["TAGR"],aln["TAGQ"]=parse_line(line)
		# except ValueError:
		# 	print "Unexpected number(",len(parse_line(line)),"/11 of elements in line",line
		# 	print "Ensure nucmer output was processed with show-coords -l -c DELTAFILE"
		# 	sys.exit(1)


		for k,v in aln.items():
			if k.startswith("TAG"):
				continue
			try:
				aln[k]=int(v)
			except :
				try:
					aln[k]=float(v)
				except:
					pass
		if aln["TAGR"] != aln["TAGQ"] :
			alignments.append(aln)
		line = results_file.readline()
	results_file.close()
	return alignments

def print_aligned_contigs(alignments,use_contig_orientation=False,reverse_mix_orientation=False,max_width=110):
	if use_contig_orientation:
		print colorstr("Reference orientation used, alignments can be weaved if on the same background color","RED")
	else:
		print colorstr("NO Reference orientation used, alignments can be weaved if \n\t* displayed with different background color and are read in oposite directions; \n\t* Or displayed with the same background color and are read in the same direction","RED")

	these_contigs = set([x['TAGQ'] for x in alignments]).union(set([x['TAGR'] for x in alignments]))
	contig_lenghts={}

	for a in alignments:
		contig_lenghts[a['TAGQ']]=int(a['LENQ'])
		contig_lenghts[a['TAGR']]=int(a['LENR'])
	largest_length = max(contig_lenghts.values())
	factor = float(max_width) / largest_length 
	align_index=0
	align_keys="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890œ∑†¥øπ∆©ƒ∂ßåΩç√∫µŒ‰ÁØ∏ÅÍÎÏÓÓÔÒÚÂ˜ı◊ıÇ"
	for c in these_contigs:
		if use_contig_orientation:
			contig_orientation=my_assembly.GRAPH.GRAPH.node.get(c,{}).get('orientation','NA')
			if reverse_mix_orientation:
				if contig_orientation=="reverse":
					contig_orientation="forward"
				elif contig_orientation=="forward":
					contig_orientation="reverse"


		else:
			contig_orientation=""
		print "\n\n-------",c,contig_orientation,contig_lenghts[c],"nt"
		totlen=int(factor*contig_lenghts[c])
		for align_index in range(len(alignments)):
			an_align=alignments[align_index]
			if an_align['TAGQ']==c:
				start=int(int(an_align['S2'])*factor)
				end=int(int(an_align['E2'])*factor)
				alen=int(int(an_align['LENQ'])*factor)
				label="_q"
			elif an_align['TAGR']==c:
				start=int(int(an_align['S1'])*factor)
				end=int(int(an_align['E1'])*factor)
				alen=int(int(an_align['LENR'])*factor)
				label="_r"
			else:
				continue
			k=align_keys[align_index]
			if end <= start:
				# start,end = totlen-start,totlen-end
				start,end = end,start
				if contig_orientation=="reverse":
					print (str(align_index)+label).ljust(12)+"-   "+"_"*(totlen-end-1)+(k*(end-start))+"_"*(start-1) # Reversed alignment on a reversed contig
				else:
					print (str(align_index)+label).ljust(12)+"-   "+"_"*(start-1)+colorstr(k*(end-start),"INVERT")+"_"*(totlen-end) # Reversed alignment on a forward contig 
			else:
				if contig_orientation=="reverse":
					print (str(align_index)+label).ljust(12)+"+   "+"_"*(totlen-end)+colorstr(k*(end-start),"INVERT")+"_"*(start-1) # Forward alignments on a reversed  contig 
				else:
					print (str(align_index)+label).ljust(12)+"+   "+"_"*(start-1)+k*(end-start)+"_"*(totlen-end) # Forward alignments on a Forward contig
				




def main():

	parser=argparse.ArgumentParser(description="Display nucmer alignments in ASCII")
	parser.add_argument('-w',dest="max_width",action="store",type=int,default=110,help="Width of text output")
	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,help="Nucmer coords file")


	args=parser.parse_args()
	alignments=mummerParser.parse_mummerFile(args.infile)
	# Filter out self alignments
	alignments=[x for x in alignments if x['TAGQ']!=x['TAGR']]
	print "Found",len(alignments)
	if len(alignments)<1:
		sys.exit(0)
	print_aligned_contigs(alignments,max_width=args.max_width)


if __name__ == "__main__":
	main()
