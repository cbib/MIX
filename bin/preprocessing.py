#! /usr/bin/python
import sys
import re
from optparse import OptionParser 

def process_contig_renaming(ass_adr, i, output_file):
	print ass_adr, " will correspond to the assembly number ", i
	ass_file = open(ass_adr, "r")
	line = ass_file.readline()
	while line : 
		if re.match("^>", line) : 
			identifiant = line.replace(">", ">A"+str(i)+"_")
			line = ass_file.readline()
			output_file.write(identifiant)
			while line and not re.match("^>", line) :
				output_file.write(line)
				line = ass_file.readline()
		else : 
			line = ass_file.readline()
	ass_file.close()

def main():
	parser = OptionParser()
	parser.add_option("-o", "--out", dest="out", help ="output file concatenating the assemblies")
	(options, args) = parser.parse_args()
	out_adr = options.out
	output_file = open(out_adr, "w")
	if sys.argv[1]=="-o":
		for i in range (3, len(sys.argv), 1) : 
			process_contig_renaming(sys.argv[i], i-2, output_file)
	else : 
		for i in range (1, len(sys.argv)-2, 1) : 
			process_contig_renaming(sys.argv[i], i, output_file)
	output_file.close()

if __name__ == "__main__":
    main()

