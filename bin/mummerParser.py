import sys

import logging
if "logger" not in globals():
	logger = logging.getLogger('mix_logger')
	logger.setLevel(logging.DEBUG)

	# while len(logger.handlers()) > 0:
	#  	logger.pop()

	# create console handler and set level to debug
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)

	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s',"%Y-%m-%d %H:%M:%S")
	# formatter = logging.Formatter('%(asctime)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)

	# add ch to logger
	logger.addHandler(ch)


logger = logging.getLogger('mix_logger')



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

def parse_mummerFile (file_adr,skip_self=True):
	##
	# @brief Read the COORD file 'file_adr' to fill and return a table containing all the alignments between two different contigs.
	# @param file_adr file containing the alignments between two assemblies (MUMmer COORD file)
	# @return a table containing all the alignments between two different contigs
	alignments=[]
	if type(file_adr)==file:
		results_file=file_adr
	else:
		results_file = open(file_adr, "r")
	# Skip header lines 
	line=results_file.readline()
	while line and ((len(line.strip())==0) or not line.split()[0].isdigit()):
		line=results_file.readline()

	while line :
		if(len(alignments) % 1000)==0:
			logger.debug("Parsed %d alignments"%(len(alignments)))
		# if ((len(line.strip())==0) or not line.split()[0].isdigit()):
		# 	line = results_file.readline()
		# 	continue

		aln={}
		try:
			aln["S1"],aln["E1"],aln["S2"],aln["E2"],aln["LEN1"],aln["LEN2"],aln["IDY"],aln["LENR"],aln["LENQ"],aln["COVR"],aln["COVQ"],aln["TAGR"],aln["TAGQ"]=parse_line(line)
		except ValueError:
			print "Unexpected number(",len(parse_line(line)),"/11 of elements in line",line
			print "Ensure nucmer output was processed with show-coords -l -c DELTAFILE"
			sys.exit(1)
		if skip_self and (aln["TAGQ"]==aln['TAGR']):
			line = results_file.readline()
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
		alignments.append(aln)
		line = results_file.readline()
	results_file.close()
	return alignments

def print_alignment(aln,outStream=sys.stdout,pretty=False):
	if pretty:
		elements=map(str,[aln["S1"],aln["E1"],"|",aln["S2"],aln["E2"],"|",aln["LEN1"],aln["LEN2"],"|",aln["IDY"],"|",aln["LENR"],aln["LENQ"],"|",aln["COVR"],aln["COVQ"],"|",aln["TAGR"],aln["TAGQ"]])
		for i in range(len(elements)):
			if elements[i]=="|":
				continue
			elements[i]=elements[i].ljust(10)
		print >>outStream," ".join(elements)
	else:
		print >>outStream,aln["S1"],aln["E1"],"|",aln["S2"],aln["E2"],"|",aln["LEN1"],aln["LEN2"],"|",aln["IDY"],"|",aln["LENR"],aln["LENQ"],"|",aln["COVR"],aln["COVQ"],"|",aln["TAGR"],aln["TAGQ"]

