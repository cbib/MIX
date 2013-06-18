import re
import scipy
from IPython.core.debugger import Tracer; debug_here = Tracer()

import os 

def parse_field(field, nb_arg):
	##
	# @return the value(s) contained in the MUMmer cell
	field=field.split(" ")
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
	tab=line.split("|")
	S1, E1=parse_field(tab[0],2)
	S2, E2=parse_field(tab[1],2)
	LEN1,LEN2=parse_field(tab[2],2)
	IDY=parse_field(tab[3],1)
	LENR, LENQ=parse_field(tab[4],2)
	COVR, COVQ=parse_field(tab[5], 2)
	tab2=parse_field(tab[6],1).rstrip().split("\t")
	TAGR,TAGQ=tab2[0],tab2[1]
	return S1,E1,S2,E2,LEN1,LEN2,IDY,LENR,LENQ,COVR,COVQ,TAGR,TAGQ

def parse_mummerFile (file_adr):
	##
	# @brief Read the COORD file 'file_adr' to fill and return a table containing all the alignments between two different contigs.
	# @param file_adr file containing the alignments between two assemblies (MUMmer COORD file)
	# @return a table containing all the alignments between two different contigs
	alignments=[]
	results_file = open(file_adr, "r")
	for i in range (0, 6, 1):
		line = results_file.readline()
	while line :
		aln={}
		aln["S1"],aln["E1"],aln["S2"],aln["E2"],aln["LEN1"],aln["LEN2"],aln["IDY"],aln["LENR"],aln["LENQ"],aln["COVR"],aln["COVQ"],aln["TAGR"],aln["TAGQ"]=parse_line(line)

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



mix_output= "~/Documents/MIX/Mix-1.1/outRhodo_APSPBB_no_removal/Mix_results_A125_C0/Mix_assembly.fasta"
soap="~/Documents/MIX/Rhodobacter_Assembly/SOAPdenovo/genome.ctg.fasta"
coords_file = "~/Documents/MIX/AllpathSoapBambus.coords"
mix_line_re = re.compile(">\d+(.*)")
mix_line_re_2 = re.compile(">\d+_A\d_\w+\.(.*)")
mix_line_re_3 = re.compile("A\d_\w+\.(.*)")
mix_id = set()
handle = open(os.path.expanduser(mix_output), "rU")
for line in handle:
	if line.startswith(">"):
		if "{" in line:
			print line

			rec=eval(mix_line_re.findall(line)[0])
			for r in rec:
			    mix_id.add(mix_line_re_3.findall(r['contig'])[0])
		else:
			mix_id.add(mix_line_re_2.findall(line)[0])
handle.close()

handle = open(os.path.expanduser(soap), "rU")
soap_id=set()
soap_length={}
for record in SeqIO.parse(handle, "fasta") :
	soap_id.add(record.id)
	soap_length[record.id]=len(record.seq)
# for line in open(os.path.expanduser(soap),'r'):
	# if line.startswith(">"):
	# soap_id.add(line.strip()[1:].split(" ")[0])

soap_id.difference(mix_id)
len(soap_id.difference(mix_id))
acc_length=0
for c in soap_id.difference(mix_id):
	print c,soap_length[c]
	acc_length+=soap_length[c]
print acc_length


#parse the coord file 
al=parse_mummerFile(os.path.expanduser(coords_file))

coord_id=set()
for a in al:
	if "SOAP" in a['TAGQ']:
		coord_id.add(mix_line_re_3.findall(a['TAGQ'])[0])
	elif "SOAP" in a['TAGR']:
		coord_id.add(mix_line_re_3.findall(a['TAGR'])[0])


coord_id.difference(mix_id)

soap_id.difference(mix_id).intersection(soap_id.difference(coord_id))

# Any conting not in coord file but in mix? 
soap_id.difference(coord_id).intersection(mix_id)


for c in soap_id.difference(mix_id):
	c="A2_SOAP."+c
	print "*"*12,c
	implicated=[a for a in al if ((a['TAGQ']==c or a['TAGR']==c)and(a['TAGQ']!=a['TAGR']))]
	if len(implicated)>1:
		continue
	for a in implicated:
		print a
		print a['TAGQ'] in mix_id,a['TAGR'] in mix_id



for a in al:
	if a['TAGQ']=="A2_SOAP.scaffold17.8":
		if a['COVQ']>80:
			print a
	elif a['TAGQ']=="A2_SOAP.scaffold17.8":
		if a['COVR']>80:
			print a

