#! /usr/bin/python
import sys
from operator import itemgetter
import copy
import scipy
import networkx as nx
import collections
import integer_set as iset
import os
import re
#from decimal import Decimal, getcontext
from optparse import OptionParser 
#from types import *
from Bio.Seq import Seq
from Bio import SeqIO

from Bio.Alphabet import generic_dna

from IPython.core.debugger import Tracer; debug_here = Tracer()

import graph as graph 

import mummerParser

def reverse_is_in_alignments (aln, alignments) :
	new_aln = {}
	new_aln["LEN1"],new_aln["LEN2"],new_aln["IDY"],new_aln["LENR"],new_aln["LENQ"],new_aln["COVR"],new_aln["COVQ"],new_aln["TAGR"],new_aln["TAGQ"] = aln["LEN2"],aln["LEN1"],aln["IDY"],aln["LENQ"],aln["LENR"],aln["COVQ"],aln["COVR"],aln["TAGQ"],aln["TAGR"]
	ref_forward = (aln["S1"] < aln["E1"])
	query_forward = (aln["S2"] < aln["E2"])
	if ref_forward and query_forward : 
		new_aln["S1"], new_aln["E1"] = aln["S2"], aln["E2"]
		new_aln["S2"], new_aln["E2"] = aln["S1"], aln["E1"]
	else : # ref_forward and not query_forward : 
		new_aln["S1"], new_aln["E1"] = aln["E2"], aln["S2"]
		new_aln["S2"], new_aln["E2"] = aln["E1"], aln["S1"]
	result = False
	for A in alignments : 
		A_same_as_new_aln = False
		if A["TAGQ"] == aln["TAGR"] and A["TAGR"] == aln["TAGQ"] :
			A_same_as_new_aln = True
			for key in A.keys() :
				if A[key] != new_aln[key] :
					A_same_as_new_aln = False
		if A_same_as_new_aln : 
			result = True
	return result

def aln_pass_thresholds(aln,aln_threshold,ctg_threshold):
	idy_threshold = 75 #98
	cov_threshold = 99 #95
	# aln_threshold = 125 
	# ctg_threshold= 0
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
	if (aln["IDY"] > idy_threshold) and threshold_aln_OK and threshold_ctg_OK:
		return True
	else:
		return False
	# mummerParser.print_alignment(aln,sys.stderr,pretty=True)
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

def parse_alignments (file_adr, aln_threshold, ctg_threshold):
	##
	# @brief Read the COORD file 'file_adr' to fill and return a table containing all the alignments between two different contigs.
	# @param file_adr file containing the alignments between two assemblies (MUMmer COORD file)
	# @return a table containing all the alignments between two different contigs
	idy_threshold = 75 #98
	cov_threshold = 95 #95
	all_alignments = mummerParser.parse_mummerFile(file_adr)
	alignments=[]
	contigs = {}
	included_contigs = []
	for aln in all_alignments:
		# Build mapping
		for c in ["R", "Q"] :
			if aln["TAG"+c] not in contigs : 
				contigs[aln["TAG"+c]] = aln["LEN"+c]


		if aln['TAGQ']==aln['TAGR']:
			continue

		# We know that the two alignments are different, 
		# Do we have contigs completely included (AKA "covered") by another ? 

		if ((aln["COVR"] > cov_threshold) or ((aln["LENR"] - aln["LEN1"]) < ctg_threshold)):
			included_contigs.append(aln["TAGR"])
		elif ((aln["COVQ"] > cov_threshold) or ((aln["LENQ"] - aln["LEN2"]) < ctg_threshold)):
			included_contigs.append(aln["TAGQ"])

		if (aln_pass_thresholds(aln,aln_threshold,ctg_threshold)) and (not reverse_is_in_alignments(aln, alignments)) and (is_extremal(aln)):
			alignments.append(aln)
	#print "list of contigs "+str(contigs),"included:",included_contigs

	return all_alignments,alignments, contigs, included_contigs

### Output Assembly Sequences #############################################################################
def get_assembler_name (contig_name):
	if "_" in contig_name :
		return contig_name.split("_")[0].replace(">", "")
	else :
		return "None"

def parse_input_contigs_file (file_adr, contigs_list) : 
	contigs_sequences = {}
	contigs_file = open(file_adr, "r")
	line = contigs_file.readline()
	while line != "":
		contig_sequence=""
		if (re.match('>', line)):
			tab1=line.rstrip().split(">")
			tab2=tab1[1].split(" ")
			contig_id=tab2[0]
			line = contigs_file.readline()
		else :
			0/0
			sys.exit("Error : the file doesn't begin with a contig id like \">id\"\n"+line)
		if (re.match('^[A,T,G,C,N,a,t,c,g,n,\*,w,k,y,r,m,s,W,K,Y,R,M,S]+' , line)):
			while (re.match('^[A,T,G,C,N,a,t,c,g,n,\*,w,k,y,r,m,s,W,K,Y,R,M,S]+' , line)):
				contig_sequence = contig_sequence+line.upper().rstrip()
				line = contigs_file.readline()
		else :
			sys.exit("Error : no sequence for the contig id :"+contig_id+"\n"+line) 
		if contig_id in contigs_list : 
			contigs_sequences[contig_id]= contig_sequence
	contigs_file.close()
	return contigs_sequences

def write_sequence (output_file, contigs_sequences, path) :
	for c in path :
		# print "simple contig "+str(c["contig"])
		if c["sens"] == "f" :
			sub_sequence = contigs_sequences[c["contig"]][c["start"]:c["end"]]
		else : 
			sub_sequence = str(Seq(contigs_sequences[c["contig"]][c["start"]:c["end"]]).reverse_complement())
		output_file.write(sub_sequence)
	output_file.write("\n")

def compute_global_length_from_each_assembly(paths, single_contigs, contigs):
	nb_contigs_from_each_assembly = {}
	nb_part_of_contigs_from_each_assembly = {}
	length_of_parts_from_each_assembly = {}
	length_of_contigs_from_each_assembly = {}
	for p in paths : 
		for contig in p : 
			assembly = get_assembler_name(contig["contig"])
			length = contig["end"]-contig["start"]
			if assembly in length_of_parts_from_each_assembly.keys() : 
				length_of_parts_from_each_assembly[assembly] += length
				nb_part_of_contigs_from_each_assembly[assembly] += 1
			else : 
				length_of_parts_from_each_assembly[assembly] = length
				nb_part_of_contigs_from_each_assembly[assembly] = 1
	for k in length_of_parts_from_each_assembly.keys() : 
		nb_contigs_from_each_assembly[k] = 0
		length_of_contigs_from_each_assembly[k] = 0
	for s in single_contigs : 
		if isinstance(s, str) : 
			length = contigs[s]
			assembly = get_assembler_name(s)
		else : 
			contig_name = s[0]["contig"]
			length = contigs[contig_name]
			assembly = get_assembler_name(contig_name)
		if assembly in length_of_contigs_from_each_assembly.keys() : 
			length_of_contigs_from_each_assembly[assembly] += length
			nb_contigs_from_each_assembly[assembly] += 1
		else : 
			length_of_contigs_from_each_assembly[assembly] = length
			nb_contigs_from_each_assembly[assembly] = 1
	return nb_contigs_from_each_assembly, nb_part_of_contigs_from_each_assembly, length_of_parts_from_each_assembly, length_of_contigs_from_each_assembly

def select_single_contigs_of_one_assembly (paths, single_contigs, contigs) :
	nb_contigs_from_each_assembly, nb_part_of_contigs_from_each_assembly, length_of_parts_from_each_assembly, length_of_contigs_from_each_assembly = compute_global_length_from_each_assembly(paths, single_contigs, contigs)
	max_len = 0
	max_ass = None
	for a in length_of_contigs_from_each_assembly.keys() : 
		length_ass = length_of_contigs_from_each_assembly[a] + length_of_parts_from_each_assembly[a]
		if length_ass > max_len : 
			max_len = length_ass
			max_ass = a
	selected_single_contigs = []
	for s in single_contigs : 
		if get_assembler_name(s) == max_ass : 
			selected_single_contigs.append([{"contig":s,"start":0,"end":contigs[s],"sens":"f"}])
	return selected_single_contigs

def select_single_contigs_of_one_assembly_byL50 (paths, single_contigs, contigs) :
	### DEBUG SAM: If we use the L50 criterion, then it's the min L50 that we should choose, not the max as Florence coded!
	### DEBUGGED 
	if len(single_contigs)<1:
		return single_contigs
	sorted_contigs_length_lists_by_assembler = {}
	total_length_by_assembler = {}
	for s in single_contigs : 
		assembler = get_assembler_name(s)
		if assembler not in sorted_contigs_length_lists_by_assembler.keys() :
			sorted_contigs_length_lists_by_assembler[assembler] = [contigs[s]]
			total_length_by_assembler[assembler] = contigs[s]
		else : 
			sorted_contigs_length_lists_by_assembler[assembler].append(contigs[s])
			total_length_by_assembler[assembler] += contigs[s]
	L50_by_assembler = {}
	N50_by_assembler = {}
	for a in sorted_contigs_length_lists_by_assembler.keys() :
		sorted_contigs_length_lists_by_assembler[a].sort(reverse=True)
		L50_by_assembler[a] = 0
		N50_by_assembler[a] = 0
		for c in range(len(sorted_contigs_length_lists_by_assembler[a])) : 
			if N50_by_assembler[a] < total_length_by_assembler[a]*50/100 :
				N50_by_assembler[a] += sorted_contigs_length_lists_by_assembler[a][c]
				if N50_by_assembler[a] >= total_length_by_assembler[a]*50/100 : 
					L50_by_assembler[a] = c+1

	best_assembler,best_L50 = sorted(L50_by_assembler.items(),key=itemgetter(1))[0] # Select the assembler with min L50 (It's the number of contigs corresponding to the N50 Length, AKA N500 count)

	print "Selected assembler",best_assembler,"with L50 of",L50_by_assembler[best_assembler], "out of", L50_by_assembler
	selected_single_contigs = []
	for s in single_contigs : 
		if get_assembler_name(s) == best_assembler : 
			selected_single_contigs.append([{"contig":s,"start":0,"end":contigs[s],"sens":"f"}])
	return selected_single_contigs

def output_fasta (contigs_adr, paths, single_contigs, contigs, output_dir, alignments) : 
	print len(single_contigs)
	contigs_list = []
	for path in paths : 
		for n in path : 
			if n not in contigs_list :
				contigs_list.append(n["contig"])
	for s in single_contigs : 
		for n in s : 
			if n not in contigs_list :
				contigs_list.append(n["contig"])	
	contigs_sequences = parse_input_contigs_file(contigs_adr, contigs_list)
	output_file_path = open(output_dir+"/Mix_assembly_path.fasta", "w")
	output_file_ctg = open(output_dir+"/Mix_assembly_ctg.fasta", "w")
	output_file_tot = open(output_dir+"/Mix_assembly.fasta", "w")
	p = 0
	for p in range (len(paths)) : 
		if len(paths[p])==0:
			continue
		output_file_path.write(">"+str(p+1)+str(paths[p])+"\n")
		write_sequence (output_file_path, contigs_sequences, paths[p])
		output_file_tot.write(">"+str(p+1)+str(paths[p])+"\n")
		write_sequence (output_file_tot, contigs_sequences, paths[p])
	output_tab = open(output_dir+"/Mix_selected_isolated_contigs.csv", "w")
	for c in range (len(single_contigs)):
		if len(single_contigs[c])==0:
			continue
		output_file_ctg.write(">"+str(p+1+c+1)+"_"+str(single_contigs[c][0]["contig"])+"\n")
		write_sequence (output_file_ctg, contigs_sequences, single_contigs[c])
		output_file_tot.write(">"+str(p+1+c+1)+"_"+str(single_contigs[c][0]["contig"])+"\n")
		write_sequence (output_file_tot, contigs_sequences, single_contigs[c])
		output_table(c, alignments, single_contigs, paths, output_tab, str(p+1+c+1)) 
	output_file_path.close()
	output_file_ctg.close()
	output_file_tot.close()
	output_tab.close()

def output_table(c, alignments, single_contigs, paths, output_tab, nb_ctg) : 
	for aln in alignments : 
		aln_found = False
		input_contig_id = single_contigs[c][0]["contig"]
		if input_contig_id == aln["TAGR"] : 
			aln_found = True
			per_cov, ctg_len = aln["COVR"], aln["LENR"]
			contig_id_in_extension = aln["TAGQ"]
			output_extension_id = find_contig_in_paths(paths, contig_id_in_extension)
		elif input_contig_id == aln["TAGQ"] :
			aln_found = True
			per_cov, ctg_len = aln["COVQ"], aln["LENQ"]
			contig_id_in_extension = aln["TAGR"]
			output_extension_id = find_contig_in_paths(paths, contig_id_in_extension)
		if aln_found and output_extension_id :
			output_tab.write(nb_ctg+"\t"+input_contig_id+"\t|\t"+contig_id_in_extension+"\t"+str(output_extension_id)+"\t|\t"+str(per_cov)+"\t"+str(ctg_len)+"\n")

def find_contig_in_paths(paths, ctg_id) : 
	for p in range(len(paths)) : 
		for c in range(len(paths[p])) :
			if paths[p][c]["contig"] == ctg_id : 
				return p
	return None			

def print_stats(paths, single_contigs, contigs):
	nb_contigs_from_each_assembly, nb_part_of_contigs_from_each_assembly, length_of_parts_from_each_assembly, length_of_contigs_from_each_assembly = compute_global_length_from_each_assembly(paths, single_contigs, contigs)
	for a in set(nb_contigs_from_each_assembly.keys()+nb_part_of_contigs_from_each_assembly.keys()) : 
		print a
		if a in nb_contigs_from_each_assembly.keys() : 
			print "\t"+str(nb_contigs_from_each_assembly[a])+" contigs entirely kept"
		else : 
			print "\t"+str(0)+" contigs entirely kept"
		if a in nb_part_of_contigs_from_each_assembly.keys() : 
			print "\t"+str(nb_part_of_contigs_from_each_assembly[a])+" contigs in extensions"
		else : 
			print "\t"+str(0)+" contigs in extensions"
		if a not in length_of_parts_from_each_assembly.keys() : 
			print "\t"+str(0+length_of_contigs_from_each_assembly[a])+" bp"
		elif a not in length_of_contigs_from_each_assembly.keys() : 
			print "\t"+str(length_of_parts_from_each_assembly[a]+0)+" bp"
		else : 
			print "\t"+str(length_of_parts_from_each_assembly[a]+length_of_contigs_from_each_assembly[a])+" bp"
		
# def remove_included_contigs_in_paths_from_selected_contigs(paths, selected_single_contigs, alignments, contig_threshold) :
# 	cov_threshold = 50 #Debug macha
# 	contigs_in_paths = []
# 	for p in paths : 
# 		for n in p : 
# 			contigs_in_paths.append(n["contig"])
# 	contigs_included = []
# 	contigs_not_included = list(selected_single_contigs)
# 	for s in selected_single_contigs : 
# 		contig_name = s[0]["contig"]
# 		if is_contig_aligned_with_a_path(cov_threshold, contig_threshold, contig_name, paths, alignments) : 
# 				contigs_included.append(s)
# 	for i in contigs_included : 
# 		if i in contigs_not_included :
# 			contigs_not_included.remove(i)
# 	return contigs_not_included

# def is_contig_aligned_with_a_path(cov_threshold, contig_threshold, contig, paths, alignments):
# 	cov_sum_by_ctg = {}
# 	# debug_here()
# 	# if contig=="A2_SOAP.scaffold8.1":
# 	# 	debug_here()
# 	for a in alignments :
# 		# if ("A2_SOAP.scaffold8.1" == a['TAGR']) or ("A2_SOAP.scaffold8.1" == a['TAGQ']):
# 		# 	debug_here()
# 		if (a["TAGR"] == contig) and (find_contig_in_paths(paths, a["TAGQ"])!=None): 
# 			if ((a["COVR"] >= cov_threshold) or ((a["LENR"] - a["LEN1"]) < contig_threshold*2)) :
# 				return True
# 			elif a["TAGQ"] in cov_sum_by_ctg.keys() : 
# 				cov_sum_by_ctg[a["TAGQ"]] += float(a["COVR"])
# 			else : 
# 				cov_sum_by_ctg[a["TAGQ"]] = float(a["COVR"])
# 		elif (a["TAGQ"] == contig) and (find_contig_in_paths(paths, a["TAGR"])!=None): 
# 			if ((a["COVQ"] >= cov_threshold) or ((a["LENQ"] - a["LEN2"]) < contig_threshold*2)) :
# 				return True
# 			elif a["TAGR"] in cov_sum_by_ctg.keys() : 
# 				cov_sum_by_ctg[a["TAGR"]] += float(a["COVQ"])
# 			else : 
# 				cov_sum_by_ctg[a["TAGR"]] = float(a["COVQ"])
# 	# debug_here()
# 	print cov_sum_by_ctg
# 	for path_nb in range(len(paths)) :
# 		for p in paths[path_nb] : 
# 			p_cov_sum = 0
# 			p_id = p["contig"]
# 			if p_id in cov_sum_by_ctg.keys() : 
# 				p_cov_sum += cov_sum_by_ctg[p_id]
# 			if (p_cov_sum > cov_threshold) : 
# 				return True
# 	return False


			


def lower_duplication(paths,selected_single_contigs,cv_thr=0.65):
	# debug_here()
	all_elements=paths+selected_single_contigs
	# all_elements=paths
	scores=scipy.zeros((len(all_elements),len(all_elements)))

	for i in range(len(all_elements)):
		for j in range(i+1,len(all_elements)):
			p1=all_elements[i]
			p2=all_elements[j]
			cv=coverage(alignments,p1,p2)
			# if cv!=(0,0):
			# 	print p1,p2,cv
			# print i,j,cv
			if cv[0]>cv_thr and cv[1]>cv_thr: #competition 
				if cv[0]>cv[1]: #p1 is more covered by p2, we add an "inclusion" edge from p1 to p2 
					scores[i,j]=int(cv[0]*100)
				else:
					scores[j,i]=int(cv[0]*100) # we do it from p2 to p1 
			elif cv[0]>cv_thr:
				scores[i,j]=int(cv[0]*100)
			elif cv[1]>cv_thr:
				scores[j,i]=int(cv[1]*100)

	D=nx.DiGraph(scores)
	# Detect elements to remove 
	somethingChanged=True
	while (somethingChanged):
		somethingChanged=False
		# print "New iter"
		for (n,d_in),(n,d_out) in zip(D.in_degree_iter(),D.out_degree_iter()):
			if (d_in==0) and (d_out >0) :
				src_ctg=[x['contig'] for x in all_elements[n]]
				src_assemblers=[x[:2] for x in src_ctg]
				# print n,src_assemblers,src_ctg
				if len(src_ctg)>1: #is a path, either covered by a path or by a contig
					# print "is path"
					# print "removing",src_ctg,"covered by",D[n],"which is",[all_elements[x] for x in D[n]]
					paths[n]=[]
					D.remove_node(n)
					somethingChanged=True
				else:
					# print "is single"
					tgt_assemblers=set()
					covered_by_a_path=False
					for tgt in D[n]:
						tgt_ctg=[x['contig'] for x in all_elements[tgt]]
						if len(tgt_ctg)>1:
							covered_by_a_path=True
						tgt_assemblers.update([x[:2] for x in tgt_ctg])
						# print "\t"*2,tgt_assemblers,tgt_ctg,D[n][tgt]
					# if  (not covered_by_a_path):
					# 	continue
					# if "A2" in tgt_assemblers:
					# 	continue
					# print "removing",src_ctg,"covered by",D[n],"which is",[all_elements[x] for x in D[n]]
					selected_single_contigs[n-len(paths)]=[]
					D.remove_node(n)
					somethingChanged=True
	return paths,selected_single_contigs




### Global processing #############################################################################
def processing (all_alignments,alignments, contigs_file, output_base, contig_threshold, aln_threshold, display_dot, display_cytoscape, contigs, included_contigs):
	##
	# @brief Construct an assembly from the alignments, and then make a graph with it to extend the contigs and write them in a new file. 
	# @param alignments		list containing all the alignments between tow different contigs
	# @param contigs_files		list of file(s) containing the contigs sequences
	# @param output_repository	repository address where the results will be written
	# @param contig_threshold	minimum contig length to consider
	# @param aln_threshold		minimum length alignment to consider
	# @param pattern			pattern to use to distinguish the contigs from the two assemblers
	
	global graph_assembly,single_contigs,selected_single_contigs,selected_single_contigs_not_included,selected_paths,selected_contigs,paths
	output_dir = output_base+"/Mix_results_A"+str(aln_threshold)+"_C"+str(contig_threshold)
	try:
		os.mkdir(output_base)
	except:
		pass

	try:
		os.mkdir(output_dir)
	except OSError:
		pass

	print "included?",("A2_scf1_33" in included_contigs)

	graph_assembly = graph.graph(alignments, contigs, included_contigs)
	initial_assembly_graph=copy.deepcopy(graph_assembly)
	single_contigs = graph_assembly.get_single_contigs(contigs)

	if display_dot : 
		graph_assembly.write_dot_graph(output_dir+"/graph1.dot")

	# Make the extensions
	paths = graph_assembly.select_extensions(contigs)
#	debug_here()

	if display_dot : 
		graph_assembly.write_dot_graph(output_dir+"/graph2.dot")
	if display_cytoscape : 
		graph_assembly.write_cytoscape_graph(output_dir+"/")

	#selected_single_contigs = select_single_contigs_of_one_assembly(paths, single_contigs, contigs)
	selected_single_contigs = select_single_contigs_of_one_assembly_byL50(paths, single_contigs, contigs)
	print "selected_single_contigs"+str(selected_single_contigs)
	print len(selected_single_contigs)
	#selected_single_contigs_not_included = remove_included_contigs_from_selected_contigs(selected_single_contigs, alignments)
	# selected_single_contigs_not_included = remove_included_contigs_in_paths_from_selected_contigs(paths, selected_single_contigs, alignments, contig_threshold)
	# selected_paths=paths
	# selected_paths = remove_paths_included_in_selected_contigs(paths, contigs, alignments)

	# Todo: Real optimization to select the assembly
	#Remove paths covered by selected contig 
	# Todo: Iterative computation, removing one element at a time 
	# Todo: Maybe consider the problem of coverage by the same assembly: should be kept anyhow?
	# Might be order dependant! 
	index_alignment(all_alignments)
	# cv_thr=0.50
	# selected_paths=[]
	# for p in paths:
	# 	covered=False
	# 	print p
	# 	for c in selected_single_contigs:
	# 		cv1,cv2=coverage(all_alignments,p,c)
	# 		if cv1>=cv_thr : 
	# 			print "should remove",p,"covered at",cv1,"by",c
	# 			covered=True
	# 			break
	# 		# if cv2>=cv_thr:
	# 		# 	print "should remove",c,"covered at",cv2,"by",p
	# 	if not covered:
	# 		selected_paths.append(p)

	# selected_contigs=[]
	# for c in selected_single_contigs:
	# 	print c
	# 	covered=False
	# 	for p in selected_paths:
	# 		cv1,cv2=coverage(all_alignments,c,p)
	# 		if cv1>=cv_thr : 
	# 			print "should remove",c,"covered at",cv1,"by",p
	# 			covered=True
	# 			break
	# 	if not covered:
	# 		selected_contigs.append(c)
	selected_paths,selected_contigs=lower_duplication(paths,selected_single_contigs)
	print "After duplication minimisation"
	print "selected_single_contigs"+str(selected_single_contigs)
	print len(selected_single_contigs)

	save_statistics_to_files(\
		initial_assembly_graph=initial_assembly_graph,\
		assembly_graph=graph_assembly,\
		all_alignments=all_alignments,\
		selected_alignements=alignments,\
		all_paths=paths,\
		selected_paths=selected_paths,\
		all_contigs=contigs,\
		included_contigs=included_contigs,\
		single_contigs=single_contigs,\
		selected_single_contigs=selected_single_contigs,\
		selected_contigs=selected_contigs,\
		output_dir=output_dir)

	output_fasta(contigs_file, selected_paths, selected_contigs, contigs, output_dir, alignments)
	#output_fasta(contigs_file, paths, selected_single_contigs, contigs, output_dir)
	# print_stats(selected_paths, selected_contigs, contigs)
	#print_stats(paths, selected_single_contigs, contigs)

def save_statistics_to_files(initial_assembly_graph,assembly_graph,all_alignments,selected_alignements,all_paths, selected_paths, all_contigs,included_contigs, selected_contigs,selected_single_contigs,single_contigs,output_dir):
	if len(alignments)<1:
		return
	global aln_threshold,contig_threshold
	nx.write_gml(initial_assembly_graph.GRAPH,output_dir+"/initial_assembly_graph.gml")
	nx.write_gml(assembly_graph.GRAPH,output_dir+"/reduced_assembly_graph.gml")
	all_alignments_file = open(output_dir+"/all_alignments.csv","w")
	a=alignments[0]
	header = sorted(a.keys())+["is_selected","is_extremal","pass_thr"]
	print >>all_alignments_file,",".join(header)
	for a in alignments:
		is_selected=("is_selected",a in selected_alignements)
		is_extremal_a = ("is_extremal",is_extremal(a))
		pass_thr = ("pass_thr",aln_pass_thresholds(a,aln_threshold,contig_threshold))
		items = a.items()
		items.sort(key=itemgetter(0))
		items +=[is_selected,is_extremal_a,pass_thr]
		print>>all_alignments_file, ",".join(map(str,[x[1] for x in items]))

	all_contigs_file = open(output_dir+"/all_contigs.csv","w")
	# Build a tally 
	contig_tally=collections.defaultdict(int)
	for a in selected_alignements:
		contig_tally[a['TAGQ']]+=1
		contig_tally[a['TAGR']]+=1

	selected_contigs_c=set()
	for xl in selected_contigs:
		selected_contigs_c.update([x['contig'] for x in xl])
	selected_single_contigs_c=set()
	for xl in selected_single_contigs:
		selected_single_contigs_c.update([x['contig'] for x in xl])

	print >>all_contigs_file,",".join(("ID","length","assembler","is_included","is_selected","is_selected_single","is_single","n_alignments","in_path","in_selected_path","in_path_idx","idx_in_selected_path","in_assembly_graph","in_final_assembly_graph"))
	for c in all_contigs.keys():
		in_path=False
		in_selected_path=False
		in_path_idx=-1
		idx_in_selected_path=-1
		for p in all_paths:
			if c in p:
				in_path=True
				break
		for i in range(len(selected_paths)):
			p=selected_paths[i]
			p.sort(key=lambda x:x['start'])
			if c in [x['contig'] for x in p]:
				in_selected_path=True
				in_path_idx=i
				idx_in_selected_path=[x['contig'] for x in p].index(c)
				break

		in_assembly_graph=len([x for x in initial_assembly_graph.GRAPH.nodes(data=True) if x[1]['contig']==c])
		in_final_assembly_graph=len([x for x in assembly_graph.GRAPH.nodes(data=True) if x[1]['contig']==c])
		items=(c,all_contigs[c],get_assembler_name(c),c in included_contigs, c in selected_contigs_c,c in selected_single_contigs_c,c in single_contigs,contig_tally[c],in_path,in_selected_path,in_path_idx,idx_in_selected_path,in_selected_path,in_final_assembly_graph)
		print >>all_contigs_file,",".join(map(str,items))
	all_alignments_file.close()
	all_contigs_file.close()


def main():
	global alignments,contigs,included_contigs,aln_threshold,contig_threshold
	parser = OptionParser()
	parser.add_option("-a", "--aln", dest="aln", help ="the file containing the alignments (.coords)")
	parser.add_option("-o", "--out", dest="out", help ="the output directory where the scaffolds will be written (it must already exist)") 
	parser.add_option("-c", "--ctg", dest="ctg", help ="the file containing all the contigs that were used in the alignment") 
	parser.add_option("-A", "--ath", dest="ath", help ="minimum length of alignment (optionnal)")
	parser.add_option("-C", "--cth", dest="cth", help ="minimum length of contig (optionnal)")
	parser.add_option("-d", "--dot", dest="dot", help ="write the graphs in dot format", default=False, action='store_true')
	parser.add_option("-g", "--graph", dest="graph", help ="write the graphs in cytoscape format", default=False, action='store_true')
	parser.add_option("-r", "--restrict-to-aligned", dest="restrict_aligned", help ="If on, restrict output to aligned coords",default=False,action="store_true")
	(options, args) = parser.parse_args()
	aln_adr = options.aln
	ctg_adr = options.ctg
	output_repository = options.out
	ath = options.ath
	cth = options.cth
	display_dot = options.dot
	display_cytoscape = options.graph
	restrict_aligned=options.restrict_aligned

	if cth :
		contig_threshold = int(cth)
		if ath :
			aln_threshold = int(ath)
		else :
			aln_threshold = 0
	else :
		if ath :
			aln_threshold = int(ath)
			contig_threshold = 0

	alignments,included_alignemnts, contigs, included_contigs = parse_alignments(aln_adr, aln_threshold, contig_threshold)
	cleaned_alignments=[x for x in included_alignemnts if ((x['TAGQ'] not in included_contigs)and(x['TAGR'] not in included_contigs))]

	# We check whether we have contigs with too many alignements. E.g. One of the contigs of CLC consisted of reads that are mapped by MIRA onto k contigs. Once aligned, we have a CLC contig aligned to k MIRA contigs, without any inclusion 
	contig_tally=collections.defaultdict(int)
	for a in cleaned_alignments:
		contig_tally[a['TAGQ']]+=1
		contig_tally[a['TAGR']]+=1
	mu=scipy.average(contig_tally.values())
	sigma=scipy.std(contig_tally.values())
	extreme_contigs = [k for k,v in contig_tally.items() if v>=3*mu]
	print "Removed contigs with too many alignments:",[contig_tally[k] for k in extreme_contigs]
	cleaned_alignments=[x for x in cleaned_alignments if ((x['TAGQ'] not in extreme_contigs)and(x['TAGR'] not in extreme_contigs))]
	included_contigs=list(set(included_contigs).union(extreme_contigs))




	# Add contigs without any alignement 
	# contigs['A2_SOAP.C769.1']=12345
	# included_contigs.append('A2_SOAP.C769.1')

	if not restrict_aligned:
		handle = open(ctg_adr, "rU")
		for record in SeqIO.parse(handle, "fasta") :
			if (record.id in included_contigs) or (record.id in contigs):
				continue
			found=False
			for a in alignments:
				if (a['TAGR']==record.id) or (a['TAGQ']==record.id):
					found=True
					break
			if found:
				continue
			contigs[record.id]=len(record.seq)
			included_contigs.append(record.id)
			print "DEBUG SAM Added conting ID",record.id
		handle.close()

	processing(alignments,cleaned_alignments, ctg_adr, output_repository, contig_threshold, aln_threshold, display_dot, display_cytoscape, contigs, included_contigs)


def assemble_path_as_intset(path):
	#Determine for each path the offset to apply to each of its contig for the weaving 
	# Build the intset indicating the coord system of the final assembly
	offsets=[]
	running_offset=0
	for contig in path:
		offsets.append(running_offset)
		running_offset+=contig['end']
	path_assembly=iset.IntSet(0)
	for i in range(len(path)):
		contig= path[i]
		coords=sorted((contig['start'],contig['end']))
		# offset 
		coords[1] = coords[1] + (offsets[i] - coords[0])
		coords[0] = offsets[i]
		path_assembly |= iset.IntSet(tuple(coords))
	return path_assembly


def index_alignment(all_alignments):
	global indexed_alignments 
	indexed_alignments={}
	for a in all_alignments:
		if a['TAGQ'] not in indexed_alignments:
			indexed_alignments[a['TAGQ']]=collections.defaultdict(list)
		if a['TAGR'] not in indexed_alignments:
			indexed_alignments[a['TAGR']]=collections.defaultdict(list)
		indexed_alignments[a['TAGQ']][a['TAGR']].append(a)
		indexed_alignments[a['TAGR']][a['TAGQ']].append(a)

def test_indexation():
	import random
	for c1 in random.sample(contigs.keys(),k=10):
		for c2 in random.sample(contigs.keys(),k=10):
			al=[a for a in alignments if set([a['TAGQ'],a['TAGR']])==set([c1,c2])]
			print al

			print indexed_alignments.get(c1,{}).get(c2,[])
			print indexed_alignments.get(c2,{}).get(c1,[])
			print "--"*12
			assert len(al)==len(indexed_alignments.get(c1,{}).get(c2,[]))


def coverage(all_alignments,contigs1,contigs2,verbose=False):

	#contigs1,contigs2 are assumed to be paths, possibly of lenght 1
	contigs1_covered=iset.IntSet()
	contigs2_covered=iset.IntSet()
	c1_c2_coverage={}
	for c1 in contigs1:
		for c2 in contigs2:
			# get alignment between them 
			# al=[a for a in all_alignments if set([a['TAGQ'],a['TAGR']])==set([c1['contig'],c2['contig']])]
			al=indexed_alignments.get(c1['contig'],{}).get(c2['contig'],[])
			if len(al)<1:
				continue
			if c1['contig'][:2]==c2['contig']:#Same assembler
				continue
			#build two iset.intsets to compute the coverages 
			coverage_c1_in_c2=iset.IntSet()
			coverage_c2_in_c1=iset.IntSet()
			if verbose:
				print len(al),"alignments found"
			for a in al:
				if a['TAGR']==c1['contig']:
					coords=sorted((a['S1'],a['E1']))
					coverage_c1_in_c2=coverage_c1_in_c2 | iset.IntSet(tuple(coords))

					coords=sorted((a['S2'],a['E2']))
					coverage_c2_in_c1=coverage_c2_in_c1 | iset.IntSet(tuple(coords))

				else:
					coords=sorted((a['S2'],a['E2']))
					coverage_c1_in_c2=coverage_c1_in_c2 | iset.IntSet(tuple(coords))

					coords=sorted((a['S1'],a['E1']))
					coverage_c2_in_c1=coverage_c2_in_c1 | iset.IntSet(tuple(coords))
			if verbose:
				print c1,c2
				print coverage_c1_in_c2,len(coverage_c1_in_c2)*1.0/c1['end']
				print coverage_c2_in_c1,len(coverage_c2_in_c1)*1.0/c2['end']
				print "--"
			c1_c2_coverage[(c1['contig'],c2['contig'])]=coverage_c1_in_c2
			c1_c2_coverage[(c2['contig'],c1['contig'])]=coverage_c2_in_c1
	if len(c1_c2_coverage)==0:
		return (0,0)

	# compute the running offsets 
	offsets_c1=[]
	running_offset=0
	for contig in contigs1:
		offsets_c1.append(running_offset)
		running_offset+=contig['end']
	offsets_c2=[]
	running_offset=0
	for contig in contigs2:
		offsets_c2.append(running_offset)
		running_offset+=contig['end']

	#Update the coverage coordinates 
	contigs1_coverage=iset.IntSet()
	contigs2_coverage=iset.IntSet()

	for c1,c2 in c1_c2_coverage.keys():
		if c1 in [c['contig'] for c in contigs1]:
			c1i = [c['contig'] for c in contigs1].index(c1)
			# c2i = [c['contig'] for x in contigs2].index(c2)
			c1lb=contigs1[c1i]['start']
			c1o=offsets_c1[c1i]-c1lb
			# c2o=offsets_c2[c2i]
			# c2lb=contigs1[c2i]['start']
			if verbose:
				print c1,c2,c1o,c1lb,c1_c2_coverage[(c1,c2)].shift(c1o).trim(lbound=c1lb)
			contigs1_coverage|=c1_c2_coverage[(c1,c2)].shift(c1o).trim(lbound=c1lb)
		else:
			c1i = [c['contig'] for c in contigs2].index(c1)
			# c2i = [c['contig'] for x in contigs2].index(c2)
			c1lb=contigs2[c1i]['start']
			c1o=offsets_c2[c1i]-c1lb
			# c2o=offsets_c2[c2i]
			# c2lb=contigs1[c2i]['start']
			if verbose:
				print c1,c2,c1o,c1lb,c1_c2_coverage[(c1,c2)].shift(c1o).trim(lbound=c1lb)
			contigs2_coverage|=c1_c2_coverage[(c1,c2)].shift(c1o).trim(lbound=c1lb)
	# debug_here()
	contigs1_assembly=assemble_path_as_intset(contigs1)
	contigs2_assembly=assemble_path_as_intset(contigs2)
	cov_of_c1=len(contigs1_coverage)*1.0/len(contigs1_assembly)
	cov_of_c2=len(contigs2_coverage)*1.0/len(contigs2_assembly)
	if verbose:
		print contigs1_assembly,contigs1_coverage,cov_of_c1
		print contigs2_assembly,contigs2_coverage,cov_of_c2
	return (cov_of_c1,cov_of_c2)

if __name__ == "__main__":
	main()
# sys.exit(0)



# #Make a competition between paths 
# import scipy
# import networkx as nx
# scores=scipy.zeros((len(selected_paths),len(selected_paths)))
# for i in range(len(selected_paths)):
# 	for j in range(i+1,len(selected_paths)):
# 		p1=selected_paths[i]
# 		p2=selected_paths[j]
# 		cv=coverage(alignments,p1,p2)
# 		print i,j,cv
# 		if cv[0]>0.5 and cv[1]>0.5: #competition 
# 			if cv[0]>cv[1]: #p1 is more covered by p2, we add an "inclusion" edge from p1 to p2 
# 				scores[i,j]=int(cv[0]*100)
# 			else:
# 				scores[j,i]=int(cv[0]*100) # we do it from p2 to p1 
# 		elif cv[0]>0.5:
# 			scores[i,j]=int(cv[0]*100)
# 		elif cv[1]>0.5:
# 			scores[j,i]=int(cv[1]*100)
# D=nx.DiGraph(scores)


# scores=scipy.zeros((len(selected_contigs),len(selected_contigs)))
# for i in range(len(selected_contigs)):
# 	for j in range(i+1,len(selected_contigs)):
# 		p1=selected_contigs[i]
# 		p2=selected_contigs[j]
# 		cv=coverage(alignments,p1,p2)
# 		print i,j,cv
# 		if cv[0]>0.5 and cv[1]>0.5: #competition 
# 			if cv[0]>cv[1]: #p1 is more covered by p2, we add an "inclusion" edge from p1 to p2 
# 				scores[i,j]=int(cv[0]*100)
# 			else:
# 				scores[j,i]=int(cv[0]*100) # we do it from p2 to p1 
# 		elif cv[0]>0.5:
# 			scores[i,j]=int(cv[0]*100)
# 		elif cv[1]>0.5:
# 			scores[j,i]=int(cv[1]*100)
# D=nx.DiGraph(scores)
# # Detect contigs to remove 
# for (n,d_in),(n,d_out) in zip(D.in_degree_iter(),D.out_degree_iter()):
# 	if d_in==0 and d_out >0 : 
# 		selected_contigs[n]=[]


# Both combined ?
# all_elements=selected_paths
# all_elements=paths
# all_elements=selected_paths+selected_contigs

# f=open("rhodo_state.cPickle","r")
# selected_paths,selected_contigs,contigs,paths,selected_single_contigs=cPickle.load(f)
# f.close()



# # Todo: remove node then reiterate? 

# # f=open("rhodo_state.cPickle","w")
# # cPickle.dump((selected_paths,selected_contigs,contigs,paths,selected_single_contigs),f)
# # f.close()

# # f=open("aerius_state.cPickle","w")
# # cPickle.dump((selected_paths,selected_contigs,contigs),f)
# # f.close()

# # Manual edition for aerius 
# # selected_paths[5]=[] #high coverage with path[6]
# # selected_paths[8]=[] #high coverage with path[7]
# # selected_paths[9]=[] #high coverage with path[7]
# # selected_paths[10]=[] #high coverage with path[7]
# # selected_paths[11]=[] #high coverage with path[6]
# # output_fasta("../aerius_AP_SOAP_BB.fasta",selected_paths,selected_contigs,contigs,"manual",alignments)
# # f=open("aerius_state.cPickle","r")
# # selected_paths,selected_contigs,contigs=cPickle.load(f)
# # f.close()

# # Manual edition for rhodo
# output_fasta("../contigs.AllpathSoapBambus.fa",selected_paths,selected_contigs,contigs,"manualRho",alignments)
# output_fasta("../contigs.AllpathSoapBambus.fa",paths,selected_contigs,contigs,"manualRho",alignments)
# output_fasta("../contigs.AllpathSoapBambus.fa",paths,selected_single_contigs,contigs,"manualRho",alignments)
# f=open("rhodo_state.cPickle","r")
# selected_paths,selected_contigs,contigs,paths,selected_single_contigs=cPickle.load(f)
# f.close()

