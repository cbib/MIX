import collections
assemblies=[x.strip() for x in open("datasets/GAGE-B/assemblies_index.txt")]

by_specie=collections.defaultdict(list)
for a in assemblies:
	specie,assembly=a.split("/")
	by_specie[specie].append(assembly)



references={
'A_hydrophila_HiSeq' : "Aeromonas_hydrophila",
'B_cereus_MiSeq' : "Bacillus_cereus",
'B_fragilis_HiSeq' : "Bacteroides_fragilis",
'M_abscessus_HiSeq' : "Mycobacterium_abscessus",
'R_sphaeroides_HiSeq' : "Rhodobacter_sphaeroides",
'S_aureus_HiSeq' : "Staphylococcus_aureus",
'V_cholerae_HiSeq' : "Vibrio_cholerae",
'X_axonopodis_HiSeq' : "Xanthomonas_axonopodis",
 }


dependencies="""temp_assemblies/%(tag)s.fasta: %(files)s"""
quasts="""result_statistics/pairwise-gageb/%(specie)s_quast: %(all_tags)s
	rm -rf result_statistics/pairwise-gageb/%(specie)s_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/%(ref)s/%(ref)s_ref.fasta -G datasets/reference/%(ref)s/%(ref)s_ref.gff  $^ """


all_targets=[]
for specie,available_assemblies in by_specie.items():
	if specie not in references:
		continue
	all_tags_for_specie=[]
	for i in range(len(available_assemblies)):
		for j in range(i+1,len(available_assemblies)):
			asbl1=available_assemblies[i]
			asbl2=available_assemblies[j]
			tag=specie+"_"+available_assemblies[i].split("_")[0]+"-"+available_assemblies[j].split("_")[0]
			print dependencies%{"tag":tag,"files":" ".join(["datasets/GAGE-B/"+specie+"/"+asbl1,"datasets/GAGE-B/"+specie+"/"+asbl2])}
			all_tags_for_specie.append(tag+"_mix")
	print quasts%{'specie':specie,\
					"all_tags":" ".join(["result_assemblies/"+x+".fasta" for x in all_tags_for_specie]),\
					"ref":references[specie]}
	print '\n'
	all_targets.append("result_statistics/pairwise-gageb/%s_quast"%(specie))
print "AllMixGAGEB:%s"%(" ".join(all_targets))
