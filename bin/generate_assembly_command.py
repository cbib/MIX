import os
targets=["MOVI", "MMC", "MSCe", "MSCd", "MSCc", "MSCb"]
template="""
%(mollicute)s_%(assemblies)s: 
	./preprocessing.py ../%(mollicute)s/%(assemblyAfname)s ../%(mollicute)s/%(assemblyBfname)s -o ../%(mollicute)s/%(mollicute)s_%(assemblies)s.fasta
	../MUMmer3.23/nucmer -p "%(mollicute)s_%(assemblies)s" --maxmatch -l 30 -banded ../%(mollicute)s/%(mollicute)s_%(assemblies)s.fasta ../%(mollicute)s/%(mollicute)s_%(assemblies)s.fasta 2>/dev/null
	../MUMmer3.23/show-coords -l -c %(mollicute)s_%(assemblies)s.delta > ../%(mollicute)s/%(mollicute)s_%(assemblies)s.coords
	rm %(mollicute)s_%(assemblies)s.delta
	$(pyinterp) Mix.py -A 125 -C 0 -o %(mollicute)s_%(assemblies)s  -a ../%(mollicute)s/%(mollicute)s_%(assemblies)s.coords -c ../%(mollicute)s/%(mollicute)s_%(assemblies)s.fasta
	mv %(mollicute)s_%(assemblies)s/Mix_results_A125_C0/Mix_assembly.fasta ../%(mollicute)s/MIX_%(mollicute)s_%(assemblies)s.fasta
"""

template="""temp_assemblies/%(mollicute)s_%(assemblies)s.fasta: datasets/Mollicutes/%(mollicute)s/%(assemblyAfname)s datasets/Mollicutes/%(mollicute)s/%(assemblyBfname)s"""

all_targets=[]
for mollicute in targets:
	these_targets=[]
	avail_assemblies=sorted([x for x in os.listdir("datasets/Mollicutes/"+mollicute) if "fa" in os.path.splitext(x)[1]])
	# print mollicute
	# print "#",avail_assemblies

	for i in range(len(avail_assemblies)):
		a1=avail_assemblies[i]
		for j in range(i+1,len(avail_assemblies)):
			a2=avail_assemblies[j]
			# print "#",i,j,a1,a2
			assemblies=""
			if "step2" in a1:
				assemblies+="MIRA"
			elif "CLC" in a1:
				assemblies+="CLC"
			elif "ABySS" in a1:
				assemblies+="AB"
			if "step2" in a2:
				assemblies+="_MIRA"
			elif "CLC" in a2:
				assemblies+="_CLC"
			elif "ABySS" in a2:
				assemblies+="_AB"
			# if assemblies!="CLC_MIRA":
			# 	continue
			print template%{'mollicute':mollicute,'assemblies':assemblies,'assemblyAfname':a1,'assemblyBfname':a2}
			all_targets.append(mollicute+"_"+assemblies)
			these_targets.append(mollicute+"_"+assemblies)
	# print "%s:%s"%(mollicute," ".join(these_targets))
# print "Allmolli:%s"%(" ".join(all_targets))
