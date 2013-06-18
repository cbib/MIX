import os
targets=["MOVI", "MMC", "MSCe", "MSCd", "MSCc", "MSCb"]
template="""
%(molicute)s_%(assemblies)s: 
	./preprocessing.py ../%(molicute)s/%(assemblyAfname)s ../%(molicute)s/%(assemblyBfname)s -o ../%(molicute)s/%(molicute)s_%(assemblies)s.fasta
	../MUMmer3.23/nucmer -p "%(molicute)s_%(assemblies)s" --maxmatch -l 30 -banded ../%(molicute)s/%(molicute)s_%(assemblies)s.fasta ../%(molicute)s/%(molicute)s_%(assemblies)s.fasta 2>/dev/null
	../MUMmer3.23/show-coords -l -c %(molicute)s_%(assemblies)s.delta > ../%(molicute)s/%(molicute)s_%(assemblies)s.coords
	rm %(molicute)s_%(assemblies)s.delta
	$(pyinterp) Mix.py -A 125 -C 0 -o %(molicute)s_%(assemblies)s  -a ../%(molicute)s/%(molicute)s_%(assemblies)s.coords -c ../%(molicute)s/%(molicute)s_%(assemblies)s.fasta
	mv %(molicute)s_%(assemblies)s/Mix_results_A125_C0/Mix_assembly.fasta ../%(molicute)s/MIX_%(molicute)s_%(assemblies)s.fasta
"""
all_targets=[]
for molicute in targets:
	these_targets=[]
	avail_assemblies=sorted([x for x in os.listdir("../"+molicute) if "fa" in os.path.splitext(x)[1]])
	# print molicute
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
			if assemblies!="CLC_MIRA":
				continue
			print template%{'molicute':molicute,'assemblies':assemblies,'assemblyAfname':a1,'assemblyBfname':a2}
			all_targets.append(molicute+"_"+assemblies)
			these_targets.append(molicute+"_"+assemblies)
	print "%s:%s"%(molicute," ".join(these_targets))
print "AllMoli:%s"%(" ".join(all_targets))
