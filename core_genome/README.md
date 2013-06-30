Computation of the core genome conservation for Mycoplasmas genomes
======================================================

* Based on orthology clusters defining the core Mycoplasma genome
The core mollicute genome consists of ortologous conserved proteins found in the clade. These proteins are clustered by similarity.

* CoreGenome 
  * CoreGenome_Mycoplasmas.fasta, file with protein sequences for the core genomes used for the assembly evaluation of 10 Mycoplasma genomes 
    It contains 170 orthologous clusters, 31 proteins per cluster. "1169_Mcapa_s03_5480" indicates that the MCAPA protein is in the cluster 1169.

  * mycoplasmas_all_core_conservation.tsv, file that conatins the core genome conservation results 
    For each orthologous cluster we select the best alignment (in terms of its e-value, over all members of the cluster) and check if this alignment covers over 50%, 85% or 99% of the protein. For each of these, if the coverage is present we put 1 in the corresponding column.


# To generate the *.hits files 
	
    /usr/share/Modules/apps/blast/2.2.25/bin/blastall -p tblastn -m 8 -dbgcode4 -d ./MOVI  -i CoreGenome_Mycoplasmas.fasta -o MOVI.hits -e 0.1

# To generate the *.hits.core files 

    for f in *hits; do perl ./clusters.perl --blast MOVI.hits --prots CoreGenome_170.fasta  > .hits.core; done


# *.hits.core 

For each orthologous cluster we select the best alignment (in terms of its e-value)
We check if this alignment covers over 50%, 85% or 99% of the protein
For each of these, if the coverage is present we put 1 in the corresponding column in the output .core file

# Example session 

    makeblastdb -dbtype 'nucl' -in ../result_assemblies/MMC_CLC_MIRA_mix.fasta -out MMC_CLC-MIRA -parse_seqids
    blastall -p tblastn -m 8 -d ./MMC_CLC-MIRA -i CoreGenome_Mycoplasmas.fasta -o MMC_CLC-MIRA.hits -e 0.1
    perl clusters.perl --blast MMC_CLC-MIRA.hits -prots CoreGenome_Mycoplasmas.fasta > MMC_CLC-MIRA.core
    # [Generate results for all genomes....]
	for f in *core; do n5=`cut -f 2 $f | grep 1 |wc -l`; n8=`cut -f 3 $f | grep 1 |wc -l`; n9=`cut -f 4 $f | grep 1 | wc -l`; echo "$f $n5 $n8 $n9"  >> Core_all; done;

	rm Core_all ; here=`pwd`; for dir in M*; do echo $dir >> Core_all; cat $dir/Core_all >> Core_all; echo >> Core_all; done

