Computation of the core genome conservation for MMC
======================================================
* Based on orthology clusters defining the core Mollicute genome
The core mollicute genome consists of ortologous conserved proteins found in the clade. These proteins are clustered by similarity.
 
* 1169_Mcapa_s03_5480 indicate that the MCAPA protein is in cluster 1169

* Colnames of *.hits files 

1. Query id
1. Subject id
1. %identity
1. alignment length
1. mismatches
1. gap openings
1. q. start
1. q. end
1. s. start
1. s. end
1. e-value
1. bit score


# To generate the *.hits files 

    /usr/share/Modules/apps/blast/2.2.25/bin/blastall -p tblastn -m 8 -dbgcode4 -d ./MOVI  -i ../CorePep.fasta -o MOVI.hits -e 0.1

# To generate the *.hits.core files 

    for f in *hits; do perl ./clusters.perl --blast MOVI.hits --prots CoreGenome_170.fasta  > .hits.core; done


# *.hits.core 
For each orthologous cluster we select the best alignment (in terms of its e-value)
We check if this alignment covers over 50%, 85% or 99% of the protein
For each of these, if the coverage is present we put 1 in the corresponding column in the output .core file

# Example

    makeblastdb -dbtype 'nucl' -in MMC_CLC-MIRA.fasta -out MMC_CLC-MIRA -parse_seqids
    blastall -p tblastn -m 8 -d ./MMC_CLC-MIRA -i CoreGenome_170.fasta -o MMC_CLC-MIRA.hits -e 0.1
    perl clusters.perl --blast MMC_CLC-MIRA.hits -prots CoreGenome_170.fasta > MMC_CLC-MIRA.core
