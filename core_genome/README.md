Computation of the core genome conservation for MMC
======================================================
* Based on orthology cluster defining the core Mollicute genome
The core mollicute genome consist of the proteins found in the clade. These proteins are clustered by similarity.
 
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
    /usr/share/Modules/apps/blast/2.2.25/bin/blastall -p tblastn -m 8 -d ./MOVI  -i ../CorePep.fasta -o MOVI.hits -e 0.1

# To generate the *.hits.core files 
    for f in *hits; do perl ./clusters.perl --blast   > .hits.core; done


# *.hits.core 
On prend le meilleur alignement parmi chaque cluster (Par E-value)
si ce meilleur alignement couvre plus de 50% de la prot correspondante, alors 1 dans la col du cluster 
170 clust de ~ 30 prot, les 30 sont assez similaires 
par ligne pourcentage du cluster couvert: , 50%, 85%, 99% couverture (trois colonnes)


# Calcul 
Sur le serveur de calcul rainman, 
Se placer dans le rep coreGenome 
    /usr/share/Modules/apps/blast/2.2.25plus/bin/makeblastdb -dbtype 'nucl' -in MMC_CLC-MIRA_gam-merge.gam.fasta -out MMC_CLC-MIRA_gam -parse_seqids
    /usr/share/Modules/apps/blast/2.2.25/bin/blastall -p tblastn -m 8 -d ./MMC_CLC-MIRA_gam  -i ~/Mycoplasmes/CorePep.fasta -o MMC_CLC-MIRA_gam.hits.tsv -e 0.1
    perl clusters.perl --blast MMC_CLC-MIRA_gam.hits.tsv > MMC_CLC-MIRA_gam.hits.core.tsv 2> foo
