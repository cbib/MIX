pyinterp:=python
	
MUMmer:=$(shell pwd)/bin/MUMmer3.23/
MIXPARAMS:=-A 500 -C 0
MIXPARAMSFOLDER:=A500_C0
QUAST:=bin/quast-2.1/quast.py --min-contig 200 --threads 8
GAGEBREF:=datasets/GAGE-B

result_statistics/gage-b-mono/R_sphaeroides_HiSeq: $(GAGEBREF)/R_sphaeroides_HiSeq/abyss_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/cabog_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/mira_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/msrca_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/sga_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/soap_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/spades_ctg.fasta $(GAGEBREF)/R_sphaeroides_HiSeq/velvet_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Rhodobacter_sphaeroides/Rhodobacter_sphaeroides_ref.fasta -G datasets/reference/Rhodobacter_sphaeroides/Rhodobacter_sphaeroides_ref.gff  $^ 

result_statistics/gage-b-mono/A_hydrophila_HiSeq: $(GAGEBREF)/A_hydrophila_HiSeq/abyss_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/cabog_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/mira_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/msrca_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/sga_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/soap_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/spades_ctg.fasta $(GAGEBREF)/A_hydrophila_HiSeq/velvet_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Aeromonas_hydrophila/Aeromonas_hydrophila_ref.fasta  -G datasets/reference/Aeromonas_hydrophila/Aeromonas_hydrophila_ref.gff $^ 

result_statistics/gage-b-mono/B_cereus_MiSeq: $(GAGEBREF)/B_cereus_MiSeq/abyss_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/cabog_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/mira_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/msrca_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/sga_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/soap_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/spades_ctg.fasta $(GAGEBREF)/B_cereus_MiSeq/velvet_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/abyss_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/cabog_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Bacillus_cereus/Bacillus_cereus_ref.fasta -G datasets/reference/Bacillus_cereus/Bacillus_cereus_ref.gff  $^ 


result_statistics/gage-b-mono/B_fragilis_HiSeq: $(GAGEBREF)/B_fragilis_HiSeq/mira_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/msrca_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/sga_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/soap_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/spades_ctg.fasta $(GAGEBREF)/B_fragilis_HiSeq/velvet_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/abyss_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/cabog_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.fasta  -G datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.gff $^ 

result_statistics/gage-b-mono/M_abscessus_HiSeq: $(GAGEBREF)/M_abscessus_HiSeq/mira_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/msrca_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/sga_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/soap_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/spades_ctg.fasta $(GAGEBREF)/M_abscessus_HiSeq/velvet_ctg.fasta $(GAGEBREF)/M_abscessus_MiSeq/abyss_ctg.fasta $(GAGEBREF)/M_abscessus_MiSeq/cabog_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Mycobacterium_abscessus/Mycobacterium_abscessus_ref.fasta  -G datasets/reference/Mycobacterium_abscessus/Mycobacterium_abscessus_ref.gff $^ 

result_statistics/gage-b-mono/S_aureus_HiSeq: $(GAGEBREF)/S_aureus_HiSeq/abyss_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/cabog_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/mira_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/msrca_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/sga_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/soap_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/spades_ctg.fasta $(GAGEBREF)/S_aureus_HiSeq/velvet_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/abyss_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/cabog_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Staphylococcus_aureus/Staphylococcus_aureus_ref.fasta  -G datasets/reference/Staphylococcus_aureus/Staphylococcus_aureus_ref.gff $^ 

result_statistics/gage-b-mono/V_cholerae_HiSeq: $(GAGEBREF)/V_cholerae_HiSeq/mira_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/msrca_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/sga_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/soap_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/spades_ctg.fasta $(GAGEBREF)/V_cholerae_HiSeq/velvet_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/abyss_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/cabog_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Vibrio_cholerae/Vibrio_cholerae_ref.fasta  -G datasets/reference/Vibrio_cholerae/Vibrio_cholerae_ref.gff $^ 

result_statistics/gage-b-mono/X_axonopodis_HiSeq: $(GAGEBREF)/X_axonopodis_HiSeq/mira_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/msrca_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/sga_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/soap_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/spades_ctg.fasta $(GAGEBREF)/X_axonopodis_HiSeq/velvet_ctg.fasta
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Xanthomonas_axonopodis/Xanthomonas_axonopodis_ref.fasta  -G datasets/reference/Xanthomonas_axonopodis/Xanthomonas_axonopodis_ref.gff $^ 


ALLGAGEBQUASTS:result_statistics/gage-b-mono/R_sphaeroides_HiSeq \
result_statistics/gage-b-mono/A_hydrophila_HiSeq \
result_statistics/gage-b-mono/B_cereus_MiSeq \
result_statistics/gage-b-mono/B_fragilis_HiSeq \
result_statistics/gage-b-mono/M_abscessus_HiSeq \
result_statistics/gage-b-mono/S_aureus_HiSeq \
result_statistics/gage-b-mono/V_cholerae_HiSeq \
result_statistics/gage-b-mono/X_axonopodis_HiSeq 