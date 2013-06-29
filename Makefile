####
# Example usage 
# make rhodo_AP_SP_mix.fasta will generate the Mixed assembly from the assemblies of AllPaths and SOAP. 
# make rhodo_quast will generate QUAST reports for a selection of Rhodbacter assemblies
###

pyinterp:= /Library/Frameworks/Python.framework/Versions/Current/bin/python
pyinterp:= /Library/Frameworks/Python.framework/Versions/Current/bin/ipython -i --
pyinterp:=python
	
MUMmer:=$(shell pwd)/bin/MUMmer/
MIXPARAMS:=-A 500 -C 0
MIXPARAMSFOLDER:=A500_C0
QUAST:=bin/quast-2.1/quast.py --min-contig 200 --threads 2
# QUAST:=bin/generate_assembly_command.py

export PATH := $(MUMmer):$(PATH)

.PHONY: clean 

clean: 
	rm -rf temp_assemblies/*
	rm -rf result_statistics/rhodo_quast result_statistics/b_cereus_quast result_statistics/B_fragilis_quast

#IMPORTANT: Indicate here how to match the concatenated Fastas with the assemblers result
# Indicate the dependancies for the concatenated files
# Using these dependancies, the next rule is fired and the files are concatenated 


# Rhodbacter sphaeroides GAGE-B 
RHODO:=datasets/GAGE-B/R_sphaeroides_HiSeq
RHODOGAM:=datasets/GAM-NGS/GAGE-B/R_sphaeroides_HiSeq/
temp_assemblies/rhodo_MS_SOAP.fasta: $(RHODO)/msrca_ctg.fasta $(RHODO)/soap_ctg.fasta
temp_assemblies/rhodo_MI_MS.fasta: $(RHODO)/mira_ctg.fasta $(RHODO)/msrca_ctg.fasta
temp_assemblies/rhodo_MI_SOAP.fasta: $(RHODO)/mira_ctg.fasta $(RHODO)/soap_ctg.fasta
temp_assemblies/rhodo_CBG_SOAP.fasta: $(RHODO)/cabog_ctg.fasta $(RHODO)/soap_ctg.fasta
temp_assemblies/rhodo_MS_MI_SOAP.fasta: $(RHODO)/msrca_ctg.fasta $(RHODO)/mira_ctg.fasta $(RHODO)/soap_ctg.fasta
temp_assemblies/rhodo_MS_CBG_SOAP.fasta: $(RHODO)/msrca_ctg.fasta $(RHODO)/cabog_ctg.fasta $(RHODO)/soap_ctg.fasta
temp_assemblies/rhodo_MS_SPD.fasta: $(RHODO)/msrca_ctg.fasta $(RHODO)/spades_ctg.fasta 

# temp_assemblies/rhodo_AP_SP.fasta: datasets/GAGE/Rhodobacter_sphaeroides/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Rhodobacter_sphaeroides/SOAPdenovo/SOAP.genome.ctg.fasta
# temp_assemblies/rhodo_AP_BB.fasta: datasets/GAGE/Rhodobacter_sphaeroides/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Rhodobacter_sphaeroides/Bambus2/bambus2.ctg.fasta
# temp_assemblies/rhodo_BB_SP.fasta: datasets/GAGE/Rhodobacter_sphaeroides/Bambus2/bambus2.ctg.fasta datasets/GAGE/Rhodobacter_sphaeroides/SOAPdenovo/SOAP.genome.ctg.fasta
# temp_assemblies/rhodo_AP_BB_SP.fasta: datasets/GAGE/Rhodobacter_sphaeroides/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Rhodobacter_sphaeroides/Bambus2/bambus2.ctg.fasta datasets/GAGE/Rhodobacter_sphaeroides/SOAPdenovo/SOAP.genome.ctg.fasta

result_statistics/rhodo_quast: \
	result_assemblies/rhodo_MS_SOAP_mix.fasta \
	result_assemblies/rhodo_MI_MS_mix.fasta \
	result_assemblies/rhodo_MI_SOAP_mix.fasta \
	result_assemblies/rhodo_CBG_SOAP_mix.fasta \
	result_assemblies/rhodo_MS_MI_SOAP_mix.fasta \
	result_assemblies/rhodo_MS_CBG_SOAP_mix.fasta	 \
	result_assemblies/rhodo_MS_SPD_mix.fasta	\
	$(RHODO)/soap_ctg.fasta \
	$(RHODO)/msrca_ctg.fasta \
	$(RHODO)/mira_ctg.fasta \
	$(RHODO)/cabog_ctg.fasta \
	$(RHODO)/spades_ctg.fasta \
	$(RHODOGAM)/GAM_mira-msrca.fasta \
	$(RHODOGAM)/GAM_msrca-soap.fasta \
	$(RHODOGAM)/GAM_soap-msrca.fasta \
	$(RHODOGAM)/GAM_mira-soap.fasta \
	$(RHODOGAM)/GAM_msrca-spades.fasta
	rm -rf result_statistics/rhodo_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Rhodobacter_sphaeroides/Rhodobacter_sphaeroides_ref.fasta -G datasets/reference/Rhodobacter_sphaeroides/Rhodobacter_sphaeroides_ref.gff  $^ 






# GAGE-B assemblies abyss_ctg.fasta cabog_ctg.fasta mira_ctg.fasta msrca_ctg.fasta sga_ctg.fasta soap_ctg.fasta spades_ctg.fasta velvet_ctg.fasta

# Best GAGE-B are 1)Masurca 2)SPADES 3)SOAP-de novo & 4) MIRA 
# temp_assemblies/aureus_AP_SP.fasta: datasets/GAGE/Staphylococcus_aureus/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/SOAPdenovo/SOAP.genome.ctg.fasta
# temp_assemblies/aureus_AP_BB.fasta: datasets/GAGE/Staphylococcus_aureus/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/Bambus2/bambus2.genome.ctg.fasta
# temp_assemblies/aureus_BB_SP.fasta: datasets/GAGE/Staphylococcus_aureus/Bambus2/bambus2.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/SOAPdenovo/SOAP.genome.ctg.fasta
# temp_assemblies/aureus_AP_MS.fasta: datasets/GAGE/Staphylococcus_aureus/Allpaths-LG/ap.genome.ctg.fasta $(AUREUS)/msrca_ctg.fasta
# temp_assemblies/aureus_BB_SP-B.fasta: datasets/GAGE/Staphylococcus_aureus/Bambus2/bambus2.genome.ctg.fasta $(AUREUS)/soap_ctg.fasta
# temp_assemblies/aureus_AP_BB_SP.fasta: datasets/GAGE/Staphylococcus_aureus/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/Bambus2/bambus2.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/SOAPdenovo/SOAP.genome.ctg.fasta
# temp_assemblies/aureus_AP_BB_MS.fasta: datasets/GAGE/Staphylococcus_aureus/Allpaths-LG/ap.genome.ctg.fasta datasets/GAGE/Staphylococcus_aureus/Bambus2/bambus2.genome.ctg.fasta $(AUREUS)/msrca_ctg.fasta
AUREUS:=datasets/GAGE-B/S_aureus_HiSeq
AUREUSGAM:=datasets/GAM-NGS/GAGE-B/S_aureus_HiSeq
temp_assemblies/aureus_MS_SPD.fasta: $(AUREUS)/msrca_ctg.fasta $(AUREUS)/spades_ctg.fasta
temp_assemblies/aureus_MS_SOAP.fasta: $(AUREUS)/msrca_ctg.fasta $(AUREUS)/soap_ctg.fasta
temp_assemblies/aureus_SPD_SOAP.fasta: $(AUREUS)/spades_ctg.fasta $(AUREUS)/soap_ctg.fasta
temp_assemblies/aureus_MS_SPD_SOAP.fasta: $(AUREUS)/msrca_ctg.fasta $(AUREUS)/spades_ctg.fasta $(AUREUS)/soap_ctg.fasta
temp_assemblies/aureus_MI_MS.fasta: $(AUREUS)/mira_ctg.fasta $(AUREUS)/msrca_ctg.fasta
temp_assemblies/aureus_MI_SOAP.fasta: $(AUREUS)/mira_ctg.fasta $(AUREUS)/soap_ctg.fasta
temp_assemblies/aureus_MI_SOAP_SPD_MS.fasta: $(AUREUS)/mira_ctg.fasta $(AUREUS)/soap_ctg.fasta $(AUREUS)/spades_ctg.fasta $(AUREUS)/msrca_ctg.fasta

result_statistics/aureus_quast: \
	result_assemblies/aureus_MS_SPD_mix.fasta \
	result_assemblies/aureus_MS_SOAP_mix.fasta \
	result_assemblies/aureus_SPD_SOAP_mix.fasta \
	result_assemblies/aureus_MS_SPD_SOAP_mix.fasta \
	result_assemblies/aureus_MI_MS_mix.fasta \
	result_assemblies/aureus_MI_SOAP_mix.fasta	 \
	result_assemblies/aureus_MI_SOAP_SPD_MS_mix.fasta\
	$(AUREUS)/soap_ctg.fasta \
	$(AUREUS)/msrca_ctg.fasta \
	$(AUREUS)/mira_ctg.fasta \
	$(AUREUS)/spades_ctg.fasta \
	$(AUREUSGAM)/GAM_mira-msrca.fasta \
	$(AUREUSGAM)/GAM_msrca-soap.fasta \
	$(AUREUSGAM)/GAM_msrca-spades.fasta \
	$(AUREUSGAM)/GAM_spades-soap.fasta \
	$(AUREUSGAM)/GAM_soap-mira.fasta
	rm -rf result_statistics/aureus_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Staphylococcus_aureus/Staphylococcus_aureus_ref.fasta  -G datasets/reference/Staphylococcus_aureus/Staphylococcus_aureus_ref.gff $^ 




# Best GAGE for B. cereus are MaSuRCA CABOG, SOAP  and MIRA 
CEREUS:=datasets/GAGE-B/B_cereus_MiSeq
CEREUSGAM:=datasets/GAM-NGS/GAGE-B/B_cereus_MiSeq/
temp_assemblies/b_cereus_AB_SP.fasta: $(CEREUS)/abyss_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_AB_MS.fasta: $(CEREUS)/abyss_ctg.fasta $(CEREUS)/msrca_ctg.fasta
temp_assemblies/b_cereus_MS_SP.fasta: $(CEREUS)/msrca_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_AB_MS_SP.fasta: $(CEREUS)/abyss_ctg.fasta $(CEREUS)/msrca_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_MS_SOAP.fasta: $(CEREUS)/msrca_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_CBG_SOAP.fasta: $(CEREUS)/cabog_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_MS_CBG_SOAP.fasta: $(CEREUS)/msrca_ctg.fasta $(CEREUS)/cabog_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_MI_MS.fasta: $(CEREUS)/mira_ctg.fasta $(CEREUS)/msrca_ctg.fasta
temp_assemblies/b_cereus_MI_SOAP.fasta: $(CEREUS)/mira_ctg.fasta $(CEREUS)/soap_ctg.fasta
temp_assemblies/b_cereus_MI_SOAP_CBG_MS.fasta: $(CEREUS)/mira_ctg.fasta $(CEREUS)/soap_ctg.fasta $(CEREUS)/cabog_ctg.fasta $(CEREUS)/msrca_ctg.fasta

result_statistics/b_cereus_quast: \
	result_assemblies/b_cereus_AB_SP_mix.fasta \
	result_assemblies/b_cereus_AB_MS_mix.fasta \
	result_assemblies/b_cereus_MS_SP_mix.fasta \
	result_assemblies/b_cereus_AB_MS_SP_mix.fasta \
	result_assemblies/b_cereus_MS_SOAP_mix.fasta \
	result_assemblies/b_cereus_CBG_SOAP_mix.fasta \
	result_assemblies/b_cereus_MS_CBG_SOAP_mix.fasta \
	result_assemblies/b_cereus_MI_MS_mix.fasta \
	result_assemblies/b_cereus_MI_SOAP_mix.fasta \
	result_assemblies/b_cereus_MI_SOAP_CBG_MS_mix.fasta \
	$(CEREUS)/abyss_ctg.fasta \
	$(CEREUS)/soap_ctg.fasta \
	$(CEREUS)/msrca_ctg.fasta \
	$(CEREUS)/spades_ctg.fasta \
	$(CEREUS)/cabog_ctg.fasta \
	$(CEREUS)/mira_ctg.fasta \
	$(CEREUSGAM)/GAM_mira-cabog.fasta\
	$(CEREUSGAM)/GAM_mira-msrca.fasta\
	$(CEREUSGAM)/GAM_msrca-cabog.fasta\
	$(CEREUSGAM)/GAM_soap-msrca.fasta
	rm -rf result_statistics/b_cereus_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Bacillus_cereus/Bacillus_cereus_ref.fasta -G datasets/reference/Bacillus_cereus/Bacillus_cereus_ref.gff  $^ 



HYDROPHILA:=datasets/GAGE-B/A_hydrophila_HiSeq
HYDROPHILAGAM:=datasets/GAM-NGS/GAGE-B/A_hydrophila_HiSeq/
temp_assemblies/a_hydrophila_AB_SP.fasta: $(HYDROPHILA)/abyss_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_AB_MS.fasta: $(HYDROPHILA)/abyss_ctg.fasta $(HYDROPHILA)/msrca_ctg.fasta
temp_assemblies/a_hydrophila_MS_SP.fasta: $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_AB_MS_SP.fasta: $(HYDROPHILA)/abyss_ctg.fasta $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_MS_SPD.fasta: $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/spades_ctg.fasta
temp_assemblies/a_hydrophila_MS_SOAP.fasta: $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_SPD_SOAP.fasta: $(HYDROPHILA)/spades_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_MS_SPD_SOAP.fasta: $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/spades_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_MI_MS.fasta: $(HYDROPHILA)/mira_ctg.fasta $(HYDROPHILA)/msrca_ctg.fasta
temp_assemblies/a_hydrophila_MI_SOAP.fasta: $(HYDROPHILA)/mira_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta
temp_assemblies/a_hydrophila_MI_SOAP_SPD_MS.fasta: $(HYDROPHILA)/mira_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta $(HYDROPHILA)/spades_ctg.fasta $(HYDROPHILA)/msrca_ctg.fasta
temp_assemblies/a_hydrophila_MS_CBG.fasta: $(HYDROPHILA)/msrca_ctg.fasta $(HYDROPHILA)/cabog_ctg.fasta
temp_assemblies/a_hydrophila_MI_CBG.fasta: $(HYDROPHILA)/mira_ctg.fasta $(HYDROPHILA)/cabog_ctg.fasta
temp_assemblies/a_hydrophila_CBG_SOAP.fasta: $(HYDROPHILA)/cabog_ctg.fasta $(HYDROPHILA)/soap_ctg.fasta

result_statistics/a_hydrophila_quast: \
	result_assemblies/a_hydrophila_AB_SP_mix.fasta \
	result_assemblies/a_hydrophila_AB_MS_mix.fasta \
	result_assemblies/a_hydrophila_MS_SP_mix.fasta \
	result_assemblies/a_hydrophila_AB_MS_SP_mix.fasta \
	result_assemblies/a_hydrophila_MS_SPD_mix.fasta \
	result_assemblies/a_hydrophila_MS_SOAP_mix.fasta \
	result_assemblies/a_hydrophila_SPD_SOAP_mix.fasta \
	result_assemblies/a_hydrophila_MS_SPD_SOAP_mix.fasta \
	result_assemblies/a_hydrophila_MI_MS_mix.fasta \
	result_assemblies/a_hydrophila_MI_SOAP_mix.fasta \
	result_assemblies/a_hydrophila_MI_SOAP_SPD_MS_mix.fasta \
	result_assemblies/a_hydrophila_MS_CBG_mix.fasta \
	result_assemblies/a_hydrophila_MI_CBG_mix.fasta \
	result_assemblies/a_hydrophila_CBG_SOAP_mix.fasta \
	$(HYDROPHILA)/abyss_ctg.fasta \
	$(HYDROPHILA)/soap_ctg.fasta \
	$(HYDROPHILA)/msrca_ctg.fasta \
	$(HYDROPHILA)/mira_ctg.fasta \
	$(HYDROPHILA)/spades_ctg.fasta \
	$(HYDROPHILAGAM)/GAM_mira-cabog.fasta \
	$(HYDROPHILAGAM)/GAM_msrca-cabog.fasta \
	$(HYDROPHILAGAM)/GAM_msrca-spades.fasta \
	$(HYDROPHILAGAM)/GAM_soap-cabog.fasta
	rm -rf result_statistics/a_hydrophila_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Aeromonas_hydrophila/Aeromonas_hydrophila_ref.fasta  -G datasets/reference/Aeromonas_hydrophila/Aeromonas_hydrophila_ref.gff $^ 


# Rules for GAGE-B
# VCHOLERAE Hiseq best are MIRA, MaSuRCA, ABySS, SOAP 
VCHOLERAE_HS:=datasets/GAGE-B/V_cholerae_HiSeq/
VCHOLERAEGAM:=datasets/GAM-NGS/GAGE-B/V_cholerae_HiSeq
temp_assemblies/V_cholerae_AB_SP.fasta: $(VCHOLERAE_HS)/abyss_ctg.fasta $(VCHOLERAE_HS)/soap_ctg.fasta
temp_assemblies/V_cholerae_AB_MS.fasta: $(VCHOLERAE_HS)/abyss_ctg.fasta $(VCHOLERAE_HS)/msrca_ctg.fasta
temp_assemblies/V_cholerae_MS_SP.fasta: $(VCHOLERAE_HS)/msrca_ctg.fasta $(VCHOLERAE_HS)/soap_ctg.fasta
temp_assemblies/V_cholerae_MS_MI.fasta: $(VCHOLERAE_HS)/msrca_ctg.fasta $(VCHOLERAE_HS)/mira_ctg.fasta
temp_assemblies/V_cholerae_AB_MS_SP.fasta: $(VCHOLERAE_HS)/abyss_ctg.fasta $(VCHOLERAE_HS)/msrca_ctg.fasta $(VCHOLERAE_HS)/soap_ctg.fasta
temp_assemblies/V_cholerae_MI_MS_SP.fasta: $(VCHOLERAE_HS)/mira_ctg.fasta $(VCHOLERAE_HS)/msrca_ctg.fasta $(VCHOLERAE_HS)/soap_ctg.fasta


result_statistics/V_cholerae_quast: \
	result_assemblies/V_cholerae_AB_SP_mix.fasta \
	result_assemblies/V_cholerae_AB_MS_mix.fasta \
	result_assemblies/V_cholerae_MS_SP_mix.fasta \
	result_assemblies/V_cholerae_MS_MI_mix.fasta \
	result_assemblies/V_cholerae_AB_MS_SP_mix.fasta \
	result_assemblies/V_cholerae_MI_MS_SP_mix.fasta \
	$(VCHOLERAE_HS)/abyss_ctg.fasta \
	$(VCHOLERAE_HS)/msrca_ctg.fasta \
	$(VCHOLERAE_HS)/mira_ctg.fasta \
	$(VCHOLERAE_HS)/soap_ctg.fasta \
	$(VCHOLERAEGAM)/GAM_abyss-msrca.fasta \
	$(VCHOLERAEGAM)/GAM_mira-msrca.fasta \
	$(VCHOLERAEGAM)/GAM_msrca-abyss.fasta
	rm -rf result_statistics/V_cholerae_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Vibrio_cholerae/Vibrio_cholerae_ref.fasta  -G datasets/reference/Vibrio_cholerae/Vibrio_cholerae_ref.gff $^ 


# M.Abscessus
MABSCESSUS:=datasets/GAGE-B/M_abscessus_HiSeq/
MABSCESSUSGAM:=datasets/GAM-NGS/GAGE-B/M_abscessus_HiSeq/
temp_assemblies/M_abscessus_AB_SP.fasta: $(MABSCESSUS)/abyss_ctg.fasta $(MABSCESSUS)/soap_ctg.fasta
temp_assemblies/M_abscessus_AB_MS.fasta: $(MABSCESSUS)/abyss_ctg.fasta $(MABSCESSUS)/msrca_ctg.fasta
temp_assemblies/M_abscessus_MS_SP.fasta: $(MABSCESSUS)/msrca_ctg.fasta $(MABSCESSUS)/soap_ctg.fasta
temp_assemblies/M_abscessus_SP_SPD.fasta: $(MABSCESSUS)/soap_ctg.fasta $(MABSCESSUS)/spades_ctg.fasta
temp_assemblies/M_abscessus_AB_MS_SP.fasta: $(MABSCESSUS)/abyss_ctg.fasta $(MABSCESSUS)/msrca_ctg.fasta $(MABSCESSUS)/soap_ctg.fasta

result_statistics/M_abscessus_quast: \
	result_assemblies/M_abscessus_AB_SP_mix.fasta \
	result_assemblies/M_abscessus_AB_MS_mix.fasta \
	result_assemblies/M_abscessus_MS_SP_mix.fasta \
	result_assemblies/M_abscessus_SP_SPD_mix.fasta \
	result_assemblies/M_abscessus_AB_MS_SP_mix.fasta \
	datasets/GAGE-B/M_abscessus_HiSeq/abyss_ctg.fasta \
	datasets/GAGE-B/M_abscessus_HiSeq/soap_ctg.fasta \
	datasets/GAGE-B/M_abscessus_HiSeq/msrca_ctg.fasta \
	$(MABSCESSUSGAM)/GAM_abyss-msrca.fasta \
	$(MABSCESSUSGAM)/GAM_msrca-abyss.fasta \
	$(MABSCESSUSGAM)/GAM_msrca-spades.fasta 
	rm -rf result_statistics/M_abscessus_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Mycobacterium_abscessus/Mycobacterium_abscessus_ref.fasta  -G datasets/reference/Mycobacterium_abscessus/Mycobacterium_abscessus_ref.gff $^ 



# X.axonopodis
XAXONOPODIS:=datasets/GAGE-B/X_axonopodis_HiSeq/
XAXONOPODISGAM:=datasets/GAM-NGS/GAGE-B/X_axonopodis_HiSeq
temp_assemblies/X_axonopodis_AB_SP.fasta: $(XAXONOPODIS)/abyss_ctg.fasta $(XAXONOPODIS)/soap_ctg.fasta
temp_assemblies/X_axonopodis_AB_MS.fasta: $(XAXONOPODIS)/abyss_ctg.fasta $(XAXONOPODIS)/msrca_ctg.fasta
temp_assemblies/X_axonopodis_MS_SP.fasta: $(XAXONOPODIS)/msrca_ctg.fasta $(XAXONOPODIS)/soap_ctg.fasta
temp_assemblies/X_axonopodis_MS_SPD.fasta: $(XAXONOPODIS)/msrca_ctg.fasta $(XAXONOPODIS)/spades_ctg.fasta
temp_assemblies/X_axonopodis_AB_MS_SP.fasta: $(XAXONOPODIS)/abyss_ctg.fasta $(XAXONOPODIS)/msrca_ctg.fasta $(XAXONOPODIS)/soap_ctg.fasta

result_statistics/X_axonopodis_quast: \
	result_assemblies/X_axonopodis_AB_SP_mix.fasta \
	result_assemblies/X_axonopodis_AB_MS_mix.fasta \
	result_assemblies/X_axonopodis_MS_SP_mix.fasta \
	result_assemblies/X_axonopodis_MS_SPD_mix.fasta \
	result_assemblies/X_axonopodis_AB_MS_SP_mix.fasta \
	$(XAXONOPODIS)/abyss_ctg.fasta \
	$(XAXONOPODIS)/soap_ctg.fasta \
	$(XAXONOPODIS)/msrca_ctg.fasta \
	$(XAXONOPODIS)/spades_ctg.fasta \
	$(XAXONOPODISGAM)/GAM_abyss-soap.fasta \
	$(XAXONOPODISGAM)/GAM_msrca-abyss.fasta \
	$(XAXONOPODISGAM)/GAM_msrca-spades.fasta 
	rm -rf result_statistics/X_axonopodis_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Xanthomonas_axonopodis/Xanthomonas_axonopodis_ref.fasta  -G datasets/reference/Xanthomonas_axonopodis/Xanthomonas_axonopodis_ref.gff $^ 




# B. Fragilis
BFRAGILIS:=datasets/GAGE-B/B_fragilis_HiSeq/
BFRAGILISGAM:=datasets/GAM-NGS/GAGE-B/B_fragilis_HiSeq
temp_assemblies/B_fragilis_AB_SP.fasta: $(BFRAGILIS)/abyss_ctg.fasta $(BFRAGILIS)/soap_ctg.fasta
temp_assemblies/B_fragilis_AB_MS.fasta: $(BFRAGILIS)/abyss_ctg.fasta $(BFRAGILIS)/msrca_ctg.fasta
temp_assemblies/B_fragilis_MS_SP.fasta: $(BFRAGILIS)/msrca_ctg.fasta $(BFRAGILIS)/soap_ctg.fasta
temp_assemblies/B_fragilis_MS_SPD.fasta: $(BFRAGILIS)/msrca_ctg.fasta $(BFRAGILIS)/spades_ctg.fasta
temp_assemblies/B_fragilis_AB_MS_SP.fasta: $(BFRAGILIS)/abyss_ctg.fasta $(BFRAGILIS)/msrca_ctg.fasta $(BFRAGILIS)/soap_ctg.fasta
temp_assemblies/B_fragilis_AB_MS_SPD.fasta: $(BFRAGILIS)/abyss_ctg.fasta $(BFRAGILIS)/msrca_ctg.fasta $(BFRAGILIS)/spades_ctg.fasta

result_statistics/B_fragilis_quast: \
	result_assemblies/B_fragilis_AB_SP_mix.fasta \
	result_assemblies/B_fragilis_AB_MS_mix.fasta \
	result_assemblies/B_fragilis_MS_SP_mix.fasta \
	result_assemblies/B_fragilis_MS_SPD_mix.fasta \
	result_assemblies/B_fragilis_AB_MS_SP_mix.fasta \
	result_assemblies/B_fragilis_AB_MS_SPD_mix.fasta \
	$(BFRAGILIS)/abyss_ctg.fasta \
	$(BFRAGILIS)/soap_ctg.fasta \
	$(BFRAGILIS)/msrca_ctg.fasta \
	$(BFRAGILIS)/spades_ctg.fasta \
	$(BFRAGILISGAM)/GAM_abyss-soap.fasta \
	$(BFRAGILISGAM)/GAM_msrca-spades.fasta \
	$(BFRAGILISGAM)/GAM_soap-msrca.fasta \
	$(BFRAGILISGAM)/GAM_spades-msrca.fasta
	rm -rf result_statistics/B_fragilis_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.fasta  -G datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.gff $^ 




include Makefiles/Makefile_mollicutes.GNUmakefile
include Makefiles/Makefile_MixGAGEB.GNUmakefile
include Makefiles/Makefile_AllGAGEB_quasts.GNUmakefile
include Makefiles/Makefile_AllGAM_over_GAGEB_quasts.GNUmakefile
include Makefiles/Makefile_mollicutesQuasts.GNUmakefile
include Makefiles/Makefile_mycoplasmas_triple_mixes.GNUmakefile




result_statistics/%_quast: result_assemblies/%_CLC_MIRA_mix.fasta result_assemblies/%_AB_CLC_mix.fasta result_assemblies/%_AB_MIRA_mix.fasta 
	$(pyinterp) $(QUAST) -o $@ $^ 

%_assembly: result_assemblies/%_AB_CLC_mix.fasta result_assemblies/%_CLC_MIRA_mix.fasta result_assemblies/%_AB_MIRA_mix.fasta
MolliStats: result_statistics/MOVI_quast result_statistics/MMC_quast result_statistics/MSCe_quast result_statistics/MSCd_quast result_statistics/MSCc_quast result_statistics/MSCb_quast result_statistics/MBVG_quast






clean_for_git:
	rm -rf result_statistics/*/contigs_reports



rhodoExpansion: temp_assemblies/rhodo_AP_BB_SP.coords
	bin/nucmer_coords_filter.py temp_assemblies/rhodo_AP_BB_SP.coords A1_119 A3_scaffold13.3 A3_scaffold13.5 > temp_assemblies/TEST_rhodoExpansion.coords
	# Display original input
	bin/display_coords.py temp_assemblies/TEST_rhodoExpansion.coords
	# Execute mix
	$(pyinterp) bin/Mix.py --graph -r -a temp_assemblies/TEST_rhodoExpansion.coords -o TEST_rhodoExpansion -c temp_assemblies/rhodo_AP_BB_SP.fasta -A 125 -C 0
	# Check the resulting assembly with nucmer 
	bin/MUMmer3.23/nucmer -p "TEST_rhodoExpansion_assembly" --maxmatch -l 30 -banded TEST_rhodoExpansion/Mix_results_A125_C0/Mix_assembly.fasta TEST_rhodoExpansion/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
	../MUMmer3.23/show-coords -l -c TEST_rhodoExpansion_assembly.delta > TEST_rhodoExpansion_assembly.coords
	cat TEST_rhodoExpansion_assembly.coords
	./display_coords.py TEST_rhodoExpansion_assembly.coords
	./compute_n50.pl TEST_rhodoExpansion/Mix_results_A125_C0/Mix_assembly.fasta




# End of configurable part

# Variables reminder: 
# $^ all pre-requisites
# $@ name of the target 
# $(@F) file part of the target 
# $(basename $(@F)) file part of the target with extension removed 

# Given k fasta file, generate a single file with concatenated and unify IDs
temp_assemblies/%.fasta: bin/preprocessing.py
	$(eval TAG:= $(basename $(@F))) 
	$(pyinterp) bin/preprocessing.py $(filter %.fasta,$^) -o temp_assemblies/$(TAG).fasta


# Indicate that coords files should not be removed upon completion
.SECONDARY:

# Given a fasta file, generate a coords file with the coordinates of the alignments
temp_assemblies/%.coords: temp_assemblies/%.fasta
	$(eval TAG:= $(basename $(@F))) 
#	cd temp_assemblies; nucmer -p $(TAG) --maxmatch -l 30 -banded $(TAG).fasta $(TAG).fasta 2>/dev/null; show-coords -l -c $(TAG).delta > $(TAG).coords
	cd temp_assemblies; nucmer -p $(TAG) --maxmatch -c 30 -l 30 -banded $(TAG).fasta $(TAG).fasta 2>$(TAG).nucmer.log.txt; show-coords -l -c $(TAG).delta > $(TAG).coords

#	rm temp_assemblies/$(TAG).delta

# Given a fasta file and a coord file of alignments, run Mix and move the resulting assembly to the result folder
result_assemblies/%_mix.fasta: temp_assemblies/%.coords temp_assemblies/%.fasta bin/Mix.py bin/graph.py bin/integer_set.py
	$(eval TAG:= $(basename $(@F))) 
	$(pyinterp) bin/Mix.py $(MIXPARAMS) -o result_assemblies/$(TAG)  -a $(filter %.coords,$^) -c $(filter %.fasta,$^) 2> temp_assemblies/$(TAG).mix.log.txt
	mv result_assemblies/$(TAG)/Mix_results_$(MIXPARAMSFOLDER)/Mix_assembly.fasta result_assemblies/$(TAG).fasta