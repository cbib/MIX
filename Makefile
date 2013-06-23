####
# Example usage 
# make rhodo_AP_SP_mix.fasta will generate the Mixed assembly from the assemblies of AllPaths and SOAP. 
# make rhodo_quast will generate QUAST reports for a selection of Rhodbacter assemblies
###

pyinterp:= /Library/Frameworks/Python.framework/Versions/Current/bin/python
pyinterp:= /Library/Frameworks/Python.framework/Versions/Current/bin/ipython -i --
pyinterp:=python
MUMmer:=$(shell pwd)/bin/MUMmer3.23/
MIXPARAMS:=-A 500 -C 0
MIXPARAMSFOLDER:=A500_C0
QUAST:=bin/quast-2.1/quast.py --min-contig 200

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
	# rm -rf result_statistics/rhodo_quast
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
	# rm -rf result_statistics/aureus_quast
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
	# rm -rf result_statistics/b_cereus_quast
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
	# rm -rf result_statistics/a_hydrophila_quast
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
	# rm -rf result_statistics/V_cholerae_quast
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
	# rm -rf result_statistics/M_abscessus_quast
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
	# rm -rf result_statistics/X_axonopodis_quast
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
	# rm -rf result_statistics/B_fragilis_quast
	$(pyinterp) $(QUAST) -o $@ -R datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.fasta  -G datasets/reference/Bacteroides_fragilis/Bacteroides_fragilis_ref.gff $^ 



# Rules for Mollicutes assemblies and statistics 
# Cf bin/generate_assembly_command.py to generate the rules 
include Makefile_mollicutes.mk
# MOLLI:=datasets/Mollicutes/
# MOLLIGAM:=datasets/GAM-NGS/Mollicutes


# temp_assemblies/MOVI_AB_CLC.fasta: $(MOLLI)/MOVI/AIP_ABySS_32-scaffolds.fasta $(MOLLI)/MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MOVI_AB_MIRA.fasta: $(MOLLI)/MOVI/AIP_ABySS_32-scaffolds.fasta $(MOLLI)/MOVI/AIP_step2_out.unpadded.fasta
# temp_assemblies/MOVI_CLC_MIRA.fasta: $(MOLLI)/MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MOVI/AIP_step2_out.unpadded.fasta

# temp_assemblies/MMC_AB_CLC.fasta: $(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MMC_AB_MIRA.fasta: $(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta
# temp_assemblies/MMC_CLC_MIRA.fasta: $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta

# result_statistics/MMC_quast: \
# 	$(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta\
# 	$(MOLLIGAM)/MMC/GAM_abyss-CLC.fasta $(MOLLIGAM)/MMC/GAM_CLC-mira.fasta $(MOLLIGAM)/MMC/GAM_mira-abyss.fasta \
# 	result_assemblies/MMC_AB_CLC_mix.fasta result_assemblies/MMC_AB_MIRA_mix.fasta result_assemblies/MMC_CLC_MIRA_mix.fasta 
# 	rm -rf result_statistics/MMC_quast
# 	$(pyinterp) $(QUAST) -o $@ $^ 

# temp_assemblies/MSCe_AB_CLC.fasta: $(MOLLI)/MSCe/AKC_ABySS_35-scaffolds.fasta $(MOLLI)/MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MSCe_AB_MIRA.fasta: $(MOLLI)/MSCe/AKC_ABySS_35-scaffolds.fasta $(MOLLI)/MSCe/AKC_step2_out.unpadded.fasta
# temp_assemblies/MSCe_CLC_MIRA.fasta: $(MOLLI)/MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCe/AKC_step2_out.unpadded.fasta

# temp_assemblies/MSCd_AB_CLC.fasta: $(MOLLI)/MSCd/AKE_ABySS_34-scaffolds.fasta $(MOLLI)/MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MSCd_AB_MIRA.fasta: $(MOLLI)/MSCd/AKE_ABySS_34-scaffolds.fasta $(MOLLI)/MSCd/AKE_step2_out.unpadded.fasta
# temp_assemblies/MSCd_CLC_MIRA.fasta: $(MOLLI)/MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCd/AKE_step2_out.unpadded.fasta

# temp_assemblies/MSCc_AB_CLC.fasta: $(MOLLI)/MSCc/AIZ_ABySS_27-scaffolds.fasta $(MOLLI)/MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MSCc_AB_MIRA.fasta: $(MOLLI)/MSCc/AIZ_ABySS_27-scaffolds.fasta $(MOLLI)/MSCc/AIZ_step2_out.unpadded.fasta
# temp_assemblies/MSCc_CLC_MIRA.fasta: $(MOLLI)/MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCc/AIZ_step2_out.unpadded.fasta

# temp_assemblies/MSCb_AB_CLC.fasta: $(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MSCb_AB_MIRA.fasta: $(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta
# temp_assemblies/MSCb_CLC_MIRA.fasta: $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta

# result_statistics/MSCb_quast: \
# 	$(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta \
# 	$(MOLLIGAM)/MSCb/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCb/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCb/GAM_mira-abyss.fasta \
# 	result_assemblies/MSCb_AB_CLC_mix.fasta result_assemblies/MSCb_AB_MIRA_mix.fasta result_assemblies/MSCb_CLC_MIRA_mix.fasta
# 	rm -rf result_statistics/MSCb_quast
# 	$(pyinterp) $(QUAST) -o $@ $^


# temp_assemblies/MBVG_AB_CLC.fasta: $(MOLLI)/MBVG/AIN_ABySS_35-scaffolds.fasta $(MOLLI)/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta
# temp_assemblies/MBVG_AB_MIRA.fasta: $(MOLLI)/MBVG/AIN_ABySS_35-scaffolds.fasta $(MOLLI)/MBVG/AIN_step2_out.unpadded.fasta
# temp_assemblies/MBVG_CLC_MIRA.fasta: $(MOLLI)/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MBVG/AIN_step2_out.unpadded.fasta

# result_statistics/MBVG_quast: \
# 	$(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta
# 	$(MOLLIGAM)/MBVG/GAM_abyss-CLC.fasta $(MOLLIGAM)/MBVG/GAM_CLC-mira.fasta $(MOLLIGAM)/MBVG/GAM_mira-abyss.fasta \
# 	result_assemblies/MBVG_AB_CLC_mix.fasta result_assemblies/MBVG_AB_MIRA_mix.fasta result_assemblies/MBVG_CLC_MIRA_mix.fasta 
# 	rm -rf result_statistics/MBVG_quast
# 	$(pyinterp) $(QUAST) -o $@ $^ 


#VPATH = src:temp_assemblies:result_assemblies

# GENERIC Rule for mollicutes








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



# AllpathSoapBambus_v2.coords: ../contigs.AllpathSoapBambus.fa
# 	cd ..
# 	MUMmer3.23/nucmer --maxmatch -l 30 -banded contigs.AllpathSoapBambus.fa contigs.AllpathSoapBambus.fa 2> /dev/null
# 	MUMmer3.23/show-coords -l -c out.delta > AllpathSoapBambus_v2.coords 

# contig_in_path: ../AllpathSoapBambus_v2.coords Mix.py
# 	# Get the stripped coords file 
# 	./nucmer_coords_filter.py ../AllpathSoapBambus_v2.coords A1_ALLPATH.151 A2_SOAP.scaffold8.1 A3_BAMBUS.scf9_0 > ../TEST_contig_in_path.coords
# 	# Display original input
# 	./display_coords.py ../TEST_contig_in_path.coords
# 	# Execute mix
# 	$(pyinterp) Mix.py -r -a ../TEST_contig_in_path.coords -o TEST_contig_in_path -c ../contigs.AllpathSoapBambus.fa -A 125 -C 0
# 	# Check the resulting assembly with nucmer 
# 	# ../MUMmer3.23/nucmer -p "TEST_contig_in_path_assembly" --maxmatch -l 30 -banded TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_contig_in_path_assembly.delta > TEST_contig_in_path_assembly.coords
# 	# ./display_coords.py TEST_contig_in_path_assembly.coords


# contig_in_path2: ../AllpathSoapBambus_v2.coords Mix.py
# 	# Get the stripped coords file 
# 	./nucmer_coords_filter.py ../AllpathSoapBambus_v2.coords A1_ALLPATH.151 A2_SOAP.scaffold8.1 A2_SOAP.scaffold16.6 A2_SOAP.scaffold17.8 A1_ALLPATH.121 A1_ALLPATH.184 A3_BAMBUS.scf9_0 > ../TEST_contig_in_path.coords
# 	# Display original input
# 	./display_coords.py ../TEST_contig_in_path.coords
# 	# Execute mix
# 	$(pyinterp) Mix.py -r -a ../TEST_contig_in_path.coords -o TEST_contig_in_path -c ../contigs.AllpathSoapBambus.fa -A 125 -C 0
# 	# Check the resulting assembly with nucmer 
# 	../MUMmer3.23/nucmer -p "TEST_contig_in_path_assembly" --maxmatch -l 30 -banded TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	../MUMmer3.23/show-coords -l -c TEST_contig_in_path_assembly.delta > TEST_contig_in_path_assembly.coords
# 	./display_coords.py TEST_contig_in_path_assembly.coords

# smallDuplication: ../AllpathSoapBambus_v2.coords Mix.py
# 	# Get the stripped coords file 
# 	./nucmer_coords_filter.py ../AllpathSoapBambus_v2.coords A2_SOAP.C1141.1 A1_ALLPATH.73 A2_SOAP.scaffold8.9 > ../TEST_contig_in_path.coords
# 	# Display original input
# 	./display_coords.py ../TEST_contig_in_path.coords
# 	# Execute mix
# 	$(pyinterp) Mix.py -r -a ../TEST_contig_in_path.coords -o TEST_contig_in_path -c ../contigs.AllpathSoapBambus.fa -A 125 -C 0
# 	# Check the resulting assembly with nucmer 
# 	../MUMmer3.23/nucmer -p "TEST_contig_in_path_assembly" --maxmatch -l 30 -banded TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta TEST_contig_in_path/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	../MUMmer3.23/show-coords -l -c TEST_contig_in_path_assembly.delta > TEST_contig_in_path_assembly.coords
# 	./display_coords.py TEST_contig_in_path_assembly.coords

# rhodoMix:
# 	$(pyinterp) Mix.py -A 125 -C 0 -o Rhodov2  -a ../AllpathSoapBambus_v2.coords -c ../contigs.AllpathSoapBambus.fa
# 	sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_rhodo_mix_A125_AP_SP_BB.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_rhodov2" --maxmatch -l 30 -banded Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_rhodov2.delta > TEST_rhodov2.coords
# 	# ./nucmer_coords_filter.py TEST_rhodov2.coords -p -s "*" 



# rhodoMix_AP_BB:
# 	./preprocessing.py ../Rhodobacter_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta  -o ../rhodo_AP_BB.fasta
# 	../MUMmer3.23/nucmer -p "rhodoMix_AP_BB" --maxmatch -l 30 -banded ../rhodo_AP_BB.fasta ../rhodo_AP_BB.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c rhodoMix_AP_BB.delta > ../rhodoMix_AP_BB.coords
# 	rm rhodoMix_AP_BB.delta

# 	$(pyinterp) Mix.py -A 500 -C 0 -o rhodoMix_AP_BB  -a ../rhodoMix_AP_BB.coords -c ../rhodo_AP_BB.fasta
# 	# export PATH=$PATH:~/Documents/MIX/MUMmer3.23/
# 	# sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta rhodoMix_AP_BB/Mix_results_A125_C0/Mix_assembly.fasta rhodoMix_AP_BB/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_rhodo_mix_A125_AP_BB.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_rhodov2" --maxmatch -l 30 -banded Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_rhodov2.delta > TEST_rhodov2.coords
# 	# ./nucmer_coords_filter.py TEST_rhodov2.coords -p -s "*" 

# rhodoMix_BB_SP:
# 	./preprocessing.py ../Rhodobacter_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta  -o ../rhodo_BB_SP.fasta
# 	../MUMmer3.23/nucmer -p "rhodoMix_BB_SP" --maxmatch -l 30 -banded ../rhodo_BB_SP.fasta ../rhodo_BB_SP.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c rhodoMix_BB_SP.delta > ../rhodoMix_BB_SP.coords
# 	rm rhodoMix_BB_SP.delta

# 	$(pyinterp) Mix.py -A 500 -C 0 -o rhodoMix_BB_SP  -a ../rhodoMix_BB_SP.coords -c ../rhodo_BB_SP.fasta
# 	# export PATH=$PATH:~/Documents/MIX/MUMmer3.23/
# 	# sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta rhodoMix_BB_SP/Mix_results_A125_C0/Mix_assembly.fasta rhodoMix_BB_SP/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_rhodo_mix_A125_BB_SP.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_rhodov2" --maxmatch -l 30 -banded Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_rhodov2.delta > TEST_rhodov2.coords
# 	# ./nucmer_coords_filter.py TEST_rhodov2.coords -p -s "*" 

# rhodoMix_AP_BB_SP:
# 	./preprocessing.py ../Rhodobacter_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Rhodobacter_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta  -o ../rhodo_AP_BB_SP.fasta
# 	../MUMmer3.23/nucmer -p "rhodoMix_AP_BB_SP" --maxmatch -l 30 -banded ../rhodo_AP_BB_SP.fasta ../rhodo_AP_BB_SP.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c rhodoMix_AP_BB_SP.delta > ../rhodoMix_AP_BB_SP.coords
# 	rm rhodoMix_AP_BB_SP.delta

# 	$(pyinterp) Mix.py -A 500 -C 0 -o rhodoMix_AP_BB_SP  -a ../rhodoMix_AP_BB_SP.coords -c ../rhodo_AP_BB_SP.fasta
# 	# export PATH=$PATH:~/Documents/MIX/MUMmer3.23/
# 	# sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta rhodoMix_BB_SP/Mix_results_A125_C0/Mix_assembly.fasta rhodoMix_BB_SP/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_rhodo_mix_A125_BB_SP.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_rhodov2" --maxmatch -l 30 -banded Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_rhodov2.delta > TEST_rhodov2.coords
# 	# ./nucmer_coords_filter.py TEST_rhodov2.coords -p -s "*" 

# rhodoMix_up3: rhodoMix_AP_BB_SP rhodoMix_BB_SP rhodoMix_AP_BB

# rhodoMix_AP_SP_BB_AB2_CG:
# 	./preprocessing.py ../Rhodobacter_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Rhodobacter_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta  ../Rhodobacter_Assembly/ABySS2/genome.ctg.fasta ../Rhodobacter_Assembly/CABOG/cabog.ctg.fasta -o ../rhodo_AP_SP_BB_AB2_CG.fasta
# 	../MUMmer3.23/nucmer -p "rhodo_AP_SP_BB_AB2_CG" --maxmatch -l 30 -banded ../rhodo_AP_SP_BB_AB2_CG.fasta ../rhodo_AP_SP_BB_AB2_CG.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c rhodo_AP_SP_BB_AB2_CG.delta > ../rhodo_AP_SP_BB_AB2_CG.coords
# 	rm rhodo_AP_SP_BB_AB2_CG.delta

# 	$(pyinterp) Mix.py -A 500 -C 0 -o rhodo_AP_SP_BB_AB2_CG  -a ../rhodo_AP_SP_BB_AB2_CG.coords -c ../rhodo_AP_SP_BB_AB2_CG.fasta
# 	#sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta rhodo_AP_SP_BB_AB2_CG/Mix_results_A125_C0/Mix_assembly.fasta rhodo_AP_SP_BB_AB2_CG/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_rhodo_mix_AP_SP_BB_AB2_CG.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_rhodov2" --maxmatch -l 30 -banded Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta Rhodov2/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_rhodov2.delta > TEST_rhodov2.coords
# 	# ./nucmer_coords_filter.py TEST_rhodov2.coords -p -s "*" 



# Aureus_AP_SOAP_BB.fasta: ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Aureus_Assembly/Bambus2/genome.ctg.fasta
# 	./preprocessing.py ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Aureus_Assembly/Bambus2/genome.ctg.fasta  -o ../Aureus_AP_SOAP_BB.fasta
# Aureus_AP_SOAP_BB.coords:../Aureus_AP_SOAP_BB.fas/gage_rhodo_mix.txt
# 	../MUMmer3.23/nucmer -p "Aureus_AP_SOAP_BB" --maxmatch -l 30 -banded ../Aureus_AP_SOAP_BB.fasta ../Aureus_AP_SOAP_BB.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c Aureus_AP_SOAP_BB.delta > ../Aureus_AP_SOAP_BB.coords
# 	rm Aureus_AP_SOAP_BB.delta

# AureusMix_AP_SP_BB:
# 	$(pyinterp) Mix.py -A 125 -C 0 -o Aureus  -a ../Aureus_AP_SOAP_BB.coords -c ../Aureus_AP_SOAP_BB.fasta
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta Aureus/Mix_results_A125_C0/Mix_assembly.fasta Aureus/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_Aureus_mix_A125_AP_SOAP_BB.txt
# 	# ../MUMmer3.23/nucmer -p "TEST_aerieus" --maxmatch -l 30 -banded Aureus/Mix_results_A125_C0/Mix_assembly.fasta Aureus/Mix_results_A125_C0/Mix_assembly.fasta 2> /dev/null
# 	# ../MUMmer3.23/show-coords -l -c TEST_aerieus.delta > TEST_aerieus.coords
# 	# ./nucmer_coords_filter.py -p -s TEST_aerieus.coords "*" 

# AureusMix_AP_BB:
# 	./preprocessing.py ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/Bambus2/BB_genome.ctg.fasta  -o ../Aureus_AP_BB.fasta
# 	../MUMmer3.23/nucmer -p "Aureus_AP_BB" --maxmatch -l 30 -banded ../Aureus_AP_BB.fasta ../Aureus_AP_BB.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c Aureus_AP_BB.delta > ../Aureus_AP_BB.coords
# 	rm Aureus_AP_BB.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o AureusMix_AP_BB  -a ../Aureus_AP_BB.coords -c ../Aureus_AP_BB.fasta
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta AureusMix_AP_BB/Mix_results_A125_C0/Mix_assembly.fasta AureusMix_AP_BB/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_Aureus_mix_A125_AP_BB.txt

# AureusMix_AP_SP:
# 	./preprocessing.py ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta  -o ../Aureus_AP_SP.fasta
# 	../MUMmer3.23/nucmer -p "Aureus_AP_SP" --maxmatch -l 30 -banded ../Aureus_AP_SP.fasta ../Aureus_AP_SP.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c Aureus_AP_SP.delta > ../Aureus_AP_SP.coords
# 	rm Aureus_AP_SP.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o AureusMix_AP_SP  -a ../Aureus_AP_SP.coords -c ../Aureus_AP_SP.fasta
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta AureusMix_AP_SP/Mix_results_A125_C0/Mix_assembly.fasta AureusMix_AP_SP/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_Aureus_mix_A125_AP_SP.txt

# AureusMix_AP_AB:
# 	./preprocessing.py ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/ABySS2/AB_genome.ctg.fasta  -o ../Aureus_AP_AB.fasta
# 	../MUMmer3.23/nucmer -p "Aureus_AP_AB" --maxmatch -l 30 -banded ../Aureus_AP_AB.fasta ../Aureus_AP_AB.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c Aureus_AP_AB.delta > ../Aureus_AP_AB.coords
# 	rm Aureus_AP_AB.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o AureusMix_AP_AB  -a ../Aureus_AP_AB.coords -c ../Aureus_AP_AB.fasta
# 	# sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta AureusMix_AP_AB/Mix_results_A125_C0/Mix_assembly.fasta AureusMix_AP_AB/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_Aureus_mix_A125_AP_AB.txt

# AureusMix_AP_AB_SP:
# 	./preprocessing.py ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/ABySS2/AB_genome.ctg.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta -o ../Aureus_AP_AB_SP.fasta
# 	../MUMmer3.23/nucmer -p "Aureus_AP_AB_SP" --maxmatch -l 30 -banded ../Aureus_AP_AB_SP.fasta ../Aureus_AP_AB_SP.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c Aureus_AP_AB_SP.delta > ../Aureus_AP_AB_SP.coords
# 	rm Aureus_AP_AB_SP.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o AureusMix_AP_AB_SP  -a ../Aureus_AP_AB_SP.coords -c ../Aureus_AP_AB_SP.fasta
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta AureusMix_AP_AB_SP/Mix_results_A125_C0/Mix_assembly.fasta AureusMix_AP_AB_SP/Mix_results_A125_C0/Mix_assembly.fasta > results/gage_Aureus_mix_A125_AP_AB.txt

# AureusBench:
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta > results/gage_Aureus_allpaths.txt
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta > results/gage_Aureus_soap.txt
# 	sh ../gage/getCorrectnessStats.sh ../Aureus_Assembly/Aureus_genome.fasta ../Aureus_Assembly/Bambus2/genome.ctg.fasta ../Aureus_Assembly/Bambus2/genome.ctg.fasta > results/gage_Aureus_bambus.txt

# rhodoBench:
# 	sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta ../Rhodobacter_Assembly/Allpaths-LG/AP_genome.ctg.fasta ../Rhodobacter_Assembly/Allpaths-LG/AP_genome.ctg.fasta > results/gage_rhodo_allpaths.txt
# 	sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta ../Rhodobacter_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta ../Rhodobacter_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta > results/gage_rhodo_soap.txt
# 	sh ../gage/getCorrectnessStats.sh ../Rhodobacter_Assembly/rho_genome.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta ../Rhodobacter_Assembly/Bambus2/bambus2.ctg.fasta > results/gage_rhodo_bambus.txt



# AIZCLCStep.fasta: 
# 	./preprocessing.py ../AIZ/AIZ_CLC_contigsCLCTrimmed_0.fasta ../AIZ/AIZ_step2_out.unpadded.fasta -o ../AIZCLCStep.fasta

# AIZCLCStep.coords: ../AIZCLCStep.fasta
# 	../MUMmer3.23/nucmer -p "AIZCLCStep" --maxmatch -l 30 -banded ../AIZCLCStep.fasta ../AIZCLCStep.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c AIZCLCStep.delta > ../AIZCLCStep.coords
# 	rm AIZCLCStep.delta

# AIZMix: 
# 	$(pyinterp) Mix.py -a ../AIZ/AIZCLSStep2.coords -o AIZ -c ../AIZCLCStep.fasta -A 125 -C 0



# # All molicutes
# MOVI_CLC_MIRA:
# 	./preprocessing.py ../MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta ../MOVI/AIP_step2_out.unpadded.fasta -o ../MOVI/MOVI_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MOVI_CLC_MIRA" --maxmatch -l 30 -banded ../MOVI/MOVI_CLC_MIRA.fasta ../MOVI/MOVI_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MOVI_CLC_MIRA.delta > ../MOVI/MOVI_CLC_MIRA.coords
# 	rm MOVI_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MOVI_CLC_MIRA  -a ../MOVI/MOVI_CLC_MIRA.coords -c ../MOVI/MOVI_CLC_MIRA.fasta
# 	mv MOVI_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MOVI/MIX_MOVI_CLC_MIRA.fasta

# MOVI:MOVI_CLC_MIRA

# MMC_CLC_MIRA:
# 	./preprocessing.py ../MMC/AIW_CLC_contigsCLCTrimmed_0.fasta ../MMC/AIW_step2_out.unpadded.fasta -o ../MMC/MMC_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MMC_CLC_MIRA" --maxmatch -l 30 -banded ../MMC/MMC_CLC_MIRA.fasta ../MMC/MMC_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MMC_CLC_MIRA.delta > ../MMC/MMC_CLC_MIRA.coords
# 	rm MMC_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MMC_CLC_MIRA  -a ../MMC/MMC_CLC_MIRA.coords -c ../MMC/MMC_CLC_MIRA.fasta
# 	mv MMC_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MMC/MIX_MMC_CLC_MIRA.fasta

# MMC:MMC_CLC_MIRA

# MSCe_CLC_MIRA:
# 	./preprocessing.py ../MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta ../MSCe/AKC_step2_out.unpadded.fasta -o ../MSCe/MSCe_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MSCe_CLC_MIRA" --maxmatch -l 30 -banded ../MSCe/MSCe_CLC_MIRA.fasta ../MSCe/MSCe_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCe_CLC_MIRA.delta > ../MSCe/MSCe_CLC_MIRA.coords
# 	rm MSCe_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCe_CLC_MIRA  -a ../MSCe/MSCe_CLC_MIRA.coords -c ../MSCe/MSCe_CLC_MIRA.fasta
# 	mv MSCe_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCe/MIX_MSCe_CLC_MIRA.fasta

# MSCe:MSCe_CLC_MIRA

# MSCd_CLC_MIRA:
# 	./preprocessing.py ../MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta ../MSCd/AKE_step2_out.unpadded.fasta -o ../MSCd/MSCd_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MSCd_CLC_MIRA" --maxmatch -l 30 -banded ../MSCd/MSCd_CLC_MIRA.fasta ../MSCd/MSCd_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCd_CLC_MIRA.delta > ../MSCd/MSCd_CLC_MIRA.coords
# 	rm MSCd_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCd_CLC_MIRA  -a ../MSCd/MSCd_CLC_MIRA.coords -c ../MSCd/MSCd_CLC_MIRA.fasta
# 	mv MSCd_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCd/MIX_MSCd_CLC_MIRA.fasta

# MSCd: MSCd_CLC_MIRA

# MSCc_CLC_MIRA:
# 	./preprocessing.py ../MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta ../MSCc/AIZ_step2_out.unpadded.fasta -o ../MSCc/MSCc_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MSCc_CLC_MIRA" --maxmatch -l 30 -banded ../MSCc/MSCc_CLC_MIRA.fasta ../MSCc/MSCc_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCc_CLC_MIRA.delta > ../MSCc/MSCc_CLC_MIRA.coords
# 	rm MSCc_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCc_CLC_MIRA  -a ../MSCc/MSCc_CLC_MIRA.coords -c ../MSCc/MSCc_CLC_MIRA.fasta
# 	mv MSCc_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCc/MIX_MSCc_CLC_MIRA.fasta

# MSCc:MSCc_CLC_MIRA

# MSCb_CLC_MIRA:
# 	./preprocessing.py ../MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta ../MSCb/AIY_step2_out.unpadded.fasta -o ../MSCb/MSCb_CLC_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MSCb_CLC_MIRA" --maxmatch -l 30 -banded ../MSCb/MSCb_CLC_MIRA.fasta ../MSCb/MSCb_CLC_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCb_CLC_MIRA.delta > ../MSCb/MSCb_CLC_MIRA.coords
# 	rm MSCb_CLC_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCb_CLC_MIRA  -a ../MSCb/MSCb_CLC_MIRA.coords -c ../MSCb/MSCb_CLC_MIRA.fasta
# 	mv MSCb_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCb/MIX_MSCb_CLC_MIRA.fasta

# MSCb_AB_CLC:
# 	./preprocessing.py ../MSCb/AIY_ABySS_29-scaffolds.fasta ../MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta -o ../MSCb/MSCb_AB_CLC.fasta
# 	../MUMmer3.23/nucmer -p "MSCb_AB_CLC" --maxmatch -l 30 -banded ../MSCb/MSCb_AB_CLC.fasta ../MSCb/MSCb_AB_CLC.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCb_AB_CLC.delta > ../MSCb/MSCb_AB_CLC.coords
# 	rm MSCb_AB_CLC.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCb_AB_CLC  -a ../MSCb/MSCb_AB_CLC.coords -c ../MSCb/MSCb_AB_CLC.fasta
# 	mv MSCb_AB_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCb/MSCb_AB_CLC.fasta


# MSCb_AB_MIRA:
# 	./preprocessing.py ../MSCb/AIY_ABySS_29-scaffolds.fasta ../MSCb/AIY_step2_out.unpadded.fasta -o ../MSCb/MSCb_AB_MIRA.fasta
# 	../MUMmer3.23/nucmer -p "MSCb_AB_MIRA" --maxmatch -l 30 -banded ../MSCb/MSCb_AB_MIRA.fasta ../MSCb/MSCb_AB_MIRA.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCb_AB_MIRA.delta > ../MSCb/MSCb_AB_MIRA.coords
# 	rm MSCb_AB_MIRA.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCb_AB_MIRA  -a ../MSCb/MSCb_AB_MIRA.coords -c ../MSCb/MSCb_AB_MIRA.fasta
# 	mv MSCb_AB_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCb/MSCb_AB_MIRA.fasta

# MSCb:MSCb_CLC_MIRA
# AllMoli: MMC_CLC_MIRA MSCe_CLC_MIRA MSCd_CLC_MIRA MSCc_CLC_MIRA




# MSCe_CLC_MIRA_no91:
# 	./preprocessing.py ../MSCe/AKC_CLC_contigsCLCTrimmed_0_sans_91.fasta ../MSCe/AKC_step2_out.unpadded.fasta -o ../MSCe/MSCe_CLC_MIRA_sans91.fasta
# 	../MUMmer3.23/nucmer -p "MSCe_CLC_MIRA_sans91" --maxmatch -l 30 -banded ../MSCe/MSCe_CLC_MIRA_sans91.fasta ../MSCe/MSCe_CLC_MIRA_sans91.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCe_CLC_MIRA_sans91.delta > ../MSCe/MSCe_CLC_MIRA_sans91.coords
# 	rm MSCe_CLC_MIRA_sans91.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCe_CLC_MIRA_sans91  -a ../MSCe/MSCe_CLC_MIRA_sans91.coords -c ../MSCe/MSCe_CLC_MIRA_sans91.fasta
# 	mv MSCe_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCe/MIX_MSCe_CLC_MIRA.fasta


# MSCe_CLC_MIRA_91_subset:
# 	./preprocessing.py ../MSCe/AKC_CLC_contigsCLCTrimmed_0_sans_91.fasta ../MSCe/AKC_step2_out.unpadded.fasta -o ../MSCe/MSCe_CLC_MIRA_sans91.fasta
# 	../MUMmer3.23/nucmer -p "MSCe_CLC_MIRA_sans91" --maxmatch -l 30 -banded ../MSCe/MSCe_CLC_MIRA_sans91.fasta ../MSCe/MSCe_CLC_MIRA_sans91.fasta 2>/dev/null
# 	../MUMmer3.23/show-coords -l -c MSCe_CLC_MIRA_sans91.delta > ../MSCe/MSCe_CLC_MIRA_sans91.coords
# 	rm MSCe_CLC_MIRA_sans91.delta
# 	$(pyinterp) Mix.py -A 125 -C 0 -o MSCe_CLC_MIRA_sans91  -a ../MSCe/MSCe_CLC_MIRA_sans91.coords -c ../MSCe/MSCe_CLC_MIRA_sans91.fasta
# 	mv MSCe_CLC_MIRA/Mix_results_A125_C0/Mix_assembly.fasta ../MSCe/MIX_MSCe_CLC_MIRA.fasta


# AureusQuast:
# 	python ../quast-2.1/quast.py -R ../Aureus_Assembly/Aureus_genome.fasta \
# 	../Mix-1.1/AureusMix_AP_BB/Mix_results_A125_C0/Mix_AP_BB_assembly.fasta \
# 	../Mix-1.1/AureusMix_AP_SP/Mix_results_A125_C0/Mix_AP_SP_assembly.fasta \
# 	../Mix-1.1/AureusMix_AP_AB/Mix_results_A125_C0/Mix_AP_AB_assembly.fasta \
# 	../Mix-1.1/AureusMix_AP_AB_SP/Mix_results_A125_C0/Mix_AP_AB_SP_assembly.fasta \
# 	../Aureus_Assembly/Allpaths-LG/AP_genome.ctg.fasta \
# 	../Aureus_Assembly/SOAPdenovo/SOAP_genome.ctg.fasta \
# 	../Aureus_Assembly/ABySS2/AB_genome.ctg.fasta \
# 	../Aureus_Assembly/Bambus2/BB_genome.ctg.fasta



# B_cereus_GAGE_eval: 
# 	sh ../gage/getCorrectnessStats.sh ../B_cereus_MiSeq/B_cereus_NC_003909.8_ref_genome.fasta ../B_cereus_MiSeq/spades_ctg.fasta ../B_cereus_MiSeq/spades_ctg.fasta > results/gage_B_B_cereus_SPADES.txt

