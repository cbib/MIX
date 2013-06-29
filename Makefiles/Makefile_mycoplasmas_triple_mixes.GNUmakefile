MOLLI:=datasets/Mycoplasmas/
QUASTMOLLI:= $(QUAST) --gene-finding --gene-thresholds 0,100,300,500,1000,1500,3000

temp_assemblies/MSCD_triple_assemblies.fasta: 	$(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_CLC.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta
temp_assemblies/MSCc_triple_assemblies.fasta: 	$(MOLLI)/MSCc/MSCc_ABySS_27.fasta $(MOLLI)/MSCc/MSCc_CLC.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta
temp_assemblies/MSCB_triple_assemblies.fasta: 	$(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_CLC.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta
temp_assemblies/MOVI_triple_assemblies.fasta: 	$(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_CLC.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta
temp_assemblies/MMC_triple_assemblies.fasta: 	$(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_CLC.fasta $(MOLLI)/MMC/MMC_MIRA.fasta
temp_assemblies/MCCP_triple_assemblies.fasta: 	$(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_CLC.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta
temp_assemblies/MBVG_triple_assemblies.fasta: 	$(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_CLC.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta
temp_assemblies/MBOVb_triple_assemblies.fasta: $(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_CLC.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta
temp_assemblies/MAUR_triple_assemblies.fasta: 	$(MOLLI)/MAUR/MAUR_ABySS_27.fasta $(MOLLI)/MAUR/MAUR_CLC.fasta $(MOLLI)/MAUR/MAUR_MIRA.fasta
temp_assemblies/MSCE_triple_assemblies.fasta: 	$(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_CLC.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta


result_statistics/Mycoplasmas-All-triples: \
	result_assemblies/MSCD_triple_assemblies_mix.fasta \
	result_assemblies/MSCc_triple_assemblies_mix.fasta \
	result_assemblies/MSCB_triple_assemblies_mix.fasta \
	result_assemblies/MOVI_triple_assemblies_mix.fasta \
	result_assemblies/MMC_triple_assemblies_mix.fasta \
	result_assemblies/MCCP_triple_assemblies_mix.fasta \
	result_assemblies/MBVG_triple_assemblies_mix.fasta \
	result_assemblies/MBOVb_triple_assemblies_mix.fasta \
	result_assemblies/MAUR_triple_assemblies_mix.fasta \
	result_assemblies/MSCE_triple_assemblies_mix.fasta 
	$(pyinterp) $(QUASTMOLLI) -o $@ $^	