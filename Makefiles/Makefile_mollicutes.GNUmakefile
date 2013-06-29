MOLLI:=datasets/Mollicutes/
MOLLIGAM:=datasets/GAM-NGS/Mollicutes

# The smallest mollicute protein has 36 aa, 108 bp. Anything below should be considered as noise
QUASTMOLLI:= $(QUAST) --gene-finding --gene-thresholds 0,100,300,500,1000,1500,3000

# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MBVG/MBVG_CLC.fasta -o result_assemblies/MBVG_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MAUR/MAUR_CLC.fasta -o result_assemblies/MAUR_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MBOVb/MBOVb_CLC.fasta -o result_assemblies/MBOVb_CLC_naive_concatenation.fasta


MolliMix: \
	result_assemblies/MAUR_AB_CLC_mix.fasta result_assemblies/MAUR_AB_MIRA_mix.fasta result_assemblies/MAUR_CLC_MIRA_mix.fasta \
	result_assemblies/MBOVb_AB_CLC_mix.fasta result_assemblies/MBOVb_AB_MIRA_mix.fasta result_assemblies/MBOVb_CLC_MIRA_mix.fasta \
	result_assemblies/MBVG_AB_CLC_mix.fasta result_assemblies/MBVG_AB_MIRA_mix.fasta result_assemblies/MBVG_CLC_MIRA_mix.fasta \
	result_assemblies/MCCP_AB_CLC_mix.fasta result_assemblies/MCCP_AB_MIRA_mix.fasta result_assemblies/MCCP_CLC_MIRA_mix.fasta \
	result_assemblies/MOVI_AB_CLC_mix.fasta result_assemblies/MOVI_AB_MIRA_mix.fasta result_assemblies/MOVI_CLC_MIRA_mix.fasta \
	result_assemblies/MMC_AB_CLC_mix.fasta result_assemblies/MMC_AB_MIRA_mix.fasta result_assemblies/MMC_CLC_MIRA_mix.fasta \
	result_assemblies/MSCe_AB_CLC_mix.fasta result_assemblies/MSCe_AB_MIRA_mix.fasta result_assemblies/MSCe_CLC_MIRA_mix.fasta \
	result_assemblies/MSCd_AB_CLC_mix.fasta result_assemblies/MSCd_AB_MIRA_mix.fasta result_assemblies/MSCd_CLC_MIRA_mix.fasta \
	result_assemblies/MSCc_AB_CLC_mix.fasta result_assemblies/MSCc_AB_MIRA_mix.fasta result_assemblies/MSCc_CLC_MIRA_mix.fasta \
	result_assemblies/MSCb_AB_CLC_mix.fasta result_assemblies/MSCb_AB_MIRA_mix.fasta result_assemblies/MSCb_CLC_MIRA_mix.fasta

MolliQuast: result_statistics/MAUR_quast result_statistics/MBOVb_quast result_statistics/MBVG_quast result_statistics/MCCP_quast result_statistics/MOVI_quast result_statistics/MMC_quast result_statistics/MSCe_quast result_statistics/MSCd_quast result_statistics/MSCc_quast result_statistics/MSCb_quast
	
temp_assemblies/MAUR_AB_CLC.fasta: $(MOLLI)/MAUR/MAUR_ABySS_27.fasta $(MOLLI)/MAUR/MAUR_CLC.fasta
temp_assemblies/MAUR_AB_MIRA.fasta: $(MOLLI)/MAUR/MAUR_ABySS_27.fasta $(MOLLI)/MAUR/MAUR_MIRA.fasta
temp_assemblies/MAUR_CLC_MIRA.fasta: $(MOLLI)/MAUR/MAUR_CLC.fasta $(MOLLI)/MAUR/MAUR_MIRA.fasta

result_statistics/MAUR_quast: \
	$(MOLLI)/MAUR/MAUR_ABySS_27.fasta $(MOLLI)/MAUR/MAUR_CLC.fasta $(MOLLI)/MAUR/MAUR_MIRA.fasta \
	$(MOLLIGAM)/MAUR/GAM_abyss-CLC.fasta $(MOLLIGAM)/MAUR/GAM_CLC-mira.fasta $(MOLLIGAM)/MAUR/GAM_mira-abyss.fasta \
	result_assemblies/MAUR_AB_CLC_mix.fasta result_assemblies/MAUR_AB_MIRA_mix.fasta result_assemblies/MAUR_CLC_MIRA_mix.fasta \
	result_assemblies/MAUR_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MAUR_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBOVb_AB_CLC.fasta: $(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_CLC.fasta
temp_assemblies/MBOVb_AB_MIRA.fasta: $(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta
temp_assemblies/MBOVb_CLC_MIRA.fasta: $(MOLLI)/MBOVb/MBOVb_CLC.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta

result_statistics/MBOVb_quast: \
	$(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_CLC.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta \
	$(MOLLIGAM)/MBOVb/GAM_abyss-CLC.fasta $(MOLLIGAM)/MBOVb/GAM_CLC-mira.fasta $(MOLLIGAM)/MBOVb/GAM_mira-abyss.fasta \
	result_assemblies/MBOVb_AB_CLC_mix.fasta result_assemblies/MBOVb_AB_MIRA_mix.fasta result_assemblies/MBOVb_CLC_MIRA_mix.fasta \
	result_assemblies/MBOVb_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MBOVb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBVG_AB_CLC.fasta: $(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_CLC.fasta
temp_assemblies/MBVG_AB_MIRA.fasta: $(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta
temp_assemblies/MBVG_CLC_MIRA.fasta: $(MOLLI)/MBVG/MBVG_CLC.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta

result_statistics/MBVG_quast: \
	$(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_CLC.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta \
	$(MOLLIGAM)/MBVG/GAM_abyss-CLC.fasta $(MOLLIGAM)/MBVG/GAM_CLC-mira.fasta $(MOLLIGAM)/MBVG/GAM_mira-abyss.fasta \
	result_assemblies/MBVG_AB_CLC_mix.fasta result_assemblies/MBVG_AB_MIRA_mix.fasta result_assemblies/MBVG_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MBVG_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MCCP_AB_CLC.fasta: $(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_CLC.fasta
temp_assemblies/MCCP_AB_MIRA.fasta: $(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta
temp_assemblies/MCCP_CLC_MIRA.fasta: $(MOLLI)/MCCP/MCCP_CLC.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta

result_statistics/MCCP_quast: \
	$(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_CLC.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta \
	$(MOLLIGAM)/MCCP/GAM_abyss-CLC.fasta $(MOLLIGAM)/MCCP/GAM_CLC-mira.fasta $(MOLLIGAM)/MCCP/GAM_mira-abyss.fasta \
	result_assemblies/MCCP_AB_CLC_mix.fasta result_assemblies/MCCP_AB_MIRA_mix.fasta result_assemblies/MCCP_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MCCP_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MOVI_AB_CLC.fasta: $(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_CLC.fasta
temp_assemblies/MOVI_AB_MIRA.fasta: $(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta
temp_assemblies/MOVI_CLC_MIRA.fasta: $(MOLLI)/MOVI/MOVI_CLC.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta

result_statistics/MOVI_quast: \
	$(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_CLC.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta \
	$(MOLLIGAM)/MOVI/GAM_abyss-CLC.fasta $(MOLLIGAM)/MOVI/GAM_CLC-mira.fasta $(MOLLIGAM)/MOVI/GAM_mira-abyss.fasta \
	result_assemblies/MOVI_AB_CLC_mix.fasta result_assemblies/MOVI_AB_MIRA_mix.fasta result_assemblies/MOVI_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MOVI_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MMC_AB_CLC.fasta: $(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_CLC.fasta
temp_assemblies/MMC_AB_MIRA.fasta: $(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_MIRA.fasta
temp_assemblies/MMC_CLC_MIRA.fasta: $(MOLLI)/MMC/MMC_CLC.fasta $(MOLLI)/MMC/MMC_MIRA.fasta

result_statistics/MMC_quast: \
	$(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_CLC.fasta $(MOLLI)/MMC/MMC_MIRA.fasta \
	$(MOLLIGAM)/MMC/GAM_abyss-CLC.fasta $(MOLLIGAM)/MMC/GAM_CLC-mira.fasta $(MOLLIGAM)/MMC/GAM_mira-abyss.fasta \
	result_assemblies/MMC_AB_CLC_mix.fasta result_assemblies/MMC_AB_MIRA_mix.fasta result_assemblies/MMC_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MMC_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCe_AB_CLC.fasta: $(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_CLC.fasta
temp_assemblies/MSCe_AB_MIRA.fasta: $(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta
temp_assemblies/MSCe_CLC_MIRA.fasta: $(MOLLI)/MSCe/MSCe_CLC.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta

result_statistics/MSCe_quast: \
	$(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_CLC.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta \
	$(MOLLIGAM)/MSCe/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCe/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCe/GAM_mira-abyss.fasta \
	result_assemblies/MSCe_AB_CLC_mix.fasta result_assemblies/MSCe_AB_MIRA_mix.fasta result_assemblies/MSCe_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCe_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCd_AB_CLC.fasta: $(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_CLC.fasta
temp_assemblies/MSCd_AB_MIRA.fasta: $(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta
temp_assemblies/MSCd_CLC_MIRA.fasta: $(MOLLI)/MSCd/MSCd_CLC.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta

result_statistics/MSCd_quast: \
	$(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_CLC.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta \
	$(MOLLIGAM)/MSCd/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCd/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCd/GAM_mira-abyss.fasta \
	result_assemblies/MSCd_AB_CLC_mix.fasta result_assemblies/MSCd_AB_MIRA_mix.fasta result_assemblies/MSCd_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCd_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCc_AB_CLC.fasta: $(MOLLI)/MSCc/MSCc_ABySS_27 $(MOLLI)/MSCc/MSCc_CLC.fasta
temp_assemblies/MSCc_AB_MIRA.fasta: $(MOLLI)/MSCc/MSCc_ABySS_27 $(MOLLI)/MSCc/MSCc_MIRA.fasta
temp_assemblies/MSCc_CLC_MIRA.fasta: $(MOLLI)/MSCc/MSCc_CLC.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta

result_statistics/MSCc_quast: \
	$(MOLLI)/MSCc/MSCc_ABySS_27 $(MOLLI)/MSCc/MSCc_CLC.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta \
	$(MOLLIGAM)/MSCc/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCc/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCc/GAM_mira-abyss.fasta \
	result_assemblies/MSCc_AB_CLC_mix.fasta result_assemblies/MSCc_AB_MIRA_mix.fasta result_assemblies/MSCc_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCc_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCb_AB_CLC.fasta: $(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_CLC.fasta
temp_assemblies/MSCb_AB_MIRA.fasta: $(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta
temp_assemblies/MSCb_CLC_MIRA.fasta: $(MOLLI)/MSCb/MSCb_CLC.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta

result_statistics/MSCb_quast: \
	$(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_CLC.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta \
	$(MOLLIGAM)/MSCb/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCb/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCb/GAM_mira-abyss.fasta \
	result_assemblies/MSCb_AB_CLC_mix.fasta result_assemblies/MSCb_AB_MIRA_mix.fasta result_assemblies/MSCb_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^