MOLLI:=datasets/Mollicutes/
MOLLIGAM:=datasets/GAM-NGS/Mollicutes

# The smallest mollicute protein has 36 aa, 108 bp. Anything below should be considered as noise
QUASTMOLLI:= $(QUAST) --gene-finding --gene-thresholds 0,100,300,500,1000,1500,3000

# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta -o result_assemblies/MBVG_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MAUR/AIM_CLC_contigsCLCTrimmed_0.fasta -o result_assemblies/MAUR_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mollicutes/MBOVb/AIQ_CLC_contigsCLCTrimmed_0.fasta -o result_assemblies/MBOVb_CLC_naive_concatenation.fasta


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


temp_assemblies/MAUR_AB_CLC.fasta: $(MOLLI)/MAUR/AIM_ABySS_27-scaffolds.fasta $(MOLLI)/MAUR/AIM_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MAUR_AB_MIRA.fasta: $(MOLLI)/MAUR/AIM_ABySS_27-scaffolds.fasta $(MOLLI)/MAUR/AIM_step2_out.unpadded.fasta
temp_assemblies/MAUR_CLC_MIRA.fasta: $(MOLLI)/MAUR/AIM_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MAUR/AIM_step2_out.unpadded.fasta

result_statistics/MAUR_quast: \
	$(MOLLI)/MAUR/AIM_ABySS_27-scaffolds.fasta $(MOLLI)/MAUR/AIM_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MAUR/AIM_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MAUR/GAM_abyss-CLC.fasta $(MOLLIGAM)/MAUR/GAM_CLC-mira.fasta $(MOLLIGAM)/MAUR/GAM_mira-abyss.fasta \
	result_assemblies/MAUR_AB_CLC_mix.fasta result_assemblies/MAUR_AB_MIRA_mix.fasta result_assemblies/MAUR_CLC_MIRA_mix.fasta \
	result_assemblies/MAUR_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MAUR_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBOVb_AB_CLC.fasta: $(MOLLI)/MBOVb/AIQ_ABySS_34-scaffolds.fasta $(MOLLI)/MBOVb/AIQ_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MBOVb_AB_MIRA.fasta: $(MOLLI)/MBOVb/AIQ_ABySS_34-scaffolds.fasta $(MOLLI)/MBOVb/AIQ_step2_out.unpadded.fasta
temp_assemblies/MBOVb_CLC_MIRA.fasta: $(MOLLI)/MBOVb/AIQ_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MBOVb/AIQ_step2_out.unpadded.fasta

result_statistics/MBOVb_quast: \
	$(MOLLI)/MBOVb/AIQ_ABySS_34-scaffolds.fasta $(MOLLI)/MBOVb/AIQ_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MBOVb/AIQ_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MBOVb/GAM_abyss-CLC.fasta $(MOLLIGAM)/MBOVb/GAM_CLC-mira.fasta $(MOLLIGAM)/MBOVb/GAM_mira-abyss.fasta \
	result_assemblies/MBOVb_AB_CLC_mix.fasta result_assemblies/MBOVb_AB_MIRA_mix.fasta result_assemblies/MBOVb_CLC_MIRA_mix.fasta \
	result_assemblies/MBOVb_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MBOVb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBVG_AB_CLC.fasta: $(MOLLI)/MBVG/AIN_ABySS_35-scaffolds.fasta $(MOLLI)/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MBVG_AB_MIRA.fasta: $(MOLLI)/MBVG/AIN_ABySS_35-scaffolds.fasta $(MOLLI)/MBVG/AIN_step2_out.unpadded.fasta
temp_assemblies/MBVG_CLC_MIRA.fasta: $(MOLLI)/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MBVG/AIN_step2_out.unpadded.fasta

result_statistics/MBVG_quast: \
	$(MOLLI)/MBVG/AIN_ABySS_35-scaffolds.fasta $(MOLLI)/MBVG/AIN_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MBVG/AIN_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MBVG/GAM_abyss-CLC.fasta $(MOLLIGAM)/MBVG/GAM_CLC-mira.fasta $(MOLLIGAM)/MBVG/GAM_mira-abyss.fasta \
	result_assemblies/MBVG_AB_CLC_mix.fasta result_assemblies/MBVG_AB_MIRA_mix.fasta result_assemblies/MBVG_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MBVG_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MCCP_AB_CLC.fasta: $(MOLLI)/MCCP/AIX_ABySS_35-scaffolds.fasta $(MOLLI)/MCCP/AIX_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MCCP_AB_MIRA.fasta: $(MOLLI)/MCCP/AIX_ABySS_35-scaffolds.fasta $(MOLLI)/MCCP/AIX_step2_out.unpadded.fasta
temp_assemblies/MCCP_CLC_MIRA.fasta: $(MOLLI)/MCCP/AIX_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MCCP/AIX_step2_out.unpadded.fasta

result_statistics/MCCP_quast: \
	$(MOLLI)/MCCP/AIX_ABySS_35-scaffolds.fasta $(MOLLI)/MCCP/AIX_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MCCP/AIX_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MCCP/GAM_abyss-CLC.fasta $(MOLLIGAM)/MCCP/GAM_CLC-mira.fasta $(MOLLIGAM)/MCCP/GAM_mira-abyss.fasta \
	result_assemblies/MCCP_AB_CLC_mix.fasta result_assemblies/MCCP_AB_MIRA_mix.fasta result_assemblies/MCCP_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MCCP_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MOVI_AB_CLC.fasta: $(MOLLI)/MOVI/AIP_ABySS_32-scaffolds.fasta $(MOLLI)/MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MOVI_AB_MIRA.fasta: $(MOLLI)/MOVI/AIP_ABySS_32-scaffolds.fasta $(MOLLI)/MOVI/AIP_step2_out.unpadded.fasta
temp_assemblies/MOVI_CLC_MIRA.fasta: $(MOLLI)/MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MOVI/AIP_step2_out.unpadded.fasta

result_statistics/MOVI_quast: \
	$(MOLLI)/MOVI/AIP_ABySS_32-scaffolds.fasta $(MOLLI)/MOVI/AIP_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MOVI/AIP_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MOVI/GAM_abyss-CLC.fasta $(MOLLIGAM)/MOVI/GAM_CLC-mira.fasta $(MOLLIGAM)/MOVI/GAM_mira-abyss.fasta \
	result_assemblies/MOVI_AB_CLC_mix.fasta result_assemblies/MOVI_AB_MIRA_mix.fasta result_assemblies/MOVI_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MOVI_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MMC_AB_CLC.fasta: $(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MMC_AB_MIRA.fasta: $(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta
temp_assemblies/MMC_CLC_MIRA.fasta: $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta

result_statistics/MMC_quast: \
	$(MOLLI)/MMC/AIW_ABySS_30-scaffolds.fasta $(MOLLI)/MMC/AIW_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MMC/AIW_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MMC/GAM_abyss-CLC.fasta $(MOLLIGAM)/MMC/GAM_CLC-mira.fasta $(MOLLIGAM)/MMC/GAM_mira-abyss.fasta \
	result_assemblies/MMC_AB_CLC_mix.fasta result_assemblies/MMC_AB_MIRA_mix.fasta result_assemblies/MMC_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MMC_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCe_AB_CLC.fasta: $(MOLLI)/MSCe/AKC_ABySS_35-scaffolds.fasta $(MOLLI)/MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MSCe_AB_MIRA.fasta: $(MOLLI)/MSCe/AKC_ABySS_35-scaffolds.fasta $(MOLLI)/MSCe/AKC_step2_out.unpadded.fasta
temp_assemblies/MSCe_CLC_MIRA.fasta: $(MOLLI)/MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCe/AKC_step2_out.unpadded.fasta

result_statistics/MSCe_quast: \
	$(MOLLI)/MSCe/AKC_ABySS_35-scaffolds.fasta $(MOLLI)/MSCe/AKC_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCe/AKC_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MSCe/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCe/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCe/GAM_mira-abyss.fasta \
	result_assemblies/MSCe_AB_CLC_mix.fasta result_assemblies/MSCe_AB_MIRA_mix.fasta result_assemblies/MSCe_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCe_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCd_AB_CLC.fasta: $(MOLLI)/MSCd/AKE_ABySS_34-scaffolds.fasta $(MOLLI)/MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MSCd_AB_MIRA.fasta: $(MOLLI)/MSCd/AKE_ABySS_34-scaffolds.fasta $(MOLLI)/MSCd/AKE_step2_out.unpadded.fasta
temp_assemblies/MSCd_CLC_MIRA.fasta: $(MOLLI)/MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCd/AKE_step2_out.unpadded.fasta

result_statistics/MSCd_quast: \
	$(MOLLI)/MSCd/AKE_ABySS_34-scaffolds.fasta $(MOLLI)/MSCd/AKE_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCd/AKE_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MSCd/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCd/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCd/GAM_mira-abyss.fasta \
	result_assemblies/MSCd_AB_CLC_mix.fasta result_assemblies/MSCd_AB_MIRA_mix.fasta result_assemblies/MSCd_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCd_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCc_AB_CLC.fasta: $(MOLLI)/MSCc/AIZ_ABySS_27-scaffolds.fasta $(MOLLI)/MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MSCc_AB_MIRA.fasta: $(MOLLI)/MSCc/AIZ_ABySS_27-scaffolds.fasta $(MOLLI)/MSCc/AIZ_step2_out.unpadded.fasta
temp_assemblies/MSCc_CLC_MIRA.fasta: $(MOLLI)/MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCc/AIZ_step2_out.unpadded.fasta

result_statistics/MSCc_quast: \
	$(MOLLI)/MSCc/AIZ_ABySS_27-scaffolds.fasta $(MOLLI)/MSCc/AIZ_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCc/AIZ_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MSCc/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCc/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCc/GAM_mira-abyss.fasta \
	result_assemblies/MSCc_AB_CLC_mix.fasta result_assemblies/MSCc_AB_MIRA_mix.fasta result_assemblies/MSCc_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCc_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCb_AB_CLC.fasta: $(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta
temp_assemblies/MSCb_AB_MIRA.fasta: $(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta
temp_assemblies/MSCb_CLC_MIRA.fasta: $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta

result_statistics/MSCb_quast: \
	$(MOLLI)/MSCb/AIY_ABySS_29-scaffolds.fasta $(MOLLI)/MSCb/AIY_CLC_contigsCLCTrimmed_0.fasta $(MOLLI)/MSCb/AIY_step2_out.unpadded.fasta \
	$(MOLLIGAM)/MSCb/GAM_abyss-CLC.fasta $(MOLLIGAM)/MSCb/GAM_CLC-mira.fasta $(MOLLIGAM)/MSCb/GAM_mira-abyss.fasta \
	result_assemblies/MSCb_AB_CLC_mix.fasta result_assemblies/MSCb_AB_MIRA_mix.fasta result_assemblies/MSCb_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^