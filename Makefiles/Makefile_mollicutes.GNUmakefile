MOLLI:=datasets/Mycoplasmas/
MOLLIGAM:=datasets/GAM-NGS/Mycoplasmas

# The smallest mollicute protein has 36 aa, 108 bp. Anything below should be considered as noise
QUASTMOLLI:= $(QUAST) --gene-finding --gene-thresholds 0,100,300,500,1000,1500,3000

# python bin/assembly_by_naive_concatenation.py datasets/Mycoplasmas/MBVG/MBVG_CLC.fasta -o result_assemblies/MBVG_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mycoplasmas/MAUR/MAUR_CLC.fasta -o result_assemblies/MAUR_CLC_naive_concatenation.fasta
# python bin/assembly_by_naive_concatenation.py datasets/Mycoplasmas/MBOVb/MBOVb_CLC.fasta -o result_assemblies/MBOVb_CLC_naive_concatenation.fasta


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
	$(MOLLIGAM)/MAUR/MAUR_ABySS-CLC_gam.fasta $(MOLLIGAM)/MAUR/MAUR_CLC-mira_gam.fasta $(MOLLIGAM)/MAUR/MAUR_mira-ABySS_gam.fasta \
	result_assemblies/MAUR_AB_CLC_mix.fasta result_assemblies/MAUR_AB_MIRA_mix.fasta result_assemblies/MAUR_CLC_MIRA_mix.fasta \
	result_assemblies/MAUR_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MAUR_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBOVb_AB_CLC.fasta: $(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_CLC.fasta
temp_assemblies/MBOVb_AB_MIRA.fasta: $(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta
temp_assemblies/MBOVb_CLC_MIRA.fasta: $(MOLLI)/MBOVb/MBOVb_CLC.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta

result_statistics/MBOVb_quast: \
	$(MOLLI)/MBOVb/MBOVb_ABySS_34.fasta $(MOLLI)/MBOVb/MBOVb_CLC.fasta $(MOLLI)/MBOVb/MBOVb_MIRA.fasta \
	$(MOLLIGAM)/MBOVb/MBOVb_ABySS-CLC_gam.fasta $(MOLLIGAM)/MBOVb/MBOVb_CLC-mira_gam.fasta $(MOLLIGAM)/MBOVb/MBOVb_mira-ABySS_gam.fasta \
	result_assemblies/MBOVb_AB_CLC_mix.fasta result_assemblies/MBOVb_AB_MIRA_mix.fasta result_assemblies/MBOVb_CLC_MIRA_mix.fasta \
	result_assemblies/MBOVb_CLC_naive_concatenation.fasta
	rm -rf result_statistics/MBOVb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MBVG_AB_CLC.fasta: $(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_CLC.fasta
temp_assemblies/MBVG_AB_MIRA.fasta: $(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta
temp_assemblies/MBVG_CLC_MIRA.fasta: $(MOLLI)/MBVG/MBVG_CLC.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta

result_statistics/MBVG_quast: \
	$(MOLLI)/MBVG/MBVG_ABySS_35.fasta $(MOLLI)/MBVG/MBVG_CLC.fasta $(MOLLI)/MBVG/MBVG_MIRA.fasta \
	$(MOLLIGAM)/MBVG/MBVG_ABySS-CLC_gam.fasta $(MOLLIGAM)/MBVG/MBVG_CLC-mira_gam.fasta $(MOLLIGAM)/MBVG/MBVG_mira-ABySS_gam.fasta \
	result_assemblies/MBVG_AB_CLC_mix.fasta result_assemblies/MBVG_AB_MIRA_mix.fasta result_assemblies/MBVG_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MBVG_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MCCP_AB_CLC.fasta: $(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_CLC.fasta
temp_assemblies/MCCP_AB_MIRA.fasta: $(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta
temp_assemblies/MCCP_CLC_MIRA.fasta: $(MOLLI)/MCCP/MCCP_CLC.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta

result_statistics/MCCP_quast: \
	$(MOLLI)/MCCP/MCCP_ABySS_35.fasta $(MOLLI)/MCCP/MCCP_CLC.fasta $(MOLLI)/MCCP/MCCP_MIRA.fasta \
	$(MOLLIGAM)/MCCP/MCCP_ABySS-CLC_gam.fasta $(MOLLIGAM)/MCCP/MCCP_CLC-mira_gam.fasta $(MOLLIGAM)/MCCP/MCCP_mira-ABySS_gam.fasta \
	result_assemblies/MCCP_AB_CLC_mix.fasta result_assemblies/MCCP_AB_MIRA_mix.fasta result_assemblies/MCCP_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MCCP_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MOVI_AB_CLC.fasta: $(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_CLC.fasta
temp_assemblies/MOVI_AB_MIRA.fasta: $(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta
temp_assemblies/MOVI_CLC_MIRA.fasta: $(MOLLI)/MOVI/MOVI_CLC.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta

result_statistics/MOVI_quast: \
	$(MOLLI)/MOVI/MOVI_ABySS_32.fasta $(MOLLI)/MOVI/MOVI_CLC.fasta $(MOLLI)/MOVI/MOVI_MIRA.fasta \
	$(MOLLIGAM)/MOVI/MOVI_ABySS-CLC_gam.fasta $(MOLLIGAM)/MOVI/MOVI_CLC-mira_gam.fasta $(MOLLIGAM)/MOVI/MOVI_mira-ABySS_gam.fasta \
	result_assemblies/MOVI_AB_CLC_mix.fasta result_assemblies/MOVI_AB_MIRA_mix.fasta result_assemblies/MOVI_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MOVI_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MMC_AB_CLC.fasta: $(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_CLC.fasta
temp_assemblies/MMC_AB_MIRA.fasta: $(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_MIRA.fasta
temp_assemblies/MMC_CLC_MIRA.fasta: $(MOLLI)/MMC/MMC_CLC.fasta $(MOLLI)/MMC/MMC_MIRA.fasta

result_statistics/MMC_quast: \
	$(MOLLI)/MMC/MMC_ABySS_30.fasta $(MOLLI)/MMC/MMC_CLC.fasta $(MOLLI)/MMC/MMC_MIRA.fasta \
	$(MOLLIGAM)/MMC/MMC_ABySS-CLC_gam.fasta $(MOLLIGAM)/MMC/MMC_CLC-mira_gam.fasta $(MOLLIGAM)/MMC/MMC_mira-ABySS_gam.fasta \
	result_assemblies/MMC_AB_CLC_mix.fasta result_assemblies/MMC_AB_MIRA_mix.fasta result_assemblies/MMC_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MMC_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCe_AB_CLC.fasta: $(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_CLC.fasta
temp_assemblies/MSCe_AB_MIRA.fasta: $(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta
temp_assemblies/MSCe_CLC_MIRA.fasta: $(MOLLI)/MSCe/MSCe_CLC.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta

result_statistics/MSCe_quast: \
	$(MOLLI)/MSCe/MSCe_ABySS_35.fasta $(MOLLI)/MSCe/MSCe_CLC.fasta $(MOLLI)/MSCe/MSCe_MIRA.fasta \
	$(MOLLIGAM)/MSCe/MSCe_ABySS-CLC_gam.fasta $(MOLLIGAM)/MSCe/MSCe_CLC-mira_gam.fasta $(MOLLIGAM)/MSCe/MSCe_mira-ABySS_gam.fasta \
	result_assemblies/MSCe_AB_CLC_mix.fasta result_assemblies/MSCe_AB_MIRA_mix.fasta result_assemblies/MSCe_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCe_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCd_AB_CLC.fasta: $(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_CLC.fasta
temp_assemblies/MSCd_AB_MIRA.fasta: $(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta
temp_assemblies/MSCd_CLC_MIRA.fasta: $(MOLLI)/MSCd/MSCd_CLC.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta

result_statistics/MSCd_quast: \
	$(MOLLI)/MSCd/MSCd_ABySS_34.fasta $(MOLLI)/MSCd/MSCd_CLC.fasta $(MOLLI)/MSCd/MSCd_MIRA.fasta \
	$(MOLLIGAM)/MSCd/MSCd_ABySS-CLC_gam.fasta $(MOLLIGAM)/MSCd/MSCd_CLC-mira_gam.fasta $(MOLLIGAM)/MSCd/MSCd_mira-ABySS_gam.fasta \
	result_assemblies/MSCd_AB_CLC_mix.fasta result_assemblies/MSCd_AB_MIRA_mix.fasta result_assemblies/MSCd_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCd_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCc_AB_CLC.fasta: $(MOLLI)/MSCc/MSCc_ABySS_27.fasta $(MOLLI)/MSCc/MSCc_CLC.fasta
temp_assemblies/MSCc_AB_MIRA.fasta: $(MOLLI)/MSCc/MSCc_ABySS_27.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta
temp_assemblies/MSCc_CLC_MIRA.fasta: $(MOLLI)/MSCc/MSCc_CLC.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta

result_statistics/MSCc_quast: \
	$(MOLLI)/MSCc/MSCc_ABySS_27.fasta $(MOLLI)/MSCc/MSCc_CLC.fasta $(MOLLI)/MSCc/MSCc_MIRA.fasta \
	$(MOLLIGAM)/MSCc/MSCc_ABySS-CLC_gam.fasta $(MOLLIGAM)/MSCc/MSCc_CLC-mira_gam.fasta $(MOLLIGAM)/MSCc/MSCc_mira-ABySS_gam.fasta \
	result_assemblies/MSCc_AB_CLC_mix.fasta result_assemblies/MSCc_AB_MIRA_mix.fasta result_assemblies/MSCc_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCc_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^

temp_assemblies/MSCb_AB_CLC.fasta: $(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_CLC.fasta
temp_assemblies/MSCb_AB_MIRA.fasta: $(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta
temp_assemblies/MSCb_CLC_MIRA.fasta: $(MOLLI)/MSCb/MSCb_CLC.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta

result_statistics/MSCb_quast: \
	$(MOLLI)/MSCb/MSCb_ABySS_29.fasta $(MOLLI)/MSCb/MSCb_CLC.fasta $(MOLLI)/MSCb/MSCb_MIRA.fasta \
	$(MOLLIGAM)/MSCb/MSCb_ABySS-CLC_gam.fasta $(MOLLIGAM)/MSCb/MSCb_CLC-mira_gam.fasta $(MOLLIGAM)/MSCb/MSCb_mira-ABySS_gam.fasta \
	result_assemblies/MSCb_AB_CLC_mix.fasta result_assemblies/MSCb_AB_MIRA_mix.fasta result_assemblies/MSCb_CLC_MIRA_mix.fasta
	rm -rf result_statistics/MSCb_quast
	$(pyinterp) $(QUASTMOLLI) -o $@ $^