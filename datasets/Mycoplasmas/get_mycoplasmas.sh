wget -N \
http://services.cbib.u-bordeaux2.fr/mix/MAUR/single-assemblies/MAUR_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MAUR/single-assemblies/MAUR_ABySS_27.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MAUR/single-assemblies/MAUR_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBOVb/single-assemblies/MBOVb_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBOVb/single-assemblies/MBOVb_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBOVb/single-assemblies/MBOVb_ABySS_34.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBVG/single-assemblies/MBVG_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBVG/single-assemblies/MBVG_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MBVG/single-assemblies/MBVG_ABySS_35.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MCCP/single-assemblies/MCCP_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MCCP/single-assemblies/MCCP_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MCCP/single-assemblies/MCCP_ABySS_35.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MMC/single-assemblies/MMC_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MMC/single-assemblies/MMC_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MMC/single-assemblies/MMC_ABySS_30.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MOVI/single-assemblies/MOVI_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MOVI/single-assemblies/MOVI_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MOVI/single-assemblies/MOVI_ABySS_32.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCb/single-assemblies/MSCb_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCb/single-assemblies/MSCb_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCb/single-assemblies/MSCb_ABySS_29.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCc/single-assemblies/MSCc_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCc/single-assemblies/MSCc_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCc/single-assemblies/MSCc_ABySS_27.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCd/single-assemblies/MSCd_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCd/single-assemblies/MSCd_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCd/single-assemblies/MSCd_ABySS_34.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCe/single-assemblies/MSCe_CLC.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCe/single-assemblies/MSCe_MIRA.fasta \
http://services.cbib.u-bordeaux2.fr/mix/MSCe/single-assemblies/MSCe_ABySS_35.fasta 

mkdir MAUR
mv ./MAUR*.fasta ./MAUR/
mkdir MBOVb
mv ./MBOVb*.fasta ./MBOVb/
mkdir MBVG
mv ./MBVG*.fasta ./MBVG/
mkdir MCCP
mv ./MCCP*.fasta ./MCCP/
mkdir MMC
mv ./MMC*.fasta ./MMC/
mkdir MOVI
mv ./MOVI*.fasta ./MOVI/
mkdir MSCb
mv ./MSCb*.fasta ./MSCb/
mkdir MSCc
mv ./MSCc*.fasta ./MSCc/
mkdir MSCd
mv ./MSCd*.fasta ./MSCd/
mkdir MSCe
mv ./MSCe*.fasta ./MSCe/