wget -N ccb.jhu.edu/gage_b/genomeAssemblies/A_hydrophila_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/B_cereus_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/B_cereus_MiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/B_fragilis_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/M_abscessus_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/M_abscessus_MiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/R_sphaeroides_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/R_sphaeroides_MiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/S_aureus_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/V_cholerae_HiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/V_cholerae_MiSeq.tar.gz \
ccb.jhu.edu/gage_b/genomeAssemblies/X_axonopodis_HiSeq.tar.gz 
for file in `ls *.tar.gz`
do
TARGET="`basename $file .tar.gz`"
mkdir $TARGET 
tar zxvf $file -C $TARGET

done
rm *.tar.gz
