wget -N http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz
mkdir Staphylococcus_aureus
tar zxvf Assembly.tgz -C Staphylococcus_aureus
mv Staphylococcus_aureus/Assembly/* Staphylococcus_aureus
rm -r Staphylococcus_aureus/Assembly
cd Staphylococcus_aureus
mv ABySS2/genome.ctg.fasta ABySS2/abyss2.genome.ctg.fasta
mv ABySS/genome.ctg.fasta ABySS/abyss.genome.ctg.fasta
mv Allpaths-LG/genome.ctg.fasta Allpaths-LG/lg.genome.ctg.fasta
mv Bambus2/genome.ctg.fasta Bambus2/bambus2.genome.ctg.fasta
mv MSR-CA/genome.ctg.fasta MSR-CA/MSR.genome.ctg.fasta
mv SGA/genome.ctg.fasta SGA/SGA.genome.ctg.fasta
mv SOAPdenovo/genome.ctg.fasta SOAPdenovo/SOAPdenovo.genome.ctg.fasta
mv Velvet/genome.ctg.fasta Velvet/velvet.genome.ctg.fasta
cd ..
rm Assembly.tgz


wget -N http://gage.cbcb.umd.edu/data/Rhodobacter_sphaeroides/Assembly.tgz
mkdir Rhodobacter_sphaeroides
tar zxvf Assembly.tgz -C Rhodobacter_sphaeroides
mv Rhodobacter_sphaeroides/Assembly/* Rhodobacter_sphaeroides
rm -r Rhodobacter_sphaeroides/Assembly
cd Rhodobacter_sphaeroides
mv ABySS2/genome.ctg.fasta ABySS2/abyss2.genome.ctg.fasta
mv ABySS/genome.ctg.fasta ABySS/abyss.genome.ctg.fasta
mv Allpaths-LG/genome.ctg.fasta Allpaths-LG/lg.genome.ctg.fasta
mv Bambus2/genome.ctg.fasta Bambus2/bambus2.genome.ctg.fasta
mv MSR-CA/genome.ctg.fasta MSR-CA/MSR.genome.ctg.fasta
mv SGA/genome.ctg.fasta SGA/SGA.genome.ctg.fasta
mv SOAPdenovo/genome.ctg.fasta SOAPdenovo/SOAPdenovo.genome.ctg.fasta
mv Velvet/genome.ctg.fasta Velvet/velvet.genome.ctg.fasta
cd ..
