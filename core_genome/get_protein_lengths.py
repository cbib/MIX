from Bio import SeqUtils
from Bio import SeqIO

from Bio.Alphabet import generic_protein

afile="mollicutes_core_genome_cluster_of_proteins.fasta"

for record in SeqIO.parse(afile, "fasta", generic_protein):
	print record.id+"\t"+str(len(record.seq))
