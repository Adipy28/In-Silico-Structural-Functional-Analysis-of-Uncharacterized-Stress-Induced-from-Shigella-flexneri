# Step 4 of the project performs Homology Analysis using NCBI BLASTp.
# The amino acid sequence of the UPF0337 family protein YjbJ (P68208) from Shigella flexneri
# is compared against the NCBI RefSeq protein database to identify homologous proteins.
# This helps in predicting possible functions and evolutionary relationships.

from Bio.Blast import NCBIWWW
from Bio import SeqIO

record = SeqIO.read("P68208.fasta", "fasta")
result_handle = NCBIWWW.qblast(
    program= "blastp",
    database= "nr",
    sequence= record.seq
)


with open("result_P68208_blast.xml", "w") as a:
    a.write(result_handle.read())

print("BLAST completed successfully!")

