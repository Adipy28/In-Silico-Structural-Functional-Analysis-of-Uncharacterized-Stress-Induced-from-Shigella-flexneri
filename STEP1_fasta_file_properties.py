#Step 1 (Biological Sequence Selection/Retrieval), I chose a stress-induced uncharacterized bacterial protein.
#A gene from Shigella flexneri, a pathogen causing human bacillary dysentery, is selected.
#This gene encodes a small 7.7 kDa protein known as YjbJ, belonging to the UPF0337 family.
#The biological function of YjbJ in S. flexneri is still not clearly understood.
#Abbreviations: UPF – Uncharacterized Protein Family; YjbJ – Stress-Induced Protein.
print("Step 1: performed by Aditya Panigrahi")
print("-"*90)
from Bio import SeqIO
record = SeqIO.read("P68208.fasta", "fasta")
print("The following protein record has Nucleotide ID:", record.id)
print(record.description)
print("The length of the sequence is:", len(record.seq), "amino acids")
print("The protein sequence is:", record.seq)
