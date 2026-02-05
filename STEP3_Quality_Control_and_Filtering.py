# Step 3 of the project performs Sequence Retrieval and Quality Control.
# The protein sequence of YjbJ (UPF0337 family protein) from Shigella flexneri (P68208)
# was retrieved from UniProt and verified for sequence quality.
# This script filters and validates the amino acid sequence to ensure it is
# biologically consistent and ready for downstream analysis.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Allowed amino acids (standard 20 + stop codon for QC check)
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY*")

# Input and output files
input_fasta = "P68208.fasta"
output_fasta = "filtered_YjbJ_P68208.fasta"

filtered_records = []

# Iterate through all sequences in FASTA (usually 1 for UniProt)
for record in SeqIO.parse(input_fasta, "fasta"):

    protein_seq = str(record.seq).upper()

    # ---- Basic QC ----
    if len(protein_seq) < 30:
        print(f"Sequence too short: {record.id} ({len(protein_seq)} aa)")
        continue

    # Check only valid amino acids
    if not set(protein_seq).issubset(VALID_AA):
        print(f"Invalid amino acid symbols found in {record.id}")
        continue

    # Check for internal stop codons
    if "*" in protein_seq[:-1]:
        print(f"Internal stop codon found in {record.id}")
        continue

    # Check if sequence ends correctly (no trailing '*')
    if protein_seq.endswith("*"):
        protein_seq = protein_seq.rstrip("*")  # remove final stop for clean FASTA

    # If all checks pass, save it
    protein_record = SeqRecord(
        record.seq,
        id=record.id,
        description="YjbJ UPF0337 family protein (QC passed)"
    )
    filtered_records.append(protein_record)

# Write filtered results
SeqIO.write(filtered_records, output_fasta, "fasta")

print(f"Saved {len(filtered_records)} high-quality protein sequence(s) to {output_fasta}.")
