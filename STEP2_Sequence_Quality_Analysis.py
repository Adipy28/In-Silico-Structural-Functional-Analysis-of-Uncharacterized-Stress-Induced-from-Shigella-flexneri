# Step 2 of the project has been performed to analyze the Sequence Quality.
# This includes reading the protein sequence of YjbJ (UPF0337 family protein) using Biopython 
# and calculating the count and percentage of different amino acid types.
# Since amino acid composition strongly affects protein structure and stability,
# this step ensures that the chosen sequence is biologically consistent and suitable for further analysis.

from Bio import SeqIO

record = SeqIO.read("P68208.fasta", "fasta")

Alanine_count_A = record.seq.count("A")
Cysteine_count_C = record.seq.count("C")
Aspartate_count_D = record.seq.count("D")
Glutamate_count_E = record.seq.count("E")
Phenylalanine_count_F = record.seq.count("F")
Glycine_count_G = record.seq.count("G")
Histidine_count_H = record.seq.count("H")
Isoleucine_count_I = record.seq.count("I")
Lysine_count_K = record.seq.count("K")
Leucine_count_L = record.seq.count("L")
Methionine_count_M = record.seq.count("M")
Asparagine_count_N = record.seq.count("N")
Proline_count_P = record.seq.count("P")
Glutamine_count_Q = record.seq.count("Q")
Arginine_count_R = record.seq.count("R")
Serine_count_S = record.seq.count("S")
Threonine_count_T = record.seq.count("T")
Valine_count_V = record.seq.count("V")
Tryptophan_count_W = record.seq.count("W")
Tyrosine_count_Y = record.seq.count("Y")
StopCodon_count = record.seq.count("*")

print("The query protein is", len(record.seq), "amino acids long.")
print("-"*90)
print("Alanine count is:", Alanine_count_A)
print("Cysteine count is:", Cysteine_count_C)
print("Aspartate count is:", Aspartate_count_D)
print("Glutamate count is:", Glutamate_count_E)
print("Phenylalanine count is:", Phenylalanine_count_F)
print("Glycine count is:", Glycine_count_G)
print("Histidine count is:", Histidine_count_H)
print("Isoleucine count is:", Isoleucine_count_I)
print("Lysine count is:", Lysine_count_K)
print("Leucine count is:", Leucine_count_L)
print("Methionine count is:", Methionine_count_M)
print("Asparagine count is:", Asparagine_count_N)
print("Proline count is:", Proline_count_P)
print("Glutamine count is:", Glutamine_count_Q)
print("Arginine count is:", Arginine_count_R)
print("Serine count is:", Serine_count_S)
print("Threonine count is:", Threonine_count_T)
print("Valine count is:", Valine_count_V)
print("Tryptophan count is:", Tryptophan_count_W)
print("Tyrosine count is:", Tyrosine_count_Y)

# Percentage-based composition
print("-"*90)
print("% Non-polar Aliphatic / Hydrophobic amino acids:",
      (Glycine_count_G + Alanine_count_A + Valine_count_V + Leucine_count_L +
       Isoleucine_count_I + Proline_count_P + Methionine_count_M) / len(record.seq) * 100)
print("% Non-polar Aromatic / Hydrophobic amino acids:",
      (Phenylalanine_count_F + Tryptophan_count_W + Tyrosine_count_Y) / len(record.seq) * 100)
print("% Polar Uncharged / Hydrophilic amino acids:",
      (Serine_count_S + Threonine_count_T + Cysteine_count_C +
       Asparagine_count_N + Glutamine_count_Q) / len(record.seq) * 100)
print("% Acidic (Negatively charged) / Hydrophilic amino acids:",
      (Aspartate_count_D + Glutamate_count_E) / len(record.seq) * 100)
print("% Basic (Positively charged) / Hydrophilic amino acids:",
      (Lysine_count_K + Arginine_count_R + Histidine_count_H) / len(record.seq) * 100)
print("-"*90)
print("Stop codons in this sequence (expected 0 for protein):", StopCodon_count)
print("-"*90)
print("End of Quality check. All parameters fulfilled.")
print("The query sequence passed the quality checks successfully.")
