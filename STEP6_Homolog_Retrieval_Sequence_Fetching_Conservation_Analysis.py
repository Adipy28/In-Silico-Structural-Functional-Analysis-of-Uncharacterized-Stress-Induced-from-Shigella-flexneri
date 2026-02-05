# Step 5 of the project: Homolog Retrieval, Sequence Fetching, and Conservation Analysis
# Using BLAST results of the YjbJ (UPF0337 family) protein from Shigella flexneri (UniProt: P68208)
# This step identifies closest homologs, fetches their sequences from NCBI,
# and analyzes conserved amino acid regions across top hits.

from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
import ssl
import os

# Fix SSL issues for NCBI access
ssl._create_default_https_context = ssl._create_unverified_context

# Step 0: Set your NCBI email ID (mandatory)
Entrez.email = "aditya.dsp05@gmail.com"  

# Step 1: Parse the BLAST results (from Step 3)
blast_xml = "blast_result.xml"

if not os.path.exists(blast_xml):
    print(f"❌ ERROR: File {blast_xml} not found.")
    print("Please ensure you have run STEP3_BlastP_Analysis.py first.")
    exit()

print("🔍 Parsing BLAST XML results for YjbJ (P68208)...\n")

with open(blast_xml) as handle:
    blast_records = list(NCBIXML.parse(handle))

blast_record = blast_records[0]  # assuming one query sequence

# Step 2: Identify closest homologs based on E-value and identity
closest_hits = []

for alignment in blast_record.alignments:
    best_hsp = max(alignment.hsps, key=lambda h: h.bits)
    identity_fraction = best_hsp.identities / best_hsp.align_length

    # Filter significant hits (E-value < 1e-5 and identity ≥ 35%)
    if best_hsp.expect < 1e-5 and identity_fraction >= 0.35:
        closest_hits.append({
            "title": alignment.title,
            "accession": alignment.accession,
            "identity": identity_fraction,
            "evalue": best_hsp.expect,
            "q_start": best_hsp.query_start,
            "q_end": best_hsp.query_end
        })

print(f"✅ Number of closest homologs found: {len(closest_hits)}\n")

print("Top 10 homologs (similarity information):")
for hit in closest_hits[:10]:
    print(f"- {hit['title']}")
    print(f"  Identity: {hit['identity']:.2%}")
    print(f"  E-value: {hit['evalue']}")
    print(f"  Query region: {hit['q_start']}-{hit['q_end']}\n")

# Step 3: Fetch sequences of top 10 homologs from NCBI
records = []
print("="*90)
print("⬇️ Fetching sequences of top homologs from NCBI (please wait)...")

for hit in closest_hits[:10]:
    try:
        handle = Entrez.efetch(
            db="protein",
            id=hit["accession"],
            rettype="fasta",
            retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        records.append(record)
        print(f"  ✓ Retrieved: {record.id}")
    except Exception as e:
        print(f"  ⚠️ Failed to fetch {hit['accession']}: {e}")

# Save sequences to FASTA
if len(records) > 0:
    SeqIO.write(records, "homologs_P68208.fasta", "fasta")
    print("\n✅ Sequences saved successfully as 'homologs_P68208.fasta'")
else:
    print("\n❌ No sequences retrieved. Check your internet connection or NCBI limits.")
    exit()

print("="*90)

# Step 4: Detect conserved regions (simple frequency-based conservation)
if len(records) > 1:
    seq_length = min(len(r.seq) for r in records)  # trim to shortest sequence
    conserved_positions = []

    for i in range(seq_length):
        column = [str(r.seq[i]) for r in records]
        most_common = Counter(column).most_common(1)[0]  # (aa, count)
        conservation = most_common[1] / len(records)
        if conservation >= 0.8:  # ≥80% conserved
            conserved_positions.append((i + 1, most_common[0]))

    print("\n⭐ Highly conserved positions (first 20 shown):")
    print(conserved_positions[:20], "...")
else:
    print("\n⚠️ Not enough homologous sequences for conserved region analysis.")

# Step 5: Evolutionary hints (organism information)
print("\n🌍 Evolutionary hints (organism info from top hits):")
for hit in closest_hits[:10]:
    print(f"- {hit['title']}")

print("\n🎯 Conservation and evolutionary analysis for YjbJ completed successfully.")
print("="*90)
