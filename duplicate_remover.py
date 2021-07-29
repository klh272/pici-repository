# This script removes duplicate (100% overlap) sequences in your fasta files
import os
from Bio import SeqIO

seen = []
records = []

for record in SeqIO.parse("PICI_results", "fasta"):  
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)

SeqIO.write(records, "PICI_results", "fasta")
