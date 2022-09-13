# This script removes duplicate (100% overlap) sequences in your fasta files
import os
from Bio import SeqIO
import argparse
import pathlib
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=pathlib.Path, default='PICI_results')
parser.add_argument("--output", type=pathlib.Path, default='PICI_results')
args = parser.parse_args()

seen = []
records = []

for record in SeqIO.parse(args.input, "fasta"):  
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)

SeqIO.write(records, args.output, "fasta")
