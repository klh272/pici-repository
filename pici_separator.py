import collections
import os
import pandas as pd
#!pip3 install Biopython
from Bio import SeqIO

handle = open("ALL_PICIs.fasta", "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"), key_function = lambda rec : rec.description.replace(" ", "__"))
handle.close()

# Identify PICI type and separate into dictionary
g_neg = dict([(key,val) for key,val in record_dict.items() if "G_neg_PICI" in key])
sapi = dict([(key,val) for key,val in record_dict.items() if "SaPI" in key])
g_pos = dict([(key,val) for key,val in record_dict.items() if "phage_satellite" in key])

# Only write to file if type exists
# G-
if bool(g_neg) == True:
  print("Writing G- PICI file...")
  with open('G_neg_PICIs.fasta', 'w') as handle:
      SeqIO.write(g_neg.values(), handle, 'fasta')
else:
  print("No G- PICIs found. No file will be outputted for G- PICIs.")

# SaPI
if bool(sapi) == True:
  print("Writing SaPIs file...")
  with open('SaPIs.fasta', 'w') as handle:
      SeqIO.write(sapi.values(), handle, 'fasta')
else:
  print("No SaPIs found. No file will be outputted for SaPIs.")

# Phage Satellite
if bool(g_pos) == True:
  print("Writing phage satellites file...")
  with open('Phage_Satellites.fasta', 'w') as handle:
      SeqIO.write(g_pos.values(), handle, 'fasta')
else:
  print("No phage satellites found. No file will be outputted for phage satellites.")

