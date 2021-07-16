# -*- coding: utf-8 -*-
"""pici_integrase_trimmer_script.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1qYyKzGVdDNsEj-h6wrKEVhVOeWSJCdNB
"""

import os
import pandas as pd
from Bio import SeqIO

# read in tBLASTn alignment output
df = pd.read_csv('tBLASTn_results.out', sep='\t', header = None)

# lists to hold headers and sequences for fasta output
name_list = []
seq_list = []
description_list = []

# read in sequences and run integrase search/trim
fasta_file = 'all.fna' 
with open(fasta_file, mode='r') as handle:
  
  for record in SeqIO.parse(handle, 'fasta'):
    for i in range(len(df)):
      if df.iloc[i,1].startswith(record.id): # ensures same sequence in multi-sequence fasta file is being used
        if df.iloc[i,0].startswith('int'): # find integrases in df
          if df.iloc[i,2] >= 90: # check it integrase identity is above 90%
            print('\nIntegrase found with identity >= 90% in', record.description)
            print('Location: {} {}'.format(df.iloc[i,8], df.iloc[i,9]))
            print('BLAST database ID: ', df.iloc[i,0])
            print('Beginning trim...')

            # get integrase start/end
            int_start = df.iloc[i,8]
            int_end = df.iloc[i,9]

            trim_low = int_start - 30000
            trim_high = int_end + 30000
            
            # fix lower limit index to 0 if it is negative
            if trim_low < 0:
              trim_low = 0
            # fix upper limit index to the length of the sequence if it is over
            if trim_high >= len(record.seq):
              trim_high = len(record.seq)
            
            # Append trim to list
            name_list.append(record.id)
            seq_list.append(record.seq[trim_low:trim_high])
            description_list.append(record.description)
            print('Trim Finished.\n')

        else:
          print('No hit')


# write trimmed sequence(s) to file
print('Writing {} trims to output file \"trimmed_file\"'.format(len(seq_list)))
trimmed_file = open("trimmed_file", "w")

for m in range(len(seq_list)):
  trimmed_file.write(">" + str(name_list[m]) + " " +str(description_list[m]) + ", integrase_" + str(m) + "\n" + str(seq_list[m]) + "\n")

trimmed_file.close()