import os
import pandas as pd
from Bio import SeqIO

# read in tBLASTn alignment output
df = pd.read_csv('tBLASTn_results.out', sep='\t', header = None)
# create dictionary for fasta file
record_dict = SeqIO.to_dict(SeqIO.parse("all.fna", "fasta"))

def range_subset(range1, range2):
    """Whether range1 is a subset of range2. 
    https://stackoverflow.com/questions/32480423/how-to-check-if-a-range-is-a-part-of-another-range-in-python-3-x
    """
    if not range1:
        return True  # empty range is subset of anything
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step
    return range1.start in range2 and range1[-1] in range2

# lists to hold headers and sequences for fasta output
name_list = []
seq_list = []
description_list = []
record_locations = []

# search BLAST results for integrases
for i in range(len(df)):
    if df.iloc[i,0].startswith('int'): # find integrases in df
      if df.iloc[i,2] >= 90: # check it integrase identity is above 90%
        print('\nIntegrase found with identity >= 90% in', record_dict[df.iloc[i,1]].description)
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
        if trim_high >= len(record_dict[df.iloc[i,1]].seq):
          trim_high = len(record.seq)
            
        # Append trim to list
        if record_dict[df.iloc[i,1]].id in name_list:
          idx = name_list.index(record_dict[df.iloc[i,1]].id)
          print(idx)
          if range_subset(range(trim_low,trim_high), record_locations[idx]) == True:
            continue
        else:
          name_list.append(record_dict[df.iloc[i,1]].id)
          seq_list.append(record_dict[df.iloc[i,1]].seq[trim_low:trim_high])
          description_list.append(record_dict[df.iloc[i,1]].description)
          record_locations.append(range(trim_low,trim_high))
        print('Trim Finished.\n')

    else:
      print('No hit')


# write trimmed sequence(s) to file
print('Writing {} trims to output file \"trimmed_file\"'.format(len(seq_list)))
trimmed_file = open("trimmed_file", "w")

for m in range(len(seq_list)):
  trimmed_file.write(">" + str(name_list[m]) + " " +str(description_list[m]) + ", integrase_" + str(m) + "\n" + str(seq_list[m]) + "\n")

trimmed_file.close()
