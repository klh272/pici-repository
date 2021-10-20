import re
import pandas as pd
from Bio import SeqIO

# read in BLAST alignment output
df = pd.read_csv('BLAST_results.out', sep='\t', header = None)

#for i in range(len(df)):
#  df.iloc[i,1] = re.split('(\d+)', df.iloc[i,1])[0] 



#


# take the max identity % (column "2") for each alignment and drop the duplicates
condensed_df = df.sort_values([2], ascending=False).sort_index()
condensed_df.reset_index(drop=True, inplace=True)


# lists to hold PICI headers and sequences for fasta output
neg_name_list = []
neg_seq_list = []
neg_desc_list = []
satellite_name_list = []
satellite_seq_list = []
satellite_desc_list = []

fasta_file = 'Phage_Satellites.fasta'
with open(fasta_file, mode='r') as handle:
  
  for record in SeqIO.parse(handle, 'fasta'):

    for i in range(len(condensed_df)):
      if condensed_df.iloc[i,0].startswith(record.id): # ensures same sequence in multi-sequence fasta file is being used
        if condensed_df.iloc[i,1].startswith('alpA'): # find alpA'ss in df

            # add PICI to list
            neg_name_list.append(record.id)
            neg_seq_list.append(record.seq)
            neg_desc_list.append(record.description)
            break
        else:
          continue

fasta_file = 'Phage_Satellites.fasta'
with open(fasta_file, mode='r') as handle:

	for record in SeqIO.parse(handle, 'fasta'):
		if str(record.seq) not in neg_seq_list:
			satellite_name_list.append(record.id)
			satellite_seq_list.append(record.seq)
			satellite_desc_list.append(record.description)

print(f"G- PICIs detected and changed: {len(neg_seq_list)}")
print(f"Phage satellites confirmed: {len(satellite_seq_list)}")

# write PICI sequence to file
neg_PICI_file = open("G_neg_PICI_reviewed", "w")

for m in range(len(neg_seq_list)):
  neg_PICI_file.write(">"  + str(neg_desc_list[m]) + "\n" + str(neg_seq_list[m]) + "\n")

neg_PICI_file.close()


phage_satellite_file = open("phage_satellite_reviewed", "w")

for n in range(len(satellite_seq_list)):
  phage_satellite_file.write(">"  + str(satellite_desc_list[n]) + "\n" + str(satellite_seq_list[n]) + "\n")

phage_satellite_file.close()

