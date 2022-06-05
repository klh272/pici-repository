import os
import re
import pandas as pd
import math

from Bio.Seq import Seq
from Bio import SeqIO
#from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

########################################################

# use both fwd and rev sequences
#def translate_seq(record):
#        #dna_seqs = [record.seq, record.seq.reverse_complement()]
#        print(record)
#        # generate all translation frames
#        aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in record)
#        print(aa_seqs)
#
#        # select the longest one
#        max_aa = max(aa_seqs, key=len)
#
#        # write new record
#        #aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
#        #SeqIO.write(aa_record, aa_fa, 'fasta')
#
#        return (max_aa)

########################################################




# Construct DB, add on if it exists
# make counter for keeping track of unique ID number
filesize = os.path.getsize("BLAST_DB.csv")
df = pd.read_csv("BLAST_DB.csv", sep=',', header = None)


if filesize == 0:
	int_id_num = 0
	alpA_id_num = 0
	sis_id_num = 0
	prirep_id_num = 0
	terS_id_num = 0
	ppi_id_num = 0
	rpp_id_num = 0
else:

	# Create list of previous DB
	prev_id = df[0].tolist()
	prev_seq = df[1].tolist()

	integrases = df[df[0].str.contains('int')]
	int_id_num = int(integrases[0].iloc[-1].split("int")[1])
	del integrases

	alpAs = df[df[0].str.contains('alpA')]
	alpA_id_num = int(alpAs[0].iloc[-1].split("alpA")[1])
	del alpAs

	siss = df[df[0].str.contains('sis')]
	sis_id_num = int(siss[0].iloc[-1].split("sis")[1])
	del siss

	prireps = df[df[0].str.contains('pri-rep')]
	prirep_id_num = int(prireps[0].iloc[-1].split("pri-rep")[1])
	del prireps

	terSs = df[df[0].str.contains('terS')]
	terS_id_num = int(terSs[0].iloc[-1].split("terS")[1])
	del terSs

	ppis = df[df[0].str.contains('ppi')]
	ppi_id_num = int(ppis[0].iloc[-1].split("ppi")[1])
	del ppis

	rpps = df[df[0].str.contains('rpp')]
	rpp_id_num = int(rpps[0].iloc[-1].split("rpp")[1])
	del rpps



# make master database lists
master_id = []
master_seq = []

# read in table of PICI/phage satellite results
tbl = pd.read_csv('RefSeq_PICI_table.tsv', sep='\t', header = None)

# read in host genomes
fasta_file = 'RefSeq_ALL_PICIs_host_genomes.fasta'
with open(fasta_file, mode='r') as handle:

	for record in SeqIO.parse(handle, 'fasta'):
		for i in range(len(tbl)):
			if tbl.iloc[i,0].startswith(record.id):

				print(record.description)

				# integrases
				int_id_num += 1
				master_id.append("int" + str(int_id_num))
				master_seq.append(record.seq[int(tbl.iloc[i,8]):int(tbl.iloc[i,9])])
				#print(translate_seq((record.seq[tbl.iloc[i,8]:tbl.iloc[i,9]], record.seq[tbl.iloc[i,8]:tbl.iloc[i,9]].reverse_complement())) + "*")
				#print(record.seq[tbl.iloc[i,8]:tbl.iloc[i,9]].translate() + "*")

				# alpAs
				if math.isnan(tbl.iloc[i,12]) == False:
					if  record.seq[int(tbl.iloc[i,12]):int(tbl.iloc[i,13])] != "":
						alpA_id_num += 1
						master_id.append("alpA" + str(alpA_id_num))
						master_seq.append(record.seq[int(tbl.iloc[i,12]):int(tbl.iloc[i,13])])

				# siss
				if math.isnan(tbl.iloc[i,16]) == False:
					if record.seq[int(tbl.iloc[i,16]):int(tbl.iloc[i,17])] != "":
						sis_id_num += 1
						master_id.append("sis" + str(sis_id_num))
						master_seq.append(record.seq[int(tbl.iloc[i,16]):int(tbl.iloc[i,17])])

				# pri-reps
				prirep_id_num += 1
				master_id.append("pri-rep" + str(prirep_id_num))
				master_seq.append(record.seq[int(tbl.iloc[i,20]):int(tbl.iloc[i,21])])

				#terS
				


# write new DB
NEW_DB_FILE = open("NEW_BLAST_DB", "w")

for l in range(len(prev_seq)):
	NEW_DB_FILE.write(">" + str(prev_id[l]) + "\n" + str(prev_seq[l]) + "\n")

for m in range(len(master_seq)):
	NEW_DB_FILE.write(">" + str(master_id[m]) + "\n" + str(master_seq[m]) + "\n")
	#print(str(master_id[m]) + "\t" + master_seq[m])

NEW_DB_FILE.close()

#				if tbl.iloc[i,8] != "NA":
#print("MASTER ID:", master_id)
#print("MASTER SEQ:", master_seq)



##########################################################################################################################################################################
