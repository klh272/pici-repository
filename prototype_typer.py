import os
import re
import pandas as pd
import pathlib
import argparse

from Bio.Seq import Seq
from Bio import SeqIO
#from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

##########################################################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--temp', type=pathlib.Path, required=True)
parser.add_argument('--aa', type=pathlib.Path, default='all.pdg.faa')
parser.add_argument('--blast', type=pathlib.Path, default='BLASTp_results.out')
parser.add_argument('--fasta', type=pathlib.Path, default='all.fna')
parser.add_argument('--output', type=pathlib.Path, default='trimmed_file')
parser.add_argument('--i', nargs='?', const=1, type=int, default=70)
parser.add_argument('--a', nargs='?', const=1, type=int, default=50)

args = parser.parse_args()
args.temp = str(args.temp)
print('Int Identity percentage:', args.i)
print('AlpA Identity percentage:', args.a)


##########################################################################################################################################################################

# read in BLAST alignment output
df = pd.read_csv(args.blast, sep='\t', header = None)
# take the max identity % (column "2") for each alignment and drop the duplicates
condensed_df = df.sort_values([2], ascending=False).drop_duplicates([0]).sort_index()
condensed_df.reset_index(drop=True, inplace=True)

##########################################################################################################################################################################

def get_locations(file, key_word):
  """Returns the locations of where the gene is found in the genome"""
  lines = open(file).readlines()

  prot_list = []
  for line in lines:
    if line[0] == '>':
        prot_list.append(line[1:])

  b = [idx for i in key_word.split() for idx, j in enumerate(prot_list) if i in j]

  result1 = re.search('(?<=\#\ )(.*?)(?=\ \#)', prot_list[b[0]])
  start = int(result1.group(1))
  result2 = re.search('(?:.*?\#\ ){2}(.*?)(?=\ \#)', prot_list[b[0]])
  end = int(result2.group(1))

  return start-1, end-1

##########################################################################################################################################################################

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

##########################################################################################################################################################################

def get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record, temp):
	"""This function determines the border (start and end coordinates) of a PICI/phage satellite by looking for paired attachement sites. If not attachment sites are 
	found then it will resort to a lower-quality trim."""

	print("Finding  border...")
	quality = ""

	# If PICI is in the forward orientation
	if orientation == "+":
		print("+ orientation")

		# Get your head and tail sequence (regions in front of int and after prirep)
		head_seq = SeqRecord(record.seq[pici_low_limit:int_location_start], id="head_seq"); 
		tail_seq = SeqRecord(record.seq[prirep_location_end:pici_high_limit], id="tail_seq");

		#Write two sequences to files

		SeqIO.write(head_seq, temp + "/head_seq", "fasta"); 
		SeqIO.write(tail_seq, temp + "/tail_seq", "fasta")

		# Run BLAST
		setupDB = "formatdb -p F -i " + temp + "/tail_seq -n " + temp + "/tail_seq -o T"
		cmd = "blastall -p blastn -d " + temp + "/tail_seq -i " + temp + "/head_seq -o " + temp + "/BLAST_for_border -m 8 -S 3 -v 4 -b 4 -a 2 -F T"
		os.system(setupDB)
		os.system(cmd)

		# Create range for the PICI general location from the integrase trimmer script
		seq_trim_start = int(record.description.split(';')[2])
		seq_trim_end = int(record.description.split(';')[3])
		trim_list = [*range(int(seq_trim_start), int(seq_trim_end)+1, 1)]

		if os.stat(temp + '/BLAST_for_border').st_size == 0:
			print("No suitable attachment site found. Using a lower-quality border.")

			final_seq_start = trim_list[pici_low_limit]
			border_start = pici_low_limit
			final_seq_end = trim_list[pici_high_limit]
			border_end = pici_high_limit
			attL_start = "NA"
			attL_end = "NA"
			attR_start = "NA"
			attR_end = "NA"
			attC = "NA"
			attC_L = "NA"
			attC_R = "NA"

			#return final_start, final_end, trim_list, border_start, border_end, attL_start, attL_end, attR_start, attR_end, attC, quality
			# taken care of: finalstart/end, borderstart/end, att stuff, quality
			############################################
		else:
			print("Potential attatchment sites found...")

			border_df = pd.read_csv(temp + '/BLAST_for_border', sep='\t', header = None)

			if (pici_low_limit + int(border_df.iloc[0,6])-1) < 0:
				final_seq_start = trim_list[pici_low_limit]
				border_start = pici_low_limit
				attL_start = "NA"
				attL_end = "NA"
				attC_L = "NA"
			else:
				final_seq_start = trim_list[pici_low_limit + int(border_df.iloc[0,6])-1]
				border_start = pici_low_limit + int(border_df.iloc[0,6])-1
				attL_start = trim_list[pici_low_limit + int(border_df.iloc[0,6])-1]
				attL_end = trim_list[pici_low_limit + int(border_df.iloc[0,7])-1]
				attC_L = record.seq[(pici_low_limit + int(border_df.iloc[0,6])-1):(pici_low_limit + int(border_df.iloc[0,7])-1)]

			if (prirep_location_end + int(border_df.iloc[0,9])-1) > pici_high_limit:
				final_seq_end = trim_list[pici_high_limit]
				border_end = pici_high_limit
				attR_start = "NA"
				attR_end = "NA"
				attC_R = "NA"
			else:
				final_seq_end = trim_list[prirep_location_end + int(border_df.iloc[0,9])-1]
				border_end = prirep_location_end + int(border_df.iloc[0,9])-1
				attR_start = trim_list[prirep_location_end + int(border_df.iloc[0,8])-1]
				attR_end = trim_list[prirep_location_end + int(border_df.iloc[0,9])-1]
				attC_R = record.seq[(prirep_location_end + int(border_df.iloc[0,8])-1):(prirep_location_end + int(border_df.iloc[0,9])-1)]

		# register attC
		if attC_L != attC_R:
			attC = str(attC_L + "/" + attC_R)
		else:
			attC = attC_L

		# determine quality
		if attC_L and attC_R == "NA":
			quality = "low"
		elif attC_L == "NA" and attC_R != "NA":
			quality = "medium"
		elif attC_L != "NA" and attC_R == "NA":
			quality = "medium"
		else:
			quality = "high"

		#print("Quality:", quality) #quality of trim
		#print("attC:", attC) #attC sequence
		#print("attC_L:", attC_L) #attL sequence
		#print("attC_R:", attC_R) #attR sequence
		#print("attL start", attL_start) #attL start
		#print("attL end", attL_end) #attL end
		#print("attR start", attR_start) #attR start
		#print("attR end", attR_end) #attR end
		#print("seq trim start:", seq_trim_start) #original trim start coordinates
		#print("seq trim end:", seq_trim_end) #original trim end coordinates
		#print("border start:", border_start) #border index for sequence start
		#print("border end:", border_end) #border index for sequence end
		#print("FINAL SEQ START:", final_seq_start) #new start coordinates
		#print("FINAL SEQ END:", final_seq_end) #new end coordinates

		#attL_start = trim_list[pici_low_limit + int(border_df.iloc[0,6])
		#attL_end = int(border_df.iloc[0,7])
		#attR_start = int(border_df.iloc[0,8])
		#attR_end = int(border_df.iloc[0,9])
		#attC = str()

		return quality, attC, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end


	# If PICI is in the backward orientation
	elif orientation == "-":

		print("- orientation")

		# Get your head and tail sequence (regions in front and after integrase)
		head_seq = SeqRecord(record.seq[int_location_end:pici_high_limit], id="head_seq")
		tail_seq = SeqRecord(record.seq[pici_low_limit:prirep_location_start], id="tail_seq")

		#Write two sequences to files
		SeqIO.write(head_seq, temp + "/head_seq", "fasta")
		SeqIO.write(tail_seq, temp + "/tail_seq", "fasta")

		# Run BLAST
		setupDB = "formatdb -p F -i "+ temp + "/tail_seq -n tail_seq -o T"
		cmd = "blastall -p blastn -d "+ temp + "/tail_seq -i "+ temp + "/head_seq -o "+ temp + "/BLAST_for_border -m 8 -S 3 -v 4 -b 4 -a 2 -F T"
		os.system(setupDB)
		os.system(cmd)

		# Create range for the PICI general location from the integrase trimmer script
		seq_trim_start = int(record.description.split(';')[2])
		seq_trim_end = int(record.description.split(';')[3])
		trim_list = [*range(int(seq_trim_start), int(seq_trim_end)+1, 1)]

		if os.stat(temp + '/BLAST_for_border').st_size == 0:
			print("No suitable attachment site found. Using a lower-quality border.")

			final_seq_start = trim_list[pici_low_limit]
			border_start = pici_low_limit
			final_seq_end = trim_list[pici_high_limit]
			border_end = pici_high_limit
			attL_start = "NA"
			attL_end = "NA"
			attR_start = "NA"
			attR_end = "NA"
			attC = "NA"
			attC_L = "NA"
			attC_R = "NA"

		else:
			print("Potential attatchment sites found...")

			border_df = pd.read_csv(temp + '/BLAST_for_border', sep='\t', header = None)

			if (pici_low_limit + int(border_df.iloc[0,8])-1) < 0:
				final_seq_start = trim_list[pici_low_limit]
				border_start = pici_low_limit
				attL_start = "NA"
				attL_end = "NA"
				attC_L = "NA"
			else:
				final_seq_start = trim_list[pici_low_limit + int(border_df.iloc[0,8])-1]
				border_start = pici_low_limit + int(border_df.iloc[0,8])-1
				attL_start = trim_list[pici_low_limit + int(border_df.iloc[0,8])-1]
				attL_end = trim_list[pici_low_limit + int(border_df.iloc[0,9])-1]
				attC_L = record.seq[(pici_low_limit + int(border_df.iloc[0,8])-1):(pici_low_limit + int(border_df.iloc[0,9])-1)]
			if (int_location_end + int(border_df.iloc[0,7])-1) > pici_high_limit:
				final_seq_end = trim_list[pici_high_limit]
				border_end = pici_high_limit
				attR_start = "NA"
				attR_end = "NA"
				attC_R = "NA"
			else:
				print("ERROR REPORT:", (int_location_end + int(border_df.iloc[0,7])-1), pici_high_limit, len(trim_list)-1)
				final_seq_end = trim_list[int_location_end + int(border_df.iloc[0,7])-1]
				border_end = int_location_end + int(border_df.iloc[0,7])-1
				attR_start = trim_list[int_location_end + int(border_df.iloc[0,6])-1]
				attR_end = trim_list[int_location_end + int(border_df.iloc[0,7])-1]
				attC_R = record.seq[(int_location_end + int(border_df.iloc[0,6])-1):(int_location_end + int(border_df.iloc[0,7])-1)]

		# register attC
		if attC_L == "":
			attC_L = "NA"
		if attC_R == "":
			attC_R = "NA"
		if attC_L != attC_R:
			attC = str(attC_L + "/" + attC_R)
		else:
			attC = attC_L

		# determine quality
		if attC_L and attC_R == "NA":
			quality = "low"
		elif attC_L == "NA" and attC_R != "NA":
			quality = "medium"
		elif attC_L != "NA" and attC_R == "NA":
			quality = "medium"
		else:
			quality = "high"

		#print("Quality:", quality)
		#print("attC:", attC)
		#print("attC_L:", attC_L)
		#print("attC_R:", attC_R)
		#print("attL start", attL_start)
		#print("attL end", attL_end)
		#print("attR start", attR_start)
		#print("attR end", attR_end)
		#print("seq trim start:", seq_trim_start)
		#print("seq trim end:", seq_trim_end)
		#print("border start:", border_start)
		#print("border end:", border_end)

		#print("FINAL SEQ START:", final_seq_start)
		#print("FINAL SEQ END:", final_seq_end)

		return quality, attC, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end

##########################################################################################################################################################################

# lists to hold PICI headers and sequences for fasta output
name_list = []
seq_list = []
PICI_type = []
desc_list = []
final_start = []
final_end = []

PICI_length_list = []
PICI_start_coord = []
PICI_end_coord = []
orientation_list = []
quality_list = []
int_ID_list = []
int_identity_list = []
int_start_coord = []
int_end_coord = []
alpA_ID_list = []
alpA_identity_list = []
alpA_start_coord = []
alpA_end_coord = []
sis_ID_list = []
sis_identity_list = []
sis_start_coord = []
sis_end_coord = []
prirep_ID_list = []
prirep_identity_list = []
prirep_start_coord = []
prirep_end_coord = []
attL_start_coord = []
attL_end_coord = []
attR_start_coord = []
attR_end_coord = []
attC_list = []
terS_list = []
rpp_list = []
ppi_list = []

master_list = []

# Read in BLAST_Db to get last unique ID num
#BLAST_file = read.table(file = 'BLAST_DB.tsv', sep = '\t', header = FALSE)
#if os.stat(BLAST_file).st_size == 0:
#	master_db_int_num = 0
#	master_db_alpA_num = 0
#	master_db_sis_num = 0
#	master_db_prirep_num = 0
#	master_db_terS_num = 0
#	master_db_ppi_num = 0
#	master_db_rpp_num = 0
#else:
#	int_idx = BLAST_file.iloc.where


#master_db_seq = ""


fasta_file = args.fasta
with open(fasta_file, mode='r') as handle:
  
  for record in SeqIO.parse(handle, 'fasta'):
    for i in range(len(condensed_df)):
      if condensed_df.iloc[i,0].startswith(record.id): # ensures same sequence in multi-sequence fasta file is being used
        if condensed_df.iloc[i,1].startswith('int'): # find integrases in df
          if condensed_df.iloc[i,2] >= args.i: # check it integrase identity is above 90%
            print('Integrase found with identity >= 70%:')
            print(condensed_df.iloc[i,0].split('||')[0])
            
            # take the protein number of the int ("0" column)
            int_prot_num = int((condensed_df.iloc[i,0]).split('_')[-1])
            
            #set lower and upper bounds for alpA, sis, and pri-rep to be found (to the left and right of int)
            alpA_low_lim_forward = int_prot_num + 1
            alpA_high_lim_forward = int_prot_num + 5 #compensating for range() function - highest number in range() not considered... essentially adding 1 to high limits to get balanced range -/+ 
            alpA_low_lim_backward = int_prot_num - 0
            alpA_high_lim_backward = int_prot_num - 4
            sis_low_lim_forward = int_prot_num + 4
            sis_high_lim_forward = int_prot_num + 13
            sis_low_lim_backward = int_prot_num - 3
            sis_high_lim_backward = int_prot_num - 12
            prirep_low_lim_forward = int_prot_num + 4
            prirep_high_lim_forward = int_prot_num + 16
            prirep_low_lim_backward = int_prot_num - 3
            prirep_high_lim_backward = int_prot_num - 15

            # To stop current iteration if PICI is found (stop_loop = True)
            stop_loop = False
            # Reset directionality
            forward_direction = False
            backward_direction = False


            #####################################################
            # FIRST search for co-localized alpA... determines G-
            ##################################################### 
            for k in range(i-4,i+4):
              # if the index is negative set it to 0
              if k < 0:
                k = 0
              #breaks loop if the index exceeds the dataframe length
              if k >= condensed_df.shape[0]:
                break
              # otherwise get the gene number for alpA
              elif condensed_df.iloc[k,0].startswith(record.id):  
                if condensed_df.iloc[k,1].startswith('alpA'):
                  if condensed_df.iloc[k,2] >= args.a: # check it alpA identity is above 50%
                    alpA_prot_num = int((condensed_df.iloc[k,0]).split('_')[-1])

                    # check to see if it is within 1-4 genes away from the integrase  
                    if alpA_prot_num in range(alpA_high_lim_backward, alpA_low_lim_backward) or alpA_prot_num in range(alpA_low_lim_forward, alpA_high_lim_forward):
                    
                      # Determine directionality of PICI
                      if alpA_prot_num > int_prot_num:
                        forward_direction = True
                      elif alpA_prot_num < int_prot_num:
                        backward_direction = True

                      # Get protein locations in host
                      # get the location of alpA and see if pri-rep is within 25 kb to left or right (50kb total)
                      int_location_start, int_location_end = get_locations(args.aa, condensed_df.iloc[i,0])
                      alpA_location_start, alpA_location_end = get_locations(args.aa, condensed_df.iloc[k,0])
                      print('\nBeginning G- PICI determination...')
                      print('Int ({}) start: {}'.format(int_prot_num, int_location_start))
                      print('Int ({}) end: {}'.format(int_prot_num, int_location_end))
                      print('alpA ({}) start: {}'.format(alpA_prot_num, alpA_location_start))
                      print('alpA ({}) end: {}'.format(alpA_prot_num, alpA_location_end))


                      # set range for pri-rep to be found from alpA
                      prirep_range_low = alpA_location_start - 25000
                      prirep_range_high = alpA_location_end + 25000

                      # fix lower limit index to 0 if it is negative
                      if prirep_range_low < 0:
                        prirep_range_low = 0

                      # fix upper limit index to the length of the sequence if it is over
                      if prirep_range_high >= len(record.seq)-1:
                        prirep_range_high = len(record.seq)-1

                      # If the forward direction is true for PICI structure
                      if forward_direction == True:
                        orientation = "+"
                        # Looking for pri-rep in the direction alpA is (forward)
                        range_to_be_discovered = range(int_location_end, prirep_range_high)

                        # check for pri-rep and if it is in the 50kb range
                        for l in range(len(condensed_df)):
                          if condensed_df.iloc[l,0].startswith(record.id):
                            if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                              prirep_location_start, prirep_location_end = get_locations(args.aa, condensed_df.iloc[l,0])
                              prirep_range = range(prirep_location_start,prirep_location_end)
                              prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                              if range_subset(prirep_range, range_to_be_discovered) == True:
                                if range_subset(range(alpA_location_start, alpA_location_end), range(int_location_end, prirep_location_start)) == True:
                                  #get PICI location using -2kb and +25kb from int
                                  pici_low_limit = int_location_start - 2000
                                  pici_high_limit = int_location_end + 25000
                                  # fix lower limit index to 0 if it is negative
                                  if pici_low_limit < 0:
                                    pici_low_limit = 0
                                  # fix upper limit index to the length of the sequence if it is over
                                  if pici_high_limit >= len(record.seq)-1:
                                    pici_high_limit = len(record.seq)-1


                                  print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                                  print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                                  print('PICI start:', pici_low_limit)
                                  print('PICI end:', pici_high_limit)
                                  print('\nPICI sequence:')
                                  print(record.seq[pici_low_limit:pici_high_limit])
                                  print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                                  print('\n')



                                  # add PICI to list
                                  trim_start = int(record.description.split(';')[2])
                                  trim_end = int(record.description.split(';')[3])
                                  trim_list = [*range(int(trim_start), int(trim_end), 1)]
                                  if pici_high_limit > len(trim_list)-1:
                                    continue
                                  else:
                                    quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record, args.temp)


                                    # Look for things like terS, rpp, and ppi in the new region
                                    # 0 = not present, 1 = present
                                    terS_present = 0
                                    ppi_present = 0
                                    rpp_present = 0

                                    #terS
                                    for r in range(len(condensed_df)):
                                      if condensed_df.iloc[r,0].startswith(record.id):
                                        if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                          terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                          terS_range = range(terS_location_start,terS_location_end)

                                          if range_subset(terS_range, range(border_start, border_end)) == True:
                                            terS_present = 1

                                    #ppi
                                    for s in range(len(condensed_df)):
                                      if condensed_df.iloc[s,0].startswith(record.id):
                                        if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                          ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                          ppi_range = range(ppi_location_start,ppi_location_end)

                                          if range_subset(ppi_range, range(border_start, border_end)) == True:
                                            ppi_present = 1

                                    #rpp
                                    for t in range(len(condensed_df)):
                                      if condensed_df.iloc[t,0].startswith(record.id):
                                        if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                          rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                          rpp_range = range(rpp_location_start,rpp_location_end)

                                          if range_subset(rpp_range, range(border_start, border_end)) == True:
                                            rpp_present = 1

                                    # Append to lists
                                    trim_start = int(record.description.split(';')[2])
                                    trim_end = int(record.description.split(';')[3])
                                    trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                                    record.description = record.description.replace(str(trim_start), str(final_seq_start))
                                    record.description = record.description.replace(str(trim_end), str(final_seq_end))

                                    #desc_list.append(record.description)
                                    #PICI_type.append("G_neg_PICI")
                                    #PICI_length_list.append(border_end - border_start)
                                    #PICI_start_coord.append(final_seq_start)
                                    #PICI_end_coord.append(final_seq_end)
                                    #orientation_list.append(orientation)
                                    #quality_list.append(quality)
                                    #int_ID_list.append(condensed_df.iloc[i,1])
                                    #int_identity_list.append(condensed_df.iloc[i,2])
                                    #int_start_coord.append(int_location_start)
                                    #int_end_coord.append(int_location_end)
                                    #alpA_ID_list.append(condensed_df.iloc[k,1])
                                    #alpA_identity_list.append(condensed_df.iloc[k,2])
                                    #alpA_start_coord.append(alpA_location_start)
                                    #alpA_end_coord.append(alpA_location_end)
                                    #sis_ID_list.append("NA")
                                    #sis_identity_list.append("NA")
                                    #sis_start_coord.append("NA")
                                    #sis_end_coord.append("NA")
                                    #prirep_ID_list.append(condensed_df.iloc[l,1])
                                    #prirep_identity_list.append(condensed_df.iloc[l,2])
                                    #prirep_start_coord.append(prirep_location_start)
                                    #prirep_end_coord.append(prirep_location_end)
                                    #attL_start_coord.append(attL_start)
                                    #attL_end_coord.append(attL_end)
                                    #attR_start_coord.append(attR_start)
                                    #attR_end_coord.append(attR_end)
                                    #attC_list.append(att_C)

                                    seq_list.append(record.seq[border_start:border_end])

                                    master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "G_neg_PICI" + ";" +
                                    str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                                    ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + str(condensed_df.iloc[k,1]) + ";" +
                                    str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(alpA_location_start)) + ";" + str(int(trim_start) + int(alpA_location_end)) + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" +
                                    ";" + str(condensed_df.iloc[l,1]) + ";" + str(condensed_df.iloc[l,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                                    ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                                    str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                                    #print("MASTER_LIST:", master_list)
                                    #name_list.append(record.id)
                                    #seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                                    #PICI_type.append("G_neg_PICI")
                                    #trim_start = int(record.description.split(';')[2])
                                    #trim_end = int(record.description.split(';')[3])
                                    #trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                    #final_start = trim_list[pici_low_limit]
                                    #final_end = trim_list[pici_high_limit-1]
                                    #record.description = record.description.replace(str(trim_start), str(final_start))
                                    #record.description = record.description.replace(str(trim_end), str(final_end))
                                    #desc_list.append(record.description)
                                    stop_loop = True


                      # If the backward direction is true for PICI structure
                      elif backward_direction == True:
                        orientation = "-"
                        # Looking for pri-rep in the direction alpA is (backward)
                        range_to_be_discovered = range(prirep_range_low, int_location_start)

                        # check for pri-rep and if it is in the 50kb range
                        for l in range(len(condensed_df)):
                          if condensed_df.iloc[l,0].startswith(record.id):
                            if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                              prirep_location_start, prirep_location_end = get_locations(args.aa, condensed_df.iloc[l,0])
                              prirep_range = range(prirep_location_start,prirep_location_end)
                              prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                              if range_subset(prirep_range, range_to_be_discovered) == True:
                                if range_subset(range(alpA_location_start, alpA_location_end), range(prirep_location_end, int_location_start)) == True:
                                  #get PICI location using -25kb and +2kb from int
                                  pici_low_limit = int_location_start - 25000
                                  pici_high_limit = int_location_end + 2000
                                  # fix lower limit index to 0 if it is negative
                                  if pici_low_limit < 0:
                                    pici_low_limit = 0
                                  # fix upper limit index to the length of the sequence if it is over
                                  if pici_high_limit >= len(record.seq)-1:
                                    pici_high_limit = len(record.seq)-1


                                  print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                                  print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                                  print('PICI start:', pici_low_limit)
                                  print('PICI end:', pici_high_limit)
                                  print('\nPICI sequence:')
                                  print(record.seq[pici_low_limit:pici_high_limit])
                                  print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                                  print('\n')


                                  ##################
                                  #get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)


                                  # add PICI to list
                                  trim_start = int(record.description.split(';')[2])
                                  trim_end = int(record.description.split(';')[3])
                                  trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                  if pici_high_limit > len(trim_list)-1:
                                    continue
                                  else:
                                    quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)

                                    # Look for things like terS, rpp, and ppi in the new region
                                    # 0 = not present, 1 = present
                                    terS_present = 0
                                    ppi_present = 0
                                    rpp_present = 0

                                    #terS
                                    for r in range(len(condensed_df)):
                                      if condensed_df.iloc[r,0].startswith(record.id):
                                        if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                          terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                          terS_range = range(terS_location_start,terS_location_end)

                                          if range_subset(terS_range, range(border_start, border_end)) == True:
                                            terS_present = 1

                                    #ppi
                                    for s in range(len(condensed_df)):
                                      if condensed_df.iloc[s,0].startswith(record.id):
                                        if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                          ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                          ppi_range = range(ppi_location_start,ppi_location_end)

                                          if range_subset(ppi_range, range(border_start, border_end)) == True:
                                            ppi_present = 1

                                    #rpp
                                    for t in range(len(condensed_df)):
                                      if condensed_df.iloc[t,0].startswith(record.id):
                                        if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                          rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                          rpp_range = range(rpp_location_start,rpp_location_end)

                                          if range_subset(rpp_range, range(border_start, border_end)) == True:
                                            rpp_present = 1

                                    # Append to lists
                                    trim_start = int(record.description.split(';')[2])
                                    trim_end = int(record.description.split(';')[3])
                                    trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                                    record.description = record.description.replace(str(trim_start), str(final_seq_start))
                                    record.description = record.description.replace(str(trim_end), str(final_seq_end))

                                    seq_list.append(record.seq[border_start:border_end])

                                    master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "G_neg_PICI" + ";" +
                                    str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                                    ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + str(condensed_df.iloc[k,1]) + ";" +
                                    str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(alpA_location_start)) + ";" + str(int(trim_start) + int(alpA_location_end)) + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" +
                                    ";" + str(condensed_df.iloc[l,1]) + ";" + str(condensed_df.iloc[l,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                                    ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                                    str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                                    #print(master_list)
                                    #name_list.append(record.id)

                                    stop_loop = True




            #######################################################
            # SECOND search for co-localized sis... determines SaPI
            ####################################################### 
            if stop_loop == False:
              for k in range(i-12,i+12):
                # if the index is negative set it to 0
                if k < 0:
                  k = 0
                #breaks loop if the index exceeds the dataframe length
                if k >= condensed_df.shape[0]:
                  break
                # otherwise get the gene number for sis  
                elif condensed_df.iloc[k,0].startswith(record.id):
                  if condensed_df.iloc[k,1].startswith('sis'):
                    sis_prot_num = int((condensed_df.iloc[k,0]).split('_')[-1])

                    # check to see if it is within 4-12 genes away from the integrase  
                    if sis_prot_num in range(sis_low_lim_forward, sis_high_lim_forward) or sis_prot_num in range(sis_high_lim_backward, sis_low_lim_backward):

                      # Reset and determine directionality of PICI
                      forward_direction = False
                      backward_direction = False
                      if sis_prot_num > int_prot_num:
                        forward_direction = True
                      elif sis_prot_num < int_prot_num:
                        backward_direction = True

                      # Get protein locations in host
                      # get the location of sis and see if pri-rep is within 25 kb to left or right (50kb total)
                      int_location_start, int_location_end = get_locations(args.aa, condensed_df.iloc[i,0])
                      sis_location_start, sis_location_end = get_locations(args.aa, condensed_df.iloc[k,0])
                      print('\nBeginning SaPI determination...')
                      print('Int ({}) start: {}'.format(int_prot_num, int_location_start))
                      print('Int ({}) end: {}'.format(int_prot_num, int_location_end))
                      print('Sis ({}) start: {}'.format(sis_prot_num, sis_location_start))
                      print('Sis ({}) end: {}'.format(sis_prot_num, sis_location_end))


                      # set range for pri-rep to be found from sis
                      prirep_range_low = sis_location_start - 25000
                      prirep_range_high = sis_location_end + 25000

                      # fix lower limit index to 0 if it is negative
                      if prirep_range_low < 0:
                        prirep_range_low = 0

                      # fix upper limit index to the length of the sequence if it is over
                      if prirep_range_high >= len(record.seq)-1:
                        prirep_range_high = len(record.seq)-1








                      # If the forward direction is true for PICI structure
                      if forward_direction == True:
                        orientation = "+"
                        # Looking for pri-rep in the direction sis is (forward)
                        range_to_be_discovered = range(int_location_end, prirep_range_high)

                        # check for pri-rep and if it is in range
                        for l in range(len(condensed_df)):
                          if condensed_df.iloc[l,0].startswith(record.id):
                            if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                              prirep_location_start, prirep_location_end = get_locations(args.aa, condensed_df.iloc[l,0])
                              prirep_range = range(prirep_location_start,prirep_location_end)
                              prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                              # Make sure structure is correct
                              if range_subset(prirep_range, range_to_be_discovered) == True:
                                if range_subset(range(sis_location_start, sis_location_end), range(int_location_end, prirep_location_start)) == True:
                                  #get PICI location using -2kb and +25kb from int
                                  pici_low_limit = int_location_start - 2000
                                  pici_high_limit = int_location_end + 25000
                                  # fix lower limit index to 0 if it is negative
                                  if pici_low_limit < 0:
                                    pici_low_limit = 0
                                  # fix upper limit index to the length of the sequence if it is over
                                  if pici_high_limit >= len(record.seq)-1:
                                    pici_high_limit = len(record.seq)-1


                                  print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                                  print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                                  print('PICI start:', pici_low_limit)
                                  print('PICI end:', pici_high_limit)
                                  print('\nPICI sequence:')
                                  print(record.seq[pici_low_limit:pici_high_limit])
                                  print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                                  print('\n')



                                  #############
                                  #get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)

                                  # add PICI to list
                                  print(record.description.split(';'))
                                  trim_start = int(record.description.split(';')[2])
                                  trim_end = int(record.description.split(';')[3])
                                  trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                  if pici_high_limit > len(trim_list)-1:
                                    continue
                                  else:
                                    quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record, args.temp)

                                    # Look for things like terS, rpp, and ppi in the new region
                                    # 0 = not present, 1 = present
                                    terS_present = 0
                                    ppi_present = 0
                                    rpp_present = 0

                                    #terS
                                    for r in range(len(condensed_df)):
                                      if condensed_df.iloc[r,0].startswith(record.id):
                                        if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                          terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                          terS_range = range(terS_location_start,terS_location_end)

                                          if range_subset(terS_range, range(border_start, border_end)) == True:
                                            terS_present = 1

                                    #ppi
                                    for s in range(len(condensed_df)):
                                      if condensed_df.iloc[s,0].startswith(record.id):
                                        if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                          ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                          ppi_range = range(ppi_location_start,ppi_location_end)

                                          if range_subset(ppi_range, range(border_start, border_end)) == True:
                                            ppi_present = 1

                                    #rpp
                                    for t in range(len(condensed_df)):
                                      if condensed_df.iloc[t,0].startswith(record.id):
                                        if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                          rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                          rpp_range = range(rpp_location_start,rpp_location_end)

                                          if range_subset(rpp_range, range(border_start, border_end)) == True:
                                            rpp_present = 1

                                    # Append to lists
                                    #name_list.append(record.id)
                                    #seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                                    #PICI_type.append("SaPI")
                                    #trim_start = int(record.description.split(';')[2])
                                    #trim_end = int(record.description.split(';')[3])
                                    #trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                    #final_start = trim_list[pici_low_limit]
                                    #final_end = trim_list[pici_high_limit-1]
                                    #record.description = record.description.replace(str(trim_start), str(final_start))
                                    #record.description = record.description.replace(str(trim_end), str(final_end))
                                    #desc_list.append(record.description)

                                    trim_start = int(record.description.split(';')[2])
                                    trim_end = int(record.description.split(';')[3])
                                    trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                                    record.description = record.description.replace(str(trim_start), str(final_seq_start))
                                    record.description = record.description.replace(str(trim_end), str(final_seq_end))

                                    seq_list.append(record.seq[border_start:border_end])

                                    master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "SaPI" + ";" +
                                    str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                                    ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + "NA" + ";" +
                                    "NA" + ";" + "NA" + ";" + "NA" + ";" + str(condensed_df.iloc[k,1]) + ";" + str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(sis_location_start)) + ";" + str(int(trim_start) + int(sis_location_end)) +
                                    ";" + str(condensed_df.iloc[l,1]) + ";" + str(condensed_df.iloc[l,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                                    ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                                    str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                                    #name_list.append(record.id)

                                    stop_loop = True


                      # If the backward direction is true for PICI structure
                      elif backward_direction == True:
                        orientation = "-"
                        # Looking for pri-rep in the direction sis is (backward)
                        range_to_be_discovered = range(prirep_range_low, int_location_start)

                        # check for pri-rep and if it is in the 50kb range
                        for l in range(len(condensed_df)):
                          if condensed_df.iloc[l,0].startswith(record.id):
                            if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                              prirep_location_start, prirep_location_end = get_locations(args.aa, condensed_df.iloc[l,0])
                              prirep_range = range(prirep_location_start,prirep_location_end)
                              prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                              # Make sure structure is correct
                              if range_subset(prirep_range, range_to_be_discovered) == True:
                                print("Range is a subset")
                                if range_subset(range(sis_location_start, sis_location_end), range(prirep_location_end, int_location_start)) == True:
                                  #get PICI location using -25kb and +2kb from int
                                  pici_low_limit = int_location_start - 25000
                                  pici_high_limit = int_location_end + 2000
                                  # fix lower limit index to 0 if it is negative
                                  if pici_low_limit < 0:
                                    pici_low_limit = 0
                                  # fix upper limit index to the length of the sequence if it is over
                                  if pici_high_limit >= len(record.seq)-1:
                                    pici_high_limit = len(record.seq)-1


                                  print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                                  print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                                  print('PICI start:', pici_low_limit)
                                  print('PICI end:', pici_high_limit)
                                  print('\nPICI sequence:')
                                  print(record.seq[pici_low_limit:pici_high_limit])
                                  print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                                  print('\n')



                                  #############
                                  #get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)

                                  # add PICI to list
                                  trim_start = int(record.description.split(';')[2])
                                  trim_end = int(record.description.split(';')[3])
                                  trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                  if pici_high_limit > len(trim_list)-1:
                                    continue
                                  else:
                                    quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)

                                    # Look for things like terS, rpp, and ppi in the new region
                                    # 0 = not present, 1 = present
                                    terS_present = 0
                                    ppi_present = 0
                                    rpp_present = 0

                                    #terS
                                    for r in range(len(condensed_df)):
                                      if condensed_df.iloc[r,0].startswith(record.id):
                                        if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                          terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                          terS_range = range(terS_location_start,terS_location_end)

                                          if range_subset(terS_range, range(border_start, border_end)) == True:
                                            terS_present = 1

                                    #ppi
                                    for s in range(len(condensed_df)):
                                      if condensed_df.iloc[s,0].startswith(record.id):
                                        if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                          ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                          ppi_range = range(ppi_location_start,ppi_location_end)

                                          if range_subset(ppi_range, range(border_start, border_end)) == True:
                                            ppi_present = 1

                                    #rpp
                                    for t in range(len(condensed_df)):
                                      if condensed_df.iloc[t,0].startswith(record.id):
                                        if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                          rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                          rpp_range = range(rpp_location_start,rpp_location_end)

                                          if range_subset(rpp_range, range(border_start, border_end)) == True:
                                            rpp_present = 1

                                    # Append to lists
                                    #name_list.append(record.id)
                                    #seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                                    #PICI_type.append("SaPI")
                                    #trim_start = int(record.description.split(';')[2])
                                    #trim_end = int(record.description.split(';')[3])
                                    #trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                                    #final_start = trim_list[pici_low_limit]
                                    #final_end = trim_list[pici_high_limit-1]
                                    #record.description = record.description.replace(str(trim_start), str(final_start))
                                    #record.description = record.description.replace(str(trim_end), str(final_end))
                                    #desc_list.append(record.description)
                                    #stop_loop = True

                                    trim_start = int(record.description.split(';')[2])
                                    trim_end = int(record.description.split(';')[3])
                                    trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                                    record.description = record.description.replace(str(trim_start), str(final_seq_start))
                                    record.description = record.description.replace(str(trim_end), str(final_seq_end))

                                    seq_list.append(record.seq[border_start:border_end])

                                    master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "SaPI" + ";" +
                                    str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                                    ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + "NA" + ";" +
                                    "NA" + ";" + "NA" + ";" + "NA" + ";" + str(condensed_df.iloc[k,1]) + ";" + str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(sis_location_start)) + ";" + str(int(trim_start) + int(sis_location_end)) +
                                    ";" + str(condensed_df.iloc[l,1]) + ";" + str(condensed_df.iloc[l,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                                    ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                                    str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                                    #name_list.append(record.id)

                                    stop_loop = True


            ######################################################################
            # THIRD search for co-localized pri-rep... determines phage satellites
            ######################################################################
            if stop_loop == False:
              for k in range(i-15,i+15):
                # if the index is negative set it to 0
                if k < 0:
                  k = 0
                #breaks loop if the index exceeds the dataframe length
                if k >= condensed_df.shape[0]:
                  break
                # otherwise get the gene number for pri-rep
                elif condensed_df.iloc[k,0].startswith(record.id):
                  if condensed_df.iloc[k,1].startswith('pri-rep'):
                    prirep_prot_num = int((condensed_df.iloc[k,0]).split('_')[-1])
                    # check to see if it is within 4-15 genes away from the integrase  
                    if prirep_prot_num in range(prirep_low_lim_forward, prirep_high_lim_forward) or prirep_prot_num in range(prirep_high_lim_backward, prirep_low_lim_backward):

                      # Resent and determine directionality of phage satellite
                      forward_direction = False
                      backward_direction = False
                      if prirep_prot_num > int_prot_num:
                        forward_direction = True
                      elif prirep_prot_num < int_prot_num:
                        backward_direction = True

                      # Get protein locations in host
                      # get the location of pri-rep and see if pri-rep is within 25 kb to left or right (50kb total)
                      int_location_start, int_location_end = get_locations(args.aa, condensed_df.iloc[i,0])
                      prirep_location_start, prirep_location_end = get_locations(args.aa, condensed_df.iloc[k,0])
                      print('\nBeginning phage satellite determination...')
                      print('Int ({}) start: {}'.format(int_prot_num, int_location_start))
                      print('Int ({}) end: {}'.format(int_prot_num, int_location_end))
                      print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                      print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))



                      # If the forward direction is true for phage satellite structure
                      if forward_direction == True:
                        orientation = "+"
                        # set range for pri-rep to be found 
                        pici_low_limit = int_location_start - 2000
                        pici_high_limit = int_location_end + 25000

                        # fix lower limit index to 0 if it is negative
                        if pici_low_limit < 0:
                          pici_low_limit = 0

                        # fix upper limit index to the length of the sequence if it is over
                        if pici_high_limit >= len(record.seq)-1:
                          pici_high_limit = len(record.seq)-1

                        print('PICI start:', pici_low_limit)
                        print('PICI end:', pici_high_limit)
                        print('\nPICI sequence:')
                        print(record.seq[pici_low_limit:pici_high_limit])
                        print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                        print('\n')




                        #############
                        #get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)


                        # add PICI to list
                        trim_start = int(record.description.split(';')[2])
                        trim_end = int(record.description.split(';')[3])
                        trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                        if pici_high_limit > len(trim_list)-1:
                          continue
                        else:
                          quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)

                          # Look for things like terS, rpp, and ppi in the new region
                          # 0 = not present, 1 = present
                          terS_present = 0
                          ppi_present = 0
                          rpp_present = 0

                          #terS
                          for r in range(len(condensed_df)):
                            if condensed_df.iloc[r,0].startswith(record.id):
                              if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                terS_range = range(terS_location_start,terS_location_end)

                                if range_subset(terS_range, range(border_start, border_end)) == True:
                                  terS_present = 1

                          #ppi
                          for s in range(len(condensed_df)):
                            if condensed_df.iloc[s,0].startswith(record.id):
                              if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                ppi_range = range(ppi_location_start,ppi_location_end)

                                if range_subset(ppi_range, range(border_start, border_end)) == True:
                                  ppi_present = 1

                          #rpp
                            for t in range(len(condensed_df)):
                              if condensed_df.iloc[t,0].startswith(record.id):
                                if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                  rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                  rpp_range = range(rpp_location_start,rpp_location_end)

                                  if range_subset(rpp_range, range(border_start, border_end)) == True:
                                    rpp_present = 1

                          # Append to lists
                          #name_list.append(record.id)
                          #seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                          #PICI_type.append("phage_satellite")
                          #trim_start = int(record.description.split(';')[2])
                          #trim_end = int(record.description.split(';')[3])
                          #trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                          #final_start = trim_list[pici_low_limit]
                          #final_end = trim_list[pici_high_limit-1]
                          #record.description = record.description.replace(str(trim_start), str(final_start))
                          #record.description = record.description.replace(str(trim_end), str(final_end))
                          #desc_list.append(record.description)
                          #stop_loop = True

                          trim_start = int(record.description.split(';')[2])
                          trim_end = int(record.description.split(';')[3])
                          trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                          record.description = record.description.replace(str(trim_start), str(final_seq_start))
                          record.description = record.description.replace(str(trim_end), str(final_seq_end))

                          seq_list.append(record.seq[border_start:border_end])

                          master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "phage_satellite" + ";" +
                          str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                          ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + "NA" + ";" +
                          "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" +
                          ";" + str(condensed_df.iloc[k,1]) + ";" + str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                          ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                          str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                          #print(master_list)

                          #name_list.append(record.id)

                          stop_loop = True

                      # If the backward direction is true for phage satellite structure
                      elif backward_direction == True:
                        orientation = "-"
                        # set range for pri-rep to be found 
                        pici_low_limit = int_location_start - 25000
                        pici_high_limit = int_location_end + 2000

                        # fix lower limit index to 0 if it is negative
                        if pici_low_limit < 0:
                          pici_low_limit = 0

                        # fix upper limit index to the length of the sequence if it is over
                        if pici_high_limit >= len(record.seq)-1:
                          pici_high_limit = len(record.seq)-1

                        print('PICI start:', pici_low_limit)
                        print('PICI end:', pici_high_limit)
                        print('\nPICI sequence:')
                        print(record.seq[pici_low_limit:pici_high_limit])
                        print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                        print('\n')




                        #############
                        #get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record)


                        # add PICI to list
                        trim_start = int(record.description.split(';')[2])
                        trim_end = int(record.description.split(';')[3])
                        trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                        if pici_high_limit > len(trim_list)-1:
                          continue
                        else:
                          quality, att_C, attL_start, attL_end, attR_start, attR_end, border_start, border_end, final_seq_start, final_seq_end = get_PICI_border(int_location_start, int_location_end, prirep_location_start, prirep_location_end, pici_low_limit, pici_high_limit, orientation, record, args.temp)

                          # Look for things like terS, rpp, and ppi in the new region
                          # 0 = not present, 1 = present
                          terS_present = 0
                          ppi_present = 0
                          rpp_present = 0

                          #terS
                          for r in range(len(condensed_df)):
                            if condensed_df.iloc[r,0].startswith(record.id):
                              if condensed_df.iloc[r,1].startswith('terS'): # find terS in df
                                terS_location_start, terS_location_end = get_locations(args.aa, condensed_df.iloc[r,0])
                                terS_range = range(terS_location_start,terS_location_end)

                                if range_subset(terS_range, range(border_start, border_end)) == True:
                                  terS_present = 1

                          #ppi
                          for s in range(len(condensed_df)):
                            if condensed_df.iloc[s,0].startswith(record.id):
                              if condensed_df.iloc[s,1].startswith('ppi'): # find terS in df
                                ppi_location_start, ppi_location_end = get_locations(args.aa, condensed_df.iloc[s,0])
                                ppi_range = range(ppi_location_start,ppi_location_end)

                                if range_subset(ppi_range, range(border_start, border_end)) == True:
                                  ppi_present = 1

                          #rpp
                            for t in range(len(condensed_df)):
                              if condensed_df.iloc[t,0].startswith(record.id):
                                if condensed_df.iloc[t,1].startswith('rpp'): # find terS in df
                                  rpp_location_start, rpp_location_end = get_locations(args.aa, condensed_df.iloc[t,0])
                                  rpp_range = range(rpp_location_start,rpp_location_end)

                                  if range_subset(rpp_range, range(border_start, border_end)) == True:
                                    rpp_present = 1

                          # Append to lists
                          #name_list.append(record.id)
                          #seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                          #PICI_type.append("phage_satellite")
                          #trim_start = int(record.description.split(';')[2])
                          #trim_end = int(record.description.split(';')[3])
                          #trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]
                          #final_start = trim_list[pici_low_limit]
                          #final_end = trim_list[pici_high_limit-1]
                          #record.description = record.description.replace(str(trim_start), str(final_start))
                          #record.description = record.description.replace(str(trim_end), str(final_end))
                          #desc_list.append(record.description)
                          #stop_loop = True

                          trim_start = int(record.description.split(';')[2])
                          trim_end = int(record.description.split(';')[3])
                          trim_list = [*range(int(trim_start), int(trim_end)+1, 1)]

                          record.description = record.description.replace(str(trim_start), str(final_seq_start))
                          record.description = record.description.replace(str(trim_end), str(final_seq_end))

                          seq_list.append(record.seq[border_start:border_end])

                          master_list.append(str(record.description.replace(str(record.id), str(record.id).split("*")[0]).split(";")[0] + ";" + "phage_satellite" + ";" +
                          str((border_end - border_start)) + ";" + str(final_seq_start) + ";" +  str(final_seq_end) + ";" + orientation + ";" + str(condensed_df.iloc[i,1]) +
                          ";" + str(condensed_df.iloc[i,2]) + ";" + str(int(trim_start) + int(int_location_start)) + ";" + str(int(trim_start) + int(int_location_end)) + ";" + "NA" + ";" +
                          "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" + ";" + "NA" +
                          ";" + str(condensed_df.iloc[k,1]) + ";" + str(condensed_df.iloc[k,2]) + ";" + str(int(trim_start) + int(prirep_location_start)) + ";" + str(int(trim_start) + int(prirep_location_end)) +
                          ";" + str(attL_start) + ";" + str(attL_end) + ";" + str(attR_start) + ";" + str(attR_end) + ";" + str(att_C) + ";" + quality + ";" +
                          str(terS_present) + ";" + str(ppi_present) + ";" + str(rpp_present)))

                          #print(master_list)

                          #name_list.append(record.id)

                          stop_loop = True















                        
                        


            
        else:
          print('No hit')



# write PICI sequence to file
PICI_file = open(args.output, "w")

print("SEQ LIST LENGTH:", len(seq_list))
print("MASTER LIST LENGTH:", len(master_list))

for m in range(len(seq_list)):
  PICI_file.write(">"  + str(master_list[m]) +  ";" + ";" + "PICI_#" + str(m) + "\n" + str(seq_list[m]) + "\n")
  print(">"  + str(master_list[m]) +  ";" + ";" + "PICI_#" + str(m))
  #PICI_file.write(">"  + str(desc_list[m].replace(name_list[m], name_list[m].split("*")[0])) +  ";" + str(PICI_type[m]) + ";" + str(m) + "\n" + str(seq_list[m]) + "\n")
  #print(">"  + str(desc_list[m].replace(name_list[m], name_list[m].split("*")[0])) +  ";" + str(PICI_type[m]) + ";" + str(m))
#desc_list[m].replace(name_list[m], name_list[m].split("_")[0]
#str(desc_list[m])
PICI_file.close()

# sequenceID, type, length, start, end, orientation, trim_quality, int_ID, int_identity, int_start_coord, int_end coord, alpA_ID, alpA_identity, alpA_start_coord, alpA_end coord, sis_ID, sis_identity, sis_start_coord, sis_end coord, prirep_ID, prirep_identity, prirep_start_coord, prirep_end_coord, terS, rpp, ppi, attL, attR, attC seq
