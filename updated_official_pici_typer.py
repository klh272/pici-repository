import re
import pandas as pd
from Bio import SeqIO

# read in BLAST alignment output
df = pd.read_csv('BLASTp_results.out', sep='\t', header = None)
# take the max identity % (column "2") for each alignment and drop the duplicates
condensed_df = df.sort_values([2], ascending=False).drop_duplicates([0]).sort_index()
condensed_df.reset_index(drop=True, inplace=True)



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

  return start, end

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



# lists to hold PICI headers and sequences for fasta output
name_list = []
seq_list = []
PICI_type = []
desc_list = []

fasta_file = 'all.fna'
with open(fasta_file, mode='r') as handle:
  
  for record in SeqIO.parse(handle, 'fasta'):

    for i in range(len(condensed_df)):
      if condensed_df.iloc[i,0].startswith(record.id): # ensures same sequence in multi-sequence fasta file is being used
        if condensed_df.iloc[i,1].startswith('int'): # find integrases in df
          if condensed_df.iloc[i,2] >= 90: # check it integrase identity is above 90%
            print('Integrase found with identity >= 90%:')
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
                  alpA_prot_num = int((condensed_df.iloc[k,0]).split('_')[-1])

                  # check to see if it is within 1-4 genes away from the integrase  
                  if alpA_prot_num in range(alpA_high_lim_backward, alpA_low_lim_backward) or alpA_prot_num in range(alpA_low_lim_forward, alpA_high_lim_forward):

                    # Get protein locations in host
                    # get the location of alpA and see if pri-rep is within 25 kb to left or right (50kb total)
                    int_location_start, int_location_end = get_locations("all.pdg.faa", condensed_df.iloc[i,0])
                    alpA_location_start, alpA_location_end = get_locations("all.pdg.faa", condensed_df.iloc[k,0])
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
                    if prirep_range_high >= len(record.seq):
                      prirep_range_high = len(record.seq)

                    range_to_be_discovered = range(prirep_range_low, prirep_range_high)

                    # check for pri-rep and if it is in the 50kb range
                    for l in range(len(condensed_df)):
                      if condensed_df.iloc[l,0].startswith(record.id):
                        if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                          prirep_location_start, prirep_location_end = get_locations("all.pdg.faa", condensed_df.iloc[l,0])
                          prirep_range = range(prirep_location_start,prirep_location_end)
                          prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                          if range_subset(prirep_range, range_to_be_discovered) == True:
                            #get PICI location using 25kb to the left and right of pri-rep
                            pici_low_limit = prirep_location_start - 25000
                            pici_high_limit = prirep_location_end + 25000
                            # fix lower limit index to 0 if it is negative
                            if pici_low_limit < 0:
                              pici_low_limit = 0
                            # fix upper limit index to the length of the sequence if it is over
                            if pici_high_limit >= len(record.seq):
                              pici_high_limit = len(record.seq)


                            print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                            print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                            print('PICI start:', pici_low_limit)
                            print('PICI end:', pici_high_limit)
                            print('\nPICI sequence:')
                            print(record.seq[pici_low_limit:pici_high_limit])
                            print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                            print('\n')

                            # add PICI to list
                            name_list.append(record.id)
                            seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                            PICI_type.append("G_neg_PICI")
                            desc_list.append(record.description)
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

                      # Get protein locations in host
                      # get the location of sis and see if pri-rep is within 25 kb to left or right (50kb total)
                      int_location_start, int_location_end = get_locations("all.pdg.faa", condensed_df.iloc[i,0])
                      sis_location_start, sis_location_end = get_locations("all.pdg.faa", condensed_df.iloc[k,0])
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
                      if prirep_range_high >= len(record.seq):
                        prirep_range_high = len(record.seq)

                      range_to_be_discovered = range(prirep_range_low, prirep_range_high)

                      # check for pri-rep and if it is in the 50kb range
                      for l in range(len(condensed_df)):
                        if condensed_df.iloc[l,0].startswith(record.id):
                          if condensed_df.iloc[l,1].startswith('pri-rep'): # find prirep in df
                            prirep_location_start, prirep_location_end = get_locations("all.pdg.faa", condensed_df.iloc[l,0])
                            prirep_range = range(prirep_location_start,prirep_location_end)
                            prirep_prot_num = int((condensed_df.iloc[l,0]).split('_')[-1])

                            if range_subset(prirep_range, range_to_be_discovered) == True:
                              #get PICI location using 25kb to the left and right of pri-rep
                              pici_low_limit = prirep_location_start - 25000
                              pici_high_limit = prirep_location_end + 25000
                              # fix lower limit index to 0 if it is negative
                              if pici_low_limit < 0:
                                pici_low_limit = 0
                              # fix upper limit index to the length of the sequence if it is over
                              if pici_high_limit >= len(record.seq):
                                pici_high_limit = len(record.seq)

                              print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                              print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))
                              print('PICI start:', pici_low_limit)
                              print('PICI end:', pici_high_limit)
                              print('\nPICI sequence:')
                              print(record.seq[pici_low_limit:pici_high_limit])
                              print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                              print('\n')

                              # add PICI to list
                              name_list.append(record.id)
                              seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                              PICI_type.append("SaPI")
                              desc_list.append(record.description)
                              stop_loop = True


            ########################################################
            # THIRD search for co-localized pri-rep... determines G+
            ######################################################## 
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

                      # Get protein locations in host
                      # get the location of pri-rep and see if pri-rep is within 25 kb to left or right (50kb total)
                      int_location_start, int_location_end = get_locations("all.pdg.faa", condensed_df.iloc[i,0])
                      prirep_location_start, prirep_location_end = get_locations("all.pdg.faa", condensed_df.iloc[k,0])
                      print('\nBeginning G+ PICI determination...')
                      print('Int ({}) start: {}'.format(int_prot_num, int_location_start))
                      print('Int ({}) end: {}'.format(int_prot_num, int_location_end))
                      print('Pri-rep ({}) start: {}'.format(prirep_prot_num, prirep_location_start))
                      print('Pri-rep ({}) end: {}'.format(prirep_prot_num, prirep_location_end))

                      # set range for pri-rep to be found 
                      pici_low_limit = prirep_location_start - 25000
                      pici_high_limit = prirep_location_end + 25000
            

                      # fix lower limit index to 0 if it is negative
                      if pici_low_limit < 0:
                        pici_low_limit = 0


                      # fix upper limit index to the length of the sequence if it is over
                      if pici_high_limit >= len(record.seq):
                        pici_high_limit = len(record.seq)

                      print('PICI start:', pici_low_limit)
                      print('PICI end:', pici_high_limit)
                      print('\nPICI sequence:')
                      print(record.seq[pici_low_limit:pici_high_limit])
                      print('\nPICI length:', len(record.seq[pici_low_limit:pici_high_limit]))
                      print('\n')

                      # add PICI to list
                      name_list.append(record.id)
                      seq_list.append(record.seq[pici_low_limit:pici_high_limit])
                      PICI_type.append("G_pos_PICI")
                      desc_list.append(record.description)
                      stop_loop = True

                        
                        


            
        else:
          print('No hit')



# write PICI sequence to file
PICI_file = open("PICI_results", "w")

for m in range(len(seq_list)):
  PICI_file.write(">"  + str(desc_list[m]) +  ";" + str(PICI_type[m]) + ";" + str(m) + "\n" + str(seq_list[m]) + "\n")

PICI_file.close()
