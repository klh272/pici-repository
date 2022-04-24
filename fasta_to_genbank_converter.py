#from Bio import SeqIO

#count = SeqIO.convert("all.fna", "fasta", "all.gb", "genbank", alphabet = IUPACUnambiguousDNA())
#print("Converted %i records" % count)


# New style
from Bio import SeqIO

# This file has a single record only
fasta_file = 'all.fna'
with open(fasta_file, mode='r') as handle:
  
  for record in SeqIO.parse(handle, 'fasta'):

    #record = SeqIO.read("all.fna", "fasta")
    record.annotations["molecule_type"] = "DNA"
SeqIO.write(record, "all.gb", "genbank")
