Download:
  1) pici_scipt.sh: wrapper script that iterates protocol over every sequence in the directory "sequences". You may need to change the paths to fit your set up.
  2) official_pici_typer.py: main algorithm to determine SaPI, G+ and G- PICIs. 
  3) pici_integrase_trimmer_script.py: trims original sequences +/-35kb around >= 90% identity integrase matches. Reduces computational time and complexity in VirSorter2
  4) BLAST_protein_db.faa: BLAST database of putative PICI-specific proteins (amino acids). This is necessary for tBLASTn and BLASTp procedures.


IMPORTANT:
Before doing anything create a child directory named "sequences". Here you will store your sequences of interest. The script should be ran in the parent directory 
directly above "sequences". Make sure there are only fasta files in the "sequences" directory.

NOTE: 
At the moment there may be direct overlaps in PICIs for a given sequence... Currently working on removing complete duplicates in the output.

The protocol will create two new folders in the parent directory (above child directory "sequences") called "results" and "PICIs". The "results" directory contains
the output from VirSorter2, as well as the results from BLASTp. There isn't much of interest hear regarding PICIs unless the user wishes to explore them. The "PICIs" 
directory holds the identified PICIs from the corresponding sequences. They are in fasta file format and contain information regarding what type of PICI they are in
their header. 
