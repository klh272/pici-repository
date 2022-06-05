#! /bin/sh

# Make master directory
mkdir pici_typer

# Set up directories
cd pici_typer
mkdir data
mkdir databases
mkdir scripts


# set up databases
wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_protein_db_new_integrases_21032022.fasta -P ./databases/putative/ -O ./databases/putative/BLAST_protein_db.faa
makeblastdb -in ./databases/putative/BLAST_protein_db.faa -input_type fasta -dbtype prot -out ./databases/putative/PICI_BLAST_DB
wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_nucleotide_db.fna.fasta -P ./databases/derived/ -O ./databases/putative/BLAST_nucleotide_db.fna
makeblastdb -in ./databases/derived/BLAST_protein_db.fna -input_type fasta -dbtype nucl -out ./databases/derived/PICI_BLAST_DB


# download scripts
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_integrase_trimmer_script.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/prototype_typer.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/phage_satellite_review.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/duplicate_remover.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/updated_pici_script.sh -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_collector.sh -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_separator.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/genome_collector.sh -P ./scripts/

chmod +x ./scripts/updated_pici_script.sh
chmod +x ./scripts/pici_collector.sh
chmod +x ./scripts/genome_collector.sh

rm setup.sh
