#! /bin/sh

# Make master directory
mkdir pici_typer

# Set up directories
cd pici_typer
mkdir data
mkdir databases
mkdir databases/putative
mkdir databases/derived
mkdir scripts


# set up databases
cd databases/putative
wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_protein_db_new_integrases_21032022.fasta -O BLAST_protein_db.faa
makeblastdb -in ./BLAST_protein_db.faa -input_type fasta -dbtype prot -out ./PICI_BLAST_DB
cd ../derived
wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_nucleotide_db.fna -O BLAST_nucleotide_db.fna
makeblastdb -in ./BLAST_nucleotide_db.fna -input_type fasta -dbtype nucl -out ./PICI_BLAST_DB

cd ../..

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

# rm setup.sh
