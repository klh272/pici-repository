#! /bin/sh

# Make master directory
mkdir pici_typer

# Set up directories
cd pici_typer
mkdir data
mkdir databases
mkdir scripts


# set up databases
wget https://raw.githubusercontent.com/klh272/pici-repository/main/BLAST_protein_db.faa -P ./databases/
makeblastdb -in ./databases/BLAST_protein_db.faa -input_type fasta -dbtype prot -out ./databases/PICI_BLAST_DB

# download scripts
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_integrase_trimmer_script.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/improved_pici_typer_script.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/phage_satellite_review.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/duplicate_remover.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/updated_pici_script.sh -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_collector.sh -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/pici_separator.py -P ./scripts/
wget https://raw.githubusercontent.com/klh272/pici-repository/main/genome_collector.sh 

chmod +x ./scripts/updated_pici_script.sh
chmod +x ./scripts/pici_collector.sh
chmod +x ./scripts/genome_collector.sh
