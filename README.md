Installation:

1) Download and install VirSorter 2
https://github.com/jiarong/VirSorter2/blob/master/README.md

2) Download and install BLAST+ (if error occurs make sure it is the latest version from https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz 
tar zxvpf ncbi-blast-2.12.0+-x64-linux.tar.gz

3) Download setup.sh script
wget https://raw.githubusercontent.com/klh272/pici-repository/main/setup.sh
chmod +x setup.sh




You now are all setup to run the PICI typer! Some errors may occur with dependencies not installed (often related to VirSorter2), so a test run is recommended. Simply install these dependencies with pip if this is the case.

In your "data" directory you will create your projects (e.g., ./data/EXAMPLE_PROJECT_NAME/).
These project directories are where you PICI, BLAST, and VirSorter2 results will be stored for that particular project.

Within your project directory (EXAMPLE_PROJECT_NAME) create a sub-directory called "sequences" where you will store all your data relevant to that project (e.g., ./data/EXAMPLE_PROJECT_NAME/sequences/EXAMPLE_DATA_NAME.fasta)

You will run the pici script ./../../scripts/updated_pici_script.sh from EXAMPLE_PROJECT_NAME.

The output should create two additional directories: "PICIs" and "results"
PICIs: Where identified PICIs are stored
results: Where VirSorter2 and BLAST results are stored



