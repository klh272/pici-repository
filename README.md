Installation:

1) Download and install Prodigal <br>
https://github.com/hyattpd/prodigal/wiki/Introduction <br>
<br>
2) Download and install BLAST+ (if error occurs make sure it is the latest version from https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)<br>
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz <br>
tar zxvpf ncbi-blast-2.12.0+-x64-linux.tar.gz <br>
<br>
3) Download setup.sh script <br>
wget https://raw.githubusercontent.com/klh272/pici-repository/main/setup.sh <br>
chmod +x setup.sh <br>
./setup.sh<br>
<br>




You now are all setup to run the PICI typer! Some errors may occur with dependencies not installed, so a test run is recommended. Simply install these dependencies with pip if this is the case.<br>
<br>
In your "data" directory you will create your projects (e.g., ./data/EXAMPLE_PROJECT_NAME/). <br>
These project directories are where you PICI, BLAST, and Prodigal results will be stored for that particular project. <br>
<br>
Within your project directory (EXAMPLE_PROJECT_NAME) create a sub-directory called "sequences" where you will store all your data relevant to that project (e.g., ./data/EXAMPLE_PROJECT_NAME/sequences/EXAMPLE_DATA_NAME.fasta) <br>
<br>
You will run the pici script ./../../scripts/updated_pici_script.sh from your project directory (EXAMPLE_PROJECT_NAME). <br>
<br>
The output should create two additional directories: "PICIs" and "results" <br>
PICIs: Where identified PICIs are stored <br>
results: Where Prodigal and BLAST results are stored <br>



