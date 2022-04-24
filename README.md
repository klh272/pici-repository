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
PICIs: Where pre-reviewed PICIs are stored <br>
results: Where Prodigal and BLAST results are stored <br>
<br>
<br>
In your project directory will be the final results, with the project's name as the prefix to the files.


<br>
<br>
<br>
Table Output:<br>
1. <b>sequenceID</b>: The record's name and description<br>
2. <b>type</b>: G_neg_PICI, SaPI, phage_satellite<br>
3. <b>length</b>: PICI length in bp<br>
4. <b>start</b>: The beginning coordinate of the PICI in the host's genome<br>
5. <b>end</b>:  The end coordinate of the PICI in the host's genome<br>
6. <b>orientation</b>: + (forward) or - (backward)<br>
7. <b>int_ID</b>: The corresponding int unique ID in the BLAST database<br>
8. <b>int_identity</b>: The identity percentage of the int hit<br>
9. <b>int_start_coord</b>: The start coordinate of the int<br>
10. <b>int_end coord</b>: The end coordinate of the int<br>
11. <b>alpA_ID</b>: The corresponding alpA unique ID in the BLAST database<br>
12. <b>alpA_identity</b>: The identity percentage of the alpA hit<br>
13. <b>alpA_start_coord</b>: The start coordinate of the alpA<br>
14. <b>alpA_end coord</b>: The end coordinate of the alpA<br>
15. <b>sis_ID</b>: The corresponding sis unique ID in the BLAST database<br>
16. <b>sis_identity</b>: The identity percentage of the sis hit<br>
17. <b>sis_start_coord</b>: The start coordinate of the sis<br>
18. <b>sis_end coord</b>: The end coordinate of the sis<br>
19. <b>prirep_ID</b>: The corresponding pri-rep unique ID in the BLAST database<br>
20. <b>prirep_identity</b>: The identity percentage of the pri-rep hit<br>
21. <b>prirep_start_coord</b>: The start coordinate of the pri-rep<br>
22. <b>prirep_end_coord</b>: The end coordinate of the pri-rep<br>
23. <b>attL_start</b>: The start coordinate of the attL<br>
24. <b>attL_end</b>: The end coordinate of the attL<br>
25. <b>attR_start</b>: The start coordinate of the attR<br>
26. <b>attR_end</b>: The end coordinate of the attR<br>
27. <b>attC</b>: the attC sequence that was predicted. This can be one sequence if attR == attL, or two sequences separated by a "/" if they are not. If attL or attR or both are not found, a "NA" will appear.<br>
28. <b>trim_quality</b>: "high" indicates that attL and attR have been predicted, "medium" indicates only one attachment site has been predicted, and "low" means none have been predicted. The trimming methods for each of the qualities differs. That is, "high" quality means the PICI coordinates are exact or close to exact and only include the PICI itself; "medium" quality means the PICI coordinates are semi-exact and use an arbitrary border up/downstream of the failed attachment site prediction with consideration of the PICI's genetic architecture; and "low" quality means and arbitrary region has been extracted around the PICI with consideration of it's genetic architecture.<br>
29. <b>terS</b>: Indicates the presence of a terS in the PICI. 0 = not present, 1 = present<br>
30. <b>ppi</b>: Indicates the presence of a ppi in the PICI. 0 = not present, 1 = present<br>
31. <b>rpp</b>: Indicates the presence of a rpp in the PICI. 0 = not present, 1 = present<br>
32. <b>PICI_#</b>: The order in which the PICI was found by the typer. "0" means it was the first PICI found in the sequence, "1" means it was the second PICI found in the sequence, etc. This does not necessarily tell you the number of PICIs in the genome.<br>
