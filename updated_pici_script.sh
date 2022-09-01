#! /bin/sh
set -x
# Set up your projects in the "data" directory
# For example: ./data/EXAMPLE_PROJECT_NAME
# Inside EXAMPLE_PROJECT_NAME is where your sequences, PICIs, and VirSorter2/BLAST results will be stored for that project
# In EXAMPLE_PROJECT_NAME set up a child directory "sequences" with all your data (./sequences/EXAMPLE_DATA.fasta)
# Make sure when you run this script you are in the parent directory of "sequences" (EXAMPLE_PROJECT_NAME)
# Make sure there is ONLY fasta files in your "sequences" directory
#test
db_nuc_path=./../../../databases/derived/BLAST_nucleotide_db.fna
db_prot_path=./../../pici_typer/databases/putative/BLAST_protein_db.faa
database=${d:=0} # 0 is the putative db, 1 is the derived db
integrase_identity=${g:=70}
alpa_identity=${a:=50}
input=${i:=.}
output=${o:=.}
dir_in=${y:=./sequences}
dir_out=${z:=.}


while getopts d:g:a:i:o:y:z: flag
do
    case "${flag}" in
        d) database=${OPTARG};;
        g) integrase_identity=${OPTARG};;
        a) alpa_identity=${OPTARG};;
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        y) dir_in=${OPTARG};;
        z) dir_out=${OPTARG};;
    esac
done


mkdir -p "${dir_in}/results"
mkdir -p "${dir_in}/tmp_PICIs"

# Wrapper script that iterates over every sequence in directory "sequences"
for f in ${dir_in}/${input}/*
do
    echo "${f}"
    # Checks if file has already been processed
	if test -e "${dir_in}/results/${f##*/}"; then
    		echo "${f##*/} has already been processed. The sequence will be skipped. If you wish you re-run this sequence remove the corresponding directory from \"result\"."
		continue
	else
		# Temporary folder to create a copy of sequence named "all.fna" so that it may be fed into pici_integrase_trimmer_script.py
        	echo "Processing ${f##*/}..."
        	TEMP=`mktemp -d`
        	cp $dir_in/${f##*/} $TEMP
#        	cd $TEMP
        	mv "${TEMP}/${f##*/}" "${TEMP}/all.fna"
        	TEMP_FASTA="${TEMP}/all.fna"
	
        	#  Run BLAST (note: the output name is the same regardless of BLAST operation)
                if [ "$database" -eq "0" ]; then
                        echo "Performing tBLASTn on ${f##*/}..."
                	tblastn -query ${db_prot_path} -subject ${TEMP_FASTA} -task tblastn -evalue 0.001 -outfmt 6 -out "${TEMP}/tBLASTn_results.out"
                elif [ "$database" -eq "1" ]; then
        		echo "Performing BLASTn on ${f##*/}..."
        		blastn -query ${db_nuc_path} -subject ${TEMP_FASTA} -task blastn -evalue 0.001 -outfmt 6 -out "${TEMP}/tBLASTn_results.out"
        	fi

        	# Run integrase trimmer (calls the script in the data
        	echo "Beginning Trim on ${f##*/}..."
        	python ./../../../scripts/pici_integrase_trimmer_script.py --i $integrase_identity
        	
        	# If no integrases >= 90% identity then stop the current iteration
        	if ! [ -s trimmed_file ]; then
        	        cd ..
        	        rm -r tmp
        	        echo "Terminating ${f##*/}: no match for integrase..."
        	        continue
        	        
        	# else continue with script
        	else
        	        # convert trimmed file back to original filename
        	        mv trimmed_file ${f##*/}
        	        cd ..

			# Run Prodigal
			echo "Running Prodigal on ${f##*/}..."
			mkdir ./tmp_PICIs/${f##*/}
			mkdir ./results/${f##*/}
			prodigal -i ./tmp/${f##*/} -a ./results/${f##*/}/all.pdg.faa -f gff -o ./results/${f##*/}/all.pdg.gff -p meta

			cp ./tmp/${f##*/} ./results/${f##*/} # this will preserve host info
			rm -r tmp
        	        
        	        # Run Blastp 
        	        echo "Performing BLASTp on ${f##*/}..."
        	        cd results/${f##*/}/
                        if [ "$database" -eq "0" ]; then
        	        	blastp -query ./all.pdg.faa -db ./../../../../databases/putative/PICI_BLAST_DB -task blastp -evalue 0.001 -outfmt 6 -out BLASTp_results.out
        	        elif [ "$database" -eq "1" ]; then
                                tblastn -query ./all.pdg.faa -db ./../../../../databases/derived/PICI_BLAST_DB -task tblastn -evalue 0.001 -outfmt 6 -out BLASTp_results.out
                        fi
        	        # Set up for PICI typer
        	        echo "Running PICI-typer script on ${f##*/}..."
        	        mkdir python
			cp ${f##*/} ./python # Preserves host info
        	        cp all.pdg.faa ./python
        	        cp BLASTp_results.out ./python
        	        cd python
			mv ${f##*/} all.fna
        	        
        	        
        	        # Run PICI typer script 
        	        python3 ./../../../../../scripts/prototype_typer.py --i $integrase_identity --a $alpa_identity
			
			# Remove duplicates
			python3 ./../../../../../scripts/duplicate_remover.py
        	        
        	        # Move PICI results to PICI directory
        	        mv ./PICI_results ./../../../tmp_PICIs/${f##*/}
        	        cd ../../../
        	        echo "Finished ${f##*/}..."
        	fi
	fi
done

echo "Compiling all PICIs into one file 'ALL_PICIs.fasta'..."
./../../scripts/pici_collector.sh

echo "Separating PICI type and phage satellites..."
python3 ./../../scripts/pici_separator.py

echo "Reviewing phage satellites..."
if [ "$database" -eq "0" ]; then
	blastx -query ./Phage_Satellites.fasta -subject ${db_nuc_path} -task blastx -evalue 0.001 -outfmt 6 -out BLAST_results.out
elif [ "$database" -eq "1" ]; then
	blastn -query ./Phage_Satellites.fasta -subject ./../../databases/derived/BLAST_nucleotide_db.fna -task blastn -evalue 0.001 -outfmt 6 -out BLAST_results.out
fi
python3 ./../../scripts/phage_satellite_review.py --a $alpa_identity
sed -i -- 's/phage_satellite/G_neg_PICI/g' ./G_neg_PICI_reviewed
cat G_neg_PICI_reviewed G_neg_PICIs.fasta SaPIs.fasta phage_satellite_reviewed > new_ALL_PICIs.fasta
mv new_ALL_PICIs.fasta ALL_PICIs.fasta
cat G_neg_PICI_reviewed G_neg_PICIs.fasta > new_G_neg_PICIs.fasta
mv new_G_neg_PICIs.fasta G_neg_PICIs.fasta
mv phage_satellite_reviewed ./Phage_Satellites.fasta
rm G_neg_PICI_reviewed


echo "Collecting host genomes..."
./../../scripts/genome_collector.sh

echo "Creating table..."
cat $(basename "$PWD")_ALL_PICIs.fasta | grep -e "^>" | sed 's/>//g' | sed 's/;/\t/g' > $(basename "$PWD")_PICI_table.tsv

echo "Done."



if [ "$output" == "." ]; then
        continue
else
        mkdir $dir_out/$output/
        mv $(basename "$PWD")_ALL_PICIs.fasta $dir_out/$output/${output}_ALL_PICIs.fasta
        mv $(basename "$PWD")_G_neg_PICIs.fasta $dir_out/$output/${output}_G_neg_PICIs.fasta
        #mv $(basename "$PWD")_SaPIs.fasta $dir_out/$output/${output}_SaPIs.fasta
        mv $(basename "$PWD")_Phage_Satellites.fasta $dir_out/$output/${output}_Phage_Satellites.fasta
        mv $(basename "$PWD")_PICI_table.tsv $dir_out/$output/${output}_PICI_table.tsv
        mv $(basename "$PWD")_ALL_PICIs_host_genomes.fasta $dir_out/$output/${output}_ALL_PICIs_host_genomes.fasta
        mv $(basename "$PWD")_G_neg_PICIs_host_genomes.fasta $dir_out/$output/${output}_G_neg_PICIs_host_genomes.fasta
        mv $(basename "$PWD")_Phage_Satellites_host_genomes.fasta $dir_out/$output/${output}_Phage_Satellites_host_genomes.fasta
        mv $(basename "$PWD")_ALL_PICIs_host_names.txt $dir_out/$output/${output}_ALL_PICIs_host_names.txt
        mv $(basename "$PWD")_G_neg_PICIs_host_names.txt $dir_out/$output/${output}_G_neg_PICIs_host_names.txt
        mv $(basename "$PWD")_Phage_Satellites_host_names.txt $dir_out/$output/${output}_Phage_Satellites_host_names.txt

fi

