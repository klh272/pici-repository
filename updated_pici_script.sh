#! /bin/sh

# Set up your projects in the "data" directory
# For example: ./data/EXAMPLE_PROJECT_NAME
# Inside EXAMPLE_PROJECT_NAME is where your sequences, PICIs, and VirSorter2/BLAST results will be stored for that project
# In EXAMPLE_PROJECT_NAME set up a child directory "sequences" with all your data (./sequences/EXAMPLE_DATA.fasta)
# Make sure when you run this script you are in the parent directory of "sequences" (EXAMPLE_PROJECT_NAME)
# Make sure there is ONLY fasta files in your "sequences" directory

mkdir results
mkdir PICIs

if test -e "./BLAST_DB.tsv"; then

	touch BLAST_DB.tsv
	echo "Creating BLAST DB file."
fi

# Wrapper script that iterates over every sequence in directory "sequences"
for f in ./sequences/*
do
	# Checks if file has already been processed
	if test -e "./results/${f##*/}"; then
    		echo "${f##*/} has already been processed. The sequence will be skipped. If you wish you re-run this sequence remove the corresponding directory from \"result\"."
		continue
	else
		# Temporary folder to create a copy of sequence named "all.fna" so that it may be fed into pici_integrase_trimmer_script.py
        	echo "Processing ${f##*/}..."
        	mkdir tmp
        	cp ./sequences/${f##*/} ./tmp
        	cd tmp
        	mv ${f##*/} all.fna
        	
        	#  Run tBLASTn 
        	echo "Performing tBLASTn on ${f##*/}..."
        	tblastn -query ./../../../databases/BLAST_protein_db.faa -subject ./all.fna -task tblastn -evalue 0.001 -outfmt 6 -out tBLASTn_results.out
        	
        	# Run integrase trimmer (calls the script in the data
        	echo "Beginning Trim on ${f##*/}..."
        	python ./../../../scripts/pici_integrase_trimmer_script.py
        	
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
			mkdir ./PICIs/${f##*/}
			mkdir ./results/${f##*/}
			prodigal -i ./tmp/${f##*/} -a ./results/${f##*/}/all.pdg.faa -f gff -o ./results/${f##*/}/all.pdg.gff -p meta

			# Run HattCI
			#echo "Running HattCI on ${f##*/}..."
			#../../HattCI/hattci.out -b -t 4 ./tmp/${f##*/} ./results/${f##*/}/hattci.out

        	        # Run VirSorter2
        	        #echo "Running VirSorter2 on ${f##*/}..."
        	        #mkdir ./PICIs/${f##*/}
        	        #virsorter run -w ./results/${f##*/} -i ./tmp/${f##*/} -j 20 all
			cp ./tmp/${f##*/} ./results/${f##*/} # this will preserve host info after running VirSorter2
			rm -r tmp
        	        
        	        # Run Blastp 
        	        echo "Performing BLASTp on ${f##*/}..."
        	        cd results/${f##*/}/
        	        blastp -query ./all.pdg.faa -db ./../../../../databases/PICI_BLAST_DB -task blastp -evalue 0.001 -outfmt 6 -out BLASTp_results.out
        	        
        	        # Set up for PICI typer
        	        echo "Running PICI-typer script on ${f##*/}..."
        	        mkdir python
			cp ${f##*/} ./python # Preserves host info
        	        cp all.pdg.faa ./python
        	        cp BLASTp_results.out ./python
			csplit ./hattci.out '/^--------------------------------------------$/' '{*}' #parses output to obtain table
			cp xx00 ./python #parsed table from HattCI
        	        cd python
			mv ${f##*/} all.fna
        	        
        	        
        	        # Run PICI typer script 
        	        python3 ./../../../../../scripts/prototype_typer.py
			
			# Remove duplicates
			python3 ./../../../../../scripts/duplicate_remover.py
        	        
        	        # Move PICI results to PICI directory
        	        mv ./PICI_results ./../../../PICIs/${f##*/}
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
blastx -query ./Phage_Satellites.fasta -subject ./../../databases/BLAST_protein_db.faa -task blastx -evalue 0.001 -outfmt 6 -out BLAST_results.out
python3 ./../../scripts/phage_satellite_review.py
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
