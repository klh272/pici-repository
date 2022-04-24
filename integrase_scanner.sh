#! /bin/sh

# Set up your projects in the "data" directory
# For example: ./data/EXAMPLE_PROJECT_NAME
# Inside EXAMPLE_PROJECT_NAME is where your sequences, PICIs, and VirSorter2/BLAST results will be stored for that project
# In EXAMPLE_PROJECT_NAME set up a child directory "sequences" with all your data (./sequences/EXAMPLE_DATA.fasta)
# Make sure when you run this script you are in the parent directory of "sequences" (EXAMPLE_PROJECT_NAME)
# Make sure there is ONLY fasta files in your "sequences" directory

mkdir results
mkdir PICIs
touch list_of_potential_hosts

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
                mv ./list_of_potential_hosts ./tmp/
                cd tmp
                mv ${f##*/} all.fna
                
                #  Run tBLASTn 
                echo "Performing tBLASTn on ${f##*/}..."
                tblastn -query ./../../../databases/BLAST_protein_db.faa -subject ./all.fna -task tblastn -evalue 0.001 -outfmt 6 -out tBLASTn_results.out

                # Run search
                echo "Performing potential host search on ${f##*/}..."
                python3 ../../../scripts/integrase_scanner.py
                mv ./list_of_potential_hosts ./..
                cd ..
                rm -r tmp
                
	fi
done
