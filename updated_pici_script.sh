#! /bin/sh

# Assumes you are in a parent directory that contains a child directory labeled 'sequences' that contains all data (genomes)!
mkdir results
mkdir PICIs

# Wrapper script that iterates over every sequence in directory "sequences"
for f in ./sequences/*
do
        echo "Processing ${f##*/}..."
        # Temporary folder to create a copy of sequence named "all.fna" so that it may be fed into pici_integrase_trimmer_script.py
        mkdir tmp
        cp ./sequences/${f##*/} ./tmp
        cd tmp
        mv ${f##*/} all.fna
        
        #  Run tBLASTn (*** add you path to BLAST_protein_db.faa ***)
        echo "Performing tBLASTn on ${f##*/}..."
        tblastn -query /PATH/TO/BLAST_protein_db.faa -subject ./all.fna -task tblastn -evalue 0.001 -outfmt 6 -out tBLASTn_results.out
        
        # Run integrase trimmer (*** add your path to pici_integrase_trimmer_script.py ***)
        echo "Beginning Trim on ${f##*/}..."
        python /PATH/TO/pici_integrase_trimmer_script.py
        
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
                
                # Run VirSorter2
                echo "Running VirSorter2 on ${f##*/}..."
                mkdir ./PICIs/${f##*/}
                virsorter run -w ./results/${f##*/} -i ./tmp/${f##*/} -j 20 all
		            cp ./tmp/${f##*/} ./results/${f##*/}/iter-0 # this will preserve host info after running VirSorter2
		            rm -r tmp
                
                # Run Blastp (*** add your path to BLAST_protein_db ***)
                echo "Performing BLASTp on ${f##*/}..."
                cd results/${f##*/}/iter-0
                blastp -query ./all.pdg.faa -db /PATH/TO/BLAST_protein_db -task blastp -evalue 0.001 -outfmt 6 -out BLASTp_results.out
                
                # Set up for PICI typer
                echo '\033[0;35m'"Running PICI-typer script on ${f##*/}..."'\033[0m'
                mkdir python
		            cp ${f##*/} ./python # Preserves host info
                cp all.pdg.faa ./python
                cp BLASTp_results.out ./python
                cd python
		            mv ${f##*/} all.fna
                
                
                # Run PICI typer script (*** add your path to updated_official_pici_typer.py ***)
                python3 /PATH/TO/updated_official_pici_typer.py
                
                # Remove temporary files and move PICI results to PICI directory
                rm all.fna all.pdg.faa BLASTp_results.out
                mv ./* ./../../../../PICIs/${f##*/}
                cd ..
                rm -r python
                cd ../../../
                echo '\033[0;35m'"Finished ${f##*/}..."'\033[0m'
        fi
done
