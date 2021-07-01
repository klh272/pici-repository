#! /bin/sh

# Assumes you are in a parent directory that contains a child directory labeled 'sequences' that contains all data (genomes)
mkdir results
mkdir PICIs


for f in ./sequences/*
do
        echo "Processing ${f##*/}..."
        virsorter run -w ./results/${f##*/}.results -i ./sequences/${f##*/} -j 64 all
        
	echo "Performing BLASTp on ${f##*/}..."
        cd results/${f##*/}.results/iter-0
        blastp -query ./all.pdg.faa -db /mibi/users/mwiddowson/BLAST_db/PICI_BLAST_DB -task blastp -evalue 0.001 -outfmt 6 -out BLAST_results.out 

        echo "Running PICI-typer script on ${f##*/}..."
        mkdir python
        cp all.fna ./python
        cp all.pdg.faa ./python
        cp BLAST_results.out ./python
        cd python
        python3 /mibi/users/mwiddowson/pici_typer_script.py
        rm all.fna all.pdg.faa BLAST_results.out
        mv ./* ./../../../../PICIs
        cd ..
        rm -r python
        cd ../../../

        echo "Finished ${f##*/}..."
done
