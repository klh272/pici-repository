#! /bin/sh

# Assumes you are in a parent directory that contains a child directory labeled 'sequences' that contains all data (genomes)!
mkdir results
mkdir PICIs

# Wrapper script that iterates over every sequence in directory "sequences"
for f in ./sequences/*
do
	echo '\033[0;35m'"Processing ${f##*/}..."'\033[0m'
	# Temporary folder to create a copy of sequence named "all.fna" so that it may be fed into pici_integrase_trimmer_script.py
	mkdir tmp
	cp ./sequences/${f##*/} ./tmp
	cd tmp
	mv ${f##*/} all.fna
	echo '\033[0;35m'"Performing tBLASTn on ${f##*/}..."'\033[0m'
	# Run tBLASTn (primarily interested in integrases for trimming)
	tblastn -query ../BLAST_protein_db.faa -subject ./all.fna -task tblastn -evalue 0.001 -outfmt 6 -out tBLASTn_results.out
	echo '\033[0;35m'"Beginning Trim on ${f##*/}..."'\033[0m'
	# Run integrase trimmer
	python ../pici_integrase_trimmer_script.py
	# If no integrases >= 90% identity then stop the current iteration
	if ! [ -s trimmed_file ]; then
		cd ..
		rm -r tmp
		echo '\033[0;31m'"Terminating ${f##*/}: no match for integrase..."'\033[0m'
		continue
	else
		mv trimmed_file ${f##*/}
		cd ..
		echo '\033[0;35m'"Running VirSorter2 on ${f##*/}..."'\033[0m'
		mkdir ./PICIs/${f##*/}
		# Run VirSorter2
		virsorter run -w ./results/${f##*/} -i ./tmp/${f##*/} -j 64 all
		rm -r tmp
		echo '\033[0;35m'"Performing BLASTp on ${f##*/}..."'\033[0m'
		cd results/${f##*/}/iter-0
		# Run BLASTp
		blastp -query ./all.pdg.faa -db ../../../PICI_BLAST_DB -task blastp -evalue 0.001 -outfmt 6 -out BLAST_results.out
		echo '\033[0;35m'"Running PICI-typer script on ${f##*/}..."'\033[0m'
		mkdir python
		cp all.fna ./python
		cp all.pdg.faa ./python
		cp BLAST_results.out ./python
		cd python
		# Run PICI typer script
		python3 ../../../../official_pici_typer.py
		rm all.fna all.pdg.faa BLAST_results.out
		mv ./* ./../../../../PICIs/${f##*/}
		cd ..
		rm -r python
		cd ../../../
		echo '\033[0;35m'"Finished ${f##*/}..."'\033[0m'
	fi
done
