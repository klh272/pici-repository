#! /bin/sh
set -x
# Set up your projects in the "data" directory
# For example: ./data/EXAMPLE_PROJECT_NAME
# Inside EXAMPLE_PROJECT_NAME is where your sequences, PICIs, and VirSorter2/BLAST results will be stored for that project
# In EXAMPLE_PROJECT_NAME set up a child directory "sequences" with all your data (./sequences/EXAMPLE_DATA.fasta)
# Make sure when you run this script you are in the parent directory of "sequences" (EXAMPLE_PROJECT_NAME)
# Make sure there is ONLY fasta files in your "sequences" directory


scripts_path=$(dirname $(readlink -f $0))
db_nuc_path=./../../../databases/derived/BLAST_nucleotide_db.fna
db_prot_path=./../../pici_typer/databases/putative/BLAST_protein_db.faa
db_blastp_path=./../../pici_typer/databases/putative/PICI_BLAST_DB
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
sequences_dir=${dir_in}/${input}/sequences
# Wrapper script that iterates over every sequence in directory "sequences"
for f in ${sequences_dir}/*
do
    echo "${f}"
    # Checks if file has already been processed
    if test -e "${dir_in}/results/${f##*/}"; then
    	echo "${f##*/} has already been processed. The sequence will be skipped. If you wish you re-run this sequence remove the corresponding directory from \"result\"."
	continue
	
	else
	    # Temporary folder to create a copy of sequence named "all.fna" so that it may be fed into pici_integrase_trimmer_script.py
	    TEMP_NAME=${f##*/}
            echo "Processing ${TEMP_NAME}..."
            TEMP=`mktemp -d`
	    echo `pwd`
	    ls ${f}
        	cp ${f} "$TEMP/all.fna"
        	TEMP_FASTA="${TEMP}/all.fna"
		echo ${TEMP} #Delete this
        	#  Run BLAST (note: the output name is the same regardless of BLAST operation)
                if [ "$database" -eq "0" ]; then
                        echo "Performing tBLASTn on ${TEMP_NAME}..."
                	tblastn -query ${db_prot_path} -subject ${TEMP_FASTA} -task tblastn -evalue 0.001 -outfmt 6 -out "${TEMP}/tBLASTn_results.out"
                elif [ "$database" -eq "1" ]; then
        		echo "Performing BLASTn on ${TEMP_NAME}..."
        		blastn -query ${db_nuc_path} -subject ${TEMP_FASTA} -task blastn -evalue 0.001 -outfmt 6 -out "${TEMP}/tBLASTn_results.out"
        	fi

        	# Run integrase trimmer (calls the script in the data
        	echo "Beginning Trim on ${TEMP_NAME}..."
        	python ${scripts_path}/pici_integrase_trimmer_script.py --i $integrase_identity --blast "${TEMP}/tBLASTn_results.out" --fasta "${TEMP}/all.fna" --output "${TEMP}/trimmed_file"
        	
        	# If no integrases >= 90% identity then stop the current iteration
        	if ! [ -s "${TEMP}/trimmed_file" ]; then

        	        #rm -r ${TEMP} 
        	        echo "Terminating ${TEMP_NAME}: no match for integrase..."
        	        continue

        	# else continue with script
        	else
        	       # convert trimmed file back to original filename
        	        cp "${TEMP}/trimmed_file" ${TEMP}/${TEMP_NAME} #WTF
        	     #   cd ..

			# Run Prodigal
			echo "Running Prodigal on ${TEMP_NAME}..."
			mkdir -p ${dir_in}/tmp_PICIs/${TEMP_NAME}
			mkdir -p ${dir_in}/results/${TEMP_NAME}
			prodigal -i ${TEMP}/${TEMP_NAME} -a ${dir_in}/results/${TEMP_NAME}/all.pdg.faa -f gff -o ${dir_in}/results/${TEMP_NAME}/all.pdg.gff -p meta
			AAFILE=${dir_in}/results/${TEMP_NAME}/all.pdg.faa
			
			cp ${TEMP}/${TEMP_NAME} ${dir_in}/results/${TEMP_NAME} # this will preserve host info
			#rm -r tmp
        	        
        	        # Run Blastp 
        	        echo "Performing BLASTp on ${f##*/}..."
        	      #  cd results/${f##*/}/
                        if [ "$database" -eq "0" ]; then
        	        	blastp -query ${dir_in}/results/${TEMP_NAME}/all.pdg.faa -db ${db_blastp_path} -task blastp -evalue 0.001 -outfmt 6 -out ${TEMP}/BLASTp_results.out
        	        elif [ "$database" -eq "1" ]; then
                                tblastn -query ${dir_in}/results/${TEMP_NAME}/all.pdg.faa -db ${db_nuc_path} -task tblastn -evalue 0.001 -outfmt 6 -out ${TEMP}/BLASTp_results.out
                        fi
        	        # Set up for PICI typer
        	        echo "Running PICI-typer script on ${f##*/}..."
        	        mkdir ${TEMP}/python
			cp ${sequences_dir}/${TEMP_NAME} ${TEMP}/python # Preserves host info
        	        cp ${AAFILE} ${TEMP}/python
        	        cp ${TEMP}/BLASTp_results.out ${TEMP}/python
        	       # cd python
			#mv ${f##*/} all.fna #WTF, think it should be covered above
        	        
        	        
        	        # Run PICI typer script 
        	        python ${scripts_path}/prototype_typer.py --i $integrase_identity --a $alpa_identity --blast "${TEMP}/tBLASTn_results.out" --fasta "${TEMP}/all.fna" --output "${TEMP}/PICI_results" --aa "${TEMP}/python/all.pdg.faa" --temp ${TEMP}
			
			# Remove duplicates
			python ${scripts_path}/duplicate_remover.py --input "${TEMP}/PICI_results" --output "${TEMP}/PICI_results"
        	        
        	        # Move PICI results to PICI directory
        	        mv "${TEMP}/PICI_results" "${dir_in}/tmp_PICIs/"
        	        #cd ../../../
        	        echo "Finished ${TEMP_NAME}..."
        	fi
	fi
done

exit 1
#Exit the typer if no PICIs were found
if [ ! "$(ls -A ${dir_in}/results)" ]; then
    echo "${dir_in}/results is  empty"
    exit 0
fi


echo "Compiling all PICIs into one file 'ALL_PICIs.fasta'..."
#${scripts_path}/scripts/pici_collector.sh
touch "${TEMP}/ALL_PICIs.fasta"

for f in ./PICIs/*
do
        cat ./PICIs/${f##*/}/PICI_results >> ALL_PICIs.fasta

done


echo "Separating PICI type and phage satellites..."
python ${scripts_path}/scripts/pici_separator.py

echo "Reviewing phage satellites..."
if [ "$database" -eq "0" ]; then
	blastx -query ./Phage_Satellites.fasta -subject ${db_prot_path} -task blastx -evalue 0.001 -outfmt 6 -out BLAST_results.out
elif [ "$database" -eq "1" ]; then
	blastn -query ./Phage_Satellites.fasta -subject ${db_nuc_path} -task blastn -evalue 0.001 -outfmt 6 -out BLAST_results.out
fi
python ${scripts_path}/scripts/phage_satellite_review.py --a $alpa_identity
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

