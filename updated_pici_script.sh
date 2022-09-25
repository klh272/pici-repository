#! /bin/sh

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
        		blastn -query ${db_nuc_path} -subject ${TEMP_FASTA} -task blastn -evalue 0.001 -outfmt 6 -out "${TEMP}/tBLASTn_results.out" #TODO: is this output name correct?
        	fi

        	# Run integrase trimmer (calls the script in the data
        	echo "Beginning Trim on ${TEMP_NAME}..."
        	python ${scripts_path}/pici_integrase_trimmer_script.py --i $integrase_identity --blast "${TEMP}/tBLASTn_results.out" --fasta "${TEMP}/all.fna" --output "${TEMP}/trimmed_file"
        	
        	# If no integrases >= 90% identity then stop the current iteration
        	if ! [ -s "${TEMP}/trimmed_file" ]; then
        	        echo "Terminating ${TEMP_NAME}: no match for integrase..."
        	        continue
        	# else continue with script
        	else
        	       # convert trimmed file back to original filename
        	        cp "${TEMP}/trimmed_file" ${TEMP}/${TEMP_NAME} #WTF
			# Run Prodigal
			echo "Running Prodigal on ${TEMP_NAME}..."
			mkdir -p ${dir_in}/tmp_PICIs/${TEMP_NAME}
			mkdir -p ${dir_in}/results/${TEMP_NAME}
			prodigal -i ${TEMP}/${TEMP_NAME} -a ${dir_in}/results/${TEMP_NAME}/all.pdg.faa -f gff -o ${dir_in}/results/${TEMP_NAME}/all.pdg.gff -p meta
			AAFILE=${dir_in}/results/${TEMP_NAME}/all.pdg.faa
			
			cp ${TEMP}/${TEMP_NAME} ${dir_in}/results/${TEMP_NAME} # this will preserve host info
        	        # Run Blastp 
        	        echo "Performing BLASTp on ${f##*/}..."
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

        	        # Run PICI typer script 
        	        python ${scripts_path}/prototype_typer.py --i $integrase_identity --a $alpa_identity --blast "${TEMP}/BLASTp_results.out" --fasta ${TEMP}/${TEMP_NAME} --output "${TEMP}/PICI_results" --aa "${TEMP}/python/all.pdg.faa" --temp ${TEMP}
			
			# Remove duplicates
			python ${scripts_path}/duplicate_remover.py --input "${TEMP}/PICI_results" --output "${TEMP}/PICI_results"
        	        
        	        # Move PICI results to PICI directory
			cp "${TEMP}/tBLASTn_results.out" "${dir_in}/tmp_PICIs/${TEMP_NAME}/"
			cp ${TEMP}/BLASTp_results.out "${dir_in}/tmp_PICIs/${TEMP_NAME}/"
        	        mv "${TEMP}/PICI_results" "${dir_in}/tmp_PICIs/${TEMP_NAME}/"
        	        echo "Finished ${TEMP_NAME}..."
        	fi
	fi
done


#Exit the typer if no PICIs were found
if [ ! "$(ls -A ${dir_in}/results)" ]; then
    echo "${dir_in}/results is  empty"
    exit 0
fi


echo "Compiling all PICIs into one file 'ALL_PICIs.fasta'..."
find ${dir_in}/tmp_PICIs/ -name PICI_results -exec cat {} \; >> ${dir_in}/results/ALL_PICIs.fasta


echo "Separating PICI type and phage satellites..."
python ${scripts_path}/pici_separator.py --input ${dir_in}/results/ALL_PICIs.fasta --output ${dir_in}/results/

echo "Reviewing phage satellites..."
if [ "$database" -eq "0" ]; then
	blastx -query ${dir_in}/results/Phage_Satellites.fasta -subject ${db_prot_path} -task blastx -evalue 0.001 -outfmt 6 -out ${TEMP}/satellite_BLAST_results.out
elif [ "$database" -eq "1" ]; then
	blastn -query ${dir_in}/results//Phage_Satellites.fasta -subject ${db_nuc_path} -task blastn -evalue 0.001 -outfmt 6 -out ${TEMP}/satellite_BLAST_results.out
fi
python ${scripts_path}/phage_satellite_review.py --a $alpa_identity --output-gram-negative ${dir_in}/results/G_neg_PICI_reviewed --output-reviewed ${dir_in}/results/phage_satellite_reviewed --blast ${TEMP}/satellite_BLAST_results.out --fasta ${dir_in}/results/Phage_Satellites.fasta
sed -i -- 's/phage_satellite/G_neg_PICI/g' ${dir_in}/results//G_neg_PICI_reviewed
cat ${dir_in}/results/G_neg_PICI_reviewed ${dir_in}/results/G_neg_PICIs.fasta ${dir_in}/results/SaPIs.fasta ${dir_in}/results/phage_satellite_reviewed > ${dir_in}/results/ALL_PICIs.fasta

echo "Creating table..."
cat ${dir_in}/results/ALL_PICIs.fasta | grep -e "^>" | sed 's/>//g' | sed 's/;/\t/g' > ${dir_in}/results/PICI_table.tsv

echo "Done."


prefix=$(basename ${dir_in})

for file in $(ls -p  ${dir_in}/results/ | grep -v /)
do
    mv ${dir_in}/results/${file} ${dir_in}/results/${prefix}_${file}
done


