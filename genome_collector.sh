touch Host_genomes.fasta

for f in ./PICIs/*
do
        if grep -q ">" ./PICIs/${f##*/}/PICI_results; then
                cat ./sequences/${f##*/} >> Host_genomes.fasta
        else
                echo "Empty PICI file for ${f##*/}."
        fi
  
done

echo "Making list of host names for ALL_PICIs.fasta. Storing in ALL_PICIs_host_names.txt..."
grep -oP '(?<=>)(.*?)(?=;)' ALL_PICIs.fasta > ALL_PICIs_host_names.txt

echo "Linearizing Host_genomes.fasta for extraction..."
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < Host_genomes.fasta > linear.fasta

echo "Extracting specific contigs from host genomes..."
grep -w -A 1 -Ff ALL_PICIs_host_names.txt linear.fasta --no-group-separator > tmp_ALL_PICIs_host_genomes.fasta

echo "Formatting back to fasta..."
tr "\t" "\n" < tmp_ALL_PICIs_host_genomes.fasta > ALL_PICIs_host_genomes.fasta

rm tmp_ALL_PICIs_host_genomes.fasta

echo "Host genome retro-extraction completed. Check ALL_PICIs_host_names.txt and ALL_PICIs_host_genomes.fasta for results."


if test -f "G_neg_PICIs.fasta"; then
        echo "G_neg_PICIs.fasta exists. Gathering genomes..."
        echo "Making list of host names for G_neg_PICIs.fasta. Storing in G_neg_PICIs_host_names.txt..."
        grep -oP '(?<=>)(.*?)(?=;)' G_neg_PICIs.fasta > G_neg_PICIs_host_names.txt

        echo "Extracting specific contigs from host genomes..."
        grep -w -A 1 -Ff G_neg_PICIs_host_names.txt linear.fasta --no-group-separator > tmp_G_neg_PICIs_host_genomes.fasta

        echo "Formatting back to fasta..."
        tr "\t" "\n" < tmp_G_neg_PICIs_host_genomes.fasta > G_neg_PICIs_host_genomes.fasta

        rm tmp_G_neg_PICIs_host_genomes.fasta

        echo "Host genome retro-extraction completed. Check G_neg_PICIs_host_names.txt and G_neg_PICIs_host_genomes.fasta for results."
else
        echo "No G_neg_PICIs.fasta detected."
fi



if test -f "SaPIs.fasta"; then
        echo "SaPIs.fasta exists. Gathering genomes..."
        echo "Making list of host names for SaPIs.fasta. Storing in SaPIs_host_names.txt..."
        grep -oP '(?<=>)(.*?)(?=;)' SaPIs.fasta > SaPIs_host_names.txt

        echo "Extracting specific contigs from host genomes..."
        grep -w -A 1 -Ff SaPIs_host_names.txt linear.fasta --no-group-separator > tmp_SaPIs_host_genomes.fasta

        echo "Formatting back to fasta..."
        tr "\t" "\n" < tmp_SaPIs_host_genomes.fasta > SaPIs_host_genomes.fasta

        rm tmp_SaPIs_host_genomes.fasta

        echo "Host genome retro-extraction completed. Check SaPIs_host_names.txt and SaPIs_host_genomes.fasta for results."
else
        echo "No SaPIs.fasta detected."
fi



if test -f "Phage_Satellites.fasta"; then
        echo "Phage_Satellites.fasta exists. Gathering genomes..."
        echo "Making list of host names for Phage_Satellites.fasta. Storing in Phage_Satellites_host_names.txt..."
        grep -oP '(?<=>)(.*?)(?=;)' Phage_Satellites.fasta > Phage_Satellites_host_names.txt

        echo "Extracting specific contigs from host genomes..."
        grep -w -A 1 -Ff Phage_Satellites_host_names.txt linear.fasta --no-group-separator > tmp_Phage_Satellites_host_genomes.fasta

        echo "Formatting back to fasta..."
        tr "\t" "\n" < tmp_Phage_Satellites_host_genomes.fasta > Phage_Satellites_host_genomes.fasta

        rm tmp_Phage_Satellites_host_genomes.fasta

        echo "Host genome retro-extraction completed. Check Phage_Satellites_host_names.txt and Phage_Satellites_host_genomes.fasta for results."
else
        echo "No Phage_Satellites.fasta detected."
fi

rm linear.fasta
rm Host_genomes.fasta

