touch ALL_PICIs.fasta

for f in ./PICIs/*
do
	cat ./PICIs/${f##*/}/PICI_results >> ALL_PICIs.fasta
  
done
