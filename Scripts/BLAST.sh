#!/bin/bash
Dir=/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/FASTA_for_BLAST/
Fastas=$Dir*.fasta

resdir=/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/BLAST_results_cl/

for f in $Fastas
do
  locus=${f/$Dir/}
  resulttag=_result

  echo "Running BLAST search for $locus..."

  outfile=$resdir$locus$resulttag

  #echo $outfile

  # take action on each file. $f store current file name
  # head $f
  blastp -query $f -db "nr" -remote \
  -outfmt "6 qseqid sseqid sacc ssciname stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -evalue 1e-3 -num_alignments 10 -out $outfile       

  mv $f /mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/FASTA_for_BLAST/completed/

done

