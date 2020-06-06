#!/bin/bash

# Script for automating BLASTp search of protein sequences against non-redundant database

Dir=/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/FASTA_for_BLAST/ #directory containing FASTA sequences to blast
Fastas=$Dir*.fasta

resdir=/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/BLAST_results_cl/ #path to directory for BLAST output

for f in $Fastas
do
  locus=${f/$Dir/}
  resulttag=_result

  echo "Running BLAST search for $locus..."

  outfile=$resdir$locus$resulttag

  # run BLASTp, -outfmt controls columns of output table
  blastp -query $f -db "nr" -remote \
  -outfmt "6 qseqid sseqid sacc ssciname stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -evalue 1e-3 -num_alignments 10 -out $outfile       

  mv $f /mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/FASTA_for_BLAST/completed/ #move FASTA of completed BLAST search to new dir

done

