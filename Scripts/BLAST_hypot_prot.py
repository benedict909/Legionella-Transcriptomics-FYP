#!/usr/bin/python3.6

from Bio.Blast import NCBIWWW
from Bio import SeqIO
import sys
import os
import shutil

directory_in_str = "/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/FASTA_for_BLAST"
directory = os.fsencode(directory_in_str)

resdir = "/mnt/c/Users/bened/Dropbox/Final_Year/Project/Data/BLAST_results/"

print("running BLASTp search on input protein(s)")

for file in os.listdir(directory):

	filename = os.fsdecode(file)
	filepath = directory_in_str + "/" +  filename

	if os.path.isfile(filepath):

		print(filename)
		
		gene = SeqIO.read(filepath, "fasta") # load sequence record from a FASTA file
		locus = filename.replace(".fasta", "")

		# Tabular Result Format 
		result_tab = NCBIWWW.qblast("blastp", "nr", gene.format("fasta"), 
			hitlist_size=10, format_type="Tabular",	expect=1e-3)

		# Open a file for the BLAST output
		output = resdir + locus + "_blast_output_tab.out"
		blast_file = open(output,"w")
		blast_file.write(result_tab.read()) # write output
		blast_file.close() # close connection  


		# Text Result Format 
		result_text = NCBIWWW.qblast("blastp", "nr", gene.format("fasta"), 
			hitlist_size=10, format_type="Text", expect=1e-3)

		# Open a file for the BLAST output
		output = resdir + locus + "_blast_output_text.out"
		blast_file = open(output,"w")
		blast_file.write(result_text.read()) # write output
		blast_file.close() # close connection 

		completed_path = directory_in_str + "/completed/" + filename
		shutil.move(filepath, completed_path) # move FASTA of completed BLAST search to a seperate directory 