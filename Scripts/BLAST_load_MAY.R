fastafiles = list.files(file.path(datadir,"BLAST_results_cl"))
blacklist = c("lpg0076.fasta_result", "lpg0364.fasta_result", "lpg0406.fasta_result", "lpg0634.fasta_result",
              "lpg1496.fasta_result","lpg1683.fasta_result", "lpg2400.fasta_result", "lpg2765.fasta_result", "lpg2867.fasta_result")
fastafiles = fastafiles[-which(fastafiles %in% blacklist)]

fasta_con_all = NA
for(eff in fastafiles){
  print(eff)
  tmp_fasta = read.table(file.path(datadir,"BLAST_results_cl",eff), sep = "\t")
  if(all(is.na(fasta_con_all))){
    fasta_con_all = tmp_fasta
  }else{
    fasta_con_all = rbind(fasta_con_all, tmp_fasta)
  }
}
colnames(fasta_con_all) = unlist(strsplit("qseqid sseqid sacc ssciname stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"," "))

newloci = c()
load(file.path(datadir, "T4SS_effectors.RData"))
effectors_loci = c(unlist(strsplit(as.character(effectors_df$L..pneumophila), " ")), newloci); length(effectors_loci)


blast_eff = fasta_con_all %>% 
  select(locus_tag = qseqid, blast_result  = stitle) %>% 
  distinct() %>% 
  filter(grepl("Dot", blast_result)) %>% 
  mutate(new = ifelse(locus_tag %in% effectors_loci, FALSE, TRUE))

newloci = unique(as.character(blast_eff$locus_tag[blast_eff$new == TRUE]))
effectors_loci = c(unlist(strsplit(as.character(effectors_df$L..pneumophila), " ")), newloci); length(effectors_loci)
