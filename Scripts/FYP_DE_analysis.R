# Final year project data analysis May 2020
# RNA-Seq DE analysis of Legionella pneumophila in response to sources of environmental stress

source('C:/Users/bened/Dropbox/Final_Year/Project/Scripts/FYP_DE_functions.R') # path to source script

install_necessary_packages() # install any currently uninstalled packages that are required for this script

library(tidyverse)
library(gdata)
library(reshape2)
library(data.table)
library(gage) # for KEGG & Gene set enrichment analysis
library(ComplexHeatmap) # for heatmap plots 
library(circlize) # for colour functions 


# Load data ---------------------------------------------------------------------------------------
setwd("~")
datadir = "../Dropbox/Final_Year/Project/Data/" # path to data directory
resdir_parent = "../Dropbox/Final_Year/Project/Results/" # path to desired parent results directory 

load(file.path(datadir,"DE_KEGG.RData"))
load(file.path(datadir,"DE_genes.RData"))
load(file.path(datadir,"KEGG_pathway_lists_concat/general_pthwys_list.RData"))
load(file.path(datadir,"KEGG_pathway_lists_concat/pthwys_list.RData"))


# Set parameters ----------------------------------------------------------------------------------

fold_filter = c(1.5, log2(1.5))[2]
p_filter = 0.05

fold_suffix = "_log2FC"
p_suffix = "_pvalue"

DE_cols = names(DE_genes) 
cond_colnames = DE_cols[grepl("log2|pvalue|change", DE_cols)]
conditions_list = unique(substr(cond_colnames,1,2))

mycol = c("firebrick3","firebrick3","blue3","blue3","blue3", "grey45")
names(mycol) = c("up", "upregulated", "downregulated","dn", "down", "non-DE")

foldchanges = cond_colnames[grepl("log2", cond_colnames)]
pcols =       cond_colnames[grepl("pvalue", cond_colnames)]
changecols =  cond_colnames[grepl("_change", cond_colnames)]

allgenes_folds = DE_genes %>% # make a matrix of just fold changes
  mutate(gene = as.character(locus_tag)) %>% 
  dplyr::select(gene, all_of(foldchanges)) %>% 
  arrange(gene) %>% 
  column_to_rownames("gene"); colnames(allgenes_folds) = conditions_list

allgenes_changes = DE_genes %>% # make a matrix of just fold changes
  mutate(gene = as.character(locus_tag)) %>% 
  dplyr::select(gene, all_of(changecols)) %>% 
  arrange(gene) %>% 
  column_to_rownames("gene"); colnames(allgenes_folds) = conditions_list

#  1. Make stats table --------------------------------------------------------------------------------

DE_table = (matrix(nrow = 13, ncol = length(conditions_list)))
colnames(DE_table) = conditions_list
row.names(DE_table) = c("n_up","n_hi_up","mean_raw_up","sd_raw_up","mean_up","sd_up","n_dn","n_hi_dn","mean_raw_dn","sd_raw_dn","mean_dn","sd_dn","total")


for(i in 1:length(conditions_list)){
  condition = conditions_list[i]; print(condition)
  fold_colname = paste0(condition,fold_suffix)
  p_colname = paste0(condition,p_suffix)
  
  up_tmp = DE_genes %>% 
    filter(!is.na(fold_colname) & !!as.name(fold_colname) >= fold_filter &
             !!as.name(p_colname) <=p_filter)
  
  DE_table["n_up", condition] = nrow(up_tmp)
  DE_table["n_hi_up", condition] = up_tmp %>% filter(!!as.name(fold_colname)>= 2) %>% nrow()
  DE_table["mean_raw_up", condition] = mean(up_tmp %>% pull(!!as.name(fold_colname)) %>% 2^.)
  DE_table["sd_raw_up", condition] = sd(up_tmp %>% pull(!!as.name(fold_colname)) %>% 2^.)
  DE_table["mean_up", condition] = mean(up_tmp %>% pull(!!as.name(fold_colname)))
  DE_table["sd_up", condition] = sd(up_tmp %>% pull(!!as.name(fold_colname)))
  
  dn_tmp = DE_genes %>% 
    filter(!is.na(fold_colname) & !!as.name(fold_colname) <= -fold_filter &
             !!as.name(p_colname) <=p_filter)
  
  DE_table["n_dn", condition] = nrow(dn_tmp)
  DE_table["n_hi_dn", condition] = dn_tmp %>% filter(!!as.name(fold_colname)<= -2) %>% nrow()
  DE_table["mean_raw_dn", condition] = -mean(dn_tmp %>% pull(!!as.name(fold_colname)) %>% abs() %>% 2^.)
  DE_table["sd_raw_dn", condition] = sd(dn_tmp %>% pull(!!as.name(fold_colname)) %>% abs() %>% 2^.)
  DE_table["mean_dn", condition] = mean(dn_tmp %>% pull(!!as.name(fold_colname)))
  DE_table["sd_dn", condition] = sd(dn_tmp %>% pull(!!as.name(fold_colname)))
  
  DE_table["total", condition] = DE_table["n_up", condition] + DE_table["n_dn", condition]
}
t(DE_table)

write.csv(t(DE_table), file = file.path(resdir_parent, "summary_table.txt"), quote = F)






#  2. Summary boxplots ----

resdir = file.path(resdir_parent, "summary_box"); if(!file.exists(resdir)) dir.create(resdir)

box_in = allgenes_folds %>%
  reshape2::melt() 

tmp = DE_genes %>% 
  dplyr::select(locus_tag, all_of(foldchanges), all_of(changecols)) %>% 
  drop_na(all_of(foldchanges)) %>% 
  reshape2::melt(id.vars = c("locus_tag"), measure.vars = c(foldchanges, changecols)) %>% 
  mutate(cond = substr(variable,1,2),
         unique_tag = paste0(cond, "_",locus_tag) ) 

change_df = tmp %>% 
  filter(!grepl("log2",variable)) %>% 
  dplyr::select(value, unique_tag)

box_in = tmp %>%
  filter(grepl("log2",variable)) %>% 
  left_join(change_df, by = "unique_tag") %>% 
  rename(change = value.y, foldchange = value.x) %>% 
  mutate(foldchange = as.numeric(foldchange),
         direction = ifelse(foldchange >= 0, "positive", "negative"),
         #log2fc = log2(abs(foldchange)),
         log2fc = foldchange,
         # log2fc = ifelse(direction == "negative", -log2fc, log2fc)
  ) %>% 
  dplyr::add_count(cond, change,name = "n_nonDE") %>% 
  mutate(n_nonDE = ifelse(change == "non-DE", n_nonDE, NA),
         condnum = as.numeric(factor(cond)))

defaultcol = hue_pal()(9); names(defaultcol) = unique(box_in$cond)
box_in$jitcol = defaultcol[match(box_in$cond,names(defaultcol))]



png(file.path(resdir, "summ_box_rb_nogrid.png"), width = 2000, height = 1400, res = 200)
ggplot() +
  geom_segment(aes(x = 0.6, y =  0.5,xend = 10.5, yend=  0.5), linetype = "dashed", colour = "grey48") + 
  geom_segment(aes(x = 0.6, y = -0.5,xend = 10.5, yend= -0.5), linetype = "dashed", colour = "grey48") + 
  geom_segment(aes(x = 0.6, y =  2.0,xend = 10.5, yend=  2.0), linetype = "dashed", colour = "grey48") + 
  geom_segment(aes(x = 0.6, y = -2.0,xend = 10.5, yend= -2.0), linetype = "dashed", colour = "grey48") + 
  
  ggpol::geom_boxjitter(filter(box_in, change == "upregulated"),
                        mapping = aes(x = condnum, y = log2fc, group = cond),
                        fill = "grey87", jitter.fill = "firebrick3",
                        jitter.shape = 21, jitter.size = 2, #show.legend = T,
                        outlier.shape = NA, width = .5, 
                        lwd = 1, fatten = 1,
                        errorbar.draw = T, errorbar.length = 0.4) +
  ggpol::geom_boxjitter(filter(box_in, change == "downregulated"),
                        mapping = aes(x = condnum, y = log2fc, group = cond), 
                        jitter.shape = 21, width = .5,
                        lwd = 1, fatten = 1,
                        fill = "grey87", jitter.fill = "blue3", #show.legend = T,
                        outlier.shape = NA, jitter.size = 2, errorbar.draw = T, errorbar.length = 0.4) +
  
  scale_y_continuous(limits = c(-8,12), expand = c(0,0),# position = "right",
                     breaks = c(seq(-8,12,2))) + 
  scale_x_continuous(breaks = c(1:9), labels = unique(box_in$cond), expand = c(0,0), limits = c(.6,10.5)) + 
  geom_text(box_in, mapping = aes(x = condnum, y = 0, label = n_nonDE),
            position = "identity", na.rm = T, check_overlap = T, size = 3.5) +
  # geom_label("text", mapping = aes(data = box_in, x = condnum, y = 0, label = n_nonDE),
  #            position = "identity", na.rm = T) + 
  xlab("Stress Condition") + ylab(expression(log[2](FC))) + 
  coord_cartesian(clip = 'off') + # 'off' keeps the labels from disappearing #xlim = c(.6, 9.25), # This focuses the x-axis on the range of interest
  
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(), panel.grid.minor.x = element_blank(),legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 15), plot.margin=unit(c(0.02,0.02,0,0),"npc"), axis.title.x = element_text(size = 15),
        #plot.margin=unit(c(0,2,0,0),"cm"),  # order = top, right, bottom, left (/never eat shredded wheat)) 
        axis.text.y = element_text(vjust = .25, size = 15), axis.title.y = element_text(hjust = 0.55, size = 15) )+
  
  annotate("text",x = 10, y =  0.00,colour = "black", label = "no. of non-DEGs") +
  annotate("text",x = 10, y =  1.25,colour = "black", label = "upregulated") +
  annotate("text",x = 10, y = -1.25,colour = "black", label = "downregulated") +
  annotate("text",x = 10, y =  2.75,colour = "black", label = "highly upreg.") +
  annotate("text",x = 10, y = -2.75,colour = "black", label = "highly downreg.") 
dev.off()


#  3. Heatshock and Fe proteins ----

resdir = file.path(resdir_parent, "HS_&_Fe"); if(!dir.exists(resdir)) dir.create(resdir)

# Heat shock
hs_genes = DE_genes$locus_tag[grepl("heat s", DE_genes$product, ignore.case = T) | grepl("heat s", DE_genes$KO_description, ignore.case = T)]

hs_df = DE_genes %>% 
  filter(locus_tag %in% hs_genes) %>% 
  mutate(gene = as.character(gene))

hs_df$gene[hs_df$product == "heat shock hsp20"] = "hsp20"
hs_df$gene[hs_df$product == "small HspC2 heat shock protein"] = "HspC2.1"
hs_df$gene[hs_df$product == "small heat shock protein HspC2"] = "HspC2.2"

hs_mat = hs_df %>% 
  column_to_rownames("gene") %>%
  select(all_of(changecols) ) %>% 
  t() %>% 
  cbind(n_up =  NA, n_dn = NA)

for(i in 1:nrow(hs_mat)){
  hs_mat[i,"n_up"] = as.numeric(length(which(hs_mat[i,] == "upregulated")))
  hs_mat[i,"n_dn"] = as.numeric(length(which(hs_mat[i,] == "downregulated")))
}
  
hs_plot = hs_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("cond") %>% 
  reshape2::melt("cond") %>% 
  filter(variable %in% c("n_up", "n_dn")) %>% 
  mutate(n_change = ifelse(variable == "n_dn", -as.numeric(value), as.numeric(value)),
         variable = ifelse(variable == "n_dn", "downregulated", "upregulated"),
         cond = substr(cond,1,2))
  
png(file.path(resdir, "hs_bars.png"), width = 600, height = 600, res = 200)
ggplot(hs_plot) + 
  geom_hline(yintercept = 0, linetype = "solid")+ 
  geom_bar(stat='identity', aes(x=factor(cond, levels = sort(unique(cond), decreasing = F)),
                                y=n_change, fill = variable), width=.5, colour = "black")  +
  scale_fill_manual(values = mycol, limits = c("upregulated", "downregulated")) +
  scale_y_continuous(limits = c(-3.25,10), expand = c(0,0), breaks = seq(-2, 10, 2),
                     labels = c(2,0,seq(2,10,2))) + 
  theme_bw() + ylab("No. of heat shock genes") + xlab("Condition")+
  theme(legend.title = element_blank(),
        legend.position = "none", 
        panel.border = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), #plot.margin=unit(c(0,.03,0,0),"npc"),
        panel.grid = element_blank())
dev.off()



fe_genes = DE_genes$locus_tag[grepl("iron", DE_genes$product, ignore.case = T) | grepl("iron", DE_genes$KO_description, ignore.case = T)]

Fe_acq_loci = c("lpg0858","lpg0124","lpg2657","lpg2800","lpg0232","lpg0024","lpg1998","lpg1285", # from Cianciotto, Future Microbiol. (2015) 10(5), 841-851 
                "lpg0746", "lpg2815","lpg1325","lpg1324","lpg1323","lpg1326","lpg2278","lpg0265","lpg2647","lpg0467","lpg2906") 

fe_functions = list("lpg0858" = "Legiobactin production",
                    "lpg0124" = "Legiobactin production",
                    "lpg2657" = "Ferrous iron uptake",
                    "lpg2800" = "Putative siderophore synthetase",
                    "lpg0232" = "Transcriptional regulator",
                    "lpg0024" = "Hemin-binding protein",
                    "lpg1998" = "HGA production",
                    "lpg1285" = "HGA degradation",
                    "lpg0746" = "Putative peptide transporter",
                    "lpg2815" = "Ferrous iron uptake",
                    "lpg1325" = "Legiobactin biosynthesis",
                    "lpg1324" = "IM export of legiobactin",
                    "lpg1323" = "IM import of legiobactin",
                    "lpg1326" = "OM receptor for legiobactin",
                    "lpg2278" = "HGA production",
                    "lpp0339" = "lpg0265 Multicopper oxidase",
                    "lpg2647" = "HGA production",
                    "lpg0467" = "Degrades transferrin",
                    "lpg2906" = "Growth on low-iron media")

fe_func_df = as.data.frame(fe_functions) %>% t() %>% as.data.frame() %>% rownames_to_column("locus_tag")


Fe_df =  DE_genes %>%
  filter(locus_tag %in% Fe_acq_loci) %>%
  select(locus_tag,gene,product,KO_description, Li_log2FC, Li_pvalue, Li_change) %>% 
  left_join(fe_func_df) %>% 
  arrange(desc(Li_log2FC))

write.csv(Fe_df, file = file.path(resdir, "Fe_acq_genes.csv"), row.names = F)


#  4. KEGG analysis May 2020 ----

resdir. = file.path(resdir_parent, "KEGG"); if(!dir.exists(resdir.)) dir.create(resdir.)

# KEGG pathway enrichment analysis
KEGG_df = DE_KEGG %>% dplyr::select(locus_tag, KEGG_pthwy_grp) %>% filter(!is.na(KEGG_pthwy_grp)) %>% distinct()
KEGG_sets = split(KEGG_df, KEGG_df$KEGG_pthwy_grp)
for(item in names(KEGG_sets)){
  currentitem = data.frame(KEGG_sets[item])
  currentitemloci = currentitem %>% pull(paste0(item,".locus_tag"))
  KEGG_sets[item] = NA
  KEGG_sets[item] = list(currentitemloci)
}

keggres = gage(allgenes_folds, gsets= KEGG_sets, same.dir=T) # df of p-values 
# View(keggres$greater); View(keggres$less); View(keggres$stats)

generaldf = reshape2::melt(general_pthwys_list) %>% 
  rename("pathway" = value, "general_pth" = L1)

keggres_dn = reshape2::melt(t(keggres$less)) %>% 
  filter(Var1 %!in% c("As", "Os", "Tm"), Var2 != "Global_and_overview_maps") %>%
  rename(cond = Var1, pathway = Var2, pval = value) %>% 
  mutate(variable = "dn")

keggres_df = reshape2::melt(t(keggres$greater)) %>%
  filter(Var1 %!in% c("As", "Os", "Tm"), Var2 != "Global_and_overview_maps") %>% 
  rename(cond = Var1, pathway = Var2, pval = value) %>%
  mutate(variable = "up") %>% 
  rbind(keggres_dn) %>% 
  arrange(cond, pathway)

input = plot_nchanges_input(paths = general_pthwys_list, input_df = DE_KEGG,
                            input_path_col = "KEGG_pthwy_grp", conditions = conditions_list) %>% 
  filter(variable != "non-DE") %>% 
  left_join(generaldf) %>% 
  filter(pathway != "Global_and_overview_maps") %>% 
  filter(cond %!in% c("As", "Os", "Tm")) %>% 
  left_join(keggres_df, by = c("pathway", "cond", "variable")) %>% 
  mutate(signif_up = ifelse(variable == "up" & pval <= 0.05, "*", NA),
         signif_dn = ifelse(variable == "dn" & pval <= 0.05, "*", NA))


input = input %>% 
  mutate(perc = prop * 100,
         pathway = gsub("_"," ", pathway),    ## edit pathway labels for plot 
         general_pth =  gsub("_"," ", general_pth),
         pathway = gsub("metabolism","met.", pathway),
         pathway = gsub("Metabolism","Met.", pathway),
         pathway = gsub("and","&", pathway),
         pathway = gsub("degradation","degrad.", pathway), 
         pathway = gsub(" prokaryotes","", pathway),
         pathway = gsub("biosynthesis","Biosynth.", pathway, ignore.case = T),
         pathway = gsub("secondary","2\u00B0", pathway),
         general_pth = gsub("Information","Info", general_pth),
         general_pth = gsub("l I","l\nI", general_pth),
         general_pth = gsub("r P","r\nP", general_pth),
         general_pth = gsub("n D","n\nD", general_pth),
         general_pth = gsub("Genetic Info P","Genetic Info.\nP", general_pth),
         general_pth = gsub("eases","ease", general_pth))

png(file = file.path(resdir.,paste0("KEGG_pathways_fig_",gsub("-","",Sys.Date()),".png")),  width = 2000, height = 2250, res = 200)
ggplot(input) +
  geom_bar(mapping = aes(x = factor(pathway, levels = sort(unique(pathway),decreasing = T)),
                         y =  perc, fill = variable,
                         group = factor(variable, levels = sort(c("up","dn"),decreasing = F))),
           stat = "identity", position = "dodge", width = 0.75, colour = "black") +
  scale_fill_manual(values = mycol, limits = c("up","dn")) +
  
  geom_text(mapping = aes(x = factor(pathway, levels = sort(unique(pathway),decreasing = T)), # add * labels for stat. sig. 
                          y = perc, label = signif_up, group = factor(variable, levels = sort(c("up","dn"),decreasing = F))),
             na.rm = T, check_overlap = T, size = 6, nudge_x = 0.15, nudge_y = 6) +
  geom_text(mapping = aes(x = factor(pathway, levels = sort(unique(pathway),decreasing = T)),
                          y = perc, label = signif_dn, group = factor(variable, levels = sort(c("up","dn"),decreasing = F))),
            na.rm = T, check_overlap = T, size = 6, nudge_x = -.25, nudge_y = 6) +
  
  #facet_grid(cond ~ factor(general_pth), scales = "free", space = "free_x") + # pathways in x axis, or:
  facet_grid(factor(general_pth)~ cond, scales = "free", space = "free_y") + coord_flip() + # pathways in y axis 
  xlab("KEGG Pathway") + ylab("Percentage of genes") + 
  scale_y_continuous(expand = expansion(mult = c(0,0), add = 0), limits = c(0,105),breaks = c(0,50,100)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -90, size = 10, hjust = 0, vjust = 0.5), 
        panel.border = element_rect(size = 1),
        axis.text.y = element_text(size = 9), axis.title.y = element_text(vjust = 0),
        axis.line = element_line(), axis.title = element_text(size = 15),
        panel.grid.minor = element_blank(),  panel.grid.major = element_line(colour = "grey80"), 
        strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 10),
        # legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 15) #uncomment for legend
        legend.position = "none" # no legend
  )
dev.off()

input %>% filter(pathway == "Cell motility") %>% dplyr::select(-c(general_pth))  %>% arrange(desc(value))
input %>% filter(pathway == "Transcription") %>% dplyr::select(-c(general_pth))  %>% arrange(desc(value))
input %>% filter(grepl("Infectious",pathway)) %>% dplyr::select(-c(general_pth)) %>% arrange(desc(value))
input %>% filter(grepl("Energy",pathway)) %>% dplyr::select(-c(general_pth)) %>% arrange(desc(value))



#  5. Clustering heatmap ----

resdir = file.path(resdir_parent,"clustering_may"); if(!file.exists(resdir)) dir.create(resdir)


# Load Dot/Icm effector locus tag names
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

blast_eff = fasta_con_all %>% 
  select(locus_tag = qseqid, blast_result  = stitle) %>% 
  distinct() %>% 
  filter(grepl("Dot", blast_result)) %>% 
  mutate(new = ifelse(locus_tag %in% effectors_loci, FALSE, TRUE))

newloci = unique(as.character(blast_eff$locus_tag[blast_eff$new == TRUE]))
effectors_loci = c(unlist(strsplit(as.character(effectors_df$L..pneumophila), " ")), newloci); length(effectors_loci)

eff_presence = rownames(allgenes_folds) %in% effectors_loci

# T4SS machinery annot
load(file.path(datadir, "T4SS_genes.RData"))
additional_loci = c("lpg0525", "lpg0472") #lphA = icmN, lpg0525 = LvgA, lpg0472 = dotV

T4SS_loci = DE_genes %>% filter(gene %in% T4SS_genes) %>% pull(locus_tag)
T4SS_loci = c(T4SS_loci, additional_loci)
T4SS_presence = as.character(rownames(allgenes_folds) %in% T4SS_loci)

# KEGG pathways annot
met_loci = DE_KEGG %>% filter(KEGG_pthwy_general_grp == "Metabolism") %>% pull(locus_tag)
met_presence = as.character(rownames(allgenes_folds) %in% met_loci)

Genet_loci = DE_KEGG %>% filter(KEGG_pthwy_general_grp == "Genetic_Information_Processing") %>% pull(locus_tag)
Genet_presence = as.character(rownames(allgenes_folds) %in% Genet_loci)

Cellul_loci = DE_KEGG %>% filter(KEGG_pthwy_general_grp == "Cellular_Processes") %>% pull(locus_tag)
Cellul_presence = as.character(rownames(allgenes_folds) %in% Cellul_loci)

enviro_loci = DE_KEGG %>% filter(KEGG_pthwy_general_grp == "Environmental_Information_Processing") %>% pull(locus_tag)
enviro_presence = as.character(rownames(allgenes_folds) %in% enviro_loci)

Diseas_loci = DE_KEGG %>% filter(KEGG_pthwy_general_grp == "Human_Diseases") %>% pull(locus_tag)
Diseas_presence = as.character(rownames(allgenes_folds) %in% Diseas_loci)


# Make row annot
row_ha = rowAnnotation(T4SS_eff. = eff_presence, T4SS = T4SS_presence, 
                       Met = met_presence, Genet = Genet_presence, Cellul = Cellul_presence,
                       Enviro = enviro_presence, Diseas = Diseas_presence,
                       col = list(T4SS_eff. = c("TRUE" = "black", "FALSE" = "white"),
                                  T4SS = c("TRUE" = "black", "FALSE" = "white"),
                                  Met = c("TRUE" = "red", "FALSE" = "white"),
                                  Genet = c("TRUE" = "blue", "FALSE" = "white"),
                                  Cellul = c("TRUE" = "chartreuse1", "FALSE" = "white"),
                                  Enviro = c("TRUE" = "darkorange1", "FALSE" = "white"),
                                  Diseas = c("TRUE" = "mediumorchid1", "FALSE" = "white")),
                       border = TRUE)

# heatmap colour function, change num range to change where colour is max
# col_fun = colorRamp2(c(-2,-0.4, 0,0.4, 2), c("blue3","white", "white","white", "firebrick3"))
# col_fun = colorRamp2(c(-2, 0, 2), c("blue3","white", "firebrick3"))
col_fun2 = colorRamp2(c(-2,-0.4, 0,0.4, 2), c("blue3","white", "white","white", "firebrick3"))


# # plot heatmap for all 9 conditions
# allgeneshc = Heatmap((as.matrix(allgenes_folds[,])), 
#                      show_row_names = FALSE, name = "Fold change", border = T, col = col_fun,
#                      clustering_distance_columns = cosine_distance, clustering_method_columns = "ward.D",
#                      clustering_distance_rows = cosine_distance, clustering_method_rows = "ward.D",
#                      right_annotation = row_ha,
#                      column_dend_height = unit(40, "mm"),
#                      row_dend_width = unit(40, "mm"), row_split = 2,
#                      # cluster_rows = F,
#                      heatmap_width = unit(0.75, "npc"))
# 
# pdf(file = file.path(resdir, "all_genes_clust_all_cond.pdf"), width = 10, height = 15)
# allgeneshc
# dev.off()

# plot heatmap for 6 most stressful conditions 
cartedechaleur = Heatmap(t(as.matrix(allgenes_folds[,c("Li","Mi","Nd","Ns","Ox","Sp")])), 
                         show_column_names = FALSE, name = "log2 FC", border = T,
                         col = col_fun2, column_title = "Genes", row_title = "Stress Condition", 
                         clustering_distance_columns = cosine_distance, clustering_method_columns = "ward.D",
                         clustering_distance_rows = cosine_distance, clustering_method_rows = "ward.D",
                         #right_annotation = row_ha,
                         column_dend_height = unit(0.2, "npc"), row_dend_width = unit(0.1, "npc"), # control heights of dendrograms relative to heatmap
                         #row_dend_width = unit(0.1, "npc"), row_split = 2,
                         column_split = 2, column_gap = unit(0.01, "npc"),
                         heatmap_width = unit(1, "npc"),
                         row_names_gp = gpar(fontsize = 20), # row_title = "Stress Condition", #row_title_side = "right",
                         raster_device = "png")

# pdf(file = file.path(resdir, "cosine_hc.pdf"), width = 10, height = 15)
png(file = file.path(resdir, "cosine_hc.png"), width = 2400, height = 1000, res = 200)
cartedechaleur # plot heatmap
decorate_heatmap_body("log2 FC", column_slice = c(1), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))}) # plot thicker borders
decorate_heatmap_body("log2 FC", column_slice = c(2), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()



# # euclidean distance clustering
# euc_hm = Heatmap((as.matrix(allgenes_folds[,c("Li","Mi","Nd","Ns","Ox","Sp")])), 
#                  show_row_names = FALSE, name = "Fold change", border = T, col = col_fun,
#                  clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D",
#                  clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D",
#                  right_annotation = row_ha,
#                  column_dend_height = unit(40, "mm"),
#                  row_dend_width = unit(40, "mm"), row_split = 3,
#                  heatmap_width = unit(0.75, "npc"))
# 
# pdf(file = file.path(resdir, "euclidean_hc.pdf"), width = 10, height = 15)
# euc_hm
# dev.off()

#  6. Clustering bar ----
resdir = file.path(resdir_parent,"clustering_may"); if(!file.exists(resdir)) dir.create(resdir)

geneord = column_order(cartedechaleur)
clust1_loci = rownames(allgenes_folds)[geneord[[1]]]
clust2_loci = rownames(allgenes_folds)[geneord[[2]]]

overall_counts = DE_KEGG %>% 
  select(locus_tag, KEGG_pthwy_grp) %>% 
  distinct() %>% 
  count(KEGG_pthwy_grp, name = "overall_n")

KEGG_conversion = DE_KEGG %>% 
  select(KEGG_pthwy_grp, KEGG_pthwy_general_grp) %>% filter(!is.na(KEGG_pthwy_grp)) %>% 
  distinct()

clust1_KEGG_counts = KEGG_df %>% 
  filter(locus_tag %in% clust1_loci) %>% 
  count(KEGG_pthwy_grp, name = "clust_n") %>% 
  as.data.frame() %>% 
  mutate(variable = "clust1")

clust_KEGG_counts = KEGG_df %>% 
  filter(locus_tag %in% clust2_loci) %>% 
  count(KEGG_pthwy_grp, name = "clust_n") %>% 
  as.data.frame() %>% 
  mutate(variable = "clust2") %>% 
  rbind(clust1_KEGG_counts) %>% 
  left_join(overall_counts) %>% 
  left_join(KEGG_conversion) %>% 
  rename(pathway = KEGG_pthwy_grp, general_pth = KEGG_pthwy_general_grp) %>% 
  mutate(prop = ifelse(variable=="clust1", -clust_n / overall_n, clust_n / overall_n)*100)#,

T4SS_counts1 = data.frame(general_pth = "Dot/Icm_effectors", variable = "clust1",
                         n = length(effectors_loci[effectors_loci%in%clust1_loci]),
                         ntot = length(effectors_loci))
T4SS_counts2 = data.frame(general_pth = "Dot/Icm_effectors", variable = "clust2",
                         n = length(effectors_loci[effectors_loci%in%clust2_loci]),
                         ntot = length(effectors_loci))
T4SS_counts3 = data.frame(general_pth = "Dot/Icm_machinery", variable = "clust1",
                          n = length(T4SS_loci[T4SS_loci%in%clust1_loci]),
                          ntot = length(T4SS_loci))
T4SS_counts4 = data.frame(general_pth = "Dot/Icm_machinery", variable = "clust2",
                          n = length(T4SS_loci[T4SS_loci%in%clust2_loci]),
                          ntot = length(T4SS_loci))
T4SS_counts = rbind(T4SS_counts1,T4SS_counts2,T4SS_counts3,T4SS_counts4)

keep_KEGG_counts = clust_KEGG_counts %>% 
  filter(pathway %in% c("Cell_motility")) %>% 
  select(general_pth =pathway,variable,n=clust_n,ntot=overall_n,prop)

general_counts = clust_KEGG_counts %>% 
  select(clust_n, general_pth, variable) %>% 
  filter(general_pth %in% c("Genetic_Information_Processing", "Metabolism")) %>% 
  group_by(general_pth, variable) %>% 
  summarise(n = sum(clust_n)) %>% ungroup() %>% 
  group_by(general_pth) %>% 
  add_tally(n, name = "ntot") %>% ungroup() %>% 
  rbind(T4SS_counts) %>% 
  mutate(prop = ifelse(variable == "clust1", -n/ntot, n/ntot)*100) %>% 
  rbind(keep_KEGG_counts) %>% 
  mutate(general_pth =  gsub("_"," ", general_pth),
         general_pth = gsub("l I","l\nI", general_pth),
         general_pth = gsub("r P","r\nP", general_pth),
         general_pth = gsub("n D","n\nD", general_pth),
         general_pth = gsub("Icm ","Icm\n", general_pth),
         general_pth = gsub("Genetic Information P","Genetic Info.\nP", general_pth),
         nclust2 = ifelse(variable == "clust2", n, NA),
         nclust1 = ifelse(variable == "clust1", n, NA),
         variable = ifelse(variable == "clust1","Group 1", "Group 2"))

bar_ord = general_counts %>% select(general_pth, prop) %>% arrange(desc(prop)) %>% pull(general_pth) %>% unique()

png(file.path(resdir, "clust_bargraphs.png"), width = 1600, height = 800, res = 200)
ggplot(general_counts) + 
  geom_hline(yintercept = 0, linetype = "solid")+ 
  geom_bar(stat='identity', aes(x=factor(general_pth, levels = bar_ord),
                                y=prop, fill = variable), width=.5, colour = "black")  +
  geom_text(mapping = (aes(x=factor(general_pth, levels = bar_ord), y = prop, label = nclust1)),
            na.rm = T, check_overlap = T, size = 5, nudge_x = 0, nudge_y = -1, hjust = 1) +
  geom_text(mapping = (aes(x=factor(general_pth, levels = bar_ord), y = prop, label = nclust2)),
            na.rm = T, check_overlap = T, size = 5, nudge_x = 0, nudge_y = 1, hjust = 0) +
  
  scale_fill_manual(values = c("Group 1"="#8B3399",
                               "Group 2"="#ADE25D"),
                    limits = c("Group 1", "Group 2")) +
  scale_y_continuous(limits = c(-100,100), expand = c(0,0), breaks = seq(-100, 100, 25), 
                     labels = c(seq(100,0,-25),seq(25,100,25)))+ # coord_cartesian(clip = 'off') +
  #facet_grid(factor(KEGG_pthwy_general_grp)~ ., scales = "free", space = "free_y") +
  coord_flip() + theme_bw() + ylab("Percentage of genes")+  xlab("Pathway")+
  theme(#legend.title = element_blank(),
        legend.position = "none", legend.title = element_blank(), legend.direction = "vertical",
        panel.border = element_blank(), axis.line = element_line(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 15), plot.margin=unit(c(0,.03,0,0),"npc"),
        panel.grid = element_blank())
dev.off()

# summary statements 
print(paste("There are", length(clust1_loci), "genes in cluster 1 (left cluster)"), quote = F)
print(paste("There are", length(clust2_loci), "genes in cluster 2 (left cluster)"), quote = F)

clust1_mean = allgenes_folds %>% filter(rownames(.) %in% clust1_loci) %>% select(-c("As", "Tm", "Os")) %>% as.matrix() %>% mean()
print(paste("The avergae log2 fold change in cluster 1 is:",round(clust1_mean,digits = 3)), quote = F)

clust2_mean = allgenes_folds %>% filter(rownames(.) %in% clust2_loci) %>% select(-c("As", "Tm", "Os")) %>% as.matrix() %>% mean()
print(paste("The avergae log2 fold change in cluster 2 is:",round(clust2_mean,digits = 3)), quote = F); 



#  7. clustering box ----
resdir = file.path(resdir_parent,"clustering_may"); if(!file.exists(resdir)) dir.create(resdir)

corinput1 = allgenes_folds %>% filter(rownames(.) %in% clust1_loci) %>% select(-c("As","Os","Tm")) %>% as.matrix() %>% as.numeric()
corinput2 = allgenes_folds %>% filter(rownames(.) %in% clust2_loci) %>% select(-c("As","Os","Tm")) %>% as.matrix() %>% as.numeric()
corinput_pre = data.frame(clust1 = corinput1) %>% 
  reshape2::melt()
corinput_pre2 = data.frame(clust2 = corinput2) %>% 
  reshape2::melt() %>% 
  rbind(corinput_pre) %>% 
  mutate(diffex = ifelse(value <= -log2(1.5), "dn", "non-DE"),
         diffex = ifelse(value >= log2(1.5), "up", diffex))
corinput = data.frame(mycol) %>% rownames_to_column("diffex") %>% filter(diffex %in% c("up","dn","non-DE")) %>% 
  right_join(corinput_pre2, by = "diffex") %>% 
  mutate(x = ifelse(variable == "clust1",1,2))

png(file.path(resdir, "clust_boxplots.png"), width = 800, height = 800, res = 200)
ggplot() +
  geom_segment(aes(x = 0.33, y =  0.5,xend = 2.77, yend=  0.5), linetype = "dashed", colour = "grey48", lwd = 1) + 
  geom_segment(aes(x = 0.33, y = -0.5,xend = 2.77, yend= -0.5), linetype = "dashed", colour = "grey48", lwd = 1) + 
  ggpol::geom_boxjitter(filter(corinput, variable == "clust1"), 
                        mapping =  aes(x = x, y = value, group = variable),# position=position_dodge(0),
                        outlier.shape = NA, jitter.shape = 21, width = .9,jitter.fill = corinput[corinput$variable == "clust1", "mycol"],
                        lwd = .8, fatten = 1, errorbar.draw = T, errorbar.length = 0.33, fill = "grey87", jitter.size = 2.5,
                        show.legend = as.logical(corinput[corinput$variable == "clust1", "mycol"] )) + 
  ggpol::geom_boxjitter(filter(corinput, variable == "clust2"), 
                        mapping =  aes(x = x, y = value, group = variable),
                        outlier.shape = NA, jitter.shape = 21, width = .9,jitter.fill = corinput[corinput$variable == "clust2", "mycol"],
                        lwd = .8, fatten = 1, errorbar.draw = T, errorbar.length = 0.33, fill = "grey87",  jitter.size = 2.5) + 
  xlab("") + ylab(expression(log[2](FC)))+
  scale_fill_manual(values = mycol, limits = unique(corinput$diffex)) + 
  scale_y_continuous(limits = c(-8,12), expand = c(0,0),breaks = c(seq(-8,12,2))) + 
  scale_x_continuous(breaks = c(1,2), labels = c("Group 1", "Group 2"), expand = c(0,0), limits = c(0.33,2.77)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), axis.line = element_line(), legend.position = "top", 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 15), plot.margin=unit(c(0.03,0,0,0),"npc"),
        axis.text.y = element_text(vjust = .25), axis.title.y = element_text(hjust = 0.55, size = 15) )
dev.off()

corinput %>% 
  count(variable, diffex) %>% 
  arrange(desc(variable), desc(diffex))

t.test(corinput1, corinput2)

#  8. Correlation analysis ----

resdir = file.path(resdir_parent, "corplots"); if(!file.exists(resdir)) dir.create(resdir)

corconds = names(allgenes_folds)[names(allgenes_folds) %!in% c("As", "Tm", "Os")]
completedcombos = c()
for(cond in corconds){
  print(cond)
  for(cond2 in corconds[corconds != cond]){
    combo = sort(c(cond,cond2))
    combo = paste0(combo[1],combo[2])
    if(combo %!in% completedcombos){
      pdf(file = file.path(resdir, paste0(cond,"_vs_",cond2,".pdf")), height = 7, width  = 7)
      myplot = ggplot(allgenes_folds, aes(x = !!as.name(cond), y = !!as.name(cond2))) +
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7)
      print(myplot)
      dev.off()
      completedcombos = c(completedcombos,  combo)
    }
  }
}


# plots for report 

png(file = file.path(resdir, paste0("Sp_vs_Nd.png")), height = 800, width  = 800, res = 200)
ggplot(allgenes_folds, aes(x = Sp, y = Nd)) +
  geom_point(colour = "#7D98A1") +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7) +
  scale_y_continuous(breaks = c(-6,-4,-2,0,2,4,6,8))+
  scale_x_continuous(breaks = c(-6,-4,-2,0,2,4,6,8,10))+
  theme_bw() + ylab("log2(FC) in nutrional downshift") + xlab("log2(FC) in stationary phase")+
  theme(axis.line = element_line(), panel.border = element_blank(), panel.grid = element_blank())
dev.off()

png(file = file.path(resdir, paste0("Ns_vs_Nd.png")), height = 1600, width  = 1600, res = 400)
ggplot(allgenes_folds, aes(x = Ns, y = Nd)) +
  geom_point(colour = "#0D1F2D") + 
  # geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7) + 
  scale_y_continuous(breaks = c(-6,-4,-2,0,2,4,6,8))+
  scale_x_continuous(breaks = c(-6,-4,-2,0,2,4,6,8,10))+
  theme_bw() + xlab("log2(FC) in nitrosative stress ") + ylab("")+
  theme(axis.line.x = element_line(), panel.border = element_blank(), panel.grid = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
dev.off()



corconds = names(allgenes_folds)[names(allgenes_folds) %!in% c("As", "Tm", "Os")]

completedcombos = c()
cor.output <- matrix(NA,nrow=length(corconds),ncol=length(corconds),dimnames=list(corconds, corconds))
p.output <-  matrix(NA,nrow=length(corconds),ncol=length(corconds),dimnames=list(corconds, corconds))

for(cond in corconds){
  print(cond)
  for(cond2 in rev(corconds[corconds != cond])){
    combo = sort(c(cond,cond2))
    combo = paste0(combo[1],combo[2])
    if(combo %!in% completedcombos){
      cortest = cor.test(allgenes_folds[,cond], allgenes_folds[,cond2])
      cor.output[cond, cond2] = cortest$estimate
      p.output[cond, cond2] = cortest$p.value
      
      completedcombos = c(completedcombos,  combo)
    }
  }
}

p.output = p.output[-length(corconds),-1]
cor.output = cor.output[-length(corconds),-1]


write.table(cor.output, file.path(resdir, "cor.output.csv"), sep = ",", row.names = T, col.names = T)
write.table(p.output, file.path(resdir, "p.output.csv"), sep = ",", row.names = T, col.names = T, quote = F)


#  9. Dot/Icm effectors ----

resdir = file.path(resdir_parent, "Dot.Icm"); if(!file.exists(resdir)) dir.create(resdir)

eff_hm = allgenes_folds %>%
  select(-c("As","Os","Tm")) %>%
  rownames_to_column(var = "locus_tag") %>%
  filter(locus_tag %in% effectors_loci) %>%
  # left_join(sigma_loci) %>%
  column_to_rownames(("locus_tag")) %>%
  as.matrix()

col_fun2 = colorRamp2(c(-2,-0.4, 0,0.4, 2), c("blue3","white", "white","white", "firebrick3"))
Heatmap(eff_hm, col = col_fun2, show_row_names = F, show_column_names = T, name = "log2(FC)",
        clustering_distance_columns = cosine_distance, clustering_distance_rows = "euclidean",
        # clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        row_split = 2, #column_split = 2
)
decorate_heatmap_body("log2(FC)", row_slice = c(1), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))}) # plot thicker borders
decorate_heatmap_body("log2(FC)", row_slice = c(2), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})

# rotate heatmap
eff_hm_t = eff_hm %>% 
  t()

eff_hm_plot = Heatmap(eff_hm_t, col = col_fun2, show_row_names = T, show_column_names = F, name = "log2(FC)",
                      clustering_distance_rows = cosine_distance, clustering_distance_columns  = "euclidean",
                      clustering_method_rows = "ward.D", clustering_method_columns =  "ward.D",
                      column_split = 2,  row_title = "Stress Condition",
                      column_title = "Dot/Icm Effectors",#column_split = 2
                      column_dend_height = unit(0.15, "npc"), row_dend_width = unit(0.05, "npc"), # control heights of dendrograms relative to heatmap
                      #row_dend_width = unit(0.1, "npc"), row_split = 2,
                      column_gap = unit(0.025, "npc"),
)

png(file.path(resdir, "effectors_hm.png"), width = 1600, height = 700, res = 200)
eff_hm_plot
decorate_heatmap_body("log2(FC)", row_slice = c(1), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))}) # plot thicker borders
decorate_heatmap_body("log2(FC)", column_slice = c(2), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()

efford = column_order(eff_hm_plot)
clust1_eff = colnames(eff_hm_t)[efford[[1]]]; length(clust1_eff)
clust2_eff = colnames(eff_hm_t)[efford[[2]]]; length(clust2_eff)

mean(eff_hm_t[,colnames(eff_hm_t) %in% clust1_eff])
mean(eff_hm_t[,colnames(eff_hm_t) %in% clust2_eff])
DE_genes %>% filter(locus_tag %in% clust2_eff) %>% filter(Nd_change == "upregulated" & Sp_change == "upregulated") %>%
  filter(Nd_log2FC >= 2 & Sp_log2FC >= 2) %>% nrow()

DE_KEGG %>% 
  filter(locus_tag %in% clust2_eff) %>% View() 

eff_top10 = eff_hm_t %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("locus_tag") %>% 
  filter(locus_tag %in% clust2_eff) %>% 
  select(-c("Ns")) %>% 
  group_by(locus_tag) %>% 
  summarise(mean = mean(c(Li, Mi, Nd, Ox, Sp))) %>% 
  ungroup() %>% 
  arrange(desc(mean)) %>% 
  mutate(order = c(1:131))
 
DE_KEGG %>% 
  filter(locus_tag %in% eff_top10) %>% View()

fasta_con = NA
for(eff in eff_top10){
  print(eff)
  tmp_fasta = read.table(file.path(datadir,"BLAST_results_cl",paste0(eff,".fasta_result")),sep = "\t")
  if(all(is.na(fasta_con))){
    fasta_con = tmp_fasta
  }else{
    fasta_con = rbind(fasta_con, tmp_fasta)
  }
}
colnames(fasta_con) = unlist(strsplit("qseqid sseqid sacc ssciname stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"," "))



eff_countdf = eff_hm_t
eff_countdf[,305] = NA

effchanges = allgenes_changes[rownames(allgenes_changes) %in% effectors_loci, c(paste0(c("Li","Mi","Nd","Ox","Sp"),"_change"))]
effchanges[,"n_up"] = NA

for(i in 1:nrow(effchanges)){
  effchanges[i,"n_up"] = length(which(effchanges[i,] == "upregulated"))
}


# 10. Dot/Icm machinery ----
resdir = file.path(resdir_parent, "Dot.Icm"); if(!file.exists(resdir)) dir.create(resdir)

doticm_df = DE_genes %>% 
  filter(locus_tag %in% T4SS_loci) %>% 
  select(locus_tag, all_of(changecols)) %>%
  select(-c("As_change", "Os_change","Tm_change")) %>% 
  column_to_rownames("locus_tag") %>% 
  t() %>% as.data.frame()



doticm_df$n_up = NA; doticm_df$n_dn = NA
for(i in 1:nrow(doticm_df)){
  doticm_df[i,"n_up"] = length(which(doticm_df[i,] == "upregulated"))
  doticm_df[i,"n_dn"] = length(which(doticm_df[i,] == "downregulated"))
}

doticm_input = doticm_df %>% rownames_to_column("locus_tag") %>% 
  reshape2::melt(value.name = "n_changed") %>% arrange(locus_tag) %>% 
  mutate(prop = n_changed / length(T4SS_loci) * 100,
         prop = ifelse(variable == "n_dn", - prop, prop),
         n_changed = ifelse(variable == "n_dn", - n_changed, n_changed),
         cond = substr(locus_tag,1,2),
         variable = ifelse(variable == "n_up","upregulated","downregulated"))

png(file.path(resdir, "machinery_bars_top.png"), width = 800, height = 500, res = 200)
ggplot(doticm_input) + 
  geom_hline(yintercept = 0, linetype = "solid")+ 
  geom_bar(stat='identity', aes(x=factor(cond, levels = sort(unique(cond), decreasing = F)),
                                y=n_changed, fill = variable), width=.5, colour = "black")  +
  # geom_text(mapping = (aes(x=factor(general_pth, levels = bar_ord), y = prop, label = nclust1)),
  #           na.rm = T, check_overlap = T, size = 5, nudge_x = 0, nudge_y = -1, hjust = 1) +
  # geom_text(mapping = (aes(x=factor(general_pth, levels = bar_ord), y = prop, label = nclust2)),
  #           na.rm = T, check_overlap = T, size = 5, nudge_x = 0, nudge_y = 1, hjust = 0) +
  scale_fill_manual(values = mycol, limits = c("upregulated", "downregulated")) +
  scale_y_continuous(limits = c(-9,18), expand = c(0,0), breaks = seq(-6, 27, 6), 
                     #labels = seq(-9,27,9)
  )+ # coord_cartesian(clip = 'off') +
  #facet_grid(factor(KEGG_pthwy_general_grp)~ ., scales = "free", space = "free_y") +
  # coord_flip() + 
  theme_bw() + ylab("No. of Dot/Icm machinery genes") + xlab("Condition")+
  theme(legend.title = element_blank(),
        legend.position = "top", 
        panel.border = element_blank(), axis.line = element_line(), axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), #plot.margin=unit(c(0,.03,0,0),"npc"),
        panel.grid = element_blank())
dev.off()

# 11. Differentiation mediators ----

resdir = file.path(resdir_parent,"differentiation_control"); if(!file.exists(resdir)) dir.create(resdir)

# sigma factors 
sigma_loci = c("lpg0477", "lpg1284","lpg1577", "lpg1762", "lpg1782", "lpg2361", "lpg2667")
names(sigma_loci) = c("sigma 54", "rpoS", "rpoE", "fleR", "fliA", "rpoD", "rpoH")
sigma_loci = data.frame(sigma_loci) %>% rownames_to_column(var = "gene") %>% rename(locus_tag = sigma_loci)


sigma_loci = c("lpg0477", "lpg1284","lpg1577", "lpg0853", "lpg2361", "lpg2667", "lpg1782")
names(sigma_loci) = c("RpoN (\u3C3 54)", "RpoS (\u3C3)", "RpoE (\u3C3)", "FleQ (\u3C3)", "RpoD (\u3C3 70)", "RpoH (\u3C3 32)", "FliA (\u3C3 28)")
sigma_loci = data.frame(sigma_loci) %>% rownames_to_column(var = "gene") %>% rename(locus_tag = sigma_loci)

sigcol_list = c("darkgreen", "deepskyblue4", "grey", "orangered1", "darkred", "goldenrod1", "deeppink4"); names(sigcol_list) =  unique(sigma_input$gene)
sigcol = data.frame(sigcol_list) %>% rownames_to_column(var = "gene")

sigma_input = allgenes_folds %>% 
  select(-c("As","Os","Tm")) %>% 
  rownames_to_column(var = "locus_tag") %>% 
  filter(locus_tag %in% sigma_loci$locus_tag) %>% 
  reshape2::melt() %>% 
  left_join(sigma_loci) %>% 
  left_join(sigcol) %>% mutate(sigcol_list = as.character(sigcol_list), gene = factor(gene))




# diff_loci = c("lpg2009", "lpg1457", "lpg2094", "lpg0537", "lpg1437", "lpg1438", "lpg2646", "lpg1912")
# names(diff_loci) = c("spoT", "relA", "crsA", "LetE", "cpxA",  "cpxR", "letA", "letS")
diff_loci = c("lpg2009", "lpg1457", "lpg2094", "lpg0537", "lpg1437", "lpg1438", "lpg2646", "lpg1912", "lpg2731","lpg2732","lpg2733","lpg2734", "lpg1762")
names(diff_loci) = c("SpoT", "RelA", "csrA", "LetE", "cpxA",  "cpxR", "LetA", "LetS", "LqsA", 
                     "LqsR" , "HdeD" , "LqsS", "FleR")

# lqsR_loci = c("lpg2731","lpg2732","lpg2733","lpg2734"); names(lqsR_loci) = c("lqsA", "lqsR" , "hdeD" , "lqsS")
 # = as.character(c(diff_loci, lqsR_loci))
diff_loci = data.frame(diff_loci) %>% rownames_to_column(var = "gene") %>% rename(locus_tag = diff_loci)


diff_input = allgenes_folds %>% 
  select(-c("As","Os","Tm", "Ns")) %>% 
  rownames_to_column(var = "locus_tag") %>% 
  filter(locus_tag %in% diff_loci$locus_tag) %>% 
  reshape2::melt() %>% 
  left_join(diff_loci)
# library(jcolors)


hm_genes = rbind(sigma_loci, diff_loci)

sigma_hm = allgenes_folds %>% 
  select(-c("As","Os","Tm")) %>%
  rownames_to_column(var = "locus_tag") %>% 
  filter(locus_tag %in% hm_genes$locus_tag) %>% 
  left_join(hm_genes) %>% 
  mutate(gene = gsub("c","C", gene)) %>% 
  column_to_rownames(var = "gene") %>% select(-c(locus_tag)) %>% 
  as.matrix() %>% t()

png(file = file.path(resdir, "corr_heatmap_t.png"), width = 2000, height = 800, res = 200)
col_fun2 = colorRamp2(c(-2,-0.4, 0,0.4, 2), c("blue3","white", "white","white", "firebrick3"))
Heatmap(sigma_hm, col = col_fun2, show_row_names = T, show_column_names = T, name = "log2(FC)",
        heatmap_width = unit(.8, "npc"), column_dend_height = unit(.15,"npc"),#row_dend_width = unit(.2, "npc"),
        heatmap_height = unit(.5,"npc"), column_title = "Life cycle regulatory genes", row_title = "Stress Condition",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(sigma_hm[i, j] < -1) grid.text(sprintf("%.1f", sigma_hm[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
          if(sigma_hm[i, j] <= -log2(1.5) & sigma_hm[i, j] > -1) grid.text(sprintf("%.1f", sigma_hm[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
          if(sigma_hm[i, j] >= log2(1.5) & sigma_hm[i, j] < 1.7) grid.text(sprintf("%.1f", sigma_hm[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
          if(sigma_hm[i, j] >= 1.7) grid.text(sprintf("%.1f", sigma_hm[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
          # if(sigma_hm[i, j] < log2(1.5) & sigma_hm[i, j] > -log2(1.5)) grid.text(sprintf("%.1f", sigma_hm[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
        },
        heatmap_legend_param = list(
          direction = "horizontal" ,
          title_position = "topcenter"
        ),
        rect_gp = gpar(col = "black", lwd = 1)
        )
decorate_heatmap_body("log2(FC)", column_slice = c(1), {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))}) # plot thicker borders
dev.off()

Legend(col_fun = col_fun, title = "foo", direction = "horizontal")



# diff correlation 

flag_fc = DE_KEGG %>% 
  filter(KEGG_pthwy_grp == "Cell_motility", gene != "fliA") %>% 
  select(locus_tag, all_of(foldchanges)) %>% 
  distinct() %>% 
  column_to_rownames("locus_tag")

flagmeans = matrix(NA, nrow = 2, ncol = 9, dimnames = list(c("mean_flag","FliA"), c(paste0(unique(substr(cond_colnames,1,2))))))
for(i in 1:ncol(flag_fc)){
  flagmeans[1,i] = mean(flag_fc[,i]) # add mean of flag genes
  cond = substr(colnames(flag_fc)[i],1,2)
  flagmeans[2,i] =  DE_genes[DE_genes$gene == "fliA", paste0(cond, "_log2FC")] # add FliA FC
}

flagmeans = flagmeans %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Stress.condition") 

hded_fc = DE_KEGG %>% 
  # filter(locus_tag %in% effectors_loci) %>%  
  filter(locus_tag %in% clust2_eff) %>%  
  select(locus_tag, all_of(foldchanges)) %>% 
  distinct() %>% 
  column_to_rownames("locus_tag")

hdedmeans = matrix(NA, nrow = 2, ncol = 9, dimnames = list(c("mean_eff","HdeD"), c(paste0(unique(substr(cond_colnames,1,2))))))
for(i in 1:ncol(hded_fc)){
  hdedmeans[1,i] = mean(hded_fc[,i]) # add mean of hded genes
  cond = substr(colnames(hded_fc)[i],1,2)
  hdedmeans[2,i] =  DE_genes[DE_genes$locus_tag == "lpg2733", paste0(cond, "_log2FC")] # add FliA FC
}

hdedmeans = hdedmeans %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Stress.condition") 
bothmeans = hdedmeans %>% left_join(flagmeans)
# bothmeansfilt = filter(bothmeans, Stress.condition %!in% c("As","Os","Tm"))
bothmeansfilt = bothmeans

condcolours = c("#FAFA00","#ff5e5b","#5438dc","#403233","#70D6FF","#5fad56","#f2c14e","#A4D747","#b4436c") # https://coolors.co/fafa00-df9a57-5438dc-403233-70d6ff-5fad56-f2c14e-a4d747-b4436c

# effectors vs. HdeD
png(file = file.path(resdir, paste0("Dot.Icm_vs_HdED.png")), height = 1600, width  = 1600, res = 400)
# png(file = file.path(resdir, paste0("Dot.Icm_vs_HdED_legend.png")), height = 1600, width  = 2000, res = 400)
ggplot(bothmeansfilt, aes(y = mean_eff, x = HdeD, colour = Stress.condition)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7) +
  geom_point(size = 5) +
  scale_colour_manual(values = condcolours) +
  # scale_y_continuous(breaks = c(-.25,0,.25,.5,.75,1, 1.25), limits = c(-.25,1.25))+
  # scale_x_continuous(breaks = c(0,.5,1, 1.5,2,2.5), limits = c(-.1,2.5))+
  theme_bw() + 
  xlab(expression(HdeD~log[2](FC)))+ ylab(expression(mean~Dot/Icm~effector~log[2](FC)))+
  labs(colour="Condition")+
  theme(axis.line = element_line(), panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.direction = "horizontal", legend.justification = c(0,0), legend.text.align = 0, legend.text = element_text(size = 15))
dev.off()
cor.test(bothmeansfilt$mean_eff, bothmeansfilt$HdeD, method = "pearson")

# flagella vs. FliA
png(file = file.path(resdir, paste0("motility_vs_FliA.png")), height = 1600, width  = 1600, res = 400)
ggplot(bothmeansfilt, aes(y = mean_flag, x = FliA, colour = Stress.condition)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7) +
  geom_point(size = 5) +
  # scale_colour_manual(values = c("forestgreen", "darkturquoise", "grey", "orangered1", "darkolivegreen3",
  #                                "goldenrod1", "magenta4","mediumblue", "maroon1")) +
  scale_colour_manual(values  = condcolours) + 
  scale_y_continuous(breaks = c(-.5,0,.5,1,1.5,2), limits = c(-.5,2))+
  scale_x_continuous(breaks = c(-.5,0,.5,1,1.5,2,2.5), limits = c(-.6,2.5))+
  theme_bw() +# ylab("mean Cell motility log2(FC)") + xlab("FliA log2(FC)")+
  xlab(expression(FliA~log[2](FC)))+ ylab(expression(mean~cell~motility~log[2](FC)))+
  theme(axis.line = element_line(), panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()
cor.test(bothmeansfilt$mean_flag, bothmeansfilt$FliA, method = "pearson")


cor.test(bothmeansfilt$mean_flag, bothmeansfilt$FliA, method = "pearson")
cor.test(bothmeansfilt$mean_flag, bothmeansfilt$HdeD, method = "pearson")
cor.test(bothmeansfilt$mean_eff, bothmeansfilt$FliA, method = "pearson")


# metabolism genes vs RpoD
met_fc = DE_KEGG %>% 
  filter(KEGG_pthwy_general_grp == "Metabolism") %>% 
  select(locus_tag, all_of(foldchanges)) %>% 
  distinct() %>% 
  column_to_rownames("locus_tag")

metmeans = matrix(NA, nrow = 2, ncol = 9, dimnames = list(c("mean_met","RpoD"), c(paste0(unique(substr(cond_colnames,1,2))))))
for(i in 1:ncol(met_fc)){
  metmeans[1,i] = mean(met_fc[,i]) # add mean of met genes
  cond = substr(colnames(met_fc)[i],1,2)
  metmeans[2,i] =  DE_genes[DE_genes$locus_tag == "lpg2361", paste0(cond, "_log2FC")] # add RpoD
}

metmeans = metmeans %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Stress.condition") #%>% 
  # filter(Stress.condition %!in% c("As","Os","Tm"))

png(file = file.path(resdir, paste0("metabolism_vs_RpoD.png")), height = 1600, width  = 1600, res = 400)
ggplot(metmeans, aes(y = mean_met, x = RpoD, colour = Stress.condition)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = "dashed", size = 0.7) +
  geom_point(size = 5) +
  scale_colour_manual(values = condcolours) +
  scale_x_continuous(breaks = seq(-2,1,.5), limits = c(-2,0.7))+
  theme_bw() +
  xlab(expression(RpoD~log[2](FC)))+ ylab(expression(mean~metabolism~log[2](FC)))+ labs("Condition") + 
  theme(axis.line = element_line(), panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()

cor.test(metmeans$mean_met, metmeans$RpoD, method = "pearson")