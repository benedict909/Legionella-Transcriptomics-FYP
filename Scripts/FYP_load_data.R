# Final year project data load & process February 2020
# Load DE data and add KEGG annotations

library(tidyverse)
library(gdata)
library(dplyr)
library(reshape2)
#library(KEGGREST)
library(gage)
library(gdata)


# Load data ---------------------------------------------------------------------------------------
setwd("~")
datadir = "../Dropbox/Final_Year/Project/Data/"
datafile = file.path(datadir,"Legionella_DE.xlsx")
resdir_parent = "../Dropbox/Final_Year/Project/Results/"

source('C:/Users/bened/Dropbox/Final_Year/Project/Scripts/FYP_DE_functions.R')

# DE_data_load = read.xls(xls = datafile, sheet = 1)
# save(DE_data_load, file = file.path(datadir,"Legionella_DE_unfiltered.RData"))

load(file.path(datadir,"Legionella_DE_rm_na.RData"))



# Prepare data ------------------------------------------------------------------------------------

DE_cols = colnames(DE_data) 
cond_colnames = DE_cols[grepl("vs", DE_cols)]
conditions_list = c()
for(n4me in cond_colnames){
  conditions_list[n4me] = unlist(strsplit(n4me, "[.]"))[1]
}
conditions_list = unique(unlist(conditions_list))
# rm_tags = c("Bs|Vic")  ## original filtering
# rm_cols = DE_cols[grep(rm_tags,DE_cols)]
# 
# DE_data = DE_data_load %>% 
#  dplyr::select(-c(rm_cols))
# 
# save(DE_data, file = file.path(datadir,"Legionella_DE_rm_na.RData"))

mycol = c("red","red","blue","blue","blue", "grey45"); names(mycol) = c("up", "upregulated", "downregulated","dn", "down", "non-DE")

fold_filter = log2(1.5)
p_filter = 0.05

fold_suffix = ".vs..Ctrl...Log.fold.change"
p_suffix = ".vs..Ctrl...FDR.p.value"

# Histograms --------------------------------------------------------------------------------------

# resdir = file.path(resdir_parent, "Histograms"); if(!dir.exists(resdir)) dir.create(resdir)
# 
# for(i in 1:length(conditions_list)){
#   condition = conditions_list[i]; print(condition)
#   fold_colname = paste0(condition,fold_suffix)
#   p_colname = paste0(condition,p_suffix)
#   
# 
#   DE_folds = DE_data %>%
#     filter(!!as.name(p_colname) < p_filter) %>% 
#     dplyr::select(fold_colname) %>% 
#     na.omit()
#    
#   
#   DE_foldslist = DE_folds[,1]
#   print(range(DE_foldslist))
#   
#   pdf(file = file.path(resdir,paste0(condition,".pdf")),width = 5, height = 4)
#   myplot = ggplot() +
#               geom_histogram(mapping = aes(x = DE_foldslist), binwidth = .1) 
#   
#   myplot_build <- ggplot_build(myplot)
#   tmp = myplot_build$data[[1]]
#   max.y = max(tmp$ymax)
#   
#   myplot = myplot +
#     geom_segment(aes(x=-fold_filter,xend=-fold_filter,y=0,yend=max.y),
#                colour = "red", linetype = 2, size  = 0.5) + 
#     geom_segment(aes(x=fold_filter,xend=fold_filter,y=0,yend=max.y),
#                  colour = "red", linetype = 2, size  = 0.5) +
#     ggtitle(condition)
#   print(myplot)
#   dev.off()
# }



# Extract DE genes --------------------------------------------------------------------------------

DE_cols_filt = colnames(DE_data)
extract_rm_tags = c("Max.group|...P.Fold.change|Bonferroni|...P.value|...Fold.change")
extract_rm_cols = c("new.locus.tag", as.character(DE_cols_filt[grep(extract_rm_tags, DE_cols_filt)]))

DE_genes = DE_data %>% 
  dplyr::select(-c(all_of(extract_rm_cols))) %>% 
  mutate(locus_tag = as.character(locus_tag))
  

for(i in 1:length(conditions_list)){
  condition = conditions_list[i]; print(condition)
  fold_colname = paste0(condition,fold_suffix)
  p_colname = paste0(condition,p_suffix)
  reg_colname = paste0(substr(condition,1,2),"_change")
  
  DE_genes[,reg_colname] = "non-DE"
  for( j in 1:nrow(DE_genes)){
    if(!is.na(DE_genes[j,fold_colname])){
      if(DE_genes[j,fold_colname] >= fold_filter 
         & DE_data[j,p_colname] <= p_filter) DE_genes[j,reg_colname] = "upregulated"
      if(DE_genes[j,fold_colname] <= -fold_filter
         & DE_data[j,p_colname] <= p_filter) DE_genes[j,reg_colname] = "downregulated"
    }
  }
  print(nrow(filter(DE_genes, !!as.name(reg_colname) == "upregulated")) + nrow(filter(DE_genes, !!as.name(reg_colname) == "downregulated")))
}

change_colnames = colnames(DE_genes)[grepl("_change", colnames(DE_genes))]

non_DE_locus_tags = DE_genes$locus_tag[all(DE_genes[change_colnames] == "non-DE")]  

for(i in 1:nrow(DE_genes)){
  if(all(DE_genes[i,change_colnames] == "non-DE")){
  non_DE_locus_tags[i] = DE_genes$locus_tag[i]
  }
}
non_DE_locus_tags = as.character(non_DE_locus_tags[!is.na(non_DE_locus_tags)])

# DE_genes = DE_genes %>% 
#   filter(locus_tag %!in% non_DE_locus_tags)

core_colnames = colnames(DE_genes)[1:6]
cond_colnames = sort(colnames(DE_genes)[grepl("Ctrl|change",colnames(DE_genes)) & !grepl(extract_rm_tags, colnames(DE_genes))])

DE_genes = DE_genes[,c(core_colnames,cond_colnames)]




# Add gene names from spreadsheets ----------------------------------------------------------------

# KO_names = read.xls(xls = file.path(datadir,"Legionella KO.xlsx"), sheet = 1) %>% 
#   dplyr::select(-c(paste0("kegg",c(12:16))), locus_tag = Locustag) %>% 
#   filter(locus_tag %in% DE_genes$locus_tag)

Patrics_names = read.xls(xls = file.path(datadir,"Patrics for Legionella.xlsx"), sheet = 1) %>% 
  #dplyr::select(KO_description = Description, locus_tag = Locustag, everything()) %>% # use this line to get patrics names too
  dplyr::select(KO_description = Description, locus_tag = Locustag) %>% 
  filter(locus_tag %in% DE_genes$locus_tag)

dim(DE_genes)
  
DE_genes = DE_genes %>% 
    #left_join(y = KO_names, by = "locus_tag") %>% use this line to add KO
    left_join(y = Patrics_names, by = "locus_tag")
dim(DE_genes)
length(unique(DE_genes$locus_tag))


newcolnames = colnames(DE_genes)[colnames(DE_genes) %!in% c(core_colnames, cond_colnames)]
DE_genes = DE_genes[, c(core_colnames, newcolnames, cond_colnames)]


for(cond in unique(substr(colnames(DE_genes),1,2))) {  ## change colnames to make them shorter
  print(cond)
  colnames(DE_genes)[substr(colnames(DE_genes),1,2) == cond & grepl("fold.change",colnames(DE_genes))] = 
     paste0(cond, "_log2FC")
  
  colnames(DE_genes)[substr(colnames(DE_genes),1,2) == cond & grepl("p.value",colnames(DE_genes))] = 
    paste0(cond, "_pvalue")
}

cond_colnames = names(DE_genes)[names(DE_genes) %!in% core_colnames & !grepl("KO", names(DE_genes))]

# mem = DE_genes
# DE_genes = mem 

DE_genes = DE_genes %>% 
  mutate(KO_description = as.character(KO_description)) %>% 
  mutate(KO_description = ifelse(locus_tag == "lpg2369","HipA", KO_description)) %>% # duplivcated gene 
  distinct() %>% # remove duplicate lines 
  filter(!is.na(!!as.name(cond_colnames[!grepl("change", cond_colnames)]))) %>%  #remove genes with NA fold changes
  mutate(locus_tag = gsub("_","", locus_tag)) %>% # remove underscore from locus tags 
  filter(KO_description != "COG2827: putative endonuclease containing a URI domain" | is.na(KO_description)) # duplicated gene

save(DE_genes, file = file.path(datadir,"DE_genes.RData"))


# KEGG add pathways / create DE_KEGG --------------------------------------------------------------

kg.lpn <- kegg.gsets( "lpn" )
kegg.gs2 <- kg.lpn$kg.sets  # take all Kg pathway gene sets 

KEGG_set_table = reshape2::melt(kegg.gs2) %>% 
  dplyr::select(locus_tag = value, KEGG_set = L1) %>% 
  filter(locus_tag %in% DE_genes$locus_tag) %>% 
  mutate_all(as.character)

# KEGG_onerow = data.frame(locus_tag = unique(KEGG_set_table$locus_tag), KEGG1 = NA, KEGG2)
for(loci in unique(KEGG_set_table$locus_tag)){
  tmp = KEGG_set_table %>% 
    filter(locus_tag == loci) %>% 
    pull(KEGG_set)
}

DE_KEGG = DE_genes %>% 
  left_join(KEGG_set_table, by = "locus_tag")

## add pathway "groups" (from "kegg_pathways_lists.R")

load(file.path(datadir,"general_pthwys_list.RData"))
load(file.path(datadir,"pthwys_list.RData"))


DE_KEGG$KEGG_pthwy_grp = NA
for(pthwy in  names(pthwys_list)){
  print(pthwy)
  for(i in 1:nrow(DE_KEGG)){
    if(DE_KEGG$KEGG_set[i] %in% unlist(pthwys_list[pthwy])){
      DE_KEGG$KEGG_pthwy_grp[i] = pthwy
    }
  }
}
head(table(DE_KEGG$KEGG_pthwy_grp))
# Amino_acid_metabolism             Biosynthesis_of_other_secondary_metabolites 
# 206                               19 
# Carbohydrate_metabolism           Cell_motility 
# 192                               44 
# Cellular_community_prokaryotes    Drug_resistance_antimicrobial 
# 33                                31 

DE_KEGG$KEGG_pthwy_general_grp = NA
for(pthwy in  names(general_pthwys_list)){
  print(pthwy)
  for(i in 1:nrow(DE_KEGG)){
    if(DE_KEGG$KEGG_pthwy_grp[i] %in% unlist(general_pthwys_list[pthwy])){
      DE_KEGG$KEGG_pthwy_general_grp[i] = pthwy
    }
  }
}
table(DE_KEGG$KEGG_pthwy_general_grp)
# Cellular_Processes                 Environmental_Information_Processing       Genetic_Information_Processing 
# 77                                 127                                        188 
# Human_Diseases                     Metabolism 
# 56                                 1856 

newcolnames = colnames(DE_KEGG)[colnames(DE_KEGG) %!in% c(core_colnames, cond_colnames)]
DE_KEGG = DE_KEGG[, c(core_colnames, newcolnames, cond_colnames)]


save(DE_KEGG, file = file.path(datadir,"DE_KEGG.RData"))








# Effectors from Burstein et al. 2015

effectors_load = read.xls(xls = file.path(datadir, "Burstein_2015_T4SS_effectors_for_R_load.xlsx"), sheet = 1) 

effectors_df = effectors_load %>% 
  filter(L..pneumophila != "-") %>% 
  filter(Effector == "+")



save(effectors_df, file = file.path(datadir, "T4SS_effectors.RData"))







