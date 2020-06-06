# Final year project functions source file March 2020

# 1. not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# 2. load file to named variable
load2object <- function (filename) 
{
  if (file.exists(filename)) 
    return(eval(parse(text = load(filename))))
  cat(paste("error! the file ", filename, 
            " was not found :("))
  NULL
}

# 3. remove cond_colnames from dataframe
condselect = function (x, rm_colnames = cond_colnames) 
{
  dplyr::select(x, -c(rm_colnames))
}

# 4. remove cond_colnames from dataframe and view dataframe
viewnocond = function (x, rm_colnames = cond_colnames) 
{
  tmp = dplyr::select(x, -c(all_of(rm_colnames)))
  View(tmp)
}

# 5. extract fold chnages input


plot_nchanges_input = function(conditions = conditions_list, paths = NA, input_df = NA, input_path_col = "KEGG_set"){
  
  newlist = c()
  
  if(typeof(paths) == "list"){
    paths = unique(input_df[,input_path_col][!is.na(input_df[,input_path_col])]) %>% sort() 
  }
  
  newlist = input_df %>% # create counts of total no. of genes in each pathway 
    dplyr::select(locus_tag, !!as.name(input_path_col)) %>% 
    rename(pathway  := !!as.name(input_path_col)) %>% 
    distinct() %>% 
    count(pathway) %>% 
    as.data.frame() %>% 
    filter(!is.na(pathway))
  
  for(condition in conditions){
    print(condition)
    fold_colname = paste0(condition,fold_suffix)
    p_colname = paste0(condition,p_suffix)
    reg_colname = paste0(condition,"_change")
    
    res = matrix(nrow = length(paths), ncol = 3, dimnames = list(paths, c("up","dn","non-DE")) )
    
    for(chemin in paths){ # count no. of genes in each direction in each condition & pathway 
      df_filt = input_df %>% 
        dplyr::select(locus_tag, !!as.name(input_path_col), !!as.name(reg_colname)) %>% 
        distinct() %>% 
        filter(!!as.name(input_path_col) == chemin) %>% 
        pull(!!as.name(reg_colname)) 
      
      res[chemin,"up"] = length(df_filt[df_filt == "upregulated"])
      res[chemin,"dn"] = length(df_filt[df_filt == "downregulated"])
      res[chemin,"non-DE"] = length(df_filt[df_filt == "non-DE"])
    }
    
    resdf = as.data.frame(res) %>% 
      rownames_to_column(var = "pathway") %>% 
      reshape2::melt() %>% 
      arrange(pathway, variable) %>% 
      left_join(newlist, by = "pathway") %>% 
      mutate(prop = value / n,
             cond =  condition)
    
    if(condition == conditions[1]) res_final = resdf
    if(condition != conditions[1]) res_final = rbind(res_final, resdf)
  }
  return(res_final)
}

# 6. plot fold number of genes up or down 

plot_nchanges = function(input = NULL, resdir = NULL, conditions_list = conditions_list,
                            y_colname = "value", x_colname = "pathway", fill_colname = "variable",
                            proportion = TRUE, 
                            counts = NULL, filename = NULL){
  if(!dir.exists(resdir)) dir.create(resdir)
  for(condition in conditions_list){
    print(condition)
    
    input_filt = input %>% 
      filter(cond == condition)
    
    n_xvar = length(unique(input_filt[,x_colname]))
    plot_width = ifelse(n_xvar >= 20, 15, 7)
    
    pdf(file = file.path(resdir,paste0(condition,filename,".pdf")), width = plot_width, height = 7)
    myplot = ggplot() +
      geom_bar(mapping = aes(x = input_filt[,x_colname], y =  input_filt[,y_colname], fill = input_filt[,fill_colname]),
               stat = "identity", position = "dodge") +
      ggtitle(paste("stress condition:",condition)) + 
      scale_fill_manual(values = mycol) +
      xlab("Pathway") + ylab("Counts") + 
      scale_y_continuous(expand = expand_scale(mult = c(0,.1), add = 0)) +
      theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5), legend.position = "none")
    print(myplot)
    
    if(proportion == TRUE){
      input_filt = input_filt %>% 
        filter(variable != "non-DE")
      myprop = ggplot() +
        geom_bar(mapping = aes(x = input_filt[,x_colname], y =  input_filt[,"prop"], fill = input_filt[,fill_colname]),
                 stat = "identity", position = "dodge") +
        scale_fill_manual(values = mycol) +
        ggtitle(paste("stress condition:",condition)) + 
        xlab("Pathway") + ylab("Proportion") + 
        scale_y_continuous(expand = expand_scale(mult = c(0,0), add = 0), limits = c(0,1)) +
        theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5), legend.position = "none")
      print(myprop)
    }
    dev.off()
  }
}

# 7. plot fold changes of genes in selection 
plot_foldchange = function(input = NULL, gene_col = NULL, fill_col = NULL, fold_col = NULL,
                           fold_filt = fold_filter, condit = condition){
  
  input = input %>% 
    mutate(gene_col = !!as.name(gene_col),
           fold_col = !!as.name(fold_col),
           fill_col = !!as.name(fill_col))
  
  max.y = max(input$fold_col)*1.1; min.y = min(input$fold_col)*1.1
  dottedlen = nrow(input) + 0.5
  
  myplot  = ggplot() +
    geom_bar(mapping = aes(x = factor(input$gene_col,levels = sort(input$gene_col)), y =  input$fold_col, 
                           fill = input$fill_col),
             stat = "identity", position = "dodge", width = 0.75) +
    ggtitle(paste("stress condition:",condit)) + 
    scale_y_continuous(expand = c(0,0), limits = c(min.y, max.y)) + 
    scale_fill_manual(values = mycol)  +
    geom_segment(aes(x=0.5,xend=dottedlen,y=fold_filt,yend=fold_filt),colour = "black", linetype = 2, size  = 0.5)+
    geom_segment(aes(x=0.5,xend=dottedlen,y=-fold_filt,yend=-fold_filt),colour = "black", linetype = 2, size  = 0.5)+
    theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5), legend.position = "none")
  
  print(myplot)
}




# # run BLAST in R 
# blastSeqKK <- function (x, database = "nr", hitListSize = "10", 
#                         filter = "L", expect = "10", program = "blastn",
#                         attempts = 10) {
#   baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
#   query <- paste("QUERY=", as.character(x), "&DATABASE=", database, 
#                  "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
#                  expect, "&PROGRAM=", program, sep = "")
#   url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
#   results <- tempfile()
#   Sys.sleep(5)
#   require(XML)
#   post <- htmlTreeParse(url0, useInternalNodes = TRUE)
#   x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
#   rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
#   rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
#                          x))
#   url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
#                   rid)
#   Sys.sleep(rtoe)
#   .tryParseResult <- function(url, attempts){
#     for (i in 1:(attempts+1)) {
#       result <- tryCatch({
#         xmlTreeParse(url, useInternalNodes=TRUE,
#                      error = xmlErrorCumulator(immediate=FALSE))
#       }, error=function(err) NULL)
#       if (!is.null(result)) return(result)
#       Sys.sleep(10)
#     }
#     stop(paste("no results after ", attempts, 
#                " attempts; please try again later", sep = ""))
#   }
#   result <- .tryParseResult(url1, attempts)
#   qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
#   hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
#   require(Biostrings)
#   res <- list()
#   for (i in seq_len(length(qseq))) {
#     res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
#                                    rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
#                                                                                           "NormalIRanges"))
#   }
# }



# cosine distance function, taken from the package Palimpsest: https://github.com/FunGeST/Palimpsest 
# method adapted from Shinde et al. (2018) Bioinformatics & Letouzé et al. (2017) Nat. Comms. 

cosine_distance = function (m) 
{
  requireNamespace("lsa", quietly = T)
  nsamp <- nrow(m)
  res <- matrix(NA, nrow = nsamp, ncol = nsamp)
  rownames(res) <- colnames(res) <- rownames(m)
  for (i in 1:nsamp) {
    for (j in 1:nsamp) {
      res[i, j] <- lsa::cosine(m[i, ], m[j, ])
    }
  }
  as.dist(1 - res)
}


# install necessary packages
install_necessary_packages = function(){
  cran_packages <- c("tidyverse", "gdata", "reshape2", "data.table", "lsa", "BiocManager", "ggpol")
  new.packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  bioconductor_packages <- c("gage", "ComplexHeatmap", "circlize")
  new_bioc_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
  if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages)
}