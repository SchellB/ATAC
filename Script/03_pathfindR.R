###############################################################################
################ Chromatin Accessibility Pathway Enrichment ###################

# 01 - Loaging data ####

# working directory
wd <- getwd()

# ATACseq normalized counts
library(readr)
library(dplyr)
library(tibble)
ATAC <- read_csv2(paste0(wd,"/Results/Normalized.Counts.csv"))
data <- ATAC %>% column_to_rownames("...1")

# ATACseq DEG loading
load(paste0(wd,"/Script/ATACseqDEG.Rdata"))

# 02 - PathfindR funciton #####

library(pathfindR)


pathenrich <- function(data,table, gene_sets="KEGG",
                       p_val_threshold=0.01,
                       adj_method="BH", wd, name){
  tab <- table %>% select(c(genes, logFC, P.Value))
  
  pdf(file = paste0(wd,"/Output/pathenrich_",gene_sets,"_",name,".pdf"))
  output <- run_pathfindR(tab, gene_sets = gene_sets,
                             p_val_threshold =p_val_threshold,
                          adj_method=adj_method)
  
  clustered <- cluster_enriched_terms(output,
                                      plot_hmap=TRUE,
                                      plot_dend=TRUE)
  
  term_gene_heatmap(output)
  term_gene_graph(output)
  UpSet_plot(output)
  score_terms(output, as.matrix(data))
  dev.off()
  
  write.csv2(output,file=paste0(wd,"/Results/pathenrich_",gene_sets,"_",name,".csv"), row.names = FALSE)
}


# 03 - Performing nerichment ####

## HR vs HD
pathenrich(data=data, table = HRvsHD, name="HRvsHD", gene_sets="KEGG",wd=wd)

## LR vs HD
pathenrich(data=data, table = LRvsHD, name="LRvsHD",gene_sets="KEGG",wd=wd)

## AML vs HD
pathenrich(data=data, table = AMLvsHD, name="AMLvsHD", gene_sets="KEGG",wd=wd)

## LR vs AML
pathenrich(data=data, table = LRvsAML, name="LRvsAML", gene_sets="KEGG",wd=wd)

## HR vs AML
pathenrich(data=data, table = HRvsAML, name="HRvsAML", gene_sets="KEGG",wd=wd)

## HR vs LR
pathenrich(data=data, table = HRvsLR, name="HRvsLR", gene_sets="KEGG",wd=wd)





