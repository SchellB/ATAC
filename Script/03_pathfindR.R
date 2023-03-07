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
load("ATAC.Rdata")

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


# 03 - Performing enrichment ####

## HR vs HD
pathenrich(data=data, table = HRvsHD_ATAC, name="HRvsHD", gene_sets="GO-BP",wd=wd)

## AML vs HD
pathenrich(data=data, table = AMLvsHD_ATAC, name="AMLvsHD", gene_sets="GO-BP",wd=wd)

## HR vs AML
pathenrich(data=data, table = HRvsAML_ATAC, name="HRvsAML", gene_sets="GO-BP",wd=wd)

## HR&AML vs HD
pathenrich(data=data, table = HR_AMLvsHD_ATAC, name="HR_AMLvsHD", gene_sets="GO-All",wd=wd)





