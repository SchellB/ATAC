################################################################################
############################### TARGETED Pathway ###############################


# 01 - Importing DAR data ####

# directory
wd <- getwd()
githubdir <- ("D:/Berenice")

# ATAC
library(readr)
ATAC <- read.csv2("Results/Normalized.Counts.csv")

# RNA
RNA <- read.csv2(paste0(githubdir,"/RNAseqNKValeria/Results/Normalized.Counts.csv"))
colnames(RNA) <- c("X",gsub(".*_","",colnames(RNA[,2:length(RNA)])))

# 02 - Find pathway of interest ####

# load package
library(KEGGREST)
library(org.Hs.eg.db)
library(dplyr)

# find all pathways
hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))

# find all genes
hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )

# find pathway description
hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)

# 03 - Natural killer cell cytotoxicity ####

# find the pathway hsa04650
genes <- hsa_kegg_anno[hsa_kegg_anno$pathway=="path:hsa04650",]

# get the expression from the genes
library(tidyr)
table_rna <- RNA[RNA$X%in%genes$symbol,]
table_rna <- table_rna %>% gather(Key, Value,-X)

# get the accessibility from the genes
table_atac <- ATAC[ATAC$X%in%genes$symbol,]
table_atac <- table_atac %>% gather(Key, Value,-X)

# merge both
table <- merge(table_atac,table_rna,by=c("X","Key"), suffixes=c("ATAC","RNA")) %>% mutate(group=case_when(
  grepl("HR", Key) ~ "HR",
  grepl("HD", Key) ~ "HD",
  grepl("sAML", Key) ~ "sAML",
  .default = "NK alone"
))

# plot
library(ggplot2)
library(viridis)

ggplot(data = table,aes(x=group,y=ValueRNA,color=ValueATAC))+
  geom_point()+theme_minimal()+
  facet_wrap(~X,scales = "free_y")+ scale_color_viridis()+
  theme(axis.text.x=element_text(angle = 60))

