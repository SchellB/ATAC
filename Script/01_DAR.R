################################################################################
############################### ATACseq Valeria ################################


# 01 - Data preprocessing

wd <- getwd()

# load the data using readxl
library(readr)
Data <- read.csv(paste0(wd,"/Data/Data.csv"))

# BiomaRt loading
library("biomaRt")

## set up data base
listMarts()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

## filter by HGNC symbols
filters <- c("entrezgene_id")

## select attributes that we want from the query
attributes <- c("hgnc_symbol", "entrezgene_id")

## perform the query
genes <- getBM(attributes=attributes, filters=filters, values=Data$Geneid, mart=ensembl)

# renames the genes
data <- merge(Data, genes, by.x = "Geneid",by.y="entrezgene_id")

# process the data
library(dplyr)
library(tidyverse)
data <- data %>% distinct(hgnc_symbol, .keep_all=TRUE)%>%drop_na() %>% column_to_rownames(var="hgnc_symbol") %>% dplyr::select(c(7:30))


# create a metadata file
Metadata <- data.frame(sample=colnames(data),
                       group=c(rep("HD",6),rep("LR",5),rep("HR",6),rep("sAML",6),"NKalone"),
                       category=c(rep("HD",6),rep("MDS",11),rep("AML",6),"NKalone"))

Metadata$group <- factor(Metadata$group, levels=c("HD","LR","HR","sAML", "NKalone"))

# load edgeR
library(edgeR)


# creating a DGEList object
DGE_list <- DGEList(counts=data, colSums(data),
                    samples=Metadata, group=Metadata$group,
                    genes=rownames(data))

# creating design matrix
design <- model.matrix(~ 0+Metadata$group)

# plot cpm against density
library(ggplot2)
ggplot(data=as.data.frame(cpm(DGE_list$counts)))+
  geom_density(aes(X23_HD.1))+
  geom_density(aes(x=X08_LR.5))+
  geom_density(aes(x=X10_HR.1))+
  geom_density(aes(x=X19_sAML.4))

# plot count distribution
ggplot(data=as.data.frame(DGE_list$counts))+
  geom_density(aes(X23_HD.1), color="blue")+
  geom_density(aes(x=X08_LR.5))+
  geom_density(aes(x=X10_HR.1), color="orange")+
  geom_density(aes(x=X19_sAML.4), color="red")+
  scale_x_continuous(limits=c(0,25))

# filtering by expression
keep <- filterByExpr(DGE_list, design, min.count = 5, min.total.count = 2,large.n = 6, min.prop=0.5)
table(keep)
DGE_list <- DGE_list[keep,]

# plot cpm against density
library(ggplot2)
ggplot(data=as.data.frame(log(cpm(DGE_list$counts))))+
  geom_density(aes(X23_HD.1), color="blue")+
  geom_density(aes(x=X08_LR.5))+
  geom_density(aes(x=X10_HR.1),color="orange")+
  geom_density(aes(x=X19_sAML.4),color="red")

# plot count distribution
ggplot(data=as.data.frame(DGE_list$counts))+
  geom_density(aes(X23_HD.1), color="blue")+
  geom_density(aes(x=X08_LR.5))+
  geom_density(aes(x=X10_HR.1), color="orange")+
  geom_density(aes(x=X19_sAML.4), color="red")+
  scale_x_continuous(limits=c(0,25))

# perform TMM normalisation
library(ggplot2)
library(reshape2)
w <- as.data.frame(DGE_list$counts)
w.plot <- melt(w)

barplot(DGE_list$samples$lib.size*1e-6, ylab="Library size (millions)")

DGE_norm <- calcNormFactors(DGE_list)
Normalized.table <- cpm(DGE_norm)
write.csv2(Normalized.table,file="Results/Normalized.Counts.csv")
z=as.data.frame(DGE_norm$counts)
z.plot <- melt(z) 
ggplot(aes(x=log(value)), colour="black", data=w.plot) + geom_density()+
  geom_density(data = z.plot, aes(x=log(value)), colour="red", linetype=2)

DGE_norm$samples
barplot(cpm(DGE_norm)*1e-6,ylab="Library size (millions)")

# plot MDS
col <- data.frame(group=Metadata$group, col=c(rep("blue",6),
                                              rep("orange",5),
                                              rep("red",6),
                                              rep("green",6),
                                              "black"))
plotMDS(DGE_norm, col=col$col)


# heatmap
logcpm <- cpm(DGE_norm, log=TRUE, normalized.lib.sizes=TRUE)
heatmap(logcpm)


# 02 - Limma Voom ####

# model voom
y <- voom(DGE_norm, design, plot = T)

# fitting the linear model
fit <- lmFit(y, design)
head(coef(fit))
#colnames(fit$coefficients) <- c("HD","LR","HR","sAML","NKalone")

# 02 - Calculate all DEGs ####

# function
calculateDEG <- function(name="AvsB",fit=fit, contrast=c(-1,1,0,0,0), wd=wd){
  AvsB <- contrasts.fit(fit, contrast=contrast)
  AvsB <- eBayes(AvsB)
  print(summary(decideTests(AvsB,method = "separate", adjust.method = "BH", p.value = 0.05,
                      lfc = 1)))
  plotMD(AvsB)
  AvsB.table <- topTable(AvsB, sort.by = "P", n = Inf)
  print(head(AvsB.table, 20))
  write.csv2(AvsB.table,file=paste0(wd,"/Results/",name,".csv"))
  return(AvsB.table)
  
}

## LR Vs HD
LRvsHD <- calculateDEG(name="LRvsHD", fit=fit, contrast=c(-1,1,0,0,0), wd=wd)

## HR Vs HD
HRvsHD<-calculateDEG(name="HRvsHD", fit=fit, contrast=c(-1,0,1,0,0), wd=wd)

## AML Vs HD
AMLvsHD<-calculateDEG(name="AMLvsHD", fit=fit, contrast=c(-1,0,0,1,0), wd=wd)

## LR Vs AML
LRvsAML<-calculateDEG(name="LRvsAML", fit=fit, contrast=c(0,1,0,-1,0), wd=wd)

## HR Vs AML
HRvsAML<-calculateDEG(name="HRvsAML", fit=fit, contrast=c(0,0,1,-1,0), wd=wd)

## HR Vs LR
HRvsLR<-calculateDEG(name="HRvsLR", fit=fit, contrast=c(0,-1,1,0,0), wd=wd)


# 03 - Venn Dagram ####

## Versus HD
# process the data
DEG <- list(LR=subset(LRvsHD$genes,subset=abs(LRvsHD$logFC)>=1&LRvsHD$P.Value<0.05),
            HR=subset(HRvsHD$genes,subset=abs(HRvsHD$logFC)>=1&HRvsHD$P.Value<0.05),
            AML=subset(AMLvsHD$genes,subset=abs(AMLvsHD$logFC)>=1&AMLvsHD$P.Value<0.05))

devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(viridis)
ggvenn(
  DEG, 
  fill_color =viridis(3),
  show_percentage = FALSE,
  stroke_size = 0.5, set_name_size = 4, show_elements = FALSE, text_size=2,label_sep = "\n")

# get output
write.csv2(DEG[["LR"]], file=paste0(wd,"/Results/SignifLRvsHD.csv"))
write.csv2(DEG[["HR"]], file=paste0(wd,"/Results/SignifHRvsHD.csv"))
write.csv2(DEG[["AML"]], file=paste0(wd,"/Results/SignifAMLvsHD.csv"))

## only patho
# process the data
DEG2 <- list(LRversusAML=subset(LRvsAML$genes,subset=abs(LRvsAML$logFC)>=1&LRvsAML$P.Val<0.05),
            HRversusAML=subset(HRvsAML$genes,subset=abs(HRvsAML$logFC)>=1&HRvsAML$P.Val<0.05),
            HRversusLR=subset(HRvsLR$genes,subset=abs(HRvsLR$logFC)>=1&HRvsLR$P.Val<0.05))

ggvenn(
  DEG2, 
  fill_color =viridis(3),
  show_percentage = FALSE,
  stroke_size = 0.5, set_name_size = 4, show_elements = FALSE, text_size=2,label_sep = "\n")

# getoutput
write.csv2(DEG[["LRversusAML"]], file=paste0(wd,"/Results/SignifLRvsAML.csv"))
write.csv2(DEG[["HRversusAML"]], file=paste0(wd,"/Results/SignifHRvsAML.csv"))
write.csv2(DEG[["HRversusLR"]], file=paste0(wd,"/Results/SignifHRvsLR.csv"))

# 04 - Volcanoplot ####

library(ggrepel)
library(viridis)

# function
volcano <- function(table=table, name="name",title="title"){
  
  pdf(paste0(wd,"/Output/Volcano_",name,".pdf"))
  print(ggplot(data=table,
         aes(x=logFC,
             y=-log10(P.Value),
             label=genes)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel(data=subset(table, abs(logFC)>1&P.Value<0.05), aes(x=logFC,
                                                                              y=-log10(P.Value),label=genes),
                    max.overlaps = 50, color="darkgrey") +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    labs(title=title,
         x="log2FC", y="-log10(Pvalue)")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text()))
  
  print(ggplot(data=table,
         aes(x=logFC,
             y=-log10(P.Value),
             label=genes)) +
    geom_point() +
    theme_minimal() +
    geom_point(data=subset(table, abs(logFC)>1&P.Value<0.05),
               aes(x=logFC, y=-log10(P.Value)),color="red")+
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    labs(title=title,
         x="log2FC", y="-log10(Pvalue)")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text()))
  
  dev.off()
  
}

## LR vs HD
volcano(table=LRvsHD,name="LRvsHD",title = "Low Risk MDS versus HD")


## HR vs HD
volcano(table=HRvsHD,name="HRvsHD",title = "High Risk MDS versus HD")


## AML vs HD
volcano(table=AMLvsHD,name="AMLvsHD",title = "sAML MDS versus HD")


## LR vs AML
volcano(table=LRvsAML,name="LRvsAML",title = "Low Risk MDS versus AML")


## HR vs AML
volcano(table=HRvsAML,name="HRvsAML",title = "High Risk MDS versus AML")


## HR vs LR
volcano(table=HRvsLR,name="HRvsLR",title = "High Risk MDS versus Low Risk MDS")

# 05 - saving data ####

save(HRvsAML,HRvsHD,HRvsLR, LRvsAML, LRvsHD, AMLvsHD, DGE_norm, file=paste0(wd, "/Script/ATACseqDEG.Rdata"))

# 06 - Heatmap ####

## data
data <-read_csv2("Results/Normalized.Counts.csv")
data <- data %>% column_to_rownames('...1')

## variables
library(viridis)
cols = list("HD"="blue",
            "LR"="green",
            "HR"="orange","sAML"="red")


## Heatmap
library(ComplexHeatmap)
library(circlize)
library(rlist)

DEG_heatmap <- function(data,table,cond,Metadata,cols,name){
  
  DEGgenes <- subset(table, abs(logFC)>1&P.Value<0.05) # DEG
  df <- data %>% dplyr::filter(rownames(data)%in% DEGgenes$genes) %>% select(contains(c(cond))) # Expr
  Meta <- Metadata %>% filter(sample%in% colnames(df))%>% arrange(group) # Metadata
  
  col_fun = colorRamp2(c(0,10, 50, 100, 150), c("#FDE725FF","#75D054FF","#55C667FF","#404688FF","#440154FF"))
  
  # Create the heatmap annotation
  ha <- HeatmapAnnotation(df=Meta$group)
  
  pdf(file=paste0(wd,"/Output/Heatmap",name,".pdf"))
  print(Heatmap(as.matrix(df), 
                name = name, #title of legend
                column_title = "Patients", row_title = "Genes",
                show_row_dend = FALSE,
                top_annotation = ha,
                col = col_fun,
                row_names_gp = gpar(fontsize = 4) # Text size for row names
  ))
  dev.off()
  
}


# HR vsHD
DEG_heatmap(data=data, Metadata = Metadata, table=HRvsHD,
            name="HRvsHD",cols=cols,cond=c("HD","HR"))

# LR vs HD
DEG_heatmap(data=data, Metadata = Metadata, table=LRvsHD,
            name="LRvsHD",cols=cols,cond=c("HD","LR"))

# AML vs HD
DEG_heatmap(data=data, Metadata = Metadata, table=AMLvsHD,
            name="AMLvsHD",cols=cols,cond=c("HD","AML"))


# LR vs AML
DEG_heatmap(data=data, Metadata = Metadata, table=LRvsAML,
            name="LRvsAML",cols=cols,cond=c("LR","AML"))

# HR vs AML
DEG_heatmap(data=data, Metadata = Metadata, table=HRvsAML,
            name="HRvsAML",cols=cols,cond=c("HR","AML"))

# HR vs LR
DEG_heatmap(data=data, Metadata = Metadata, table=HRvsLR,
            name="HRvsLR",cols=cols,cond=c("HR","LR"))

