################################################################################
############################### ATACseq Valeria ################################


# 01 - Regions ####

wd <- getwd()

# load the data using readxl
library(readxl)
data_regions <- read_excel("Data/01W5Inserm_ATAC_mergedregs.xlsx", 
                                          col_types = c("numeric", "text", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "text", "text", 
                                                        "text", "text", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric"))

# 02 - Genes ####

# Import Data
library(readxl)
data_genes <- read_excel("Data/01W5Inserm_ATAC_genes.xlsx", 
                                     col_types = c("text", "text", "text", 
                                                   "numeric", "numeric", "numeric", 
                                                   "text", "numeric", "text", "text", 
                                                   "text", "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric", "text", "numeric", "text", 
                                                   "text", "text", "text", "numeric", 
                                                   "text"))

# 03 - Calculated DEG ###
# process the data
library(dplyr)
library(tidyverse)
data <- data_genes  %>% column_to_rownames(var="Gene Name") %>% dplyr::select(c(30:46))
patients <- gsub("_ATAC.*","",colnames(data))
colnames(data) <- gsub(".*[[:digit:]][[:digit:]]_INSERM_","",patients)


# create a metadata file
Metadata <- data.frame(sample=colnames(data),
                       group=c(rep("HR",6),rep("AML",5),rep("HD",5),"NKalone"))

Metadata$group <- factor(Metadata$group, levels=c("HD","HR","AML", "NKalone"))

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
  geom_density(aes(x=`HR-1`))+
  geom_density(aes(x=`sAML-1`))
  #geom_density(aes(X=`HD-6`))

# plot count distribution
ggplot(data=as.data.frame(DGE_list$counts))+
  geom_density(aes(x=log(`HR-1`)))+
  geom_density(aes(x=log(`sAML-1`)))

# filtering by expression
keep <- filterByExpr(DGE_list, design, min.count = 5, min.total.count = 2,large.n = 6, min.prop=0.5)
table(keep)
DGE_list <- DGE_list[keep,]

# plot cpm against density
library(ggplot2)
ggplot(data=as.data.frame(log(cpm(DGE_list$counts))))+
  geom_density(aes(x=`HR-1`))+
  geom_density(aes(x=`sAML-1`))

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
col <- data.frame(group=Metadata$group, col=c(rep("orange",6),
                                              rep("red",5),
                                              rep("blue",5),
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


## HR Vs HD
HRvsHD_ATAC<-calculateDEG(name="HRvsHD", fit=fit, contrast=c(-1,1,0,0), wd=wd)

## AML Vs HD
AMLvsHD_ATAC<-calculateDEG(name="AMLvsHD", fit=fit, contrast=c(-1,0,1,0), wd=wd)

## HR Vs AML
HRvsAML_ATAC<-calculateDEG(name="HRvsAML", fit=fit, contrast=c(0,1,-1,0), wd=wd)

## HR&AML Vs HD
HR_AMLvsHD_ATAC<-calculateDEG(name="HR_AMLvsHD", fit=fit, contrast=c(-1,0.5,0.5,0), wd=wd)


# save DEG data
save(HRvsHD_ATAC,AMLvsHD_ATAC,HRvsAML_ATAC,HR_AMLvsHD_ATAC, file="ATAC.Rdata")

# 03 - Venn Dagram ####

## Versus HD
# process the data
DEG <- list(HR=subset(HRvsHD$genes,subset=abs(HRvsHD$logFC)>=1&HRvsHD$P.Value<0.05),
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
write.csv2(DEG[["HR"]], file=paste0(wd,"/Results/SignifHRvsHD.csv"))
write.csv2(DEG[["AML"]], file=paste0(wd,"/Results/SignifAMLvsHD.csv"))


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

## HR vs HD
volcano(table=HRvsHD,name="HRvsHD",title = "High Risk MDS versus HD")


## AML vs HD
volcano(table=AMLvsHD,name="AMLvsHD",title = "sAML versus HD")

## HR vs AML
volcano(table=HRvsAML,name="HRvsAML",title = "High Risk MDS versus AML")


## HR&AML vs HD
volcano(table=HR_AMLvsHD,name="HR&AMLvsHD",title = "High Risk MDS & AML versus HD")

# 05 - saving data ####

save(HRvsAML,HRvsHD,HRvsLR, LRvsAML, LRvsHD, AMLvsHD, DGE_norm, file=paste0(wd, "/Script/ATACseqDEG.Rdata"))

# 06 - Heatmap ####

## data
data <-read_csv2("Results/Normalized.Counts.csv")
data <- data %>% column_to_rownames('...1')

## variables
library(viridis)
cols = list("HD"="blue",
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

# AML vs HD
DEG_heatmap(data=data, Metadata = Metadata, table=AMLvsHD,
            name="AMLvsHD",cols=cols,cond=c("HD","AML"))

# HR vs AML
DEG_heatmap(data=data, Metadata = Metadata, table=HRvsAML,
            name="HRvsAML",cols=cols,cond=c("HR","AML"))

# HR & AML vs HD
DEG_heatmap(data=data, Metadata = Metadata, table=HR_AMLvsHD,
            name="HR&AMLvsHD",cols=cols,cond=c("HR","AML","LR"))








