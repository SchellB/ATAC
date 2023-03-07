###############################################################################
############################### Dim reduction #################################


# 01 - Starting ####

# working directory
wd <- getwd()

# load data
library(readr)
library(dplyr)
library(tidyverse)
ATAC <- read_csv2("Results/Normalized.Counts.csv")
data <- ATAC %>% column_to_rownames("...1")

# metadata
Metadata <- data.frame(sample=colnames(data),
                       group=c(rep("HR",6),rep("sAML",5),rep("HD",5),"NKalone"))


# correlation plot
library("corrplot")
cor.data <- round(cor(as.matrix(data)),2)

corrplot(cor.data, type="upper", order="hclust", 
         tl.col="black", tl.srt=45)

# install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(data, histogram=TRUE, pch=19)

# 02 - PCA ####
# PCA
library(FactoMineR)
library(factoextra)
res.pca <- PCA(t(data))

eigenvalues <- res.pca$eig

barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")

# variable plot
fviz_pca_var(res.pca, col.var="contrib")

# individual plot
fviz_pca_ind(res.pca, col.ind="cos2") +
  scale_color_gradient2(low="blue", mid="white", 
                        high="red", midpoint=0.50)

# select genes with high cos2
gene.contrib.table <- as.data.frame(res.pca$var$contrib)
genes.interest <- top_n(gene.contrib.table,100,Dim.1)

# do pca on these genes
interest.data <- subset(data, subset=rownames(data)%in%rownames(genes.interest))

# 03 - ploting genes of interest ####
library(tidyr)
long_df <- interest.data %>% rownames_to_column("genes")%>% gather(Key, Value,-genes)

library(viridis)
ggplot(data = long_df, aes(y=genes,x=Value))+
  geom_point()+
  scale_color_manual(values = viridis(3), labels=Metadata$group, guide = "legend")+
  theme_minimal()

# 03 - UMAP ####

#packages
library(plotly) 
library(umap)

# compute umap
umap <- umap(t(data))

layout <-umap[["layout"]] 

layout <- data.frame(layout) 

layout$group <- c(rep("HR",6),rep("sAML",5),rep("HD",5),"NK")


fig <- plot_ly(layout, x = ~X1, y = ~X2, color = ~layout$group, colors = c('#EF553B','#636EFA',
                                                                           "orange",'#00CC96',"black"),
               type = 'scatter', mode = 'markers')%>%  
  
  layout(
    
    plot_bgcolor = "#e5ecf6",
    
    legend=list(title=list(text='groups')))


fig
