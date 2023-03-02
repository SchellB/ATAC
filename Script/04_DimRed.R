###############################################################################
############################### Dim reduction #################################


# 01 - Starting ####

# load data
library(readr)
library(dplyr)
ATAC <- read_csv2("Results/Normalized.Counts.csv")
data <- ATAC %>% column_to_rownames("genes")

# metadata
Metadata <- data.frame(sample=colnames(data),
                       group=c(rep("HD",6),rep("LR",5),rep("HR",6),rep("sAML",6),"NKalone"),
                       category=c(rep("HD",6),rep("MDS",11),rep("AML",6),"NKalone"))


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
res.pca <- PCA(data)

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
gene.contrib.table <- as.data.frame(res.pca$ind$contrib)
genes.interest <- rownames(subset(gene.contrib.table, subset = gene.contrib.table$Dim.1>0.5))

# do pca on these genes
interest.data <- subset(data, subset=rownames(data)%in%genes.interest)

new.pca <- PCA(t(interest.data))


# look variance explained
fviz_screeplot(new.pca, ncp=10)

# variable plot
fviz_pca_var(new.pca, col.var="contrib")

# individual plot
library(viridis)
fviz_pca_ind(new.pca, col.ind=Metadata$group, palette=viridis(5), geom="text")

# biplot
fviz_pca_biplot(new.pca, col.ind=Metadata$group, palette=viridis(5), geom="point",
                alpha.var=0.2, repel=TRUE)+
  theme_minimal()


# 03 - UMAP ####

#packages
library(plotly) 
library(umap)

# compute umap
umap <- umap(t(data))

layout <-umap[["layout"]] 

layout <- data.frame(layout) 

layout$group <- c(rep("HD",6),rep("LR",5),rep("HR",6),rep("AML",6),"NK")


fig <- plot_ly(layout, x = ~X1, y = ~X2, color = ~layout$group, colors = c('#EF553B','#636EFA',
                                                                           "orange",'#00CC96',"black"),
               type = 'scatter', mode = 'markers')%>%  
  
  layout(
    
    plot_bgcolor = "#e5ecf6",
    
    legend=list(title=list(text='groups')))


fig
