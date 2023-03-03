###############################################################################
########################### Correlation ATAC - RNA seq ########################

# 01 - Importing data ####

# directory
wd <- getwd()

# ATAC data
load("ATAC.Rdata")

# RNA data
load("C:/Users/beren/Documents/RNAseqNKValeria/RNAseqTables.Rdata")


# 02 - plot ####

# merge data HR vs HD
data <- merge(HRvsHD, HRvsHD_ATAC, by="genes",suffixes = c(".RNA",".ATAC"))

# plot HR vs HD

library(ggplot2)

ggplot(data = data)+
  geom_point(aes(x = logFC.RNA, y=logFC.ATAC, color=P.Value.RNA))+theme_minimal()


ggplot(data = data)+
  geom_point(aes(x = AveExpr.RNA, y=AveExpr.ATAC, color=P.Value.RNA))+theme_minimal()


# merge data AML vs HD
data <- merge(AMLvsHD, AMLvsHD_ATAC, by="genes",suffixes = c(".RNA",".ATAC"))

# plot AML vs HD

library(ggplot2)

ggplot(data = data)+
  geom_point(aes(x = logFC.RNA, y=logFC.ATAC, color=P.Value.RNA))+theme_minimal()


ggplot(data = data)+
  geom_point(aes(x = AveExpr.RNA, y=AveExpr.ATAC, color=P.Value.RNA))+theme_minimal()


# subset signif

ggplot(data = subset(data, subset=abs(data$logFC.RNA)&data$P.Value.RNA<0.05))+
  geom_point(aes(x = logFC.RNA, y=logFC.ATAC, color=P.Value.RNA))+theme_minimal()




