###############################################################################
########################### Correlation ATAC - RNA seq ########################

# 01 - Importing DAR data ####

# directory
wd <- getwd()

# ATAC data
load("ATAC.Rdata")

# RNA data
#load("C:/Users/beren/Documents/RNAseqNKValeria/RNAseqTables.Rdata")
load("~/Documents/Github/RNAseqNKValeria/RNAseqTables.Rdata")

# 02 - plot ####
library(ggplot2)
library(ggrepel)

# funciton

CommonDEG <- function(a,b, name="name",title="title") {
  data <- merge(a, b, by="genes",suffixes = c(".RNA",".ATAC"))

  pdf(paste0(wd,"/Output/CommonVolcano_",name,".pdf"))  
  print(ggplot(data = data)+
    geom_point(aes(x = logFC.RNA, y=-log10(P.Value.RNA)),color="blue",alpha=0.3)+
    geom_point(aes(x = logFC.ATAC, y=-log10(P.Value.ATAC)),color="red",alpha=0.3)+
    theme_minimal()+
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    geom_text_repel(data=subset(data,
                                abs(logFC.RNA)>1&P.Value.RNA<0.05&abs(logFC.ATAC)>1&P.Value.ATAC<0.05),
                    aes(x=logFC.RNA,
                        y=-log10(P.Value.RNA),label=genes),
                    max.overlaps = 50, color="black") +
    labs(title=title,
         x="log2FC", y="-log10(Pvalue)")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text()))
  dev.off()
  
  write.csv2(data,file=paste0(wd,"/Results/",name,".csv"))
  
}


# HR vs HD
CommonDEG(a=HRvsHD,b=HRvsHD_ATAC,name="commonHRvsHD",title ="HR vs HD ATACseq and RNAseq" )

# AML vs HD
CommonDEG(a=AMLvsHD,b=AMLvsHD_ATAC,name="commonAMLvsHD",title ="AML vs HD ATACseq and RNAseq" )


# 03 - Plotting Raw Expression and Opening ####

# ATAC data
library(readr)
library(dplyr)
ATAC <- read.csv2("Results/Normalized.Counts.csv")
data2 <- ATAC %>% group_by(X) %>% mutate(
  meanmetHD= mean(HD.1,HD.2,HD.3,HD.4,HD.6),
  meanmetHR= mean(HR.1,HR.2,HR.3,HR.4,HR.5,HR.6),
  meanmetAML= mean(sAML.1,sAML.2,sAML.4,sAML.5,sAML.6))

# RNA data
#RNA <- read.csv2("C:/Users/beren/Documents/RNAseqNKValeria/Results/Normalized.Counts.csv")
RNA <- read.csv2("~/Documents/Github/RNAseqNKValeria/Results/Normalized.Counts.csv")
data1 <- RNA %>% group_by(X) %>% mutate(
  meanexpHD= (HD.1_RNA+HD.2_RNA+HD.3_RNA+HD.4_RNA+HD.5_RNA+HD.6_RNA),
  meanexpLR= mean(LR.1_RNA+ LR.2_RNA+LR.4_RNA+LR.5_RNA+LR.6_RNA),
  meanexpHR= mean(HR.1_RNA+HR.2_RNA+HR.3_RNA+HR.4_RNA+HR.5_RNA+HR.6_RNA),
  meanexpAML= mean(sAML.1_RNA+sAML.2_RNA+sAML.3_RNA+sAML.4_RNA+sAML.5_RNA+sAML.6_RNA)
)

# expression
data <- merge(data1, data2, by="X",suffixes = c(".RNA",".ATAC"))


# colors
cols <- data.frame("HD"="blue","HR"="orange","sAML"="red")

# plot
ggplot(data = data)+
  geom_point(aes(x = HD.1_RNA, y=HD.1, color="HD"))+
  geom_point(aes(x = HD.2_RNA, y=HD.2, color="HD"))+
  geom_point(aes(x = HD.3_RNA, y=HD.3, color="HD"))+
  geom_point(aes(x = HD.4_RNA, y=HD.4, color="HD"))+
  geom_point(aes(x = HD.6_RNA, y=HD.6, color="HD"))+
  geom_smooth(aes(x=meanexpHD,y=meanmetHD,color="HD"))+
  geom_point(aes(x = HR.1_RNA, y=HR.1, color="HR"))+
  geom_point(aes(x = HR.2_RNA, y=HR.2, color="HR"))+
  geom_point(aes(x = HR.3_RNA, y=HR.3, color="HR"))+
  geom_point(aes(x = HR.4_RNA, y=HR.4, color="HR"))+
  geom_point(aes(x = HR.5_RNA, y=HR.5, color="HR"))+
  geom_point(aes(x = HR.6_RNA, y=HR.6, color="HR"))+
  geom_smooth(aes(x=meanexpHR,y=meanmetHR,color="HR"))+
  geom_point(aes(x = sAML.1_RNA, y=sAML.1, color="sAML"))+
  geom_point(aes(x = sAML.2_RNA, y=sAML.2, color="sAML"))+
  geom_point(aes(x = sAML.4_RNA, y=sAML.4, color="sAML"))+
  geom_point(aes(x = sAML.5_RNA, y=sAML.5, color="sAML"))+
  geom_point(aes(x = sAML.6_RNA, y=sAML.6, color="sAML"))+
  geom_smooth(aes(x=meanexpAML,y=meanmetAML,color="sAML"))+
  theme_minimal()+
  scale_x_continuous(limits=c(0,5000))+
  scale_y_continuous(limits=c(0,100))+
  labs(y="chromatin accessibility", x="gene expression")


# 04 - DEG --> opening ####

data <- subset(ATAC, X%in% HRvsHD[abs(HRvsHD$logFC)>=1&HRvsHD$P.Value<=0.05,"genes"])
library(tidyr)
long_df <- data %>% gather(Key, Value,-X) %>% mutate(group=case_when(
  grepl("HR", Key) ~ "HR",
  grepl("HD", Key) ~ "HD",
  grepl("sAML", Key) ~ "sAML",
  .default = "NK alone"
))

data <- merge(long_df,HRvsHD, by.x="X",by.y = "genes")

library(viridis)
ggplot(data = data, aes(y=X,x=Value, color=group))+
  geom_point()+
  scale_color_manual(values = viridis(3), labels=Mgroup, guide = "legend")+
  theme_minimal()

ggplot(data = data)+ facet_wrap(~X)+
  geom_point(aes(x = logFC, y=-log10(P.Value)),color="blue",alpha=0.3)+
  theme_minimal()+
  labs(title=title,
       x="log2FC", y="-log10(Pvalue)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text())

