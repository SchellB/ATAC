###############################################################################
########################### Correlation ATAC - RNA seq ########################

# 01 - Importing DAR data ####

# directory
wd <- getwd()

# ATAC data
load("ATAC.Rdata")

# RNA data
load("C:/Users/beren/Documents/RNAseqNKValeria/RNAseqTables.Rdata")


# 02 - plot ####
library(ggplot2)

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
RNA <- read.csv2("C:/Users/beren/Documents/RNAseqNKValeria/Results/Normalized.Counts.csv")
data1 <- RNA %>% group_by(X) %>% mutate(
  meanexpHD= (X23_HD.1+X24_HD.2+X25_HD.3+X01_HD.4+X02_HD.5+X03_HD.6),
  meanexpLR= mean(X04_LR.1+X05_LR.2+X07_LR.4+X08_LR.5+X09_LR.6),
  meanexpHR= mean(X10_HR.1+X11_HR.2+X12_HR.3+X13_HR.4+X14_HR.5+X15_HR.6),
  meanexpAML= mean(X16_sAML.1+X17_sAML.2+X18_sAML.3+X19_sAML.4+X20_sAML.5+X21_sAML.6)
)

# expression
data <- merge(data1, data2, by="X",suffixes = c(".RNA",".ATAC"))


# colors
cols <- data.frame("HD"="blue","HR"="orange","sAML"="red")

# plot
ggplot(data = data)+
  geom_point(aes(x = X23_HD.1, y=HD.1, color="HD"))+
  geom_point(aes(x = X24_HD.2, y=HD.2, color="HD"))+
  geom_point(aes(x = X25_HD.3, y=HD.3, color="HD"))+
  geom_point(aes(x = X01_HD.4, y=HD.4, color="HD"))+
  geom_point(aes(x = X03_HD.6, y=HD.6, color="HD"))+
  geom_smooth(aes(x=meanexpHD,y=meanmetHD,color="HD"))+
  geom_point(aes(x = X10_HR.1, y=HR.1, color="HR"))+
  geom_point(aes(x = X11_HR.2, y=HR.2, color="HR"))+
  geom_point(aes(x = X12_HR.3, y=HR.3, color="HR"))+
  geom_point(aes(x = X13_HR.4, y=HR.4, color="HR"))+
  geom_point(aes(x = X14_HR.5, y=HR.5, color="HR"))+
  geom_point(aes(x = X15_HR.6, y=HR.6, color="HR"))+
  geom_smooth(aes(x=meanexpHR,y=meanmetHR,color="HR"))+
  geom_point(aes(x = X16_sAML.1, y=sAML.1, color="sAML"))+
  geom_point(aes(x = X17_sAML.2, y=sAML.2, color="sAML"))+
  geom_point(aes(x = X19_sAML.4, y=sAML.4, color="sAML"))+
  geom_point(aes(x = X20_sAML.5, y=sAML.5, color="sAML"))+
  geom_point(aes(x = X21_sAML.6, y=sAML.6, color="sAML"))+
  geom_smooth(aes(x=meanexpAML,y=meanmetAML,color="sAML"))+
  theme_minimal()+
  scale_x_continuous(limits=c(0,5000))




