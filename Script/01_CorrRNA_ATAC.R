###############################################################################
########################### Correlation ATAC - RNA seq ########################

# 01 - Importing data ####

# directory
wd <- getwd()

# ATACseq
library(readr)
ATAC <- read_csv2("Results/Normalized.Counts.csv")
colnames(ATAC) <- c("genes",colnames(ATAC[,2:25]))

# RNAseq
RNA <- read_csv2("C:/Users/beren/Documents/RNAseqNKValeria/Results/Normalized.Counts.csv")
colnames(RNA) <- c("genes",colnames(RNA[,2:25]))

# calculate mean per gene
## RNA
library(dplyr)
data1 <- RNA %>% group_by(genes) %>% mutate(
  meanexpHD= ((`01_HD-4`+`02_HD-5`+`03_HD-6`+`23_HD-1`+`24_HD-2`+`25_HD-3`)/6),
  meanexpLR= ((`04_LR-1`+`05_LR-2`+`07_LR-4` +`08_LR-5` + `09_LR-6`)/5),
  meanexpHR= ((`10_HR-1`+`11_HR-2`+`12_HR-3`+`13_HR-4`+`14_HR-5`+`15_HR-6`)/6),
  meanexpAML= ((`16_sAML-1`+`17_sAML-2`+`18_sAML-3`+`19_sAML-4`+`20_sAML-5`+`21_sAML-6`)/6)
)

## ATAC
data2 <- ATAC %>% group_by(genes) %>% mutate(
  meanmetHD= (X23_HD.1+X24_HD.2+X25_HD.3+X01_HD.4+X02_HD.5+X03_HD.6)/6,
  meanmetLR= mean(X04_LR.1+X05_LR.2+X07_LR.4+X08_LR.5+X09_LR.6)/5,
  meanmetHR= mean(X10_HR.1+X11_HR.2+X12_HR.3+X13_HR.4+X14_HR.5+X15_HR.6)/6,
  meanmetAML= mean(X16_sAML.1+X17_sAML.2+X18_sAML.3+X19_sAML.4+X20_sAML.5+X21_sAML.6)/6
)

# merge both files
data <- inner_join(data2[,c("genes","meanmetHD","meanmetHR","meanmetLR","meanmetAML")],
                   data1[,c("genes","meanexpHD", "meanexpLR", "meanexpHR","meanexpAML")], by="genes")

# 02 - Correlation plot ####
library(ggplot2)
colors <- c("HD"="blue","LR"="green","HR"="orange","AML"="red")

ggplot(data=data)+geom_point(aes(x=meanmetHD,y=meanexpHD,color="HD"))+
  geom_smooth(aes(x=meanmetHD,y=meanexpHD, color="HD"))+
  geom_point(aes(x=meanmetLR,y=meanexpLR,color="LR"))+
  geom_smooth(aes(x=meanmetLR,y=meanexpLR,color="LR"))+
  geom_point(aes(x=meanmetHR,y=meanexpHR,color="HR"))+
  geom_smooth(aes(x=meanmetHR,y=meanexpHR,color="HR"))+
geom_point(aes(x=meanmetAML,y=meanexpAML,color="AML"))+
  geom_smooth(aes(x=meanmetAML,y=meanexpAML,color="AML"))+
  scale_x_continuous(limits=c(500,5000))+
  scale_color_manual(values=colors,
                     breaks = c("HD","LR","HR","AML"),
                     labels = c("Healthy Donor","Low Risk MDS","High Risk MDS","sAML"))



# 03 - DEG  data ####

## ATAC DEG HR vs HD

HRvsHD_ATAC <- read_csv2("Results/HRvsHD.csv")

## RNA DEG HR vs HD

#HRvsHD_RNA <- read_csv2("~/Documents/Github/RNAseqNKValeria/Results/HRvsHD.csv")
HRvsHD_RNA <- read_csv2("C:/Users/beren/Documents/RNAseqNKValeria/Results/HRvsHD.csv")


## pulling the data together
library(dplyr)
dt <- inner_join(HRvsHD_ATAC, HRvsHD_RNA,by="genes",suffix = c(".ATAC", ".RNA"))

## comparison plot
library(viridis)
library(ggrepel)

ggplot(data=dt) + geom_point(aes(x=logFC.ATAC, y=logFC.RNA, color=P.Value.ATAC))+
  scale_color_viridis_b()+
  theme_minimal() +
  geom_text_repel(data=subset(dt, abs(logFC.ATAC)>1&abs(logFC.RNA)>1&(P.Value.ATAC<=0.05|P.Value.RNA<=0.05)),
                  aes(x=logFC.ATAC, y=logFC.RNA,label=genes),
                  max.overlaps = 50, color="darkgrey")+
  labs(title="Correlation between ATACseq and RNAseq for HR vs HD")

## select common genes

common_genes <- dt %>% filter(abs(dt$logFC.ATAC)>=1&abs(dt$logFC.RNA)>=1&dt$P.Value.ATAC<=0.05&dt$P.Value.RNA<=0.05) %>% select("genes")


# 04 - Function to automatize the script ####

which_deg <- function(name){
  DEG_ATAC <- read_csv2(paste0("Results/",name,".csv"))
  #HRvsHD_RNA <- read_csv2("~/Documents/Github/RNAseqNKValeria/Results/HRvsHD.csv")
  DEG_RNA <- read_csv2(paste0("C:/Users/beren/Documents/RNAseqNKValeria/Results/",name,".csv"))
  ## pulling the data together
  library(dplyr)
  dt <- inner_join(DEG_ATAC, DEG_RNA,by="genes",suffix = c(".ATAC", ".RNA"))
  #return(dt)
  common_genes <- dt %>% filter(abs(dt$logFC.ATAC)>=1&abs(dt$logFC.RNA)>=1&dt$P.Value.ATAC<=0.05&dt$P.Value.RNA<=0.05)
  return(common_genes)
}

LRvsHD_corrDEG <- which_deg(name="LRvsHD")
AMLvsHD_corrDEG <- which_deg(name="AMLvsHD")
LRvsAML_corrDEG <- which_deg(name="LRvsAML")
HRvsAML_corrDEG <- which_deg(name="HRvsAML")
HRvsLR_corrDEG <- which_deg(name="HRvsLR")



