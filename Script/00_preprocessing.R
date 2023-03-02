################################################################################
######################### Pre processing of ATACseq Data #######################

# 01 - Packages ####

install.packages(c("knitr","rmdformats", "DT","magrittr"))
library(ggplot2,tydir, devtools)

BiocManager::install(c("Rsamtools",
                       "GenomicAlignements",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "soGGi","rtracklayer","ChIPQC",
                       "ChIPseeker","rGREAT"))

BiocManager::install(c("tracktables","clusterProfiler","org.Mm.eg.db",
                       "MotifDb","Biostrings","BSgenome.Hsapiens.UCSC.hg19"))
library(limma)
library(DESeq2)

devtools::install_github('ThomasCarroll/soGGi')

# 02 - QC ####

# Check proportion of mapped reads

## load BAM
sortedBAM <- "C:/Users/beren/OneDrive/Documents/Projets/RNA&ATACseqValeria/BAM_ATAC_44507/15_0DLB_01W5INSERM_HR-1_ATAC_hg38_i219.bam"

## check number of mapped reads
library(Rsubread)
pmapped <- propmapped(sortedBAM)
pmapped # 100% mapped !!

## check distribution of mapped reads
library(Rsamtools)
library(ggplot2)
library(magrittr)

idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
  geom_bar(stat = "identity") + coord_flip()


# Reading mapped reads

library(GenomicAlignments)

atacReads <- readGAlignmentPairs(sortedBAM,
                                 param = ScanBamParam(mapqFilter = 1,
                                                      flag = scanBamFlag(isPaired = TRUE,
                                                                         isProperPair = TRUE),
                                                      what = c("qname",
                                                               "mapq",
                                                               "isize"),
                                                      which = GRanges("chr20", IRanges(1, 63025520))))

length(atacReads)
atacReads

# Insert sizes

## retrieving insert sizes
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)

## plotting insert sizes
library(magrittr)
library(dplyr)
library(ggplot2)
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                     Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
  geom_line()

fragLenPlot + theme_bw()




