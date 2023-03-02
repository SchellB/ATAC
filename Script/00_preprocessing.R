################################################################################
######################### Pre processing of ATACseq Data #######################

# 01 - Packages ####

install.packages(c("knitr","rmdformats", "DT","magrittr"))
library(ggplot2,tydir, devtools)

BiocManager::install(c("Rsamtools",
                       "GenomicAlignements",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
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

fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()

fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
                                                                                               437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw()

fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
                                                                             247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()

# TSS sites

## find TSS in human genome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), fix = "start", 1)
TSSs

## plotting ATACseq signal over TSS
library(soGGi)

# Nucleosome free
nucFree <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
                      forceFragment = 50)

# Mononucleosome
monoNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
                      forceFragment = 80)

# Dinucleosome
diNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                    format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
                    forceFragment = 160)

# nucFree_gL <- nucFree monoNuc_gL <- monoNuc diNuc_gL <- diNuc
# save(monoNuc_gL,nucFree_gL,diNuc_gL,file='ATAC_Data/ATAC_RData/gL_soGGiResults.RData')

## plots
library(soGGi)
plotRegion(nucFree, outliers = 0.01)

# isolate by insert size

atacReads_Open <- atacReads[insertSizes < 100, ]
atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]

# create BAM files
openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM)
monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM)
diNucBam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM)

library(rtracklayer)
export(atacReads_Open, openRegionBam, format = "bam")
export(atacReads_MonoNuc, monoNucBam, format = "bam")
# export(atacReads_Open,diNucBam,format = 'bam')


