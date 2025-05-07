#######################################
#ANNOTATION REGARDING THE GENIC POSITION
#######################################
#2025/01/23
#Uxue Lazkano

###Summary#####
# This script is used to assign a unique annotation to each CpG.

#R in lamarr
R
.libPaths("/home/CICBIOGUNE/ulazcano/R_lamarr/epigenomics")


##################
# Libraries
##################
library(methylKit)
library(annotatr)
library(bsseq)
library(plyranges)
library(UCSCRepeatMasker)
library(dplyr)
library(genomation)
library(GenomicFeatures)
library(ChIPseeker)
library(data.table)

setwd("W:/ulazcano/AC80_RRBS/05_Methylation_extraction/")
dir_infiles <- "W:/ulazcano/AC80_RRBS/05_Methylation_extraction/"

setwd("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction")
dir_infiles <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/"


#Read coverage files (containing cytosines only from my samples)
#List all files TIENE QUE ESTAR TODO AS LIST
file_list<-list.files(path=dir_infiles, pattern = "*_bismark_bt2_pe.bismark.cov.gz", all.files=FALSE,
                      full.names=FALSE)
path_list <- as.list(paste0(dir_infiles,'/', unlist(file_list)))

sample_name <- strsplit(file_list, "_R1_val_1_bismark_bt2_pe.bismark.cov.gz")

##########################################
# Load data (.cov) from MspI + BS [S01,S04,S07]
##########################################
# Select S01, S04, S07
mspI_BS_paths <- Filter(function(x) any(grepl("S01|S04|S07", x)), path_list)
mspI_BS_sample_name <- sapply(mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspI_BS_sample_name <- as.list(mspI_BS_sample_name)

mspI_BS_BismarkCoverage<- methRead(mspI_BS_paths, sample.id = mspI_BS_sample_name,  
                                   assembly = "hg38", mincov =0, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspI_BS_BismarkCoverage_meth<- methylKit::unite(mspI_BS_BismarkCoverage)
mspI_BS_df <- as.data.frame(mspI_BS_BismarkCoverage_meth)
mspI_BS_CpGs <- mspI_BS_df[mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspI_BS_CpGs)
# 9203822      13
rm(mspI_BS_BismarkCoverage_meth)
rm(mspI_BS_BismarkCoverage)
rm(mspI_BS_df)
#To GRanges for annotation
mspI_BS_CpGs_gr <- as(mspI_BS_CpGs,"GRanges")
seqlevelsStyle(mspI_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-TaqaI + BS [S02,S05,S08]
##########################################
# Select S02, S05, S08
mspItaqaI_BS_paths <- Filter(function(x) any(grepl("S02|S05|S08", x)), path_list)
mspItaqaI_BS_sample_name <- sapply(mspItaqaI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspItaqaI_BS_sample_name <- as.list(mspItaqaI_BS_sample_name)

mspItaqaI_BS_BismarkCoverage<- methRead(mspItaqaI_BS_paths, sample.id = mspItaqaI_BS_sample_name,  
                                   assembly = "hg38", mincov =0, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspItaqaI_BS_BismarkCoverage_meth<- methylKit::unite(mspItaqaI_BS_BismarkCoverage)
mspItaqaI_BS_df <- as.data.frame(mspItaqaI_BS_BismarkCoverage_meth)
mspItaqaI_BS_CpGs <- mspItaqaI_BS_df[mspItaqaI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspItaqaI_BS_CpGs)
# 9580989      13
rm(mspItaqaI_BS_BismarkCoverage_meth)
rm(mspItaqaI_BS_BismarkCoverage)
rm(mspItaqaI_BS_df)
#To GRanges for annotation
mspItaqaI_BS_CpGs_gr <- as(mspItaqaI_BS_CpGs,"GRanges")
seqlevelsStyle(mspItaqaI_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-HaeIII + BS [S03,S06,S09]
##########################################
# Select S03, S06, S09
mspIhaeIII_BS_paths <- Filter(function(x) any(grepl("S03|S06|S09", x)), path_list)
mspIhaeIII_BS_sample_name <- sapply(mspIhaeIII_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspIhaeIII_BS_sample_name <- as.list(mspIhaeIII_BS_sample_name)

mspIhaeIII_BS_BismarkCoverage<- methRead(mspIhaeIII_BS_paths, sample.id = mspIhaeIII_BS_sample_name,  
                                        assembly = "hg38", mincov =0, 
                                        treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspIhaeIII_BS_BismarkCoverage_meth<- methylKit::unite(mspIhaeIII_BS_BismarkCoverage)
mspIhaeIII_BS_df <- as.data.frame(mspIhaeIII_BS_BismarkCoverage_meth)
mspIhaeIII_BS_CpGs <- mspIhaeIII_BS_df[mspIhaeIII_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspIhaeIII_BS_CpGs)
# 13341644       13

rm(mspIhaeIII_BS_BismarkCoverage_meth)
rm(mspIhaeIII_BS_BismarkCoverage)
rm(mspIhaeIII_BS_df)
#To GRanges for annotation
mspIhaeIII_BS_CpGs_gr <- as(mspIhaeIII_BS_CpGs,"GRanges")
seqlevelsStyle(mspIhaeIII_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI + enzymatic [S10,S13,S16]
##########################################
# Select S10,S13,S16

mspI_enzym_paths <- Filter(function(x) any(grepl("S10|S13|S16", x)), path_list)
mspI_enzym_sample_name <- sapply(mspI_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspI_enzym_sample_name <- as.list(mspI_enzym_sample_name)

mspI_enzym_BismarkCoverage<- methRead(mspI_enzym_paths, sample.id = mspI_enzym_sample_name,  
                                   assembly = "hg38", mincov =0, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspI_enzym_BismarkCoverage_meth<- methylKit::unite(mspI_enzym_BismarkCoverage)
mspI_enzym_df <- as.data.frame(mspI_enzym_BismarkCoverage_meth)
mspI_enzym_CpGs <- mspI_enzym_df[mspI_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspI_enzym_CpGs)
# 12270483       13
rm(mspI_enzym_BismarkCoverage_meth)
rm(mspI_enzym_BismarkCoverage)
rm(mspI_enzym_df)
#To GRanges for annotation
mspI_enzym_CpGs_gr <- as(mspI_enzym_CpGs,"GRanges")
seqlevelsStyle(mspI_enzym_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-TaqaI + enzymatic [S11,S14,S17]
##########################################
# Select S11,S14,S17

mspItaqaI_enzym_paths <- Filter(function(x) any(grepl("S11|S14|S17", x)), path_list)
mspItaqaI_enzym_sample_name <- sapply(mspItaqaI_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspItaqaI_enzym_sample_name <- as.list(mspItaqaI_enzym_sample_name)

mspItaqaI_enzym_BismarkCoverage<- methRead(mspItaqaI_enzym_paths, sample.id = mspItaqaI_enzym_sample_name,  
                                      assembly = "hg38", mincov =0, 
                                      treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspItaqaI_enzym_BismarkCoverage_meth<- methylKit::unite(mspItaqaI_enzym_BismarkCoverage)
mspItaqaI_enzym_df <- as.data.frame(mspItaqaI_enzym_BismarkCoverage_meth)
mspItaqaI_enzym_CpGs <- mspItaqaI_enzym_df[mspItaqaI_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspItaqaI_enzym_CpGs)
# 14491039       13
rm(mspItaqaI_enzym_BismarkCoverage_meth)
rm(mspItaqaI_enzym_BismarkCoverage)
rm(mspItaqaI_enzym_df)
#To GRanges for annotation
mspItaqaI_enzym_CpGs_gr <- as(mspItaqaI_enzym_CpGs,"GRanges")
seqlevelsStyle(mspItaqaI_enzym_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-haeIII + enzymatic [S12,S15,S18]
##########################################
# Select S12,S15,S18
mspIhaeIII_enzym_paths <- Filter(function(x) any(grepl("S12|S15|S18", x)), path_list)
mspIhaeIII_enzym_sample_name <- sapply(mspIhaeIII_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspIhaeIII_enzym_sample_name <- as.list(mspIhaeIII_enzym_sample_name)

mspIhaeIII_enzym_BismarkCoverage<- methRead(mspIhaeIII_enzym_paths, sample.id = mspIhaeIII_enzym_sample_name,  
                                           assembly = "hg38", mincov =0, 
                                           treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspIhaeIII_enzym_BismarkCoverage_meth<- methylKit::unite(mspIhaeIII_enzym_BismarkCoverage)
mspIhaeIII_enzym_df <- as.data.frame(mspIhaeIII_enzym_BismarkCoverage_meth)
mspIhaeIII_enzym_CpGs <- mspIhaeIII_enzym_df[mspIhaeIII_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspIhaeIII_enzym_CpGs)
# 18251289       13
rm(mspIhaeIII_enzym_BismarkCoverage_meth)
rm(mspIhaeIII_enzym_BismarkCoverage)
rm(mspIhaeIII_enzym_df)
#To GRanges for annotation
mspIhaeIII_enzym_CpGs_gr <- as(mspIhaeIII_enzym_CpGs,"GRanges")
seqlevelsStyle(mspIhaeIII_enzym_CpGs_gr) = "UCSC"  # necessary

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/AC80_RRBS.RData")
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/AC80_RRBS.RData")
