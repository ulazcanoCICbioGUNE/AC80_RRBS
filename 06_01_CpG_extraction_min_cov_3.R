#######################################
#ANNOTATION REGARDING THE GENIC POSITION
#######################################
#2025/01/23
#Uxue Lazkano

###Summary#####
# This script is used to assign a unique annotation to each CpG.

#R in lamarr
R
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")


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

setwd("W:/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/")
dir_infiles <- "W:/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/"

setwd("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/")
dir_infiles <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/"


#Read coverage files (containing cytosines only from my samples)
#List all files TIENE QUE ESTAR TODO AS LIST
file_list<-list.files(path=dir_infiles, pattern = ".CpG_report.merged_CpG_evidence.cov", all.files=FALSE,
                      full.names=FALSE)
path_list <- as.list(paste0(dir_infiles,'/', unlist(file_list)))

sample_name <- strsplit(file_list, ".CpG_report.merged_CpG_evidence.cov")

##########################################
# Load data (.cov) from MspI + BS [S01,S04,S07]
##########################################
# Select S01
S01_mspI_BS_paths <- Filter(function(x) any(grepl("S01", x)), path_list)
S01_mspI_BS_sample_name <- sapply(S01_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S01_mspI_BS_sample_name <- as.list(S01_mspI_BS_sample_name)

S01_mspI_BS_BismarkCoverage<- methRead(S01_mspI_BS_paths, sample.id = S01_mspI_BS_sample_name,  
                                   assembly = "hg38",  mincov = 3, 
                                   treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S01_mspI_BS_df <- as.data.frame(S01_mspI_BS_BismarkCoverage)

S01_mspI_BS_CpGs <- S01_mspI_BS_df[S01_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S01_mspI_BS_CpGs)
# 6138838      13
rm(S01_mspI_BS_BismarkCoverage)
rm(S01_mspI_BS_df)
#To GRanges for annotation
S01_mspI_BS_CpGs_gr <- as(S01_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S01_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S04
S04_mspI_BS_paths <- Filter(function(x) any(grepl("S04", x)), path_list)
S04_mspI_BS_sample_name <- sapply(S04_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S04_mspI_BS_sample_name <- as.list(S04_mspI_BS_sample_name)

S04_mspI_BS_BismarkCoverage<- methRead(S04_mspI_BS_paths, sample.id = S04_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S04_mspI_BS_df <- as.data.frame(S04_mspI_BS_BismarkCoverage)

S04_mspI_BS_CpGs <- S04_mspI_BS_df[S04_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S04_mspI_BS_CpGs)
# 7158730      13
rm(S04_mspI_BS_BismarkCoverage)
rm(S04_mspI_BS_df)
#To GRanges for annotation
S04_mspI_BS_CpGs_gr <- as(S04_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S04_mspI_BS_CpGs_gr) = "UCSC"  # necessary


# Select S07
S07_mspI_BS_paths <- Filter(function(x) any(grepl("S07", x)), path_list)
S07_mspI_BS_sample_name <- sapply(S07_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S07_mspI_BS_sample_name <- as.list(S07_mspI_BS_sample_name)

S07_mspI_BS_BismarkCoverage<- methRead(S07_mspI_BS_paths, sample.id = S07_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S07_mspI_BS_df <- as.data.frame(S07_mspI_BS_BismarkCoverage)

S07_mspI_BS_CpGs <- S07_mspI_BS_df[S07_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S07_mspI_BS_CpGs)
# 7696699      13
rm(S07_mspI_BS_BismarkCoverage)
rm(S07_mspI_BS_df)
#To GRanges for annotation
S07_mspI_BS_CpGs_gr <- as(S07_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S07_mspI_BS_CpGs_gr) = "UCSC"  # necessary


# Select S01, S04, S07
mspI_BS_paths <- Filter(function(x) any(grepl("S01|S04|S07", x)), path_list)
mspI_BS_sample_name <- sapply(mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspI_BS_sample_name <- as.list(mspI_BS_sample_name)

mspI_BS_BismarkCoverage<- methRead(mspI_BS_paths, sample.id = mspI_BS_sample_name,  
                                   assembly = "hg38", mincov =3, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspI_BS_BismarkCoverage_meth<- methylKit::unite(mspI_BS_BismarkCoverage)
mspI_BS_df <- as.data.frame(mspI_BS_BismarkCoverage_meth)

mspI_BS_CpGs <- mspI_BS_df[mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspI_BS_CpGs)
# 4316467      13
rm(mspI_BS_BismarkCoverage)
rm(mspI_BS_df)
#To GRanges for annotation
mspI_BS_CpGs_gr <- as(mspI_BS_CpGs,"GRanges")
seqlevelsStyle(mspI_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-TaqaI + BS [S02,S05,S08]
##########################################
# Select S02
S02_mspI_BS_paths <- Filter(function(x) any(grepl("S02", x)), path_list)
S02_mspI_BS_sample_name <- sapply(S02_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S02_mspI_BS_sample_name <- as.list(S02_mspI_BS_sample_name)

S02_mspI_BS_BismarkCoverage<- methRead(S02_mspI_BS_paths, sample.id = S02_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S02_mspI_BS_df <- as.data.frame(S02_mspI_BS_BismarkCoverage)

S02_mspI_BS_CpGs <- S02_mspI_BS_df[S02_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S02_mspI_BS_CpGs)
# 7210733      13
rm(S02_mspI_BS_BismarkCoverage)
rm(S02_mspI_BS_df)
#To GRanges for annotation
S02_mspI_BS_CpGs_gr <- as(S02_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S02_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S05
S05_mspI_BS_paths <- Filter(function(x) any(grepl("S05", x)), path_list)
S05_mspI_BS_sample_name <- sapply(S05_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S05_mspI_BS_sample_name <- as.list(S05_mspI_BS_sample_name)

S05_mspI_BS_BismarkCoverage<- methRead(S05_mspI_BS_paths, sample.id = S05_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S05_mspI_BS_df <- as.data.frame(S05_mspI_BS_BismarkCoverage)

S05_mspI_BS_CpGs <- S05_mspI_BS_df[S05_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S05_mspI_BS_CpGs)
# 7627119      13
rm(S05_mspI_BS_BismarkCoverage)
rm(S05_mspI_BS_df)
#To GRanges for annotation
S05_mspI_BS_CpGs_gr <- as(S05_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S05_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S08
S08_mspI_BS_paths <- Filter(function(x) any(grepl("S08", x)), path_list)
S08_mspI_BS_sample_name <- sapply(S08_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S08_mspI_BS_sample_name <- as.list(S08_mspI_BS_sample_name)

S08_mspI_BS_BismarkCoverage<- methRead(S08_mspI_BS_paths, sample.id = S08_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S08_mspI_BS_df <- as.data.frame(S08_mspI_BS_BismarkCoverage)

S08_mspI_BS_CpGs <- S08_mspI_BS_df[S08_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S08_mspI_BS_CpGs)
# 7472889      13
rm(S08_mspI_BS_BismarkCoverage)
rm(S08_mspI_BS_df)
#To GRanges for annotation
S08_mspI_BS_CpGs_gr <- as(S08_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S08_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S02, S05, S08
mspItaqaI_BS_paths <- Filter(function(x) any(grepl("S02|S05|S08", x)), path_list)
mspItaqaI_BS_sample_name <- sapply(mspItaqaI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspItaqaI_BS_sample_name <- as.list(mspItaqaI_BS_sample_name)

mspItaqaI_BS_BismarkCoverage<- methRead(mspItaqaI_BS_paths, sample.id = mspItaqaI_BS_sample_name,  
                                   assembly = "hg38", mincov =3, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspItaqaI_BS_BismarkCoverage_meth<- methylKit::unite(mspItaqaI_BS_BismarkCoverage)
mspItaqaI_BS_df <- as.data.frame(mspItaqaI_BS_BismarkCoverage_meth)
mspItaqaI_BS_CpGs <- mspItaqaI_BS_df[mspItaqaI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspItaqaI_BS_CpGs)
# 4467404      13
rm(mspItaqaI_BS_BismarkCoverage_meth)
rm(mspItaqaI_BS_BismarkCoverage)
rm(mspItaqaI_BS_df)
#To GRanges for annotation
mspItaqaI_BS_CpGs_gr <- as(mspItaqaI_BS_CpGs,"GRanges")
seqlevelsStyle(mspItaqaI_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-HaeIII + BS [S03,S06,S09]
##########################################
# Select S03
S03_mspI_BS_paths <- Filter(function(x) any(grepl("S03", x)), path_list)
S03_mspI_BS_sample_name <- sapply(S03_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S03_mspI_BS_sample_name <- as.list(S03_mspI_BS_sample_name)

S03_mspI_BS_BismarkCoverage<- methRead(S03_mspI_BS_paths, sample.id = S03_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S03_mspI_BS_df <- as.data.frame(S03_mspI_BS_BismarkCoverage)

S03_mspI_BS_CpGs <- S03_mspI_BS_df[S03_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S03_mspI_BS_CpGs)
# 9577772      13
rm(S03_mspI_BS_BismarkCoverage)
rm(S03_mspI_BS_df)
#To GRanges for annotation
S03_mspI_BS_CpGs_gr <- as(S03_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S03_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S06
S06_mspI_BS_paths <- Filter(function(x) any(grepl("S06", x)), path_list)
S06_mspI_BS_sample_name <- sapply(S06_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S06_mspI_BS_sample_name <- as.list(S06_mspI_BS_sample_name)

S06_mspI_BS_BismarkCoverage<- methRead(S06_mspI_BS_paths, sample.id = S06_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S06_mspI_BS_df <- as.data.frame(S06_mspI_BS_BismarkCoverage)

S06_mspI_BS_CpGs <- S06_mspI_BS_df[S06_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S06_mspI_BS_CpGs)
# 9523707      13
rm(S06_mspI_BS_BismarkCoverage)
rm(S06_mspI_BS_df)
#To GRanges for annotation
S06_mspI_BS_CpGs_gr <- as(S06_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S06_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S09
S09_mspI_BS_paths <- Filter(function(x) any(grepl("S09", x)), path_list)
S09_mspI_BS_sample_name <- sapply(S09_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S09_mspI_BS_sample_name <- as.list(S09_mspI_BS_sample_name)

S09_mspI_BS_BismarkCoverage<- methRead(S09_mspI_BS_paths, sample.id = S09_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S09_mspI_BS_df <- as.data.frame(S09_mspI_BS_BismarkCoverage)

S09_mspI_BS_CpGs <- S09_mspI_BS_df[S09_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S09_mspI_BS_CpGs)
# 10321580      13
rm(S09_mspI_BS_BismarkCoverage)
rm(S09_mspI_BS_df)
#To GRanges for annotation
S09_mspI_BS_CpGs_gr <- as(S09_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S09_mspI_BS_CpGs_gr) = "UCSC"  # necessary


# Select S03, S09, S09
mspIhaeIII_BS_paths <- Filter(function(x) any(grepl("S03|S06|S09", x)), path_list)
mspIhaeIII_BS_sample_name <- sapply(mspIhaeIII_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspIhaeIII_BS_sample_name <- as.list(mspIhaeIII_BS_sample_name)

mspIhaeIII_BS_BismarkCoverage<- methRead(mspIhaeIII_BS_paths, sample.id = mspIhaeIII_BS_sample_name,  
                                        assembly = "hg38", mincov =3, 
                                        treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspIhaeIII_BS_BismarkCoverage_meth<- methylKit::unite(mspIhaeIII_BS_BismarkCoverage)
mspIhaeIII_BS_df <- as.data.frame(mspIhaeIII_BS_BismarkCoverage_meth)
mspIhaeIII_BS_CpGs <- mspIhaeIII_BS_df[mspIhaeIII_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspIhaeIII_BS_CpGs)
# 7173275       13

rm(mspIhaeIII_BS_BismarkCoverage_meth)
rm(mspIhaeIII_BS_BismarkCoverage)
rm(mspIhaeIII_BS_df)
#To GRanges for annotation
mspIhaeIII_BS_CpGs_gr <- as(mspIhaeIII_BS_CpGs,"GRanges")
seqlevelsStyle(mspIhaeIII_BS_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI + enzymatic [S10,S13,S16]
##########################################

# Select S10
S10_mspI_BS_paths <- Filter(function(x) any(grepl("S10", x)), path_list)
S10_mspI_BS_sample_name <- sapply(S10_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S10_mspI_BS_sample_name <- as.list(S10_mspI_BS_sample_name)

S10_mspI_BS_BismarkCoverage<- methRead(S10_mspI_BS_paths, sample.id = S10_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S10_mspI_BS_df <- as.data.frame(S10_mspI_BS_BismarkCoverage)

S10_mspI_BS_CpGs <- S10_mspI_BS_df[S10_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S10_mspI_BS_CpGs)
# 8998878      13
rm(S10_mspI_BS_BismarkCoverage)
rm(S10_mspI_BS_df)
#To GRanges for annotation
S10_mspI_BS_CpGs_gr <- as(S10_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S10_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S13
S13_mspI_BS_paths <- Filter(function(x) any(grepl("S13", x)), path_list)
S13_mspI_BS_sample_name <- sapply(S13_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S13_mspI_BS_sample_name <- as.list(S13_mspI_BS_sample_name)

S13_mspI_BS_BismarkCoverage<- methRead(S13_mspI_BS_paths, sample.id = S13_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S13_mspI_BS_df <- as.data.frame(S13_mspI_BS_BismarkCoverage)

S13_mspI_BS_CpGs <- S13_mspI_BS_df[S13_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S13_mspI_BS_CpGs)
# 8098828      13
rm(S13_mspI_BS_BismarkCoverage)
rm(S13_mspI_BS_df)
#To GRanges for annotation
S13_mspI_BS_CpGs_gr <- as(S13_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S13_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S16
S16_mspI_BS_paths <- Filter(function(x) any(grepl("S16", x)), path_list)
S16_mspI_BS_sample_name <- sapply(S16_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S16_mspI_BS_sample_name <- as.list(S16_mspI_BS_sample_name)

S16_mspI_BS_BismarkCoverage<- methRead(S16_mspI_BS_paths, sample.id = S16_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S16_mspI_BS_df <- as.data.frame(S16_mspI_BS_BismarkCoverage)

S16_mspI_BS_CpGs <- S16_mspI_BS_df[S16_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S16_mspI_BS_CpGs)
# 7662379      13
rm(S16_mspI_BS_BismarkCoverage)
rm(S16_mspI_BS_df)
#To GRanges for annotation
S16_mspI_BS_CpGs_gr <- as(S16_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S16_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S10,S13,S16
mspI_enzym_paths <- Filter(function(x) any(grepl("S10|S13|S16", x)), path_list)
mspI_enzym_sample_name <- sapply(mspI_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspI_enzym_sample_name <- as.list(mspI_enzym_sample_name)

mspI_enzym_BismarkCoverage<- methRead(mspI_enzym_paths, sample.id = mspI_enzym_sample_name,  
                                   assembly = "hg38", mincov = 3, 
                                   treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspI_enzym_BismarkCoverage_meth<- methylKit::unite(mspI_enzym_BismarkCoverage)
mspI_enzym_df <- as.data.frame(mspI_enzym_BismarkCoverage_meth)
mspI_enzym_CpGs <- mspI_enzym_df[mspI_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspI_enzym_CpGs)
# 5382033       13
rm(mspI_enzym_BismarkCoverage_meth)
rm(mspI_enzym_BismarkCoverage)
rm(mspI_enzym_df)
#To GRanges for annotation
mspI_enzym_CpGs_gr <- as(mspI_enzym_CpGs,"GRanges")
seqlevelsStyle(mspI_enzym_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-TaqaI + enzymatic [S11,S14,S17]
##########################################
# Select S11
S11_mspI_BS_paths <- Filter(function(x) any(grepl("S11", x)), path_list)
S11_mspI_BS_sample_name <- sapply(S11_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S11_mspI_BS_sample_name <- as.list(S11_mspI_BS_sample_name)

S11_mspI_BS_BismarkCoverage<- methRead(S11_mspI_BS_paths, sample.id = S11_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S11_mspI_BS_df <- as.data.frame(S11_mspI_BS_BismarkCoverage)

S11_mspI_BS_CpGs <- S11_mspI_BS_df[S11_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S11_mspI_BS_CpGs)
# 10437592      13
rm(S11_mspI_BS_BismarkCoverage)
rm(S11_mspI_BS_df)
#To GRanges for annotation
S11_mspI_BS_CpGs_gr <- as(S11_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S11_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S14
S14_mspI_BS_paths <- Filter(function(x) any(grepl("S14", x)), path_list)
S14_mspI_BS_sample_name <- sapply(S14_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S14_mspI_BS_sample_name <- as.list(S14_mspI_BS_sample_name)

S14_mspI_BS_BismarkCoverage<- methRead(S14_mspI_BS_paths, sample.id = S14_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S14_mspI_BS_df <- as.data.frame(S14_mspI_BS_BismarkCoverage)

S14_mspI_BS_CpGs <- S14_mspI_BS_df[S14_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S14_mspI_BS_CpGs)
# 10138771      13
rm(S14_mspI_BS_BismarkCoverage)
rm(S14_mspI_BS_df)
#To GRanges for annotation
S14_mspI_BS_CpGs_gr <- as(S14_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S14_mspI_BS_CpGs_gr) = "UCSC"  # necessary


# Select S17
S17_mspI_BS_paths <- Filter(function(x) any(grepl("S17", x)), path_list)
S17_mspI_BS_sample_name <- sapply(S17_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S17_mspI_BS_sample_name <- as.list(S17_mspI_BS_sample_name)

S17_mspI_BS_BismarkCoverage<- methRead(S17_mspI_BS_paths, sample.id = S17_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S17_mspI_BS_df <- as.data.frame(S17_mspI_BS_BismarkCoverage)

S17_mspI_BS_CpGs <- S17_mspI_BS_df[S17_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S17_mspI_BS_CpGs)
#8918220       13
rm(S17_mspI_BS_BismarkCoverage)
rm(S17_mspI_BS_df)
#To GRanges for annotation
S17_mspI_BS_CpGs_gr <- as(S17_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S17_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S11,S14,S17
mspItaqaI_enzym_paths <- Filter(function(x) any(grepl("S11|S14|S17", x)), path_list)
mspItaqaI_enzym_sample_name <- sapply(mspItaqaI_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspItaqaI_enzym_sample_name <- as.list(mspItaqaI_enzym_sample_name)

mspItaqaI_enzym_BismarkCoverage<- methRead(mspItaqaI_enzym_paths, sample.id = mspItaqaI_enzym_sample_name,  
                                      assembly = "hg38", mincov =3, 
                                      treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspItaqaI_enzym_BismarkCoverage_meth<- methylKit::unite(mspItaqaI_enzym_BismarkCoverage)
mspItaqaI_enzym_df <- as.data.frame(mspItaqaI_enzym_BismarkCoverage_meth)
mspItaqaI_enzym_CpGs <- mspItaqaI_enzym_df[mspItaqaI_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspItaqaI_enzym_CpGs)
# 6331083       13
rm(mspItaqaI_enzym_BismarkCoverage_meth)
rm(mspItaqaI_enzym_BismarkCoverage)
rm(mspItaqaI_enzym_df)
#To GRanges for annotation
mspItaqaI_enzym_CpGs_gr <- as(mspItaqaI_enzym_CpGs,"GRanges")
seqlevelsStyle(mspItaqaI_enzym_CpGs_gr) = "UCSC"  # necessary

##########################################
# Load data (.cov) from MspI-haeIII + enzymatic [S12,S15,S18]
##########################################
# Select S12
S12_mspI_BS_paths <- Filter(function(x) any(grepl("S12", x)), path_list)
S12_mspI_BS_sample_name <- sapply(S12_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S12_mspI_BS_sample_name <- as.list(S12_mspI_BS_sample_name)

S12_mspI_BS_BismarkCoverage<- methRead(S12_mspI_BS_paths, sample.id = S12_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S12_mspI_BS_df <- as.data.frame(S12_mspI_BS_BismarkCoverage)

S12_mspI_BS_CpGs <- S12_mspI_BS_df[S12_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S12_mspI_BS_CpGs)
# 11518757      13
rm(S12_mspI_BS_BismarkCoverage)
rm(S12_mspI_BS_df)
#To GRanges for annotation
S12_mspI_BS_CpGs_gr <- as(S12_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S12_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S15
S15_mspI_BS_paths <- Filter(function(x) any(grepl("S15", x)), path_list)
S15_mspI_BS_sample_name <- sapply(S15_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S15_mspI_BS_sample_name <- as.list(S15_mspI_BS_sample_name)

S15_mspI_BS_BismarkCoverage<- methRead(S15_mspI_BS_paths, sample.id = S15_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S15_mspI_BS_df <- as.data.frame(S15_mspI_BS_BismarkCoverage)

S15_mspI_BS_CpGs <- S15_mspI_BS_df[S15_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S15_mspI_BS_CpGs)
# 11329817      13
rm(S15_mspI_BS_BismarkCoverage)
rm(S15_mspI_BS_df)
#To GRanges for annotation
S15_mspI_BS_CpGs_gr <- as(S15_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S15_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S18
S18_mspI_BS_paths <- Filter(function(x) any(grepl("S18", x)), path_list)
S18_mspI_BS_sample_name <- sapply(S18_mspI_BS_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
S18_mspI_BS_sample_name <- as.list(S18_mspI_BS_sample_name)

S18_mspI_BS_BismarkCoverage<- methRead(S18_mspI_BS_paths, sample.id = S18_mspI_BS_sample_name,  
                                       assembly = "hg38",  mincov = 3, 
                                       treatment = c(1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

S18_mspI_BS_df <- as.data.frame(S18_mspI_BS_BismarkCoverage)

S18_mspI_BS_CpGs <- S18_mspI_BS_df[S18_mspI_BS_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(S18_mspI_BS_CpGs)
# 10479629      13
rm(S18_mspI_BS_BismarkCoverage)
rm(S18_mspI_BS_df)
#To GRanges for annotation
S18_mspI_BS_CpGs_gr <- as(S18_mspI_BS_CpGs,"GRanges")
seqlevelsStyle(S18_mspI_BS_CpGs_gr) = "UCSC"  # necessary

# Select S12,S15,S18
mspIhaeIII_enzym_paths <- Filter(function(x) any(grepl("S12|S15|S18", x)), path_list)
mspIhaeIII_enzym_sample_name <- sapply(mspIhaeIII_enzym_paths, function(x) {
  sub(".*(/)(S[0-9]{2})_.*", "\\2", x)
})
mspIhaeIII_enzym_sample_name <- as.list(mspIhaeIII_enzym_sample_name)

mspIhaeIII_enzym_BismarkCoverage<- methRead(mspIhaeIII_enzym_paths, sample.id = mspIhaeIII_enzym_sample_name,  
                                           assembly = "hg38", mincov = 3, 
                                           treatment = c(1,1,1), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

mspIhaeIII_enzym_BismarkCoverage_meth<- methylKit::unite(mspIhaeIII_enzym_BismarkCoverage)
mspIhaeIII_enzym_df <- as.data.frame(mspIhaeIII_enzym_BismarkCoverage_meth)
mspIhaeIII_enzym_CpGs <- mspIhaeIII_enzym_df[mspIhaeIII_enzym_df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(mspIhaeIII_enzym_CpGs)
# 8936424       13
rm(mspIhaeIII_enzym_BismarkCoverage_meth)
rm(mspIhaeIII_enzym_BismarkCoverage)
rm(mspIhaeIII_enzym_df)
#To GRanges for annotation
mspIhaeIII_enzym_CpGs_gr <- as(mspIhaeIII_enzym_CpGs,"GRanges")
seqlevelsStyle(mspIhaeIII_enzym_CpGs_gr) = "UCSC"  # necessary

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/CpGs_cov3.RData")
