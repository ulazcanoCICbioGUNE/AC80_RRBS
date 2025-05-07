#######################################
#ANNOTATION REGARDING THE GENIC POSITION
#######################################
#2025/01/23
#Uxue Lazkano

###Summary#####
# This script is used to assign a unique annotation to each CpG.

#R in lamarr
R
#.libPaths("/home/CICBIOGUNE/ulazcano/R_lamarr/epigenomics")
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")


##################
# Libraries
##################
options(repos = c(CRAN = "https://cran.r-project.org"))

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

setwd("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/")
dir_infiles <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/"

#########################################################
#Load extracted CpGs and annotated CpGs
########################################################
#Samples CpG data
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/AC80_RRBS.RData")
load("W:/ulazcano/AC80_RRBS/AC80_RRBS.RData")

#Genome wide CpG anotation
cpgs_df_annotated <- fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/cpgAnotation20250313.csv")
cpgs_df_annotated <- fread("W:/ulazcano/AC80_RRBS/06_Annotation/cpgAnotation20250313.csv")

###########
# mspI_BS
###########
head(cpgs_df_annotated_prueba)
mspI_BS_CpGs <- getData(mspI_BS_CpGs)
mspI_BS_CpGs$index <- paste0("chr", mspI_BS_CpGs$chr, "_", mspI_BS_CpGs$start)

merged_mspI_BS <- merge(mspI_BS_CpGs, cpgs_df_annotated, by = "index")
write.csv(merged_mspI_BS, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspI_BS_geneAnnotation20250507.csv")

##############
# mspItaqaI_BS
##############
head(cpgs_df_annotated_prueba)
mspI_BS_CpGs <- getData(mspItaqaI_BS_CpGs)
mspI_BS_CpGs$index <- paste0("chr", mspItaqaI_BS_CpGs$chr, "_", mspItaqaI_BS_CpGs$start)

merged_mspItaqaI_BS <- merge(mspItaqaI_BS_CpGs, cpgs_df_annotated, by = "index")
write.csv(merged_mspItaqaI_BS, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspItaqaI_BS_geneAnnotation20250507.csv")

#################
# mspIhaeIII_BS
#################
head(cpgs_df_annotated_prueba)
mspIhaeIII_BS_CpGs <- getData(mspIhaeIII_BS_CpGs)
mspIhaeIII_BS_CpGs$index <- paste0("chr", mspIhaeIII_BS_CpGs$chr, "_", mspIhaeIII_BS_CpGs$start)

merged_mspIhaeIII_BS <- merge(mspIhaeIII_BS_CpGs, cpgs_df_annotated, by = "index")
write.csv(mspIhaeIII_BS, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspIhaeIII_BS_geneAnnotation20250507.csv")

###########
# mspI_Enzym
###########
head(cpgs_df_annotated_prueba)
mspI_Enzym_CpGs <- getData(mspI_Enzym_CpGs)
mspI_Enzym_CpGs$index <- paste0("chr", mspI_Enzym_CpGs$chr, "_", mspI_Enzym_CpGs$start)

merged_mspI_Enzym <- merge(mspI_Enzym_CpGs, cpgs_df_annotated, by = "index")
write.csv(merged_mspI_Enzym, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspI_Enzym_geneAnnotation20250507.csv")

##############
# mspItaqaI_Enzym
##############
head(cpgs_df_annotated_prueba)
mspI_Enzym_CpGs <- getData(mspItaqaI_Enzym_CpGs)
mspI_Enzym_CpGs$index <- paste0("chr", mspItaqaI_Enzym_CpGs$chr, "_", mspItaqaI_Enzym_CpGs$start)

merged_mspItaqaI_Enzym <- merge(mspItaqaI_Enzym_CpGs, cpgs_df_annotated, by = "index")
write.csv(merged_mspItaqaI_Enzym, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspItaqaI_Enzym_geneAnnotation20250507.csv")

#################
# mspIhaeIII_Enzym
#################
head(cpgs_df_annotated_prueba)
mspIhaeIII_Enzym_CpGs <- getData(mspIhaeIII_Enzym_CpGs)
mspIhaeIII_Enzym_CpGs$index <- paste0("chr", mspIhaeIII_Enzym_CpGs$chr, "_", mspIhaeIII_Enzym_CpGs$start)

merged_mspIhaeIII_Enzym <- merge(mspIhaeIII_Enzym_CpGs, cpgs_df_annotated, by = "index")
write.csv(mspIhaeIII_Enzym, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/mspIhaeIII_Enzym_geneAnnotation20250507.csv")


