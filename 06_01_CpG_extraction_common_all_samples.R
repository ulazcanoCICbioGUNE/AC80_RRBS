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
#library(bsseq)
#library(plyranges)
#library(UCSCRepeatMasker)
library(dplyr)
#library(genomation)
library(GenomicFeatures)
library(ChIPseeker)
library(data.table)
library(pheatmap)
library(matrixStats)  # for rowVars


setwd("W:/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/")
dir_infiles <- "W:/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/"

setwd("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/")
dir_infiles <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/05_Methylation_extraction/cov2cyt/"


#Read coverage files (containing cytosines only from my samples)
#List all files TIENE QUE ESTAR TODO AS LIST
file_list<-list.files(path=dir_infiles, pattern = ".CpG_report.merged_CpG_evidence.cov", all.files=FALSE,
                      full.names=FALSE)
path_list <- as.list(paste0(dir_infiles,'/', unlist(file_list)))

sample_name <- as.list(strsplit(file_list, ".CpG_report.merged_CpG_evidence.cov"))

############################################
# Load data (.cov) from all samples together
############################################
# Select all samples

BismarkCoverage<- methRead(path_list, sample.id = sample_name,  
                                   assembly = "hg38", mincov =1, 
                                   treatment = rep(1, 18), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

BismarkCoverage_meth<- methylKit::unite(BismarkCoverage)
df <- as.data.frame(BismarkCoverage_meth)
CpGs <- df[df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(CpGs)
# 3529935      13
rm(BismarkCoverage_meth)
rm(BismarkCoverage)
rm(df)
#To GRanges for annotation
CpGs_gr <- as(CpGs,"GRanges")
seqlevelsStyle(CpGs_gr) = "UCSC"  # necessary
#dim 3529935

# CALCULATE FRACTIONAL METHYLATION
# Loop through all 18 samples 
for (i in 1:18) {
  # Build column names dynamically
  numCs_col <- paste0("numCs", i)
  numTs_col <- paste0("numTs", i)
  frac_methyl_col <- paste0("frac_methyl", i)
  
  # Compute fractional methylation and add as a new metadata column
  mcols(CpGs_gr)[[frac_methyl_col]] <- mcols(CpGs_gr)[[numCs_col]] / 
    (mcols(CpGs_gr)[[numCs_col]] + mcols(CpGs_gr)[[numTs_col]])
}

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov1_3529935.RData")
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov1_3529935.RData")
#load("W:/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov1_3529935.RData")

CpGs_mincov1 <- as.data.frame(CpGs_gr)
CpGs_mincov1$index <- paste(CpGs_mincov1$seqnames,CpGs_mincov1$start, sep="_")


load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov10_1238930.RData")
#load("W:/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov10_1238930.RData")
CpGs_mincov10 <- as.data.frame(CpGs_gr)
CpGs_mincov10$index <- paste(CpGs_mincov10$seqnames,CpGs_mincov10$start, sep="_")

head(CpGs_mincov10)

CpGs_mincov1_common_CpG <- CpGs_mincov1[CpGs_mincov1$index %in% CpGs_mincov10$index, ]
head(CpGs_mincov1_common_CpG)
#Plot of the results to see how do fractional methylation values correlate in the same CpGs between the two techniques

# Separate the samples in the two grous regarding their treatment
group1_cols <- paste0("frac_methyl", 1:9)
group2_cols <- paste0("frac_methyl", 10:18)

# Extract the data matrices
group1 <- as.matrix(CpGs_mincov1_common_CpG[, group1_cols])
group2 <- as.matrix(CpGs_mincov1_common_CpG[, group2_cols])

# Compute pairwise correlation: each column in group1 vs each in group2
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs")
# Set custom row and column names
rownames(cor_matrix) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(cor_matrix) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")
# Plot heatmap and save to PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/pairwise_correlation_common_CpGs_cov1_3529935.pdf", width = 8, height = 6)
pheatmap(cor_matrix,
         main = "Pairwise Correlation: BS vs Enzym",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Extract methylation fraction columns for both groups
group1 <- CpGs_mincov1_common_CpG[, paste0("frac_methyl", 1:9)]
group2 <- CpGs_mincov1_common_CpG[, paste0("frac_methyl", 10:18)]

# Calculate variance per CpG for group1 (can also do for combined if preferred)
vars <- rowVars(as.matrix(group1))

# Get indices for top 500 and bottom 500 by variance
top_idx <- order(vars, decreasing = TRUE)[1:500]
bottom_idx <- order(vars, decreasing = FALSE)[1:500]

# Combine indices
subset_idx <- c(top_idx, bottom_idx)

# Subset data matrices
group1_sub <- as.matrix(group1)[subset_idx, ]
group2_sub <- as.matrix(group2)[subset_idx, ]

# Manually rename columns with samples name , enzyme and treatment
colnames(group1_sub) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(group2_sub) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")
# Save group1 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/BS_heatmap_common_CpGs_cov1_3529935.pdf", width = 8, height = 12)
pheatmap(group1_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - BS", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

# Save group2 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/Enzym_heatmap_common_CpGs_cov1_3529935.pdf", width = 8, height = 12)
pheatmap(group2_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - Enzym", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()


############################################
# Load data (.cov) from all samples together
############################################
# Select all samples
BismarkCoverage<- methRead(path_list, sample.id = sample_name,  
                           assembly = "hg38", mincov =3, 
                           treatment = rep(1, 18), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

BismarkCoverage_meth<- methylKit::unite(BismarkCoverage)
df <- as.data.frame(BismarkCoverage_meth)
CpGs <- df[df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(CpGs)
# 2242268      13
rm(BismarkCoverage_meth)
rm(BismarkCoverage)
rm(df)
#To GRanges for annotation
CpGs_gr <- as(CpGs,"GRanges")
seqlevelsStyle(CpGs_gr) = "UCSC"  # necessary
#dim 2242268

# CALCULATE FRACTIONAL METHYLATION
# Loop through all 18 samples 
for (i in 1:18) {
  # Build column names dynamically
  numCs_col <- paste0("numCs", i)
  numTs_col <- paste0("numTs", i)
  frac_methyl_col <- paste0("frac_methyl", i)
  
  # Compute fractional methylation and add as a new metadata column
  mcols(CpGs_gr)[[frac_methyl_col]] <- mcols(CpGs_gr)[[numCs_col]] / 
    (mcols(CpGs_gr)[[numCs_col]] + mcols(CpGs_gr)[[numTs_col]])
}

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov3_2242268.RData")
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov3_2242268.RData")
CpGs_mincov3 <- as.data.frame(CpGs_gr)
CpGs_mincov3$index <- paste(CpGs_mincov3$seqnames,CpGs_mincov3$start, sep="_")


CpGs_mincov3_common_CpG <- CpGs_mincov3[CpGs_mincov3$index %in% CpGs_mincov10$index, ]
head(CpGs_mincov3_common_CpG)


# Select sample columns
group1_cols <- paste0("frac_methyl", 1:9)
group2_cols <- paste0("frac_methyl", 10:18)

# Extract the data matrices
group1 <- as.matrix(CpGs_mincov3_common_CpG[, group1_cols])
group2 <- as.matrix(CpGs_mincov3_common_CpG[, group2_cols])

# Compute pairwise correlation: each column in group1 vs each in group2
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs")

# Set custom row and column names
rownames(cor_matrix) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(cor_matrix) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")

# Plot heatmap and save to PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/pairwise_correlation_common_CpGs_cov3_2242268.pdf", width = 8, height = 6)
pheatmap(cor_matrix,
         main = "Pairwise Correlation: BS vs Enzym",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Extract methylation fraction columns for both groups
group1 <- CpGs_mincov3_common_CpG[, paste0("frac_methyl", 1:9)]
group2 <- CpGs_mincov3_common_CpG[, paste0("frac_methyl", 10:18)]

# Calculate variance per CpG for group1 (can also do for combined if preferred)
vars <- rowVars(as.matrix(group1))

# Get indices for top 500 and bottom 500 by variance
top_idx <- order(vars, decreasing = TRUE)[1:500]
bottom_idx <- order(vars, decreasing = FALSE)[1:500]

# Combine indices
subset_idx <- c(top_idx, bottom_idx)

# Subset data matrices
group1_sub <- as.matrix(group1)[subset_idx, ]
group2_sub <- as.matrix(group2)[subset_idx, ]

# Manually rename columns with samples name , enzyme and treatment
colnames(group1_sub) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(group2_sub) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")
# Save group1 heatmap as PDF

# Save group1 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/BS_heatmap_common_CpGs_cov3_2242268.pdf", width = 8, height = 12)
pheatmap(group1_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - BS", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

# Save group2 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/Enzym_heatmap_common_CpGs_cov3_2242268.pdf", width = 8, height = 12)
pheatmap(group2_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - Enzym", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

############################################
# Load data (.cov) from all samples together
############################################
# Select all samples
BismarkCoverage<- methRead(path_list, sample.id = sample_name,  
                           assembly = "hg38", mincov =5, 
                           treatment = rep(1, 18), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

BismarkCoverage_meth<- methylKit::unite(BismarkCoverage)
df <- as.data.frame(BismarkCoverage_meth)
CpGs <- df[df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(CpGs)
# 1789918      13
rm(BismarkCoverage_meth)
rm(BismarkCoverage)
rm(df)
#To GRanges for annotation
CpGs_gr <- as(CpGs,"GRanges")
seqlevelsStyle(CpGs_gr) = "UCSC"  # necessary
#dim 1789918

# CALCULATE FRACTIONAL METHYLATION
# Loop through all 18 samples 
for (i in 1:18) {
  # Build column names dynamically
  numCs_col <- paste0("numCs", i)
  numTs_col <- paste0("numTs", i)
  frac_methyl_col <- paste0("frac_methyl", i)
  
  # Compute fractional methylation and add as a new metadata column
  mcols(CpGs_gr)[[frac_methyl_col]] <- mcols(CpGs_gr)[[numCs_col]] / 
    (mcols(CpGs_gr)[[numCs_col]] + mcols(CpGs_gr)[[numTs_col]])
}

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov5_1789918.RData")
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov5_1789918.RData")
CpGs_mincov5 <- as.data.frame(CpGs_gr)
CpGs_mincov5$index <- paste(CpGs_mincov5$seqnames,CpGs_mincov5$start, sep="_")

CpGs_mincov5_common_CpG <- CpGs_mincov5[CpGs_mincov5$index %in% CpGs_mincov10$index, ]
head(CpGs_mincov5_common_CpG)

# Select sample columns
group1_cols <- paste0("frac_methyl", 1:9)
group2_cols <- paste0("frac_methyl", 10:18)

# Extract the data matrices
group1 <- as.matrix(CpGs_mincov5_common_CpG[, group1_cols])
group2 <- as.matrix(CpGs_mincov5_common_CpG[, group2_cols])

# Compute pairwise correlation: each column in group1 vs each in group2
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs")
# Set custom row and column names
rownames(cor_matrix) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(cor_matrix) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")
# Plot heatmap and save to PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/pairwise_correlation_common_CpGs_cov5_1789918.pdf", width = 8, height = 6)
pheatmap(cor_matrix,
         main = "Pairwise Correlation: BS vs Enzym",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Extract methylation fraction columns for both groups
group1 <- CpGs_mincov5_common_CpG[, paste0("frac_methyl", 1:9)]
group2 <- CpGs_mincov5_common_CpG[, paste0("frac_methyl", 10:18)]

# Calculate variance per CpG for group1 (can also do for combined if preferred)
vars <- rowVars(as.matrix(group1))

# Get indices for top 500 and bottom 500 by variance
top_idx <- order(vars, decreasing = TRUE)[1:500]
bottom_idx <- order(vars, decreasing = FALSE)[1:500]

# Combine indices
subset_idx <- c(top_idx, bottom_idx)

# Subset data matrices
group1_sub <- as.matrix(group1)[subset_idx, ]
group2_sub <- as.matrix(group2)[subset_idx, ]

# Manually rename columns with samples name , enzyme and treatment
colnames(group1_sub) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(group2_sub) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")


# Save group1 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/BS_heatmap_common_CpGs_cov5_1789918.pdf", width = 8, height = 12)
pheatmap(group1_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - BS", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

# Save group2 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/Enzym_heatmap_common_CpGs_cov5_1789918.pdf", width = 8, height = 12)
pheatmap(group2_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - Enzym", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

############################################
# Load data (.cov) from all samples together
############################################
# Select all samples
BismarkCoverage<- methRead(path_list, sample.id = sample_name,  
                           assembly = "hg38", mincov =10, 
                           treatment = rep(1, 18), context = "CpG", resolution ="base", pipeline = "bismarkCoverage")

BismarkCoverage_meth<- methylKit::unite(BismarkCoverage)
df <- as.data.frame(BismarkCoverage_meth)
CpGs <- df[df$chr%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"), ]
dim(CpGs)
# 1238930      13
rm(BismarkCoverage_meth)
rm(BismarkCoverage)
rm(df)
#To GRanges for annotation
CpGs_gr <- as(CpGs,"GRanges")
seqlevelsStyle(CpGs_gr) = "UCSC"  # necessary
#dim 1238930

# CALCULATE FRACTIONAL METHYLATION
# Loop through all 18 samples 
for (i in 1:18) {
  # Build column names dynamically
  numCs_col <- paste0("numCs", i)
  numTs_col <- paste0("numTs", i)
  frac_methyl_col <- paste0("frac_methyl", i)
  
  # Compute fractional methylation and add as a new metadata column
  mcols(CpGs_gr)[[frac_methyl_col]] <- mcols(CpGs_gr)[[numCs_col]] / 
    (mcols(CpGs_gr)[[numCs_col]] + mcols(CpGs_gr)[[numTs_col]])
}

save.image(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov10_1238930.RData")
load("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/common_CpGs_cov10_1238930.RData")
CpGs <- as.data.frame(CpGs_gr)


# Select sample columns
group1_cols <- paste0("frac_methyl", 1:9)
group2_cols <- paste0("frac_methyl", 10:18)

# Extract the data matrices
group1 <- as.matrix(CpGs[, group1_cols])
group2 <- as.matrix(CpGs[, group2_cols])

# Compute pairwise correlation: each column in group1 vs each in group2
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs")
# Set custom row and column names
rownames(cor_matrix) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(cor_matrix) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")
# Plot heatmap and save to PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/pairwise_correlation_common_CpGs_cov10_1238930.pdf", width = 8, height = 6)
pheatmap(cor_matrix,
         main = "Pairwise Correlation: BS vs Enzym ",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Extract methylation fraction columns for both groups
group1 <- CpGs[, paste0("frac_methyl", 1:9)]
group2 <- CpGs[, paste0("frac_methyl", 10:18)]

# Calculate variance per CpG for group1 (can also do for combined if preferred)
vars <- rowVars(as.matrix(group1))

# Get indices for top 500 and bottom 500 by variance
top_idx <- order(vars, decreasing = TRUE)[1:500]
bottom_idx <- order(vars, decreasing = FALSE)[1:500]

# Combine indices
subset_idx <- c(top_idx, bottom_idx)

# Subset data matrices
group1_sub <- as.matrix(group1)[subset_idx, ]
group2_sub <- as.matrix(group2)[subset_idx, ]

# Manually rename columns with samples name , enzyme and treatment
colnames(group1_sub) <- c("S_22_MspI_BS", "S_22_TaqaI_BS", "S_22_HaeIII_BS", "S_27_MspI_BS", "S_27_TaqaI_BS", 
                          "S_27_HaeIII_BS", "S_31_MspI_BS", "S_31_TaqaI_BS", "S_31_HaeIII_BS")

colnames(group2_sub) <- c("S_22_MspI_Enzym", "S_22_TaqaI_Enzym", "S_22_HaeIII_Enzym", "S_27_MspI_Enzym", "S_27_TaqaI_Enzym", 
                          "S_27_HaeIII_Enzym", "S_31_MspI_Enzym", "S_31_TaqaI_Enzym", "S_31_HaeIII_Enzym")

# Save group1 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/BS_heatmap_common_CpGs_cov10_1238930.pdf", width = 8, height = 12)
pheatmap(group1_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - BS", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

# Save group2 heatmap as PDF
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/Enzym_heatmap_common_CpGs_cov10_1238930.pdf", width = 8, height = 12)
pheatmap(group2_sub, 
         main = "Top 500 and Bottom 500 Variable CpGs - Enzym", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE)
dev.off()

