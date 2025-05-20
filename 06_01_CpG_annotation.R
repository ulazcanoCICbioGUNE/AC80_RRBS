/fastdata/GPArkaitz_fastdata/R-4.4.1/bin/R
.libPaths("/fastdata/GPArkaitz_fastdata/R_fastdata/epigenomics")

#options(repos = c(CRAN = "https://cran.rstudio.com/"))  # Set a CRAN mirror
#install.packages("tidyverse", dependencies = TRUE)
#################################################################################
#                 ANNOTATING genome wide CpGs
#################################################################################

# Summary
#---------
#This script is used to annotate all the CpGs found in the human genome regarding their position
#to closest gen. 

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
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(readxl)
# Get the reference genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Find all CpG sites in genome
findCpG <- function(seqname) {
  seq <- genome[[seqname]]  # Get the sequence for a given chromosome
  cpgs <- matchPattern("CG", seq)  # Find "CG" motifs
  data.frame(
    chrom = seqname,
    start = start(cpgs),
    end = end(cpgs)
  )
}

# Run for all autosomes + X, Y
chroms <- paste0("chr", c(1:22, "X", "Y"))
#Create the data for all the CG dinucleotides in 
cpg_sites <- do.call(rbind, lapply(chroms, findCpG))
dim(cpg_sites)
#29401360        3

#To GRanges for annotation
cpgs_gr <- as(cpg_sites,"GRanges")
seqlevelsStyle(cpgs_gr) = "UCSC"  

cpgs_df <- as.data.frame(cpgs_gr)
cpgs_df$index <- paste(cpgs_df$seqnames,cpgs_df$start, sep="_")
length(unique(cpgs_df$index))
#29401360
all_indexes <- as.list(cpgs_df$index)

# Select annotations for intersection with regions in annotrt()
annots <- c("hg38_basicgenes")#, "hg38_genes_intergenic")

# Build the annotations (a single GRanges object)
annotations <- build_annotations(genome = 'hg38', annotations = annots)

#Annotate using annotrt() package the full GR
cpgs_gr_annotated <- annotate_regions(
  regions = cpgs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)
cpgs_df_annotated <- data.frame(cpgs_gr_annotated)
head(cpgs_df_annotated)
dim(cpgs_df_annotated)
# 136540315        15
#Create index column to have a unique idebntifier per CpG
cpgs_df_annotated$index <- paste(cpgs_df_annotated$seqnames,cpgs_df_annotated$start, sep="_")
annot_indexes <- as.list((unique(cpgs_df_annotated$index)))
length(annot_indexes)
#20434042
# Difference between the 29,401,360  and 20,434,042 are the ones not annotated, previously as Distal intergenic.
# Save them for the final data 
# Find the indexes in all_indexes that are not in annot_indexes
length(all_indexes)
#29401360
length(annot_indexes)
#20434042
diff_index <- setdiff(all_indexes, annot_indexes)
length(diff_index)
#8967318

# Convert diff_index to a vector if needed
diff_index <- unlist(diff_index)
# Extract rows from taqaII_df where the index is in diff_index
diff_rows <- cpgs_df[cpgs_df$index %in% diff_index, ]
# Add "Distal integrenic" category
diff_rows$prio_index <- "Distal_intergenic"

#As we do not have the biotype of each transcript I dowloaded on 2024/12/18 the ensembl database
#Add ensmbl info about biotype of each transcript
#library(biomaRt)
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#biotype_data <- getBM(
#  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "transcript_biotype"),
#  mart = ensembl
#)
#dim(biotype_data)

biotype_data <- read.csv("/vols/GPArkaitz_bigdata/ulazcano/AC80/06_ANNOTATION/biotype_data.csv")
# Create a new column in df by removing everything after the "." in annot.tx_id to easier classification
cpgs_df_annotated$ensembl_transcript_id <- sub("\\..*", "", cpgs_df_annotated$annot.tx_id)
#Create a new column with the clean genic location info
cpgs_df_annotated$type_prom <- sapply(strsplit(as.character(cpgs_df_annotated$annot.id), ":"), `[`, 1)

# Perform the join to add biotype info from ensembl to all the annotated transcripts
cpgs_df_annotated_V2 <- cpgs_df_annotated %>%
  left_join(biotype_data, by = c("ensembl_transcript_id" = "ensembl_transcript_id"))

# Create a table with the info retrieved from Ensembl and manually clasificate in three groups
biotype_data_grouped <- read.csv("/vols/GPArkaitz_bigdata/ulazcano/AC80/06_ANNOTATION/biotype_data_table.csv")
#assign group "type" from biotype data processed by me
cpgs_df_annotated_V2 <- cpgs_df_annotated_V2 %>%
  left_join(biotype_data_grouped, by = c("transcript_biotype" = "Var1"))
dim(cpgs_df_annotated_V2)
# 133389018    35
length(unique(cpgs_df_annotated_V2$index ))
#20434042

#How many trascripts with no info in type_prom.x and type_prom.y
print(paste("NA in 'type_prom.x' ", sum(is.na(cpgs_df_annotated_V2$type_prom.x))))
#0
print(paste("NA in 'type_prom.y' ", sum(is.na(cpgs_df_annotated_V2$type_prom.y))))
#809102

#na_rows <- taqaI_df_annotated_V2[is.na(taqaI_df_annotated_V2$type_prom.y), ]
#write.csv(na_rows, "/vols/GPArkaitz_bigdata/ulazcano/na_rows.csv")

#Eliminate this NA 
cpgs_df_annotated_V2_cleaned <- cpgs_df_annotated_V2[!is.na(cpgs_df_annotated_V2$type_prom.y), ]
dim(cpgs_df_annotated_V2_cleaned)
#136477098 
length(unique(cpgs_df_annotated_V2_cleaned$index))
#20466004

#First create a new category to use in the priorization including the info in the two lists
cpgs_df_annotated_V2_cleaned$prio_index <- paste(cpgs_df_annotated_V2_cleaned$type_prom.y,cpgs_df_annotated_V2_cleaned$type_prom.x, sep="_")
head(cpgs_df_annotated_V2_cleaned$prio_index)
# Create priorization list
v_type <- c("Gene", "Non_coding","Pseudogene" )
v_promoter <- c("promoter", "1to5kb","5UTR", "intron","exon", "3UTR")
# Generate all combinations
prio_list <- paste(rep(v_type, each = length(v_promoter)), v_promoter, sep = "_")

# Order the data frame by the priority order
cpgs_df_annotated_V2_cleaned$prio_index  <- factor(cpgs_df_annotated_V2_cleaned$prio_index , levels = prio_list)
cpgs_df_annotated_V2_ordered <- cpgs_df_annotated_V2_cleaned[order(cpgs_df_annotated_V2_cleaned$prio_index ), ]

# Remove duplicate CpG sites, keeping the highest priority annotation
cpgs_df_annotated_V2_out <- cpgs_df_annotated_V2_ordered[!duplicated(cpgs_df_annotated_V2_ordered$index), ]
dim(cpgs_df_annotated_V2_out)
#20466004

# Add missing columns to diff_rows with NA values
missing_cols <- setdiff(names(cpgs_df_annotated_V2_out), names(diff_rows))
diff_rows[missing_cols] <- NA

# Reorder columns of diff_rows to match cpgs_df_annotated_V2_out
diff_rows <- diff_rows[names(cpgs_df_annotated_V2_out)]
# Combine the datasets
cpgs_df_annotated__V3 <- rbind(cpgs_df_annotated_V2_out, diff_rows)
dim(cpgs_df_annotated__V3)
#29384467      37

write.csv(cpgs_df_annotated__V3, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/cpgAnotation20250313.csv")
#cp /fastdata/GPArkaitz_fastdata/ulazcano/AC80/06_ANNOTATION/genic_annotation_taqaI_9043564_20250124.csv /vols/GPArkaitz/bigdata/ulazcano/AC80/06_ANNOTATION/
cpgs_df_annotated__V3 <- fread("/fastdata/GPArkaitz_fastdata/ulazcano/AC80/cpgAnotation20250313.csv")




