#############################################
#ANNOTATION USING CHROMHMM DATA
#############################################
#2025/01/23
#Uxue Lazkano

###Summary#####
# Annotation of my CpGs to ChromHMM data regardin 18 state model coming from 98 genomes of 
# Road Epigenomics data. website: https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
# Ref paper: https://www.nature.com/articles/nature14248


#R in lamarr
R
#.libPaths("/home/CICBIOGUNE/ulazcano/R_lamarr/epigenomics")
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")


##################
# Libraries
##################
library(methylKit)
library(data.table)
library(annotatr)
#library(bsseq)
#library(plyranges)
library(UCSCRepeatMasker)
library(dplyr)
library(genomation)
library(GenomicFeatures)
library(ChIPseeker)
library(readxl)
library(rtracklayer)
library(GenomicRanges)   # For working with genomic ranges and overlaps (e.g., findOverlaps, queryHits, subjectHits)
library(rtracklayer)     # For importing BED files (e.g., import() function in the safe_import function)
library(magrittr)        # For the pipe operator (%>%)
library(IRanges) 

###############################################################################
######################
# Import ChromHMM data
######################
#Read excel containing information of the clasification of the 18 states (made manually by Uxue)
roadmap_metadata <- read_excel("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/Annotation_files/ChromHMM/Road_epigenomics.xlsx")

state_clasification <- read_excel("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/Annotation_files/ChromHMM/18state_metadata.xlsx")
state_clasification$name <- paste(state_clasification$`STATE NO.`, state_clasification$MNEMONIC, sep = "_")

#Metadata from ChromHMM data available 
metadata_roadmap<-read_excel("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/Annotation_files/ChromHMM/Metadata_summary.xlsx")
cell_types <- metadata_roadmap$`Epigenome ID (EID)`
## Dowload systematically all the data from web
# Base URL for ChromHMM data
base_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"

# Function to safely import BED files from the url, handling both warnings and 
# errors because in the case of the 18 state model some of the codes don't have
# all the epigenetic mark info.
safe_import <- function(file_url) {
  result <- NULL
  tryCatch(
    {
      # Capture warnings as errors
      withCallingHandlers(
        result <- import(BEDFile(file_url)),
        warning = function(w) {
          message("Warning while importing file:", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      message("Error importing file:", conditionMessage(e))
      result <- NULL
    }
  )
  return(result)
}


##########################################
# Annotation of CpGs in all genome
##########################################

#Genome wide CpGs extracted and annotated in 06_01_CpG_annotation
cpgs_df_annotated <- fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/cpgAnotation20250313.csv")
dim(cpgs_df_annotated)
#29392464       28


# To GRanges
cpgs_gr <-  makeGRangesFromDataFrame(cpgs_df_annotated,
                                        keep.extra.columns=TRUE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames"),
                                        start.field="start",
                                        end.field=c("end"),
                                        starts.in.df.are.0based=FALSE)
#dim(mspI_BS_gr)
#7780133
seqlevelsStyle(cpgs_gr) = "UCSC"  # necessary


#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(cpgs_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(cpgs_gr))) {
      mcols(cpgs_gr)[[cell_type]] <- rep(NA, length(cpgs_gr))
      mcols(cpgs_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(cpgs_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(cpgs_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(cpgs_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}

#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
cpgs_df <- as.data.frame(mcols(cpgs_gr))

head(cpgs_df)
# Define the mapping of categories to numeric values

category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
cpgs_df <- create_numeric_columns(cpgs_df)
head(cpgs_df)

# !!!!!!!!!!!!!!!!!!!!!!!!
#De aqui pabajo revisar


############################################################
#### Create a column with the minimun value for each CpG
############################################################

# Select only columns that end with "_num"
num_cols <- grep("_num$", names(haeIII_df), value = TRUE)
haeIII_df <- as.data.frame(haeIII_df)

# Create a new column "min_cat" with the minimum value across these columns
haeIII_df$min_cat <- apply(haeIII_df[num_cols], 1, min)

# Print the updated dataframe
head(haeIII_df)

############################################################
#### Create a column to mark as enhancer and another with 
#the proportion from the total of enhancer states
############################################################
haeIII_df$enhancer <- ifelse(haeIII_df$min_cat == 1, 1, 0)
table(haeIII_df$enhancer)

#This part is to create the proportion column 
# Select only columns that end with "_num"
num_cols <- grep("_num$", names(haeIII_df), value = TRUE)

# Calculate the proportion of columns with value 1
haeIII_df$enhancer_prop <- apply(haeIII_df[num_cols], 1, function(x) sum(x == 1, na.rm = TRUE) / sum(!is.na(x)))

# Print the updated dataframe
head(haeIII_df)



############################################################
#### Create a column to mark the active cpgs  and another with 
#the proportion from the total of active states
############################################################
# Create "active" column: 1 if min_cat is in (1,2,4,5), otherwise 0
haeIII_df$active <- ifelse(haeIII_df$min_cat %in% c(1, 4, 5), 1, 0)

# Calculate "active_prop": proportion of _num columns with values in (1,2,4,5)
haeIII_df$active_prop <- apply(haeIII_df[num_cols], 1, function(x) sum(x %in% c(1, 4, 5), na.rm = TRUE) / sum(!is.na(x)))

# Print the updated dataframe
head(haeIII_df)

############################################################
#### Create a column to mark the inactive cpgs  and another with 
#the proportion from the total of inactive states
############################################################
# Create "active" column: 1 if min_cat is in (1,2,4,5), otherwise 0
haeIII_df$inactive <- ifelse(haeIII_df$min_cat %in% c(3, 6, 7), 1, 0)

# Calculate "active_prop": proportion of _num columns with values in (1,2,4,5)
haeIII_df$inactive_prop <- apply(haeIII_df[num_cols], 1, function(x) sum(x %in% c(3, 6, 7), na.rm = TRUE) / sum(!is.na(x)))

# Print the updated dataframe
head(haeIII_df)

############################################################
#### Create a column to mark the inactive cpgs  and another with 
#the proportion from the total of inactive states
############################################################
# Create "active" column: 1 if min_cat is in (1,2,4,5), otherwise 0
haeIII_df$bivalent <- ifelse(haeIII_df$min_cat %in% c(2), 1, 0)

# Calculate "active_prop": proportion of _num columns with values in (1,2,4,5)
haeIII_df$bivalent_prop <- apply(haeIII_df[num_cols], 1, function(x) sum(x %in% c(2), na.rm = TRUE) / sum(!is.na(x)))

# Print the updated dataframe
head(haeIII_df)
write.csv(haeIII_df,"/fastdata/GPArkaitz_fastdata/ulazcano/AC80/ChomHMM_processed_haeIII_20250307.csv" )





# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(cpgs_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_all_cpgs_20250526.csv")



##########################################
# Annotation of  mspI_BS
##########################################
mspI_BS <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspI_BS_7780133_20250519.csv"))
any(is.na(mspI_BS$seqnames))
any(is.na(mspI_BS$start))
any(is.na(mspI_BS$end))

# To GRanges
mspI_BS_gr <-  makeGRangesFromDataFrame(mspI_BS,
                                        keep.extra.columns=TRUE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames"),
                                        start.field="start",
                                        end.field=c("end"),
                                        starts.in.df.are.0based=FALSE)
#dim(mspI_BS_gr)
#7780133
seqlevelsStyle(mspI_BS_gr) = "UCSC"  # necessary

###### Merge annotations of all CpGs with index column
mspI_BS_CpGs$index <- paste0("chr", mspI_BS_CpGs$chr, "_", mspI_BS_CpGs$start)

merged_mspI_BS <- merge(mspI_BS_CpGs, cpgs_df_annotated, by = "index", all.x=T)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspI_BS_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_mspI_BS_7780133_20250519.csv")

##########################################
# Annotation of  mspItaqaI_BS
##########################################
mspItaqaI_BS <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspItaqaI_BS_8519467_20250519.csv"))
any(is.na(mspItaqaI_BS$seqnames))
any(is.na(mspItaqaI_BS$start))
any(is.na(mspItaqaI_BS$end))

# To GRanges
mspItaqaI_BS_gr <-  makeGRangesFromDataFrame(mspItaqaI_BS,
                                        keep.extra.columns=TRUE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames"),
                                        start.field="start",
                                        end.field=c("end"),
                                        starts.in.df.are.0based=FALSE)
#dim(mspItaqaI_BS_gr)
#7780133
seqlevelsStyle(mspItaqaI_BS_gr) = "UCSC"  # necessary

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(mspItaqaI_BS_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(mspItaqaI_BS_gr))) {
      mcols(mspItaqaI_BS_gr)[[cell_type]] <- rep(NA, length(mspItaqaI_BS_gr))
      mcols(mspItaqaI_BS_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(mspItaqaI_BS_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(mspItaqaI_BS_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(mspItaqaI_BS_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
mspItaqaI_BS_df <- as.data.frame(mcols(mspItaqaI_BS_gr))

head(mspItaqaI_BS_df)
# Define the mapping of categories to numeric values
#category_mapping <- c(
#  "Transcription" = 1,
#  "Promoter" = 2,
#  "Enhancer" = 3,
#  "Bivalent" = 4,
#  "Quiescent" = 5,
#  "Polycomb-Repressed" = 6,
#  "Repetitive-Het" = 7
#)
category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
mspItaqaI_BS_df <- create_numeric_columns(mspItaqaI_BS_df)
head(mspItaqaI_BS_df)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspItaqaI_BS_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/chromhmm_annotation_mspItaqaI_BS_8519467_20250519.csv")


##########################################
# Annotation of  mspIhaeIII_BS
##########################################
mspIhaeIII_BS <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspIhaeIII_BS_10719257_20250519.csv"))
any(is.na(mspIhaeIII_BS$seqnames))
any(is.na(mspIhaeIII_BS$start))
any(is.na(mspIhaeIII_BS$end))

# To GRanges
mspIhaeIII_BS_gr <-  makeGRangesFromDataFrame(mspIhaeIII_BS,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("seqnames"),
                                             start.field="start",
                                             end.field=c("end"),
                                             starts.in.df.are.0based=FALSE)
#dim(mspIhaeIII_BS_gr)
#7780133
seqlevelsStyle(mspIhaeIII_BS_gr) = "UCSC"  # necessary

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(mspIhaeIII_BS_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(mspIhaeIII_BS_gr))) {
      mcols(mspIhaeIII_BS_gr)[[cell_type]] <- rep(NA, length(mspIhaeIII_BS_gr))
      mcols(mspIhaeIII_BS_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(mspIhaeIII_BS_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(mspIhaeIII_BS_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(mspIhaeIII_BS_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
mspIhaeIII_BS_df <- as.data.frame(mcols(mspIhaeIII_BS_gr))

head(mspIhaeIII_BS_df)
# Define the mapping of categories to numeric values
#category_mapping <- c(
#  "Transcription" = 1,
#  "Promoter" = 2,
#  "Enhancer" = 3,
#  "Bivalent" = 4,
#  "Quiescent" = 5,
#  "Polycomb-Repressed" = 6,
#  "Repetitive-Het" = 7
#)
category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
mspIhaeIII_BS_df <- create_numeric_columns(mspIhaeIII_BS_df)
head(mspIhaeIII_BS_df)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspIhaeIII_BS_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_mspIhaeIII_BS_10719257_20250521.csv")

##########################################
# Annotation of mspI_enzym
##########################################
mspI_enzym <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspI_Enzym_9532433_20250519.csv"))
any(is.na(mspI_enzym$seqnames))
any(is.na(mspI_enzym$start))
any(is.na(mspI_enzym$end))

# To GRanges
mspI_enzym_gr <-  makeGRangesFromDataFrame(mspI_enzym,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field=c("seqnames"),
                                              start.field="start",
                                              end.field=c("end"),
                                              starts.in.df.are.0based=FALSE)
#dim(mspI_enzym_gr)
#7780133
seqlevelsStyle(mspI_enzym_gr) = "UCSC"  # necessary

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(mspI_enzym_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(mspI_enzym_gr))) {
      mcols(mspI_enzym_gr)[[cell_type]] <- rep(NA, length(mspI_enzym_gr))
      mcols(mspI_enzym_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(mspI_enzym_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(mspI_enzym_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(mspI_enzym_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
mspI_enzym_df <- as.data.frame(mcols(mspI_enzym_gr))

head(mspI_enzym_df)
# Define the mapping of categories to numeric values
#category_mapping <- c(
#  "Transcription" = 1,
#  "Promoter" = 2,
#  "Enhancer" = 3,
#  "Bivalent" = 4,
#  "Quiescent" = 5,
#  "Polycomb-Repressed" = 6,
#  "Repetitive-Het" = 7
#)
category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
mspI_enzym_df <- create_numeric_columns(mspI_enzym_df)
head(mspI_enzym_df)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspI_enzym_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_mspI_Enzym_9532433_20250521.csv")


##########################################
# Annotation of mspItaqaI_enzym
##########################################
mspItaqaI_enzym <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspItaqaI_enzym_9532433_20250519.csv"))
any(is.na(mspItaqaI_enzym$seqnames))
any(is.na(mspItaqaI_enzym$start))
any(is.na(mspItaqaI_enzym$end))

# To GRanges
mspItaqaI_enzym_gr <-  makeGRangesFromDataFrame(mspItaqaI_enzym,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("seqnames"),
                                           start.field="start",
                                           end.field=c("end"),
                                           starts.in.df.are.0based=FALSE)
#dim(mspItaqaI_enzym_gr)
#7780133
seqlevelsStyle(mspItaqaI_enzym_gr) = "UCSC"  # necessary

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(mspItaqaI_enzym_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(mspItaqaI_enzym_gr))) {
      mcols(mspItaqaI_enzym_gr)[[cell_type]] <- rep(NA, length(mspItaqaI_enzym_gr))
      mcols(mspItaqaI_enzym_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(mspItaqaI_enzym_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(mspItaqaI_enzym_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(mspItaqaI_enzym_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
mspItaqaI_enzym_df <- as.data.frame(mcols(mspItaqaI_enzym_gr))

head(mspItaqaI_enzym_df)
# Define the mapping of categories to numeric values
#category_mapping <- c(
#  "Transcription" = 1,
#  "Promoter" = 2,
#  "Enhancer" = 3,
#  "Bivalent" = 4,
#  "Quiescent" = 5,
#  "Polycomb-Repressed" = 6,
#  "Repetitive-Het" = 7
#)
category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
mspItaqaI_enzym_df <- create_numeric_columns(mspItaqaI_enzym_df)
head(mspItaqaI_enzym_df)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspItaqaI_enzym_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_mspItaqaI_enzym_9532433_20250521.csv")


##########################################
# Annotation of mspIhaeIII_enzym
##########################################
mspIhaeIII_enzym <- as.data.frame(fread("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/cpgIsland_annotation_mspIhaeIII_Enzym_11900623_20250519.csv"))
any(is.na(mspIhaeIII_enzym$seqnames))
any(is.na(mspIhaeIII_enzym$start))
any(is.na(mspIhaeIII_enzym$end))

# To GRanges
mspIhaeIII_enzym_gr <-  makeGRangesFromDataFrame(mspIhaeIII_enzym,
                                                keep.extra.columns=TRUE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c("seqnames"),
                                                start.field="start",
                                                end.field=c("end"),
                                                starts.in.df.are.0based=FALSE)
#dim(mspIhaeIII_enzym_gr)
#7780133
seqlevelsStyle(mspIhaeIII_enzym_gr) = "UCSC"  # necessary

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(mspIhaeIII_enzym_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(mspIhaeIII_enzym_gr))) {
      mcols(mspIhaeIII_enzym_gr)[[cell_type]] <- rep(NA, length(mspIhaeIII_enzym_gr))
      mcols(mspIhaeIII_enzym_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(mspIhaeIII_enzym_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(mspIhaeIII_enzym_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(mspIhaeIII_enzym_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
mspIhaeIII_enzym_df <- as.data.frame(mcols(mspIhaeIII_enzym_gr))

head(mspIhaeIII_enzym_df)
# Define the mapping of categories to numeric values
#category_mapping <- c(
#  "Transcription" = 1,
#  "Promoter" = 2,
#  "Enhancer" = 3,
#  "Bivalent" = 4,
#  "Quiescent" = 5,
#  "Polycomb-Repressed" = 6,
#  "Repetitive-Het" = 7
#)
category_mapping <- c(
  "Enhancer" = 1,
  "Bivalent" = 2,
  "Polycomb-Repressed" = 3,
  "Transcription" = 4,
  "Promoter" = 5,
  "Quiescent" = 6,
  "Repetitive-Het" = 7
)
# Function to create new numeric columns based on category columns
create_numeric_columns <- function(df) {
  for (col in names(df)) {
    if (grepl("_cat$", col)) {
      new_col <- sub("_cat$", "_num", col)
      df[[new_col]] <- as.numeric(factor(df[[col]], levels = names(category_mapping), labels = category_mapping))
    }
  }
  return(df)
}

# Apply the function to your dataframe
mspIhaeIII_enzym_df <- create_numeric_columns(mspIhaeIII_enzym_df)
head(mspIhaeIII_enzym_df)

# Nomeclatura artxibuak izendatzeko _annotation_haeIII_20232718_20250124.csv
write.csv(mspIhaeIII_enzym_df, "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_ANNOTATION/chromhmm_annotation_mspIhaeIII_Enzym_11900623_202505.csv")


################################################################################
#                        EPIC -  Data loading
################################################################################
EPICV2<-read.csv("W:/ulazcano/AC83/Manifest/MethylationEPICV2_manifest.csv", row.names=1)
#EPICV2<-read.csv("/vols/GPArkaitz_bigdata/ulazcano/AC83/Manifest/MethylationEPICV2_manifest.csv", row.names=1)

#Remove rows with no MAPINFO
EPICV2_clean <- EPICV2[!is.na(EPICV2$MAPINFO),]

EPICV2_gr <-  makeGRangesFromDataFrame(EPICV2_clean,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqinfo=NULL,
                                            seqnames.field=c("CHR"),
                                            start.field="MAPINFO",
                                            end.field="MAPINFO",
                                            starts.in.df.are.0based=FALSE)


#Read excel containing information of the clasification of the 18 states (made manually by Uxue)
state_clasification <- read_excel("C:/Users/ulazcano/CIC bioGUNE/Arkaitz group - Documentos/Individual folders/Uxue Lazkano/PhD/3.Projects/AC80_RRBS_Test/Annotation/ChromHMM/18states/18state_metadata.xlsx")
state_clasification$name <- paste(state_clasification$`STATE NO.`, state_clasification$MNEMONIC, sep = "_")

#Metadata from ChromHMM data available 
metadata_roadmap<-read_excel("C:/Users/ulazcano/CIC bioGUNE/Arkaitz group - Documentos/Individual folders/Uxue Lazkano/PhD/3.Projects/AC80_RRBS_Test/Annotation/Roadmap Epigenomics/Metadata/Metadata_summary.xlsx")
cell_types <- metadata_roadmap$`Epigenome ID (EID)`
## Dowload systematically all the data from web
# Base URL for ChromHMM data
base_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"

# Function to safely import BED files from the url, handling both warnings and 
# errors because in the case of the 18 state model some of the codes don't have
# all the epigenetic mark info.
safe_import <- function(file_url) {
  result <- NULL
  tryCatch(
    {
      # Capture warnings as errors
      withCallingHandlers(
        result <- import(BEDFile(file_url)),
        warning = function(w) {
          message("Warning while importing file:", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      message("Error importing file:", conditionMessage(e))
      result <- NULL
    }
  )
  return(result)
}

#######################################################################
# This is a modification of the for loop so we have for each CpG (rows) and 
# and each tissue (columns) which is the chromatin state of a given CpG

#cell_types <- c("E017" ,"E002" ,"E008", "E001")
for (cell_type in cell_types) {
  message(paste("Processing:", cell_type))
  
  # Construct the URL for the specific cell type
  file_url <- paste0(base_url, cell_type, "_18_core_K27ac_hg38lift_mnemonics.bed.gz")
  
  # Safely attempt to import the BED file
  mnemonics <- safe_import(file_url)
  
  if (!is.null(mnemonics)) {
    message(paste("Annotating overlaps for:", cell_type))
    
    # Find overlaps
    overlaps <- findOverlaps(EPICV2_gr, mnemonics)
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Initialize or update the cell_type-specific column
    if (!cell_type %in% colnames(mcols(EPICV2_gr))) {
      mcols(EPICV2_gr)[[cell_type]] <- rep(NA, length(EPICV2_gr))
      mcols(EPICV2_gr)[[paste0(cell_type, "_cat")]] <- rep(NA, length(EPICV2_gr))
    }
    
    # Loop through each category in state_classification$name
    for (category in unique(state_clasification$name)) {
      message(paste("Processing category:", category))
      
      # Filter for the current category
      category_rows <- state_clasification[state_clasification$name == category, ]
      category_indices <- which(mnemonics$name %in% category_rows$name)
      
      # Identify overlaps corresponding to the current category
      category_overlap_indices <- which(subject_hits %in% category_indices)
      category_query_hits <- query_hits[category_overlap_indices]
      
      # Assign the category name to the cell_type-specific column for the matching positions
      mcols(EPICV2_gr)[[cell_type]][category_query_hits] <- category
      category_group <- state_clasification[state_clasification$name == category, ]$Functional
      mcols(EPICV2_gr)[[paste0(cell_type, "_cat")]][category_query_hits] <- category_group
    }
  } else {
    message(paste("Skipped cell type:", cell_type))
  }
}


#Sum total counts that need to be the same in all
# Extract metadata columns as a data frame first
metadata_df <- as.data.frame(mcols(EPICV2_gr))

# Sum the last 7 columns and create a new column 'Total_counts'
EPICV2_gr$Total_counts <- rowSums(metadata_df[, c("Promoter", "Transcription", "Enhancer", 
                                                                 "Repetitive.Het", "Bivalent", "Polycomb.Repressed", "Quiescent")], 
                                                 na.rm = TRUE)
not98_annotated_AC80_haeIII_gr <- EPICV2_gr[EPICV2_gr$Total_counts != 98]


saveRDS(EPICV2_gr, "W:/ulazcano/AC80/07_ChromHMM/ChomHMM_annotated_EPIC_937055_20250124.rds")
annotated_AC80_EPICV2_gr <- readRDS("W:/ulazcano/AC80/07_ChromHMM/ChomHMM_annotated_AC80_EPICV2_gr.rds")

# List of columns to evaluate
columns_to_check <- c("Promoter", "Transcription", "Enhancer", "Repetitive-Het", 
                      "Bivalent", "Polycomb-Repressed", "Quiescent")

# Compute the major category
annotated_AC80_EPICV2_gr$major_cat <- apply(
  as.data.frame(annotated_AC80_EPICV2_gr@elementMetadata[, columns_to_check]), 
  1, 
  function(x) columns_to_check[which.max(x)]
)

# Check the result
head(annotated_AC80_EPICV2_gr)

summary_table <- table(annotated_AC80_EPICV2_gr$major_cat)

# Convert summary_table to data frame
EPICV2_summary_df <- as.data.frame(summary_table)
colnames(EPICV2_summary_df) <- c("Category", "Count")

# Create the bar plot using ggplot2
ggplot(EPICV2_summary_df, aes(x = Category, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Major Categories", x = "Category", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

