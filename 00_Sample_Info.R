
################################################################################
#                               METADATA                          
################################################################################


# Summary
# --------
#
# Generate a sample information file also refereed as metadata. The sample 
# information file is used since the beginning of the analysis to generate 
# different scripts and in the RRBS analysis which is based 
# in the project information. Allowing the execution of the scripts with minor 
# adjustments.
# 
# This file must be MANUALLY adjusted for the different projects. 

#######CAMBIAR SEGÚN YO LO HAGA, HABRÁ UNA VERSIÓN PARA SRA Y OTRA PARA CUANDO ME MANDE ANA LOS DATOS ######
# Content
#   - Named Sample_info.csv
#   - Must have at least, the following variables:
#     - Sample: Variable with the samples names. MUST NAME Sample
#     - Treatment: Specific treatment per sample. Must be a factor variable, can 
#       be named as you prefer but must be consistent for the following analysis


# Folder
# Input: W:/DATA_shared/Sequencing_data_folder{ACXX}/
# Output: W:/PersonalFolder/Project folder


################################################################################
#                                 LOAD DATA                         
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "AC80"


# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
files_rocky <- "/vols/GPArkaitz_bigdata/ulazcano/AC80/"
path <- "W:/ulazcano/AC80/"

# Local directory with Git folders 
local_dir <- "C:/Users/ulazcano/CIC bioGUNE/Arkaitz group - Documentos/Individual folders/Uxue Lazkano/PhD/3.Projects/AC80_RRBS_Test/"


# Input directory
# Must be the folder where the library preparation pdf is found and folder 
# with the fastq are included 
files <- "W:/DATA_shared/AC-80_RRBS_Test/"
#files <- paste(path,project,"_fastqs","/FASTQ/",sep="" )
# Select input directory
dir_in <- files

# First and last sample found in the Library Preparation pdf table which is 
# usually in page 5, corresponding to the firs column which can be Library ID or GAP ID
# This is a key step to generate a data frame based on the table from the pdf
first_sample <- "AC-80_23"
last_sample <- "AC-80_32"

# Condition
trt <- "Enzyme"

# Contrast level order 
lvl_order <- c("TaqaI", "HaeIII")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


################################################################################
#                         PROJECT CONFIGURATION
################################################################################
## Create folders BigData folders
# Create project folder in BigData
#dir.create(file.path(path, project))
#dir_out <- paste(path, project, sep = "")
dir_out <- "W:/ulazcano/AC80-prueba/"
setwd(dir_out)

# Create log file folder
dir.create(file.path(dir_out, "log"))
dir_log <- paste(dir_out, "/log", sep = "")


################################################################################
#                     LOAD LIBRARIES AND FUNCTIONS                           
################################################################################
# Load libraries
source(paste(local_dir, "Scripts", "/libraries_RREMseq.R", sep = ""))
# Load functions 
source(paste(local_dir, "Scripts", "/function_pdf_to_tab.R", sep = ""))

################################################################################
#                         SAMPLE INFORMATION FILE
################################################################################
#In this case we do not have the proper sample infor file from GAP, need to revise with them the delivery for subsequent analysis
#I'm going to manually add the information of the samples

# List of samples names
samples <- unique(gsub("_1.fastq.gz","", list.files(path = paste(dir_in, "FASTQs/", sep=""), pattern = "_1.fastq.gz")))


# PDF file with the library preparation report
# This file should be in the project DATA_shared folder 
file <- pdf_text(paste(dir_in, list.files(path = dir_in, pattern = "Library_Preparation"), sep = ""))


################################################################################
#                             SAMPLE INFORMATION                           
################################################################################
# In the case of AC80 sample names in "Library prep" and  in fastq files does not match
#This need to be corrected in following project but for the moment I manually crate the table
# Run the function
data <- pdf_to_tab(file, first_sample, last_sample)
#print(head(data))
# Create the enzyme column
enzyme <- c(rep("TaqaI", 3), rep("HaeIII", 3))

# Combine into a data frame
data <- data.frame(Sample = samples, Enzyme = enzyme)

# Number of samples
n <- length(samples)

# Verification
ifelse(nrow(data) == n, print("The number of samples from the data table in the pdf match the samples in the folder"), paste("ERROR: Samples do not match", "CHANGE THE NUMBER OF ROWS REMOVED IN THE FUNCTION", sep = "\n."))


# Save data as sample information 
sample_info <- data
#######################################################################
#                             SAVE DATA                         
#######################################################################

# Save sample information file, it is saved in dir_out
write.csv(sample_info, file = paste("Sample_info.csv", sep = ""), row.names = FALSE)

#######################################################################
#                            LOG FILE                        
#######################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project_name <- project
log_data$condition <- trt
log_data$condition_order <- paste0(lvl_order, collapse =",")
log_data$path <- dir_out
log_data$filedir <- dir_in
log_data$filedirRocky <- files_rocky

write.table(as.data.frame(log_data), paste(dir_log, "/0_Sample_info_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



