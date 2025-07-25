################################################################################
#                    TRIMMING FASTQ WITH TRIMGALORE!
################################################################################
# Author: Uxue Lazcano
# Date of creation: 2024/11/21

# Summary
# ---------
#Trim Galore is a tool designed to preprocess sequencing data by trimming adapter sequences and low-quality bases.
# In this project --rrbs mode is going to be used so and additional step is going to be performed.
# Key functions of TrimGalore!
# 1. Specific trimming for reduce representation data (default MspI digestion)
# As recommended by Vazyme (respues A Aransay 2025/04/11 por Teams) For the RRBS, you need to trim 12 nucleotides 
# at the 3' end of R1 reads (post-QC) and at the 5' end of R2 reads (post-QC)
# 2. Adapter trimming
# It removes Illumina adapter sequences that may remain in the reads, ensuring these do not interfere with downstream analysis.
# 3. Trimming Low-Quality Bases
# Bases below a specified quality threshold (default is usually Phred score 20) at the ends of the reads are trimmed to 
# improve the overall data quality.
# 4. Length Filtering
# After trimming, very short reads (default < 20 bp) are discarded as they are unlikely to align uniquely to the reference genome.
# 5. Single- and Paired-End Support
# In paired-end mode, it ensures that both reads of a pair are trimmed consistently, 
#and if one read becomes too short after trimming, the pair is discarded.
# 6. Directionality of libraries
# Handling Directionality in RRBS Libraries:
# RRBS libraries are often strand-specific due to the bisulfite conversion process. Trim Galore accounts for this by recognizing 
# the strand-specific nature of RRBS data and ensuring that only reads aligned to the correct strand orientation are processed correctly, 
# which is essential for accurate methylation calling.


################################################################################
#                         TRIMGALORE VERSION
################################################################################

# Conda version is the latest in March 2022
# Conda 4.12.0
#
# Condaconda is located in
# /opt/ohpc/pub/apps/anaconda3/cic-env
# 
# Trim Galore version version 0.6.10 Last update: 02 02 2023, with Python 3.10.10 in environment iRRBS_prueba
# 
# Link: https://anaconda.org/bioconda/trim-galore

################################################################################
#                             PIPELINE
################################################################################
### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
# Load R libraries 
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")


#-------------------------------------------------------------------------------------------------------------------------------------------------------
### General Project ###

# Project Name 
project <- "AC80_RRBS"

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "/fastdata/GPArkaitz_fastdata/ulazcano/"
files_rocky <- "/vols/GPArkaitz_bigdata/DATA_shared/"

# Date of the log file 0_Sample_info_XXX.log
log_folder <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/log/"
logdate <- "20250401"


### Process information ###
partition <- "FAST"
time <- c("06:00:00")
memory <- c("10")
cpu <- 4


### Trim Galore! parameters ###

#   --rrbs | It triggers specialized trimming optimized for MspI-digested libraries, 
#           including removing known bias at fragment ends and handling overhangs.
#   --non_directional | For non directional libraries, not my case (RREMseq library)
#   --paired | Indicates that the input files are paired-end reads.
#   -o | output directory
#   -- illumina | It trims Illumina adapter sequences by default, not need to specify (as in Cutadapt)
#   -q | by default is set to 20
#   -l (--lenght) | by default is set to 20bp 

# For this analysis I leave all the default parameters as they are
# q <- 20
# l <- 20

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Output directories 
dir_out <- paste(path, project, sep = "")

# Load log file 
logfile <- read.table(paste(log_folder, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory
#dir_infiles <- paste(logfile$filedirRocky, "FASTQs", sep = "")
dir_infiles <- paste(files_rocky,"AC-80_RRBS/TEST_Vazyme_EMvsBS/FASTQs", sep="")

# Create output directory
dir.create(file.path(dir_out,"02_TRIMMED"))
dir_outfiles <- paste(dir_out,"/02_TRIMMED",sep='')
setwd(dir_outfiles)


### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path = dir_infiles, pattern = ".fastq.gz")
# Filter sample name 
samples_names <- gsub(".fastq.gz", "", samples)
samples_names <- unique(gsub("_.*", "", samples_names))

###############################################################
#                      TRIMMING
############################################################### 


### Generate PBS ###
for (i in 1:length(samples_names)) {
#for (i in 1:1) {
  print(samples_names[i])
  # Job name
  job_name <- paste(samples_names[i],"_Trimmed",sep='');

  # Input file
  input1 <- paste(dir_infiles, "/", samples_names[i],"_R1.fastq.gz", sep="")
  input2 <- paste(dir_infiles, "/", samples_names[i],"_R2.fastq.gz", sep="")
  
  # Output file: Adapter trimmed, q<10 and  m>30
  output1 <- paste(dir_outfiles, "/", samples_names[i],"R1_trmd.fastq.gz", sep="")
  output2 <- paste(dir_outfiles, "/", samples_names[i],"R2_trmd.fastq.gz", sep="")

  # Cutadapt command
  #The Illimuna universal adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCTATCAATCTCGTA
  # and all the overepresented seuqences in the 01 FATSQC
  command <- paste("trim_galore --illumina  --clip_R2 12 --three_prime_clip_R1 12 --paired", input1, input2, "-o", dir_outfiles, sep= " ")
  
  # SBATCH File
  filename <- paste(job_name,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    c("#SBATCH  --export=ALL"),
    paste("#SBATCH --job-name=",job_name,sep=''),
    paste("#SBATCH --partition=",partition,sep=''),
    paste("#SBATCH --cpus-per-task=",cpu,sep=''),
    paste("#SBATCH --time=",time,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",job_name,".out",sep=''),
    paste("#SBATCH -e ",job_name,".err",sep=''),
    c(paste("cd ", dir_infiles)),
    c("source /opt/ohpc/pub/apps/anaconda3/cic-env"),
    c("conda activate /vols/GPArkaitz_bigdata/ulazcano/envs/eRREMseq"),
    c(command),
    file=filename,sep = "\n",append=F)
  
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
  }


#######################################################################
#                            LOG FILE                        
#######################################################################


# Save log file information
logdate <- format(Sys.time(), "%Y%m%d")
log_data <- c()
log_data$Date <- Sys.time()
log_data$project <- project
log_data$outputdir <- dir_outfiles
log_data$inputdir <- dir_infiles

log_data$command <- command

write.table(as.data.frame(log_data), paste(log_folder, "2_Trimming_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")


q()


