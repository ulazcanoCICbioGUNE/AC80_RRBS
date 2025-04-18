################################################################################
#                           FASTQC RAW FASTQ              
################################################################################

# Summary
# ---------

# Quality control of the raw fastq files

# Folder
# Input: W:/DATA_shared/Sequencing_name/
# Output: Project_folder/01_FASTQC


################################################################################
#                             FASTQC VERSION
################################################################################

# FastQC version is the latest in 2023: 
# FastQC v0.12.1
#
# FastQC is located in: 
# /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/fastqc 

# Link: https://github.com/s-andrews/FastQC

################################################################################
#                             PIPELINE
################################################################################
#

### Access R ###
source /opt/ohpc/pub/apps/R/R-4.2.1/cic-R
R
# Load R libraries 
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")

#-------------------------------------------------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project_name <- "AC80_RRBS"


# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "/vols/GPArkaitz_bigdata/ulazcano/"
files_rocky <- "/vols/GPArkaitz_bigdata/DATA_shared/"

# Date of the log file 0_Sample_info_XXX.txt
logdate <- "20250401"

### Process information ###
cluster <- "FAST"
walltime <- c("00:25:00") 
cpu <- 1
memory <- c("3")
#-------------------------------------------------------------------------------------------------------------------------------------------------------


# Load log file 
logfile <- read.table(paste(path, project_name, "/log/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

# Input directory 
dir_infiles <- paste(files_rocky,"AC-80_RRBS/TEST_Vazyme_EMvsBS/FASTQs", sep="")

# Create output directory
dir_outfiles <- paste(path, project_name, sep = "")
dir.create(file.path(dir_outfiles, "01_FASTQC"))
dir_outfiles <- paste(dir_outfiles, "/01_FASTQC",sep='')
setwd(dir_outfiles)

### FASTQs ###
# Instead of using the project name, here we use 
# the files names
# List fastq.gz files
samples <- list.files(path = dir_infiles, pattern = ".fastq.gz")
# Filter sample name 
samples_names <- gsub(".fastq.gz", "", samples)


### Generate PBS ###
for (i in 1:length(samples)) {
  # Job name
  name.job <- paste(samples_names[i],"_FastQC",sep='');
  # Command
  command <- paste("/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/FastQC/fastqc", samples[i], "-o ", dir_outfiles)
  # SBATCH File
  filename <- paste(name.job,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    c("#SBATCH  --export=ALL"),
    paste("#SBATCH --job-name=",name.job,sep=''),
    paste("#SBATCH --partition=",cluster,sep=''),
    paste("#SBATCH --cpus-per-task=",cpu,sep=''),
    paste("#SBATCH --time=",walltime,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",name.job,".out",sep=''),
    paste("#SBATCH -e ",name.job,".err",sep=''),
    c(paste("cd ", dir_infiles)),
    c(command),
    file=filename,sep = "\n",append=F)
  # Run PBS
  system(paste("sbatch",filename,sep=' '));
  Sys.sleep(3)
}

q()

