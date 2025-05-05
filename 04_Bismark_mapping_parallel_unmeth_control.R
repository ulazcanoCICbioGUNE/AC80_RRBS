
#################################################################################
#                            MAPPING WITH Bismark 
#################################################################################

# Summary
#---------



#################################################################################
#                           Bismark VERSION
#################################################################################
# Bismark Version: v0.24.2
# Copyright 2010-23 Felix Krueger, Altos Bioinformatics
# https://github.com/FelixKrueger/Bismark
#
#
################################################################################
#                         PARAMETERS SELECTION
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
#path <- "W:/ulazcano/"

# Log file date from 0_Sample_info_XXX.log
log_folder <- "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/log/"
logdate <- "20250401"

### Process information ###
partition <- "FAST"
time <- c("00:40:00")
memory <- c("2")
node <- 1
cpu <- 18
ram <- as.numeric(memory)*10^9

#-------------------------------------------------------------------------------------------------------------------------------------------------------
dir_infiles <- paste(path, project, "/02_TRIMMED", sep = "" )

# Output directory
dir_out <- paste(path, project, "/04_Bismark", sep = "")
# Create output directory
dir.create(file.path(dir_out,"/Unmeth_control"))
dir_outfiles <- paste(dir_out,"/Unmeth_control",sep='')
# Set directory
setwd(dir_outfiles)

# Genome index
indexfolder <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Ref_genome/RREMseq_controls/Unmeth/"
### FASTQs ###
# List  fastq files
pattern = "_R1_val_1.fq.gz"
pattern2 = "_R2_val_2.fq.gz"

### FASTQs ###
# List  fastq.gz files
samples <- list.files(path=dir_infiles, pattern = pattern)
# Filter sample name 
samples_names <- gsub(pattern,"", samples)


###############################################################
#                           Bismark
###############################################################

### Generate PBS ###
for (i in 1:length(samples_names)) {
  # Job name
  job_name <- paste(samples_names[i],"_Bismark_alignment_unmeth_control",sep='');
  
  
  # Input files
  input1 <- paste(samples_names[i], pattern, sep="")
  input2 <- paste(samples_names[i], pattern2, sep="")
  
  # Output files in this directory
  output <- paste(dir_outfiles, "/", job_name, "/", sep="")
  
  # Command
  command <- paste("bismark --genome /vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Ref_genome/RREMseq_controls/Unmeth/ -1 ", path, project, "/02_TRIMMED/" ,input1,
                   " -2 ",path, project, "/02_TRIMMED/", input2, " --multicore 18", sep='')
  
  
  # SBATCH File
  filename <- paste(job_name,".sh",sep='');
  cat(
    c("#!/bin/sh"),
    c("#SBATCH  --export=ALL"),
    paste("#SBATCH --job-name=",job_name,sep=''),
    paste("#SBATCH --partition=",partition,sep=''),
    c("#SBATCH --ntasks=1"),
    c("#SBATCH --exclude=gn[00,01,02]"),
        paste("#SBATCH --nodes=", node,sep =''), 
    paste("#SBATCH --cpus-per-task=",cpu, sep=''),
    paste("#SBATCH --time=",time,sep=''),
    paste("#SBATCH --mem=",memory,"GB",sep=''),
    paste("#SBATCH -o ",job_name,".out",sep=''),
    paste("#SBATCH -e ",job_name,".err",sep=''),
    c("source /opt/ohpc/pub/apps/anaconda3/cic-env"),
    c("conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/Bismark"),
    c(paste("cd ", dir_outfiles)), 
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
logfile$Date <- Sys.time()
logfile$Trimming <- trmd

write.table(as.data.frame(logfile), paste(path, project_name, "/log/4_Bismark_align_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



q()