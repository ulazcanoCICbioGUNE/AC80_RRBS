
#################################################################################
#                            Methylation extraction Bismark 
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
# path <- "W:/ulazcano/"
path_bigdata <- "/vols/GPArkaitz_bigdata/ulazcano/"
# Trimmed fastq/ raw fastqs
trmd <- TRUE

# Log file date from 0_Sample_info_XXX.log
log_folder <- "/vols/GPArkaitz_bigdata/ulazcano/AC80/log/"
logdate <- "20241121"

# Load log file 
logfile <- read.table(paste(log_folder, "/0_Sample_info_", logdate, ".log", sep = ""), header = TRUE)

### Process information ###
partition <- "FAST"
time <- c("08:00:00")
memory <- c("15")
node <- 1
cpu <- 18
ram <- as.numeric(memory)*10^9

#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Specie
specie <- "Human"

#Input file directory, mapped files
dir_infiles <- paste(path_bigdata, project, "/04_Bismark/Meth_control/", sep="")

# Create output directory
dir.create(file.path(path_bigdata, project,"/05_Methylation_extraction/Meth_control/"))
dir_outfiles <- paste(path_bigdata, project,"/05_Methylation_extraction/Meth_control",sep='')
#LUEGO BORRAR
#dir_outfiles <- "/vols/GPArkaitz_bigdata/ulazcano/"
# Set directory
setwd(dir_outfiles)


# Genome index
indexfolder <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Ref_genome/RREMseq_controls/Meth/"

### FASTQs ###
# List  fastq.gz files
pattern = "_R1_val_1_bismark_bt2_pe.bam"
samples <- list.files(path=dir_infiles, pattern = pattern)
# Filter sample name 
samples_names <- gsub(pattern,"", samples)

###############################################################
#          Bismark Methylation Extractor
###############################################################

### Generate PBS ###

for (i in 1:length(samples_names)) {
#for (i in 1:1) {
  # Job name
  job_name <- paste("Meth_",samples_names[i],"_BISMARK_methylation_extractor",sep='');
  
  # Input files
  input1 <- paste(samples_names[i], pattern, sep="")
  
  # Output files in this directory
  output <- paste(dir_outfiles, "/", job_name, "/", sep="")

  # Command
  command <- paste("bismark_methylation_extractor -p --no_overlap --report --comprehensive --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder ",indexfolder,  " -o ", dir_outfiles, " ",dir_infiles, samples_names[i],"_R1_val_1_bismark_bt2_pe.bam  --multicore 18", sep="")

  
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
    paste("#SBATCH --cpus-per-task=",cpu,
          sep=''),
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
logfile$command_methextract <- command

write.table(as.data.frame(logfile), paste(log_folder, "/5_Bismark_methExtractor_", logdate, ".log", sep = ""), row.names = FALSE, eol = "\r")



q()

