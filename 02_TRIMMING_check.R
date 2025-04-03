################################################################################
#                             TRIMMED FASTQ SUMMARY                            
################################################################################
#Uxue Lazkano (Based in a script of Maria Ponce https://github.com/AC-lab-pipelines/DEG_Reference/blob/main/Scripts/02_trimming_result_check.R)
#2025/03/03 

# Summary
# ---------

# After running Trimgalore!, we want to evaluate the trimming and phred score base quality removal.
# We select the total number of reads processed, the reads with adapter content, 
# reads too short, quality trimmed and the total reads written.

# Read .err file of read 2 that contains the info of removed sequences and extract the desired data



################################################################################
#                         LOAD FILES AND DATA                           
################################################################################
library(dplyr)

# -----------------------------------------------------------------------------------------------------------------
### General Project ###
# Project Name 
project_name <- "AC80_RRBS"

# File path
path <- "/fastdata/GPArkaitz_fastdata/ulazcano/"
#path <- "/vols/GPArkaitz_bigdata/ulazcano/"

# -----------------------------------------------------------------------------------------------------------------

# Output and input directories SSH 
input_dir <- paste(path, project_name, "/02_TRIMMED/", sep = "" )
setwd(input_dir)

# List fastq.zip files
samples <- list.files(pattern = ".fq.gz")
# Filter sample name 
#samples_names <- gsub("_val_.*", "", samples)
#samples_names <- unique(sub("^(S_\\d+)_.*", "\\1", samples))
samples_names <- unique(gsub("_R[12].*", "", samples))

# Define patterns (escaping special characters)
pattern <- c(
  "Total number of sequences analysed:",
  "Number of sequence pairs removed because at least one read was shorter ",
  "Sequences were truncated to a varying degree because of deteriorating ",
  "RRBS reads trimmed by additional 2 bp when adapter contamination was detected:"
)



# -- Adaptador --
# Total read pairs processed
# Reads with adapters
# Reads written (passing filters)
#
# -- seq menores 30 --
# Total basepairs processed
# 
# -- seq calidad mayor 10 --
# Quality-trimmed
# Total written (filtered)


# Final matrix
def_file <- data.frame()


################################################################################
#                            PROCESS                           
################################################################################


for (i in 1:length(samples_names)){
  # Sample name
  file_name <- paste(samples_names[i], "_Trimmed.err", sep="")
  
  # Read .out file
  # out_file <- read.fwf(file_name)
  # out_data <- data.frame(read.delim(file = file_name, header = FALSE))
  #out_data <- readLines(file_name)
  out_data <- tail(readLines(file_name), 24)
  # Select presenting the information of interest
  out_data <- data.frame(x = out_data[grepl(paste(pattern, collapse = "|"), out_data)])
  
  
  # Split the strings into two columns using ": " as separator
  #   First column. Variable of interest corresponding to pattern
  #   Second column. Value of the variable
  #out_file <- data.frame(do.call("rbind",strsplit(as.character(out_data[,1]), ": ", fixed = TRUE)))
  out_file <- data.frame(do.call("rbind", lapply(as.character(out_data[,1]), function(x) {
    parts <- strsplit(x, ":\\s*(?=[^:]+$)", perl = TRUE)[[1]]
    return(parts)
  })))
  
  ## Key step for the final output
  # Transpose matrix
  # Output
  #   CPU Efficiency  Job Wall-clock  Memory Efficiency
  #      Value 1          Value 2           Value 3
  out_data1 <- setNames(data.frame(t(out_file[,-1])), out_file[,1])
  # out_data2 <- out_data1[,-c(3,7)]
  colnames(out_data1) <- c("Phred score <20", "RRBS reads trimmed", "Total reads", "Length < 20bps")
  
  # Aggregate all the information together of the corresponding file 
  out_data3 <- cbind(Sample=paste(samples_names[i]), out_data1)
  
  # Reorder the columns
  out_data3 <- out_data3[, c("Sample", "Total reads", "Length < 20bps", "RRBS reads trimmed", "Phred score <20")]
  
  # Aggregate all the information from every file
  def_file <- rbind(def_file, out_data3)
}

################################################################################
#                   SAVED DATA: Performance               
################################################################################
###ERREPASAU HAUUUUUUUUUUUUUUUU
# Function to remove content inside parentheses
remove_parentheses <- function(x) gsub(" \\(.*?\\)", "", x)

# Apply function to all columns except Sample and Total reads
def_file <- def_file %>%
  mutate(across(-c(Sample, `Total reads`), remove_parentheses))

# Convert numeric columns to proper type
def_file <- def_file %>%
  mutate(across(-c(Sample), as.numeric))
# Create 'Trimmed Reads' column
def_file$`Trimmed Reads` <- as.numeric(def_file$`Total reads`) - as.numeric(def_file$`Length < 20bps`)

# Add a new column with the percentage of reads that passed trimming
def_file <- def_file %>%
  mutate(`Trimmed Reads %` = (`Trimmed Reads` / `Total reads`) * 100)

# Print the updated dataframe
print(def_file)

write.csv(def_file,file = paste(input_dir ,"/", project, "_Trimming_check.csv", sep = ""), row.names=FALSE)