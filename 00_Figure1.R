library(tidyverse)
library(dplyr)

#########################################################################################
# Figure 1B. lollipop plot to check nº of fragments or CpGs recovered with each technique
#########################################################################################
#Lollipop plot by total reads

data_raw <- tribble(
  ~Enzyme_Type,    ~Treatment,  ~Sample_ID, ~Total_Reads,
  "MspI",          "Bisulfite", "S01",      53883217,
  "MspI",          "Bisulfite", "S04",      68340105,
  "MspI",          "Bisulfite", "S07",      83244373,
  "MspI",          "Enzymatic", "S10",      103565480,
  "MspI",          "Enzymatic", "S13",      87664425,
  "MspI",          "Enzymatic", "S16",      89082837,
  "MspI + TaqαI",  "Bisulfite", "S02",      70149496,
  "MspI + TaqαI",  "Bisulfite", "S05",      77184022,
  "MspI + TaqαI",  "Bisulfite", "S08",      77425132,
  "MspI + TaqαI",  "Enzymatic", "S11",      106995493,
  "MspI + TaqαI",  "Enzymatic", "S14",      103848112,
  "MspI + TaqαI",  "Enzymatic", "S17",      91887167, 
  "MspI + HaeIII", "Bisulfite", "S03",      80432295,
  "MspI + HaeIII", "Bisulfite", "S06",      79690265,
  "MspI + HaeIII", "Bisulfite", "S09",      87022675,
  "MspI + HaeIII", "Enzymatic", "S12",      101439442,
  "MspI + HaeIII", "Enzymatic", "S15",      100495232,
  "MspI + HaeIII", "Enzymatic", "S18",      81900855
)

# 2. Process the data
# Convert Total_Reads to millions and set Enzyme_Type as an ordered factor
data_processed <- data_raw %>%
  mutate(Reads_in_Millions = Total_Reads / 1000000) %>%
  mutate(Enzyme_Type = factor(Enzyme_Type, levels = c("MspI + HaeIII", "MspI + TaqαI", "MspI"))) %>%
  # Create a unique identifier for each combination to potentially stagger points slightly if needed
  # For a basic lollipop, we'll map directly, but this can be useful for avoiding overplotting
  mutate(Plot_Y_Axis = paste(Enzyme_Type, Treatment, sep = " - "))

# 3. Create the Dot Plot (Lollipop without the lines)
dot_plot_reads <- ggplot(data_processed, aes(x = Reads_in_Millions, y = Enzyme_Type, color = Treatment)) +
  # geom_segment() is removed as requested
  geom_point(size = 4, alpha = 0.9) + # Only the points remain
  scale_color_manual(values = c("Bisulfite" = "goldenrod", "Enzymatic" = "darkorchid4")) + # Updated colors
  labs(
    title = "Total Reads per Enzyme Digestion and Treatment Type",
    x = "Total Reads (Millions)",
    y = NULL,
    color = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(), # Still remove horizontal grid lines for clarity
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "gray80"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# 4. Save the plot as a PDF file with specified dimensions
ggsave("W:/ulazcano/AC80_RRBS/Fig1/total_reads_dot_plot.pdf",dot_plot_reads,
  width = 8,height = 4)


# Lollipop plot of total CpG mincov=1

# 1. Create the new data frame from your second table (only 'Total CpGs' columns)
data_new <- tribble(
  ~Enzyme_Type,    ~Treatment,  ~Sample_ID, ~Total_CpGs,
  "MspI",          "Bisulfite", "S01",      11777671,
  "MspI",          "Bisulfite", "S04",      13554362,
  "MspI",          "Bisulfite", "S07",      13749902,
  "MspI",          "Enzymatic", "S10",      13832279,
  "MspI",          "Enzymatic", "S13",      13428418,
  "MspI",          "Enzymatic", "S16",      12662225,
  "MspI + TaqαI",  "Bisulfite", "S02",      13856783,
  "MspI + TaqαI",  "Bisulfite", "S05",      14509435,
  "MspI + TaqαI",  "Bisulfite", "S08",      14326337,
  "MspI + TaqαI",  "Enzymatic", "S11",      16788610,
  "MspI + TaqαI",  "Enzymatic", "S14",      16852355,
  "MspI + TaqαI",  "Enzymatic", "S17",      15301635,
  "MspI + HaeIII", "Bisulfite", "S03",      15464879,
  "MspI + HaeIII", "Bisulfite", "S06",      15505657,
  "MspI + HaeIII", "Bisulfite", "S09",      15789285,
  "MspI + HaeIII", "Enzymatic", "S12",      16800483,
  "MspI + HaeIII", "Enzymatic", "S15",      16569745,
  "MspI + HaeIII", "Enzymatic", "S18",      15131428
) 


# 2. Process the data (Total_CpGs to millions, and factor order)
data_processed_new <- data_new %>%
  mutate(CpGs_in_Millions = Total_CpGs / 1000000) %>%
  mutate(Enzyme_Type = factor(Enzyme_Type, levels = c("MspI + HaeIII", "MspI + TaqαI", "MspI")))


# 3. Create the Dot Plot with the new data
dot_plot_cpgs <- ggplot(data_processed_new, aes(x = CpGs_in_Millions, y = Enzyme_Type, color = Treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Bisulfite" = "darkgoldenrod", "Enzymatic" = "darkorchid4")) + 
  labs(
    title = "Total CpGs per Enzyme Digestion and Treatment Type", 
    x = "Total CpGs (Millions)", 
    y = NULL,
    color = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "gray80"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# Print the plot (optional)
print(dot_plot_cpgs)


# 4. Save the plot as a PDF file with specified dimensions
ggsave("W:/ulazcano/AC80_RRBS/Fig1/total_cpgs_dot_plot.pdf",dot_plot_cpgs,
       width = 8,height = 4, unit="in")

################################################
# Figure 1C. Upset plot to check for common CpGs
################################################

#R in lamarr
R
.libPaths("/vols/GPArkaitz_bigdata/DATA_shared/Rocky_R/epigenomics_Rocky")

# Load mincov=1 data to do the upset plot
load(file = "/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/06_Annotation/CpGs_cov1.RData")
load(file = "W:/ulazcano/AC80_RRBS/06_Annotation/CpGs_cov1.RData")
rm(list = ls(pattern = "^S"))
save.image("W:/ulazcano/AC80_RRBS/06_Annotation/6conditions_CpGs_cov1.RData")
library("UpSetR")
library("GenomicRanges")

# Now, extract the unique 'index' values for each condition and put them in a list

list_of_cpg_sets <- list(
  MspI_Bisulfite = unique(paste0(
    as.character(seqnames(mspI_BS_CpGs_gr)),
    ":",
    start(mspI_BS_CpGs_gr),
    ":",
    as.character(strand(mspI_BS_CpGs_gr))
  )),
  MspI_Enzymatic = unique(paste0(
    as.character(seqnames(mspI_enzym_CpGs_gr)),
    ":",
    start(mspI_enzym_CpGs_gr),
    ":",
    as.character(strand(mspI_enzym_CpGs_gr))
  )),
  MspI_TaqI_Bisulfite = unique(paste0(
    as.character(seqnames(mspItaqaI_BS_CpGs_gr)),
    ":",
    start(mspItaqaI_BS_CpGs_gr),
    ":",
    as.character(strand(mspItaqaI_BS_CpGs_gr))
  )),
  MspI_TaqI_Enzymatic = unique(paste0(
    as.character(seqnames(mspItaqaI_enzym_CpGs_gr)),
    ":",
    start(mspItaqaI_enzym_CpGs_gr),
    ":",
    as.character(strand(mspItaqaI_enzym_CpGs_gr))
  )),
  MspI_HaeIII_Bisulfite = unique(paste0(
    as.character(seqnames(mspIhaeIII_BS_CpGs_gr)),
    ":",
    start(mspIhaeIII_BS_CpGs_gr),
    ":",
    as.character(strand(mspIhaeIII_BS_CpGs_gr))
  )),
  MspI_HaeIII_Enzymatic = unique(paste0(
    as.character(seqnames(mspIhaeIII_enzym_CpGs_gr)),
    ":",
    start(mspIhaeIII_enzym_CpGs_gr),
    ":",
    as.character(strand(mspIhaeIII_enzym_CpGs_gr))
  ))
)
# --- Step 2: Generate the UpSet plot ---
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/Fig1/1C_upset_plot_common_cpgs_V4.pdf", width = 4, height = 4)
pdf("W:/ulazcano/AC80_RRBS/Fig1/1C_upset_plot_common_cpgs_V7.pdf", width = 4, height = 4.5)
upset(fromList(list_of_cpg_sets), nsets = 6, nintersects = 20, order.by = "freq", decreasing = T,
      mainbar.y.label = "Number of Common CpGs in Intersection",
      sets.x.label = "Total Number of Common CpGs Per Condition",
      text.scale = c(0.8, 0.7, 0.5, 0.6, 1, 0.3),
      point.size = 2, line.size = 0.5)
dev.off()

################################################
# Figure 1D. CpGs vs depth
################################################
all_granges_list <- list(
  MspI_Bisulfite = mspI_BS_CpGs_gr,
  MspI_Enzymatic = mspI_enzym_CpGs_gr,
  MspI_TaqI_Bisulfite = mspItaqaI_BS_CpGs_gr,
  MspI_TaqI_Enzymatic = mspItaqaI_enzym_CpGs_gr,
  MspI_HaeIII_Bisulfite = mspIhaeIII_BS_CpGs_gr,
  MspI_HaeIII_Enzymatic = mspIhaeIII_enzym_CpGs_gr
)

# --- Step 1: Extract and prepare coverage data from all GRanges objects ---
# Function to extract coverage data from a single GRanges object
extract_and_pivot_coverage <- function(gr_object, condition_name) {
  # Get all metadata columns
  mcols_df <- as.data.frame(mcols(gr_object))
  
  # Filter for columns that start with "coverage"
  coverage_cols_df <- mcols_df %>%
    select(starts_with("coverage"))
  
  if (ncol(coverage_cols_df) == 0) {
    warning(paste("No 'coverage' columns found in GRanges object for condition:", condition_name))
    return(NULL) # Return NULL if no coverage data
  }
  
  # Pivot to long format: each row is a CpG from a specific replicate
  df_long <- coverage_cols_df %>%
    pivot_longer(
      cols = everything(), # Select all coverage columns
      names_to = "Replicate_ID", # New column for replicate name (e.g., "coverage1")
      values_to = "Coverage_Value" # New column for the coverage value
    ) %>%
    mutate(
      Condition = condition_name, # Add condition name
      # Optionally, extract just the number for replicate ID if needed for plotting replicates separately
      Replicate_Number = as.integer(str_extract(Replicate_ID, "\\d+"))
    )
  return(df_long)
}

# Apply the function to all GRanges objects in the list and combine
all_coverage_data <- map_dfr(names(all_granges_list), ~ {
  extract_and_pivot_coverage(all_granges_list[[.x]], .x)
})

# Filter out rows with NA coverage values if any, and ensure Coverage_Value is numeric
all_coverage_data <- all_coverage_data %>%
  filter(!is.na(Coverage_Value)) %>%
  mutate(Coverage_Value = as.numeric(Coverage_Value))

# --- Step 2: Plot the distribution of CpGs with coverage ---
# We'll use geom_freqpoly to show the frequency distribution as lines.
# Each line will represent a condition, combining its replicates.
pdf("/vols/GPArkaitz_bigdata/ulazcano/AC80_RRBS/Fig1/1D_cpg_coverage_distribution.pdf", width = 8, height = 4)
pdf("W:/ulazcano/AC80_RRBS/Fig1/1D_cpg_coverage_distribution_V6.pdf", width = 4, height = 4)
ggplot(all_coverage_data, aes(x = Coverage_Value, color = Condition)) +
   geom_freqpoly(binwidth = 1, linewidth = 0.8, alpha = 0.7) +
   labs(
     title = "Distribution of CpG Coverage Depth per Condition",
     x = "Coverage Depth",
     y = "Number of CpGs (Frequency)",
     color = "Condition"
   ) +
   scale_color_brewer(palette = "Set1") +
   theme_classic(base_size = 8) +
   theme(
     plot.title = element_text(hjust = 0, face = "bold"),
     legend.position = "right",
     legend.title = element_text(size = 4), # Smaller size for the legend title
     legend.text = element_text(size = 4),
     legend.spacing.y = unit(-0.8, "cm")
   ) +
   coord_cartesian(xlim = c(0, 100)) # Adjust max_coverage as needed
 dev.off()
 
 ########Same plot but each tratment in one line
 
 # Define the specific conditions you want to combine for each line
 # Replace these with the actual names from your 'Condition' column
 bisulfite_conditions <- c("MspI_Bisulfite", "MspI_TaqI_Bisulfite", "MspI_HaeIII_Enzymatic")
  enzymatic_conditions <- c("MspI_Enzymatic", "MspI_TaqI_Enzymatic", "MspI_HaeIII_Enzymatic")
 
 # Filter the data and create a new grouping variable
 filtered_coverage_data <- all_coverage_data %>%
   filter(Condition %in% enzymatic_conditions | Condition %in% bisulfite_conditions) %>%
   mutate(
     Line_Type = case_when(
       Condition %in% enzymatic_conditions ~ "Enzymatic",
       Condition %in% bisulfite_conditions ~ "Bisulfite",
       TRUE ~ "Other" # This ensures any non-matching conditions are handled, though filtered out
     )
   )
 
 pdf("W:/ulazcano/AC80_RRBS/Fig1/1D_cpg_coverage_distribution_V7.pdf", width = 4, height = 4)
 ggplot(filtered_coverage_data, aes(x = Coverage_Value, color = Line_Type)) +
   geom_freqpoly(binwidth = 1, linewidth = 0.8, alpha = 0.7) +
   scale_color_manual(values = c("Bisulfite" = "darkgoldenrod", "Enzymatic" = "darkorchid4")) + 
   labs(
     title = "Distribution of CpG Coverage Depth per Condition",
     x = "Coverage Depth",
     y = "Number of CpGs (Frequency)",
     color = "Condition Type"
   ) +
   #scale_color_brewer(palette = "Set1") +
   theme_classic(base_size = 8) +
   theme(
     plot.title = element_text(hjust = 0, face = "bold"),
     legend.position = "right",
     legend.title = element_text(size = 4),
     legend.text = element_text(size = 4),
     legend.spacing.y = unit(-0.8, "cm")
   ) +
   coord_cartesian(xlim = c(0, 100))
 dev.off()

 ########Same plot but each tratment in one line
 
 # Define your desired colors for each condition
 # You'll need to list all your unique 'Condition' names and assign a color to each.
 # For demonstration, let's use some example condition names.
 # Replace these with the actual unique values from your 'Condition' column.
 custom_colors <- c(
   "MspI_TaqI_Bisulfite" = "#FFFF00",  # Bright Yellow
   "MspI_HaeIII_Enzymatic" = "#FFD700",  # Gold Yellow
   "MspI_TaqI_Enzymatic" = "#800080",  # Purple
   "MspI_Enzymatic" = "#BA55D3",  # Medium Orchid (a lighter purple)
   "MspI_Bisulfite" = "#FFEC8B", # Light Yellow
   "MspI_HaeIII_Enzymatic" = "#9400D3" # Dark Violet (a darker purple)
 )
 
 
 # Assuming all_coverage_data is already loaded
 
 # Define your desired color palette.
 # We'll create a function to generate a set of yellow and purple colors
 # based on how many unique conditions you have.
 # You might need to adjust the numbers of colors if you have many conditions.
 generate_yellow_purple_colors <- function(n_colors) {
   # These are good starting points for yellow and purple hues
   yellow_hues <- seq(from = 40, to = 60, length.out = ceiling(n_colors / 2)) # Yellow range
   purple_hues <- seq(from = 270, to = 300, length.out = floor(n_colors / 2)) # Purple range
   
   # Combine and convert to hex codes. You can play with `c` (chroma) and `l` (luminance)
   # to get different shades and saturation.
   colors <- c(
     hcl(h = yellow_hues, c = 70, l = 75),  # Lighter, more vibrant yellows
     hcl(h = purple_hues, c = 70, l = 50)   # Deeper, richer purples
   )
   
   # If n_colors is odd, we might have one extra yellow or purple.
   # Just take exactly n_colors.
   head(colors, n_colors)
 }
 
 # Get the unique conditions from your data
 unique_conditions <- unique(all_coverage_data$Condition)
 num_unique_conditions <- length(unique_conditions)
 
 # Generate colors for each unique condition
 # You can manually assign if you prefer specific colors for specific conditions
 # otherwise, this will automatically assign from the generated palette.
 # If you have specific conditions you want to always be yellow and others always purple,
 # a `case_when` for `scale_color_manual` is more robust, but this is for simplicity.
 condition_colors <- setNames(generate_yellow_purple_colors(num_unique_conditions), unique_conditions)
 
 pdf("W:/ulazcano/AC80_RRBS/Fig1/1D_cpg_coverage_distribution_V8.pdf", width = 4, height = 4)
 ggplot(all_coverage_data, aes(x = Coverage_Value, color = Condition)) +
   geom_freqpoly(binwidth = 1, linewidth = 0.8, alpha = 0.7) +
   labs(
     title = "Distribution of CpG Coverage Depth per Condition",
     x = "Coverage Depth",
     y = "Number of CpGs (Frequency)",
     color = "Condition" # Legend will show each unique condition
   ) +
   # Apply the generated colors to your conditions
   scale_color_manual(values = condition_colors) +
   theme_classic(base_size = 8) +
   theme(
     plot.title = element_text(hjust = 0, face = "bold"),
     legend.position = "right",
     legend.title = element_text(size = 4),
     legend.text = element_text(size = 4),
     legend.spacing.y = unit(-0.8, "cm")
   ) +
   coord_cartesian(xlim = c(0, 100))
 dev.off()
 
 ################################################
 # Figure 1E. Insert sie distribution
 ################################################
 #First need to extract insert size for each samples
 #This is the code in bash
 #source /opt/ohpc/pub/apps/samtools/samtools-1.19.2/cic-env
 #samtools view S01_R1_val_1_bismark_bt2_pe.bam | awk '{print $9}' > insert_sizes_S01.txt

 #Calculate CG content for each read, do not need to do the previous step because in this 
 #fragment size is together wit CG content so I can do all the plots from the same file
 samtools view S04_R1_val_1_bismark_bt2_pe.bam |
 awk '{read_id=$1; frag_len=$9; seq=$10; cg=gsub(/CG/, "", seq); print read_id"\t"frag_len"\t"cg}' > S04_CG_content.txt
 
 awk '{read_id=$1; seq=$10; len=length(seq); cg=gsub(/CG/, "", seq); print read_id"\t"len"\t"cg}' > S01_CG_content.txt
 
 
 S01 <- as.data.table(fread("W:/ulazcano/AC80_RRBS/Fig1/S01_CG_content.txt"))
 S04 <- as.data.table(fread("W:/ulazcano/AC80_RRBS/Fig1/S04_CG_content.txt"))
 
 