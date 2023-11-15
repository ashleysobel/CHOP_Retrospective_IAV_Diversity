# This script imports and organizes the sequences from the Mid-Atlantic States for additional analysis of the outlier "CHOA-101"


# Clear environment to get rid of old files and get path -------------------
rm(list = ls())
R_path <- getwd()

# Create file path for permutation analysis
Outlier_analysis_path <- file.path(R_path,"Outlier_analysis")

# Set segment names
segment_name <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")

# Load packages
library(tidyverse) #use for str_replace_all, str_replace
library(phylotools) # use for read.fasta
library(dplyr) # use for select
library(stringi) # use for stri_extract_last
library(ggplot2) # Use for generating the coverage plots
library(seqinr) # used for write.fasta, translate
library(tidysq)
library(patchwork)
library(purrr)
library(Biostrings)
library(DECIPHER)
library(boot)

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("tidyverse", "phylotools", "dplyr","stringi","seqinr","ggplot2","tidysq","patchwork","purrr","Biostrings","DECIPHER","boot")

# Initialize an empty character vector to hold citations
all_citations <- character(0)

# Loop through each package and get its citation in BibTeX format
for (pkg in package_list) {
  citation_info <- utils::toBibtex(citation(pkg))
  all_citations <- c(all_citations, citation_info)
}

# Combine all citations into a single string
all_citations_text <- paste(all_citations, collapse = "\n\n")

# Write the combined citations to a .bib file
writeLines(all_citations_text, con = file.path(Outlier_analysis_path,"outlier_citations.bib"))


# Load the fasta files 
st <- 1

## Define functions ---
# Function selects 2 random sequences from dataframe 
Select_Random_Seqs <- function(MidAtlantic_seq){
  # Create a vector of the unique sequence names
  MidAtlantic_names <- unique(MidAtlantic_seq$Name )
  
  # Randomly select 2 names from MidAtlantic_names
  selected_names <- sample(MidAtlantic_names, 2)
  
  # Extract all entries in MidAtlantic_seq where Name matches one of the selected names
  Seq_1 <- MidAtlantic_seq %>% filter(Name == selected_names[1]) 
  Seq_2 <- MidAtlantic_seq %>% filter(Name == selected_names[2])
  
  # Check if both Seq_1 and Seq_2 have 8 rows each
  if(nrow(Seq_1) != 8 | nrow(Seq_2) != 8){
    stop("Error: Seq_1 or Seq_2 does not have 8 entries")
  }
  
  # Sort Seq_1 and Seq_2 according to the defined seg_order
  Seq_1_sorted <- Seq_1 %>% 
    arrange(factor(Seg, levels = segment_name))
  
  Seq_2_sorted <- Seq_2 %>% 
    arrange(factor(Seg, levels = segment_name))
  
  Seqs <- list(Seq_1_sorted,Seq_2_sorted)
  return(Seqs)
}

# Function to align random sequences and return the number of errors
Get_Mismatches <- function(Seq_1,Seq_2){ 
  # Initialize an empty list to store the alignments
  aligned_sequences_decipher <- list()
  mismatch_counts <- numeric(8)
  
  # Loop through the sorted sequences to align them
  sg <- 1
  for(sg in 1:8){
    
    # Create a DNAStringSet object for the sequences to be aligned
    seqs <- DNAStringSet(c(Seq_1$Seq[sg], Seq_2$Seq[sg]))
    
    # Perform the alignment
    alignment <- invisible(AlignSeqs(seqs,verbose = FALSE))
    
    # Store the alignment in the list
    aligned_sequences_decipher[[sg]] <- alignment
    
    # Count the mismatches
    seq1_aligned <- as.character(alignment[1])
    # seq1_trimmed <- TrimSeq(seq1_aligned)
    
    seq2_aligned <- as.character(alignment[2])
    # seq2_trimmed <- TrimSeq(seq2_aligned)
    
    # Create vectors from the aligned strings
    vec1 <- unlist(strsplit(seq1_aligned, ''))
    vec2 <- unlist(strsplit(seq2_aligned, ''))
    
    # Filter out positions where either sequence has a gap
    non_gap_positions <- which(vec1 != '-' & vec2 != '-')
    
    # Count mismatches at non-gap positions
    mismatches <- sum(vec1[non_gap_positions] != vec2[non_gap_positions])
    
    mismatch_counts[sg] <- mismatches
  }
  
  # Names from Seq_1 and Seq_2
  seq1_name <- Seq_1$Name[1]
  seq2_name <- Seq_2$Name[1]
  
  # Create a new row as a data frame
  new_row <- data.frame("Seq1_name" = seq1_name, 
                        "Seq2_name" = seq2_name,
                        "PB2" = mismatch_counts[1], 
                        "PB1" = mismatch_counts[2], 
                        "PA" = mismatch_counts[3], 
                        "HA" = mismatch_counts[4], 
                        "NP" = mismatch_counts[5], 
                        "NA" = mismatch_counts[6], 
                        "MP" = mismatch_counts[7], 
                        "NS" = mismatch_counts[8])
  
  return(new_row)
}

## Load sequences -----
Seasonal_seq <- data.frame()

fasta_name_this <- file.path(Outlier_analysis_path,paste0("gisaid_epiflu_sequence_complete.fasta"))

fasta_name_this
tmp_fasta <- phylotools::read.fasta(fasta_name_this) # Phylotools package

# Step 1: Generate a list by splitting the seq.name
list_seq_name <- strsplit(as.character(tmp_fasta$seq.name), split = ";")

# Step 2: Extract individual components from the list
Name <- sapply(list_seq_name, `[`, 1)
EpiID <- sapply(list_seq_name, `[`, 2)
Seg <- sapply(list_seq_name, `[`, 3)
Date <- sapply(list_seq_name, function(x) x[length(x)])

# Step 3: Create a new dataframe
new_fasta_df <- data.frame(Name, Seg, EpiID, Date, Seq = tmp_fasta$seq.text)
Seasonal_seq <- rbind(Seasonal_seq,new_fasta_df)


# Run Code -----
#Create a dataframe to hold the sequence mismatches
mismatch_df <- data.frame("Seq1_name"=character(),"Seq2_name"=character(),"PB2"=double(),"PB1"=double(),"PA"=double(),"HA"=double(),"NP"=double(),"NA"=double(),"MP"=double(),"NS"=double())

for (rep in 1:100000){
  print(rep)
  Seqs <- Select_Random_Seqs(Seasonal_seq) # Function returns 2 randomly selected sorted sequences
  Seq_1 <- Seqs[[1]]
  Seq_2 <- Seqs[[2]]

  Mismatches_this <- Get_Mismatches(Seq_1, Seq_2) # Function determines the number of mismatches per segment

  # Optionally, you can combine mismatch counts into a data frame
  mismatch_df <- rbind(mismatch_df,Mismatches_this)
  mismatch_df
}

write.csv(x = mismatch_df,file = file.path(Outlier_analysis_path,"mismatch_df.csv"))
mismatch_df <- read.csv(file.path(Outlier_analysis_path,"mismatch_df.csv"))
colnames(mismatch_df)[colnames(mismatch_df) == "NA."] <- "NA_"
mismatch_df$X <- NULL

# Column names for gene segments
seg_names <- c("PB2", "PB1", "PA", "HA", "NP", "NA_", "MP", "NS")

# Initialize a vector to store the indices of rows to be removed
rows_to_remove <- c()

# Loop over each gene segment by index
i <- 1
for (i in 1:length(seg_names)) {
  # Assign the current segment name to a variable for debugging
  seg <- seg_names[i]
  seg
  # Calculate mean and standard deviation
  seg_mean <- mean(mismatch_df[, seg], na.rm = TRUE)
  seg_sd <- sd(mismatch_df[, seg], na.rm = TRUE)
  print(c(seg,seg_mean,seg_sd))
  hist(mismatch_df[,seg])
  
  # Calculate the threshold for removal
  removal_threshold <- 4 * seg_sd + seg_mean
  print(removal_threshold)

  # Find rows where the value for the segment is greater than the threshold
  rows_to_remove_seg <- which(mismatch_df[, seg] > removal_threshold)
  length(rows_to_remove_seg)
  
  # Append these row indices to the overall vector of rows to remove
  rows_to_remove <- c(rows_to_remove, rows_to_remove_seg)
}

# Remove duplicate row indices
rows_to_remove <- unique(rows_to_remove)

# Remove rows from dataframe
filtered_mismatch_df <- mismatch_df[-rows_to_remove, ]

# Reorder columns to have HA as the first column
ordered_cols <- c("HA", "PB2", "PB1", "PA", "NP", "NA_", "MP", "NS")
filtered_mismatch_df <- filtered_mismatch_df[, c("Seq1_name", "Seq2_name", ordered_cols)]

# Reshape the data from wide to long format
long_df <- tidyr::gather(filtered_mismatch_df, Gene, Mutations, HA:NS)

# Specify the order of the Gene levels
long_df$Gene <- factor(long_df$Gene, levels = ordered_cols)


# Create a custom labeller function
custom_labeller <- function(x) {
  ifelse(x == "NA_", "NA", x)
}


# Create the histograms
ggplot(long_df, aes(x = Mutations)) +
  geom_histogram(binwidth = 1, fill = "#333333", alpha = 0.7) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free", switch = "y", labeller = as_labeller(custom_labeller)) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0,15000)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  ) +
  xlab("Number of Mutations") +
  ylab("Frequency") 


## Determine distributions for # of HA mutations -----
# Extract rows where the number of HA mismatches equals 21
HA_21_df <- dplyr::filter(filtered_mismatch_df, HA == 21)
HA_21_df <- HA_21_df[, c("Seq1_name", "Seq2_name", ordered_cols)]
long_21_df <- tidyr::gather(HA_21_df, Gene, Mutations, ordered_cols)

# Create the histograms
ggplot(long_21_df, aes(x = Mutations)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free", switch = "y") +
  coord_cartesian(xlim = c(0, 40)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  ) +
  xlab("Number of Mutations") +
  ylab("Frequency")

# Initialize list to store CI and mean results
CI_results <- list()

# Loop through each column by index for debugging
for (i in 3:ncol(HA_21_df)) {
  # Extract the specific column
  col_data <- HA_21_df[, i]
  
  # Calculate 2.5% and 97.5% quantiles
  lower_bound <- quantile(col_data, 0.025)
  upper_bound <- quantile(col_data, 0.975)
  
  # Calculate mean
  col_mean <- mean(col_data)
  
  # Store in list in the desired order
  CI_results[[colnames(HA_21_df)[i]]] <- c("Lower" = lower_bound, "Mean" = col_mean, "Upper" = upper_bound)
}

# Convert list to dataframe for easier viewing
CI_df <- do.call(rbind, CI_results)
print(CI_df)

# Filter the original dataframe to include only rows where HA mutations are between 17 and 23
HA_17_21_df <- subset(filtered_mismatch_df, HA >= 17 & HA <= 21)

# Reshape the data from wide to long format
long_17_21_df <- tidyr::gather(HA_17_21_df, Gene, Mutations, HA:NS)

# Specify the order of the Gene levels
long_17_21_df$Gene <- factor(long_17_21_df$Gene, levels = ordered_cols)
long_17_21_df

# Initialize list to store CI and mean results
CI_results_17_21 <- list()

# Loop through each column by index for debugging
for (i in 3:ncol(HA_17_21_df)) {
  # Extract the specific column
  col_data <- HA_17_21_df[, i]
  col_data
  
  # Calculate 2.5% and 97.5% quantiles
  lower_bound <- quantile(col_data, 0.025)
  upper_bound <- quantile(col_data, 0.975)
  lower_bound
  upper_bound
  # Calculate mean
  col_mean <- mean(col_data)
  col_mean
  # Store in list in the desired order
  CI_results_17_21[[colnames(HA_17_21_df)[i]]] <- c("Lower 2.5%" = lower_bound, "Mean" = col_mean, "Upper 97.5%" = upper_bound)
}
CI_results_17_21


# Create the histograms
ggplot(long_17_21_df, aes(x = Mutations)) +
  geom_histogram(binwidth = 1, fill = "#333333", alpha = 0.7) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free", switch = "y") +
  coord_cartesian(xlim = c(0, 40), ylim = c(0,600)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  ) +
  xlab("Number of Mutations") +
  ylab("Frequency") + 
  ggtitle("Mutation distribution of Influenza Sequences wher nMutations(HA) ranges between 17-21")


# Define the desired order
desired_order <- c("HA", "PB2", "PB1", "PA", "NP", "NA_", "MP", "NS")


Muts_101 <- data.frame(Segment = desired_order, NMut = c(21,2,8,1,9,3,7,1))
Muts_101

# Combine the long_17_21_df with the reference lines from Muts_101
combined_df <- left_join(long_17_21_df, Muts_101, by = c("Gene" = "Segment"))


# Re-factorize the Gene column with the desired order
combined_df$Gene <- factor(combined_df$Gene, levels = desired_order)

# Convert list to dataframe for easier viewing
CI_df_17_21 <- do.call(rbind, CI_results_17_21)
CI_df_17_21

# Convert the matrix to a dataframe
CI_df_17_21 <- as.data.frame(CI_df_17_21)
# Rename the columns
colnames(CI_df_17_21) <- c("Lower", "Mean", "Upper")
# Add Gene column
CI_df_17_21$Gene <- rownames(CI_df_17_21)
# Update the Gene name for NA
#CI_df_17_21$Gene <- gsub("NA", "NA_", CI_df_17_21$Gene)



# Join CI_df_17_21 to combined_df
combined_df <- left_join(combined_df, CI_df_17_21, by = "Gene")


# Create a custom labeller function
custom_labeller <- function(x) {
  ifelse(x == "NA_", "NA", x)
}

# Generate initial version of plot with the observed iSNVs shown in red

p <- ggplot(combined_df, aes(x = Mutations)) +
  geom_histogram(binwidth = 1, fill = "#333333", alpha = 0.7) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free", switch = "y", labeller = as_labeller(custom_labeller)) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0,600)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  ) +
  xlab("Number of Mutations") +
  ylab("Frequency") +
  geom_vline(data = unique(combined_df[,c("Gene", "NMut")]), aes(xintercept = NMut), color = "red", size = 1)

# Show the plot
print(p)


# Generate a second version of the plot where the 95% CI for all gene segments other than HA are shown with a dashed line

# Generated modified dataframe where the 95% CI values for the HA gene semgnet are removed
mod_combined_df <- combined_df
mod_combined_df <- mod_combined_df %>%
  mutate(
    Lower = ifelse(Gene == "HA", NaN, Lower),
    Upper = ifelse(Gene == "HA", NaN, Upper)
  )

mod_combined_df$Gene <- factor(mod_combined_df$Gene, levels = desired_order)


# Create the updated plot
p <- ggplot(mod_combined_df, aes(x = Mutations)) +
  geom_histogram(binwidth = 1, fill = "#333333", alpha = 0.7) +
  geom_vline(data = unique(mod_combined_df[, c("Gene", "NMut")]), aes(xintercept = NMut), color = "red", size = 1) +
  geom_vline(data = unique(mod_combined_df[, c("Gene", "Lower")]), aes(xintercept = Lower), color = "black", size = 0.5, linetype = "dashed" ) +
  geom_vline(data = unique(mod_combined_df[, c("Gene", "Upper")]), aes(xintercept = Upper), color = "black", size = 0.5, linetype = "dashed" ) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free", switch = "y", labeller = as_labeller(function(x) ifelse(x == "NA_", "NA", x))) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 600)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, size = 14),  # Slightly increased font size for facet labels
    text = element_text(size = 12),  # Slightly increased base font size
    panel.spacing.y = unit(1, "cm")  # Increased space between vertically adjacent plots
  ) +
  xlab("Number of Mutations") +
  ylab("Frequency") 

# Show the plot
print(p)

