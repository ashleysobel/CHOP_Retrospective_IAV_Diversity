
# Prepare the run environment ------

# Clear environment 
rm(list = ls())

# Get path to working directory 
R_path <- getwd()
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")

# Load packages --------
library(data.table) # for the data.table function
library(ggplot2) # for ggplot, geom_point, geom_smooth, scale_x_continuous, scale_y_continuous, theme_classic, theme, element_rect, element_blank, element_text, ggsave
library(dplyr) # for mutate, case_when, group_by, summarise, pivot_wider
library(tidyr) # for bind_rows
library(cowplot)
library(reshape2)
library(MASS)
library(genbankr)
library(data.table)
library(trackViewer)
library(GenomicRanges)
library(svglite) # for SetEPS
library(Biostrings)



# Define functions ---------

# The get_pileup function processes sample-specific pileup data: it loads the
# appropriate pileup file, partitions the data based on sequence names, then
# reshapes each partition by replacing non-standard nucleotides, aggregating
# counts, and calculating the depth at each position. It returns a list of these
# reshaped data tables.
get_pileup <- function(sample){
  # Extract relevant data from the sample input.
  sample_name <- sample$Sample 
  seq_date <- sample$SeqDate
  
  # Construct the path to the pileup data files.
  pileup_path <- file.path(R_path,paste0(sample$SeqDate,"_PipelineOutput"),"Pileup")
  
  # Search for the appropriate pileup file in the given directory.
  pileup_file <- list.files(path=pileup_path,pattern = sample$Sample,full.names = TRUE)
  
  # Load the pileup file as a data frame.
  working_sample_pileup <- read.csv(pileup_file)
  working_sample_pileup$X <- NULL
  working_sample_pileup <- data.table::data.table(working_sample_pileup)
  
  # Create a new column 'seqnum' that contains the last character of each entry in the 'seqnames' colum
  working_sample_pileup[, seqnum := substr(seqnames, nchar(seqnames), nchar(seqnames))]
  
  # Partition the pileup data into 8 separate data tables based on 'seqnum'.
  PB2_working_pileup <- working_sample_pileup[seqnum == "1"]
  PB1_working_pileup <- working_sample_pileup[seqnum == "2"]
  PA_working_pileup <- working_sample_pileup[seqnum == "3"]
  HA_working_pileup <- working_sample_pileup[seqnum == "4"]
  NP_working_pileup <- working_sample_pileup[seqnum == "5"]
  NA_working_pileup <- working_sample_pileup[seqnum == "6"]
  MP_working_pileup <- working_sample_pileup[seqnum == "7"]
  NS_working_pileup <- working_sample_pileup[seqnum == "8"]
  
  # Combine the separated data tables into a list for easy processing.
  pileup_list <- list(PB2_working_pileup, PB1_working_pileup, PA_working_pileup, 
                      HA_working_pileup, NP_working_pileup, NA_working_pileup, 
                      MP_working_pileup, NS_working_pileup)
  
  # Apply a function over each data table in the list that reshapes the data,
  # replaces non-standard nucleotides with "N", and calculates the depth at each
  # position.
  pileup_list_reshaped <- lapply(pileup_list, function(pileup_df) {
    pileup_df %>%
      dplyr::mutate(nucleotide = dplyr::case_when(
        nucleotide %in% c("A", "T", "G", "C") ~ nucleotide,
        TRUE ~ "N"
      )) %>%
      dplyr::group_by(pos, nucleotide) %>%
      dplyr::summarise(count = sum(count)) %>%
      pivot_wider(names_from = nucleotide, values_from = count, values_fill = 0) %>%
      dplyr::mutate(Depth = A + T + G + C)
  })
  # Return the reshaped pileup data.
  return(pileup_list_reshaped)
}

# The get_PairwiseDistance function takes a list of pileup data frames, combines
# them into a single tibble, calculates a pairwise distance measure D_l for each
# position based on the counts of each nucleotide and depth, and returns the
# updated tibble.
get_PairwiseDistance <- function(pileup_tmp){
  # Combine all data frames in the input list into a single tibble.
  pileup_tibble <- dplyr::bind_rows(pileup_tmp)
  
  # Calculate Di for each nucleotide 'A', 'T', 'G', 'C'.
  pileup_tibble <- pileup_tibble %>%
    dplyr::mutate(
      D_lA = A * (A - 1),
      D_lT = T * (T - 1),
      D_lG = G * (G - 1),
      D_lC = C * (C - 1)
    )
  
  
  # Calculate D_l, the proportion of mismatches at each position.
  pileup_tibble <- pileup_tibble %>%
    dplyr::mutate(
      D_matches = D_lA + D_lT + D_lG + D_lC
    )
  
  
  # Calculate D_l, the proportion of mismatches at each position.
  pileup_tibble <- pileup_tibble %>%
    dplyr::mutate(
      D_l = (Depth * (Depth - 1) - D_matches) / (Depth * (Depth - 1))
    )
  
  # Return the updated tibble
  return(pileup_tibble)
}

# The get_ExpectedDiversity function calculates the expected diversity ED_l at
# each position in a given diversity tibble, based on the frequencies of each
# nucleotide at that position, and returns the updated tibble.
get_ExpectedDiversity <- function(diversity_tibble){
  diversity_tibble <- diversity_tibble %>%
    dplyr::mutate(
      pG2 = (G/Depth)^2,
      pT2 = (T/Depth)^2,
      pA2 = (A/Depth)^2,
      pC2 = (C/Depth)^2
    )
  diversity_tibble <- diversity_tibble %>%
    dplyr::mutate(
      pi_l = pG2 + pT2 + pA2 + pC2) # pi_l nulcotide frequencies at site l
  
  diversity_tibble <- diversity_tibble %>%
    dplyr::mutate(
      ED_l = 1 - (pG2 + pT2 + pA2 + pC2)^2) # ED_l = expected diversity at site l
  return(diversity_tibble)
}


# Evaluate mutation distribution -------
generate_mut_df <- function(summary_df) {
  # Convert the table to a data frame
  mut_df <- as.data.frame.table(summary_df)
  
  names(mut_df) <- c("seg_num", "iSNV")
  
  # Convert segment numbers to numeric
  mut_df$seg_num <- as.numeric(as.character(mut_df$seg_num))
  
  # Add missing segments
  all_segments <- 1:8
  missing_segments <- setdiff(all_segments, mut_df$seg_num)
  
  # Create data frame for missing segments
  missing_data <- data.frame(seg_num = missing_segments, iSNV = rep(0, length(missing_segments)))
  
  # Combine the original and missing data
  mut_df <- rbind(mut_df, missing_data)
  
  # Order the data frame by seg_num
  mut_df <- mut_df[order(mut_df$seg_num), ]
  mut_df
  # Return the data frame
  return(mut_df)
}

# Define parameters ----
Analysis_date <- "81523"

# Load data ---------------------
# Load finalized parameters from "Compile_data.R")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
CHOA_RunParam <- read.csv( paste0(Compiled_OutPut_path,"/CHOA_RunParam.csv"))
CHOA_RunParam
min_var <- CHOA_RunParam$min_var
#min_var_replicate <- CHOA_RunParam$min_var_replicate
min_var_replicate <- CHOA_RunParam$min_var_replicate
min_cov <- CHOA_RunParam$min_cov

# Load mutations
CHOA_Muts_FINAL <- read.csv(paste0(Compiled_OutPut_path,"/FINAL_MutationOutput_minvar3.csv"))

# Load merged data
merged_data <- read.csv(paste0(Compiled_OutPut_path,"/FINAL_CHOA_merged_data.csv"))

# Load sample data table 
sample_DataTable <- read.csv(paste0(Compiled_OutPut_path,"/","FinalMuts-CHOA_SampleDataTable.csv"))

# Calculate summary statistics (mean/variance) for iSNVs
mean_Syn_iSNV <- mean(sample_DataTable$Synonymous_iSNV)
std_Syn_iSNV <- sd(sample_DataTable$Synonymous_iSNV)
mean_nonSyn_iSNV <- mean(sample_DataTable$nonSynonymous_iSNV)
std_nonSyn_iSNV <- sd(sample_DataTable$nonSynonymous_iSNV)

print(mean_Syn_iSNV)
print(std_Syn_iSNV)
print(mean_nonSyn_iSNV)
print(std_nonSyn_iSNV)

# Calculate diversity -----------

get_Diversity_SampleSet <- function(Original_Samples) {
  # Create an empty dataframe to store L_segment and D_sample for each segment
  div_sample <- data.frame()
  
  # Create an empty list to store pi_s and Epi_s for each sample
  Div_s_list <- list()
  os <- 1
  # Loop over each original sample
  for (os in 1:nrow(Original_Samples)){
    
    sample_tmp <- Original_Samples[os,]
    sample_pileup <- get_pileup(sample_tmp)
    
    # Create temporary vectors to store L_segment and D_sample for each segment
    L_segments <- vector()
    D_segments <- vector()
    ED_segments <- vector()
    ls <- 1
    for (ls in 1:length(sample_pileup)){
      pileup_tmp <- sample_pileup[ls]
      diversity_tibble <- get_PairwiseDistance(pileup_tmp)
      diversity_tibble <- get_ExpectedDiversity(diversity_tibble)
      # Compute L_segment and D_sample for the current segment and store in the temporary vectors
      L_segments <- c(L_segments, sum(!is.na(diversity_tibble$D_l))) # summing only non-NA D_l values
      D_segments <- c(D_segments, sum(diversity_tibble$D_l, na.rm = TRUE)) # summing D_l values excluding NA's
      ED_segments <- c(ED_segments,sum(diversity_tibble$ED_l, na.rm = TRUE)) # summing ED_l values excluding NA's
    }
    
    # Save L_segment and D_sample for each segment in the div_sample dataframe
    div_sample <- rbind(div_sample, data.frame(L_segments = L_segments, D_segments = D_segments,ED_segments = ED_segments))
    
    # Calculate pi_s and Epi_s for the current sample and store in a Diversity list
    pi_s <- sum(D_segments) / sum(L_segments) # pi_s = pi (nucleotide diversity) for segment i
    Epi_S <- sum(ED_segments) / sum(L_segments)  # Epi_s = expected value of pi (uses expectation value for D_i from Zhao & Illingworth 2019)
    Div_s_list[[os]] <- cbind(sample_tmp$ID, pi_s, Epi_S)
    
  }
  return(Div_s_list)
}

# # Get Sample Diversity for the original sample
# Original_diversity <- get_Diversity_SampleSet(Original_Samples)
# Replicate_diversity <- get_Diversity_SampleSet(Replicate_Samples)
# 
# write.csv(Original_diversity,"Original_diversity.csv")
# write.csv(Replicate_diversity,"Replicate_diversity.csv")
# 
# Original_diversity <- read.csv("Original_diversity.csv")
# Original_diversity
# Replicate_diversity <- read.csv("Replicate_diversity.csv")
# # Change output for pi_s and Epi_s into dataframes: 
# 
# # Bind all list elements into a data frame
# original_diversity_df <- as.data.frame(do.call(rbind, Original_diversity))
# replicate_diversity_df <- as.data.frame(do.call(rbind, Replicate_diversity))
# 
# # Convert all elements to appropriate data types
# original_diversity_df[,1] <- as.character(original_diversity_df[,1])
# original_diversity_df[,2] <- as.numeric(original_diversity_df[,2])
# original_diversity_df[,3] <- as.numeric(original_diversity_df[,3])
# 
# replicate_diversity_df[,1] <- as.character(replicate_diversity_df[,1])
# replicate_diversity_df[,2] <- as.numeric(replicate_diversity_df[,2])
# replicate_diversity_df[,3] <- as.numeric(replicate_diversity_df[,3])
# 
# # Set column names
# colnames(original_diversity_df) <- c("ID", "Opi_s","OEpi_s") # Opi_s = original pi for sample s
# colnames(replicate_diversity_df) <-c("ID", "Rpi_s","REpi_s") # Rpi_s = replicate pi for sample s
# 
# # Merge original and replicate data into a single dataframe 
# Diversity_data <- merge(original_diversity_df,replicate_diversity_df, by = "ID")
# 
# # Add columns for the mean pi_s and Epi_s for the original and replicate samples
# Diversity_data <- Diversity_data %>%
#   dplyr::mutate(
#     Avg_pi_s = (Opi_s + Rpi_s) / 2,
#     Avg_Epi_s = (OEpi_s + REpi_s) / 2
#   )
# 
# 
# # Plot the replicate (later sampling date) diversity and expected diversity for each subject 
# ggplot(Diversity_data, aes(x=ID, y=Rpi_s)) +
#   geom_col(fill="steelblue") +
#   theme_minimal() +
#   labs(x="ID", y="pi", title="Column plot of pi for each ID") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# ggplot(Diversity_data, aes(x=ID, y=REpi_s)) +
#   geom_col(fill="steelblue") +
#   theme_minimal() +
#   labs(x="ID", y="pi", title="Column plot of pi for each ID") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # Merge diversity dataframes (pi_s_df and Epi_s_df) with merged_data
# merged_data <- merge(merged_data, Diversity_data, by.x = "SubjectID", by.y = "ID")

# Generate a table describing nucleotide diversity per sample for # iSNVs,
# nucleotide diversity (pi) and expected nucleotide diversity (Epi)
Diversity_dataframe <- data.frame(ID = merged_data$SubjectID, Age=merged_data$Age_mo/12, Vaccine=merged_data$Flu_vaccine, Strain = merged_data$Strain, CT = merged_data$CT, All_iSNV = merged_data$all_iSNV_present, 
                                  Syn_iSNV = merged_data$Syn_iSNV_present, nonSyn_iSNV = merged_data$nonSyn_iSNV_present, nonSyn_HA_iSNV = merged_data$nonSyn_HA_present)

# Check to see if there is a difference between H1N1 and H3N2 iSNVs by fitting a negative binonial regression model
# For each iSNV category, fit a Negative Binomial regression model
model_syn = MASS::glm.nb(Syn_iSNV ~ Strain, data = Diversity_dataframe)
model_nonSyn = MASS::glm.nb(nonSyn_iSNV ~ Strain, data = Diversity_dataframe)

# Now conduct an analysis of deviance on each model
deviance_syn = anova(model_syn, test = "Chisq")
deviance_nonSyn = anova(model_nonSyn, test = "Chisq")

# Account for multiple comparisons using the bonferroni method
pval_corrected_syn = min(1, deviance_syn$`Pr(>Chi)`[2] * 2)
pval_corrected_nonSyn = min(1, deviance_nonSyn$`Pr(>Chi)`[2] * 2)


# Print the results with and without accounting for multiple comparisons
print(deviance_syn)
print(deviance_nonSyn)
print(pval_corrected_syn)
print(pval_corrected_nonSyn)

# Reshape data to organize by iSNV type
long_df <- reshape2::melt(Diversity_dataframe, id.vars = c("ID", "Age", "Vaccine", "Strain", "CT", "All_iSNV", "nonSyn_HA_iSNV"), variable.name = "iSNV_Type", value.name = "Count")

# Generate plot for the distribution of iSNV and boxplot for the association of
# strain and iSNV
iSNV_distribution <- ggplot2::ggplot(long_df[long_df$iSNV_Type %in% c("Syn_iSNV", "nonSyn_iSNV"), ], aes(x = Count, fill = iSNV_Type)) +
  geom_histogram(position = "identity", binwidth = 1) +
  labs(x = "iSNV count", y = "Frequency", fill = "iSNV Type") +
  theme_minimal() + 
  theme(text = element_text(family = "Arial", size = 12),
        axis.text = element_text(size = 10))

# Save the plot
ggplot2::ggsave(file.path(Figure_OutPut_path,"iSNV_distribution.eps"), plot = iSNV_distribution, width = 6, height = 4, device = "eps")

print(iSNV_distribution)

# Creating a new variable for interaction of Strain and iSNV_Type
long_df$Strain_iSNV_Type <- interaction(long_df$iSNV_Type, long_df$Strain, sep = ", ")

# Ordering factor levels
long_df$Strain_iSNV_Type <- factor(long_df$Strain_iSNV_Type, levels = c("Syn_iSNV, H1N1", "Syn_iSNV, H3N2", "nonSyn_iSNV, H1N1", "nonSyn_iSNV, H3N2"))

# Create the boxplot
Strain_iSNV_association <- ggplot2::ggplot(long_df[long_df$iSNV_Type %in% c("Syn_iSNV", "nonSyn_iSNV"), ], aes(x = Strain_iSNV_Type, y = Count, fill = Strain_iSNV_Type)) +
  geom_boxplot() +
  labs(x = "iSNV Type and Strain", y = "Count", fill = "iSNV Type and Strain") +
  theme(text = element_text(family = "Arial", size = 12),
        axis.text = element_text(size = 10)) + 
  theme_minimal()

# Print the plot and export
print(Strain_iSNV_association)
ggplot2::ggsave(file.path(Figure_OutPut_path,"Strain_iSNV_association.eps"), plot = Strain_iSNV_association, width = 6, height = 4, device = "eps")


# Categories the Diversity information by age and vaccine status
Diversity_dataframe <- Diversity_dataframe %>%
  dplyr::mutate(
    AgeCategory = factor(
      dplyr::case_when(
        Age >= 0.5 & Age < 5   ~ "0 - 4",
        Age >= 5 & Age < 12    ~ "5 - 11",
        Age >= 12 & Age < 19   ~ "12 - 18",
        TRUE                   ~ "Other"
      ), 
      levels = c("0 - 4", "5 - 11", "12 - 18", "Other")),
    Vaccine = ifelse(Vaccine == "Unknown", "No", Vaccine)
  )

Diversity_dataframe$Group <- factor(interaction(Diversity_dataframe$Vaccine, Diversity_dataframe$AgeCategory), 
                                    levels = c("No.0 - 4", "Yes.0 - 4", "No.5 - 11", "Yes.5 - 11", "No.12 - 18", "Yes.12 - 18", "No.Other", "Yes.Other"))
write.csv(Diversity_dataframe,"CHOP_DiversityDataframe.csv")
Diversity_dataframe <- read.csv("CHOP_DiversityDataframe.csv")
Diversity_dataframe$X <- NULL
print(Diversity_dataframe)




# Association plots ------
# Assess association between clinical metadata and # iSNV:
association_data <- data.frame(SubjectID = merged_data$SubjectID,CT = merged_data$CT,
                               Age_yo = merged_data$Age_mo/12, Vaccine = merged_data$Flu_vaccine,
                               Symptom_dates = merged_data$Symptom_days, PMCA = merged_data$PMCA,
                               Setting = merged_data$Setting, Oseltamivir=merged_data$Prior_tamiflu, Admit=merged_data$Admit,
                               Syn_iSNV = merged_data$Syn_iSNV_present, nonSyn_iSNV = merged_data$nonSyn_iSNV_present)
# Replace "Unknown" in "Vaccine" column with "No"
association_data$Vaccine[association_data$Vaccine == "Unknown"] <- "No"


## Assess associations:

# custom jitter function
custom_jitter <- function(width = 0.5) {
  function(x) {
    x + runif(length(x), min = -width/2, max = width/2)
  }
}

# Association plots ------
# Assess association between clinical metadata and # iSNV:
association_data <- data.frame(SubjectID = merged_data$SubjectID,CT = merged_data$CT,
                               Age_yo = merged_data$Age_mo/12, Vaccine = merged_data$Flu_vaccine,
                               Symptom_dates = merged_data$Symptom_days, PMCA = merged_data$PMCA,
                               Setting = merged_data$Setting, Oseltamivir=merged_data$Prior_tamiflu, Admit=merged_data$Admit,
                               Syn_iSNV = merged_data$Syn_iSNV_present, nonSyn_iSNV = merged_data$nonSyn_iSNV_present)
# Replace "Unknown" in "Vaccine" column with "No"
association_data$Vaccine[association_data$Vaccine == "Unknown"] <- "No"


## Assess associations

# Fit a Negative Binomial GLM
glm_model_Syn_iSNV <- MASS::glm.nb(Syn_iSNV ~ Vaccine + PMCA + Setting + Admit, 
                                   data = association_data)
glm_model_nonSyn_iSNV <- MASS::glm.nb(nonSyn_iSNV ~ Vaccine + PMCA + Setting + Admit, 
                                      data = association_data)
# Summary of the model
summary(glm_model_Syn_iSNV)
summary(glm_model_nonSyn_iSNV)

# Adjust P values for multiple comparisons
# Extract p-values from the models
pvals_Syn_iSNV <- summary(glm_model_Syn_iSNV)$coefficients[, "Pr(>|z|)"]
pvals_Syn_iSNV
pvals_nonSyn_iSNV <- summary(glm_model_nonSyn_iSNV)$coefficients[, "Pr(>|z|)"]
pvals_nonSyn_iSNV

# Apply Bonferroni correction
adjusted_pvals_Syn_iSNV  <- p.adjust(pvals_Syn_iSNV, method = "bonferroni")
adjusted_pvals_Syn_iSNV

adjusted_pvals_nonSyn_iSNV  <- p.adjust(pvals_nonSyn_iSNV, method = "bonferroni")
adjusted_pvals_nonSyn_iSNV

## Age vs iSNV -----

# Filter the data for regression

# Calculate r^2 for Sym_iSNV vs Age_yo
model_Syn_iSNV <- lm(Syn_iSNV ~ Age_yo, data = association_data)
summary_Syn_iSNV <- summary(model_Syn_iSNV)
r2_Syn_iSNV <- summary_Syn_iSNV$r.squared


# Calculate r^2 for nonSyn_iSNV vs Age_yo
model_nonSyn_iSNV <- lm(nonSyn_iSNV ~ Age_yo, data = association_data)
summary_nonSyn_iSNV <- summary(model_nonSyn_iSNV)
r2_nonSyn_iSNV <- summary_nonSyn_iSNV$r.squared
print(paste0("R-squared for nonSyn_iSNV vs Age_yo: ", r2_nonSyn_iSNV))

# Create scatter plot
age_iSNV_plot <- ggplot(association_data, aes(Age_yo)) +
  geom_point(data = association_data,
             aes(y = Syn_iSNV), color = "black", fill = "black", shape = 21, size = 1.25) +
  geom_point(data = association_data,
             aes(y = nonSyn_iSNV), color = "#258786", fill = "#258786", shape = 21, size = 1.25) +
  geom_smooth(data = association_data, method = "lm", aes(y = Syn_iSNV),
              linetype = "solid", color = "black", se = FALSE) +
  geom_smooth(data = association_data, method = "lm", aes(y = nonSyn_iSNV),
              linetype = "solid", color = "#258786", se = FALSE) +
  scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) + # change the range of the x-axis 
  scale_y_continuous(limits = c(-1, 40), expand = c(0, 0)) + # change the range of the y-axis
  labs(x = "Age (years)", y = "iSNV per Sample") + # label the x and y axis
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"), # set the global font size and type
        axis.title = element_text(size=12)) # set the font size for axis titles

# Show the plot
print(age_iSNV_plot)

# Save the plot
age_iSNV_name <- paste0("Age_iSNV_plot.eps")
ggplot2::ggsave(file.path(Figure_OutPut_path,age_iSNV_name), plot = age_iSNV_plot, width = 5, height = 3.5, device = "eps")

# Extract the regression model and calculate r^2 value
lm_model_syn_Age <- lm(Syn_iSNV ~ Age_yo, data = association_data)
r_squared_syn_Age <- summary(lm_model_syn_Age)$r.squared
print(lm_model_syn_Age)
print(r_squared_syn_Age)

lm_model_nonsyn_Age <- lm(nonSyn_iSNV ~ Age_yo, data = association_data)
r_squared_nonsyn_Age <- summary(lm_model_nonsyn_Age)$r.squared
print(lm_model_nonsyn_Age)
print(r_squared_nonsyn_Age)


## Generate plot for iSNV vs vaccine status ----

# Gather the data into long format
association_data_long <- tidyr::pivot_longer(association_data, 
                                             cols = c(Syn_iSNV, nonSyn_iSNV),
                                             names_to = "iSNV_type",
                                             values_to = "iSNV_values")

# Create a new column for color based on SubjectID
association_data_long$color <- ifelse(association_data_long$iSNV_type == "Syn_iSNV", "black", "#258786")

# Update levels of interaction to ensure correct order on x-axis
association_data_long$interaction_level <- interaction(association_data_long$iSNV_type, 
                                                       association_data_long$Vaccine, 
                                                       sep = ", ")
levels <- c("Syn_iSNV, No", "Syn_iSNV, Yes", "nonSyn_iSNV, No", "nonSyn_iSNV, Yes")
association_data_long$interaction_level <- factor(association_data_long$interaction_level, levels = levels)

# Generate plot
iSNV_vaccine_plot <- ggplot(association_data_long, aes(x = interaction_level, 
                                                       y = iSNV_values, 
                                                       color = iSNV_type,
                                                       shape = iSNV_type)) +
  geom_point(position = position_dodge2(width = 0.5), size = 1.25) +
  scale_color_manual(values = c("Syn_iSNV" = "black", "nonSyn_iSNV" = "#258786")) +
  scale_shape_manual(values = c("Syn_iSNV" = 16, "nonSyn_iSNV" = 16)) +
  labs(x = "iSNV Type, Vaccine Status", y = "iSNV per Sample") +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position="none", 
        panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"), # set the global font size and type
        axis.title = element_text(size=12)) + # set the font size for axis titles
  theme_bw() + theme(legend.position="none", panel.grid.minor = element_blank()) 

# Show the plot
print(iSNV_vaccine_plot)
vaccine_iSNV_name <- paste0("vaccine_iSNV_plot.eps")
ggplot2::ggsave(file.path(Figure_OutPut_path,vaccine_iSNV_name), plot = iSNV_vaccine_plot, width = 3.5, height = 3.5, device = "eps")


## Symptom Dates vs iSNV

# Decrease Symptom_dates by 1 day
association_data$Symptom_days <- association_data$Symptom_dates - 1

# Calculate r^2 for Syn_iSNV vs Symptom_days
model_Syn_iSNV <- lm(Syn_iSNV ~ Symptom_days, data = association_data)
summary_Syn_iSNV <- summary(model_Syn_iSNV)
r2_Syn_iSNV <- summary_Syn_iSNV$r.squared
print(paste0("R-squared for Syn_iSNV vs Symptom_days: ", r2_Syn_iSNV))

# Calculate r^2 for nonSyn_iSNV vs Symptom_days
model_nonSyn_iSNV <- lm(nonSyn_iSNV ~ Symptom_days, data = association_data)
summary_nonSyn_iSNV <- summary(model_nonSyn_iSNV)
r2_nonSyn_iSNV <- summary_nonSyn_iSNV$r.squared
print(paste0("R-squared for nonSyn_iSNV vs Symptom_days: ", r2_nonSyn_iSNV))

# Create scatter plot
SymptomDays_iSNV_plot <- ggplot(association_data, aes(Symptom_days)) +
  geom_point(data = association_data,
             aes(y = Syn_iSNV), color = "black", fill = "black", shape = 21, size = 1.25, position = position_dodge2(width = 0.5)) +
  geom_point(data = association_data,
             aes(y = nonSyn_iSNV), color = "#258786", fill = "#258786", shape = 21, size = 1.25, position = position_dodge2(width = 0.5)) +
  geom_smooth(data = association_data, method = "lm", aes(y = Syn_iSNV),
              linetype = "solid", color = "black", se = FALSE) +
  geom_smooth(data = association_data, method = "lm", aes(y = nonSyn_iSNV),
              linetype = "solid", color = "#258786", se = FALSE) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(x = "Days of Symptoms", y = "iSNV per Sample") + # label the x and y axis
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"), # set the global font size and type
        axis.title = element_text(size=12)) # set the font size for axis titles

# Show the plot
print(SymptomDays_iSNV_plot)
SymptomDays_iSNV_name <- paste0("SymptomDays_iSNV_plot.eps")

ggplot2::ggsave(file.path(Figure_OutPut_path,SymptomDays_iSNV_name), plot = SymptomDays_iSNV_plot, width = 5, height = 3.5, device = "eps")


## Generate plot for iSNV vs PMCA
# Gather the data into long format
association_data_long <- tidyr::pivot_longer(association_data, 
                                             cols = c(Syn_iSNV, nonSyn_iSNV),
                                             names_to = "iSNV_type",
                                             values_to = "iSNV_values")

# Create a new column for color based on SubjectID
association_data_long$color <- ifelse(association_data_long$iSNV_type == "Syn_iSNV", "black", "#258786")

# Update levels of interaction to ensure correct order on x-axis
association_data_long$interaction_level <- interaction(association_data_long$iSNV_type, 
                                                       association_data_long$PMCA, 
                                                       sep = ", ")
levels <- c("Syn_iSNV, N-CD", "Syn_iSNV, NC-CD", "Syn_iSNV, C-CD", "nonSyn_iSNV, N-CD", "nonSyn_iSNV, NC-CD","nonSyn_iSNV, C-CD")
association_data_long$interaction_level <- factor(association_data_long$interaction_level, levels = levels)

# Generate plot
iSNV_PMCA_plot <- ggplot(association_data_long, aes(x = interaction_level, 
                                                    y = iSNV_values, 
                                                    color = iSNV_type)) +
  geom_point(position = position_dodge2(width = 0.5), size = 1.25) +
  scale_color_manual(values = c("Syn_iSNV" = "black", "nonSyn_iSNV" = "#258786")) +
  labs(x = "iSNV Type, PMCA", y = "iSNV per sample") +
  scale_y_continuous(limits = c(0, 40)) +
  theme_bw() +
  theme(legend.position="none", 
        panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"), # set the global font size and type
        axis.title = element_text(size=12)) # set the font size for axis titles

# Show the plot
print(iSNV_PMCA_plot)
vaccine_iSNV_name <- paste0("iSNV_PMCA_plot.eps")
ggplot2::ggsave(file.path(Figure_OutPut_path,vaccine_iSNV_name), plot = iSNV_PMCA_plot, width = 3.5, height = 3.5, device = "eps")



# Arrange the final plots adjacent to each other
final_plot <- cowplot::plot_grid(age_iSNV_plot, iSNV_vaccine_plot, SymptomDays_iSNV_plot, iSNV_PMCA_plot, ncol = 2, nrow = 2, align = "v", rel_widths = c(1.5, 1))

print(final_plot)


# Save the plot
ggplot2::ggsave(file.path(Figure_OutPut_path,"final_plot.eps"), plot = final_plot, width = 8, height = 8, device = "eps")

merged_data <- merged_data[merged_data$SubjectID !%in% ]

## Create table of cohort distribution based on the categories established for the violin plots
Cohort_factors <- data.frame(SubjectID = merged_data$SubjectID, 
                             CT = merged_data$CT, Age = merged_data$Age_mo/12,
                             Symptoms = merged_data$Symptom_days,
                             Vaccine = merged_data$Flu_vaccine,
                             PMCA = merged_data$PMCA, CLD=merged_data$CLD, Asthma=merged_data$Asthma, IC=merged_data$Immunosuppression)

# Replace 2 with "No" and 3 with "Yes" in CLD, Asthma, and IC columns
Cohort_factors <- Cohort_factors %>% 
  dplyr::mutate(
    CLD = dplyr::case_when(
      CLD == 2 ~ "No",
      CLD == 3 ~ "Yes",
      TRUE ~ as.character(CLD)  # keeps the original value if it's not 2 or 3
    ),
    Asthma = dplyr::case_when(
      Asthma == 2 ~ "No",
      Asthma == 3 ~ "Yes",
      TRUE ~ as.character(Asthma)  # keeps the original value if it's not 2 or 3
    ),
    IC = dplyr::case_when(
      IC == 2 ~ "No",
      IC == 3 ~ "Yes",
      TRUE ~ as.character(IC)  # keeps the original value if it's not 2 or 3
    )
  )

# Create LungDisease column
Cohort_factors <- Cohort_factors %>% 
  dplyr::mutate(
    LungDisease = dplyr::case_when(
      CLD == "Yes" ~ "Yes",
      Asthma == "Yes" ~ "Yes",
      TRUE ~ "No"
    )
  )

# Categories the Diversity information by age and vaccine status
Cohort_factors <- Cohort_factors %>%
  dplyr::mutate(
    AgeCategory = factor(
      dplyr::case_when(
        Age >= 0.5 & Age < 5   ~ "0 - 4",
        Age >= 5 & Age < 12    ~ "5 - 11",
        Age >= 12 & Age < 19   ~ "12 - 18",
        TRUE                   ~ "Other"
      ), 
      levels = c("0 - 4", "5 - 11", "12 - 18", "Other")
    )
  )

Cohort_factors$Group <- factor(interaction(Cohort_factors$AgeCategory), 
                               levels = c("0 - 4", "5 - 11", "12 - 18"))

Cohort_factors$Vaccine <- ifelse(Cohort_factors$Vaccine == "Unknown", "No", Cohort_factors$Vaccine)

# Generate summary table for Cohort_factors dataframe by AgeCategory and calculates summary statistics for each group, including counts and proportions for specific conditions, and mean and standard deviation for CT and Age values.
DemoTab1 <- Cohort_factors %>% 
  dplyr::group_by(AgeCategory) %>% 
  dplyr::summarise(
    Entries = n(),
    Mean_CT = mean(CT, na.rm = TRUE),
    SD_CT = sd(CT, na.rm = TRUE),
    Mean_Symptoms= mean(Symptoms),
    SD_Symptoms = sd(Symptoms),
    Mean_Age = mean(Age, na.rm = TRUE),
    SD_Age = sd(Age, na.rm = TRUE),
    Vaccine_Yes = sum(Vaccine == "Yes"),
    Prop_Vaccine_Yes = mean(Vaccine == "Yes")*100,
    Vaccine_No = sum(Vaccine == "No"),
    Prop_Vaccine_No = mean(Vaccine == "No")*100,
    PMCA_N_CD = sum(PMCA == "N-CD"),
    Prop_N_CD = mean(PMCA == "N-CD")*100,
    PMCA_NC_CD = sum(PMCA == "NC-CD"),
    Prop_NC_CD = mean(PMCA == "NC-CD")*100,
    PMCA_C_CD = sum(PMCA == "C-CD"),
    Prop_C_CD = mean(PMCA == "C-CD")*100,
    LungDisease_Yes = sum(LungDisease == "Yes"),
    Prop_LungDisease_Yes = mean(LungDisease == "Yes")*100,
    IC_Yes = sum(IC == "Yes"),
    Prop_IC_Yes = mean(IC == "Yes")*100
  )

# Merge the relevant columns first
DemoTab1<- DemoTab1 %>%
  dplyr::mutate(
    CT = paste(round(Mean_CT, 1), " +/- ", round(SD_CT, 1), sep=""), Symptoms=paste(round(Mean_Symptoms,1), "+/-", round(SD_Symptoms,1), sep=""),
    Age = paste(round(Mean_Age, 1), " +/- ", round(SD_Age, 1), sep=""),
    Vaccine_Yes = paste(Vaccine_Yes, " (", round(Prop_Vaccine_Yes), "%)", sep=""),
    Vaccine_No = paste(Vaccine_No, " (", round(Prop_Vaccine_No), "%)", sep=""),
    LungDisease_Yes = paste(LungDisease_Yes, " (", round(Prop_LungDisease_Yes), "%)", sep=""),
    IC_Yes = paste(IC_Yes, " (", round(Prop_IC_Yes), "%)", sep=""),
    "N-CD" = paste(PMCA_N_CD, " (", round(Prop_N_CD), "%)", sep=""),
    "NC-CD" = paste(PMCA_NC_CD, " (", round(Prop_NC_CD), "%)", sep=""),
    "C-CD" = paste(PMCA_C_CD, " (", round(Prop_C_CD), "%)", sep="")
  )

# Create a new dataframe containing only the merged columns
CHOP_Demographics_Table_1 <- data.frame(Cat=DemoTab1$AgeCategory, N=DemoTab1$Entries,
                                        Age=DemoTab1$Age, CT=DemoTab1$CT, Symptoms = DemoTab1$Symptoms,
                                        No=DemoTab1$Vaccine_No,Yes=DemoTab1$Vaccine_Yes,
                                        'N-CD'=DemoTab1$`N-CD`,'NC-CD'=DemoTab1$`NC-CD`,
                                        'C-CD'=DemoTab1$`C-CD`,CLD=DemoTab1$LungDisease_Yes,IC=DemoTab1$IC_Yes)


write.csv(x = CHOP_Demographics_Table_1,file.path(Compiled_OutPut_path,"CHOP_Demographics_Table_1.csv"))

