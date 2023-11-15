# Prepare the run environment ------

# Clear environment 
rm(list = ls())

# Load packages
library(genbankr)
library(data.table)
library(trackViewer)
library(GenomicRanges)
library(svglite) # for SetEPS
library(Biostrings)

# Analysis date
Analysis_date <- Sys.Date()

# Get path to working directory 
R_path <- getwd()
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
Reference_file_path <- file.path(R_path,"Reference_files/")
SegmentSpecific_Muts_path <- file.path(R_path,"PostProcessing_OutPut","SegmentSpecific_Muts")

# Load Reference information
Muts_VPC <- read.csv(paste0(Compiled_OutPut_path,"/FINAL_MutationOutput_minvar3.csv"))
IAV_length <- read.csv(file.path(Reference_file_path,"IAV_segment_length.csv"))
gene_segment_names <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")

# Define functions --------

# Evaluate mutation distribution
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

# Function to perform chi-square tests
perform_chi_square_test <- function(df, category,gene_segment) {
  seg_index <- which(df$seg_name == gene_segment)
  non_seg_index <- which(df$seg_name != gene_segment)
  observed <- c(df[seg_index, category], sum(df[non_seg_index, category]))
  expected <- c(df[seg_index, paste0("expected_iSNV_", category)], sum(df[non_seg_index, paste0("expected_iSNV_", category)]))
  expected
  
  # Handling zero values in observed or expected
  valid_entries <- (observed >= 2) & (expected > 0)
  observed <- observed[valid_entries]
  expected <- expected[valid_entries]
  
  if(length(observed) >= 2) {
    chisq_test <- chisq.test(observed, p = expected / sum(expected))
    return(chisq_test$p.value)
  } else {
    return(NA) # Return NA if not enough valid entries for chi-squared test
  }
}

GetMutations <- function(Mutations, seg_num, genbank_ref, gene_segment_names){
  
  Muts_VPC <- Mutations
  segment_muts <- Muts_VPC[Muts_VPC$segment_num == seg_num,]
  
  nonSynonymous_seg_muts <- segment_muts[segment_muts$effect_any == "Nonsynonymous",]
  Synonymous_seg_muts <- segment_muts[segment_muts$effect_any == "Synonymous",]
  
  if(nrow(Synonymous_seg_muts) > 0){
    snp_freq_synonymous <- table(Synonymous_seg_muts$position)
    snp_gr_synonymous <- GenomicRanges::GRanges(
      gene_segment_names[seg_num], 
      IRanges(as.integer(names(snp_freq_synonymous)), width=1, names=rep("", length(snp_freq_synonymous))),
      color = rep("black", length(snp_freq_synonymous)), 
      score = as.integer(snp_freq_synonymous)
    )
  } else {
    snp_gr_synonymous <- GRanges()
  }
  
  if(nrow(nonSynonymous_seg_muts) > 0){
    snp_freq_nonsynonymous <- table(nonSynonymous_seg_muts$position)
    snp_gr_nonsynonymous <- GenomicRanges::GRanges(
      gene_segment_names[seg_num], 
      IRanges(as.integer(names(snp_freq_nonsynonymous)), width=1, names=rep("", length(snp_freq_nonsynonymous))),
      color = rep("#258786", length(snp_freq_nonsynonymous)), 
      score = as.integer(snp_freq_nonsynonymous)
    )
  } else {
    snp_gr_nonsynonymous <- GRanges()
  }
  
  snp_gr <- c(snp_gr_synonymous, snp_gr_nonsynonymous)
  
  return(snp_gr)
}

Get_chisq_iSNV_distribution <- function(all_iSNV_segment_summary, nonSynon_iSNV_segment_summary, IAV_length, gene_segment) {
  # Generate datatable with the mutation distribution per gene segment:
  iSNV_seg <- as.data.frame(all_iSNV_segment_summary); iSNV_seg
  nonSyn_seg <- as.data.frame(nonSynon_iSNV_segment_summary); nonSyn_seg
  iSNV_per_seg <- data.frame(seg_num=iSNV_seg$seg_num,seg_name=IAV_length$seg_name,seg_length=IAV_length$seg_length, all=iSNV_seg$iSNV,syn= iSNV_seg$iSNV - nonSyn_seg$iSNV, nonSyn = nonSyn_seg$iSNV); iSNV_per_seg
  iSNV_per_seg
  
  # Calculate the total number of mutations and the total length of all segments
  total_length <- sum(iSNV_per_seg$seg_length, na.rm = TRUE); total_length
  total_all <- sum(iSNV_per_seg$all, na.rm = TRUE); total_all
  total_syn <- sum(iSNV_per_seg$syn, na.rm = TRUE); total_syn
  total_nonSyn <- sum(iSNV_per_seg$nonSyn, na.rm = TRUE); total_nonSyn
  
  # Calculate expected mutations for each segment
  iSNV_per_seg$expected_iSNV_all <- (iSNV_per_seg$seg_length / total_length) * total_all
  iSNV_per_seg$expected_iSNV_syn <- (iSNV_per_seg$seg_length / total_length) * total_syn
  iSNV_per_seg$expected_iSNV_nonSyn <- (iSNV_per_seg$seg_length / total_length) * total_nonSyn
  
  # Perform chi-square tests and print p-values
  p_value_all <- perform_chi_square_test(iSNV_per_seg, "all",gene_segment)
  p_value_syn <- perform_chi_square_test(iSNV_per_seg, "syn",gene_segment)
  p_value_nonSyn <- perform_chi_square_test(iSNV_per_seg, "nonSyn",gene_segment)
  
  p_values <- data.frame(class = c("all","syn","nonSyn"),p_value=c(p_value_all,p_value_syn,p_value_nonSyn)); p_values
  return(p_values)
}



GetFeatures <- function(seg_num, genbank_ref, gene_features_palette, gene_segment_names){
  # Isolate reference file for target segment
  Ref_seg <- genbank_ref[[seg_num]]
  
  # Get features to annotate 
  gene_features <- Ref_seg@genes
  
  # Update the gene_features object so the seqnames match "entire_gene"
  seqlevels(gene_features) <- gene_segment_names[seg_num]
  seqnames(gene_features) <- gene_segment_names[seg_num]
  
  # Assign colors from the gene_features_palette to each feature
  feature_colors <- gene_features_palette[1:length(gene_features)]
  
  # Update fill and color attributes of gene features
  gene_features$fill <- feature_colors
  
  # Assign different heights to each feature
  # You can adjust the offset (0.02 in this example) as needed
  gene_features$height <- seq(from = 0.02, by = 0.02, length.out = length(gene_features))
  
  # Create a GRanges object for the entire segment
  entire_gene <- GRanges(
    seqnames = gene_segment_names[seg_num],
    ranges = IRanges::IRanges(start = 1, end = nchar(Ref_seg@sequence), names = paste0(gene_segment_names[seg_num]," gene")),
    fill = "#CCCCCC",  # a light gray color
    height = 0.02
  )
  
  # Separate the gene features into individual GRanges objects
  unique_genes <- unique(gene_features$gene_id)
  all_gene_features <- list()
  
  for(i in seq_along(unique_genes)) {
    temp_gene_features <- gene_features[gene_features$gene_id == unique_genes[i]]
    temp_gene_features$type <- paste0("annotation_track_", i)
    all_gene_features[[i]] <- temp_gene_features
  }
  
  # Combine all features into a single GRanges object
  all_features <- c(entire_gene, do.call(c, all_gene_features))
  
  # Set a default height for all features if it's NA
  all_features$height[is.na(all_features$height)] <- 0.02
  return(all_features)
}

Create_dandelion <- function(Muts_VPC, Strain, genbank_ref, gene_segment_names, gene_features_palette) {
  Muts <- Muts_VPC[Muts_VPC$Strain == Strain,]
  gs <- 2
  for (gs in 1:8){
    snp_gr <- GetMutations(Mutations = Muts, seg_num = gs, genbank_ref = genbank_ref, gene_segment_names )
    all_features <- GetFeatures(gs, genbank_ref = genbank_ref, gene_features_palette, gene_segment_names)
    
    seg_length <- nchar(genbank_ref[[gs]]@sequence)
    
    increments <- 500
    last_num <- ceiling(seg_length / increments) * increments
    dandy_x <- seq(0, last_num, by = increments)
    dandy_x[1] <- 1 
    dandy_x[length(dandy_x)] <- seg_length 
    
    par(mar=c(0.5, 0,0,.5))
    setEPS()
    postscript(file.path(Figure_OutPut_path,paste0("CHOP_", Strain, "_dandelion_", gene_segment_names[gs], ".eps")))
    tmp_dandelion <- trackViewer::dandelion.plot(snp_gr, all_features, ranges = GRanges(gene_segment_names[gs], IRanges::IRanges(start = 1, end = seg_length)), type="circle", newpage = TRUE, xaxis = dandy_x, cex = 0.8, ylab = "", maxgaps = 1/90)
    dev.off()
  }
}

calculate_mutations <- function(Muts_VPC, strain) {
  # Filter the data for the given strain
  Muts_VPC_strain <- Muts_VPC[Muts_VPC$Strain == strain,]
  
  # Determine the number of mutations per gene segment for all data and generate dataframe
  # Get table for synonymous mutations 
  synonymous_VPC <- Muts_VPC_strain[Muts_VPC_strain$effect_any == "Synonymous",]
  syn_iSNV_segment_summary <- table(synonymous_VPC$segment_num)
  
  Nonsynonymous_VPC <- Muts_VPC_strain[Muts_VPC_strain$effect_any == "Nonsynonymous",]
  nonSynon_iSNV_segment_summary <- table(Nonsynonymous_VPC$segment_num)
  
  syn_iSNV_segment_summary <- generate_mut_df(syn_iSNV_segment_summary)
  nonSynon_iSNV_segment_summary <- generate_mut_df(nonSynon_iSNV_segment_summary)
  
  all_syn <- sum(syn_iSNV_segment_summary$iSNV)
  all_nonsyn <- sum(nonSynon_iSNV_segment_summary$iSNV)
  
  # Determine the number of mutations per gene segment for all data excluding outliers CHOA-101 and CHOA-117
  outliers <- c("CHOA-101", "CHOA-117")
  Muts_sO <- Muts_VPC_strain[!Muts_VPC_strain$Sample_ID_mut %in% outliers,] # Muts_sO means muts sans outliers
  
  Syn_data <- Muts_sO[Muts_sO$effect_any == "Synonymous",]
  Syn_iSNV_s0_summary <- table(Syn_data$segment_num)
  nonSyn_data <- Muts_sO[Muts_sO$effect_any == "Nonsynonymous",]
  nonSynon_iSNV_s0_summary <- table(nonSyn_data$segment_num)
  
  syn_iSNV_sO_summary <- generate_mut_df(Syn_iSNV_s0_summary)
  nonSynon_iSNV_s0_summary <- generate_mut_df(nonSynon_iSNV_s0_summary)
  
  # Return the calculated summaries
  return(list(
    syn_iSNV_segment_summary = syn_iSNV_segment_summary,
    nonSynon_iSNV_segment_summary = nonSynon_iSNV_segment_summary,
    all_syn = all_syn,
    all_nonsyn = all_nonsyn,
    syn_iSNV_sO_summary = syn_iSNV_sO_summary,
    nonSynon_iSNV_s0_summary = nonSynon_iSNV_s0_summary
  ))
}

# Extracting the results for H3N2
results_H3N2 <- calculate_mutations(Muts_VPC, "H3N2")
syn_iSNV_segment_summary_H3N2 <- results_H3N2$syn_iSNV_segment_summary
nonSynon_iSNV_segment_summary_H3N2 <- results_H3N2$nonSynon_iSNV_segment_summary

all_iSNV_segment_summary_H3N2 <- syn_iSNV_segment_summary_H3N2
all_iSNV_segment_summary_H3N2$iSNV <- syn_iSNV_segment_summary_H3N2$iSNV + nonSynon_iSNV_segment_summary_H3N2$iSNV
all_syn_H3N2 <- results_H3N2$all_syn
all_nonsyn_H3N2 <- results_H3N2$all_nonsyn
syn_iSNV_sO_summary_H3N2 <- results_H3N2$syn_iSNV_sO_summary
nonSynon_iSNV_s0_summary_H3N2 <- results_H3N2$nonSynon_iSNV_s0_summary
all_iSNV_sO_summary_H3N2 <- syn_iSNV_sO_summary_H3N2
all_iSNV_sO_summary_H3N2$iSNV <- syn_iSNV_sO_summary_H3N2$iSNV + nonSynon_iSNV_segment_summary_H3N2$iSNV

all_syn_H3N2/8
all_nonsyn_H3N2/8
syn_iSNV_segment_summary_H3N2
nonSynon_iSNV_segment_summary_H3N2

# Extracting the results for H1N1
results_H1N1 <- calculate_mutations(Muts_VPC, "H1N1")
syn_iSNV_segment_summary_H1N1 <- results_H1N1$syn_iSNV_segment_summary
nonSynon_iSNV_segment_summary_H1N1 <- results_H1N1$nonSynon_iSNV_segment_summary
all_syn_H1N1 <- results_H1N1$all_syn
all_nonsyn_H1N1 <- results_H1N1$all_nonsyn
syn_iSNV_sO_summary_H1N1 <- results_H1N1$syn_iSNV_sO_summary
nonSynon_iSNV_s0_summary_H1N1 <- results_H1N1$nonSynon_iSNV_s0_summary

all_iSNV_segment_summary_H1N1 <- syn_iSNV_segment_summary_H1N1
all_iSNV_segment_summary_H1N1$iSNV <- syn_iSNV_segment_summary_H1N1$iSNV + nonSynon_iSNV_segment_summary_H1N1$iSNV
all_iSNV_segment_summary_H1N1

all_iSNV_sO_summary_H1N1 <- syn_iSNV_sO_summary_H1N1
all_iSNV_sO_summary_H1N1$iSNV <- syn_iSNV_sO_summary_H1N1$iSNV + nonSynon_iSNV_segment_summary_H1N1$iSNV

all_syn_H1N1/8
all_nonsyn_H1N1/8
syn_iSNV_segment_summary_H1N1
nonSynon_iSNV_segment_summary_H1N1

# Check for mutation distribution for "CHOA-117"
Muts_117_Syn <- Muts_VPC[Muts_VPC$Sample_ID_mut == "CHOA-117" & Muts_VPC$effect_any == "Synonymous",]
Muts_117_nonSyn <-Muts_VPC[Muts_VPC$Sample_ID_mut == "CHOA-117" & Muts_VPC$effect_any == "Nonsynonymous",]
Syn_iSNV_C117_summary <- table(Muts_117_Syn$segment_num)
nonSynon_iSNV_C117_summary <- table(Muts_117_nonSyn$segment_num)

Syn_iSNV_C117_summary
nonSynon_iSNV_C117_summary
Syn_iSNV_C117_summary <- generate_mut_df(Syn_iSNV_C117_summary)
nonSynon_iSNV_C117_summary <- generate_mut_df(nonSynon_iSNV_C117_summary)


# Compare HA iSNVs for all subjects, all subjects sans outliers, CHOA-101 and CHOA-117
p_values_all_sub_H3N2 <- Get_chisq_iSNV_distribution(all_iSNV_segment_summary_H3N2, nonSynon_iSNV_segment_summary_H3N2, IAV_length,"HA"); p_values_all_sub_H3N2
p_values_sO_H3N2 <- Get_chisq_iSNV_distribution(syn_iSNV_sO_summary_H3N2, nonSynon_iSNV_s0_summary_H3N2,IAV_length,"HA"); p_values_sO_H3N2

p_values_all_sub_H1N1 <- Get_chisq_iSNV_distribution(syn_iSNV_segment_summary_H1N1, nonSynon_iSNV_segment_summary_H1N1, IAV_length,"HA"); p_values_all_sub_H1N1
p_values_sO_H1N1 <- Get_chisq_iSNV_distribution(syn_iSNV_sO_summary_H1N1, nonSynon_iSNV_s0_summary_H1N1,IAV_length,"HA"); p_values_sO_H1N1

p_values_NA_sub_H1N1 <- Get_chisq_iSNV_distribution(syn_iSNV_segment_summary_H1N1, nonSynon_iSNV_segment_summary_H1N1, IAV_length,"NA"); p_values_NA_sub_H1N1

p_values_C101 <- Get_chisq_iSNV_distribution(Syn_iSNV_C101_summary,nonSynon_iSNV_C101_summary,IAV_length,"HA"); p_values_C101
p_values_C117 <- Get_chisq_iSNV_distribution(Syn_iSNV_C117_summary,nonSynon_iSNV_C117_summary,IAV_length,"HA"); p_values_C117


# Load your genbank file 
test <- genbankr::readGenBank(file.path(Reference_file_path,"A_Washington_17_2016_H3N2.gb"))

# Set color palette
gene_features_palette <- c("#91d7e2", "#3e93c0", "#40d3b0")

# We load NS gene segment reference separately 
test2 <- genbankr::readGenBank(file.path(Reference_file_path,"H3N2_ref_NS.gb"))
test[[8]] <- test2
H3N2_ref_genbank <- test

# Load genbank file for H1H1
test <- genbankr::readGenBank(file.path(Reference_file_path,"A_Michigan_45_2015_H1N1_ref.gb")) 
test2 <- genbankr::readGenBank(file.path(Reference_file_path,"H1N1_ref_NS.gb"))
test[[8]] <- test2
H1N1_ref_genbank <- test


Strain = "H3N2"
genbank_ref <- H3N2_ref_genbank

Create_dandelion(Muts_VPC, "H3N2", H3N2_ref_genbank, gene_segment_names, gene_features_palette)

Create_dandelion(Muts_VPC, "H1N1", H1N1_ref_genbank, gene_segment_names, gene_features_palette)

# Extract the HA non-synonymous mutations and save separately for H1N1 and H3N2
H3N2_HA_Muts <- Muts_VPC[Muts_VPC$Strain == "H3N2" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 4,]
file_name <- paste0(Analysis_date,"-H3N2_HA_Muts.csv")
write.csv(H3N2_HA_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H3N2_NA_Muts <- Muts_VPC[Muts_VPC$Strain == "H3N2" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 6,]
file_name <- paste0(Analysis_date,"-H3N2_NA_Muts.csv")
write.csv(H3N2_NA_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H1N1_HA_Muts <- Muts_VPC[Muts_VPC$Strain == "H1N1" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 4,]
file_name <- paste0(Analysis_date,"-H1N1_HA_Muts.csv")
write.csv(H1N1_HA_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H1N1_NA_Muts <- Muts_VPC[Muts_VPC$Strain == "H1N1" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 6,]
file_name <- paste0(Analysis_date,"-H1N1_NA_Muts.csv")
write.csv(H1N1_NA_Muts,file.path(SegmentSpecific_Muts_path,file_name))


# Extract Mutations for potential T Cell Epitopes
H3N2_NP_Muts <- Muts_VPC[Muts_VPC$Strain == "H3N2" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 5,]
file_name <- paste0(Analysis_date,"-H3N2_NP_Muts.csv")
write.csv(H3N2_NP_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H3N2_PA_Muts <- Muts_VPC[Muts_VPC$Strain == "H3N2" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 3,]
file_name <- paste0(Analysis_date,"-H3N2_PA_Muts.csv")
write.csv(H3N2_PA_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H3N2_PB1_Muts <- Muts_VPC[Muts_VPC$Strain == "H3N2" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 2,]
file_name <- paste0(Analysis_date,"-H3N2_PB1_Muts.csv")
write.csv(H3N2_PB1_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H1N1_NP_Muts <- Muts_VPC[Muts_VPC$Strain == "H1N1" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 5,]
file_name <- paste0(Analysis_date,"-H1N1_NP_Muts.csv")
write.csv(H1N1_NP_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H1N1_PA_Muts <- Muts_VPC[Muts_VPC$Strain == "H1N1" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 3,]
file_name <- paste0(Analysis_date,"-H1N1_PA_Muts.csv")
write.csv(H1N1_PA_Muts,file.path(SegmentSpecific_Muts_path,file_name))

H1N1_PB1_Muts <- Muts_VPC[Muts_VPC$Strain == "H1N1" & Muts_VPC$effect_any == "Nonsynonymous" & Muts_VPC$segment_num == 2,]
file_name <- paste0(Analysis_date,"-H1N1_PB1_Muts.csv")
write.csv(H1N1_PB1_Muts,file.path(SegmentSpecific_Muts_path,file_name))

# Calculate probability of iSNV distribution within HA gene segment
nMut_H3_nonSyn <- 31
HA_mat_length <- 550
Antigenic_site_residues <- 62
Antigenic_site_muts<- 13 
NonAntigenic_site_muts <- 31 - 13
Expected_antigenic_site_muts <- nMut_H3_nonSyn/HA_mat_length * Antigenic_site_residues
Expected_antigenic_site_muts
Expected_NonAntigenic_site_muts <- nMut_H3_nonSyn/HA_mat_length * (HA_mat_length - Antigenic_site_muts)
Expected_NonAntigenic_site_muts


# Expected probability of a mutation occurring at an antigenic site by chance
p_expected <- Antigenic_site_residues / HA_mat_length

# Perform the binomial test
result <- binom.test(Antigenic_site_muts, nMut_H3_nonSyn, p_expected, alternative = "greater")

# Show results
print(result)
