

# Prepare environment -----------------------------------------------------

rm(list = ls())
R_path <- getwd()

# Set options for run
RunDate <-
  c(33122, 50422, 61722, 91222, 22323, 31023, 50123) # These are the run dates of the sample tracker to include
segment_name <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
H3N2_Ref <- read.csv("Reference_files/H3N2_Ref_A_2016.csv")
H1N1_Ref <- read.csv("Reference_files/H1N1_Ref_A_2015.csv")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
dir.create(Compiled_OutPut_path)
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
dir.create(Figure_OutPut_path)

# Load Packages -----------------------------------------------------
# Load required packages
# 'stringr' and 'stringi' for string manipulation, 
# 'phylotools' for reading FASTA files, 'dplyr' for data manipulation,
# 'gridExtra' for arranging grid graphical objects, 
# 'seqinr' for biological sequences retrieval and analysis,
# 'GenomicRanges' for representing and manipulating genomic intervals and variables,
# 'ShortRead' for input, quality assessment, manipulation and output of high-throughput sequencing data,
# 'Biostrings' for string objects representing biological sequences, 
# 'ggplot2' for creating elegant data visualisations,
# 'data.table' for fast aggregation of large data

library(stringr) 
library(stringi) 
library(phylotools) 
library(dplyr)
library(gridExtra) 
library(seqinr)
library(GenomicRanges)
library(ShortRead)
library(Biostrings)
library(ggplot2) 
library(data.table)


# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("stringr","stringi", "phylotools", "dplyr","gridExtra", "seqinr","ShortRead","Biostrings","ggplot2","data.table")

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
writeLines(all_citations_text, con = file.path(Compiled_OutPut_path,"Compile_data_citations.bib"))



# Define Functions -----------------------------------------------------

writeFasta <- function(data, filename) {
  fastaLines = c()
  for (rowNum in 1:nrow(data)) {
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

GetStringDiff <-
  function(a,
           b,
           exclude = c("-", "?"),
           ignore.case = TRUE,
           show.excluded = FALSE) {
    if (nchar(a) != nchar(b))
      stop("Lengths of input strings differ. Please check your input.")
    if (ignore.case)
    {
      a <- toupper(a)
      b <- toupper(b)
    }
    split_seqs <- strsplit(c(a, b), split = "")
    only.diff <- (split_seqs[[1]] != split_seqs[[2]])
    only.diff[(split_seqs[[1]] %in% exclude) |
                (split_seqs[[2]] %in% exclude)] <- NA
    diff.info <- data.frame(which(is.na(only.diff) | only.diff),
                            split_seqs[[1]][only.diff], split_seqs[[2]][only.diff])
    names(diff.info) <- c("position", "poly.seq.a", "poly.seq.b")
    if (!show.excluded)
      diff.info <- na.omit(diff.info)
    diff.info
  }


# Compile SampleTrackers -----------------------------------------------------

SampleTracker_Compiled <- data.frame()
sd <- 1
for (sd in 1:length(RunDate)) {
  # Generate name of sample tracker
  tracker_name <-
    paste("PostProcessing_Output/",RunDate[sd], "_SampleTracker_CHOA.csv", sep = "")
  tracker_tmp <- read.csv(tracker_name)
  tracker_tmp <-
    tracker_tmp[-c(1)] # Get rid of the first index column
  # Add sample tracker to compiled tracker list
  SampleTracker_Compiled <-
    rbind(SampleTracker_Compiled, tracker_tmp)
  rownames(SampleTracker_Compiled) <- NULL
}
colnames(SampleTracker_Compiled) <-
  c("Sample",
    "Strain",
    "SeqDate",
    "PB2_accession",
    "propN",
    "QC_Cov")

if (length(unique(SampleTracker_Compiled$Sample)) != nrow(SampleTracker_Compiled)) {
  stop("There are duplicate samples IDs that need to be fixed before we proceed.")
}

# Save the compiled sample trackers
SampleTracker_Compiled_path <- file.path(Compiled_OutPut_path,"SampleTracker_Compiled_ALL.csv")
write.csv(SampleTracker_Compiled, SampleTracker_Compiled_path)

# Here we filter for samples passing QC for coverage
SampleTracker_Compiled <-
  SampleTracker_Compiled[SampleTracker_Compiled$propN < 5, ]
passingSamples <- SampleTracker_Compiled$Sample
SampleTracker_Compiled_path <- file.path(Compiled_OutPut_path,"SampleTracker_Compiled.csv")
write.csv(SampleTracker_Compiled,SampleTracker_Compiled_path)


# Compile Mutation Output -----------------------------------------------------

MutationOutput_Compiled <- data.frame()
md <- 1
# Compile Mutations
for (md in 1:length(RunDate)) {
  mutation_name <-
    paste0("PostProcessing_Output/MutationOutput_", RunDate[md], ".csv")
  mutation_tmp <- read.csv(mutation_name)
  MutationOutput_Compiled <-
    rbind(MutationOutput_Compiled, mutation_tmp)
}
MutationOutput_Compiled$X <- NULL

# Filter mutations to include only those samples with low propN's
MutationOutput_Compiled <-
  MutationOutput_Compiled[MutationOutput_Compiled$subject %in% passingSamples, ]
unique(MutationOutput_Compiled$subject)

# Create Fasta Copies -----------------------------------------------------

# Create a directory for the masked consensus fasta files
dir.create("PostProcessing_OutPut/MaskedConsensus_CHOP")
rd <- 1
for (rd in 1:length(RunDate)) {
  run <- RunDate[rd]
  # Set the source and destination directories
  src_dir <-
    paste0(R_path, "/Pipeline_Output_Files/",RunDate[rd],"_PipelineOutput/MaskedConsensus")
  dest_dir <-
    paste0(R_path, "/PostProcessing_OutPut/MaskedConsensus_CHOP/", run, "_MaskedConsensus")

  # Create the destination directory
  dir.create(dest_dir)
  
  # List all files in the source directory
  files <- list.files(src_dir)
  # Copy files from source to destination directory
  file.copy(file.path(src_dir, files),
            file.path(dest_dir, files),
            overwrite = TRUE)
}


# Fix Mutations -----------------------------------------------------

MutationOutput_Compiled$Modified <- 0
rownames(MutationOutput_Compiled) <- NULL
to_delete <- data.frame()
to_add <- data.frame()
# Working script to fix reference/variant consensus issues

target_muts <-
  MutationOutput_Compiled[MutationOutput_Compiled$var_freq > 0.5 &
                            MutationOutput_Compiled$variant_level == "Intra", ]
nrow(target_muts)
# Loop through the list of mutations and update the mutation and consensus sequence
MutationOutput_Unmodified <- MutationOutput_Compiled

MutationOutput_Unmodified_path <- paste0(Compiled_OutPut_path,"/MutationOutput_Unmodified.csv")
write.csv(MutationOutput_Unmodified, MutationOutput_Unmodified_path)
mut <- 1
for (mut in 1:nrow(target_muts)) {
  print(nrow(target_muts) - mut)
  mut_tmp <- target_muts[mut, ]
  row_number <- as.numeric(rownames(mut_tmp))
  # Extract the current segment_num, position, and var_nt
  segment_num <- as.integer(mut_tmp$segment_num)
  position <- as.integer(mut_tmp$position)
  var_nt <- mut_tmp$var_nt
  ref_nt <- mut_tmp$ref_nt
  
  sample_tmp <-
    SampleTracker_Compiled[SampleTracker_Compiled$Sample %in% mut_tmp$subject, ]
  
  # Load corresponding reference
  if (sample_tmp$Strain == "H3N2") {
    Ref <- H3N2_Ref
  } else if (sample_tmp$Strain == "H1N1") {
    Ref <- H1N1_Ref
  } else {
    stop("Non-existing strain!!!")
  }
  
  # Load the fasta file for the sample
  masked_path <-
    paste0(R_path,
           "/PostProcessing_Output/MaskedConsensus_CHOP/",
           sample_tmp$SeqDate,
           "_MaskedConsensus")
  masked_path
  target_fasta <-
    intersect(
      list.files(path = masked_path, pattern = mut_tmp$subject),
      list.files(masked_path, pattern = ".fasta")
    )
  tmp_fasta <-
    phylotools::read.fasta(file.path(masked_path, target_fasta))
  
  # Get the sequence corresponding to the segment_num
  seq <- tmp_fasta$seq.text[segment_num]
  
  
  # Determine the current nucleotide at the position
  target_nucleotide <-
    substring(tmp_fasta$seq.text[segment_num], position, position)
  target_nucleotide
  
  # Make sure the nucleotide isn't "N", if it is the the nucleotide is masked anyways
  if (target_nucleotide != "N") {
    # Make sure the target_nucleotide matches the current ref_NT
    if (target_nucleotide == ref_nt) {
      # Replace the nucleotide at the specified position
      seq_modified <- stringr::str_sub(seq, 1, position - 1) %>%
        stringr::str_c(var_nt) %>%
        stringr::str_c(stringr::str_sub(seq, position + 1))
      
      # Update the sequence in tmp_fasta
      tmp_fasta$seq.text[segment_num] <- seq_modified
      
      # Save the updated fasta file
      updated_fasta <-
        dplyr::data_frame(name = tmp_fasta$seq.name, seq = tmp_fasta$seq.text) #This is weird formatting, but it's the only way it works
      writeFasta(data = updated_fasta, file.path(masked_path, target_fasta))
      rm(updated_fasta, seq_modified)
    } else{
      stop("The reference nt's don't match!")
    }
  }
  
  # Check if the updated mutation meets the criteria for a SNV
  # Check consensus reference nucleotide
  ref_seq <- Ref$NT[segment_num]
  ref_seq <- gsub("U", "T", ref_seq)
  consensus_ref_nt <-
    substring(ref_seq, position, position)
  consensus_ref_nt
  
  # See if there's an existing population level mutation for the current mutation
  pop_tmp <-
    MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                              MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                              MutationOutput_Compiled$position %in% mut_tmp$position &
                              MutationOutput_Compiled$variant_level == "Pop", ]
  pop_tmp
  if (nrow(pop_tmp) == 0) {
    # There is no population level mutation for this currently, so we reclassify the current mutation
    if (consensus_ref_nt != mut_tmp$var_nt) {
      # This is a consensus level mutation, we re-classify it as such
      # Update the mutation entry
      mut_updated <- mut_tmp
      mut_updated$variant_level <- "Pop"
      mut_updated$Modified <- 1
      MutationOutput_Compiled[row_number, ] <- mut_updated
      rm(mut_updated)
    } else if (consensus_ref_nt == mut_tmp$var_nt) {
      # This is an incorrect classification of reference/variant. We've already fixed the fasta file, so we can just delete this mutation - we'll do this at the end though so the numbering doesn't get messed up
      row_number <- as.numeric(rownames(mut_tmp))
      to_delete <- rbind(to_delete, row_number)
    }
  }
  
  # We must also see if this mutation qualifies as an intrahost mutation that doesn't already exist
  swap_vf <- 1 - mut_tmp$var_freq
  if (swap_vf >= 0.01) {
    intra_tmp <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                                MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                                MutationOutput_Compiled$position %in% mut_tmp$position &
                                MutationOutput_Compiled$var_nt == mut_tmp$var_nt &
                                MutationOutput_Compiled$variant_level == "Intra", ]
    depth_check <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                                MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                                MutationOutput_Compiled$position %in% mut_tmp$position &
                                MutationOutput_Compiled$variant_level == "Intra", ]
    if (nrow(intra_tmp) >= 1 && nrow(pop_tmp) == 1) {
      # There was a population level mutation already, but we still need to make this an in intrahost mutation. We can sub it out here since we don't need to add a mutation
      update_tmp <- mut_tmp
      update_tmp$ref_nt <- mut_tmp$var_nt
      update_tmp$var_nt <- mut_tmp$ref_nt
      update_tmp$var_freq <- 1 - mut_tmp$var_freq
      update_tmp$aa_ref <- mut_tmp$aa_sub
      update_tmp$aa_sub <- mut_tmp$aa_ref
      update_tmp$aa_ref2 <- mut_tmp$aa_sub2
      update_tmp$aa_sub2 <- mut_tmp$aa_ref2
      update_tmp$variant_level <- "Intra"
      if (nrow(depth_check) > 1) {
        update_tmp$var_depth <- NaN
      } else {
        update_tmp$var_depth <- mut_tmp$depth - mut_tmp$var_depth
      }
      update_tmp$Modified <- 1
      MutationOutput_Compiled[row_number, ] <- update_tmp
      rm(update_tmp)
    }
    else if (nrow(intra_tmp) == 1 && nrow(pop_tmp) == 0) {
      stop("check me out!")
    } else if (nrow(intra_tmp) == 0) {
      # Side note - this will be negative because we've already changed the initial mutation from intra --> pop
      #This qualifies as an intrahost variant as well, so we need to generate one
      update_tmp <- mut_tmp
      update_tmp$ref_nt <- mut_tmp$var_nt
      update_tmp$var_nt <- mut_tmp$ref_nt
      update_tmp$var_freq <- 1 - mut_tmp$var_freq
      update_tmp$aa_ref <- mut_tmp$aa_sub
      update_tmp$aa_sub <- mut_tmp$aa_ref
      update_tmp$aa_ref2 <- mut_tmp$aa_sub2
      update_tmp$aa_sub2 <- mut_tmp$aa_ref2
      update_tmp$variant_level <- "Intra"
      if (nrow(depth_check) > 1) {
        update_tmp$var_depth <- NaN
      } else {
        update_tmp$var_depth <- mut_tmp$depth - mut_tmp$var_depth
      }
      update_tmp$Modified <- 1
      update_tmp
      to_add <- rbind(to_add, update_tmp)
      rm(update_tmp)
    } else {
      stop("Check me out!")
    }
  }
  rm(mut_tmp)
}

if (nrow(to_delete) > 0) {
  # Remove to_delete rows
  colnames(to_delete) <- "Row"
  
  # Sort to_delete$Row in descending order
  to_delete_sorted <- sort(to_delete$Row, decreasing = TRUE)
  
  # Remove the rows specified in to_delete$Row from MutationOutput_Compiled
  for (row in to_delete_sorted) {
    MutationOutput_Compiled <- MutationOutput_Compiled[-row,]
  }
}


# Add the to_add rows
MutationOutput_Compiled <- rbind(MutationOutput_Compiled, to_add)
MutationOutput_Compiled <-
  MutationOutput_Compiled[order(MutationOutput_Compiled$position), ]
MutationOutput_Compiled <-
  
  MutationOutput_Compiled[order(MutationOutput_Compiled$segment_num), ]
MutationOutput_Compiled <-
  MutationOutput_Compiled[order(MutationOutput_Compiled$subject), ]
Mutation_Output_Compiled_path <- paste0(Compiled_OutPut_path,"/MutationOutput_Compiled.csv")
write.csv(MutationOutput_Compiled, Mutation_Output_Compiled_path)

# Generate Consensus Sequences -----------------------------------------------------
# Split the sample names into ID (CHOA-XXX), replicate number (D1 ...), and version number (V...)
Sample_ID <- stringr::str_sub(SampleTracker_Compiled$Sample, 1, 8)
uID <- unique(Sample_ID)

get_majority_consensus <- function(seq_text) {
  # Create a matrix with the counts of each nucleotide for each position
  counts <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(seq_text), baseOnly = TRUE)
  # Extract the majority nucleotides using the max count for each column
  majority_nucleotides <- apply(counts, 2, function(x) names(which.max(x)))
  majority_nucleotides <- gsub("other", "N", majority_nucleotides)
  # Combine the majority nucleotides into a single consensus sequence
  consensus_seq <- paste(majority_nucleotides, collapse = "")
  return(consensus_seq)
}

dir.create(file.path(R_path,"PostProcessing_OutPut/CompiledOutPut/Consensus_seqs"))
Generate_consensus_sequences <- function(uID, SampleTracker_Compiled, R_path, segment_name) {
  mismatches <-
    data.frame(
      Sample = character(),
      RunDate = character(),
      SegNum = double(),
      Position = double(),
      NT_1 = character(),
      NT_2 = character()
    )
  PassingSample_Matches <- data.frame()
  
  for (ns in 1:length(uID)) {
    IDtmp <- uID[ns]
    # Pull all passing versions of this sample
    passing_samples <-
      SampleTracker_Compiled[grep(IDtmp, SampleTracker_Compiled$Sample), ]
    CHOA_PB2 <- data.frame(seq.name = character(), seq.text = character())
    CHOA_PB1 <- data.frame(seq.name = character(), seq.text = character())
    CHOA_PA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_HA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NP <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_MP <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NS <- data.frame(seq.name = character(), seq.text = character())
    
    # Now we determine the consensus sequence for each sample by comparing all masked fasta sequences and choosing the majority nucleotide
    ps <- 1
    for (ps in 1:nrow(passing_samples)) {
      # Get name of masked consensus sequence
      passing_tmp <- passing_samples[ps, ]
      masked_path <-
        paste0(R_path,
               "/PostProcessing_Output/MaskedConsensus_CHOP/",
               passing_tmp$SeqDate,
               "_MaskedConsensus")
      target_fasta <-
        list.files(path = masked_path, pattern = passing_tmp$Sample)
      tmp_fasta <-
        phylotools::read.fasta(file.path(masked_path, target_fasta))
      tmp_fasta$seq.name
      CHOA_PB2 <- rbind(CHOA_PB2, tmp_fasta[1, ])
      CHOA_PB1 <- rbind(CHOA_PB1, tmp_fasta[2, ])
      CHOA_PA <- rbind(CHOA_PA, tmp_fasta[3, ])
      CHOA_HA <- rbind(CHOA_HA, tmp_fasta[4, ])
      CHOA_NP <- rbind(CHOA_NP, tmp_fasta[5, ])
      CHOA_NA <- rbind(CHOA_NA, tmp_fasta[6, ])
      CHOA_MP <- rbind(CHOA_MP, tmp_fasta[7, ])
      CHOA_NS <- rbind(CHOA_NS, tmp_fasta[8, ])
    }
    
    PB2_consensus <-
      get_majority_consensus(CHOA_PB2$seq.text)
    PB1_consensus <-
      get_majority_consensus(CHOA_PB1$seq.text)
    PA_consensus <-
      get_majority_consensus(CHOA_PA$seq.text)
    HA_consensus <-
      get_majority_consensus(CHOA_HA$seq.text)
    NP_consensus <-
      get_majority_consensus(CHOA_NP$seq.text)
    NA_consensus <-
      get_majority_consensus(CHOA_NA$seq.text)
    MP_consensus <-
      get_majority_consensus(CHOA_MP$seq.text)
    NS_consensus <-
      get_majority_consensus(CHOA_NS$seq.text)
    
    new_name <- paste0(IDtmp, "_consensus.fasta")
    new_frame <- data.frame(seq.name = character(), seq.text = character())
    seq_names <-
      c(
        paste0(IDtmp, "_PB2"),
        paste0(IDtmp, "_PB1"),
        paste0(IDtmp, "_PA"),
        paste0(IDtmp, "_HA"),
        paste0(IDtmp, "_NP"),
        paste0(IDtmp, "_NA"),
        paste0(IDtmp, "_MP"),
        paste0(IDtmp, "_NS")
      )
    seq_names
    seqs <-
      rbind(
        PB2_consensus,
        PB1_consensus,
        PA_consensus,
        HA_consensus,
        NP_consensus,
        NA_consensus,
        MP_consensus,
        NS_consensus
      )
    rownames(seqs) <- NULL
    new_frame <- data.frame(seq.name <- seq_names, seq.text <- seqs)
    colnames(new_frame) <- c("seq.name", "seq.text")
    
    new_fasta_df <-
      dplyr::data_frame(name = new_frame$seq.name, seq = new_frame$seq.text) #This is weird formatting, but it's the only way it works
    consensus_seq_path <- file.path(R_path,"PostProcessing_OutPut/CompiledOutPut/Consensus_seqs",new_name)
    consensus_seq_path
    writeFasta(data = new_fasta_df, filename = consensus_seq_path) # This may actually be seqinr, if there's an error, repeat
    segment_match <-
      matrix(data = NaN,
             nrow = nrow(passing_samples),
             ncol = 10)
    colnames(segment_match) <- c(segment_name, 'All', "nMismatch")
    segment_match <- as.data.frame(segment_match)
    segment_match$nMismatch <- as.numeric(segment_match$nMismatch)
    pps <- 1
    for (pps in 1:nrow(passing_samples)) {
      # Get name of masked consensus sequence
      passing_tmp <- passing_samples[pps, ]
      passing_tmp
      Sample <- passing_tmp$Sample
      masked_path <-
        paste0(R_path,
               "/PostProcessing_Output/MaskedConsensus_CHOP/",
               passing_tmp$SeqDate,
               "_MaskedConsensus")
      target_fasta <-
        list.files(path = masked_path, pattern = passing_tmp$Sample)
      tmp_fasta <-
        phylotools::read.fasta(file.path(masked_path, target_fasta))
      for (seg in 1:8) {
        if (nchar(tmp_fasta$seq.text[seg]) == nchar(new_fasta_df$seq[seg])) {
          setDiff <-
            GetStringDiff(a = tmp_fasta$seq.text[seg],
                          b = new_fasta_df$seq[seg],
                          exclude = "N")
          if (nrow(setDiff) == 0) {
            segment_match[pps, seg] <- "Y"
          } else {
            segment_match[pps, seg] <- "N"
            Position <- setDiff[1]
            
            NT_1 <- setDiff[2]
            NT_2 <- setDiff[3]
            # Add new mismatches to tabke of mismatches
            mismatches_this <-
              cbind(Sample, seg, Position, NT_1, NT_2)
            mismatches_this
            mismatches <- rbind(mismatches, mismatches_this)
            
          }
        } else{
          segment_match[pps, seg] <- "N"
        }
      }
      
      mismatches_this <- NULL
      if ('N' %in% segment_match[pps, ]) {
        segment_match[pps, 9] = "N"
      } else {
        segment_match[pps, 9] = 'Y'
      }
      nMismatch <- sum(mismatches$Sample == Sample)
      segment_match[pps, 10] <- nMismatch
      nMismatch <- NULL
    }
    PassingSample_Matches_tmp <-
      cbind(passing_samples$Sample, segment_match)
    PassingSample_Matches <-
      rbind(PassingSample_Matches, PassingSample_Matches_tmp)
    PassingSample_Matches_tmp
    
  }
  return(PassingSample_Matches)
}

PassingSample_Matches <- Generate_consensus_sequences(uID, SampleTracker_Compiled, R_path, segment_name)

# See if there are any significant mismatches
Mismatched_Samples <-
  PassingSample_Matches[PassingSample_Matches$nMismatch > 100, ]
Mismatched_Samples

# Save mismatched samples (may be able to salvage later if we can find their pair)
Mismatched_Samples_path <- paste0(Compiled_OutPut_path,"/CHOA_Mismatched_Samples.csv")

write.csv(Mismatched_Samples,Mismatched_Samples_path)


# Lots of mismatches can screw with consensus generation. We remove these samples from the SampleTracker_Compiled and MutationOutput_Compiled then repeat the analysis
SampleTracker_Compiled <-
  SampleTracker_Compiled[!SampleTracker_Compiled$Sample %in% Mismatched_Samples$`passing_samples$Sample`, ]

MutationOutput_Compiled <-
  MutationOutput_Compiled[!MutationOutput_Compiled$subject %in% Mismatched_Samples$`passing_samples$Sample`, ]


# Repeat consensus generation --------------------------------------------------

# Split the sample names into ID (CHOA-XXX), replicate number (D1 ...), and version number (V...)
Sample_ID_revised <- stringr::str_sub(SampleTracker_Compiled$Sample, 1, 8)
uID_revised <- unique(Sample_ID_revised)

# We repeat the consensus sequence generation, excludingly highly mismatched samples
PassingSample_Matches_Revised <- Generate_consensus_sequences(uID_revised, SampleTracker_Compiled, R_path, segment_name)

# Sort Passing_Compiled and PassingSample_Matches dataframe by sample ID
Passing_CompiledSamples <-
  SampleTracker_Compiled[order(SampleTracker_Compiled$Sample), ]
PassingSample_Matches_Revised <-
  PassingSample_Matches_Revised[order(PassingSample_Matches_Revised$`passing_samples$Sample`), ]
Passing_CompiledSamples
PassingSample_Matches_Revised

# Split the sample ID into components
SampleID_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 1, 8)
SampleID_only
Replicates_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 10, 11)
Replicates_only
Versions_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 12, 13)
Versions_only


# Remove row indices from the PassingSample_Matches ad Passing_CompiledSamples
rownames(Passing_CompiledSamples) <-
  NULL
head(Passing_CompiledSamples)
rownames(PassingSample_Matches) <- NULL
head(PassingSample_Matches)
nrow(Passing_CompiledSamples)
nrow(PassingSample_Matches)

PassingSample_QC <-
  data.frame(
    Sample = Passing_CompiledSamples$Sample,
    ID = SampleID_only,
    Replicate = Replicates_only,
    Version = Versions_only,
    Strain = Passing_CompiledSamples$Strain,
    SeqDate = Passing_CompiledSamples$SeqDate,
    PB2_accession = Passing_CompiledSamples$PB2_accession,
    propN = Passing_CompiledSamples$propN,
    QC_Cov = Passing_CompiledSamples$QC_Cov
  )
head(PassingSample_QC)

if (!all(PassingSample_Matches_Revised$`passing_samples$Sample` == PassingSample_QC$Sample)) {
  stop("The rows are not equivalent")
}

tmp_combo <- cbind(PassingSample_QC, PassingSample_Matches_Revised)
tmp_combo

# All samples match, so now we drop the V1 (repeated name) column
CHOA_SampleQC <- tmp_combo
head(CHOA_SampleQC)
dim(CHOA_SampleQC)

write.csv(CHOA_SampleQC, paste0(Compiled_OutPut_path,"/CHOA_SampleQC.csv"))

# Extract Samples that pass all QC metrics unambiguously
passing_CHOASamples <- CHOA_SampleQC[CHOA_SampleQC$nMismatch == 0, ]
head(passing_CHOASamples)
uID_pass <- unique(passing_CHOASamples$ID)
uID_pass

# We also save those samples with mismatches as they may still meet quality control standards with manual inspection
mismatch_CHOASamples <- CHOA_SampleQC[CHOA_SampleQC$nMismatch > 0, ]
write.csv(mismatch_CHOASamples, paste0(Compiled_OutPut_path,"/mismatch_CHOASamples.csv"))
mismatch_CHOASamples <- read.csv(paste0(Compiled_OutPut_path,"/mismatch_CHOASamples.csv"))

View(mismatch_CHOASamples)

# Select Replicates -----------------------------------------------------
sample_mat <- data.frame()
selected_samples <- data.frame()
not_replicate_samples <- data.frame()
rs <- 1
for (rs in 1:length(uID_pass)) {
  uID_tmp <- uID_pass[rs]
  uID_tmp
  uPassing <-
    passing_CHOASamples[passing_CHOASamples$ID == uID_tmp, ]
  uPassing <-
    uPassing[order(uPassing$propN, decreasing = FALSE), ]
  uPassing
  # Check for duplicates
  if (any(duplicated(uPassing$Sample))) {
    stop("Error: Duplicated samples found in 'Sample' column of 'uPassing'")
  }
  
  distinctPassing <-
    uPassing %>% distinct(uPassing$ID, uPassing$Replicate, .keep_all = TRUE) # Get the rows with the unique replicates
  distinctPassing
  distinctPassing <-
    distinctPassing[-c(21, 22)] # Get rid of the two sorting columns
  distinctPassing
  
  # Select samples for further mutational analysis. We prioritize samples that
  # are techinical replicates and use the two samples with the lowest propN. If
  # 2 technical replicates are not available, we expand selection to different
  # versions
  if (nrow(uPassing) == 1) {
    # There are no replicates
    next
  }
  if (nrow(distinctPassing) >= 2) {
    selected_samples_tmp <- distinctPassing[1:2, ]
    selected_samples_tmp
    print(selected_samples_tmp)
    selected_samples <- rbind(selected_samples, selected_samples_tmp)
  } else {
    uPassing <-
      uPassing[order(uPassing$propN, decreasing = FALSE), ]
    uPassing
    not_replicate_samples_tmp <-
      uPassing[1:2, ]
    not_replicate_samples_tmp
    not_replicate_samples <-
      rbind(not_replicate_samples, not_replicate_samples_tmp)
  }
}
selected_samples
CHOA_Selected_Samples <- selected_samples
CHOA_Selected_Samples
head(CHOA_Selected_Samples)
write.csv(CHOA_Selected_Samples, paste0(Compiled_OutPut_path,"/CHOA_Selected_Samples.csv"))

CHOA_NotReplicated_Samples <- not_replicate_samples
not_replicate_samples
write.csv(CHOA_NotReplicated_Samples,paste0(Compiled_OutPut_path,"/CHOA_NotReplicated_Samples.csv"))


# Identify Samples with mismatches that can still be considered ----------

# Determine which samples have not been included in selected sample sets
NotSelected_Samples <-
  dplyr::setdiff(CHOA_SampleQC$ID, CHOA_Selected_Samples$ID)

# A total of 44 samples were sequenced but not included in the final sample set
# due to the absence of a paired sample replicate. This does not include samples
# with propN > 5%.

# Extract unmatched samples
Unmatched_Samples <-
  CHOA_SampleQC[CHOA_SampleQC$ID %in% NotSelected_Samples, ]

write.csv(Unmatched_Samples, paste0(Compiled_OutPut_path,"/CHOA_UnMatched_Samples.csv"))

# We inspect the unmatched samples manually to see if the samples may still be
# included. The reason for automated exclusion is either: 1) no sequence
# replicate (which can't be fixed) or the presence of mismatches between the
# consensus sequences. We can consider consensus mismatches if they occur within
# 20bp of the segment ends, represent a mismatch between a called and masked
# nucleotide, are in a region of poor sequence quality, or occur positions with
# high variation. The final set of selected samples can be found in the CSV
# file: CHOA_Selected_Samples_revised.csv. See below for a summary of the
# samples we have chosen to include following manual inspection and the
# reasoning.


# CHOA-051_D1V1 and CHOA-051_D3V1. There was a single mismatch between the
# consensus sequences for these samples at NS 306. Inspection of the
# CHOA-051_D3V1 isolate showed a high variant (49%) iSNV at NS 306 with the
# variant allele matching the consensus for D1V1, so we will use these samples.

# CHOA-128_D1V1 and CHOA-128_D2V1. These can be kept as there was a single
# mimatch at PB2 813. Both samples had high frequency variants at that site.


# Fix inter/intrahost variants --------------------------------------------------------

# Read in the csv file that contains selected samples
CHOA_Selected_Samples <-
  read.csv(paste0(Compiled_OutPut_path,"/CHOA_Selected_Samples_Revised.csv"))

# Create a list of samples from the data frame
sample_list <- CHOA_Selected_Samples$Sample

# Initialize an empty data frame to store the mutation output
CHOA_MutationOutput <- data.frame()

# Loop over the sample list
for (nv in 1:length(sample_list)) {
  # Print the number of remaining iterations to the console
  length(sample_list) - nv
  
  # Get the current sample from the list
  sample_tmp <- sample_list[nv]
  
  # If the current sample is not a NaN value
  if (is.nan(sample_tmp) == FALSE) {
    # Extract the mutations present in that sample from the compiled mutation output data frame
    sample_muts <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject == sample_tmp, ]
    
    # Add the extracted mutations to the mutation output data frame for selected samples
    CHOA_MutationOutput <- rbind(CHOA_MutationOutput, sample_muts)
    
    # Remove the temporary data frame to free up memory
    rm(sample_muts)
  }
}


# Extract the subjects from the mutation output data frame and store them in a new variable
subject_mut <- CHOA_MutationOutput$subject

# Extract various parts of the subject field for further analysis

# Extract the sample ID (first 8 characters of the subject field)
Sample_ID_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 1, 8)

# Extract the replicate number (characters 10 and 11 of the subject field)
Replicates_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 10, 11)

# Extract the version number (characters 12 and 13 of the subject field)
Versions_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 12, 13)

# Combine the extracted fields with the original mutation output data frame
CHOA_mutations <-
  cbind(Sample_ID_mut,
        Replicates_mut,
        Versions_mut,
        CHOA_MutationOutput)

# Fix variants that have NA as a frequency. These were hard to calculate for positions with multiple alleles present

# Remove the row names from the mutation data frame
rownames(CHOA_mutations) <- NULL

# Find the rows with NA in the var_depth column
to_fix <- CHOA_mutations[is.na(CHOA_mutations$var_depth), ]

# Get the row indices of the rows to fix
indices <- as.numeric(rownames(to_fix))

# Calculate var_depth by multiplying var_freq by depth and round it off
to_fix$var_depth <- round(to_fix$var_freq * to_fix$depth)

# Print the data frame with fixed var_depth values to the console
to_fix

if (isEmpty(to_fix) != TRUE){
  # Initialize a loop counter
  tf <- 1
  
  # Loop over the rows to fix
  for (tf in 1:nrow(to_fix)) {
    # Get the target index from the indices list
    ind_target <- indices[tf]
    
    
    # Replace the row in the original data frame with the fixed row
    CHOA_mutations[ind_target, ] <- to_fix[tf, ]
  }
}


# Select only the rows with 'Intra' in the variant_level column
CHOA_intrahost <-
  CHOA_mutations[CHOA_mutations$variant_level == "Intra", ]

# Print the number of rows in the intra-host data frame to the console
nrow(CHOA_intrahost)

# Remove the row names from the intra-host data frame
rownames(CHOA_intrahost) <- NULL

write.csv(CHOA_intrahost,paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))

# Find samples to exclude ------------------------------------------------------

# Load all mutations
CHOA_intrahost <- read.csv(paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))

# Filter samples into originals vs replicates 
Original_Samples <- CHOA_Selected_Samples %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup()
# Generate a list of samples tagged as "original"
osID <- Original_Samples$Sample

# Filter the second occurrence of each ID and store in Replicate_Samples
Replicate_Samples <- CHOA_Selected_Samples %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup()
# Generate a list of samples tagged as "replicate"
rsID <- Replicate_Samples$Sample


# Create a data frame to store information about each replicate mutation
replicate_mut_info <-
  data.frame(
    ReplicateSubject = character(),
    Presence = character(),
    RepFreq = double(),
    RepDepth = double()
  )


# Filter mutations for min_var = 10% to exclude samples w/poor repeatability,
# this is before we establish variant calling thresholds so we just want to make
# sure that variants called with high frequency in one replicate are present in
# the second replicate. The lowest frequency considered is 2%. 

maxFreq <- 0.03
replicateFreq <- 0.03
min_cov <- 100
Check_ReplicateMut <- function(CHOA_targets, mut_tmp, CHOA_Selected_Samples, pileup_path) {
  # Remove mutations from the same sample
  paired_mutation <-
    CHOA_targets[CHOA_targets$subject != mut_tmp$subject, ]
  
  # Find the subject's replicate
  replicate_subject <-
    CHOA_Selected_Samples[CHOA_Selected_Samples$ID == mut_tmp$Sample_ID_mut &
                            CHOA_Selected_Samples$Sample != mut_tmp$subject, ]$Sample
  
  # Filter the paired mutation data frame to only include rows where the subject is the replicate subject
  paired_mutation <-
    paired_mutation[paired_mutation$subject == replicate_subject, ]
  
  # Continue to filter the paired mutation data frame based on segment number and position
  paired_mutation <-
    paired_mutation[paired_mutation$segment_num == mut_tmp$segment_num, ]
  
  paired_mutation <-
    paired_mutation[paired_mutation$position == mut_tmp$position, ]
  
  # If there is more than one paired mutation
  if (nrow(paired_mutation) > 1) {
    # Filter the paired mutation data frame to only include rows where the variant nucleotide matches the current mutation's variant nucleotide
    paired_mutation <-
      paired_mutation[paired_mutation$var_nt == mut_tmp$var_nt, ]
    
    print(paired_mutation)
  }
  
  # If there is exactly one paired mutation
  if (nrow(paired_mutation) == 1) {
    # mutation is present as replicate, huzzah!, the mutation is present as a replicate
    ReplicateSubject <- replicate_subject
    Presence <- "Y"
    RepFreq <- paired_mutation$var_freq
    RepDepth <- paired_mutation$depth
  } else {
    # Now we see if the replicate is present but not called
    IDtmp <- mut_tmp$Sample_ID_mut
    IDtmp
    test <-
      CHOA_Selected_Samples[CHOA_Selected_Samples$ID == IDtmp, ]
    test
    if (is.na(test$ID[2] == TRUE)) {
      # There is no replicate mutation because there's no replicate
      ReplicateSubject <- "No Rep"
      Presence <- "No Rep"
      RepFreq <- NaN
      RepDepth <- NaN
    } else {
      # There is a replicate but the mutation doesn't seem to be in that replicate. Now we need to assess the pileup file
      Presence <- "N"
      replicate <-
        test[test$Sample != mut_tmp$subject, ]
      replicate # Determine the replicate sample info
      ReplicateSubject <- replicate_subject
      target_file <-
        list.files(path = pileup_path, pattern = replicate$Sample)
      target_file
      pileup_tmp <- read.csv(file.path(pileup_path, target_file))
      pileup_tmp <- data.table::data.table(pileup_tmp)
      Rep_reads <-
        pileup_tmp[SegNum == mut_tmp$segment_num & pos == mut_tmp$position]
      Rep_reads
      RepDepth <- sum(Rep_reads$count)
      RepDepth
      Var_reads <-
        Rep_reads[nucleotide == mut_tmp$var_nt]
      Var_reads
      if (nrow(Var_reads) == 0) {
        RepFreq <- 0
      } else {
        RepFreq <- Var_reads$count / RepDepth
      }
    }
  }
  replicate_mut_info_this <-
    data.frame(ReplicateSubject=ReplicateSubject, Presence=Presence, RepFreq=RepFreq, RepDepth=RepDepth)
  return(replicate_mut_info_this)
}
CHOA_targets <-
  CHOA_intrahost[CHOA_intrahost$var_freq >= maxFreq, ]
pileup_path <-
  paste0(R_path, "/PostProcessing_OutPut/Pileup_renamed/")
Selected_IDs <- unique(CHOA_Selected_Samples$ID)
Exclude_IDs <- c() # We haven't identifed samples to exclude yet

nm <- 1
# Check for replicates 
for (nm in 1:nrow(CHOA_targets)) {
  print(nrow(CHOA_targets) - nm)
  
  # Get the current mutation from the data frame
  mut_tmp <- CHOA_targets[nm, ]
  replicate_mut_info_this <- Check_ReplicateMut(CHOA_targets, mut_tmp, CHOA_Selected_Samples, pileup_path)
  
  # Add the replicate mutation information to the data frame
  replicate_mut_info <- rbind(replicate_mut_info,replicate_mut_info_this)
}

# Update formatting of replicate_mut_info
# Convert RepFreq and RepDepth to numeric
replicate_mut_info$RepFreq <-
  as.numeric(replicate_mut_info$RepFreq)
replicate_mut_info$RepDepth <-
  as.numeric(replicate_mut_info$RepDepth)
# Bind the targets and replicate mutation info into a single data frame
CHOA_HiFreq_Muts <- as.data.table(cbind(CHOA_targets, replicate_mut_info))
head(CHOA_HiFreq_Muts)
CHOA_HiFreq_Muts$X <- NULL # Remove index column 

# Filter High Frequency Mutations
CHOA_HiFreq_Muts <- CHOA_HiFreq_Muts[depth >= min_cov & RepDepth >= min_cov]

Check_Original_Replicates <- function(CHOA_Muts,osID,rsID,min_var_rep){
  # Split the data into original and replicate mutations.
  Original_Muts <- CHOA_Muts[CHOA_Muts$subject %in% osID, ]
  Replicate_Muts <- CHOA_Muts[CHOA_Muts$subject %in% rsID, ]
  # Create a new column in the original mutations dataframe, initialized with "N"
  Original_Muts[, Presence_corrected := "N"]
  # Mark the presence of mutations in the original data where RepFreq is greater than min_var_rep
  Original_Muts$Presence_corrected <- ifelse(Original_Muts$RepFreq >= min_var_rep, "Y", "N")
  return(Original_Muts)
}



Get_Bad_IDs <- function(CHOA_Muts, osID, rsID,min_var_rep,ID_version){
  # Create an empty dataframe to store the results
  SampleReproducibility_data <- data.frame(
    Basic_ID = character(),
    SampleID = character(),
    ReplicateID = character(),
    niSNV = integer(),
    nSharediSNV = integer(),
    pShared = double(),
    stringsAsFactors = FALSE
  )
  
  Original_HiFreq_Muts <- Check_Original_Replicates(CHOA_Muts,osID,rsID,min_var_rep=min_var_rep)
  
  # Loop through the list of original sample IDs (osID)
  ix <- 1
  for (ix in 1:length(osID)) {
    SampleID <- osID[ix]
    ReplicateID <- rsID[ix]
    Basic_ID <- sub("_.*", "", SampleID)
    
    # Filter the mutations for the current osID entry that have "Y" in the "Presence_corrected" column
    mutations <- Original_HiFreq_Muts[subject == SampleID, ]
    
    # Calculate the total number of mutations and the number of shared mutations
    niSNV <- nrow(mutations)
    nSharediSNV <- sum(mutations$Presence_corrected == "Y")
    
    # Calculate the proportion of shared mutations
    pShared <- nSharediSNV / niSNV
    
    # Add the current sample's data to the SampleReproducibility_data dataframe
    SampleReproducibility_data <- rbind(SampleReproducibility_data, 
                                        data.frame(Basic_ID, SampleID, ReplicateID, niSNV, nSharediSNV, pShared))
  }
  # Create a bar graph showing the proportion of unshared mutations for each sample
  graph <-
    ggplot(SampleReproducibility_data, aes(x = Basic_ID, y = 1 - pShared)) +
    geom_col(fill = "darkgrey") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("ID") +
    ylab("p(Not Shared)") +
    ggtitle("Sample Reproducibility Data")
  
  ggsave(filename = paste0(Figure_OutPut_path,"/",ID_version,"reproducibility_plot.pdf"), plot = graph)
  
  # Identify the samples where the proportion of shared mutations is less than 0.5
  bad_OG_samples <- SampleReproducibility_data %>%
    filter(!is.na(pShared) & pShared < 0.5 & niSNV > 1)
  
  # Extract the basic IDs of these samples
  Bad_IDs <- bad_OG_samples$Basic_ID
  
  return(Bad_IDs)
}
Bad_IDs_1 <- Get_Bad_IDs(CHOA_HiFreq_Muts,osID, rsID,min_var_rep=0.03,ID_version = "OG")
Bad_IDs_2 <- Get_Bad_IDs(CHOA_HiFreq_Muts,rsID, osID,min_var_rep=0.03,ID_version = "Rep")

Bad_IDs_1
Bad_IDs_2
# Create a list of all bad sample IDs from both rounds of analysis
Bad_IDs_all <- unique(c(Bad_IDs_1, Bad_IDs_2))
Bad_IDs_all

# Sort the list in ascending order
Bad_IDs_all <- Bad_IDs_all[order(Bad_IDs_all)]
Bad_IDs_all
# Write the list to a csv file
write.csv(Bad_IDs_all, paste0(Compiled_OutPut_path,"/CHOA_Bad_IDs_All.csv"))



# Save the finalized sample set by removing Bad_IDs_All from CHOA_Selected_Samples
CHOA_SampleSet_FINAL <- CHOA_Selected_Samples[!CHOA_Selected_Samples$ID %in% Bad_IDs_all,]
write.csv(CHOA_SampleSet_FINAL, paste0(Compiled_OutPut_path,"/CHOA_SampleSet_FINAL.csv"))


# Determine threshold criteria for variant calling ------------------

# Read in data 
Bad_IDs_all <- read.csv(paste0(Compiled_OutPut_path,"/CHOA_Bad_IDs_All.csv"))
Bad_IDs_all$X <- NULL
CHOA_intrahost <- fread(paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))
CHOA_intrahost$V1 <- NULL

# Remove bad IDS from CHOA_intrahost
Muts <- CHOA_intrahost[!Sample_ID_mut %in% Bad_IDs_all$x]
replicate_mut_info <- data.frame()
nm <- 1

Plot_muts <- function(iMuts, min_cov, min_var, min_var_rep) {
  plot_muts <- iMuts[, .(Sample_ID_mut, var_freq, RepFreq)]
  
  # Fit a linear model
  fit <- lm(RepFreq ~ var_freq, data = plot_muts)
  
  # Compute R-squared value
  r_squared <- summary(fit)$r.squared
  
  # Create scatterplot with color-coded points
  var_comp <- ggplot(plot_muts, aes(x = var_freq, y = RepFreq, color = Sample_ID_mut)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + # Add linear model
    annotate("text", x = min(plot_muts$var_freq), y = max(plot_muts$RepFreq), label = paste("R^2 = ", round(r_squared, 2)), hjust = 0, vjust = 1) + # Add R-squared value
    labs(x = "iSNV frequency", y = "Replicate iSNV frequency", color = "Sample") +
    scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) + 
    ggtitle(paste("Replicate Comparison: min_cov =", min_cov, ", min_var =", min_var, ", min_var_rep =", min_var_rep)) + # Add title
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)) # Center the title
  
  print(var_comp)
  
  # Export plot as PDF
  scatter_name <- paste0(Figure_OutPut_path,"/ReplicateFreq_Comparison_MinVar", min_var * 100, "_F1.pdf")
  ggsave(scatter_name, plot = var_comp, width = 10, height = 8)
  
  return(r_squared)
}

# Check for replicates 
for (nm in 1:nrow(Muts)) {
  print(nrow(Muts) - nm)
  
  # Get the current mutation from the data frame
  mut_tmp <- Muts[nm, ]
  replicate_mut_info_this <- Check_ReplicateMut(Muts, mut_tmp, CHOA_Selected_Samples, pileup_path)
  
  
  # Add the replicate mutation information to the data frame
  replicate_mut_info <- rbind(replicate_mut_info,replicate_mut_info_this)
}

# Convert RepFreq and RepDepth to numeric
replicate_mut_info$RepFreq <-
  as.numeric(replicate_mut_info$RepFreq)
replicate_mut_info$RepDepth <-
  as.numeric(replicate_mut_info$RepDepth)

# Bind the targets and replicate mutation info into a single data frame
CHOA_All_Muts <- as.data.table(cbind(Muts, replicate_mut_info))
head(CHOA_All_Muts)
CHOA_All_Muts$X <- NULL # Remove index column 


# Generate vector of variant calling thresholds
min_var_values <- seq(0.01, 0.1, by = 0.005)

# Create an empty list to hold the results
Replicate_Concordance <- data.frame(min_var_freq = double(),min_rep_freq = double(),min_cov=double(),pShared = double(),r_squared = double())
mvf <- 1
for (mvf in 1:length(min_var_values)){
  min_var <- min_var_values[mvf]
  min_var_replicate <- min_var
  min_var
  # Filter mutations 
  Muts_this <- CHOA_All_Muts[var_freq >= min_var ]
  Muts_mod <- Check_Original_Replicates(Muts_this,osID,rsID,min_var_rep=min_var_replicate)
  Muts_mod <- Muts_mod[depth >= min_cov & RepDepth >= min_cov]
  
  Muts_mod2 <- Check_Original_Replicates(Muts_this,rsID,osID,min_var_rep=min_var_replicate)
  Muts_mod2 <- Muts_mod2[depth >= min_cov & RepDepth >= min_cov]
  Muts_all <- rbind(Muts_mod,Muts_mod2)
  
  # Muts_all <-  distinct(Muts_all, Sample_ID_mut, position, segment_num, var_nt, .keep_all = TRUE)
  # dim(Muts_all)
  # head(Muts_all)
  iMuts_all <- Muts_all[Presence_corrected == "Y"]
  
  # Calculate the total number of mutations and the number of shared mutations
  niSNV <- nrow(Muts_all)
  nSharediSNV <- nrow(iMuts_all)
  
  # Calculate the proportion of shared mutations
  pShared <- nSharediSNV / niSNV
  
  r_squared <- Plot_muts(iMuts_all,min_cov,min_var,min_var_replicate)
  Replicate_Concordance_tmp <- data.frame(min_var_freq = min_var,min_rep_freq = min_var_replicate,min_cov=min_cov,pShared = pShared,r_squared = r_squared)
  Replicate_Concordance_tmp
  Replicate_Concordance <- rbind(Replicate_Concordance,Replicate_Concordance_tmp)
  rm(Replicate_Concordance_tmp,Muts_this, iMuts_all,r_squared)
}

# Plot Replicate concordance 
var_concordance <- ggplot2::ggplot(Replicate_Concordance, aes(x = min_var_freq, y = pShared)) +
  geom_point(alpha = 1, color = "black") +
  labs(x = "Minimum iSNV frequency", y = "Proportion of Variants Shared") +
  scale_x_continuous(limits = c(0.005, 0.105), breaks = seq(0.01, 0.1, 0.01), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0.0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

print(var_concordance)

scatter_name <- paste0(Figure_OutPut_path,"/ReplicateConcordance_minVar_pShared.eps")
ggplot2::ggsave(scatter_name, plot = var_concordance, width = 6, height = 5)

# Save replicate concordance
write.csv(Replicate_Concordance,paste0(Compiled_OutPut_path,"/CHOA_Replicate_Concordance.csv"))


# Set chosen parameter
min_var_final <- 0.03
min_var_replicate_final <- 0.03
min_cov_final <- 100

# Filter the data using the final parameter choice 
Muts_this <- CHOA_All_Muts[var_freq >= min_var_final & RepFreq >= min_var_replicate_final]
CHOA_Muts_1 <- Check_Original_Replicates(Muts_this,osID,rsID,min_var_rep=min_var_replicate_final)
CHOA_Muts_1 <- CHOA_Muts_1[depth >= min_cov & RepDepth >= min_cov]
CHOA_Muts_1[, avg_var_freq := (var_freq + RepFreq) / 2]
CHOA_Muts_1[, avg_cov := (depth + RepDepth)/2 ]


CHOA_Muts_2 <- Check_Original_Replicates(Muts_this,rsID,osID,min_var_rep=min_var_replicate_final)
CHOA_Muts_2 <- CHOA_Muts_2[depth >= min_cov & RepDepth >= min_cov]
CHOA_Muts_2[, avg_var_freq := (var_freq + RepFreq) / 2]
CHOA_Muts_2[, avg_cov := (depth + RepDepth)/2 ]

# Combine the mutations identified in each replicates
CHOA_Muts_FINAL <- rbind(CHOA_Muts_1,CHOA_Muts_2)
dim(CHOA_Muts_1)
dim(CHOA_Muts_2)
dim(CHOA_Muts_FINAL)

# Keep only distinct rows based on the specified columns
CHOA_Muts_FINAL_FILTERED <- distinct(CHOA_Muts_FINAL, Sample_ID_mut, position, segment_num, var_nt, .keep_all = TRUE)

# Find the non-distinct rows by comparing with the distinct rows
CHOA_Muts_All_dup <- anti_join(CHOA_Muts_FINAL, CHOA_Muts_FINAL_FILTERED)
dim(CHOA_Muts_All_dup)

# Add mutation effect
CHOA_Muts_FINAL <- CHOA_Muts_FINAL_FILTERED
CHOA_Muts_FINAL$effect_any <- ifelse(CHOA_Muts_FINAL$effect == "Nonsynonymous" | CHOA_Muts_FINAL$effect2 == "Nonsynonymous", "Nonsynonymous", "Synonymous")
write.csv(CHOA_Muts_FINAL,paste0(Compiled_OutPut_path,"/FINAL_MutationOutput_minvar3.csv"))

# Calculate summary statistics for identified iSNVs
mean_iSNV_frequency <- mean(CHOA_Muts_FINAL$avg_var_freq)
mean_iSNV_frequency

std_iSNV_frequency <- sd(CHOA_Muts_FINAL$avg_var_freq)
std_iSNV_frequency


# Generate plot for variant frequency comparison
Plot_muts_pub_arrangmement <- function(iMuts, Analysis_date) {
  # Select relevant columns from the input data frame.
  plot_muts <- iMuts[, .(Sample_ID_mut, var_freq, RepFreq)]
  
  # Fit a linear model predicting RepFreq from var_freq.
  fit <- lm(RepFreq ~ var_freq, data = plot_muts)
  
  # Calculate the R-squared value from the linear model fit to represent the goodness of fit.
  r_squared <- summary(fit)$r.squared
  
  # Create a scatter plot of var_freq against RepFreq using ggplot2.
  var_comp <- ggplot2::ggplot(plot_muts, aes(x = var_freq, y = RepFreq)) +
    geom_point(alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + 
    labs(x = "iSNV frequency", y = "Replicate iSNV frequency") +
    scale_x_continuous(limits = c(0.02, 0.505), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0.02, 0.505), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) + 
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10))
  # Display the created scatter plot.
  print(var_comp)
  
  # Define the filename for the saved plot
  scatter_name <- paste0(Analysis_date,"-Replicate_Frequency_Comparison_VPC.pdf")
  scatter_eps <- paste0(Analysis_date,"-Replicate_Frequency_Comparison_VPC.eps")
  # Save the scatter plot as pdf
  ggplot2::ggsave(file.path(Figure_OutPut_path,scatter_name),plot = var_comp, width = 6, height = 5)
  
  # Save the scatter plot as eps
  ggplot2::ggsave(file.path(Figure_OutPut_path,scatter_eps), plot = var_comp, width = 6, height = 5, device = "eps")
  
  # Return the R-squared value of the linear model.
  return(var_comp)
}

replicate_concordance_plot <- Plot_muts_pub_arrangmement(iMuts = CHOA_Muts_FINAL, Analysis_date = "81123")


# Save the finalized sample set by removing Bad_IDs_All from CHOA_Selected_Samples
CHOA_SampleSet_FINAL <- CHOA_Selected_Samples[!CHOA_Selected_Samples$ID %in% Bad_IDs_all,]
write.csv(CHOA_SampleSet_FINAL, paste0(Compiled_OutPut_path,"/CHOA_Muts_FINAL.csv"))

# Save parameters chosen run parameters -----
# Generate a dataframe to store parameters called CHOA_RunParam
CHOA_RunParam <- data.frame(
  min_var = min_var_final,
  min_var_replicate = min_var_replicate_final,
  min_cov = min_cov_final
)
write.csv(CHOA_RunParam, file = paste0(Compiled_OutPut_path,"/CHOA_RunParam.csv"), row.names = FALSE)

# Incorporate metadata --------------
CHOA_metadata <- read.csv("CLEAN_CHOA_metadata.csv")
CHOA_metadata$X <- NULL

merged_data<- merge(CHOA_metadata, Original_Samples, by.x = "SubjectID", by.y = "ID")
merged_data$Subject_Original <- Original_Samples$Sample
merged_data$Subject_Replicate <- Replicate_Samples$Sample


# Generate mutation frequency plot -------

# Create a histogram of avg_var_freq from Muts_VPC
iSNV_freq_hist <- ggplot2::ggplot(CHOA_Muts_FINAL, aes(x = avg_var_freq)) +
  geom_histogram(breaks = seq(0, 0.5, by = 0.01), fill = "white", color = "black") +
  labs(x = "iSNV frequency", y = "Count") +
  scale_x_continuous(limits = c(0, 0.51), breaks = seq(0, 0.5, by = 0.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 62), breaks = seq(0, 60, by = 20), expand = c(0, 0)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Display the plot
print(iSNV_freq_hist)

# Save the plot in eps format
ggplot2::ggsave(filename = paste0(Figure_OutPut_path,"/iSNV_frequency_distribution.eps"), plot = iSNV_freq_hist, device = "eps", width = 6, height = 5, units = "in")


# Get non-synonymous mutations for the Muts_VPC table for HA and NA and save them as csv files
NonSynonymous_Muts_HA <- CHOA_Muts_FINAL[effect_any == "Nonsynonymous" & segment_num == 4]
NonSynonymous_Muts_HA
write.csv(NonSynonymous_Muts_HA,paste0(Compiled_OutPut_path,"/NonSynonymous_Muts_HA.csv"))

NonSynonymous_Muts_NA <- CHOA_Muts_FINAL[effect_any == "Nonsynonymous" & segment_num == 6]
NonSynonymous_Muts_NA
write.csv(NonSynonymous_Muts_HA,paste0(Compiled_OutPut_path,"/NonSynonymous_Muts_NA.csv"))

# Determine the number of iSNVs per sample ----------------------
merged_data$all_iSNV_present <- 0
merged_data$nonSyn_iSNV_present <- 0
merged_data$Syn_iSNV_present <- 0
merged_data$nonSyn_HA_present <- 0
merged_data$nonSyn_NA_present <- 0

si <- 1

for (si in 1:length(merged_data$SubjectID)) {
  # count number of matching entries in Var_Present
  sample_id <- merged_data$SubjectID[si]
  sampleSNV_present <- CHOA_Muts_FINAL[CHOA_Muts_FINAL$Sample_ID_mut == sample_id,]
  sampleSNV_nonSynonymous <- CHOA_Muts_FINAL[CHOA_Muts_FINAL$Sample_ID_mut == sample_id & CHOA_Muts_FINAL$effect_any == "Nonsynonymous",]
  
  merged_data$all_iSNV_present[si] <- nrow(sampleSNV_present)
  merged_data$nonSyn_iSNV_present[si] <- nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous",])
  merged_data$Syn_iSNV_present[si] <- nrow(sampleSNV_present[sampleSNV_present$effect_any == "Synonymous",])
  merged_data$nonSyn_HA_present[si] <-  nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous" & sampleSNV_present$segment_num == 4,])
  merged_data$nonSyn_NA_present[si] <-  nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous" & sampleSNV_present$segment_num == 6,])
  merged_data$mean_varFreq_all[si] <- mean(sampleSNV_present$avg_var_freq)
  merged_data$mean_varFreq_nonSyn[si] <- mean(sampleSNV_nonSynonymous$avg_var_freq)
}


# Generate association plot for all iSNVs and CT
# Create scatter plot
CT_iSNV_plot <- ggplot2::ggplot(merged_data, aes(x = CT, y = all_iSNV_present)) +
  geom_point(color = "black", fill = "black", shape = 21, size = 1.25) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_x_continuous(limits = c(15, 30), breaks = seq(15, 30, by = 5), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 82), breaks = seq(0, 80, by = 10), expand = c(0, 0)) +
  labs(x = "CT", y = "iSNV per Sample") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"),
        axis.title = element_text(size=12))


# Display the plot
print(CT_iSNV_plot)

lm_model <- lm(all_iSNV_present ~ CT, data = merged_data)
r_squared <- summary(lm_model)$r.squared

# Print the r^2 value
print(paste("R-squared value:", r_squared))

# Save the plot in eps format
ggplot2::ggsave(filename = paste0(Figure_OutPut_path,"/iSNV_CT_association.eps"), plot = CT_iSNV_plot, device = "eps", width = 6, height = 5, units = "in")

# Arrange the components of the QC figure
final_QC_plot <- cowplot::plot_grid(replicate_concordance_plot, iSNV_freq_hist, CT_iSNV_plot, ncol = 1, nrow = 3, align = "v", rel_widths = c(1, 1))
print(final_QC_plot)

# Save the plot
ggplot2::ggsave(paste0(Figure_OutPut_path,"/final_QC_plot.eps"), plot = final_QC_plot, width = 3.5, height = 7, device = "eps")


# Generate final merged data set ------------
write.csv(merged_data,paste0(Compiled_OutPut_path,"/FINAL_CHOA_merged_data.csv"))
merged_data <- read.csv(file.path("PostProcessing_OutPut/CompiledOutPut/FINAL_CHOA_merged_data.csv"))

# Generate sample specific data table
sample_DataTable <- data.frame(SubjectID = merged_data$SubjectID,all_iSNV = merged_data$all_iSNV_present, 
                               mean_VarFreq = merged_data$mean_varFreq_all,
                               nonSynonymous_iSNV = merged_data$nonSyn_iSNV_present,
                               Synonymous_iSNV = merged_data$Syn_iSNV_present,
                               mean_NonSynon_VarFreq = merged_data$mean_varFreq_nonSyn,
                               nonSynoymous_HA = merged_data$nonSyn_HA_present,
                               nonSynonymous_NA = merged_data$nonSyn_NA_present)
sample_DataTable <- sample_DataTable[order(sample_DataTable$all_iSNV, decreasing = TRUE),]
row.names(sample_DataTable) <- NULL
Name_datatable <- paste0(Compiled_OutPut_path,"/","FinalMuts-CHOA_SampleDataTable.csv")
write.csv(sample_DataTable,Name_datatable)


