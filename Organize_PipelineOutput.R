

# Clear all variables
rm(list = ls()) 

# Get pathway to working directory
R_path <- getwd()
dir.create(file.path(R_path,"Pipeline_Output_Files"))
dir.create(file.path(R_path,"PostProcessing_OutPut"))
# Prepare environment -----------------------------

# Load packages
library("tidyverse")
library("Rsamtools")
library("Biostrings")
library("GenomicRanges")
library("GenomicAlignments")
library("dplyr")
library("seqinr")

# Define functions -----
create_dir <- function(parent_path, dir_name) { # Helper function to create directories and return their paths.
  full_path <- file.path(parent_path, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path)
  }
  return(full_path)
}

# Function to copy or move files from one directory to another
copy_or_move_files <- function(from_files, to_dir, move=FALSE) {
  lapply(from_files, function(file) {
    to_file <- file.path(to_dir, basename(file))
    if (move) {
      tryCatch(
        file.rename(from = file, to = to_file),
        error = function(e) cat("Error moving file:", file, "Error message:", e$message, "\n")
      )
    } else {
      tryCatch(
        file.copy(from = file, to = to_file, overwrite = TRUE),
        error = function(e) cat("Error copying file:", file, "Error message:", e$message, "\n")
      )
    }
  })
  cat(paste0(ifelse(move, "Moved ", "Copied "), length(from_files), " files to ", to_dir, "\n"))
}

# Pipeline_organisation_run_options -----------------------------

RunDate <- c(33122,50422,61722,91222,22323,31023,50123) # These are the run dates of the sample tracker to include
RD <- 1
min_Qual <- 30 #Minimum PHRED score
min_Map <- 40  #Minimum mapping quality
min_Coverage <- 100 #minimum coverage 

for (RD in 2:length(RunDate)){
  Pipeline_date <- RunDate[RD]
  Pipeline_data_path <- paste0("/Volumes/Elements/FINALIZED_CHOA_Run/IAV_Run_",Pipeline_date,"/",Pipeline_date,"_NGS_FINAL") # CHOP path
  
  # list all directories in the specified path
  dir_list <- list.dirs(path = Pipeline_data_path, recursive = TRUE)
  dir_list
  # identify directories containing the pattern "ivar"
  iVar_dirs <- grep(pattern = "ivar", x = dir_list, value = TRUE)
  iVar_dirs
  
  # Creating iVar directory and moving iVar_dirs
  new_iVar_path <- create_dir(Pipeline_data_path,"iVar")
  copy_or_move_files(from_files = iVar_dirs, to_dir = new_iVar_path, move = TRUE)
  
  # Paths for the 'sampleOutputs' and 'iVar'
  new_sample_path <- file.path(Pipeline_data_path,"sampleOutputs")
  new_intrahost_path <- file.path(Pipeline_data_path,"iVar")
  
  # Create main directories for Pipeline Output
  new_dir_R_path <- create_dir(R_path, paste0("Pipeline_Output_Files/",Pipeline_date,"_PipelineOutput"))
  
  new_consensus_path <- create_dir(new_dir_R_path,"Consensus")
  new_basecov_path <- create_dir(new_dir_R_path,"BaseCov")
  new_Variants_path <- create_dir(new_dir_R_path,"Variants")
  
  # Create subdirectories under Variants
  new_variants_intrahost_path <- create_dir(new_Variants_path,"Intrahost")
  new_variants_consensus_path <- create_dir(new_Variants_path,"Consensus")

  
  # Create main directory for pileup files
  new_pileup_path <- create_dir(new_dir_R_path,"Pileup")
  
  # Get lists of files to move and make sure the contents are equivalent -----------------------------
  all_intrahost_passing_variants <- list.files(path=new_intrahost_path,pattern="passingVariants.csv",recursive = TRUE)
  IPV_name <- str_extract(all_intrahost_passing_variants, "CHOA-\\d+_D\\dV\\d_S\\d+")
  IPV_name <- IPV_name[!is.na(IPV_name)]
  
  # Create expression to extract matching files
  pattern <- paste(IPV_name,collapse="|")
  intrahost_passing_variants <- all_intrahost_passing_variants[str_detect(all_intrahost_passing_variants,pattern)]
  
  all_consensus_sequences <- list.files(path=new_sample_path,pattern="consensus_sequence.fasta",recursive = TRUE) 
  consensus_sequences <- all_consensus_sequences[str_detect(all_consensus_sequences,pattern)]; consensus_sequences
  
  all_base_cov <- list.files(path=new_sample_path,pattern="basecov.txt",recursive = TRUE)
  base_cov <- all_base_cov[str_detect(all_base_cov, pattern)]
  
  all_consensus_passing_variants <- list.files(path=new_sample_path,pattern="passingVariants.csv",recursive = TRUE)
  consensus_passing_variants <- all_consensus_passing_variants[str_detect(all_consensus_passing_variants,pattern)]
  
  length_vec <- c(length(IPV_name),length(intrahost_passing_variants),length(consensus_passing_variants),length(consensus_sequences),length(base_cov))
  all_lengths <- all(length_vec == length_vec[1]); print(all_lengths)
  
  if (all_lengths == FALSE){
    stop("There are a different number of files")
  }
  
  # Move files -----------------------------
  
  # Generate the full paths for consensus files and copy those files
  consensus_full_paths <- file.path(new_sample_path, all_consensus_sequences)
  copy_or_move_files(from_files = consensus_full_paths, to_dir = new_consensus_path, move = FALSE)
  
  # Create the full paths for base_cov files and copy those files
  base_cov_files <- file.path(new_sample_path, base_cov)
  copy_or_move_files(from_files = base_cov_files, to_dir = new_basecov_path, move = FALSE)
  
  # Generate paths and copy consensus variants
  consensus_variants_files <- file.path(new_sample_path,consensus_passing_variants)
  copy_or_move_files(from_files = consensus_variants_files, to_dir = new_variants_consensus_path, move = FALSE)
  
  # Generate paths and copy intrahost variants
  intrahost_variants_files <- file.path(new_iVar_path,intrahost_passing_variants)
  copy_or_move_files(from_files = intrahost_variants_files, to_dir = new_variants_intrahost_path, move = FALSE)
  
  
  # iVar_bam -----------------------------
  # Create directory for the bam/bai files within the iVar folder for the run. These are big files, so don't want them in the R folder
  new_bam_path_iVar <- create_dir(new_iVar_path,"BamFiles")
  
  # Get a list of bam files in the iVar directory
  bam_files <- list.files(path=new_iVar_path,pattern=".bam",recursive = TRUE)
  bam_full_paths <- file.path(new_iVar_path, bam_files)
  copy_or_move_files(from_files = bam_full_paths, to_dir = new_bam_path_iVar, move = TRUE)
  
  
  # Consensus level bam files -----------------------------
  new_bam_path_consensus <- create_dir(new_sample_path,"BamFiles")
  
  # Get a list of base coverage files in the sample output directory
  bam_files_consensus <- list.files(path=new_sample_path,pattern=".bam",recursive = TRUE)
  bam_consensus_full_paths <- file.path(new_sample_path, bam_files_consensus)
  copy_or_move_files(from_files = bam_consensus_full_paths, to_dir = new_bam_path_consensus, move = TRUE)
  
  #  Generate readcounts from pileup files -----------------------------
  # Define pileup parameters
  mypileup_param <- PileupParam(min_base_quality = min_Qual, min_mapq = min_Map, min_nucleotide_depth = 1,max_depth = 1000000, distinguish_strands = FALSE, distinguish_nucleotides = TRUE, include_deletions = FALSE, include_insertions = FALSE, ignore_query_Ns = TRUE)
  
  # Extract the bam and index files; big_list contains bams and bais, which are sorted in the next 2 steps
  big_list <- list.files(path = new_bam_path_consensus, pattern = ".sorted.bam") 
  bam_list <- big_list[!grepl(".bai",big_list)]
  index_list <- big_list[grepl(".bai",big_list)]
  
  # Identify the consensus sequences for the bam files, which are already copied to a new directory in R
  consensus_list <- list.files(path = new_consensus_path, pattern = "consensus_sequence.fasta")
  
  bl <- 1
  for (bl in 1:length(bam_list)){
    # Generate name for pileup file
    name_bam <- unlist(strsplit(bam_list[bl],"_")); name_bam
    pileup_name <- paste(name_bam[1],"_",name_bam[2],"_Pileup.csv",sep="")
    print(pileup_name)
    # Generate pileup file
    bamfile <- bam_list[bl]
    indexfile <- index_list[bl]
    bs <- BamFile(bamfile)
    pileupfile <- pileup(file=paste(new_bam_path_consensus,"/",bam_list[bl],sep=""),index=paste(new_bam_path_consensus,"/",index_list[bl],sep=""),pileupParam = mypileup_param)
    
    # Save pileup file in allocated folder
    pileup_path <- file.path(new_pileup_path,pileup_name)
    write.csv(x = pileupfile,file=pileup_path)
  }
  
  # Create file for renamed pileup files in the main R directory
  pileup_renamed_path <- create_dir(R_path,"PostProcessing_OutPut/Pileup_renamed")
  pileup_list <- list.files(path=new_pileup_path,pattern="Pileup.csv",recursive=TRUE)
  nt <- 1
  for (nt in 1:length(pileup_list)){
    pileup_name <- pileup_list[nt];
    pileup_file <- read.csv(file.path(new_pileup_path,pileup_name))
    head(pileup_file)
    
    # Remove index column
    pileup_file$X <- NULL
    SegNum <- as.double(str_sub(pileup_file$seqnames,12,12)); head(SegNum) 
    
    # Add segment number
    new_pileup <- data.frame(SegNum=SegNum,seqnames=pileup_file$seqnames,pos=pileup_file$pos,nucleotide=pileup_file$nucleotide,count=pileup_file$count) 
    rownames(new_pileup) <- NULL
    new_path <- file.path(pileup_renamed_path,pileup_name)
    write.csv(new_pileup,new_path)
  }
  
}
