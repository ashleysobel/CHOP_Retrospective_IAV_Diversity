
# Clear environment to get rid of old files and get path -------------------
rm(list = ls())
R_path <- getwd()


# Load Reference Key (matches PB2 accession # with strain). Make sure this has been updated if your run contains different reference files
Reference_key <- read.csv("Reference_files/Reference_key.csv")

# Set segment names
segment_name <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")

# Load packages
# library(stringr) #use for str_replace_all, str_replace
library(phylotools) # use for read.fasta
library(dplyr) # use for select
library(stringi) # use for stri_extract_last
library(ggplot2) # Use for generating the coverage plots
# library(gridExtra) # Use for grid.arrange
# library(grid) # used for textGrob
library(seqinr) # used for write.fasta, translate
#library(GenomicRanges)
#library(ShortRead) 
#library(ape)
library(tidysq)
#library(devtools)
#library(scales)
#library(cowplot)
library(patchwork)
library(purrr)

# Set QC parameters -----------------------


# Enter the csv file containing the references you want to use. For the CHOP
# 2017-2018 retrospective study, we want to use H3N2 2022 and H1N1 2019.
H3N2_Ref <- read.csv("Reference_files/H3N2_Ref_A_2016.csv")
H1N1_Ref <- read.csv("Reference_files/H1N1_Ref_A_2015.csv")

H3N2_Ref_consensus <- phylotools::read.fasta("Reference_files/A_Washington_17_2016_H3N2.fasta")
H1N1_Ref_consensus <- phylotools::read.fasta("Reference_files/A_Michigan_45_2015_H1N1_ref.fasta")

min_coverage <- 100 # This is the minimum coverage we will accept for QC purposes
RunDate <- c(33122,50422,61722,91222,22323, 31023, 50123) # These are the run dates of the sample tracker to include

options(warn=1)
# Define functions ------------------------------------------

GetFasta <- function(PipelineOutput_path) {
  path_fasta <- paste0(PipelineOutput_path,"/Consensus"); path_fasta
  fasta_names <- list.files(path_fasta, pattern = "\\.fasta"); fasta_names
  fasta_info <- list(path_fasta=path_fasta,fasta_names=fasta_names); fasta_info
}


GetCov <- function(R_path, Pipeline_date) {
  path_basecov <- paste0(PipelineOutput_path,"/BaseCov"); path_basecov
  base_cov <- list.files(path_basecov, pattern = "basecov.txt"); base_cov
  cov_info <- list(path_basecov=path_basecov,basecov_names=base_cov)
}

GetCovPlot_path <- function(R_path,Pipeline_date,sample_name) {
  CovPlot_path <- paste0(CovPlot_path,"/",sample_name,"_Coverage.pdf")
} 

# Checks for differences between two strings
GetStringDiff<- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE){
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case)
  {
    a <- toupper(a)
    b <- toupper(b)
  }
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
  ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
  if(!show.excluded) diff.info<-na.omit(diff.info)
  diff.info
}

# Function to plot the read coverage for target sample, organized by segment
Plot_Segment_Coverage <- function(sorted_base_cov, segment_index, segment_name) {
  segment <- subset(sorted_base_cov, sorted_base_cov$Segment == unique(sorted_base_cov$Segment)[segment_index])
  segment_Position <- as.numeric(segment$Position)
  segment_Coverage <- as.numeric(segment$Coverage)
  segment_df <- data.frame(cbind(segment_Position,segment_Coverage))
  segment_plot <- ggplot(segment_df, aes(x=segment_Position, y=segment_Coverage)) + 
    geom_area() +
#    scale_y_log10() + # add this line to set y-axis to log scale
    ggtitle(segment_name[segment_index]) + 
    theme(plot.title=element_text(hjust=0.5)) +
    labs(x = "Position",y = "Coverage") 
  return(segment_plot)
}

# Function to write fasta file (will be used for creating fasta file for masked consensus sequence)
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

Characterize_mutations <- function(stmp, subject, Strain, Ref) {
  sample_mat <- list()
  mutation_output_tmp <- data.frame()
  # Filter mutations for substitutions only 
  subs <- stmp[stmp$TYPE == "SUB", ]
  nSub <- nrow(subs) # Number of substitutions
  print(nSub)
  
  if (nSub != 0){
    mt <- 1
    # if (mt %% 10 == 0){
    #   print(mt)
    # }
    # Evaluate each mutation
    for (mt in 1:nSub){
      sub_tmp <- subs[mt,]; sub_tmp
      accession <- strsplit(sub_tmp$CHROM,'-')[[1]][1]
      segment_num <- as.numeric(strsplit(sub_tmp$CHROM,'_')[[1]][3])
      position <- as.numeric(sub_tmp$POS)
      ref_nt <- sub_tmp$REF
      var_nt <- sub_tmp$ALT
      var_freq <- as.numeric(sub_tmp$AF)
      var_qual <- as.numeric(sub_tmp$QUAL)
      # The variant output from the intrahost and consensus level variant callers is slightly different w/respect to mapping quality (MQ), we must update the intrahost variant caller MQM to read as MQ
      if (is.null(sub_tmp$MQ)){
        sub_tmp$MQ <- sub_tmp$MQM
      }
      var_map <- as.numeric(sub_tmp$MQ)
      var_type <- sub_tmp$TYPE
      if (var_type != "SUB"){
        stop("There is an error!!!!!") # We're only considering substitutions for now, so this would be bad
      }
      depth <- as.numeric(sub_tmp$DP)
      var_depth <- as.numeric(sub_tmp$AD)
      var_info <- sub_tmp$DP4
      mut <- data.frame(subject,Strain,segment_num,position,ref_nt,var_nt,var_freq,var_qual,var_map,var_type,depth,var_depth,var_info)
      mut
      
      # Extract reference segment 
      Ref_seg <- Ref[Ref$Segment == mut$segment_num, ]
      mut$position
      c(mut$position,Ref_seg$Start,Ref_seg$Stop)
      
      if (mut$position < Ref_seg$Start) {
        # We define variant that occur outside the coding regions to be: 1)
        # noncoding, 2) resulting in a synonymous mutation, and have no relevant
        # amino acid position/identity
        location <- "NonCoding"
        effect <- NaN
        aa_pos <- NaN
        aa_sub  <- NaN
        aa_ref <- NaN
        
      } else if (mut$position > Ref_seg$Stop) {
        location <- "NonCoding"
        effect <- NaN
        aa_pos <- NaN
        aa_sub  <- NaN
        aa_ref <- NaN
      } else {
        # This variant occurs inside the coding region.
        location <- Ref_seg$Protein # List location as the protein
        seq_ref <- Ref_seg$NT # Get the sequence of the reference
        
        ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
        ref_CDS <- bite(ref_nt_seg, indices = Ref_seg$Start:Ref_seg$Stop) # Extract the CDS 
        ref_AA <- c2s(seqinr::translate(s2c(as.character(ref_CDS))))
        
        seq_sub <- seq_ref # Create the sample sequence
        seq_sub
        substr(seq_sub,mut$position,mut$position) <- mut$var_nt
        
        sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
        
        sub_CDS <- bite(sub_nt, indices = Ref_seg$Start:Ref_seg$Stop) # Extract the CDS
        sub_AA <- c2s(seqinr::translate(s2c(as.character(sub_CDS))))
        
        aa_pos <- floor((mut$position - Ref_seg$Start)/3) + 1 # Get aa residue (Thanks to Andrew for code!)
        aa_ref <- substr(x = ref_AA, start = aa_pos, stop = aa_pos) # Reference amino acid
        aa_sub <- substr(x = sub_AA, start = aa_pos, stop = aa_pos) # Subject amino acid
        c(aa_pos, aa_ref, aa_sub)
        if (aa_ref == aa_sub){
          effect <- "Synonymous"
        } else if (aa_ref != aa_sub){
          effect <- "Nonsynonymous"
        }
      }  
      mutation_effect1 <- data.frame(location,effect,aa_pos,aa_ref,aa_sub); mutation_effect1
      
      # Check for a second protein on the gene segment
      if (Ref_seg$Protein2 ==""){
        # There is a single protein in that gene segment, the rest is NaN
        location2 <- NaN; aa_pos2 <- NaN; effect2 <- NaN; aa_ref2 <- NaN; aa_sub2 <- NaN;
      } else {
        intervals <- c(Ref_seg$Start2_1,Ref_seg$Stop2_1,Ref_seg$Start2_2,Ref_seg$Stop2_2); intervals
        mut$position
        if (mut$position >= intervals[1] && mut$position <= intervals[2]){
          # Variant is within a coding region for the 2nd protein
          location2 <- Ref_seg$Protein2
          seq_ref <- Ref_seg$NT # Get the sequence of the reference
          seq_sub <- seq_ref # Create the sample sequence
          position_check_seq <- seq_sub
          substr(seq_sub,mut$position,mut$position) <- mut$var_nt
          substr(position_check_seq,mut$position,mut$position) <- 'Q'
          
          
          if (is.na(intervals[3]) == FALSE){
            seq_part1 <- substr(position_check_seq,intervals[1],intervals[2])
            seq_part2 <- substr(position_check_seq,intervals[3],intervals[4])
            seq_whole <- paste(seq_part1,seq_part2,sep=""); seq_whole
            position_new <- unlist(gregexpr('Q', seq_whole)); position_new
            
            ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
            ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
            ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2))))
            sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
            sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4]))
          } else {
            seq_part1 <- substr(position_check_seq,intervals[1],intervals[2])
            seq_whole <- seq_part1; seq_whole
            position_new <- unlist(gregexpr('Q', seq_whole)); position_new
            
            ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
            ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2])) # Extract the CDS using the indices
            ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2))))
            sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
            sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2]))
          }
          sub_AA2 <- c2s(seqinr::translate(s2c(as.character(sub_CDS2))))
          aa_pos2 <- floor((position_new - 1)/3) + 1 # Get aa residue (Thanks to Andrew for code!)
          aa_ref2 <- substr(x = ref_AA2, start = aa_pos2, stop = aa_pos2) # Reference amino acid
          aa_sub2 <- substr(x = sub_AA2, start = aa_pos2, stop = aa_pos2) # Subject amino acid
          if (ref_AA2 == sub_AA2){
            effect2 <- "Synonymous"
          } else if (ref_AA2 != sub_AA2){
            effect2 <- "Nonsynonymous"
          }
        } else if (is.na(intervals[3]) == FALSE  && mut$position >= intervals[3] && mut$position <= intervals[4]){
          # Variant is within a coding region for the 2nd protein
          location2 <- Ref_seg$Protein2
          seq_ref <- Ref_seg$NT # Get the sequence of the reference
          seq_sub <- seq_ref # Create the sample sequence
          position_check_seq <- seq_sub
          substr(seq_sub,mut$position,mut$position) <- mut$var_nt
          substr(position_check_seq,mut$position,mut$position) <- 'Q'
          seq_part1 <- substr(position_check_seq,intervals[1],intervals[2])
          seq_part2 <- substr(position_check_seq,intervals[3],intervals[4])
          seq_whole <- paste(seq_part1,seq_part2,sep=""); seq_whole
          position_new <- unlist(gregexpr('Q', seq_whole)); position_new
          
          ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
          ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
          ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2))))
          sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
          sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4]))
          sub_AA2 <- c2s(seqinr::translate(s2c(as.character(sub_CDS2))))
          aa_pos2 <- floor((position_new - 1)/3) + 1 # Get aa residue (Thanks to Andrew for code!)
          aa_ref2 <- substr(x = ref_AA2, start = aa_pos2, stop = aa_pos2) # Reference amino acid
          aa_sub2 <- substr(x = sub_AA2, start = aa_pos2, stop = aa_pos2) # Subject amino acid
          if (ref_AA2 == sub_AA2){
            effect2 <- "Synonymous"
          } else if (ref_AA2 != sub_AA2){
            effect2 <- "Nonsynonymous"
          }
        }
        else {
          location2 <- "NonCoding"
          effect2 <- NaN
          aa_pos2 <- NaN
          aa_sub2 <- NaN
          aa_ref2 <- NaN
        }
      }
      mutation_effect2 <- data.frame(location2,effect2,aa_pos2,aa_ref2,aa_sub2); mutation_effect2
      mutation_mat <- data.frame(mut,mutation_effect1,mutation_effect2)
      sample_mat[[mt]] <- mutation_mat
    }
  }
  # Combine all the data frames in the list
  mutation_output_tmp <- do.call(rbind, sample_mat)
  return(mutation_output_tmp)
}

create_dir <- function(parent_path, dir_name) { # Helper function to create directories and return their paths.
  full_path <- file.path(parent_path, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path)
  }
  return(full_path)
}

# Run QC ------------------------------------------

RunDate <- c(33122,50422,61722,91222,22323,31023,50123) # These are the run dates of the sample tracker to include
RD <- 1
min_Qual <- 30 #Minimum PHRED score
min_Map <- 40  #Minimum mapping quality
min_Coverage <- 100 #minimum coverage 

for (RD in 1:length(RunDate)){

    Pipeline_date <- RunDate[RD]
    
    # Generate data structures to hold the QC output ------------------------------
    PipelineOutput_path <- paste0(R_path,"/Pipeline_Output_Files/",Pipeline_date,"_PipelineOutput"); PipelineOutput_path
    CovPlot_path <- create_dir(PipelineOutput_path, "CovPlot")
    MaskedConsensus_path <- create_dir(PipelineOutput_path, "MaskedConsensus")
    
    # Define empty vectors to hold your data
    prop_N <- numeric()
    sample_list <- character()
    pass_Coverage <- character()
    sequencing_date <- character()
    sequence_strain <- character()
    PB2_accession_list <- character()
    
    
    # Load the fasta and coverage information from this run 
    fasta_info <- GetFasta(PipelineOutput_path)
    cov_info <- GetCov(PipelineOutput_path)
    
    # Get list of the samples with both fasta files AND base_cov files
    fasta_info$fasta_names
    cov_info$basecov_names
    
    # Extracting the conserved part of the name from cov_info$basecov_names
    cov_conserved_names <- sub("_basecov\\.txt$", "", cov_info$basecov_names)
    cov_conserved_names
    
    # Extracting the conserved part of the name from fasta_info$fasta_names
    fasta_conserved_names <- sub("_consensus_sequence\\.fasta$", "", fasta_info$fasta_names)
    fasta_conserved_names

    # Finding the entries that exist in both vectors
    not_common_entries <- !(fasta_conserved_names %in% cov_conserved_names)
    
    # Remove entries from fasta_info that are not in common
    fasta_info$fasta_names <- fasta_info$fasta_names[!not_common_entries]

    # Rename Fasta files
    n <- 1
    for (n in 1:length(fasta_info$fasta_names)){
      parts <- strsplit(fasta_info$fasta_names[n],"_")[[1]]
      parts
      sample_list[n] <- paste0(parts[[1]],"_",parts[[2]])
      sequencing_date[n] <- as.character(as.integer(Pipeline_date))
    }
    
    # Generate Sample Tracker
    ind <- 1
    for (ind in 1:length(sample_list)){
      # Read in fasta -----------------------------------------------------------
      fasta_name_path <- paste0(fasta_info$path_fasta,"/",fasta_info$fasta_names[ind])
      tmp_fasta <- phylotools::read.fasta(fasta_name_path) # Phylotools package
      
      # Get the order that the fasta file should be in. This is done by extacting
      # the final character in the accession #, which corresponds to the segment
      # number
      qzar_last <- as.numeric(stringr::str_sub(tmp_fasta$seq.name, start = -1, end = -1))
      tmp_fasta_new <- cbind(tmp_fasta,qzar_last)
      tmp_qzar <- tmp_fasta_new[order(tmp_fasta_new$qzar_last),]
      tmp_fasta <- data.frame(seq.name=tmp_qzar$seq.name,seq.text=tmp_qzar$seq.text)
      if (tmp_fasta$seq.name[1] == "CY087752_1_4"){
        # This sequence mapped to H6N1, so we need to get rid
        to_delete <- sample_list[ind]
        file.remove(list.files(pattern=to_delete, recursive=TRUE))
        next
        
      }
      
      # Create fasta file with new names ------------------------------------------------
      # We set the name for the new fasta file  
      sample_name <- sample_list[ind]
      fasta_name_string <- paste0(MaskedConsensus_path,"/",sample_name,"_Masked.fasta")
      # Rename the individual gene segments
      old_name <- tmp_fasta$seq.name
      new_name <- paste(sample_name,segment_name,sep = "_") # Generate new names
      ref2 = data.frame(old_name,new_name)
      new_fasta <- paste(sample_name,"_",Pipeline_date,".fasta",sep = "")
      
      
      # Get base coverage -------------------------------------------------------
      base_cov_path <- paste(cov_info$path_basecov,"/",cov_info$basecov_names[ind],sep="")
      base_cov <- read.delim(base_cov_path)
      colnames(base_cov) <- c("Segment","Position","Coverage")
      sorted_base_cov <- base_cov[order(stringi::stri_extract_last(base_cov$Segment,regex="\\w")),]
      
      # Convert the coverage and position columns to numeric so they can be manipulated
      base_cov$Coverage <- as.numeric(base_cov$Coverage)
      base_cov$Position <- as.numeric(base_cov$Position) 
      
      # Need to add +1 to the position column b/c it starts with 0
      base_cov$Position <- base_cov$Position + 1
      
      # Save the updated Base Coverage file with corrected position as a CSV file
      new_base_cov_path <- gsub('basecov.txt','updated_basecov.csv',base_cov_path)
      write.csv(base_cov,new_base_cov_path)
      
      # Extract the base coverage and determine which positions have insufficient coverage to be masked
      coverage <- as.double(base_cov$Coverage)
      coverage_fail <- which(coverage < min_coverage) # Identify the positions with insufficient coverage, termed min_coverage
      # Update the consensus sequence so the sites where the coverage is < min_cov, are masked with "N"
      segx <- 1
      for (segx in 1:8){
        seg_accession <- tmp_fasta$seq.name[segx]
        seg_base_cov <- base_cov[base_cov$Segment == seg_accession,]
        seg_cov <- seg_base_cov$Coverage; head(seg_cov) # Get the coverage for a segment
        
        bad_sites <- which(seg_cov < min_coverage); bad_sites # Determine if there are any sites below the minimum coverage
        if (length(bad_sites) != 0){
          fasta_seg <- tmp_fasta[tmp_fasta$seq.name == seg_accession,]; head(fasta_seg) # Get fasta sequence for current segment
          seg_seq <- fasta_seg$seq.text; head(seg_seq)
          
          for (i in bad_sites){
            str_sub(seg_seq,start = i, end = i) <-"N"   
          }
          tmp_fasta$seq.text[segx] <- seg_seq  # Update fasta sequence with the masked version
        } 
      }
      tmp_fasta$seq.name <- ref2$new_name 
      new_fasta_df <- dplyr::data_frame(name=tmp_fasta$seq.name,seq=tmp_fasta$seq.text) #This is weird formatting, but it's the only way it works
      writeFasta(data = new_fasta_df,filename = fasta_name_string)
      
      N_prop <- length(coverage_fail)/length(coverage) # get the proportion of bases with insufficient coverage and save as a vector, indexed by sample
      prop_N[ind] <- N_prop*100
      if (prop_N[ind] > 5){
        pass_Coverage[ind] <- "No"
      } else {
        pass_Coverage[ind] <- "Yes"
      }
      print(prop_N[ind])
      accession <- unique(sorted_base_cov[,1])
      
      PB2_accession <- unlist(strsplit(accession[1],'_'))[1]
      accession_index <- Reference_key[Reference_key$Accession_PB2 == PB2_accession,]
      if (grepl(PB2_accession,accession_index$Accession_PB2,fixed=TRUE) != TRUE){
        stop("Need to update accession list")
      }
      
      if (is_empty(accession_index$Accession_PB2) == FALSE){
        PB2_accession_list[ind] <- accession_index$Accession_PB2
        sequence_strain[ind] <- accession_index$Strain
      } else {
        PB2_accession_list[ind] <- NaN
        sequence_strain[ind] <- NaN
      }
      # Here we generate a table of coverage positions indexed by segment and nucleotide position. 
      for (ns in 1:length(segment_name)){
        tmp_cov <- sorted_base_cov[sorted_base_cov[,1]==accession[1],]
        sorted_base_cov <- data.frame(lapply(sorted_base_cov, function(x){
          gsub(accession[ns],segment_name[ns],x)
        }))
      }
      
      PB2_plot <- Plot_Segment_Coverage(sorted_base_cov,1,segment_name)
      PB1_plot <- Plot_Segment_Coverage(sorted_base_cov,2,segment_name)
      PA_plot <- Plot_Segment_Coverage(sorted_base_cov,3,segment_name)
      HA_plot <- Plot_Segment_Coverage(sorted_base_cov,4,segment_name)
      NP_plot <- Plot_Segment_Coverage(sorted_base_cov,5,segment_name)
      NA_plot <- Plot_Segment_Coverage(sorted_base_cov,6,segment_name)
      MP_plot <- Plot_Segment_Coverage(sorted_base_cov,7,segment_name)
      NS_plot <- Plot_Segment_Coverage(sorted_base_cov,8,segment_name)

      sample_output <- GetCovPlot_path(R_path, Pipeline_date, sample_name)

      pdf(file = sample_output,   # The directory you want to save the file in
          width = 8, # The width of the plot in inches
          height = 10) # The height of the plot in inches

      samplelabel <- paste("Coverage plots for", sample_name, sep = " ")

      # Arrange the plots using patchwork
      arranged_plot <- (PB2_plot | PB1_plot ) /
        (PA_plot | HA_plot ) /
        ( NP_plot | NA_plot) /
        (MP_plot | NS_plot) +
        patchwork::plot_annotation(title = samplelabel, theme = theme(plot.title = element_text(hjust = 0.5)))

      # Print the arranged plot to the device
      print(arranged_plot)

      dev.off()
    }
    
    
    # Generate sample tracker containing all samples
    SampleTracker_all <- data.frame(sample_list,sequence_strain,sequencing_date,PB2_accession_list,prop_N,pass_Coverage)
    colnames(SampleTracker_all) <- c("Sample","Strain","Sequencing Date","Ref_Accession_PB2","%N","Pass QC_Cov")
    rownames(SampleTracker_all) <- NULL
    
    # Save sample tracker as csv file
    tracker_name <- paste0("PostProcessing_OutPut/",as.character(as.integer(Pipeline_date)),"_SampleTracker_all.csv")
    write.csv(x = SampleTracker_all,file = tracker_name)
    
    # Generate a sample tracker for just the CHOA samples
    SampleTracker_CHOA <- SampleTracker_all[str_detect(SampleTracker_all$Sample, "CHOA"),]
    SampleTracker_CHOA <- SampleTracker_CHOA[is.na(SampleTracker_CHOA$Strain) == FALSE,]
    CHOA_tracker_name <- paste0("PostProcessing_OutPut/",as.character(as.integer(Pipeline_date)),"_SampleTracker_CHOA.csv")
    write.csv(x = SampleTracker_CHOA,file = CHOA_tracker_name)
    
    
    # Characterize sample variants -----------------------------
    H3N2_Ref$NT <- gsub("U", "T", H3N2_Ref$NT) # Convert to DNA
    H3N2_Ref$Segment <- as.numeric(H3N2_Ref$Segment)
    H1N1_Ref$NT <- gsub("U", "T", H1N1_Ref$NT) # Convert to DNA
    H1N1_Ref$Segment <- as.numeric(H1N1_Ref$Segment)
    
    # assign path to variants identified in Pipeline
    consensus_variant_path <- file.path(PipelineOutput_path,"Variants/Consensus")
    intrahost_variant_path <- file.path(PipelineOutput_path,"Variants/Intrahost")
    
    samples <- SampleTracker_CHOA$Sample
    N_sample <- length(samples)
    
    mutation_output <- data.frame() # Empty dataframe for holding mutation info 
    segment_name_num <- data.frame(SegNum <- as.double(c(1:8)),Segment <- c("PB2","PB1","PA","HA","NP","NA","MP","NS"))
    colnames(segment_name_num) <- c("SegNum","Segment")
    
    
    ### Run variant classification loops
    new_reference_intra_path <- paste0(PipelineOutput_path,"/Reference/Intrahost"); new_reference_intra_path
    new_reference_pop_path <- paste0(PipelineOutput_path,"/Reference/Consensus"); new_reference_pop_path
    smp <- 1
    for (smp in 1:N_sample){
      # Get information for this sample
      sample_tmp <- samples[smp]
      subject_parts <- unlist(strsplit(sample_tmp, "_"))
      subject <- paste(subject_parts[1],subject_parts[2],sep="_") # Get subject name
      sample_info <- SampleTracker_CHOA[SampleTracker_CHOA$Sample == subject,] # Get the subject information from the sample tracker
      Strain <- sample_info$Strain # Get subject's strain
      
      # Ensure there are no duplicate entries, this would be an error
      if(nrow(sample_info) > 1){
        print("There is a duplicate sample, this is an error")
        {stop(TRUE)}
      }
      Strain <- sample_info$Strain # Get subject's strain
      
      # Load the masked consensus file for the current subject
      name_path <- file.path(MaskedConsensus_path,paste0(sample_info$Sample,"_Masked.fasta")); name_path
      
      # Generate reference file for the current sample -------------------------------------------------
      if (is.na(Strain)){
        print("There is no matching strain.")
        {stop(TRUE)}
      }
      if (Strain == "H3N2"){
        Ref <- H3N2_Ref
        consensus_sequence_pop <- H3N2_Ref_consensus
        consensus_sequence_tmp <- tidysq::read_fasta(name_path,alphabet="dna")  
        consensus_sequence_intrahost <- tibble(seq.name = consensus_sequence_tmp$name,seq.text = consensus_sequence_tmp$sq)
      } else if (Strain == "H1N1"){
        Ref <- H1N1_Ref
        consensus_sequence_pop <- H1N1_Ref_consensus
        consensus_sequence_tmp <- tidysq::read_fasta(name_path,alphabet="dna")  
        consensus_sequence_intrahost <- tibble(seq.name = consensus_sequence_tmp$name,seq.text = consensus_sequence_tmp$sq)
      }
      
      Generate_Reference <- function(Ref, consensus_sequence, segment_name_num) {
        Ref$NT <- NaN; Ref$AA <- NaN; Ref$AA2 <- NaN #Clear original sequences, need to replace w/consensus
        
        # Load the consensus sequence and use it to update the reference
        i <- 1
        for (i in 1:8){
          consensus_tmp <- consensus_sequence[i,]
          consensus_tmp
          # Convert sequence to character vector
          consensus_char <- as.character(consensus_tmp$seq.text)
          consensus_char
          segment_tmp <- as.double(segment_name_num$SegNum[i]); segment_tmp
          Ref$NT[segment_tmp] <- consensus_char #Update NT sequence
          
          Ref_seg <- Ref[segment_tmp,]; Ref_seg
          ref_nt <- as.sq(consensus_char, alphabet = "dna_bsc")
          
          # Update first protein sequence 
          ref_CDS <- bite(ref_nt, indices = Ref_seg$Start:Ref_seg$Stop); ref_CDS # Extract the CDS using the indices
          ref_AA <- c2s(seqinr::translate(s2c(as.character(ref_CDS))))
          Ref$AA[segment_tmp] <- ref_AA
          
          #If there's a second protein, update the second AA sequence
          if (Ref_seg$Protein2 != ""){
            # if (i == 3 & smp == 37){
            #   browser()
            # }
            intervals <- c(Ref_seg$Start2_1,Ref_seg$Stop2_1,Ref_seg$Start2_2,Ref_seg$Stop2_2); intervals  
          if (is.na(intervals[3])){
              # There is a single interval
              ref_CDS_2 <- bite(ref_nt, indices = c(intervals[1]:intervals[2]))
            } else{
              # There are multiple starts
              ref_CDS_2 <- bite(ref_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
            }
            Ref$AA2[segment_tmp] <- c2s(seqinr::translate(s2c(as.character(ref_CDS_2))))
          } 
        }
        return(Ref)
      }
      
      Ref_pop <- Generate_Reference(Ref, consensus_sequence_pop, segment_name_num)
      Ref_intrahost <- Generate_Reference(Ref, consensus_sequence_intrahost, segment_name_num)
      
      # Save sample's reference file as a csv
      Sample_Ref_name_pop <- paste(Strain,'_Ref','_',subject,"_pop",sep=""); Sample_Ref_name_pop
      Sample_Ref_name_intra <- paste(Strain,'_Ref','_',subject,"_intra",sep=""); Sample_Ref_name_intra
      
      
      write.csv(x=Ref_pop,file=paste(new_reference_pop_path,'/',Sample_Ref_name_pop,'.csv',sep=""))
      write.csv(x=Ref_intrahost,file=paste(new_reference_intra_path,'/',Sample_Ref_name_intra,'.csv',sep=""))
      
      
      # Evaluate mutations in this sample ---------------------------------------
      
      # Pull mutation info
      stmp_path_consensus <- list.files(path=consensus_variant_path,pattern = samples[smp]); stmp_path_consensus
      stmp_path_intrahost <- list.files(path=intrahost_variant_path,pattern = samples[smp]); stmp_path_intrahost
      
      stmp_consensus <- read.csv(file.path(consensus_variant_path,stmp_path_consensus))
      original_colclasses <- lapply(stmp_consensus,class)
      original_colclasses$ALT <- "character"
      stmp_consensus <- read.csv(file.path(consensus_variant_path,stmp_path_consensus),colClasses = original_colclasses)
      
      stmp_intra <- read.csv(file.path(intrahost_variant_path,stmp_path_intrahost)); stmp_intra 
      original_colclasses <- lapply(stmp_intra,class)
      original_colclasses$ALT <- "character"
      stmp_intra <- read.csv(file.path(intrahost_variant_path,stmp_path_intrahost),colClasses = original_colclasses)
      
      # Characterize the variants at both the intrahost and consensus level, then concatenate the output  
      mutation_output_pop <- Characterize_mutations(stmp_consensus,subject, Strain, Ref_pop)
      mutation_output_pop$variant_level <- "Pop"
      
      mutation_output_intrahost <- Characterize_mutations(stmp_intra,subject, Strain, Ref_intrahost)
      
      if (is_empty(mutation_output_intrahost) == FALSE){
        mutation_output_intrahost$variant_level <- "Intra"
        mutation_output_tmp <- rbind(mutation_output_pop,mutation_output_intrahost)
      } else {
        mutation_output_tmp <- mutation_output_pop
      }
      mutation_output_tmp <- mutation_output_tmp[
        with(mutation_output_tmp, order(segment_num,position)),
      ]
      mutation_output <- rbind(mutation_output, mutation_output_tmp)
    }
    
    ### Save mutation output
    
    mutation_ouptut_name <-
      paste("MutationOutput", '_', Pipeline_date, '.csv', sep = "")
    mutation_ouptut_name
    write.csv(mutation_output, mutation_ouptut_name)

}
