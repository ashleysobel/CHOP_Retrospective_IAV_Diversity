
# Prepare environment -----

# Clear environment to get rid of old files 
rm(list = ls())

# Get path to working directory
R_path <- getwd()
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
Pipeline_OutPut_path <- file.path(R_path,"Pipeline_Output_Files")
Reference_path <- file.path(R_path,"Reference_files")

# Load packages ----
library(ggplot2)
library(stringi)

# Load reference/data files ----
# Load Reference Key (matches PB2 accession # with strain). Make sure this has been updated if your run contains different reference files
Reference_key <- read.csv(file.path(Reference_path,"Reference_key.csv"))

# Set segment names
segment_name <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")


# Load final sample set
sample_set <- read.csv(file.path(Compiled_OutPut_path,"CHOA_SampleSet_FINAL.csv"))
Final_Samples <- sample_set$Sample

# Generate base coverge plots ------------
GetCov_expanded <- function(Pipeline_OutPut_path, Pipeline_date) {
  path_basecov <- file.path(Pipeline_OutPut_path,"BaseCov_Expanded",paste0(as.character(Pipeline_date),"-BaseCov"))
  path_basecov
  base_cov <- list.files(path_basecov, pattern = "basecov.txt"); base_cov
  cov_info <- list(path_basecov=path_basecov,basecov_names=base_cov)
}


# Enter the csv file containing the references you want to use. For the CHOP
# 2017-2018 retrospective study, we want to use H3N2 2022 and H1N1 2019.
H3N2_Ref <- read.csv(file.path(Reference_path,"H3N2_Ref_A_2016.csv"))
H1N1_Ref <- read.csv(file.path(Reference_path,"H1N1_Ref_A_2015.csv"))


RunDate <- c(33122,50422,61722,91222,22323,31023,50123) # These are the run dates of the sample tracker to include
RD <- 1
for (RD in 1:length(RunDate)){
  
  Pipeline_date <- RunDate[RD]
  PipelineOutput_path <- file.path(Pipeline_OutPut_path,paste0(Pipeline_date,"_PipelineOutput")); PipelineOutput_path
  cov_info <- GetCov_expanded(Pipeline_OutPut_path,Pipeline_date)
  cov_info
  
  # Identify elements in cov_info$basecov_names that do not contain "CHOA". 
  # Extract these into a character vector called "Negative_controls"
  Negative_controls <- cov_info$basecov_names[!grepl("CHOA", cov_info$basecov_names)]
  
  if (length(Negative_controls != 0)){
    Negative_control_names <- c()
    for (nc in 1:length(Negative_controls)){
      Negative_control_names[nc] <- paste0("Negative_",as.character(nc))
    }
  }
  
  # 2. Identify elements in cov_info$basecov_names where the first portion of the name (components before "_S") are contained in the Final_Samples vector.
  split_basecov_names <- sapply(strsplit(cov_info$basecov_names, "_S"), `[`, 1)
  
  Samples <- cov_info$basecov_names[split_basecov_names %in% Final_Samples]

  # 3. For the samples contained in "Samples" generate a new vector called "Sample_names" which contains just the first portion of the name (components before the "_S").
  Sample_names <- sapply(strsplit(Samples, "_S"), `[`, 1)
  Sample_names <- gsub(pattern = "CHOA-",replacement = "C-",x = Sample_names)
  
  # Add the negative controls
  if (length(Negative_controls != 0)){
    Samples <- c(Samples,Negative_controls)
    Sample_names <- c(Sample_names,Negative_control_names)
  }
  
  
  # Create a data frame to store all samples data
  all_samples_df <- data.frame()
  
  st <- 1
  # Loop through Samples
  for (st in 1:length(Sample_names)){
    sample_tmp <- Samples[st]
    base_cov_path <- file.path(cov_info$path_basecov,sample_tmp); base_cov_path
    base_cov <- read.delim(base_cov_path)
    colnames(base_cov) <- c("Segment","Position","Coverage")
    sorted_base_cov <- base_cov[order(stringi::stri_extract_last(base_cov$Segment,regex="\\w")),]
    sorted_base_cov$index <- 1:nrow(sorted_base_cov) 

    # Define the breaks for the grouping
    breaks <- seq(1, max(sorted_base_cov$index), by = 100)
    
    # Generate the grouping variable
    sorted_base_cov$group <- cut(sorted_base_cov$index, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    
    # Calculate the average coverage by group
    Avg_base_cov <- aggregate(Coverage ~ group, data = sorted_base_cov, FUN = mean)
    
    # Use the lower limit of the group range as the index
    Avg_base_cov$Index <- breaks[Avg_base_cov$group]
    
    # Rename the Coverage column to AvgCov
    colnames(Avg_base_cov)[colnames(Avg_base_cov) == "Coverage"] <- "AvgCov"
    
    # Replace NA values with 0
    Avg_base_cov[is.na(Avg_base_cov)] <- 0
    
    # Remove the group column
    Avg_base_cov$group <- NULL
    Avg_base_cov
    Avg_base_cov$Index <- Avg_base_cov$Index - 1
    Avg_base_cov$Sample <- Sample_names[st]
  
    all_samples_df <- rbind(all_samples_df, Avg_base_cov)
  }

  # create a new column which checks if Sample contains "Negative"
  all_samples_df$is_negative <- grepl("Negative", all_samples_df$Sample)
  
  CovPlot <- ggplot(all_samples_df, aes(x=Index, y=log10(AvgCov+1))) + 
    geom_line(aes(colour = ifelse(is_negative, "Negative", Sample))) +
    scale_color_manual(values = c("black", rainbow(length(unique(all_samples_df$Sample[all_samples_df$is_negative == FALSE])))), 
                       name = "Sample", 
                       breaks = c("Negative", unique(all_samples_df$Sample[all_samples_df$is_negative == FALSE]))) +
    scale_y_continuous(breaks = 0:5, labels = 10^(0:5)) + # specify breaks and labels
    labs(x="Genome Position", y="Coverage") +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.position = "none") + # hide the legend
    scale_x_continuous(breaks = seq(0, 15000, by = 5000)) + # specify x-axis breaks every 5000
    coord_cartesian(ylim = c(0, 5), xlim = c(0, 15000)) # set the limits of x and y axes
  print(CovPlot)
  
  
  
  Coverage_file_name <- paste0(as.character(Pipeline_date),"_SampleDepth.eps")
  ggplot2::ggsave(file.path(Figure_OutPut_path,Coverage_file_name), plot = CovPlot, width = 7, height = 8, device = "eps")
  
}
