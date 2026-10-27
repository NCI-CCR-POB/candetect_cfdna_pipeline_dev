#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
healthy_med = args[1]
frag_file = args[2]
sampled_file = args[3]

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
})

# Function to log messages
log_message <- function(msg) {
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg))
}

tryCatch({
  # Log start of processing
  log_message("Starting GC sampling process")
  
  # Validate input files exist
  if (!file.exists(healthy_med)) {
    stop(paste("Healthy median file not found:", healthy_med))
  }
  if (!file.exists(frag_file)) {
    stop(paste("Fragment file not found:", frag_file))
  }
  
  # Load and validate healthy median data
  log_message("Loading healthy median data")
  healthy_fract <- readRDS(healthy_med)
  if (!all(c("gc_strata", "med_frag_fract") %in% names(healthy_fract))) {
    stop("Healthy median data does not contain required columns: gc_strata, med_frag_fract")
  }
  
  # Load and validate fragment data
  log_message("Loading fragment data")
  frag_file_data <- read.table(frag_file, sep = '\t', header = FALSE)
  if (ncol(frag_file_data) < 5) {
    stop("Fragment file must have at least 5 columns: chr, start, end, gc_raw, len")
  }
  
  # Process fragment data
  log_message("Processing fragment data")
  frag_bed <- frag_file_data
  names(frag_bed) <- c("chr", "start", "end", "gc_raw", "len")
  
  # Validate numeric columns
  if (!all(sapply(frag_bed[c("gc_raw", "len")], is.numeric))) {
    stop("gc_raw and len columns must be numeric")
  }
  
  # Process fragments
  frag <- frag_bed %>%
    mutate(gc_strata = round(gc_raw, 2)) %>%
    left_join(healthy_fract, by = "gc_strata")
  
  # Log data summary
  log_message(paste("Total fragments:", nrow(frag)))
  log_message(paste("Unique GC strata:", length(unique(frag$gc_strata))))
  
  # Check for missing values after join
  if (any(is.na(frag$med_frag_fract))) {
    log_message("Warning: Some GC strata do not have corresponding healthy median values")
  }
  
  # Determine sampling parameters
  log_message("Calculating sampling parameters")
  stratatotake <- frag$gc_strata[which.max(frag$med_frag_fract)]
  fragsinmaxstrata <- length(which(frag$gc_strata == stratatotake))
  fragstotake <- round(fragsinmaxstrata/stratatotake)
  
  log_message(paste("Sampling from GC strata:", stratatotake))
  log_message(paste("Fragments in max strata:", fragsinmaxstrata))
  log_message(paste("Target fragments to take:", fragstotake))
  
  # Perform sampling
  log_message("Performing sampling")
  sampled <- frag %>%
    filter(!is.na(med_frag_fract)) %>%
    slice_sample(., n = nrow(.), weight_by = med_frag_fract, replace = TRUE) %>%
    select(chr, start, end, len, gc_strata)
  
  # Validate sampled data
  if (nrow(sampled) == 0) {
    stop("Sampling resulted in zero fragments")
  }
  
  # Write output
  log_message("Writing sampled data")
  write.table(sampled, 
              file = sampled_file,
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  
  log_message("GC sampling completed successfully")
  
}, error = function(e) {
  log_message(paste("Error:", e$message))
  quit(status = 1)
})
