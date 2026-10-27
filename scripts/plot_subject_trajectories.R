#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(glue)
})

# Parse command line arguments
option_list <- list(
  make_option("--input", type="character", help="Input combined dataset file"),
  make_option("--output", type="character", help="Output plot file"),
  make_option("--cohort1", type="character", help="First cohort"),
  make_option("--cohort2", type="character", help="Second cohort"),
  make_option("--normal_control", type="character", help="Normal control cohort"),
  make_option("--log", type="character", help="Log file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

#opt = list(
#  input = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/combined_dataset_health_v_MPNST.csv",
#  data_dir = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new",
#  cohort1 = "healthy",
# # cohort2 = "MPNST",
#  normal_control = "healthy",
#  output = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/subject_trajectories_healthy_v_MPNST.png",
#  log = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/logs/combine_datasets/plot_subject_trajectories_healthy_v_MPNST.log"
# )

# Initialize logging
log_message <- function(msg) {
  write(msg, file=opt$log, append=TRUE)
}

log_message("Starting plot_subject_trajectories.R script")
log_message(glue("Parameters: {paste(names(opt), opt, sep='=', collapse=', ')}"))

# Read combined dataset
log_message("Reading combined dataset")
combined_df <- read_csv(opt$input, show_col_types = FALSE)

# Filter for relevant cohorts
combined_df <- combined_df %>%
  filter(cohort %in% c(opt$cohort1, opt$cohort2))

# Extract subject IDs from library names (assuming format SUBJECT_TIMEPOINT)
combined_df <- combined_df %>%
  mutate(
    subject = str_extract(library, "^[^_]+"),
    timepoint = as.numeric(str_extract(library, "\\d+$"))
  )

# Get list of subjects with multiple timepoints
subjects_with_trajectories <- combined_df %>%
  group_by(subject) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pull(subject) %>%
  unique()

log_message(glue("Found {length(subjects_with_trajectories)} subjects with multiple timepoints"))

if (length(subjects_with_trajectories) == 0) {
  log_message("Warning: No subjects with multiple timepoints found")
  # Create empty plot with message
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No subjects with multiple timepoints found") +
    theme_void()
} else {
  # Prepare data for plotting
  plot_data <- combined_df %>%
    filter(subject %in% subjects_with_trajectories)
  
  # Convert all numeric columns from "missing" to NA
  numeric_cols <- c("Tumor_Fraction", "lt_ctdna_burden", "mds", "nmf", "binwise_prob")
  fprofile_cols <- names(combined_df)[str_detect(names(combined_df), "^fprof")]
  all_numeric_cols <- c(numeric_cols, fprofile_cols)
  
  plot_data <- plot_data %>%
    mutate(across(all_of(all_numeric_cols), function(x) {
      suppressWarnings(as.numeric(ifelse(x == "missing", NA, x)))
    }))
  
  # Create trajectory plots for each metric
  base_metrics <- c("Tumor_Fraction", "lt_ctdna_burden", "mds", "nmf", "binwise_prob")
  base_labels <- c("Tumor Fraction", "FRAGLE ctDNA Burden", "Motif Diversity Score", 
                  "NMF Score", "Bin-wise Probability")
  
  if (length(fprofile_cols) > 0) {
    metrics <- c(base_metrics, fprofile_cols)
    metric_labels <- c(base_labels, 
                      paste("F-profile", str_extract(fprofile_cols, "\\d+")))
    log_message(glue("Including F-profile metrics: {paste(fprofile_cols, collapse=', ')}"))
  } else {
    metrics <- base_metrics
    metric_labels <- base_labels
    log_message("No F-profile metrics found in dataset")
  }
  
  plot_list <- map2(metrics, metric_labels, function(metric, label) {
    ggplot(plot_data, aes(x = timepoint, y = .data[[metric]], group = subject, color = subject)) +
      geom_line() +
      geom_point() +
      labs(x = "Timepoint", y = label) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  # Combine plots into a single figure
  p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
}

# Save plot
log_message("Saving plot")
ggsave(opt$output, p, width = 12, height = 10)
log_message("Finished plotting subject trajectories") 