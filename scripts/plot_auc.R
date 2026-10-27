#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(glue)
  library(pROC)
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
#  cohort2 = "MPNST",
#  normal_control = "healthy",
#  output = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/plot_auc_comparison_healthy_v_MPNST.png",
#  log = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/logs/combine_datasets/plot_auc_comparison_healthy_v_MPNST.log"
#)

# Initialize logging
log_message <- function(msg) {
  write(msg, file=opt$log, append=TRUE)
}

log_message("Starting plot_auc.R script")
log_message(glue("Parameters: {paste(names(opt), opt, sep='=', collapse=', ')}"))

# Read combined dataset
log_message("Reading combined dataset")
combined_df <- read_csv(opt$input, show_col_types = FALSE)

# Filter for relevant cohorts
combined_df <- combined_df %>%
  filter(cohort %in% c(opt$cohort1, opt$cohort2))

# Function to calculate ROC and AUC
calculate_roc_auc <- function(data, metric, cohort1, cohort2) {
  # Convert metric values from character to numeric, handling "missing" values
  data[[metric]] <- suppressWarnings(as.numeric(data[[metric]]))
  
  # Remove rows with NA values
  data <- data %>% filter(!is.na(.data[[metric]]))
  
  # Calculate ROC curve with higher scores indicating higher probability of cohort2.
  tryCatch({
    roc_obj <- roc(
      response = data$cohort,
      predictor = data[[metric]],
      levels = c(cohort1, cohort2),
      direction = "<",
      quiet = TRUE
    )
    list(
      auc = auc(roc_obj),
      roc = roc_obj
    )
  }, error = function(e) {
    log_message(glue("Error calculating ROC for {metric}: {e$message}"))
    NULL
  })
}

# Identify all numeric columns that could be metrics
base_metrics <- c(
  "Tumor_Fraction" = "CNA Tumor Fraction",
  "lt_ctdna_burden" = "FRAGLE Tumor Fraction",
  "mds" = "Motif Diversity Score",
  "nmf_prob" = "NMF Probability",
  "binwise_prob" = "Bin-wise Probability"
)

# Add F-profile metrics if present
fprofile_cols <- names(combined_df)[str_detect(names(combined_df), "^fprof")]
if (length(fprofile_cols) > 0) {
  fprofile_metrics <- setNames(
    paste("F-profile", str_extract(fprofile_cols, "\\d+")),
    fprofile_cols
  )
  metrics <- c(base_metrics, fprofile_metrics)
  log_message(glue("Found F-profile columns: {paste(fprofile_cols, collapse=', ')}"))
} else {
  metrics <- base_metrics
  log_message("No F-profile columns found in dataset")
}

log_message(glue("Analyzing metrics: {paste(names(metrics), collapse=', ')}"))

# Calculate ROC curves and AUCs for each metric
roc_results <- list()
auc_data <- tibble(
  metric = character(),
  metric_label = character(),
  auc = double()
)

for (metric_name in names(metrics)) {
  log_message(glue("Calculating ROC/AUC for {metric_name}"))
  
  result <- calculate_roc_auc(combined_df, metric_name, opt$cohort1, opt$cohort2)
  
  if (!is.null(result)) {
    roc_results[[metric_name]] <- result$roc
    auc_data <- bind_rows(auc_data, tibble(
      metric = metric_name,
      metric_label = metrics[[metric_name]],
      auc = as.double(result$auc)
    ))
  }
}

write_csv(
  auc_data,
  file = gsub("\\.png$", ".csv", opt$output),
  quote = "none"
)

roc_output <- ifelse(
  str_detect(basename(opt$output), "^auc_comparison_"),
  file.path(dirname(opt$output), str_replace(basename(opt$output), "^auc_comparison_", "roc_comparison_")),
  gsub("\\.png$", "_roc.png", opt$output)
)

if (length(roc_results) > 0 && nrow(auc_data) > 0) {
  roc_plot_data <- purrr::map_dfr(names(roc_results), function(metric_name) {
    roc_obj <- roc_results[[metric_name]]
    label <- auc_data$metric_label[auc_data$metric == metric_name][1]
    auc_value <- auc_data$auc[auc_data$metric == metric_name][1]
    tibble(
      metric = metric_name,
      metric_label = glue("{label} (AUC={round(auc_value, 3)})"),
      fpr = 1 - as.numeric(roc_obj$specificities),
      tpr = as.numeric(roc_obj$sensitivities),
      point_order = seq_along(roc_obj$sensitivities)
    )
  })

  roc_plot <- ggplot(roc_plot_data, aes(x = fpr, y = tpr, color = metric_label, group = metric_label)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray55") +
    geom_path(size = 0.9) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      title = glue("ROC Curves: {opt$cohort1} vs {opt$cohort2}"),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )

  ggsave(roc_output, roc_plot, width = 10, height = 7)
  log_message(glue("ROC comparison plot saved successfully: {roc_output}"))
} else {
  log_message("Warning: No valid ROC curves for ROC comparison plot")
  roc_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No valid ROC curves for analysis") +
    theme_void()
  ggsave(roc_output, roc_plot, width = 10, height = 7)
}

# Create AUC comparison plot
if (nrow(auc_data) > 0) {
  # Order metrics by AUC value
  plot_data <- auc_data %>%
    arrange(desc(auc)) %>%
    mutate(metric_label = factor(metric_label, levels = metric_label))
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = metric_label, y = auc)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = round(auc, 3)), vjust = -0.5) +
    labs(
      title = glue("AUC Comparison: {opt$cohort1} vs {opt$cohort2}"),
      x = "Metric",
      y = "AUC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    ylim(0, max(plot_data$auc) * 1.1) # Add some space for labels

  # Save plot
  ggsave(opt$output, p, width = 10, height = 6)
  log_message("AUC comparison plot saved successfully")
} else {
  log_message("Warning: No valid metrics for AUC analysis")
  # Create empty plot with message
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No valid metrics for AUC analysis") +
    theme_void()
  ggsave(opt$output, p, width = 10, height = 6)
}

log_message("Finished plotting AUC comparison") 
