#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggbeeswarm)
  library(glue)
  library(patchwork)
})

# Parse command line arguments
option_list <- list(
  make_option("--input", type="character", help="Path to combined dataset CSV"),
  make_option("--output_dir", type="character", help="Directory to save plots"),
  make_option("--cohort1", type="character", help="First cohort name"),
  make_option("--cohort2", type="character", help="Second cohort name"),
  make_option("--normal_control", type="character", help="Name of control cohort"),
  make_option("--log", type="character", help="Path to log file")
)

opt <- parse_args(OptionParser(option_list=option_list))

#opt = list(
#  input = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/combined_dataset_health_v_MPNST.csv",
#  output_dir = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs",
#  cohort1 = "healthy",
#  cohort2 = "MPNST",
#  normal_control = "healthy",
#  output = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/plot_auc_comparison_healthy_v_MPNST.png",
#  log = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/logs/combine_datasets/plot_auc_comparison_healthy_v_MPNST.log"
#)

# Initialize logging
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message <- glue("{timestamp} - {msg}")
  write(message, file=opt$log, append=TRUE)
  cat(message, "\n")
}

log_message("Starting plot_metric_distributions.R")
log_message(glue("Parameters: input={opt$input}, output_dir={opt$output_dir}, ",
                 "cohort1={opt$cohort1}, cohort2={opt$cohort2}, ",
                 "normal_control={opt$normal_control}"))

# Read the combined dataset
data <- read_csv(opt$input)
log_message(glue("Read {nrow(data)} rows from {opt$input}"))

# Filter for relevant cohorts
data <- data %>%
  filter(cohort %in% c(opt$cohort1, opt$cohort2))
log_message(glue("Filtered to {nrow(data)} rows for cohorts {opt$cohort1} and {opt$cohort2}"))

# Get all numeric columns (excluding 'library' and 'cohort')
numeric_cols <- names(data)[sapply(data, is.numeric)]
log_message(glue("Found {length(numeric_cols)} numeric columns: {paste(numeric_cols, collapse=', ')}"))

# Function to create distribution plot for a metric
create_distribution_plot <- function(data, metric_col) {
  # Get clean metric name for plot title
  metric_name <- str_replace_all(metric_col, "_", " ") %>%
    str_to_title()
  
  # Perform Wilcoxon test
  test_result <- wilcox.test(
    data %>% filter(cohort == opt$cohort2) %>% pull(!!sym(metric_col)),
    data %>% filter(cohort == opt$cohort1) %>% pull(!!sym(metric_col))
  )
  p_value <- test_result$p.value
  
  # Create plot
  p <- ggplot(data, aes(x = cohort, y = !!sym(metric_col))) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white") +
    geom_quasirandom(width = 0.25, size = 1, alpha = 0.6) +
    labs(
      title = metric_name,
      subtitle = glue("p = {format.pval(p_value, digits=3)}"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# Create plots for each numeric column
plot_list <- map(numeric_cols, ~create_distribution_plot(data, .))
names(plot_list) <- numeric_cols

# Calculate optimal grid dimensions
n_plots <- length(numeric_cols)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

# Combine plots in a grid
combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = glue("{opt$cohort2} vs {opt$cohort1} Distribution Comparison"),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

# Save combined plot
output_file <- file.path(opt$output_dir, "all_distributions.png")
ggsave(output_file, combined_plot,
       width = 3 * n_cols,  # 3 inches per plot width
       height = 3 * n_rows, # 3 inches per plot height
       dpi = 300)
log_message(glue("Saved combined plot to {output_file}"))

log_message("Finished creating distribution plots") 
