#!/usr/bin/env Rscript

# Script to plot overlapping histograms of NMF components

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(readr)
  library(purrr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--model"), type="character", help="Path to the NMF model pickle file"),
  make_option(c("--freq_csv"), type="character", help="Path to the frequency histogram CSV file"),
  make_option(c("--output"), type="character", help="Path to output plot file"),
  make_option(c("--cohort1"), type="character", help="First cohort name"),
  make_option(c("--cohort2"), type="character", help="Second cohort name"),
  make_option(c("--n_components"), type="integer", help="Number of NMF components")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to read pickle file (assuming it's a Python pickle file)
# Note: This requires the reticulate package to read Python pickle files
read_pickle <- function(file_path) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The reticulate package is required to read pickle files")
  }
  reticulate::py_load_object(file_path)
}

# Function to plot NMF components
plot_nmf_components <- function(model_path, freq_csv, output_path, cohort1, cohort2, n_components) {
  # Load the NMF model
  model <- read_pickle(model_path)
  
  # Load the frequency data
  freq_df <- read_csv(freq_csv)
  
  # Get the components from the model
  components <- model$named_steps$nmf$components_
  
  # Create a data frame for plotting
  plot_data <- data.frame()
  
  # Prepare data for each component
  for (i in 1:n_components) {
    component_data <- data.frame(
      fragment_size = as.numeric(colnames(freq_df)),
      weight = components[i,],
      component = paste0("Component ", i)
    )
    plot_data <- rbind(plot_data, component_data)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = fragment_size, y = weight, color = component)) +
    geom_line(size = 1) +
    facet_wrap(~component, ncol = 1) +
    labs(
      title = "NMF Components",
      x = "Fragment Size (bp)",
      y = "Component Weight"
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid = element_line(color = "gray90"),
      legend.position = "none"
    ) +
    scale_color_brewer(palette = "Set1")
  
  # Save the plot
  ggsave(output_path, p, width = 12, height = 4 * n_components, dpi = 300)
}

# Main execution
if (!interactive()) {
  # Create output directory if it doesn't exist
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  
  # Generate the plot
  plot_nmf_components(
    opt$model,
    opt$freq_csv,
    opt$output,
    opt$cohort1,
    opt$cohort2,
    opt$n_components
  )
} 