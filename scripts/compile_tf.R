#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compile_tf.R <ichor_dir> <output_file>")
}

ichor_dir <- args[1]
output_file <- args[2]

# Required packages
required_packages <- c("dplyr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Find all params files
param_files <- list.files(ichor_dir, pattern = "\\.params.txt$", recursive = TRUE, full.names = TRUE)
cat(sprintf("Found %d parameter files\n", length(param_files)))

# Function to extract sample ID from path
extract_sample_id <- function(path) {
  # If file is in a subdirectory, use subdirectory name
  if (dirname(path) != ichor_dir) {
    return(basename(dirname(path)))
  }
  # Otherwise use filename without extension
  return(sub("\\.params\\.txt$", "", basename(path)))
}

# Read and process each params file
results <- lapply(param_files, function(file) {
  cat(sprintf("Processing: %s\n", file))
  
  # Read first line of file (header and data)
  df <- tryCatch({
    read.table(file, header = TRUE, nrows = 1, sep = "\t")
  }, error = function(e) {
    cat(sprintf("Error reading %s: %s\n", file, e$message))
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  # Extract sample ID
  sample_id <- extract_sample_id(file)
  
  # Create result row
  data.frame(
    Sample_ID = sample_id,
    Tumor_Fraction = df$Tumor.Fraction,
    Ploidy = df$Ploidy,
    stringsAsFactors = FALSE
  )
})

# Combine all results
results_df <- do.call(rbind, results[!sapply(results, is.null)])

# Print summary
cat("\nSummary of results:\n")
print(results_df)

# Save to CSV
write.csv(results_df, output_file, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", output_file)) 