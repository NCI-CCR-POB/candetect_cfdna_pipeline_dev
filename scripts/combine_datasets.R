#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(glue)
})

# Parse command line arguments
option_list <- list(
  make_option("--config", type="character", help="Path to samples config file"),
  make_option("--data_dir", type="character", help="Data directory path"),
  make_option("--cohort1", type="character", help="First cohort"),
  make_option("--cohort2", type="character", help="Second cohort"),
  make_option("--normal_control", type="character", help="Normal control cohort"),
  make_option("--output", type="character", help="Output file path"),
  make_option("--status", type="character", help="Status file path"),
  make_option("--log", type="character", help="Log file path"),
  make_option("--tf", type="character", help="Tumor fraction file"),
  make_option("--fragle", type="character", help="FRAGLE data file"),
  make_option("--mds", type="character", help="Motif diversity scores file"),
  make_option("--fprofile_auc", type="character", help="F-profile AUC file"),
  make_option("--fprofiles", type="character", help="F-profiles file"),
  make_option("--nmf", type="character", help="NMF AUCs file"),
  make_option("--binwise", type="character", help="Bin-wise probabilities file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Initialize logging
log_message <- function(msg) {
  write(msg, file=opt$log, append=TRUE)
}

log_message("Starting combine_datasets.R script")
log_message(glue("Parameters: {paste(names(opt), opt, sep='=', collapse=', ')}"))

# Initialize status tracking
status <- list(
  tf = FALSE,
  fragle = FALSE,
  mds = FALSE,
  fprofile = FALSE,
  nmf = FALSE,
  binwise = FALSE
)
log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ...)
}
# Function to read and merge data
read_and_merge_data <- function(file_path, id_col, value_cols, id_col_combined="library") {
  tryCatch({
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Clean up ID column
    df[[id_col]] <- case_when(
      id_col == "Sample_ID" ~ str_trim(str_to_upper(df[[id_col]])),
      id_col == "library" ~ str_to_upper(str_replace(df[[id_col]], "_motifs", "")),
      id_col == "sample_id" ~ str_to_upper(df[[id_col]]),
      TRUE ~ df[[id_col]]
    )
    
    df %>% select(all_of(c(id_col, value_cols)))
  }, error = function(e) {
    log_message(glue("Error reading {file_path}: {e$message}"))
    tibble(!!id_col := character(), !!!setNames(rep(list(numeric()), length(value_cols)), value_cols))
  })
}


# Read samples config
log_message("Reading samples config file")
samples_df <- read_tsv(opt$config, show_col_types = FALSE) %>%
  mutate(library = str_to_upper(library))

# Initialize combined dataframe
combined_df <- samples_df %>% select(library, cohort)
log_message(glue("Initialized combined dataframe with {nrow(combined_df)} samples"))

# 1. Merge Tumor Fraction data
if (!is.null(opt$tf) && file.exists(opt$tf)) {
  log_message("Processing Tumor Fraction data")
  tf_df <- read_and_merge_data(opt$tf, "Sample_ID", c("Tumor_Fraction"))
  combined_df <- combined_df %>%
    left_join(tf_df, by = c("library" = "Sample_ID"))
  status$tf <- TRUE
} else {
  combined_df$Tumor_Fraction <- "missing"
  log_message("Warning: Tumor Fraction file not found")
}

# 2. Merge FRAGLE data
if (!is.null(opt$fragle) && file.exists(opt$fragle)) {
  log_message("Processing FRAGLE data")
  fragle_df <- read_and_merge_data(opt$fragle, "sample_id", c("lt_ctdna_burden"))
  combined_df <- combined_df %>%
    left_join(fragle_df, by = c("library" = "sample_id"))
  status$fragle <- TRUE
} else {
  combined_df$lt_ctdna_burden <- "missing"
  log_message("Warning: FRAGLE file not found")
}

# 3. Merge motif diversity scores
if (!is.null(opt$mds) && file.exists(opt$mds)) {
  log_message("Processing motif diversity scores")
  mds_df <- read_and_merge_data(opt$mds, "library", c("mds"))
  combined_df <- combined_df %>%
    left_join(mds_df, by = "library")
  status$mds <- TRUE
} else {
  combined_df$mds <- "missing"
  log_message("Warning: Motif diversity scores file not found")
}

# 4. Get best f-profiles and merge
if (!is.null(opt$fprofile_auc) && !is.null(opt$fprofiles) && 
    file.exists(opt$fprofile_auc) && file.exists(opt$fprofiles)) {
  log_message("Processing f-profiles")
  
  # Read and get top 3 F-profiles by AUC
  fprofile_auc <- read_csv(opt$fprofile_auc, show_col_types = FALSE)
  top_3_profiles <- fprofile_auc %>%
    arrange(desc(AUC)) %>%
    slice_head(n = 3)
  
  log_message(glue("Top 3 F-profiles by AUC: {paste(top_3_profiles$fprofile, collapse=', ')}"))
  
  # Read F-profiles data
  fprofiles_df <- read_tsv(opt$fprofiles, show_col_types = FALSE)
  
  # First column should be the profile names
  profile_col <- names(fprofiles_df)[1]
  
  # Filter for top 3 profiles and transpose
  top_profiles_df <- fprofiles_df %>%
    filter(.data[[profile_col]] %in% top_3_profiles$fprofile) %>%
    column_to_rownames(profile_col) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("library") %>%
    mutate(library = str_to_upper(str_replace(library, "_motifs", "")))
  
  # Rename columns to match profile names
  names(top_profiles_df)[-1] <- top_3_profiles$fprofile
  
  # Join with combined dataset
  combined_df <- combined_df %>%
    left_join(top_profiles_df, by = "library")
  
  status$fprofile <- TRUE
  
  log_message(glue("Added F-profile columns: {paste(top_3_profiles$fprofile, collapse=', ')}"))
} else {
  log_message("Warning: F-profile files not found")
  # Don't add placeholder columns since we don't know the profile names
}

# 5. Get best NMF components and merge
if (!is.null(opt$nmf) && file.exists(opt$nmf)) {
  log_message("Processing NMF data")
  aucs_df <- read_csv(opt$nmf, show_col_types = FALSE)
  
  best_n_components <- aucs_df %>%
    filter(cohort1 == opt$cohort1, cohort2 == opt$cohort2) %>%
    arrange(desc(med_score)) %>%
    slice_head(n = 1) %>%
    pull(n_components)
  
  log_message(glue("Best number of components for {opt$cohort1} vs {opt$cohort2}: {best_n_components}"))
  
  nmf_scores_files <- list.files(
    path = file.path(opt$data_dir, "analysis/cfdna_wgs/nmf/scores"),
    pattern = glue("{opt$cohort1}_(v|vs)_{opt$cohort2}\\.n={best_n_components}\\.csv"),
    full.names = TRUE
  )
  
  if (length(nmf_scores_files) > 0) {
    log_message(glue("Using NMF scores file: {nmf_scores_files[1]}"))
    nmf_scores_df <- read_and_merge_data(nmf_scores_files[1], "library", c("med_score"))
    combined_df <- combined_df %>%
      left_join(nmf_scores_df, by = "library") %>%
      rename(nmf_prob = med_score)
    status$nmf <- TRUE
  } else {
    log_message(glue("Warning: No NMF scores file found for n_components={best_n_components}"))
    combined_df$nmf <- "missing"
  }
} else {
  log_message("Warning: NMF AUCs file not found")
  combined_df$nmf <- "missing"
}

# 6. Merge bin-wise fragmentomics data
if (!is.null(opt$binwise) && file.exists(opt$binwise)) {
  log_message("Processing bin-wise fragmentomics data")
  
  # Read binwise data
  binwise_df <- read_csv(opt$binwise, show_col_types = FALSE)
  
  # Calculate mean probability across all repeats
  binwise_df <- binwise_df %>%
    pivot_longer(cols = starts_with("repeat_"), names_to = "repeat", values_to = "prob") %>%
    group_by(library) %>%
    summarize(binwise_prob = median(prob)) %>%
    mutate(library = str_to_upper(library))  # Ensure consistent case
  
  # Join with combined dataset
  combined_df <- combined_df %>%
    left_join(binwise_df, by = "library")
  
  status$binwise <- TRUE
  log_message(glue("Added binwise probabilities (mean across {n_distinct(binwise_df$library)} samples)"))
} else {
  log_message(glue("Warning: Bin-wise probabilities file not found for {opt$cohort1} vs {opt$cohort2}"))
  combined_df$binwise_prob <- "missing"
}

# Replace NA with "missing"
combined_df <- combined_df %>%
  mutate(across(everything(), ~if_else(is.na(.), "missing", as.character(.))))

# Save status information
status_text <- c(
  "Dataset combination status:",
  paste(names(status), 
        ifelse(unlist(status), "Included", "Missing"), 
        sep = ": ")
)
write_lines(status_text, opt$status)

# Save combined dataset
write_csv(combined_df, opt$output)
log_message("Finished combining datasets") 