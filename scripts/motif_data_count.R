#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
motif_tsv_dir = args[1]
output_counts = args[2]
output_fractions = args[3]

# Load necessary packages
library(tidyverse)

# Function to read a file and format it for merging with counts
read_motif_file_counts <- function(file) {
    df <- read_tsv(file, col_names = c("motif", "count"))
    file_name <- tools::file_path_sans_ext(basename(file))
    df <- df %>% rename(!!file_name := count)
    return(df)
}

# Function to read a file and format it for merging with fractions
read_motif_file_fractions <- function(file) {
    df <- read_tsv(file, col_names = c("motif", "count"))
    total_count <- sum(df$count)
    df <- df %>% mutate(fraction = count / total_count)
    file_name <- tools::file_path_sans_ext(basename(file))
    df <- df %>% select(motif, fraction) %>% rename(!!file_name := fraction)
    return(df)
}

# List of files
files <- list.files(path = motif_tsv_dir, pattern = "*.tsv", full.names = TRUE)

# Read and merge all files for counts
motif_data_counts <- files %>%
    map(read_motif_file_counts) %>%
    reduce(full_join, by = "motif") %>%
    na.omit(.)

# Read and merge all files for fractions
motif_data_fractions <- files %>%
    map(read_motif_file_fractions) %>%
    reduce(full_join, by = "motif") %>%
    na.omit(.)

# Save the resulting matrices
write_tsv(motif_data_counts, output_counts)
write_tsv(motif_data_fractions, output_fractions)