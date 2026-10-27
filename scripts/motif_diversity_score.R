#!/usr/bin/env Rscript

#Parse arguments
args = commandArgs(trailingOnly = TRUE)
motif_frac_combo_tsv = args[1]
annotation_table = args[2]
MDS_R_Out = args[3]
CSVout = args[4]

# Load necessary packages
library(tidyverse)
library(dplyr)

# Load data frames
annotation <- read_tsv(annotation_table) %>% 
  select(library, cohort) %>%
  mutate(library = as.character(library))  # Convert library to character type
annotation

motifs <- read_tsv(motif_frac_combo_tsv)
motifs

#Pivot data frame
motifs_long <- pivot_longer(motifs, cols = !motif, names_to = "library", values_to = "fraction") %>%
    mutate(library = as.character(library)) %>%  # Convert library to character type
    left_join(annotation, by = "library") %>%
    select(motif, library, fraction, cohort)
motifs_long

# Generate motif diversity score (mds) data.frame
mds <- motifs_long %>%
    mutate(fraction = ifelse(fraction == 0, 1e-10, fraction)) %>% # Add a small constant to avoid log(0)
    mutate(mds = -fraction * log(fraction) / log(256)) %>%
    group_by(library) %>%
    summarize(mds = sum(mds))
mds

#Save outputs
save.image(file = MDS_R_Out)
write.csv(mds, CSVout)