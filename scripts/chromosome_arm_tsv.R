#!/usr/bin/env Rscript

###########################
###   Define Chr Arms   ###
###########################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
cytobands_tsv = args[1]
arms_tsv = args[2]

#To run stand-alone, update paths below
#cytobands_tsv = "/Volumes/ShernLiquidBx/ref/cytoBand.txt"
#arms_tsv = "/Volumes/ShernLiquidBx/ref/hg38_chr_arms.tsv"

# Load required packages, data, and functions
library(tidyverse)
cytobands = read_tsv(cytobands_tsv, col_names = c("chr","start","end","band","stain"))

# Define non-acrocentric arms and their genomic coordinates
(autosomes = c(paste0("chr", seq(1:22))))
(noacro_chroms = c(paste0("chr", seq(1:12)), paste0("chr", seq(16, 20, by =1))))
(acro_chroms = autosomes[!autosomes %in% noacro_chroms])

noacro =
  cytobands %>%
  mutate(arm = substr(band, 1, 1)) %>%
  group_by(chr, arm) %>%
  summarize(armstart = min(start),
            armend = max(end)) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q"))
noacro %>% write_tsv(arms_tsv)
