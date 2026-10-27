#!/usr/bin/env Rscript

#########1#########2#########3#########4#########5#########6#########7#########8
# Mathios, et al. 2021. https://doi.org/10.1038/s41467-021-24994-w
# "In addition to fragmentation profiles, we also computed z-scores for
#  chromosomal arms and a genome-wide summary of the overall cfDNA fragmentation
#  as previously described23,51 with the modifications indicated below. We only
#  analyzed the 39 chromosomal arms that were not acrocentric. Z-scores for each
#  of the 39 autosomal arms were obtained by centering and scaling the total
#  GC-adjusted fragment count for each arm by the mean and standard deviation of
#  the corresponding arm-specific counts in the 54 non-cancer samples used as a
#  reference set."

# Load required packages and data 
library(tidyverse)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cytoband_file <- args[1]
sample_list_file <- args[2]
control_term <- args[3]
frag_counts_file <- args[4]
output_file <- args[5]

existing <- if (length(args) >= 6 && !is.na(args[6]) && nzchar(args[6])) args[6] else "no"
has_existing <- existing != "no"

# Read input files
cytoband = read_tsv(cytoband_file, 
                    col_names = c("chr","start","end","band","stain"),
                    col_types = cols(
                      chr = col_character(),
                      start = col_double(),
                      end = col_double(),
                      band = col_character(),
                      stain = col_character()
                    ))

sample_list <- read_tsv(sample_list_file)
frag_counts <- read_tsv(frag_counts_file) %>%
  filter(library %in% sample_list$library)
length(unique(frag_counts$library))

# Data cleaning and validation
# Remove any rows with NA/NaN/Inf values
frag_counts_clean <- frag_counts %>%
  filter(!is.na(count) & !is.infinite(count) & !is.nan(count))

# Check if we lost any significant data
if(nrow(frag_counts_clean) < nrow(frag_counts)) {
  warning(paste("Removed", nrow(frag_counts) - nrow(frag_counts_clean), 
                "rows with NA/NaN/Inf values"))
}

# Continue with the rest of the analysis using frag_counts_clean

# Get healthy control libraries
healthy_libs =
  sample_list %>%
  filter(cohort == control_term) %>%
  pull(library)
healthy_libs

# Define non-acrocentric arms and their genomic coordinates
(autosomes = c(paste0("chr", seq(1:22))))
(noacro_chroms = c(paste0("chr", seq(1:12)), paste0("chr", seq(16, 20, by =1))))
(acro_chroms = autosomes[!autosomes %in% noacro_chroms])

noacro = cytoband %>%
  mutate(arm = substr(band, 1, 1)) %>%
  group_by(chr, arm) %>%
  summarize(armstart = min(start),
            armend = max(end)) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  pivot_wider(names_from = arm, values_from = c(armstart,armend)) 

# LOESS correction function
gcCorrectLoess <- function(counts.mult, gc) {
  lower <- 0
  upper <- quantile(counts.mult, .99, na.rm=TRUE)
  trend.points <- counts.mult > lower & counts.mult < upper
  trend.counts <- counts.mult[trend.points]
  trend.gc <- gc[trend.points]
  num.sample.points <- min(length(trend.counts), 10E3L)
  samp <- sample(1:length(trend.counts), num.sample.points)
  include <- c(which(gc==min(gc)), which(gc==max(gc)))
  trend.counts <- c(counts.mult[include], trend.counts[samp])
  trend.gc <- c(gc[include], trend.gc[samp])
  initial.trend <- loess(trend.counts ~ trend.gc, span=0.75)
  i <- seq(min(gc, na.rm=TRUE), max(gc, na.rm=TRUE), by = 0.001)
  final.trend <- loess(predict(initial.trend, i) ~ i, span=0.75)
  counts.pred <- predict(final.trend, gc)
  return(counts.pred)
}

# Use non-acrocentric arm positions to count fragments per arm
arm_counts =
  frag_counts_clean %>% 
  # Get total counts, short and long, per position
  group_by(library,chr,start,end,gc) %>% 
  summarize(cov=sum(count)) %>%
  # Annotate arms and filter to non-acro arms
  filter(chr %in% autosomes) %>% 
  left_join(noacro, by = "chr") %>%
  mutate(arm = ifelse(end < armstart_q, "p", "q")) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  subset(select=c(library,chr,start,end,gc,arm,cov))

# Corrections over bins
arm_counts_DT <- setDT(arm_counts)
arm_counts_DT[,counts.mult:=log2((get("cov") + 1))]
arm_counts_DT[,loess.pred:=gcCorrectLoess(counts.mult, gc),by="library"]
arm_counts_DT[,adjusted:=counts.mult - loess.pred]
arm_counts_DT[,adjusted.cent:=adjusted - median(adjusted, na.rm=TRUE), by="library"]

# Aggregate bins to arms
arms <- arm_counts_DT[,.(armmean=mean(adjusted.cent, na.rm=TRUE)),
             by=.(library, chr, arm)] 
arms <- arms[,armmean:=armmean - mean(armmean, na.rm = TRUE)]

arms_tib <- as_tibble(arms)

# Get mean and sd of healthy libs
if (!has_existing) {
  print("no existing healthy mean and sd")
  healthy_arms =
    arms_tib %>%
    filter(library %in% healthy_libs) %>%
    group_by(chr,arm) %>%
    summarize(mean = mean(armmean, na.rm = T), sd = sd(armmean, na.rm=T)) 
  
  output_file1 = gsub(basename(output_file),"binwise_armz_healthy_mean_sd.csv",output_file)
  write.csv(healthy_arms, file = output_file1, row.names=F,quote = F)
} else {
  print("using existing healthy mean and sd")
  healthy_arms = read.csv(existing,header=T)
}
head(healthy_arms)

# z-score against the healthy cohort
arm_z =  
  arms_tib %>%
  left_join(healthy_arms, by = c("chr","arm")) %>%
  mutate(z = (armmean-mean)/sd) %>%
  select(!c("mean", "sd"))

arm_z2 =
  arm_z %>%
  select(-armmean) %>%
  unite(chr_arm,chr:arm,sep="_") %>%
  spread(chr_arm,z) %>%
  column_to_rownames('library')

# Add data validation before loess smoothing
if(any(is.na(frag_counts_clean$count)) || 
   any(is.infinite(frag_counts_clean$count)) || 
   any(is.nan(frag_counts_clean$count))) {
  stop("Data still contains NA/NaN/Inf values after cleaning")
}

write.csv(arm_z2, file = output_file, row.names=T)
