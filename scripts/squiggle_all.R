#!/usr/bin/env Rscript

####################################
###   Genome Fragment Profiles   ###
####################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
ratios_tsv = args[1] #"~/mpnst/analysis/cfdna-wgs/frag/ratios.tsv"
train_tsv = args[2] #"~/mpnst/analysis/cfdna-wgs/lists/cfdna_wgs_training.tsv"
libs_rds = args[3] #"~/mpnst/data-model/libraries_full.rds"
arms_tsv = args[4] #"~/mpnst/ref/hg38_chr_arms.tsv"
#bins_tsv = args[2]#~/mpnst/ref/keep_5mb.bed"

setwd("/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs")
ratios_tsv = "/data/CCRCCDI/analysis/ccrtegs65/nf1_fragmentomics_data_11-3-23/inputs/ratios.tsv" #"~/mpnst/analysis/cfdna-wgs/frag/ratios.tsv"
train_tsv = "/data/CCRCCDI/analysis/ccrtegs65/nf1_fragmentomics_data_11-3-23/inputs/cfdna_wgs_training.fixed.tsv" #"~/mpnst/analysis/cfdna-wgs/lists/cfdna_wgs_training.tsv"
libs_rds =  "/data/CCRCCDI/analysis/ccrtegs65/nf1_fragmentomics_data_11-3-23/inputs/libraries_full.rds" #"~/mpnst/data-model/libraries_full.rds"
arms_tsv = "/data/ShernLiquidBx/ref/hg38_chr_arms.tsv" #"~/mpnst/ref/hg38_chr_arms.tsv"
bins_tsv = "/data/ShernLiquidBx/ref/keep_5mb.bed" #args[2]#~/mpnst/ref/keep_5mb.bed"

# Load required packages, data, and functions
library(tidyverse)
ratios = read_tsv(ratios_tsv)
head(ratios)
bins = read_tsv(bins_tsv, col_names = c("chr","start","end","gc")) %>% mutate(bin = paste0(chr, "_", start, "_", end))
arms = read_tsv(arms_tsv)
arms
train = read_tsv(train_tsv)
libraries_full = readRDS(libs_rds)
colnames(train)
libs =
  train %>%
  left_join(libraries_full,by="library") %>%
  mutate(dx_simple = ifelse(dx_simple == "mpnst" & washout == "y", "pretx",
                     ifelse(dx_simple == "mpnst" & washout == "n", "ontx", dx_simple))) %>%
  select(library, dx_simple) %>%
  mutate(dx_simple = factor(dx_simple, levels = c("healthy", "plexiform", "anf", "pretx", "ontx")))
libs
colnames(arms )
arms_wider <- arms %>%
  pivot_wider(
    names_from = arm,
    values_from = c(armstart, armend),
    names_sep = "_"
  )
arms_wider

merged = bins %>% right_join(arms_wider, by = "chr") %>%
  mutate(bin_arm = case_when(
           start < armstart_q ~ "p",
           TRUE ~ "q"
         )) %>%
  mutate(chr_arm = paste0(chr, "_", bin_arm))

# define the order of levels
level_order <- paste0("chr", rep(1:22, each = 2), "_", rep(c("p", "q"), times = 22))

# convert chr_arm to an ordered factor
merged <- merged %>%
  mutate(chr_arm = factor(chr_arm, levels = level_order, ordered = TRUE))

# Remove specified levels
merged <- merged %>%
  filter(!(chr_arm %in% c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p")))

new_labels <- c("plexiform" = "Plexiform", "healthy" = "Healthy",
                "anf" = "ANF",
                "pretx" = "Pre-Tx MPNST",
                "ontx" = "On-Tx MPNST",
                "OTHER" = "OTHER",
                "chr1_p" = "1p",
                "chr1_q" = "1q",
                "chr2_p" = "2p",
                "chr2_q" = "2q",
                "chr3_p" = "3p",
                "chr3_q" = "3q",
                "chr4_p" = "4p",
                "chr4_q" = "4q",
                "chr5_p" = "5p",
                "chr5_q" = "5q",
                "chr6_p" = "6p",
                "chr6_q" = "6q",
                "chr7_p" = "7p",
                "chr7_q" = "7q",
                "chr8_p" = "8p",
                "chr8_q" = "8q",
                "chr9_p" = "9p",
                "chr9_q" = "9q",
                "chr10_p" = "10p",
                "chr10_q" = "11q",
                "chr11_p" = "11p",
                "chr11_q" = "11q",
                "chr12_p" = "12p",
                "chr12_q" = "12q",
                "chr13_p" = "13p",
                "chr13_q" = "13q",
                "chr14_p" = "14p",
                "chr14_q" = "14q",
                "chr15_p" = "15p",
                "chr15_q" = "15q",
                "chr16_p" = "16p",
                "chr16_q" = "16q",
                "chr17_p" = "17p",
                "chr17_q" = "17q",
                "chr18_p" = "18p",
                "chr18_q" = "18q",
                "chr19_p" = "19p",
                "chr19_q" = "19q",
                "chr20_p" = "20p",
                "chr20_q" = "20q",
                "chr21_p" = "21p",
                "chr21_q" = "21q",
                "chr22_p" = "22p",
                "chr22_q" = "22q")
colnames(ratios)[1] = "library"
#libs2 = libs %>%
#  mutate(dx_simple = as.character(dx_simple),
#         dx_simple = ifelse(dx_simple == "pretx", dx_simple,"OTHER"))
plot =
  ratios %>%
  right_join(libs, by = "library") %>%
  left_join(merged, by = c("chr","start","end")) %>%
  filter(!is.na((chr_arm))) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line(color = "#4C4E52") +
  coord_cartesian(ylim = c(-0.1, .1)) +
  scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  facet_grid(rows = vars(dx_simple), cols = vars(chr_arm), scales = "free_x", switch = "x", labeller = as_labeller(new_labels))
ggsave("test_ori.pdf", plot, width = 20, height = 8)

plot
library(cowplot)

# Prepare plot
p <- plot + theme(strip.text.y = element_blank())

# Save the plot

pdf("squiggle_ori.pdf", width = 15, height = 8)
plot
dev.off()

library(cowplot)

# Assume 'p' is your plot
p =
  ggdraw(p) +
  draw_text(text = "Healthy", x = .2, y = .95, size = 18, color = "black", hjust = 0) +
  draw_text(text = "Plexiform", x = .2, y = .76, size = 18, color = "black", hjust = 0) +
  draw_text(text = "ANF", x = .2, y = .57, size = 18, color = "black", hjust = 0) +
  draw_text(text = "Pre-Tx MPNST", x = .2, y = .39, size = 18, color = "black", hjust = 0) +
  draw_text(text = "On-Tx MPNST", x = .2, y = .19, size = 18, color = "black", hjust = 0)
ggsave("test2.pdf", p, width = 20, height = 8)

p
+
  draw_text(label = "Plexiform", x = 1, y = 0.8, hjust = 1, vjust = 1, size = 14, color = "black") +
  draw_text(label = "ANF", x = 1, y = 0.6, hjust = 1, vjust = 1, size = 14, color = "black") +
  draw_text(label = "Pre-Tx MPNST", x = 1, y = 0.4, hjust = 1, vjust = 1, size = 14, color = "black") +
  draw_text(label = "On-Tx MPNST", x = 1, y = 0.2, hjust = 1, vjust = 1, size = 14, color = "black")

merged

#########1#########2#########3#########4#########5#########6#########7#########8

bins = bins

ratios %>% ggplot

ratios %>%
  left_join(bins, by = c("chr","start","end")) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line() +
  coord_cartesian(ylim = c(-0.06, .06)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) %>%
  dplyr::filter(!is.na(arm)) %>%
  mutate(chr_arm = paste0(chr, "_",arm)) %>%
  mutate(chr_arm = factor(chr_arm, levels = non_acro_factor)) %>%
  dplyr::filter(library %in% train$library) %>%
  left_join(libs, by = "library") %>%
  mutate(bin = paste0(chr,"_",start,"_",end)) %>%
  mutate(dx = ifelse(is.na(washout), dx_simple,
              ifelse(washout == "y", "pretx",
                     "ontx"))) %>%
  mutate(dx = factor(dx, levels = c("healthy", "plexiform", "anf", "pretx", "ontx"))) %>%
  select(chr, library, dx, bin, ratio.centered, chr_arm, start) %>%
  arrange(chr_arm, start) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line() +
  coord_cartesian(ylim = c(-0.06, .06)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + facet_grid(rows = vars(dx), cols = vars(chr_arm), scales = "free_x")

#########1#########2#########3#########4#########5#########6#########7#########8

bins = read_tsv("/tmp/bins.bed", col_names = c("chr","start","end","arm"))

(chrs <- factor(paste0("chr", 1:22), levels = paste0("chr", 1:22)))

(chr_arms <- c(paste0(chrs, "_p ", chrs, "_q")))
chr_arms = unlist(strsplit(chr_arms, " "))
(non_acro_arms <- chr_arms[!chr_arms %in% c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p")])
(non_acro_factor <- factor(non_acro_arms, levels = non_acro_arms))


ratios %>%
  left_join(bins, by = c("chr","start","end")) %>%
  count(chr, arm) %>% print(n = Inf)

=%>%
  left_join(bins, by = c("chr","start","end")) %>%
  dplyr::filter(!is.na(arm)) %>%
  mutate(chr_arm = paste0(chr, "_",arm)) %>%
  mutate(chr_arm = factor(chr_arm, levels = non_acro_factor)) %>%
  dplyr::filter(library %in% train$library) %>%
  left_join(libs, by = "library") %>%
  mutate(bin = paste0(chr,"_",start,"_",end)) %>%
  mutate(dx = ifelse(is.na(washout), dx_simple,
              ifelse(washout == "y", "pretx",
                     "ontx"))) %>%
  mutate(dx = factor(dx, levels = c("healthy", "plexiform", "anf", "pretx", "ontx"))) %>%
  select(chr, library, dx, bin, ratio.centered, chr_arm, start) %>%
  arrange(chr_arm, start) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line() +
  coord_cartesian(ylim = c(-0.06, .06)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + facet_grid(rows = vars(dx), cols = vars(chr_arm), scales = "free_x")

modrat



arms = read_tsv(arms_tsv)
train = read_tsv(train_tsv)
libs = readRDS(libs_rds)

#########1#########2#########3#########4#########5#########6#########7#########8

arms %>%
  mutate(start = armstart, end = armend) %>%
  select(chr, start, end, arm) %>%
  write.table("arms.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#########1#########2#########3#########4#########5#########6#########7#########8


options(scipen = 9999)
modarms=
  arms %>%
  unite("arm_info", arm:armend, sep = "_") %>%
  mutate(arm = substr(arm_info, 1, 1)) %>%
  spread(arm, arm_info) %>%
  mutate(pstart = ifelse(!is.na(p), 0, NA)) %>%
  mutate(pend = gsub("^.*_","",p)) %>%
  mutate(qstart = gsub("^.*q_(\\d+)_(\\d+)$", "\\1", q)) %>%
  mutate(qend = gsub("^.*_","",q)) %>%
  select(chr, pstart, pend, qstart, qend)
modarms
sort(unique(paste0(bins$chr,"_",bins$arm)))
chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
         "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

non_acro_factor= c("chr1_p","chr1_q","chr2_p","chr2_q","chr3_p","chr3_q","chr4_p","chr4_q","chr5_p","chr5_q",
                   "chr6_p","chr6_q","chr7_p","chr7_q","chr8_p","chr8_q","chr9_p","chr9_q","chr10_p","chr10_q",
                   "chr11_p","chr11_q","chr12_p","chr12_q", "chr13_p","chr13_q","chr14_p","chr14_q","chr15_p","chr15_q",
                   "chr16_p","chr16_q","chr17_p","chr17_q","chr18_p","chr18_q","chr19_p","chr19_q","chr20_p","chr20_q",
                   "chr21_p","chr21_q","chr22_p","chr22_q")
ratios %>% mutate(chr = factor(chr, levels = chrs)) %>%
  filter(!is.na(chr)) %>%
  left_join(modarms, by = "chr") %>%
  mutate(arm = ifelse(end < qstart, "p", "q")) %>%
  count(chr, arm) %>% arrange(chr, arm) %>% print(n = Inf) %>%

  mutate(chr_arm = paste0(chr, "_",arm)) %>%
  mutate(chr_arm = factor(chr_arm, levels = non_acro_factor)) %>%
  filter(library %in% train$library) %>%
  left_join(libs, by = "library") %>%
  mutate(bin = paste0(chr,"_",start,"_",end)) %>%
  mutate(dx = ifelse(is.na(washout), dx_simple,
              ifelse(washout == "y", "pretx",
                     "ontx"))) %>%
  mutate(dx = factor(dx, levels = c("healthy", "plexiform", "anf", "pretx", "ontx"))) %>%
  select(chr, library, dx, bin, ratio.centered, chr_arm, start) %>%
  arrange(chr_arm, start) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line() +
  coord_cartesian(ylim = c(-0.03, .03)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + facet_grid(rows = vars(dx), cols = vars(chr), scales = "free_x")

modrat


df_now = ratios %>% filter(chr %in% chrs) %>%
  
  left_join(modarms, by = "chr") %>%
  mutate(arm = ifelse(end < qstart, "p", "q")) %>%
  mutate(chr_arm = paste0(chr, "_",arm)) %>%
  mutate(chr_arm = factor(chr_arm, levels = non_acro_factor)) %>%
  filter(library %in% train$library) %>%
  left_join(libs, by = "library") %>%
  mutate(chr = factor(chr, levels = chrs)) %>% 
  mutate(bin = paste0(chr,"_",start,"_",end)) %>%
  mutate(dx = ifelse(is.na(washout), dx_simple,
                     ifelse(washout == "y", "pretx",
                            "ontx"))) %>%
  mutate(dx = factor(dx, levels = c("healthy", "plexiform", "anf", "pretx", "ontx"))) %>%
  select(chr, library, dx, bin, ratio.centered, chr_arm, start) %>%
  arrange(chr_arm, start)

modrat =
  ratios %>% filter(chr %in% chrs)  %>% 
  left_join(modarms, by = "chr") %>%
  mutate(arm = ifelse(end < qstart, "p", "q")) %>%
  mutate(chr_arm = paste0(chr, "_",arm)) %>%
  mutate(chr_arm = factor(chr_arm, levels = non_acro_factor)) %>%
  filter(library %in% train$library) %>%
  left_join(libs, by = "library") %>%
  mutate(chr = factor(chr, levels = chrs)) %>%
  mutate(bin = paste0(chr,"_",start,"_",end)) %>%
  mutate(dx = ifelse(is.na(washout), dx_simple,
              ifelse(washout == "y", "pretx",
                     "ontx"))) %>%
  mutate(dx = factor(dx, levels = c("healthy", "plexiform", "anf", "pretx", "ontx"))) %>%
  select(chr, library, dx, bin, ratio.centered, chr_arm, start) %>%
  arrange(chr_arm, start) %>%
  ggplot(., aes(x = bin, y = ratio.centered, group = library)) + geom_line() +
  coord_cartesian(ylim = c(-0.03, .03)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + facet_grid(rows = vars(dx), cols = vars(chr), scales = "free_x")

modrat

frags = read.table("/data/CCRCCDI/analysis/ccrtegs65/nf1_fragmentomics_data_11-3-23/inputs/frag_counts.tsv", header = T)


names(frags) = c("library_id", "length", "chr", "start", "end","gc", "count")
head(frags)

test =
  frags %>%
  pivot_wider(names_from = length, values_from = count) %>%
  mutate(fract = short/long) %>%
  select(library_id, chr, start, end, fract) %>%
  pivot_wider(names_from = library_id, values_from = fract) %>%
  mutate_at(vars(starts_with("lib")), ~(scale(.) %>% as.vector))

#https://stackoverflow.com/questions/15215457/standardize-data-columns-in-r

tmp_dx = data.frame(library_id = c("lib060","lib061","lib065","lib067","lib182","lib175","lib129"),
                    current_dx = c(rep("healthy", 4), rep("mpnst", 3)))

test2 = test %>% pivot_longer(
                 cols = starts_with("lib"),
                 names_to = "library_id",
                 values_to = "fract"
                 ) %>%
  left_join(tmp_dx, by = "library_id")

ggplot(test2, aes(x = start, y = fract, group = library_id, color = current_dx)) + geom_smooth(span = .1, se = F) + facet_grid(~chr, scales = "free")

test3 = test %>% pivot_longer(
                 cols = starts_with("lib"),
                 names_to = "library_id",
                 values_to = "fract"
                 ) %>%
  left_join(tmp_dx, by = "library_id") %>%
  filter(chr %in% c("chr17"))


plt1 =
ggplot(test3, aes(x = start, y = fract, group = library_id, color = current_dx)) + geom_smooth(span = .2, se = F) + geom_point() + facet_grid(~chr, scales = "free")
plt1

ggsave("/tmp/test.png", plt1)

head(test)


mat = test %>% select(starts_with("lib"))
mat = as.matrix(mat)
colnames(mat) = tmp_dx$current_dx
rownames(mat) = paste(test$chr,test$start,test$end, sep = "_")

mat = t(mat)
head(mat)

pca = prcomp(mat)

summary(pca)$importance

# > 90% of variance in first 5 components

head(pca$sdev)

head(pca$rotation)

test = pca$rotation[,1:5]

head(test)

test2 = as_tibble(test)
test2
# Get principle component 1 & 2 values
(pve_pc1=round(100*summary(pca)$importance[2,1]))
(pve_pc2=round(100*summary(pca)$importance[2,2]))

row.names(pca$x)

library(ggrepel)

as_tibble(as.data.frame(pca$x))

dx = data.frame(library_id = row.names(pca$x))

dx = dx %>% left_join(libraries, by = "library_id") %>% left_join(specimens, by = "specimen_id") %>% left_join(subjects, by = "subject_id") %>% select(library_id, current_dx) %>% filter(!is.na(current_dx))

mat = pca$x
mat=mat[dx$library_id,]
row.names(mat) = dx$current_dx

pca_plot = as.data.frame(mat) %>%
  as_tibble() %>%
  mutate(current_dx = row.names(mat)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = current_dx)) + geom_point() + xlab("PC1 23% var") + ylab("PC2 7% var")

pca_plot

ggsave("/tmp/test.png", pca_plot, width = 4, height = 4)
+
  theme_cowplot() +
  xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
  #geom_density_2d_filled(aes(group = tx), alpha = 0.3) +
  geom_text_repel()
pca_plot


pca_plot = as_tibble(as.data.frame(pca$x)) %>%
  mutate(sample_id = as.factor(row.names(pca$x))) %>%
  left_join(sampleTable, by = "sample_id") %>%
  mutate(label = sample_id) %>%
  mutate(label = ifelse(sample_id == "s_6", "s_6", NA)) %>%
  ggplot(., aes(x = PC1, y = PC2, label = label, group = tx)) +
  geom_point(aes(color = strain, shape = tx), size = 4) +
  theme_cowplot() +
  xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
  ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
  #geom_density_2d_filled(aes(group = tx), alpha = 0.3) +
  geom_text_repel()
pca_plot

head(mat)
         )
test %>% mutate_each_(list(~scale(.) %>% as.vector), vars = c("lib060","lib061"))




dat2 <- dat %>% mutate_at(c("y", "z"), ~(scale(.) %>% as.vector))


at %>% mutate_each_(list(~scale(.) %>% as.vector), vars=c("y","z"))

head(test)
