# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)
library(pROC)
library(optparse)
library(tidyverse)
library(glue)
library(scales)

option_list <- list(
  make_option("--config", type="character", help="Path to samples config file"),
  make_option("--input", type="character", help="Path to samples config file"),
  make_option("--data_dir", type="character", help="Data directory path"),
  make_option("--cohort1", type="character", help="First cohort"),
  make_option("--cohort2", type="character", help="Second cohort"),
  make_option("--normal_control", type="character", help="Normal control cohort"),
  make_option("--output_dir", type="character", help="Output file path"),
  make_option("--status", type="character", help="Status file path"),
  make_option("--log", type="character", help="Log file path"),
  make_option("--ichor_tf", type="character", help="Tumor fraction file"),
  make_option("--nmf_score", type="character", help="NMF AUCs file"),
  make_option("--binwise_score", type="character", help="Bin-wise probabilities file")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt = list(
#  config = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/input/input_libraries_nf1_training_set.tsv",
#  input = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/combined_dataset_healthy_v_MPNST.csv",
#  data_dir = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new",
#  cohort1 = "healthy",
#  cohort2 = "MPNST",
#  normal_control = "healthy",
#  output_dir = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/plots/voting",
#  log = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/logs/combine_datasets/plot_voting_healthy_v_MPNST.log",
#  nmf_bestauc = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/nmf/aucs/best_med_aucs.csv",
#  binwise_score = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/bin_wise_fragmentomics/scores/healthy_v_MPNST.csv",
#  ichor_tf = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/ichorCNA/tf.csv"
#)

dir.create(opt$output_dir, showWarnings = F)

# load sample info
log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ...)
}
log_message("Reading samples config file")
samples_df <- read_tsv(opt$config, show_col_types = FALSE) %>%
  mutate(library = str_to_upper(library))

# load binwise_score
binwise_df <- read_csv(opt$binwise_score, show_col_types = FALSE) %>%
  mutate(library = str_to_upper(library))

# load nmf scores
best_n_components =  sub(".*\\.n=([0-9]+)\\.csv$", "\\1", opt$nmf_score)
best_n_components <- as.integer(best_n_components)
best_n_components 

nmf_df = read_csv(opt$nmf_score, show_col_types = FALSE) %>%
  mutate(library = str_to_upper(library))
head(nmf_df)

# load ichore tf
tf_df = read_csv(opt$ichor_tf, show_col_types = FALSE)
head(tf_df)
tf_df <- tf_df %>%
  mutate(
    ichor_pred_class = case_when(
      Tumor_Fraction > 0.03 ~ "1",
      TRUE ~ "0"
    )
  )


#combine score
# Read combined dataset
log_message("Reading combined dataset")
combined_df <- read_csv(opt$input, show_col_types = FALSE)

# Filter for relevant cohorts
combined_df <- combined_df 
head(combined_df)

combined_df = combined_df %>%
  left_join(binwise_df[,c("library","med_score_cutoff","med_score_pred_class")],by="library") %>%
  rename(med_score_pred_class_binwise = med_score_pred_class) %>%
  rename(med_score_cutoff_binwise = med_score_cutoff)

combined_df = combined_df %>%
  left_join(nmf_df[,c("library","med_score_cutoff","med_score_pred_class")],by="library") %>%
  rename(med_score_pred_class_nmf = med_score_pred_class) %>%
  rename(med_score_cutoff_nmf = med_score_cutoff)
combined_df$best_n_components_NMF = best_n_components

combined_df = combined_df %>%
  mutate(
    ichor_pred_class = case_when(
      Tumor_Fraction > 0.03 ~ "1",
      TRUE ~ "0"
    )
  ) %>%
  mutate(ichor_cutoff = 0.03) %>%
  mutate(
    num_positive_calls = rowSums(
      across(
        c(
          ichor_pred_class,
          med_score_pred_class_nmf,
          med_score_pred_class_binwise
        ),
        ~ as.numeric(.x)
      ),
      na.rm = TRUE
    )
  )


if (all(grepl("^LIB", combined_df$library))) {
  combined_df <- combined_df %>%
    mutate(
      cohort = factor(cohort, levels = unique(cohort)),
      lib_num = as.numeric(gsub("LIB", "", library))
    ) %>%
    arrange(cohort, lib_num) %>%
    select(-lib_num)
} else {
  combined_df <- combined_df %>%
    mutate(
      cohort = factor(cohort, levels = unique(combined_df$cohort))
    ) %>%
    arrange(cohort, library)
}

colnames(combined_df)
heatmap_data = as.data.frame(combined_df[,c("library","ichor_pred_class", "med_score_pred_class_nmf"  , "med_score_pred_class_binwise")])
row.names(heatmap_data) = combined_df$library
head(heatmap_data)
colnames(heatmap_data) = c("library","Size Selected ichorCNA", "Fragment size NMF"  , "Binwise Fragmentomics")

score_data = as.data.frame(combined_df[,c("library","Tumor_Fraction", "nmf_prob", "binwise_prob")])
row.names(score_data) = combined_df$library
colnames(score_data) = c("library","Size Selected ichorCNA", "Fragment size NMF"  , "Binwise Fragmentomics")

# Get cohort in the same order as heatmap_data rows
cohort_anno <- combined_df$cohort
names(cohort_anno) <- combined_df$library

# Define colors for cohort annotation
healthy_col <- "gray85"
pn_col      <- "#ffce00"
an_col      <- "#D55E00"
mpnst_col   <- "#5A4FCF"
cohort_colors <- c("MPNST" = mpnst_col, "healthy" = healthy_col, "PN" = pn_col, "AN" = an_col)
color_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#BAAC65", "#86B875",
                   "#F0E442", "#0072B2", "#CC79A7", "#4CB9CC", "#96AAE1",
                   "#F8766D", "#39B600", "#9590FF", "#E76BF3")

cohorts <- unique(combined_df$cohort)
n_cohorts <- length(cohorts)

if (n_cohorts > length(color_palette)) {
  cohort_ori_colors <- hue_pal()(n_cohorts)
} else {
  cohort_ori_colors <- color_palette[1:n_cohorts]
}

names(cohort_ori_colors) <- cohorts

# overwrite colors for known cohort names
matched_cohorts <- intersect(names(cohort_ori_colors), names(cohort_colors))
cohort_ori_colors[matched_cohorts] <- cohort_colors[matched_cohorts]


# Create the cohort annotation
cohort_ha <- rowAnnotation(
  Cohort = cohort_anno,
  col = list(Cohort = cohort_ori_colors),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10),
  width = unit(0.1, "cm")
)


# Create color mapping
col_fun <- colorRamp2(c(0, 1), c(cohort_colors[opt$cohort1], cohort_colors[opt$cohort2]))

plot(cohort_ha)


## prepare AUCs for each metric in the heatmap

#cohort_vec <- combined_df$cohort
#names(cohort_vec) <- combined_df$library
#cohort_vec <- cohort_vec[row.names(heatmap_data)]
#roc_call_columns <- c(
#  "ichor_pred_class",
#  "med_score_pred_class_nmf"  ,
#  "med_score_pred_class_binwise"
#
# 
#  
#)

#auc_values <- sapply(roc_call_columns, function(col) {
#  metric_vals <- combined_df[[col]]
#  metric_numeric <- ifelse(metric_vals == 1, 1, 0)
#  valid_idx <- !is.na(metric_numeric) & !is.na(cohort_vec)
#  if (length(unique(cohort_vec[valid_idx])) == 2) {
#    roc_obj <- roc(cohort_vec[valid_idx], metric_numeric[valid_idx], levels = c(opt$cohort1, opt$cohort2), direction = "<")
#    as.numeric(auc(roc_obj))
#  } else {
#    NA
#  }
#})
#auc_values
#names(auc_values) <- c("Size Selected ichorCNA", "Fragment size NMF"  , "Binwise Fragmentomics")
#auc_values


# Define colors for the barplot based on the midpoint

positive_call = combined_df$num_positive_calls
names(positive_call) = combined_df$library
midpoint <- (length(colnames(heatmap_data))-1) / 2
print(midpoint)

fbar_colors <- ifelse(positive_call> midpoint, cohort_colors[opt$cohort2], cohort_colors[opt$cohort1])

# Create the barplot annotation
ha <- rowAnnotation(
  "Number of positive Predictions" = anno_barplot(
    positive_call,
    width = unit(2, "cm"),
    gp = gpar(fill = fbar_colors, col = "black"),
    axis = TRUE,
    axis_param = list(
      at = c(0, midpoint, max(positive_call)),
      labels = c(0, midpoint, max(positive_call)),
      gp = gpar(fontsize = 8)
    ),
    border = FALSE,
    legend_param = list(
      title = paste0("Number of\npositive Predictions\n",opt$cohort1,"_v_",opt$cohort2)
    ),
    panel_fun = function(index, value, ...) {
      grid.lines(
        x = c(0, 1),
        y = rep(midpoint / max(positive_call), 2),
        gp = gpar(lty = 2, col = "black", lwd = 1)
      )
    }
  ),
  annotation_name_side = "bottom",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)
#plot(ha)

# Define AUC labels
#auc_labels <- sprintf("%.2f", auc_values)

format_score_label <- function(value) {
  if (is.na(value)) {
    return("")
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (is.na(numeric_value)) {
    return(as.character(value))
  }
  if (abs(numeric_value) > 0 && abs(numeric_value) < 0.001) {
    return(sprintf("%.1e", numeric_value))
  }
  label <- sub("\\.?0+$", "", sprintf("%.3f", numeric_value))
  if (label == "-0") {
    label <- "0"
  }
  label
}

make_score_cell_fun <- function(score_matrix, call_matrix, fontsize = 5) {
  score_matrix <- as.matrix(score_matrix)
  call_matrix <- as.matrix(call_matrix)
  function(j, i, x, y, width, height, fill) {
    label <- format_score_label(score_matrix[i, j])
    if (label == "") {
      return(NULL)
    }
    call_value <- suppressWarnings(as.numeric(call_matrix[i, j]))
    text_color <- if (!is.na(call_value) && call_value >= 0.5) "white" else "black"
    grid.text(
      label,
      x = x,
      y = y,
      gp = gpar(fontsize = fontsize, col = text_color)
    )
  }
}

# Create the main heatmap
heatmap_matrix = as.matrix(heatmap_data[-1])
heatmap_matrix_label <- ifelse(
  heatmap_matrix == 1,
  opt$cohort2,
  opt$cohort1
)

heatmap_matrix <- as.matrix(heatmap_data[-1])
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) <- heatmap_data[[1]]

score_matrix <- as.matrix(score_data[-1])
score_matrix <- apply(score_matrix, 2, function(x) suppressWarnings(as.numeric(x)))
rownames(score_matrix) <- score_data[[1]]
colnames(score_matrix) <- colnames(heatmap_matrix)

dimnames(heatmap_matrix_label) <- dimnames(heatmap_matrix)
col_fun <- colorRamp2(c(0, 1), c(cohort_colors[opt$cohort1], cohort_colors[opt$cohort2]))
n_samples <- nrow(combined_df)
score_label_fontsize <- ifelse(n_samples < 20, 7, 5)

row_name_width <- max_text_width(
  row.names(heatmap_matrix),
  gp = gpar(fontsize = 8)
) + unit(3, "mm")

ht_main <- Heatmap(
  matrix = heatmap_matrix,
  name = paste0("Prediction: \n", opt$cohort1, "_v_",opt$cohort2),
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 1),
  width = unit(ncol(heatmap_matrix) * 24, "mm"),
  column_gap = unit(4, "mm"),
  column_names_rot = 45,
  cell_fun = make_score_cell_fun(
    score_matrix,
    heatmap_matrix,
    fontsize = score_label_fontsize
  )
)

# Combine the heatmaps: main heatmap, then cohort and barplot row annotations
ht_list <- ht_main + cohort_ha + ha
ht_list

#set auto height
height_per_sample <- 8 / 56
plot_height <- n_samples * height_per_sample
plot_height
if (plot_height < 4) {
  plot_height = 4
}

#plot_width <- 7
#plot_padding <- unit(c(2, 2, 2, 2), "mm")
plot_width <- 10
plot_padding <- unit(c(5, 8, 5, 12), "mm")
# Draw the heatmap
pdf(paste0( opt$output_dir,"/roc_calls_heatmap_",opt$cohort1, "_v_", opt$cohort2,".pdf"), width = plot_width, height = plot_height)
draw(ht_list, annotation_legend_side = "right", padding = plot_padding)
dev.off()

#png(paste0( opt$output_dir,"/roc_calls_heatmap_",opt$cohort1, "_v_", opt$cohort2,".png"), width = 1200, height = plot_height * 100, res = 100)
#draw(ht_list, annotation_legend_side = "right")
#dev.off()

library(ragg)

agg_png(
  filename = paste0( opt$output_dir,"/roc_calls_heatmap_",opt$cohort1, "_v_", opt$cohort2,".png"),
  width = plot_width * 300,
  height = plot_height * 300,
  units = "px",
  res = 300
)

draw(ht_list, annotation_legend_side = "right", padding = plot_padding)
dev.off()

# Save to CSV
#combined_df = combined_df[,c("library","cohort","Tumor_Fraction","lt_ctdna_burden","mds","fprof5","fprof6","fprof2","nmf_prob","binwise_prob","med_score_cutoff_binwise","med_score_pred_class_binwise","med_score_cutoff_nmf","med_score_pred_class_nmf","best_n_components_NMF","ichor_pred_class","ichor_cutoff","num_positive_calls")]

write.csv(combined_df, file = paste0( opt$output_dir,"/roc_calls_heatmap_",opt$cohort1, "_v_", opt$cohort2,"_majority_vote.csv"), row.names = FALSE)

#save.image("img/heatmap.RData")
