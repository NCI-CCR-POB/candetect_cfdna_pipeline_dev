# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)
library(pROC)


setwd("/Volumes/ncipobrainbow/ShernLab/Taylor/Research/projects/in_progress/SARC031/figures")
setwd("/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training")
# Define color variables
healthy_col <- "#ffce00"  
mpnst_col <- "#5A4FCF"    
load("/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/test/heatmap.RData")

##### Load data #####
combined_sample_data <- read.csv("results/combined_sample_data.csv")

# Load subj_key_sample_recist_pet.csv for revised_sample_id key
subj_key <- read.csv("data/subj_key_sample_recist_pet.csv")
subj_key_map <- setNames(as.character(subj_key$revised_sample_id), as.character(subj_key$sample_id))

#### Mapping ####
# Map sample_id to revised_sample_id for cleaned row labels
revised_labels <- as.character(subj_key_map[as.character(combined_sample_data$sample_id)])
# If NA, fallback to original sample_id
revised_labels[is.na(revised_labels)] <- as.character(combined_sample_data$sample_id[is.na(revised_labels)])
# Now, change e.g. 015-006_CY0 to 015-006 Cycle 0
cleaned_samplenames <- sapply(revised_labels, function(x) {
  if (grepl("^(\\d{3}-\\d{3,4})_CY(\\d+)$", x)) {
    m <- regmatches(x, regexec("^(\\d{3}-\\d{3,4})_CY(\\d+)$", x))[[1]]
    paste(m[2], "Cycle", m[3])
  } else if (grepl("^HEALTHY", x, ignore.case = TRUE)) {
    num <- stringr::str_extract(x, "\\d+")
    paste("Healthy Control", num)
  } else {
    x
  }
})
row.names(combined_sample_data) <- cleaned_samplenames

# Define the columns for the heatmap
roc_call_columns <- c(
  "ichor_tf_3_call",
  "ichor_roc_call",
  "fragle_roc_call",
  "mds_roc_call",
  "nmf_roc_call",
  "binwise_roc_call",
  "fprofiles_fprof1_roc_call",
  "fprofiles_fprof3_roc_call",
  "fprofiles_fprof5_roc_call"
)

# Create a mapping for the calls
call_mapping <- c("mpnst" = 1, "healthy" = 0)

#### Create the heatmap data for all metrics #####
heatmap_data <- combined_sample_data
heatmap_data <- heatmap_data[roc_call_columns]
heatmap_data[] <- lapply(heatmap_data, function(x) call_mapping[x])

cleaned_columnnames <- c(
  "ichorCNA",
  "Size Selected ichorCNA",
  "Fragle ctDNA Estimate",
  "Motif Diversity Score",
  "Fragment size NMF",
  "Binwise Fragmentomics",
  "F-profile 1",
  "F-profile 3",
  "F-profile 5"
)
names(heatmap_data)<-cleaned_columnnames

# Get mpnst_calls in the same order as heatmap_data rows
mpnst_calls <- combined_sample_data$mpnst_calls
names(mpnst_calls) <- cleaned_samplenames
mpnst_calls <- mpnst_calls[row.names(heatmap_data)]

# Get cohort in the same order as heatmap_data rows
cohort_anno <- combined_sample_data$cohort
names(cohort_anno) <- cleaned_samplenames
cohort_anno <- cohort_anno[row.names(heatmap_data)]

# Define colors for cohort annotation
cohort_colors <- c("mpnst" = mpnst_col, "healthy" = healthy_col)

# Create the cohort annotation
cohort_ha <- rowAnnotation(
  Cohort = cohort_anno,
  col = list(Cohort = cohort_colors),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10),
  width = unit(0.1, "cm")
)

# Create color mapping
col_fun <- colorRamp2(c(0, 1), c(healthy_col, mpnst_col))

# Calculate AUCs for each metric in the heatmap
cohort_vec <- combined_sample_data$cohort
names(cohort_vec) <- cleaned_samplenames
cohort_vec <- cohort_vec[row.names(heatmap_data)]

auc_values <- sapply(roc_call_columns, function(col) {
  metric_vals <- combined_sample_data[[col]]
  metric_numeric <- ifelse(metric_vals == "mpnst", 1, 0)
  valid_idx <- !is.na(metric_numeric) & !is.na(cohort_vec)
  if (length(unique(cohort_vec[valid_idx])) == 2) {
    roc_obj <- roc(cohort_vec[valid_idx], metric_numeric[valid_idx], levels = c("healthy", "mpnst"), direction = "<")
    as.numeric(auc(roc_obj))
  } else {
    NA
  }
})
names(auc_values) <- cleaned_columnnames

# Define colors for the barplot based on the midpoint
midpoint <- length(roc_call_columns) / 2
bar_colors <- ifelse(mpnst_calls > midpoint, mpnst_col, healthy_col)

# Identify healthy and non-healthy sample names
is_healthy <- grepl("^Healthy Control \\d+$", cleaned_samplenames)

# Extract the numeric part for sorting
healthy_nums <- as.numeric(sub("Healthy Control ", "", cleaned_samplenames[is_healthy]))

# Get the order for healthy controls
healthy_order <- order(healthy_nums, na.last = NA)

# Indices for healthy and non-healthy
healthy_indices <- which(is_healthy)[healthy_order]
non_healthy_indices <- which(!is_healthy)

# Combine: non-healthy first, then healthy controls in numeric order (or reverse if you prefer)
final_order <- c(non_healthy_indices, healthy_indices)

# Reorder everything
heatmap_data <- heatmap_data[final_order, , drop = FALSE]
mpnst_calls <- mpnst_calls[final_order]
cohort_anno <- cohort_anno[final_order]
cleaned_samplenames <- cleaned_samplenames[final_order]
row.names(heatmap_data) <- cleaned_samplenames

# Set the names of cohort_anno to match the new row names of heatmap_data
names(cohort_anno) <- row.names(heatmap_data)
names(mpnst_calls) <- row.names(heatmap_data)

# Create the barplot annotation
ha <- rowAnnotation(
  "Number of MPNST Predictions" = anno_barplot(
    mpnst_calls,
    width = unit(2, "cm"),
    gp = gpar(fill = bar_colors, col = "black"),
    axis = TRUE,
    axis_param = list(
      at = c(0, midpoint, max(mpnst_calls)),
      labels = c(0, midpoint, max(mpnst_calls)),
      gp = gpar(fontsize = 8)
    ),
    border = FALSE,
    legend_param = list(
      title = "Number of\nMPNST Predictions"
    ),
    panel_fun = function(index, value, ...) {
      grid.lines(
        x = c(0, 1),
        y = rep(midpoint / max(mpnst_calls), 2),
        gp = gpar(lty = 2, col = "black", lwd = 1)
      )
    }
  ),
  annotation_name_side = "bottom",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

# Define AUC labels
auc_labels <- sprintf("%.2f", auc_values)

# Create the main heatmap
ht_main <- Heatmap(
  matrix = as.matrix(heatmap_data),
  name = "Prediction",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 1),
  width = unit(ncol(heatmap_data) * 5, "mm"),
  column_gap = unit(4, "mm"),
  column_names_rot = 45,
  top_annotation = HeatmapAnnotation(
    AUC = anno_barplot(
      auc_values,
      ylim = c(0, 1),
      border = FALSE,
      gp = gpar(fill = "#123456"),
      axis_param = list(
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        gp = gpar(fontsize = 8)
      ),
      panel_fun = function(index, value, ...) {
        n <- length(value)
        for (i in seq_len(n)) {
          grid.text(
            auc_labels[i],
            x = (i - 0.5) / n,
            y = min(1, value[i]) + 0.08,
            gp = gpar(fontsize = 9)
          )
        }
      }
    ),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 10)
  )
)

# Combine the heatmaps: main heatmap, then cohort and barplot row annotations
ht_list <- ht_main + cohort_ha + ha
ht_list

# Draw the heatmap
pdf("img/roc_calls_heatmap.pdf", width = 12, height = 8)
draw(ht_list, annotation_legend_side = "right")
dev.off()

png("img/roc_calls_heatmap.png", width = 1200, height = 800, res = 100)
draw(ht_list, annotation_legend_side = "right")
dev.off()

####REDUCED METRICS HEATMAP####

# Find the top AUC fprofile_roc_call
fprofile_cols <- c("fprofiles_fprof1_roc_call", "fprofiles_fprof3_roc_call", "fprofiles_fprof5_roc_call")
fprofile_aucs <- sapply(fprofile_cols, function(col) {
  metric_vals <- combined_sample_data[[col]]
  metric_numeric <- ifelse(metric_vals == "mpnst", 1, 0)
  valid_idx <- !is.na(metric_numeric) & !is.na(cohort_vec)
  if (length(unique(cohort_vec[valid_idx])) == 2) {
    roc_obj <- roc(cohort_vec[valid_idx], metric_numeric[valid_idx], levels = c("healthy", "mpnst"), direction = "<")
    as.numeric(auc(roc_obj))
  } else {
    NA
  }
})
top_fprofile_col <- fprofile_cols[which.max(fprofile_aucs)]

# Reduce to only columns used in Sundby et al, 2024
reduced_cols <- c("ichor_roc_call", "nmf_roc_call", "binwise_roc_call", top_fprofile_col)
reduced_colnames <- c(
  "Size Selected ichorCNA",
  "Fragment size NMF",
  "Binwise Fragmentomics",
  "F-profile")

reduced_heatmap_data <- combined_sample_data[reduced_cols]
reduced_heatmap_data[] <- lapply(reduced_heatmap_data, function(x) call_mapping[x])
names(reduced_heatmap_data) <- reduced_colnames

# Set row names for reduced_heatmap_data to cleaned_samplenames
row.names(reduced_heatmap_data) <- cleaned_samplenames

# Identify healthy and non-healthy sample names for reduced heatmap
is_healthy_reduced <- grepl("^Healthy Control \\d+$", cleaned_samplenames)
healthy_nums_reduced <- as.numeric(sub("Healthy Control ", "", cleaned_samplenames[is_healthy_reduced]))
healthy_order_reduced <- order(healthy_nums_reduced, na.last = NA)
healthy_indices_reduced <- which(is_healthy_reduced)[healthy_order_reduced]
non_healthy_indices_reduced <- which(!is_healthy_reduced)
final_order_reduced <- c(non_healthy_indices_reduced, healthy_indices_reduced)

# Subset everything using final_order_reduced
reduced_heatmap_data <- reduced_heatmap_data[final_order_reduced, , drop = FALSE]
reduced_mpnst_calls <- rowSums(reduced_heatmap_data == 1)
names(reduced_mpnst_calls) <- row.names(reduced_heatmap_data)
# If you use cohort_anno in the reduced plot, also subset:
# reduced_cohort_anno <- cohort_anno[final_order_reduced]

#Barplot annotation for reduced heatmap
reduced_midpoint <- length(reduced_cols) / 2
reduced_bar_colors <- ifelse(reduced_mpnst_calls >= reduced_midpoint, mpnst_col, healthy_col)
reduced_ha <- rowAnnotation(
  "Number of MPNST Predictions" = anno_barplot(
    reduced_mpnst_calls,
    width = unit(2, "cm"),
    gp = gpar(fill = reduced_bar_colors, col = "black"),
    axis = TRUE,
    axis_param = list(
      at = c(0, reduced_midpoint, max(reduced_mpnst_calls)),
      labels = c(0, reduced_midpoint, max(reduced_mpnst_calls)),
      gp = gpar(fontsize = 8)
    ),
    border = FALSE,
    legend_param = list(
      title = "Number of\nMPNST Predictions"
    ),
    panel_fun = function(index, value, ...) {
      grid.lines(
        x = c(0, 1),
        y = rep(reduced_midpoint / max(reduced_mpnst_calls), 2),
        gp = gpar(lty = 2, col = "black", lwd = 1)
      )
    }
  ),
  annotation_name_side = "bottom",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

#Cohort annotation (reuse cohort_ha, which is already in correct order)

# Calculate AUCs for reduced metrics
reduced_auc_values <- sapply(reduced_cols, function(col) {
  metric_vals <- combined_sample_data[[col]]
  metric_numeric <- ifelse(metric_vals == "mpnst", 1, 0)
  valid_idx <- !is.na(metric_numeric) & !is.na(cohort_vec)
  if (length(unique(cohort_vec[valid_idx])) == 2) {
    roc_obj <- roc(cohort_vec[valid_idx], metric_numeric[valid_idx], levels = c("healthy", "mpnst"), direction = "<")
    as.numeric(auc(roc_obj))
  } else {
    NA
  }
})
names(reduced_auc_values) <- reduced_colnames
reduced_auc_labels <- sprintf("%.2f", reduced_auc_values)

# Add AUC barplot annotation to reduced heatmap
top_reduced_auc_ha <- HeatmapAnnotation(
  AUC = anno_barplot(
    reduced_auc_values,
    ylim = c(0, 1),
    border = FALSE,
    gp = gpar(fill = "#123456"),
    axis_param = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      gp = gpar(fontsize = 8)
    ),
    panel_fun = function(index, value, ...) {
      n <- length(value)
      for (i in seq_len(n)) {
        grid.text(
          reduced_auc_labels[i],
          x = (i - 0.5) / n,
          y = min(1, value[i]) + 0.08,
          gp = gpar(fontsize = 9)
        )
      }
    }
  ),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10)
)

# Calculate AUC for combined reduced metrics
combined_pred <- ifelse(reduced_mpnst_calls >= reduced_midpoint, "mpnst", "healthy")
combined_pred_factor <- factor(combined_pred, levels = c("healthy", "mpnst"))
cohort_true <- cohort_anno # already in correct order
valid_idx_combined <- !is.na(combined_pred_factor) & !is.na(cohort_true)
if (length(unique(cohort_true[valid_idx_combined])) == 2) {
  combined_roc <- roc(cohort_true[valid_idx_combined], as.numeric(combined_pred_factor[valid_idx_combined] == "mpnst"), levels = c("healthy", "mpnst"), direction = "<")
  combined_auc <- as.numeric(auc(combined_roc))
  cat(sprintf("AUC for combined reduced metrics: %.3f\n", combined_auc))
} else {
  cat("AUC for combined reduced metrics: NA (not enough class diversity)\n")
}

#Create the reduced heatmap with top annotation
reduced_ht_main <- Heatmap(
  matrix = as.matrix(reduced_heatmap_data),
  name = "Prediction",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 1),
  width = unit(ncol(reduced_heatmap_data) * 5, "mm"),
  column_gap = unit(4, "mm"),
  column_names_rot = 45,
  top_annotation = top_reduced_auc_ha
)

#Combine and plot
reduced_ht_list <- reduced_ht_main + cohort_ha + reduced_ha
reduced_ht_list 

pdf("img/roc_calls_heatmap_reduced.pdf", width = 8, height = 8)
draw(reduced_ht_list, annotation_legend_side = "right")
dev.off()

png("img/roc_calls_heatmap_reduced.png", width = 800, height = 800, res = 100)
draw(reduced_ht_list, annotation_legend_side = "right")
dev.off()

# Add majority_decision column to reduced data
majority_decision <- ifelse(reduced_mpnst_calls >= 2, "mpnst", "healthy")
reduced_output_df <- as.data.frame(reduced_heatmap_data)
reduced_output_df$majority_decision <- majority_decision
reduced_output_df$sample_name <- row.names(reduced_output_df)
# Reorder columns to have sample_name first
reduced_output_df <- reduced_output_df[, c("sample_name", setdiff(names(reduced_output_df), "sample_name"))]

# Extract subject from sample_name (e.g., EIT Cycle 0 and EIT Cycle 3 are both subject ELA)
reduced_output_df$subject <- sapply(reduced_output_df$sample_name, function(x) {
  # If matches pattern like 'EIT Cycle 0', extract 'EIT' and map to 'ELA'
  if (grepl("^EIT Cycle", x)) {
    return("ELA")
  } else if (grepl("^Healthy Control", x)) {
    return(x) # treat each healthy as its own subject
  } else {
    # Otherwise, take the first word (subject ID)
    strsplit(x, " ")[[1]][1]
  }
})

# Add subject_majority_vote: if any majority_decision for subject is mpnst, then mpnst, else healthy
subject_majority <- reduced_output_df %>%
  group_by(subject) %>%
  mutate(subject_majority_vote = ifelse(any(majority_decision == "mpnst"), "mpnst", "healthy")) %>%
  ungroup()

# Save to CSV
write.csv(subject_majority, file = "results/reduced_metrics_subject_majority_vote.csv", row.names = FALSE)

save.image("img/heatmap.RData")

