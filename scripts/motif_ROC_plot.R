#!/usr/bin/env Rscript

#Comparison individual F-profiles ability to distinguish two cohorts (one versus one or one versus all). 10 repeats of 5-fold CV

#Parse arguments
args = commandArgs(trailingOnly = TRUE)
fprofiles_tsv = args[1]
libraries_tsv = args[2]
plot_out = args[3]
csv_out = args[4]


#fprofiles_tsv = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/features/fprofiles.tsv"
#libraries_tsv = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/input/input_libraries_nf1_training_set.tsv"
#plot_out = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_ROC_healthy_v_MPNST.pdf"
#csv_out = "/data/CCRCCDI/analysis/ccrtegs65/Analysis_NF1_training/nf1_training_set_new/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_AUC_healthy_v_MPNST.csv"


# Load necessary packages
library(readr)
library(dplyr)
library(pROC)
library(ggplot2)
library(knitr)
library(tidyr)
library(cowplot)
library(caret)

#get cohort
fname <- basename(plot_out)
fname
comparison =gsub(".pdf","", gsub("fprofile_ROC_","", fname))
comparison
cohort_1 = strsplit(comparison, "_v_")[[1]][1]
cohort_2 = strsplit(comparison, "_v_")[[1]][2]
cohorts_now = c(cohort_1 , cohort_2)
cohorts_now

# Load data frames, join data frames
fprofiles <- read_tsv(fprofiles_tsv)

# load libraries with only the needed cohorts
libraries <- read_tsv(libraries_tsv) %>%
  filter(cohort %in% cohorts_now)
dim(libraries)
unique(libraries$cohort)
colnames(fprofiles)[1]<-"fprofile"
colnames(fprofiles)[-1] <- sub("_motifs$", "", colnames(fprofiles)[-1]) #Note the fprofiles_tsv appends _motifs to library names. Pre-processing to remove
fprofiles_long <- pivot_longer(fprofiles, cols = -1, names_to = "library", values_to = "pct")
data <- fprofiles_long %>%
  filter(library %in% libraries$library) %>%
  left_join(libraries, by = "library") %>% na.omit(.)
dim(data)
#Define fprofile names and cohorts
fprofiles_list <- unique(fprofiles$fprofile)
unique(data$cohort)

cohorts <- unique(data$cohort)
if (length(cohorts) != 2) {
  stop("This script is designed for the comparison of two cohorts")
}

# Initialize an empty data frame for the results table
roc_results <- data.frame(fprofile = character(), AUC = numeric(), J_cutpoint = numeric(), stringsAsFactors = FALSE)

# Create base plot for combined ROC
roc_plot <- ggplot() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme_half_open()

# Get output directory from plot_out
out_dir <- dirname(plot_out)
dir.create(out_dir, showWarnings = F)
# Loop through each fprofile to calculate and plot ROC

for (fprofile_name in fprofiles_list) {
  
  subset_data <- data %>%
    filter(fprofile == !!fprofile_name) %>%  
    select(pct, cohort)
  
  # Calculate the main ROC curve for plotting
  main_roc <- roc(subset_data$cohort, subset_data$pct, levels = cohorts)
  
  # Set up cross-validation: 10 repeats of 5-fold CV
  train_control <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # Perform CV and collect AUC values
  auc_values <- numeric()
  
  tryCatch({
    for (i in 1:10) {
      folds <- createFolds(subset_data$cohort, k = 5, returnTrain = TRUE)
      for (fold in folds) {
        train <- subset_data[fold, ]
        test <- subset_data[-fold, ]
        
        fold_roc <- roc(test$cohort, test$pct, levels = cohorts)
        auc_values <- c(auc_values, auc(fold_roc))
      }
    }
    
    mean_auc <- mean(auc_values)
    sd_auc <- sd(auc_values)
    
    roc_results <- rbind(roc_results, data.frame(fprofile = fprofile_name, AUC = mean_auc, AUC_sd = sd_auc))
    
    # Add ROC curve to combined plot
    roc_coords <- coords(main_roc, "all")
    roc_plot <- roc_plot + 
      geom_path(data = data.frame(specificity = 1 - roc_coords$specificity,
                                 sensitivity = roc_coords$sensitivity,
                                 fprofile = fprofile_name),
                aes(x = specificity, y = sensitivity, color = fprofile),
                size = 1.2)
    
    # Create individual ROC plot for this F-profile
    individual_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      geom_path(data = data.frame(specificity = 1 - roc_coords$specificity,
                                 sensitivity = roc_coords$sensitivity),
                aes(x = specificity, y = sensitivity),
                size = 1.2, color = rainbow(length(fprofiles_list))[which(fprofiles_list == fprofile_name)]) +
      labs(x = "False Positive Rate", 
           y = "True Positive Rate",
           title = paste0(fprofile_name, " ROC Curve\nAUC: ", round(mean_auc, 3), 
                         " ± ", round(sd_auc, 3))) +
      theme_half_open()
    
    # Save individual plot
    individual_plot_file <- file.path(out_dir, paste0(fprofile_name, "_ROC_",comparison,".pdf"))
    ggsave(individual_plot_file, plot = individual_plot, width = 8, height = 6)
    
  }, error = function(e) {
    message(paste("Error in ROC computation for fprofile", fprofile_name, ":", e$message))
  })
}

# Add color scale and legend to combined plot
roc_plot <- roc_plot + 
  scale_color_manual(name = "F-profile", 
                    values = rainbow(length(fprofiles_list)), 
                    labels = sapply(1:length(fprofiles_list), function(i) {
                      fprofile <- fprofiles_list[i]
                      auc_val <- round(roc_results$AUC[roc_results$fprofile == fprofile], 2)
                      paste(fprofile, "AUC:", auc_val)
                    })) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

# Print results table
kable(roc_results, caption = "AUC and J-cutpoint for Each Fprofile")
roc_plot
# Save outputs
write.csv(roc_results, csv_out, row.names = FALSE)
ggsave(plot_out, plot = roc_plot, width = 10, height = 8)
