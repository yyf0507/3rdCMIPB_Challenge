library(dplyr)
library(tidyr)
library(tibble)

base_dir = "C:/CMI-PB/data/"
dir_RDS_objects <- paste0(base_dir, "output/")

# `codebase.R` installs required packages and all house keeping functions
# source(paste0(base_dir, "codebase.R"))

data_type <- 1 # 1: raw data
# data_type <- 2 # 2: normalized data
# data_type <- 3 # 3: batchCorrected data

# ==============================================================================

for (data_type in c(1, 2, 3)){
  master_database_data <- readRDS(paste0(dir_RDS_objects, "master_harmonized_data.RDS"))
  
  train_subject_specimen <- master_database_data$training$subject_specimen
  train_subject_specimen$age <- round(as.numeric(difftime(as.Date(train_subject_specimen$date_of_boost), as.Date(train_subject_specimen$year_of_birth), units = 'days') / 365), 2)
  train_subject_specimen <- train_subject_specimen %>%
    mutate(across(everything(), as.character))
  train_subject_specimen_long <- pivot_longer(train_subject_specimen, cols = -subject_id, names_to = 'col', values_to = 'value')
  train_subject_specimen_long <- train_subject_specimen[!duplicated(train_subject_specimen[['subject_id']]), c('subject_id', 'infancy_vac', 'biological_sex', 'ethnicity', 'race', 'age', 'dataset')]
  train_subject_specimen_long <- pivot_longer(train_subject_specimen_long, cols = -subject_id, names_to = 'col', values_to = 'value')
  train_subject_specimen_long$col <- paste(train_subject_specimen_long$col, '(subject_info)')
  
  challenge_subject_specimen <- master_database_data$challenge$subject_specimen
  challenge_subject_specimen$age <- round(as.numeric(difftime(as.Date(challenge_subject_specimen$date_of_boost), as.Date(challenge_subject_specimen$year_of_birth), units = 'days') / 365), 2)
  challenge_subject_specimen <- challenge_subject_specimen %>%
    mutate(across(everything(), as.character))
  challenge_subject_specimen_long <- challenge_subject_specimen[!duplicated(challenge_subject_specimen[['subject_id']]), c('subject_id', 'infancy_vac', 'biological_sex', 'ethnicity', 'race', 'age')]
  challenge_subject_specimen_long <- pivot_longer(challenge_subject_specimen_long, cols = -subject_id, names_to = 'col', values_to = 'value')
  challenge_subject_specimen_long$col <- paste(challenge_subject_specimen_long$col, '(subject_info)')
  
  master_database_data <- readRDS(paste0(dir_RDS_objects, "master_allData_batchCorrected.RDS"))
  
  all_plasma_ab_titer_long <- data.frame(master_database_data$plasma_ab_titer[[data_type]], check.names = FALSE)
  all_plasma_ab_titer_long$col <- paste(row.names(all_plasma_ab_titer_long), '(plasma_ab_titer)')
  all_plasma_ab_titer_long <- pivot_longer(all_plasma_ab_titer_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_plasma_cytokine_concentrations_by_olink_long <- data.frame(master_database_data$plasma_cytokine_concentrations_by_olink[[data_type]], check.names = FALSE)
  all_plasma_cytokine_concentrations_by_olink_long$col <- paste(row.names(all_plasma_cytokine_concentrations_by_olink_long), '(plasma_cytokine_concentrations_by_olink)')
  all_plasma_cytokine_concentrations_by_olink_long <- pivot_longer(all_plasma_cytokine_concentrations_by_olink_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_plasma_cytokine_concentrations_by_legendplex_long <- data.frame(master_database_data$plasma_cytokine_concentrations_by_legendplex[[ifelse(data_type > 1, 2, 1)]], check.names = FALSE) # only raw and normalized
  all_plasma_cytokine_concentrations_by_legendplex_long$col <- paste(row.names(all_plasma_cytokine_concentrations_by_legendplex_long), '(plasma_cytokine_concentrations_by_legendplex)')
  all_plasma_cytokine_concentrations_by_legendplex_long <- pivot_longer(all_plasma_cytokine_concentrations_by_legendplex_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_pbmc_cell_frequency_long <- data.frame(master_database_data$pbmc_cell_frequency[[data_type]], check.names = FALSE)
  all_pbmc_cell_frequency_long$col <- paste(row.names(all_pbmc_cell_frequency_long), '(pbmc_cell_frequency)')
  all_pbmc_cell_frequency_long <- pivot_longer(all_pbmc_cell_frequency_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_pbmc_gene_expression_long <- data.frame(master_database_data$pbmc_gene_expression$tpm[[ifelse(data_type > 1, 2, 1)]], check.names = FALSE) # raw count and tpm, only raw and batch corrected
  all_pbmc_gene_expression_long$col <- paste(row.names(all_pbmc_gene_expression_long), '(pbmc_gene_expression)')
  all_pbmc_gene_expression_long <- pivot_longer(all_pbmc_gene_expression_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_t_cell_polarization_long <- data.frame(master_database_data$t_cell_polarization[[1]], check.names = FALSE) # only raw
  all_t_cell_polarization_long$col <- paste(row.names(all_t_cell_polarization_long), '(t_cell_polarization)')
  all_t_cell_polarization_long <- pivot_longer(all_t_cell_polarization_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_t_cell_activation_long <- data.frame(master_database_data$t_cell_activation[[1]], check.names = FALSE) # only raw
  all_t_cell_activation_long$col <- paste(row.names(all_t_cell_activation_long), '(t_cell_activation)')
  all_t_cell_activation_long <- pivot_longer(all_t_cell_activation_long, cols = -col, names_to = 'specimen_id', values_to = 'value')
  
  all_data <- bind_rows(list(all_plasma_ab_titer_long, all_plasma_cytokine_concentrations_by_olink_long, all_plasma_cytokine_concentrations_by_legendplex_long,
                             all_pbmc_cell_frequency_long, all_pbmc_gene_expression_long, all_t_cell_polarization_long, all_t_cell_activation_long))
  all_data$value <- as.character(all_data$value)
  
  
  # ==========================================================
  
  
  train_data <- all_data[all_data$specimen_id %in% unique(train_subject_specimen$specimen_id), ]
  train_data <- left_join(train_data, train_subject_specimen[, c('specimen_id', 'timepoint')], by = 'specimen_id')
  train_data$col <- paste(train_data$col, train_data$timepoint)
  train_data <- left_join(train_data, train_subject_specimen[, c('specimen_id', 'subject_id')], by = 'specimen_id')
  train_data <- train_data[!colnames(train_data) %in% c('specimen_id', 'timepoint') ]
  
  train_data <- bind_rows(train_data, train_subject_specimen_long)
  train_data <- pivot_wider(train_data, names_from = subject_id, values_from = value)
  train_data <- train_data %>% 
    column_to_rownames(var = 'col')
  
  challenge_data <- all_data[all_data$specimen_id %in% unique(challenge_subject_specimen$specimen_id), ]
  challenge_data <- left_join(challenge_data, challenge_subject_specimen[, c('specimen_id', 'timepoint')], by = 'specimen_id')
  challenge_data$col <- paste(challenge_data$col, challenge_data$timepoint)
  challenge_data <- left_join(challenge_data, challenge_subject_specimen[, c('specimen_id', 'subject_id')], by = 'specimen_id')
  challenge_data <- challenge_data[!colnames(challenge_data) %in% c('specimen_id', 'timepoint') ]
  
  challenge_data <- bind_rows(challenge_data, challenge_subject_specimen_long)
  challenge_data <- pivot_wider(challenge_data, names_from = subject_id, values_from = value)
  challenge_data <- challenge_data %>%
    column_to_rownames(var = 'col')
  
  rm(master_database_data, all_plasma_ab_titer_long, all_plasma_cytokine_concentrations_by_olink_long, all_plasma_cytokine_concentrations_by_legendplex_long, 
     all_pbmc_cell_frequency_long, all_pbmc_gene_expression_long, all_t_cell_polarization_long, all_t_cell_activation_long, all_data, 
     train_subject_specimen, challenge_subject_specimen, train_subject_specimen_long, challenge_subject_specimen_long)
  
  if (data_type == 1){
    saveRDS(train_data, paste0(dir_RDS_objects, 'train_data raw.RDS'))
    saveRDS(challenge_data, paste0(dir_RDS_objects, 'challenge_data raw.RDS'))
  } else if (data_type == 2){
    saveRDS(train_data, paste0(dir_RDS_objects, 'train_data normalized.RDS'))
    saveRDS(challenge_data, paste0(dir_RDS_objects, 'challenge_data normalized.RDS'))
  } else {
    saveRDS(train_data, paste0(dir_RDS_objects, 'train_data batchCorrected.RDS'))
    saveRDS(challenge_data, paste0(dir_RDS_objects, 'challenge_data batchCorrected.RDS'))
  }
  
}


