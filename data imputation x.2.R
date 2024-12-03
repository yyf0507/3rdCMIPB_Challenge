# library(mice)
library(softImpute)
library(fastDummies)
library(dplyr)
library(stats)

seed = 1e2
set.seed(seed)

dir_RDS_objects <-  "C:/CMI-PB/data/output/"

extract_numbers <- function(x) {
  numbers <- gsub(".*?(-?\\d+)$", "\\1", x)
  if (grepl("-?\\d+", numbers)){
    return(as.numeric(numbers))
  } 
  else {
    return(0)
  }
}

y_names1 <- c('IgG_PT (plasma_ab_titer) 14', 'Monocytes (pbmc_cell_frequency) 1', 'ENSG00000277632.1 (pbmc_gene_expression) 3')
y_names2 <- c('IgG_PT (plasma_ab_titer) 0', 'Monocytes (pbmc_cell_frequency) 0', 'ENSG00000277632.1 (pbmc_gene_expression) 0')
y_names3 <- c('IgG_PT (plasma_ab_titer) 14-0 fc', 'Monocytes (pbmc_cell_frequency) 1-0 fc', 'ENSG00000277632.1 (pbmc_gene_expression) 3-0 fc')

for (task in c(1, 2, 3)){
  for (data_model in c(1, 2, 3, 4)){
    for (data_type in c('raw', 'normalized', 'batchCorrected')){
      train_x <- readRDS(sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2020-2022.RDS', task + 0.1, data_model, data_type))
      challenge_x <- readRDS(sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2023.RDS', task + 0.1, data_model, data_type))
      
      train_data <- readRDS(sprintf('%strain_data %s.RDS', dir_RDS_objects, data_type))
      train_data[train_data == 'NaN'] = NA
      rownames(train_data) <- gsub(' -15', ' -14', rownames(train_data))
      
      train_y <- sapply(train_data[y_names1[task], rownames(train_x)], as.numeric) / sapply(train_data[y_names2[task], rownames(train_x)], as.numeric)
      train_y <- t(t(train_y[!is.na(train_y)]))
      colnames(train_y) <- c(y_names3[task])
      
      common_subjects <- intersect(rownames(train_x), rownames(train_y))
      train_x <- train_x[common_subjects, ]
      train_y <- as.matrix(train_y[common_subjects, ])
      colnames(train_y) <- c(y_names3[task])
      
      train_dataset <- train_data['dataset (subject_info)', common_subjects]
      
      # train_x[, paste('reciprocal', y_names2[task])] = 1 / train_x[, y_names2[task]]
      # challenge_x[, paste('reciprocal', y_names2[task])] = 1 / challenge_x[, y_names2[task]]
      
      saveRDS(train_x, sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2020-2022.RDS', task + 0.2, data_model, data_type))
      saveRDS(challenge_x, sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2023.RDS', task + 0.2, data_model, data_type))
      saveRDS(train_y, sprintf('C:/CMI-PB/models/task %.1f/data %i/y %s 2020-2022.RDS', task + 0.2, data_model, data_type))
      
      # ====================================================================================================================
      
      train_x <- bind_cols(train_x, t(train_dataset))
      
      train_x_val <- train_x[train_x[, 'dataset (subject_info)'] %in% c('2020_dataset', '2021_dataset'), colnames(train_x)[!colnames(train_x) %in% 'dataset (subject_info)']]
      train_y_val <- train_y[train_x[, 'dataset (subject_info)'] %in% c('2020_dataset', '2021_dataset')]
      train_y_val <- as.matrix(train_y_val)
      colnames(train_y_val) <- c(y_names3[task])
      
      challenge_x_val <- train_x[train_x[, 'dataset (subject_info)'] == '2022_dataset', colnames(train_x)[!colnames(train_x) %in% 'dataset (subject_info)']]
      challenge_y_val <- train_y[train_x[, 'dataset (subject_info)'] == '2022_dataset']
      challenge_y_val <- as.matrix(challenge_y_val)
      colnames(challenge_y_val) <- c(y_names3[task])
      
      saveRDS(train_x_val, sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2020-2021.RDS', task + 0.2, data_model, data_type))
      saveRDS(challenge_x_val, sprintf('C:/CMI-PB/models/task %.1f/data %i/x %s 2022.RDS', task + 0.2, data_model, data_type))
      saveRDS(train_y_val, sprintf('C:/CMI-PB/models/task %.1f/data %i/y %s 2020-2021.RDS', task + 0.2, data_model, data_type))
      saveRDS(challenge_y_val, sprintf('C:/CMI-PB/models/task %.1f/data %i/y %s 2022.RDS', task + 0.2, data_model, data_type))
    }
  }
}

