# library(mice)
library(softImpute)
library(fastDummies)
library(dplyr)
library(stats)

seed = 1e2
set.seed(seed)

data_model <- 4

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

y_names <- c('IgG_PT (plasma_ab_titer) 14', 'Monocytes (pbmc_cell_frequency) 1', 'ENSG00000277632.1 (pbmc_gene_expression) 3')
data_source_selected <- c('subject_info')

for (task in c(1, 2, 3)){
  for (data_type in c('raw', 'normalized', 'batchCorrected')){
    train_data <- readRDS(sprintf('%strain_data %s.RDS', dir_RDS_objects, data_type))
    train_data[train_data == 'NaN'] = NA
    rownames(train_data) <- gsub(' -15', ' -14', rownames(train_data))
    
    train_y <- train_data[y_names[task], ]
    train_y <- t(train_y[, !is.na(train_y[1, ])])
    
    challenge_data <- readRDS(sprintf('%schallenge_data %s.RDS', dir_RDS_objects, data_type))
    challenge_x <- challenge_data[sapply(rownames(challenge_data), extract_numbers) <= 0, ]
    
    # train_x <- train_data[(sapply(rownames(train_data), extract_numbers) <= 0) & (sapply(rownames(train_data), extract_numbers) > -10) & (rownames(train_data) %in% rownames(challenge_data)), rownames(train_y)]
    train_x <- train_data[(sapply(rownames(train_data), extract_numbers) <= 0) & (rownames(train_data) %in% rownames(challenge_data)), rownames(train_y)]
    
    for (data_source in c('subject_info', 'ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                         'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
      
      train_x1 <- as.matrix(train_x[grepl(data_source, rownames(train_x)), ])
      for (t in c(-30, -14, 0)){
        train_x2 <- train_x1[sapply(rownames(train_x1), extract_numbers) == t, ]
        if (length(train_x2) == 0){
          next
        }
        
        if (mean(is.na(train_x2)) <= 0.5){
          data_source_selected <- c(data_source_selected, sprintf('%s\\) %i', data_source, t))
        }
        print(sprintf('%s on day %i: %i    missing proportion: %.4f', data_source, t, nrow(train_x2), mean(is.na(train_x2))))
      }
    }
    
    data_source_selected <- paste0(data_source_selected, collapse = '|')
    train_x <- train_x[grep(data_source_selected, rownames(train_x)), ]
    challenge_x <- challenge_x[grep(data_source_selected, rownames(challenge_x)), ]
    
    train_dataset <- train_data['dataset (subject_info)', rownames(train_y)]
    
    rm(train_data, challenge_data, train_x1, train_x2)
    
    # ====================================================================================================================
    
    train_x <- as.data.frame(t(train_x))
    train_subjects <- rownames(train_x)
    
    challenge_x <- as.data.frame(t(challenge_x))
    challenge_subjects <- rownames(challenge_x)
    
    print(sum(rownames(train_x) == rownames(train_y)) / length(train_y))
    
    train_y <- apply(train_y, 2, as.numeric)
    rownames(train_y) <- rownames(train_x)
    
    x <- bind_rows(train_x, challenge_x)
    x <- dummy_cols(
      x,
      select_columns = c('infancy_vac (subject_info)', 'biological_sex (subject_info)', 'ethnicity (subject_info)', 'race (subject_info)'),
      remove_most_frequent_dummy = FALSE,
      remove_selected_columns = TRUE 
    )
    
    x <- as.data.frame(apply(x, 2, as.numeric))
    rownames(x) <- c(train_subjects, challenge_subjects)
    x_col <- colnames(x)
    
    train_x <- x[train_subjects, ]
    challenge_x <- x[challenge_subjects, ]
    
    datanames <- colnames(train_x)
    
    rm(x)
    
    # ====================================================================================================================
    
    train_x_input <- as.matrix(train_x)
    # train_x_input <- bind_cols(train_x_input, train_y)
    rownames(train_x_input) <- rownames(train_x)
    colnames(train_x_input) <- make.names(1: ncol(train_x_input))
    # train_x_input <- biScale(train_x_input, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
    
    # train_x_imp <- softImpute(train_x_input, rank.max = ceiling(min(dim(train_x_input)) / 2), lambda = 0.1)
    train_x_imp <- softImpute(train_x_input, rank.max = 2, lambda = 0.1, maxit = 1000)
    # train_x_imputed <- as.data.frame(complete(train_x_input, train_x_imp)[, -ncol(train_x_input)])
    train_x_imputed <- as.data.frame(complete(train_x_input, train_x_imp))
    colnames(train_x_imputed) <- datanames
    
    
    
    challenge_x_input <- as.matrix(challenge_x)
    rownames(challenge_x_input) <- rownames(challenge_x)
    colnames(challenge_x_input) <- make.names(1: ncol(challenge_x_input))
    # challenge_x_input <- biScale(challenge_x_input, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
    
    # challenge_x_imp <- softImpute(challenge_x_input, rank.max = ceiling(min(dim(challenge_x_input)) / 2), lambda = 0.1)
    challenge_x_imp <- softImpute(challenge_x_input, rank.max = 2, lambda = 0.1, maxit = 1000)
    challenge_x_imputed <- as.data.frame(complete(challenge_x_input, challenge_x_imp))
    colnames(challenge_x_imputed) <- datanames
    
    
    
    saveRDS(train_x_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2020-2022.RDS', task, data_model, data_type))
    saveRDS(challenge_x_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2023.RDS', task, data_model, data_type))
    saveRDS(train_y, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2020-2022.RDS', task, data_model, data_type))
    
    rm(train_x_input, challenge_x_input, train_x_imp, challenge_x_imp)
    # ====================================================================================================================
    
    train_x <- bind_cols(train_x, t(train_dataset))
    
    train_x_val <- train_x[train_x[, 'dataset (subject_info)'] %in% c('2020_dataset', '2021_dataset'), colnames(train_x)[!colnames(train_x) %in% 'dataset (subject_info)']]
    train_y_val <- train_y[train_x[, 'dataset (subject_info)'] %in% c('2020_dataset', '2021_dataset'), ]
    train_y_val <- as.matrix(train_y_val)
    colnames(train_y_val) <- c(y_names[task])
    
    challenge_x_val <- train_x[train_x[, 'dataset (subject_info)'] == '2022_dataset', colnames(train_x)[!colnames(train_x) %in% 'dataset (subject_info)']]
    challenge_y_val <- train_y[train_x[, 'dataset (subject_info)'] == '2022_dataset', ]
    challenge_y_val <- as.matrix(challenge_y_val)
    colnames(challenge_y_val) <- c(y_names[task])
    
    
    
    train_x_val_input <- as.matrix(train_x_val)
    # train_x_val_input <- bind_cols(train_x_val_input, train_y_val)
    colnames(train_x_val_input) <- make.names(1: ncol(train_x_val_input))
    # train_x_val_input <- biScale(train_x_val_input, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
    
    # train_x_val_imp <- softImpute(train_x_val_input, rank.max = ceiling(min(dim(train_x_val_input)) / 2), lambda = 0.1)
    train_x_val_imp <- softImpute(train_x_val_input, rank.max = 2, lambda = 0.1, maxit = 1000)
    # train_x_val_imputed <- as.data.frame(complete(train_x_val_input, train_x_val_imp)[, -ncol(train_x_val_input)])
    train_x_val_imputed <- as.data.frame(complete(train_x_val_input, train_x_val_imp))
    colnames(train_x_val_imputed) <- datanames
    
    challenge_x_val_input <- as.matrix(challenge_x_val)
    colnames(challenge_x_val_input) <- make.names(1: ncol(challenge_x_val_input))
    # challenge_x_val_input <- biScale(challenge_x_val_input, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
    
    # challenge_x_val_imp <- softImpute(challenge_x_val_input, rank.max = ceiling(min(dim(challenge_x_val_input)) / 2), lambda = 0.1)
    challenge_x_val_imp <- softImpute(challenge_x_val_input, rank.max = 2, lambda = 0.1, maxit = 1000)
    challenge_x_val_imputed <- as.data.frame(complete(challenge_x_val_input, challenge_x_val_imp))
    colnames(challenge_x_val_imputed) <- datanames
    
    
    
    saveRDS(train_x_val_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_x_val_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2022.RDS', task, data_model, data_type))
    saveRDS(train_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2022.RDS', task, data_model, data_type))
    
    rm(train_x_val_input, challenge_x_val_input, train_x_val_imp, challenge_x_val_imp)
  }
}



