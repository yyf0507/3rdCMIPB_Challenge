# library(mice)
library(softImpute)
library(fastDummies)
library(dplyr)
library(stats)

seed = 1e2
set.seed(seed)

data_model <- 3

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
    
    for (data_source in c('plasma_ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
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
    
    train_x_imputed <- data.frame(train_x[, grepl('subject_info', datanames)], check.names = FALSE)
    train_x_input <- train_x
    
    for (data_source in c('ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
      for (t in c(-30, -14, 0)){
        selected_cols <- datanames[(grepl(data_source, datanames)) & (sapply(datanames, extract_numbers) == t)]
        
        if (length(selected_cols) == 0){
          print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, 0))
          next
        }
        
        train_x_input_1 <- as.matrix(train_x_input[, selected_cols])
        
        # train_x_input_1 <- bind_cols(train_x_input_1, train_y)
        # train_x_input_1 <- bind_cols(train_x_input_1, 0)
        # train_x_input_1 <- biScale(train_x_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)[, -ncol(train_x_input_1)]
        
        # train_x_input_1 <- biScale(train_x_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
        
        print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, ncol(train_x_input_1)))
        colnames(train_x_input_1) <- make.names(1: ncol(train_x_input_1))
    
        train_x_imp <- softImpute(train_x_input_1, rank.max = 2, lambda = 0.1, maxit = 1000)
        
        train_x_completed <- complete(train_x_input_1, train_x_imp)
        # train_x_completed <- complete(train_x_input_1, train_x_imp)[, -ncol(train_x_input_1)]
        
        colnames(train_x_completed) <- selected_cols
        train_x_imputed <- cbind(train_x_imputed, train_x_completed)
      }
      
    }
    
    
    
    challenge_x_imputed <- data.frame(challenge_x[, grepl('subject_info', datanames)], check.names = FALSE)
    challenge_x_input <- challenge_x
    
    for (data_source in c('ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
      for (t in c(-30, -14, 0)){
        selected_cols <- datanames[(grepl(data_source, datanames)) & (sapply(datanames, extract_numbers) == t)]
        
        if (length(selected_cols) == 0){
          print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, 0))
          next
        }
        
        challenge_x_input_1 <- challenge_x_input[, selected_cols]
        # filled <- as.data.frame(rowMeans(challenge_x_input_1))
        # filled[is.na(filled)] <- colMeans(filled, na.rm = TRUE)
        # challenge_x_input_1 <- as.matrix(bind_cols(challenge_x_input_1, filled))
        challenge_x_input_1 <- as.matrix(challenge_x_input_1)
        # challenge_x_input_1 <- biScale(challenge_x_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)[, -ncol(challenge_x_input_1)]
        
        print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, ncol(challenge_x_input_1)))
        colnames(challenge_x_input_1) <- make.names(1: ncol(challenge_x_input_1))
    
        challenge_x_imp <- softImpute(challenge_x_input_1, rank.max = 2, lambda = 0.1, maxit = 1000)
        challenge_x_completed <- complete(challenge_x_input_1, challenge_x_imp)
        colnames(challenge_x_completed) <- selected_cols
        challenge_x_imputed <- cbind(challenge_x_imputed, challenge_x_completed)
      }
      
    }
    
    
    
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
    
    
    train_x_val_imputed <- data.frame(train_x_val[, grepl('subject_info', datanames)], check.names = FALSE)
    train_x_val_input <- train_x_val
    
    for (data_source in c('ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
      for (t in c(-30, -14, 0)){
        selected_cols <- datanames[(grepl(data_source, datanames)) & (sapply(datanames, extract_numbers) == t)]
        
        if (length(selected_cols) == 0){
          print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, 0))
          next
        }
        
        train_x_val_input_1 <- as.matrix(train_x_val_input[, selected_cols])
        
        # train_x_val_input_1 <- bind_cols(train_x_val_input_1, 0)
        # train_x_val_input_1 <- biScale(train_x_val_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)[, -ncol(train_x_val_input_1)]
        # train_x_val_input_1 <- bind_cols(train_x_val_input_1, train_y_val)
        # train_x_val_input_1 <- biScale(train_x_val_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
        
        print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, ncol(train_x_val_input_1)))
        colnames(train_x_val_input_1) <- make.names(1: ncol(train_x_val_input_1))
        
        train_x_val_imp <- softImpute(train_x_val_input_1, rank.max = 2, lambda = 0.1, maxit = 1000)
        
        train_x_val_completed <- complete(train_x_val_input_1, train_x_val_imp)
        # train_x_val_completed <- complete(train_x_val_input_1, train_x_val_imp)[, -ncol(train_x_val_input_1)]
        
        colnames(train_x_val_completed) <- selected_cols
        train_x_val_imputed <- cbind(train_x_val_imputed, train_x_val_completed)
      }
      
    }
    
    
    
    challenge_x_val_imputed <- data.frame(challenge_x_val[, grepl('subject_info', datanames)], check.names = FALSE)
    challenge_x_val_input <- challenge_x_val
    
    for (data_source in c('ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
      for (t in c(-30, -14, 0)){
        selected_cols <- datanames[(grepl(data_source, datanames)) & (sapply(datanames, extract_numbers) == t)]
        
        if (length(selected_cols) == 0){
          print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, 0))
          next
        }
        
        challenge_x_val_input_1 <- as.matrix(challenge_x_val_input[, selected_cols])
        # filled <- as.data.frame(rowMeans(challenge_x_val_input_1))
        # filled[is.na(filled)] <- colMeans(filled, na.rm = TRUE)
        # challenge_x_val_input_1 <- as.matrix(bind_cols(challenge_x_val_input_1, filled))
        # challenge_x_val_input_1 <- biScale(challenge_x_val_input_1, maxit = 1000, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
        
        print(sprintf('%s on day %i imputation...   no. of columns: %s', data_source, t, ncol(challenge_x_val_input_1)))
        colnames(challenge_x_val_input_1) <- make.names(1: ncol(challenge_x_val_input_1))
        
        # ini <- mice(challenge_x_input_1, maxit = 0, seed = 100)
        
        # fx <- flux(challenge_x_input_1)
        # outlist1 <- row.names(fx)[fx$outflux < 0.5]
        # outlist2 <- as.character(ini$loggedEvents[, "out"])
        # 
        # outlist <- unique(c(outlist1, outlist2))
        # challenge_x_filtered <- challenge_x_input_1[, !names(challenge_x_input_1) %in% outlist]
        # pred <- quickpred(challenge_x_filtered, mincor = 0.0, method = 'pearson')
        
        # challenge_x_imp <- mice(challenge_x_filtered, m = 5, maxit = 10, pred = pred, seed = 100)
        # challenge_x_imp <- mice(challenge_x_input_1, m = 5, maxit = 10, seed = 100, method = 'pmm')
        
        # challenge_x_completed <- complete(challenge_x_imp, 1)
        
        challenge_x_val_imp <- softImpute(challenge_x_val_input_1, rank.max = 2, lambda = 0.1, maxit = 1000)
        # challenge_x_val_completed <- complete(challenge_x_val_input_1, challenge_x_val_imp)[, -ncol(challenge_x_val_input_1)]
        challenge_x_val_completed <- complete(challenge_x_val_input_1, challenge_x_val_imp)
        colnames(challenge_x_val_completed) <- selected_cols
        challenge_x_val_imputed <- cbind(challenge_x_val_imputed, challenge_x_val_completed)
      }
      
    }
    
    
    
    saveRDS(train_x_val_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_x_val_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2022.RDS', task, data_model, data_type))
    saveRDS(train_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2022.RDS', task, data_model, data_type))
    
    rm(train_x_val_input, challenge_x_val_input, train_x_val_imp, challenge_x_val_imp)
  }
}