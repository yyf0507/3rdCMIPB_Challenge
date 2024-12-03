# library(mice)
library(softImpute)
library(fastDummies)
library(dplyr)
library(stats)

seed = 1e2
set.seed(seed)

data_model <- 2

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

for (task in c(1, 2, 3)){
  for (data_type in c('raw', 'normalized', 'batchCorrected')){
    train_data <- readRDS(sprintf('%strain_data %s.RDS', dir_RDS_objects, data_type))
    train_data[train_data == 'NaN'] = NA
    rownames(train_data) <- gsub(' -15', ' -14', rownames(train_data))
    
    train_y <- train_data[y_names[task], ]
    train_y <- t(train_y[, !is.na(train_y[1, ])])
    
    challenge_data <- readRDS(sprintf('%schallenge_data %s.RDS', dir_RDS_objects, data_type))
    challenge_x <- challenge_data[sapply(rownames(challenge_data), extract_numbers) <= 0, ]
    
    train_x <- train_data[(sapply(rownames(train_data), extract_numbers) <= 0) & (sapply(rownames(train_data), extract_numbers) > -10) & (rownames(train_data) %in% rownames(challenge_data)), rownames(train_y)]
    train_x <- train_x[, colSums(is.na(train_x)) / nrow(train_x) <= 0]
    
    train_y <- train_y[colnames(train_x), ]
    
    for (data_source in c('subject_info', 'ab_titer', 'cell_frequency', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')){
    
      train_x_selected <- as.matrix(train_x[grepl(data_source, rownames(train_x)), ])
      print(sprintf('%s: %i    missing proportion: %.4f', data_source, nrow(train_x_selected), mean(is.na(train_x_selected))))
    }
    
    common_cols <- intersect(rownames(train_x), rownames(challenge_x))
    train_dataset <- train_data['dataset (subject_info)', names(train_y)]
    train_x <- train_x[common_cols, ]
    challenge_x <- challenge_x[common_cols, ]
    
    rm(train_data, challenge_data, common_cols)
    
    # ====================================================================================================================
    train_x <- as.data.frame(t(train_x))
    train_ind <- rownames(train_x)
    
    challenge_x <- as.data.frame(t(challenge_x))
    challenge_ind <- rownames(challenge_x)
    
    print(sum(rownames(train_x) == names(train_y)) / length(train_y))
    
    train_y <- apply(as.matrix(train_y), 2, as.numeric)
    rownames(train_y) <- rownames(train_x)
    
    x <- bind_rows(train_x, challenge_x)
    x <- dummy_cols(
      x,
      select_columns = c('infancy_vac (subject_info)', 'biological_sex (subject_info)', 'ethnicity (subject_info)', 'race (subject_info)'),
      remove_most_frequent_dummy = FALSE,
      remove_selected_columns = TRUE 
    )
    
    x <- as.data.frame(apply(x, 2, as.numeric))
    rownames(x) <- c(train_ind, challenge_ind)
    x_col <- colnames(x)
    
    train_x <- x[train_ind, ]
    challenge_x <- x[challenge_ind, ]
    
    rm(x)
    
    # ====================================================================================================================
    
    # colnames(challenge_x) <- make.names(1: ncol(challenge_x))
    # ini <- mice(challenge_x, maxit = 0, seed = 100)
    # 
    # fx <- flux(challenge_x)
    # outlist1 <- row.names(fx)[fx$outflux < 0.5]
    # outlist2 <- as.character(ini$loggedEvents[, "out"])
    # 
    # outlist <- unique(c(outlist1, outlist2))
    # challenge_x2 <- challenge_x[, !names(challenge_x) %in% outlist]
    # pred <- quickpred(challenge_x2, mincor = 0.0, method = 'pearson')
    # 
    # challenge_x_imp <- mice(challenge_x2, m = 5, maxit = 10, pred = pred, seed = 100)
    
    challenge_x_input <- as.matrix(challenge_x)
    challenge_x_input <- biScale(challenge_x_input, maxit = 100, col.scale = TRUE, row.scale = FALSE, trace = TRUE)
    challenge_x_imp <- softImpute(challenge_x_input, rank.max = 2, lambda = 0.1, maxit = 500)
    challenge_x_imputed <- as.data.frame(complete(challenge_x, challenge_x_imp))
    
    saveRDS(train_x, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2020-2022.RDS', task, data_model, data_type))
    saveRDS(challenge_x_imputed, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2023.RDS', task, data_model, data_type))
    saveRDS(train_y, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2020-2022.RDS', task, data_model, data_type))
    
    
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
    
    
    
    
    # train_data <- readRDS(sprintf('%strain_data %s.RDS', dir_RDS_objects, data_type))
    # train_data <- train_data[, train_data['dataset (subject_info)', ] %in% c('2020_dataset', '2021_dataset', '2022_dataset')]
    # rownames(train_data) <- gsub(' -15', ' -14', rownames(train_data))
    # 
    # train_y_val <- train_data['ENSG00000277632.1 (pbmc_gene_expression) 3', train_data['dataset (subject_info)', ] %in% c('2020_dataset', '2021_dataset')]
    # train_y_val <- t(train_y_val[, t(!is.na(train_y_val)) ])
    # 
    # challenge_y_val <- train_data['ENSG00000277632.1 (pbmc_gene_expression) 3', train_data['dataset (subject_info)', ] == '2022_dataset']
    # challenge_y_val <- t(challenge_y_val[, t(!is.na(challenge_y_val)) ])
    
    
    saveRDS(train_x_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_x_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/x %s 2022.RDS', task, data_model, data_type))
    saveRDS(train_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2020-2021.RDS', task, data_model, data_type))
    saveRDS(challenge_y_val, sprintf('C:/CMI-PB/models/task %i.1/data %i/y %s 2022.RDS', task, data_model, data_type))
  
  }
}
