library(dplyr)
# library(stats)
library(MOFA2)

seed = 1e2
set.seed(seed)

extract_numbers <- function(x) {
  numbers <- gsub(".*?(-?\\d+)$", "\\1", x)
  if (grepl("-?\\d+", numbers)){
    return(as.numeric(numbers))
  } 
  else {
    return(0)
  }
}

L1penalty_linear_train <- function(x, y, n_factors) {
  library(glmnet)
  
  formula <- as.formula(paste('y ~ ', paste(colnames(x), collapse = '+')))
  penalty_factor <- c(rep(0, n_factors), rep(1, ncol(x) - n_factors))
  model <- glmnet(x, y, alpha = 1, lambda = 0.2, penalty.factor = penalty_factor)
  
  return(model)
}

tasks <- c('1.1', '1.2', '2.1', '2.2', '3.1', '3.2')
data_types <- c('raw', 'normalized', 'batchCorrected')
y_names <- c('IgG-PT-D14-titer-Rank', 'IgG-PT-D14-FC-Rank', 'Monocytes-D1-Rank', 'Monocytes-D1-FC-Rank', 'CCL3-D3-Rank', 'CCL3-D3-FC-Rank')
y_names_baseline <- c('IgG_PT (plasma_ab_titer) 0', 'Monocytes (pbmc_cell_frequency) 0', 'ENSG00000277632.1 (pbmc_gene_expression) 0')

data_model <- 1
data_type <- 'raw'
task <- '1.2'

for (data_model in c(1, 2, 3, 4)){
  
  performance <- data.frame(matrix(nrow = length(data_types) * 2, ncol = length(tasks)))
  colnames(performance) <- tasks
  rownames(performance) <- c(paste('2020-2021 train', data_types), paste('2022 test', data_types))
  
  for (data_type in c('raw', 'normalized', 'batchCorrected')){
    
    output <- readRDS("C:/CMI-PB/data/output template.RDS")
    
    for (task in tasks){
      if (substr(task, 1, 1) == '1'){
        data_sources <- c('ab_titer', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation', 'cell_frequency')
      }
      else if (substr(task, 1, 1) == '2'){
        data_sources <- c('plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation', 'cell_frequency')
      }
      else {
        data_sources <- c('plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation')
      }
      
      dir_RDS_objects <-  sprintf("C:/CMI-PB/models/task %s/data %i/", task, data_model)
      
      x_2020_2022 <- readRDS(sprintf('%sx %s 2020-2022.RDS', dir_RDS_objects, data_type))
      y_2020_2022 <- readRDS(sprintf('%sy %s 2020-2022.RDS', dir_RDS_objects, data_type))
      
      x_2020_2021 <- readRDS(sprintf('%sx %s 2020-2021.RDS', dir_RDS_objects, data_type))
      y_2020_2021 <- readRDS(sprintf('%sy %s 2020-2021.RDS', dir_RDS_objects, data_type))
      
      x_2022 <- readRDS(sprintf('%sx %s 2022.RDS', dir_RDS_objects, data_type))
      y_2022 <- readRDS(sprintf('%sy %s 2022.RDS', dir_RDS_objects, data_type))
      
      x_2023 <- readRDS(sprintf('%sx %s 2023.RDS', dir_RDS_objects, data_type))
      
      if (data_type == 'raw' & data_model %in% c(1, 2)){
        for (data_source in data_sources){
          x_2020_2022[, grepl(data_source, colnames(x_2020_2022))] <- log(x_2020_2022[, grepl(data_source, colnames(x_2020_2022))])
          x_2020_2021[, grepl(data_source, colnames(x_2020_2021))] <- log(x_2020_2021[, grepl(data_source, colnames(x_2020_2021))])
          x_2022[, grepl(data_source, colnames(x_2022))] <- log(x_2022[, grepl(data_source, colnames(x_2022))])
          x_2023[, grepl(data_source, colnames(x_2023))] <- log(x_2023[, grepl(data_source, colnames(x_2023))])
        }
      }
      
      if (data_type == 'raw' & substr(task, 3, 3) == '1'){
        y_2020_2022 <- log(y_2020_2022)
        y_2020_2021 <- log(y_2020_2021)
        y_2022 <- log(y_2022)
      }
      
      x_col <- colnames(x_2020_2022)
      
      # ====================================================================================================================

      colnames(x_2020_2022) <- paste0('feature_', 1: ncol(x_2020_2022))
      colnames(x_2020_2021) <- paste0('feature_', 1: ncol(x_2020_2021))
      colnames(x_2022) <- paste0('feature_', 1: ncol(x_2022))
      colnames(x_2023) <- paste0('feature_', 1: ncol(x_2023))
      colnames(y_2020_2022) <- 'y'
      colnames(y_2020_2021) <- 'y'
      colnames(y_2022) <- 'y'
      
      data_source_selected <- c()
      x_2020_2022_list <- list()
      
      for (data_source in data_sources){
        
        x_2020_2022_selected <- scale(as.matrix(x_2020_2022[, grepl(data_source, x_col)]))
        rownames(x_2020_2022_selected) <- paste0('sample_', 1: nrow(x_2020_2022_selected))
        
        if(ncol(x_2020_2022_selected) == 0){
          next
        }
        else{
          data_source_selected <- c(data_source_selected, data_source)
        }
        
        x_2020_2022_list[[data_source]] <- t(x_2020_2022_selected)
        
      }
      
      names(x_2020_2022_list) <- paste0('view_', 1: length(x_2020_2022_list))
      MOFAobject <- create_mofa(x_2020_2022_list)
      

      model_opts <- get_default_model_options(MOFAobject)
      model_opts$num_factors <- 6
      MOFAobject <- prepare_mofa(
        object = MOFAobject,
        model_options = model_opts
      )

      
      MOFAobject.trained <- run_mofa(MOFAobject, sprintf("C:/CMI-PB/models/MOFA/MOFA models/model x_2020_2022 %s %s %i.hdf5", task, data_type, data_model), use_basilisk = TRUE)
      
      rm(MOFAobject, MOFAobject.trained, model_opts)
      
      # =============================================================================================
      
      model <- load_model(sprintf("C:/CMI-PB/models/MOFA/MOFA models/model x_2020_2022 %s %s %i.hdf5", task, data_type, data_model))
      factors <- get_factors(model)[[1]]
      rownames(factors) <- rownames(x_2020_2022)
      
      train_factors <- factors[rownames(x_2020_2021), ]
      test_factors <- factors[rownames(x_2022), ]
      
      if (substr(task, 3, 3) == '1'){
        inv_x_2020_2021 <- as.data.frame(x_2020_2021[, colnames(x_2020_2021)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2020_2021) <- 'y0'
        rownames(inv_x_2020_2021) <- rownames(x_2020_2021)
        inv_x_2022 <- as.data.frame(x_2022[, colnames(x_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2022) <- 'y0'
        rownames(inv_x_2022) <- rownames(x_2022)
        
        train_factors <- bind_cols(train_factors, inv_x_2020_2021)
        test_factors <- bind_cols(test_factors, inv_x_2022)
        
        n_factors <- ncol(factors) + 1
      }
      else{
        inv_x_2020_2021 <- as.data.frame(1 / exp(x_2020_2021[, colnames(x_2020_2021)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2020_2021) <- 'reciprocal.y'
        rownames(inv_x_2020_2021) <- rownames(x_2020_2021)
        inv_x_2022 <- as.data.frame(1 / exp(x_2022[, colnames(x_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2022) <- 'reciprocal.y'
        rownames(inv_x_2022) <- rownames(x_2022)
        
        train_factors <- bind_cols(train_factors, inv_x_2020_2021)
        test_factors <- bind_cols(test_factors, inv_x_2022)
        
        n_factors <- ncol(factors) + 1
      }
      
      x_2020_2021_input <- bind_cols(train_factors, x_2020_2021[, grepl('subject_info', x_col)])
      x_2022_input <- bind_cols(test_factors, x_2022[, grepl('subject_info', x_col)])
      
      model_val <- L1penalty_linear_train(x_2020_2021_input, y_2020_2021, n_factors = n_factors)
      
      y_2020_2021_pred <- stats::predict(model_val, newx = as.matrix(x_2020_2021_input))
      y_2022_pred <- stats::predict(model_val, newx = as.matrix(x_2022_input))
      
      train_cor <- cor.test(y_2020_2021, y_2020_2021_pred, method = 'spearman')
      test_cor <- cor.test(y_2022, y_2022_pred, method = 'spearman')
      print(sprintf('%s %s %s', task, data_model, data_type))
      print(sprintf('2020-2021 train spearman cor: %.4f', train_cor$estimate))
      print(sprintf('2022 test spearman cor: %.4f', test_cor$estimate))
      
      performance[c(paste('2020-2021 train', data_type), paste('2022 test', data_type)), task] <- c(ifelse(train_cor$p.value < 0.05, train_cor$estimate, NA), 
                                                                                                    ifelse(test_cor$p.value < 0.05, test_cor$estimate, NA))
      plot
      
      
      # ================================================================================================================
      
      x_2020_2023 <- bind_rows(x_2020_2022, x_2023)
      x_2020_2023_list <- list()
      
      for (data_source in data_sources){
        
        x_2020_2023_selected <- scale(as.matrix(x_2020_2023[, grepl(data_source, x_col)]))
        rownames(x_2020_2023_selected) <- paste0('sample_', 1: nrow(x_2020_2023_selected))
        
        if(ncol(x_2020_2023_selected) == 0){
          next
        }
        else{
          data_source_selected <- c(data_source_selected, data_source)
        }
        
        x_2020_2023_list[[data_source]] <- t(x_2020_2023_selected)
        
      }
      
      names(x_2020_2023_list) <- paste0('view_', 1: length(x_2020_2023_list))
      MOFAobject <- create_mofa(x_2020_2023_list)
      
      if (substr(task, 1, 1) == '5'){
        model_opts <- get_default_model_options(MOFAobject)
        MOFAobject <- prepare_mofa(
          object = MOFAobject,
          model_options = model_opts
        )
      }
      
      else{
        model_opts <- get_default_model_options(MOFAobject)
        model_opts$num_factors <- 6
        MOFAobject <- prepare_mofa(
          object = MOFAobject,
          model_options = model_opts
        )
      }
      
      MOFAobject.trained <- run_mofa(MOFAobject, sprintf("C:/CMI-PB/models/MOFA/MOFA models/model x_2020_2023 %s %s %i.hdf5", task, data_type, data_model), use_basilisk = TRUE)
      
      rm(MOFAobject, MOFAobject.trained, model_opts)
      
      # =================================================================================================================
      
      model <- load_model(sprintf("C:/CMI-PB/models/MOFA/MOFA models/model x_2020_2023 %s %s %i.hdf5", task, data_type, data_model))
      factors <- get_factors(model)[[1]]
      rownames(factors) <- rownames(x_2020_2023)
      
      train_factors <- factors[rownames(x_2020_2022), ]
      test_factors <- factors[rownames(x_2023), ]
      
      if (substr(task, 3, 3) == '1'){
        inv_x_2020_2022 <- as.data.frame(x_2020_2022[, colnames(x_2020_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2020_2022) <- 'y0'
        rownames(inv_x_2020_2022) <- rownames(x_2020_2022)
        inv_x_2023 <- as.data.frame(x_2023[, colnames(x_2023)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2023) <- 'y0'
        rownames(inv_x_2023) <- rownames(x_2023)
        
        train_factors <- bind_cols(train_factors, inv_x_2020_2022)
        test_factors <- bind_cols(test_factors, inv_x_2023)
        
        n_factors <- ncol(factors) + 1
      }
      else{
        inv_x_2020_2022 <- as.data.frame(1 / exp(x_2020_2022[, colnames(x_2020_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2020_2022) <- 'reciprocal.y'
        rownames(inv_x_2020_2022) <- rownames(x_2020_2022)
        inv_x_2023 <- as.data.frame(1 / exp(x_2023[, colnames(x_2023)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2023) <- 'reciprocal.y'
        rownames(inv_x_2023) <- rownames(x_2023)
        
        train_factors <- bind_cols(train_factors, inv_x_2020_2022)
        test_factors <- bind_cols(test_factors, inv_x_2023)
        
        n_factors <- ncol(factors) + 1
      }
      
      x_2020_2022_input <- bind_cols(train_factors, x_2020_2022[, grepl('subject_info', x_col)])
      x_2023_input <- bind_cols(test_factors, x_2023[, grepl('subject_info', x_col)])
      
      model_predict <- L1penalty_linear_train(x_2020_2022_input, y_2020_2022, n_factors = n_factors)
      
      y_2020_2022_pred <- stats::predict(model_predict, newx = as.matrix(x_2020_2022_input))
      
      train_cor <- cor(y_2020_2022, y_2020_2022_pred, method = 'spearman')
      print(sprintf('2020-2022 train spearman cor: %.4f', cor(y_2020_2022, y_2020_2022_pred, method = 'spearman')))
      
      challenge_y_pred <- rank(-stats::predict(model_predict, newx = as.matrix(x_2023_input)))

      output[, y_names[which(tasks == task)]] <- challenge_y_pred

    }

    write.csv(output[order(output$SubjectID), ], sprintf('C:/CMI-PB/models/MOFA/prediction results/challenge_y predicted from data %i %s.csv', data_model, data_type), row.names = FALSE)
    
  }
  
  write.csv(performance, sprintf("C:/CMI-PB/models/MOFA/MOFA performance %i.csv", data_model))
  
}


