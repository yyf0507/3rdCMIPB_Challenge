library(dplyr)
library(stats)
# library(MOFA2)

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

tasks <- c('1.1', '1.2', '2.1', '2.2', '3.1', '3.2')
data_types <- c('raw', 'normalized', 'batchCorrected')
y_names <- c('IgG-PT-D14-titer-Rank', 'IgG-PT-D14-FC-Rank', 'Monocytes-D1-Rank', 'Monocytes-D1-FC-Rank', 'CCL3-D3-Rank', 'CCL3-D3-FC-Rank')
y_names_baseline <- c('IgG_PT (plasma_ab_titer) 0', 'Monocytes (pbmc_cell_frequency) 0', 'ENSG00000277632.1 (pbmc_gene_expression) 0')

data_model <- 1
data_type <- 'raw'
task <- '2.1'

for (data_model in c(1, 3, 4)){
  
  performance <- data.frame(matrix(nrow = length(data_types) * 2, ncol = length(tasks)))
  colnames(performance) <- tasks
  rownames(performance) <- c(paste('2020-2021 train', data_types), paste('2022 test', data_types))
  
  for (data_type in c('raw', 'normalized', 'batchCorrected')){
    
    output <- readRDS("C:/CMI-PB/data/output template.RDS")
    
    for (task in tasks){
      if (substr(task, 1, 1) == '1'){
        data_sources <- c('subject_info', 'ab_titer', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation', 'cell_frequency')
      }
      else if (substr(task, 1, 1) == '2'){
        data_sources <- c('subject_info', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
                          'pbmc_gene_expression', 't_cell_polarization', 't_cell_activation', 'cell_frequency')
      }
      else {
        data_sources <- c('subject_info', 'plasma_cytokine_concentrations_by_olink', 'plasma_cytokine_concentrations_by_legendplex',
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
          if (data_source != 'subject_info'){
            x_2020_2022[, grepl(data_source, colnames(x_2020_2022))] <- log(x_2020_2022[, grepl(data_source, colnames(x_2020_2022))])
            x_2020_2021[, grepl(data_source, colnames(x_2020_2021))] <- log(x_2020_2021[, grepl(data_source, colnames(x_2020_2021))])
            x_2022[, grepl(data_source, colnames(x_2022))] <- log(x_2022[, grepl(data_source, colnames(x_2022))])
            x_2023[, grepl(data_source, colnames(x_2023))] <- log(x_2023[, grepl(data_source, colnames(x_2023))])
          }
        }
      }
      
      if (data_type == 'raw' & substr(task, 3, 3) == '1'){
        y_2020_2022 <- log(y_2020_2022)
        y_2020_2021 <- log(y_2020_2021)
        y_2022 <- log(y_2022)
      }
      
      x_col <- colnames(x_2020_2022)
      
      # ====================================================================================================================
      
      library(glmnet)
      
      selected_var <- paste0('X', which(x_col %in% x_col[grepl('subject_info', x_col)]))
      
      colnames(x_2020_2022) <- make.names(1: ncol(x_2020_2022))
      colnames(x_2020_2021) <- make.names(1: ncol(x_2020_2021))
      colnames(x_2022) <- make.names(1: ncol(x_2022))
      colnames(x_2023) <- make.names(1: ncol(x_2023))
      colnames(y_2020_2022) <- 'y'
      colnames(y_2020_2021) <- 'y'
      colnames(y_2022) <- 'y'
      
      for (data_source in data_sources){
        
        x_2020_2021_selected <- as.matrix(x_2020_2021[, grepl(data_source, x_col)])
        
        if(ncol(x_2020_2021_selected) == 0){
          next
        }
        
        lasso_model <- glmnet(x_2020_2021_selected, y_2020_2021, alpha = 1)
        cv_lasso <- cv.glmnet(x_2020_2021_selected, y_2020_2021, alpha = 1, nfolds = 5)
        best_lambda <- cv_lasso$lambda.min
        
        # print(sprintf('%s ncol: %i, best_lambda: %f', data_source, ncol(x_2020_2021_selected), best_lambda))
        
        coef_lasso <- as.matrix(coef(cv_lasso, s = 'lambda.min'))[-1, , drop = FALSE]
        
        selected_var <- c(selected_var, rownames(coef_lasso)[which(coef_lasso != 0)])
        
      }
      
      # print(x_col[as.numeric(gsub('X', '', selected_var))])
      print(sprintf('no. of var: %i', length(selected_var)))
      
      if (substr(task, 3, 3) == '1'){
        inv_x_2020_2021 <- as.data.frame(x_2020_2021[, colnames(x_2020_2021)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2020_2021) <- 'y0'
        rownames(inv_x_2020_2021) <- rownames(x_2020_2021)
        inv_x_2022 <- as.data.frame(x_2022[, colnames(x_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2022) <- 'y0'
        rownames(inv_x_2022) <- rownames(x_2022)
        
        train <- bind_cols(subset(x_2020_2021, select = selected_var), inv_x_2020_2021, y_2020_2021)
        test <- bind_cols(subset(x_2022, select = selected_var), inv_x_2022, y_2022)
        
        selected_var <- c(selected_var, 'y0')
      }
      else {
        inv_x_2020_2021 <- as.data.frame(1 / exp(x_2020_2021[, colnames(x_2020_2021)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2020_2021) <- 'reciprocal.y'
        rownames(inv_x_2020_2021) <- rownames(x_2020_2021)
        inv_x_2022 <- as.data.frame(1 / exp(x_2022[, colnames(x_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2022) <- 'reciprocal.y'
        rownames(inv_x_2022) <- rownames(x_2022)
        
        train <- bind_cols(subset(x_2020_2021, select = selected_var), inv_x_2020_2021, y_2020_2021)
        test <- bind_cols(subset(x_2022, select = selected_var), inv_x_2022, y_2022)
        
        selected_var <- c(selected_var, 'reciprocal.y')
      }
      
      formula <- as.formula(paste('y ~ ', paste(selected_var, collapse = '+')))
      model_val <- lm(formula = formula, data = train)
      
      y_2020_2021_pred <- stats::predict(model_val, newdata = train)
      y_2022_pred <- stats::predict(model_val, newdata = test)
      
      train_cor <- cor.test(y_2020_2021, y_2020_2021_pred, method = 'spearman')
      test_cor <- cor.test(y_2022, y_2022_pred, method = 'spearman')
      print(sprintf('%s %s %s', task, data_model, data_type))
      print(sprintf('2020-2021 train spearman cor: %.4f', train_cor$estimate))
      print(sprintf('2022 test spearman cor: %.4f', test_cor$estimate))
      
      performance[c(paste('2020-2021 train', data_type), paste('2022 test', data_type)), task] <- c(ifelse(train_cor$p.value < 0.05, train_cor$estimate, NA), 
                                                                                                    ifelse(test_cor$p.value < 0.05, test_cor$estimate, NA))
      
      
      
      # ================================================================================================================
      
      selected_var <- paste0('X', which(x_col %in% x_col[grepl('subject_info', x_col)]))
      
      for (data_source in data_sources){

        x_2020_2022_selected <- as.matrix(x_2020_2022[, grepl(data_source, x_col)])

        if(ncol(x_2020_2022_selected) == 0){
          next
        }

        lasso_model <- glmnet(x_2020_2022_selected, y_2020_2022, alpha = 1)
        cv_lasso <- cv.glmnet(x_2020_2022_selected, y_2020_2022, alpha = 1, nfolds = 5)
        best_lambda <- cv_lasso$lambda.min

        # print(sprintf('%s ncol: %i, best_lambda: %f', data_source, ncol(x_2020_2022_selected), best_lambda))

        coef_lasso <- as.matrix(coef(cv_lasso, s = 'lambda.min'))[-1, , drop = FALSE]

        selected_var <- c(selected_var, rownames(coef_lasso)[which(coef_lasso != 0)])

      }

      # print(x_col[as.numeric(gsub('X', '', selected_var))])
      print(sprintf('no. of var: %i', length(selected_var)))
      
      if (substr(task, 3, 3) == '1'){
        inv_x_2020_2022 <- as.data.frame(x_2020_2022[, colnames(x_2020_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2020_2022) <- 'y0'
        rownames(inv_x_2020_2022) <- rownames(x_2020_2022)
        inv_x_2023 <- as.data.frame(x_2023[, colnames(x_2023)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]])
        colnames(inv_x_2023) <- 'y0'
        rownames(inv_x_2023) <- rownames(x_2023)
        
        train <- bind_cols(subset(x_2020_2022, select = selected_var), inv_x_2020_2022, y_2020_2022)
        test <- bind_cols(subset(x_2023, select = selected_var), inv_x_2023)
        
        selected_var <- c(selected_var, 'y0')
      }
      else {
        inv_x_2020_2022 <- as.data.frame(1 / exp(x_2020_2022[, colnames(x_2020_2022)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2020_2022) <- 'reciprocal.y'
        rownames(inv_x_2020_2022) <- rownames(x_2020_2022)
        inv_x_2023 <- as.data.frame(1 / exp(x_2023[, colnames(x_2023)[which(x_col == y_names_baseline[as.numeric(substr(task, 1, 1))])]]))
        colnames(inv_x_2023) <- 'reciprocal.y'
        rownames(inv_x_2023) <- rownames(x_2023)
        
        train <- bind_cols(subset(x_2020_2022, select = selected_var), inv_x_2020_2022, y_2020_2022)
        test <- bind_cols(subset(x_2023, select = selected_var), inv_x_2023)
        
        selected_var <- c(selected_var, 'reciprocal.y')
      }
      
      formula <- as.formula(paste('y ~ ', paste(selected_var, collapse = '+')))
      model <- lm(formula = formula, data = train)
      
      y_2020_2022_pred <- stats::predict(model, newdata = train)
      train_cor <- cor(y_2020_2022, y_2020_2022_pred, method = 'spearman')
      print(sprintf('2020-2021 train spearman cor: %.4f', train_cor))
      
      challenge_y_pred <- rank(-stats::predict(model, newdata = test))
      output[, y_names[which(tasks == task)]] <- challenge_y_pred

    }
    
    write.csv(output[order(output$SubjectID), ], sprintf('C:/CMI-PB/models/baseline/prediction results/challenge_y predicted from data %i %s.csv', data_model, data_type), row.names = FALSE)
    
  }

  write.csv(performance, sprintf("C:/CMI-PB/models/baseline/baseline performance %i.csv", data_model))
  
}


