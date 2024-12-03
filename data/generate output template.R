
task <- '1.1'
data_model <- 1
data_type <- 'raw'

dir_RDS_objects <-  sprintf("C:/CMI-PB/models/task %s/data %i/", task, data_model)

x_2020_2022 <- readRDS(sprintf('%sx %s 2020-2022.RDS', dir_RDS_objects, data_type))
y_2020_2022 <- readRDS(sprintf('%sy %s 2020-2022.RDS', dir_RDS_objects, data_type))

x_2020_2021 <- readRDS(sprintf('%sx %s 2020-2021.RDS', dir_RDS_objects, data_type))
y_2020_2021 <- readRDS(sprintf('%sy %s 2020-2021.RDS', dir_RDS_objects, data_type))

x_2022 <- readRDS(sprintf('%sx %s 2022.RDS', dir_RDS_objects, data_type))
y_2022 <- readRDS(sprintf('%sy %s 2022.RDS', dir_RDS_objects, data_type))

x_2023 <- readRDS(sprintf('%sx %s 2023.RDS', dir_RDS_objects, data_type))

if (data_type == 'raw' & data_model %in% c(1, 2)){
  x_2020_2022[, grepl('gene_expression', colnames(x_2020_2022))] <- log(x_2020_2022[, grepl('gene_expression', colnames(x_2020_2022))])
  x_2020_2021[, grepl('gene_expression', colnames(x_2020_2021))] <- log(x_2020_2021[, grepl('gene_expression', colnames(x_2020_2021))])
  x_2022[, grepl('gene_expression', colnames(x_2022))] <- log(x_2022[, grepl('gene_expression', colnames(x_2022))])
  x_2023[, grepl('gene_expression', colnames(x_2023))] <- log(x_2023[, grepl('gene_expression', colnames(x_2023))])
}

if (task %in% c('1.1', '2.1', '3.1')){
  y_2020_2022 <- log(y_2020_2022)
  y_2020_2021 <- log(y_2020_2021)
  y_2022 <- log(y_2022)
}

x_col <- colnames(x_2020_2022)

output <- data.frame(SubjectID = rownames(x_2023), Age = floor(x_2023[, 'age (subject_info)']),
                     BiologicalSexAtBirth = ifelse(x_2023[, 'biological_sex (subject_info)_Male'] == 1, 'Male', 'Female'),
                     VaccinePrimingStatus = ifelse(x_2023[, 'infancy_vac (subject_info)_wP'] == 1, 'wP', 'aP'),
                     'IgG-PT-D14-titer-Rank' = 0,
                     'IgG-PT-D14-FC-Rank' = 0,
                     'Monocytes-D1-Rank' = 0,
                     'Monocytes-D1-FC-Rank' = 0,
                     'CCL3-D3-Rank' = 0, 
                     'CCL3-D3-FC-Rank' = 0,
                     check.names = FALSE)



saveRDS(output, "C:/CMI-PB/data/output template.RDS")
