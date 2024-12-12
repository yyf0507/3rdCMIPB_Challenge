A. Converting data to the format: columns <- variables (var_name, day, omics_type), rows <- subject_id

B. Data imputation
  There are 4 types of data imputations:
    1. Drop variables that contain NA;
    2. Drop subjects that contain NA
    3. Select subdata (eg: ab_titer on day 0) that contain less than 50% NA, use softimpute to impute those subdata.
    4. Select subdata that contain less than 50% NA, use softimpute to impute the whole data.

C. Prediction
  gene_expression ---determine---> cell frequency ---determine---> ab_titer, so the relationships between dependent variable and independent variable are

  y                  : x
  ab_titer prediction: all
  cell_frequency     : all except for ab_titer
  gene_expression    : all except ab_titer and cell_frequency

  We built 2 models for each task:
  1. Baseline model: use LASSO to select variables for each omic_types; combine selected variables with demographic information; use L1 penalized linear regression to do prediction (demographic information are L1-penalized). 
  2. MOFA: calculate 6 factors to summarize all omics data; combine factors with demographic data; use L1 penalized linear regression to do prediction (demographic information are L1-penalized).

  note: for task x.1, append corresponding y on day 0 as extra selected variable/factor;
        for task x.2, append corresponding 1/y on day 0 as extra selected variable/factor;
        log-transformed x that are positive; 
        log-transformed y for task x.1
