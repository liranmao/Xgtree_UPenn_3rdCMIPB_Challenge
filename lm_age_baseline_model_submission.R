set.seed(123)

#################
# preprocess test datq
#################
# Convert categorical variables to factors
input_mat_test$biological_sex <- as.factor(input_mat_test$biological_sex)
input_mat_test$infancy_vac <- as.factor(input_mat_test$infancy_vac)

colSums(is.na(input_mat_test))

# Use model.matrix to create dummy variables for categorical data
input_mat_test_processed <- model.matrix(~ . -1, data = input_mat_test)
input_mat_test_processed <- as.data.frame(input_mat_test_processed)
input_mat_test_imputed <- input_mat_test_processed[, -1]


# Initialize lists to store models and predictions
models <- list()
predictions <- list()

# Loop through each task in task_mat
for (task_col in colnames(task_mat_imputed)[-1]) {
  # Prepare the target variable
  target <- task_mat_imputed[[task_col]]
  
  model <- train(
    x = input_mat_imputed,
    y = target,
    method = "xgbTree",         # XGBoost Tree
    metric = "RMSE",            # Use RMSE as the evaluation metric
    tuneLength = 10             # Let caret tune hyperparameters with 10 candidate models
  )
  
  # Store the model
  models[[task_col]] <- model
  
  # Predict on test data
  preds <- predict(model, newdata = input_mat_test_imputed)
  
  # Store predictions
  predictions[[task_col]] <- preds
}

# Combine predictions into a data frame
predictions_df <- as.data.frame(predictions)
predictions_df <- cbind(subject_id = input_mat_test_processed$subject_id, predictions_df)

# Save predictions to CSV
write.csv(predictions_df, "xgboost_predictions_no_cv.csv", row.names = FALSE)


############# prepare for submission
# Extract subject IDs from the input matrix
input_subject_ids <- input_mat_test$subject_id

# Step 1: Identify missing subjects
missing_subjects <- setdiff(input_subject_ids, predictions_df$subject_id)

# Step 2: Create new rows for missing subjects with NA predictions
missing_rows <- data.frame(
  subject_id = missing_subjects,
  Monocytes_D1 = NA,
  CCL3_D3 = NA,
  IgG_PT_D14 = NA,
  Monocytes_D1_FC = NA,
  CCL3_D3_FC = NA,
  IgG_PT_D14_FC = NA
)

# Step 3: Append missing rows to predictions_df
predictions_df_extended <- bind_rows(predictions_df, missing_rows)

# Step 4: Impute missing values using median
predictions_df_imputed <- predictions_df_extended %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# Step 5: Rank the prediction results for each task
# Get the number of rows (subjects)
n_subjects <- nrow(predictions_df_imputed)

# Initialize a dataframe for ranks with the same structure as the numeric columns of predictions_df
rank_predict <- as.data.frame(matrix(NA, nrow = n_subjects, ncol = ncol(predictions_df_extended) - 1))
colnames(rank_predict) <- colnames(predictions_df_extended)[-1]

# Compute ranks while ensuring no overlapping
for (i in 1:ncol(rank_predict)) {
  # Identify non-NA and NA indices
  non_na_indices <- which(!is.na(predictions_df_extended[, i + 1]))
  na_indices <- which(is.na(predictions_df_extended[, i + 1]))
  
  # Total ranks
  total_ranks <- seq(1, n_subjects)
  
  # Assign sequential median-based ranks to missing values
  n_missing <- length(na_indices)
  if (n_missing > 0) {
    start_rank <- floor(median(total_ranks)) - floor(n_missing / 2)
    missing_ranks <- seq(start_rank, length.out = n_missing)
    rank_predict[na_indices, i] <- missing_ranks
  }
  
  # Assign remaining ranks to non-NA values based on their order
  remaining_ranks <- setdiff(total_ranks, rank_predict[na_indices, i])
  rank_predict[non_na_indices, i] <- rank(-predictions_df_extended[non_na_indices, i + 1], ties.method = "min")
  
  # Ensure non-NA values get remaining ranks
  rank_predict[non_na_indices, i] <- remaining_ranks[rank(-predictions_df_imputed[non_na_indices, i + 1], ties.method = "min")]
}

# Add the subject_id column back to the rank dataframe
rank_predict <- cbind(subject_id = predictions_df_imputed$subject_id, rank_predict)

rownames(rank_predict) <- rank_predict$subject_id
rank_predict <- rank_predict[match(input_mat_test$subject_id, rank_predict$subject_id), ]
all(rank_predict$SubjectID == input_mat_test$subject_id)
rank_predict$Age <- input_mat_test$age_at_boost
rank_predict$VaccinePrimingStatus <- input_mat_test$infancy_vac
rank_predict$BiologicalSexAtBirth <- input_mat_test$biological_sex


colnames(rank_predict)[colnames(rank_predict)=="Monocytes_D1"] <- "2.1) Monocytes-D1-Rank"
colnames(rank_predict)[colnames(rank_predict)=="Monocytes_D1_FC"] <- "2.2) Monocytes-D1-FC-Rank" 
colnames(rank_predict)[colnames(rank_predict)=="CCL3_D3"] <- "3.1) CCL3-D3-Rank" 
colnames(rank_predict)[colnames(rank_predict)=="CCL3_D3_FC"] <- "3.2) CCL3-D3-FC-Rank" 
colnames(rank_predict)[colnames(rank_predict)=="IgG_PT_D14"] <- "1.1) IgG-PT-D14-titer-Rank" 
colnames(rank_predict)[colnames(rank_predict)=="IgG_PT_D14_FC"] <- "1.2) IgG-PT-D14-FC-Rank"    

colnames(rank_predict)[1] <- 'SubjectID'
rank_predict_final <- rank_predict[, c("SubjectID", "Age", "BiologicalSexAtBirth", "VaccinePrimingStatus", "1.1) IgG-PT-D14-titer-Rank" , "1.2) IgG-PT-D14-FC-Rank","2.1) Monocytes-D1-Rank", "2.2) Monocytes-D1-FC-Rank" ,"3.1) CCL3-D3-Rank" , "3.2) CCL3-D3-FC-Rank" )]

rank_predict_final$`4.1) IFNG/IL5-Polarization-D30-Rank` <- sample(1:nrow(rank_predict_final))

example <-  `3rdChallengeSubmissionTemplate_10032024.(2)`
rank_predict_final$Age <- example$Age
write.table(rank_predict_final,'xgtree_2023_preds_11212024.tsv',sep='\t', row.names = FALSE)



######################
# combine with lexi's model to get a best one
######################
lexi <- read.delim("~/Documents/1_courses_upenn/CMI-PB/3_compare_all_models/3rdChallengeSubmission_LiMaoNguyenTan.tsv")

# > print(best_models)
# # A tibble: 6 × 4
# Task            Model        Correlation Source
# <chr>           <chr>              <dbl> <chr> 
#   1 CCL3_D3         lexi               0.710 lexi  
# 2 CCL3_D3_FC      XGBoost Tree       0.610 lm    
# 3 IgG_PT_D14      XGBoost Tree       0.571 lm    
# 4 IgG_PT_D14_FC   XGBoost Tree       0.912 lm    
# 5 Monocytes_D1    lexi               0.751 lexi  
# 6 Monocytes_D1_FC lexi              -0.549 lexi 

all(rank_predict_final$SubjectID == lexi$SubjectID)

rank_predict_final$`1.1) IgG-PT-D14-titer-Rank` <- lexi$X1.1..IgG.PT.D14.titer.Rank
rank_predict_final$`2.1) Monocytes-D1-Rank` <- lexi$X2.1..Monocytes.D1.Rank
rank_predict_final$`2.2) Monocytes-D1-FC-Rank` <- lexi$X2.2..Monocytes.D1.FC.Rank

rank_predict_lm <- rank_predict[, c("SubjectID", "Age", "BiologicalSexAtBirth", "VaccinePrimingStatus", "1.1) IgG-PT-D14-titer-Rank" , "1.2) IgG-PT-D14-FC-Rank","2.1) Monocytes-D1-Rank", "2.2) Monocytes-D1-FC-Rank" ,"3.1) CCL3-D3-Rank" , "3.2) CCL3-D3-FC-Rank" )]

for (i in 5:ncol(rank_predict_lm)){
  print(colnames(rank_predict_lm)[i])
  print(cor(rank_predict_lm[,i], lexi[,i]))
}

rownames(predictions_df_imputed) <- predictions_df_imputed$subject_id
rank_predict_lm_ordered <- rank_predict_lm[match(predictions_df_imputed$subject_id, rank_predict_lm$SubjectID), ]

all(predictions_df_imputed$subject_id == rank_predict_lm_ordered$SubjectID)

cor(rank_predict_lm_ordered$`1.1) IgG-PT-D14-titer-Rank`, predictions_df_imputed$IgG_PT_D14)




write.table(rank_predict_final,'../3_compare_all_models/combineall_2023_preds_11212024.tsv',sep='\t', row.names = FALSE)
