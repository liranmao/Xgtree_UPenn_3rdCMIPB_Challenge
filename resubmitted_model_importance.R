###### importance attribution
# SHAP values
library(SHAPforxgboost)

# Create lists to store importance results
importance_caret <- list()
importance_xgb <- list()
shap_values <- list()

###### models come from our xgboost training step in 
# for (task_col in colnames(task_mat_imputed)[-1]) {
#   # Prepare the target variable
#   target <- task_mat_imputed[[task_col]]
#   
#   model <- train(
#     x = input_mat_imputed,
#     y = target,
#     method = "xgbTree",         # XGBoost Tree
#     metric = "RMSE",            # Use RMSE as the evaluation metric
#     tuneLength = 10             # Let caret tune hyperparameters with 10 candidate models
#   )
#   
#   # Store the model
#   models[[task_col]] <- model
#   
#   # Predict on test data
#   preds <- predict(model, newdata = input_mat_test_imputed)
#   
#   # Store predictions
#   predictions[[task_col]] <- preds
# }


# Loop through models
for (task_col in colnames(task_mat_imputed)[-1]) {
  # Caret importance
  importance_caret[[task_col]] <- varImp(models[[task_col]])
  
  # XGBoost importance
  xgb_model <- models[[task_col]]$finalModel
  importance_xgb[[task_col]] <- xgb.importance(model = xgb_model)
  
  # Convert input data to matrix
  X_mat <- as.matrix(input_mat_imputed)
  
  # Calculate SHAP values
  shap_values[[task_col]] <- shap.values(
    xgb_model = xgb_model,
    X_train = X_mat
  )$shap_score
}

# Plot all results
library(gridExtra)
library(ggplot2)

# Reshape data to have one importance score per task
prepare_task_importance <- function(importance_caret, importance_xgb, shap_values) {
  tasks <- names(importance_caret)
  
  all_dfs <- lapply(tasks, function(task) {
    # Normalize each method's scores
    caret_scores <- importance_caret[[task]]$importance$Overall
    xgb_scores <- importance_xgb[[task]]$Gain * 100
    shap_scores <- colMeans(abs(shap_values[[task]])) * 100
    
    # Create dataframes with normalized scores
    caret_imp <- data.frame(
      Feature = rownames(importance_caret[[task]]$importance),
      Importance = caret_scores / sum(caret_scores),
      Method = "Caret",
      Task = task
    )
    
    xgb_imp <- data.frame(
      Feature = importance_xgb[[task]]$Feature,
      Importance = xgb_scores / sum(xgb_scores),
      Method = "XGBoost",
      Task = task
    )
    
    shap_imp <- data.frame(
      Feature = colnames(shap_values[[task]]),
      Importance = shap_scores / sum(shap_scores),
      Method = "SHAP", 
      Task = task
    )
    
    rbind(caret_imp, xgb_imp, shap_imp)
  })
  
  do.call(rbind, all_dfs)
}

importance_df <- prepare_task_importance(importance_caret, importance_xgb, shap_values)


# Create mapping for task names
task_mapping <- c(
  "IgG_PT_D14" = "1.1) IgG-PT-D14-titer-Rank",
  "IgG_PT_D14_FC" = "1.2) IgG-PT-D14-FC-Rank",
  "Monocytes_D1" = "2.1) Monocytes-D1-Rank",
  "Monocytes_D1_FC" = "2.2) Monocytes-D1-FC-Rank", 
  "CCL3_D3" = "3.1) CCL3-D3-Rank",
  "CCL3_D3_FC" = "3.2) CCL3-D3-FC-Rank"
)

importance_df$Task <- factor(task_mapping[importance_df$Task], 
                             levels = unname(task_mapping))

ggplot(importance_df, aes(x=Feature, y=Importance*100, fill=Method)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Task) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Normalized Importance Score (%)", x = "Feature") +
  scale_fill_brewer(palette = "Set2")
