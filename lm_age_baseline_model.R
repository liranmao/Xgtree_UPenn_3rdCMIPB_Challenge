library(verification)
library(pROC)
library(discretefit)
library(dplyr)
library(tidyverse)
library(corrplot)
library(data.table)
library(readr)
library(kableExtra)
library(glmnet)
library(here)
library(formattable)
library(ComplexHeatmap)
library(GetoptLong)
library(ellipse)
library(GGally)
library(caret)
library(randomForest)
library(e1071)
library(gbm)
library(xgboost)
source('https://raw.githubusercontent.com/akonstodata/mcia_mbpca/main/R/MCIA_mbpca_extra.R')

setwd("~/Documents/1_courses_upenn/CMI-PB/data_yt_code/")
################### 
# define variables
###################
df_source <- readRDS("master_allData_batchCorrected.RDS")

subject_2020 = read_tsv("2020LD_subject.tsv") 
subject_2021 = read_tsv("2021LD_subject.tsv") 
subject_2022 = read_tsv("2022LD_subject.tsv") 
subject_2023 = read_tsv("2023BD_subject.tsv")
subject = bind_rows(subject_2020, subject_2021, subject_2022, subject_2023)

DATASET <- c("2020_dataset", "2021_dataset", "2022_dataset")
TEST_DATASET = "2023_dataset"
TIMEPOINTS <- c(0, 1, 3, 14)

META_COLS <- c("specimen_id", "subject_id", "timepoint", "dataset", 
               "biological_sex", "infancy_vac", "age_at_boost")
ABTITER_COLS <- c("IgG_PT")
RNA_COLS <- c("CCL3")
CELL_COLS <- c("Monocytes")

DEMOGRAPHY_COLS <- c("age_at_boost", "biological_sex", "infancy_vac")

TASK_COLS <- c("Monocytes_D1", "CCL3_D3", "IgG_PT_D14")
TASKS_BASELINES <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
BASE_COLS <- c("Monocytes_D0", "IgG_PT_D0", "CCL3_D0")

LOG_TRANS_COLS <- c("CCL3_D0", "IgG_PT_D0",
                    "CCL3_D3",  "IgG_PT_D14")

################### 
# define variables
###################
# Obtain task df
metaDf <- data.frame(df_source[["subject_specimen"]]) %>% merge(subject)
metaDf["age_at_boost"] <- as.numeric(round(difftime(metaDf$date_of_boost, metaDf$year_of_birth,units="weeks")/52, 2))
metaDf_sel <- metaDf[, META_COLS] %>%
  data.frame()

abtiterDf <- df_source[["plasma_ab_titer"]]$batchCorrected_data %>%
  t() %>%
  data.frame() 

abtiterDf$specimen_id <- as.numeric(rownames(abtiterDf))
abtiterDf_sel <- data.frame(abtiterDf[, c("specimen_id", ABTITER_COLS)])

# Add log2 (11/10/2024)
rnaDf <- df_source[["pbmc_gene_expression"]]$tpm$batchCorrected_data %>%
  t() %>%
  data.frame() %>%
  mutate_all(~ log2(. + 1))

rnaDf$specimen_id <- as.numeric(rownames(rnaDf))
tasks_seq <- c('ENSG00000277632')
for (i in 1:length(tasks_seq)){
  rnaDf_sel <- data.frame(rnaDf %>% rename_at(vars(starts_with(tasks_seq[i])), ~RNA_COLS[i]))
}
rnaDf_sel <- data.frame(rnaDf_sel[, c("specimen_id", RNA_COLS)])


cellDf <- df_source[["pbmc_cell_frequency"]]$batchCorrected_data %>%
  t() %>%
  data.frame()  

cellDf$specimen_id <- as.numeric(rownames(cellDf))
cellDf_sel <- data.frame(cellDf[, c("specimen_id", CELL_COLS)])

list_df <- list(metaDf_sel, cellDf_sel, abtiterDf_sel, rnaDf_sel)
df_merge <- list_df %>% reduce(full_join, by="specimen_id")
df_merge <- df_merge[df_merge$timepoint %in% TIMEPOINTS, ]

df_pivot <- df_merge[, names(df_merge)!="specimen_id"] %>%
  pivot_wider(id_cols=c("subject_id", "dataset", "biological_sex",
                        "infancy_vac", "age_at_boost"),
              names_from = timepoint,
              values_from = all_of(c(CELL_COLS, RNA_COLS, ABTITER_COLS)),
              names_sep = "_D")

# df_pivot <- df_pivot[df_pivot$dataset %in% DATASET, ]
df_pivot <- data.frame(df_pivot %>%
                         mutate(across(everything(),  ~ case_when(.x >=0 ~ .x))))



targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")

fc_cols <- paste(targetY, "FC", sep="_")

# ranked_cols <- paste(c("age_at_boost", targetX, targetY, fc_cols), "Rank", sep="_")

df <- df_pivot[, c("subject_id", "dataset", DEMOGRAPHY_COLS, targetX, targetY)]

# rankingFunction <- function(x) {
#   as.numeric(rank(-x, ties.method = "min", na.last = "keep"))
# }

targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")

# Add log2 (11/10/2024)
df[,"Monocytes_D1_FC"] <- log2(df[, "Monocytes_D1"] / df[, "Monocytes_D0"])
df[,"CCL3_D3_FC"] <- log2(df[, "CCL3_D3"] / df[, "CCL3_D0"])
df[,"IgG_PT_D14_FC"] <- log2(df[, "IgG_PT_D14"] / df[, "IgG_PT_D0"])

# df[, ranked_cols] <- apply(df[, c("age_at_boost", targetX, targetY, fc_cols)],
#                            2, rankingFunction)
df <- data.frame(df)


# # obtain all time points for training data
# specimen_t0_train = df %>% filter(dataset %in% DATASET) %>% select(subject_id,dataset,age_at_boost,biological_sex,infancy_vac,Monocytes_D0, CCL3_D0, IgG_PT_D0)
# 
# # obtain all time points for test
# specimen_t0_test = df %>% filter(dataset == TEST_DATASET) %>% select(subject_id,dataset,age_at_boost,biological_sex,infancy_vac,Monocytes_D0, CCL3_D0, IgG_PT_D0)


tasks_ab<-c('IgG.PT','IgG.FHA','IgG.PRN','IgG1.PT','IgG1.FHA','IgG1.PRN','IgG4.PT','IgG4.FHA','IgG4.PRN')
tasks_cytof<-c('Monocytes','ASCs..Plasmablasts.','CD4Tcells')
tasks_seq<-c('ENSG00000277632','ENSG00000136244','ENSG00000100906','ENSG00000229807')
inputX <- c('age_at_boost', 'biological_sex', 'infancy_vac', 'Monocytes_D0',  'CCL3_D0',  'IgG_PT_D0')

# Add metadata
rownames(df) = df$subject_id
task_mat = df[df$dataset  %in% DATASET, c('subject_id', targetY, fc_cols)]
input_mat <- df[df$dataset  %in% DATASET, c('subject_id', inputX)]
input_mat_test <- df[df$dataset  == TEST_DATASET, c('subject_id', inputX)]

dim(task_mat)
dim(input_mat)
dim(input_mat_test)

colSums(is.na(task_mat))
colSums(is.na(input_mat))
colSums(is.na(input_mat_test))

#################
# preprocess
#################
# Convert categorical variables to factors
input_mat$biological_sex <- as.factor(input_mat$biological_sex)
input_mat$infancy_vac <- as.factor(input_mat$infancy_vac)

colSums(is.na(input_mat))

# Use model.matrix to create dummy variables for categorical data
input_mat_processed <- model.matrix(~ . -1, data = input_mat)
input_mat_processed <- as.data.frame(input_mat_processed)
input_mat_imputed <- input_mat_processed[, -1]

# Check for missing values
colSums(is.na(input_mat_processed))
dim(input_mat_processed)

# Repeat the same steps for task_mat
colSums(is.na(task_mat))
task_mat_imputed <- na.omit(task_mat)
task_mat_imputed <-  task_mat_imputed[match(input_mat_processed$subject_id, task_mat_imputed$subject_id), ]

all(input_mat_processed$subject_id == task_mat_imputed$subject_id)


#########################
# Define LOOCV train control with predictions saved
train_control <- trainControl(method = "LOOCV", savePredictions = 'final')

# Define a list of models with their methods and tuning grids
model_list <- list(
  list(name = 'Linear Regression', method = 'lm', tuneGrid = NULL),
  list(name = 'Lasso Regression', method = 'glmnet', tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 10))),
  list(name = 'Ridge Regression', method = 'glmnet', tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 0.1, length = 10))),
  list(name = 'Elastic Net', method = 'glmnet', tuneGrid = expand.grid(alpha = seq(0.1, 0.9, length = 5), lambda = seq(0.001, 0.1, length = 10))),
  list(name = 'Random Forest', method = 'rf', tuneGrid = NULL),
  list(name = 'SVM Linear', method = 'svmLinear', tuneGrid = NULL),
  list(name = 'SVM Radial', method = 'svmRadial', tuneGrid = NULL),
  list(name = 'Gradient Boosting', method = 'gbm', tuneGrid = NULL),
  list(name = 'XGBoost Linear', method = 'xgbLinear', tuneGrid = NULL),
  list(name = 'XGBoost Tree', method = 'xgbTree', tuneGrid = NULL)
)

# Initialize results data frame
results_df <- data.frame()

# Loop through each task
for (i in 2:ncol(task_mat_imputed)) {
  y <- task_mat_imputed[, i]
  task_name <- colnames(task_mat_imputed)[i]
  
  # Loop through each model
  for (model_info in model_list) {
    model_name <- model_info$name
    model_method <- model_info$method
    model_tuneGrid <- model_info$tuneGrid
    
    set.seed(123)
    # Try to train the model, handle exceptions
    try({
      model <- train(
        x = input_mat_imputed, y = y,
        method = model_method,
        trControl = train_control,
        tuneGrid = model_tuneGrid,
        preProcess = NULL,
        metric = 'RMSE'
      )
      
      # Get predictions and observations
      preds <- model$pred$pred
      obs <- model$pred$obs
      
      # Align predictions and observations if necessary
      if (nrow(model$pred) != length(y)) {
        preds_ordered <- preds[order(model$pred$rowIndex)]
        obs_ordered <- obs[order(model$pred$rowIndex)]
      } else {
        preds_ordered <- preds
        obs_ordered <- obs
      }
      
      # Compute correlation
      cor_value <- cor(preds_ordered, obs_ordered)
      
      # Store the results
      temp_result <- data.frame(Task = task_name, Model = model_name, Correlation = cor_value)
      results_df <- rbind(results_df, temp_result)
    }, silent = TRUE)
  }
}



# Print the results
print(results_df)

write.csv(results_df,"../combine_all_models/results_df_lm.csv", row.names = FALSE)


# Visualize the results using ggplot2
ggplot(results_df, aes(x = Model, y = Correlation, fill = Task)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = 'Model Performance (Correlation) by Task',
    y = 'Correlation between Predicted and True Values',
    x = 'Model'
  )

# Summarize results by model to find the best performing models
summary_df <- results_df %>%
  group_by(Model) %>%
  summarise(Mean_Correlation = mean(Correlation, na.rm = TRUE)) %>%
  arrange(desc(Mean_Correlation))

print(summary_df)

